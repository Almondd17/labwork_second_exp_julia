using CSV
using DataFrames
using Statistics
using Printf
using Plots

#consts instrument uncertainties
const g = 9.81
const M = 0.01662              # kg
const delta_M = 0.00001        # kg (0.01 g)
const h_s = 0.225              # m
const delta_h_s = 0.0005       # m (0.5 mm)
const delta_h_b = 0.0005       # m (0.5 mm)
const delta_D = 0.0005         # m (0.5 mm)

const r_ball = 0.0070          # m
const delta_r_ball = 0.0005    # m (0.5 mm)
const w = 0.0080               # m (track width)
const delta_w = 0.0005         # m (0.5 mm)

df = CSV.read(joinpath(@__DIR__, "data.csv"), DataFrame)

h_b = df.h_b_cm ./ 100.0 #get the h_b column in meters

d_cols = names(df, r"^D\d+") #array of all distances
D_matrix = Matrix(df[:,d_cols])./100.0 #convert distances to meters
D_sq_matrix = D_matrix .^ 2 ##square each value

##sort in ascending order
idx = sortperm(h_b)
h_b = h_b[idx]
D_matrix = D_matrix[idx, :]
D_sq_matrix = D_sq_matrix[idx, :]

# 5. Compute row-wise mean and total uncertainty of <D^2> (handling missing cells safely)
n_rows = length(h_b)
D_avg = zeros(n_rows)
D_sq_avg = zeros(n_rows)
delta_D_sq = zeros(n_rows)

for i in 1:n_rows
    valid_D = collect(skipmissing(D_matrix[i, :]))
    valid_D2 = collect(skipmissing(D_sq_matrix[i, :]))
    N = length(valid_D2)
    
    ##create average for each h_b
    D_avg[i] = mean(valid_D)
    D_sq_avg[i] = mean(valid_D2)
    
    ##deviation over sqrt(N) for the statistical error of the D^2 
    stat_error = std(valid_D2) / sqrt(N)
    ruler_error = 2.0 * D_avg[i] * delta_D ##ruler D^2 ERROR
    delta_D_sq[i] = sqrt(stat_error^2 + ruler_error^2) ##total error
end

#constant factor k = Mg/(4h_s) and its uncertainty
k = (M * g) / (4 * h_s)
delta_k = k * sqrt((delta_M / M)^2 + (delta_h_s / h_s)^2)

#Ek and its uncertainty for each h_b
Ek = k .* D_sq_avg
delta_Ek = sqrt.((D_sq_avg .* delta_k).^2 .+ (k .* delta_D_sq).^2)

println(@sprintf("Constant factor k = %.6f ± %.6f N/m\n", k, delta_k))
println("h_b [m]  |  <D^2> [m^2]       |  Ek^(cm) [J]")
println("-"^52)
for i in 1:n_rows ##print needed data
    println(@sprintf("%.4f   |  %.6f ± %.6f  |  %.6f ± %.6f",
        h_b[i], D_sq_avg[i], delta_D_sq[i], Ek[i], delta_Ek[i]))
end

#E_spin calculations: calculate r spin and its uncertainty
r_spin = sqrt(r_ball^2 - (w / 2)^2)
delta_r_spin = sqrt(((r_ball / r_spin) * delta_r_ball)^2 + ((-w / (4 * r_spin)) * delta_w)^2)
spin_factor = (2.0 / 5.0) * (r_ball / r_spin)^2 #constant factor for E_spin calculation

E_spin = spin_factor .* Ek##get each E_spin
delta_E_spin = E_spin .* sqrt.(##E spin error as calculated in docs
    (delta_r_ball / r_ball)^2 +
    (2.0 * delta_r_spin / r_spin)^2 .+
    (delta_Ek ./ Ek).^2
)

#data for the second graph: Ek/U vs U
U = M .* g .* h_b ##total energy U
##Ek/U results + error
delta_U = U .* sqrt.((delta_M / M)^2 .+ (delta_h_b ./ h_b).^2)
Ek_over_U = Ek ./ U
delta_Ek_over_U = Ek_over_U .* sqrt.((delta_Ek ./ Ek).^2 .+ (delta_U ./ U).^2)

spin_over_U = (E_spin ./ U) .* 100.0 ##get the spin/U percentage for the lab report
delta_spin_over_U = spin_over_U .* sqrt.(
    (delta_r_ball / r_ball)^2 +
    (2.0 * delta_r_spin / r_spin)^2 .+
    (delta_D_sq ./ D_sq_avg).^2 .+
    (delta_h_s / h_s)^2 .+
    (delta_h_b ./ h_b).^2
)

avg_spin_over_U = mean(spin_over_U)##e spin percentage average and error
delta_avg_spin_over_U = sqrt(sum(delta_spin_over_U .^ 2)) / length(spin_over_U)

println("\nh_b [m]  |  E_spin/U [%] ± uncertainty [%]")
println("-"^65)
for i in 1:n_rows
    println(h_b[i], "   |  ", spin_over_U[i], " ± ", delta_spin_over_U[i])
end

println("\nAverage E_spin/U [%]: ", avg_spin_over_U, " ± ", delta_avg_spin_over_U)

#Ek vs h_b
p1 = scatter(
    h_b,
    Ek,
    xerror = fill(delta_h_b, n_rows),
    yerror = delta_Ek,
    label  = "Experimental Eₖ⁽ᶜᵐ⁾",
    xlabel = "h_b [m]",
    ylabel = "Eₖ⁽ᶜᵐ⁾ [J]",
    title  = "Kinetic Energy Eₖ⁽ᶜᵐ⁾ vs Release height",
    legend = :topleft,
    markersize = 2,
    markerstrokewidth = 0.5,
    msc = :black,
    color = :dodgerblue,
    grid = true,
    dpi = 300,
)

display(p1)
savefig(p1, "Ek_vs_hb.png")

#Ek/U vs U
p2 = scatter(
    U,
    Ek_over_U,
    xerror = delta_U,
    yerror = delta_Ek_over_U,
    label  = "Experimental Eₖ/U",
    xlabel = "U [J]",
    ylabel = "Eₖ/U",
    title  = "Eₖ/U vs Potential Energy U",
    legend = :topleft,
    markersize = 2,
    markerstrokewidth = 0.5,
    msc = :black,
    color = :dodgerblue,
    grid = true,
    dpi = 300,
)

display(p2)
savefig(p2, joinpath(@__DIR__, "Ek_over_U_vs_U.png"))

#save results to csv
results_df = DataFrame(
    h_b_m = h_b,
    U_J = U,
    D_sq_avg_m2 = D_sq_avg,
    delta_D_sq_m2 = delta_D_sq,
    Ek_J = Ek,
    delta_Ek_J = delta_Ek,
    E_spin_J = E_spin,
    delta_E_spin_J = delta_E_spin,
    E_spin_over_U_pct = spin_over_U,
    delta_E_spin_over_U_pct = delta_spin_over_U
)

CSV.write("ek_results.csv", results_df)