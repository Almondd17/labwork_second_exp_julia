using CSV
using DataFrames
using Statistics
using Printf
using Plots

# Constants & instrument uncertainties
const g = 9.81
const M = 0.01662              # kg
const delta_M = 0.00001        # kg (0.01 g)
const h_s = 0.225              # m
const delta_h_s = 0.0005       # m (0.5 mm)
const delta_h_b = 0.0005       # m (0.5 mm)
const delta_D = 0.0005         # m (0.5 mm)

const r_ball = 0.0070          # m
const delta_r_ball = 0.0005    # m
const w = 0.0090               # m (track groove width)
const delta_w = 0.0005         # m

df = CSV.read("data.csv", DataFrame)

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

# 9. Plot Ek^(cm) as a function of h_b (Graph 1)
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

# 10. Save results to CSV for easy inspection and Graph 2
results_df = DataFrame(
    h_b_m = h_b,
    D_sq_avg_m2 = D_sq_avg,
    delta_D_sq_m2 = delta_D_sq,
    Ek_J = Ek,
    delta_Ek_J = delta_Ek
)
CSV.write("ek_results.csv", results_df)