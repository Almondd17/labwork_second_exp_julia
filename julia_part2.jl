using CSV
using DataFrames
using Statistics
using Plots

# 1. Load data
df = CSV.read("data.csv", DataFrame)

# 2. Extract h_b and convert from cm to meters
h_b = df.h_b_cm ./ 100.0

# 3. Extract distance columns (D1_cm through D10), convert to meters, and square
d_cols = names(df, r"^D\d+")
D_matrix = df[:, d_cols] ./ 100.0
D_sq_matrix = D_matrix .^ 2

# 4. Compute row-wise average of D^2 (ignoring missing cells)
D_sq_avg = [mean(skipmissing(row)) for row in eachrow(D_sq_matrix)]

# 5. Sort data in ascending order of h_b
idx = sortperm(h_b)
h_b = h_b[idx]
D_sq_avg = D_sq_avg[idx]

# 6. Plot ⟨D²⟩ as a function of h_b
scatter(
    h_b,
    D_sq_avg,
    label  = "Average D²",
    xlabel = "h_b [m]",
    ylabel = "⟨D²⟩ [m²]",
    title  = "⟨D²⟩ as a function of h_b",
    legend = :topleft,
    grid   = true
)

savefig("D2_vs_hb.png")