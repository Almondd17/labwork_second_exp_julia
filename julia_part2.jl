##the balls height hb is positioned at the cm
using DelimitedFiles
using Plots

# Read the CSV file, skipping the header row
raw, header = readdlm("data.csv", ',', Any; header=true)

# Extract columns 2 (h_b_cm) and 3 (D1_cm) as Float64
h_b = Float64.(raw[:, 2])
D1  = Float64.(raw[:, 3])

# Sort by release height h_b in ascending order
order = sortperm(h_b)
h_b = h_b[order]
D1  = D1[order]
D2 = D1.^2
# Plot D(h_b)
scatter(
    h_b,
    D2,##d ^ 2
    xlabel = "Release Height h_b [cm]",
    ylabel = "Horizontal Distance D [cm]",
    title = "Week 2 Preliminary Graph: D(h_b)",
    label = "Measured Data (D1)",
    markersize = 6,
    legend = :topleft,
    grid = true
)

savefig("preliminary_D_vs_hb.png")
display(current())