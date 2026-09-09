using Plots, CSV, DataFrames
include("CRTBP_dynamics.jl")

# Coordinates for L1 and L2
const L1 = 0.836915
const L2 = 1.155682
const dt = 0.00001

orbit_data = CSV.read("periodic_orbits.csv", DataFrame)

"""
Uses Euler algorithm.\n
- f is the function\n
- y0 is the initial state vector
"""
function ni(f, y0, tspan, mu, dt)
    # Initialize timestep values
    times = collect(tspan[1]:dt:tspan[2])
    N = length(times)
    # Initialize vector storage
    state_vectors = zeros(6, N)
    state_vectors[:, 1] = y0
    dy = zeros(6)

    for i in 1:(N-1)
        f(dy, state_vectors[:, i], mu, times[i]) # Evaluate current derivates
        state_vectors[:, i+1] = state_vectors[:, i] + dt*dy # Euler step
    end
    return state_vectors
end



# Plotting
plt = plot(
    title = "Cislunar Orbits",
    legend = :topright,
)

# Plot static points, using row one's mass ratio for reference
ref_mu = orbit_data[1, Symbol("Mass ratio")]
scatter!(plt, [ -ref_mu ], [0], [0], color=:blue, label="Earth")
scatter!(plt, [1 - ref_mu], [0], [0], color=:gray, label="Moon")
scatter!(plt, [L1], [0], [0], color=:red, label="L1")
scatter!(plt, [L2], [0], [0], color=:red, label="L2")

# Get orbit data from CSV
for row in eachrow(orbit_data)
    y0 = [
        row[Symbol("x0 (LU) ")],
        row[Symbol("y0 (LU) ")],
        row[Symbol("z0 (LU) ")],
        row[Symbol("vx0 (LU/TU) ")],
        row[Symbol("vy0 (LU/TU) ")],
        row[Symbol("vz0 (LU/TU) ")]
        ]

    local states = ni(CR3BP!, y0, (0.0, row[Symbol("Period (TU) ")]), row[Symbol("Mass ratio")], dt)

    plot!(plt, states[1, :], states[2, :], states[3, :], label=row[Symbol("Id ")], lw=1.5)
end

display(plt)