using Plots, CSV, DataFrames, DifferentialEquations
plotlyjs()
include("CRTBP_dynamics.jl")
orbit_data = CSV.read("periodic_orbits.csv", DataFrame)

# CONSTANTS
const L1 = 0.836915
const L2 = 1.155682
const dt = 0.0001
const PERIODS = 4


# PLOTTING
plt_orbit = plot(
    title = "Cislunar Orbits",
    xlabel = "X (LU)",
    ylabel = "Y (LU)",
    zlabel = "Z (LU)"
)
plt_jacobi = plot(
    title = "Jacobi Constant Over Time",
    xlabel = "Time (TU)",
    ylabel = "Jacobi Constant (LU^2/TU^2)"
)
# Plot static points, using row one's mass ratio (mu) for reference
ref_mu = orbit_data[1, Symbol("Mass ratio")]
#scatter!(plt, [ -ref_mu ], [0], [0], color=:blue, label="Earth")
scatter!(plt_orbit, [1 - ref_mu], [0], [0], color=:gray, label="Moon")
scatter!(plt_orbit, [L1], [0], [0], color=:red, label="L1")
scatter!(plt_orbit, [L2], [0], [0], color=:red, label="L2")

colors = palette(:tab10)


#NUMERICAL INTEGRATION
# Define function for calculating jacobi given a state vector and mu
function calc_jacobi(y::Vector{Float64}, mu::Float64)
    r1 = sqrt((y[1]+mu)^2 + y[2]^2 + y[3]^2) # Distance to the earth 
    r2 = sqrt((y[1]+mu-1)^2 + y[2]^2 + y[3]^2) # Distance to the moon 
    return (y[1]^2) + (y[2]^2) + (2*(1-mu)/r1) + (2*mu/r2) - ((y[4]^2)+(y[5]^2)+(y[6]^2))
end
# Definite function for doing euler steps
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
i=1
for row in eachrow(orbit_data)
    global i
    orbit_name = row[Symbol("Id ")]
    mu = row[Symbol("Mass ratio")]
    y0 = [
        row[Symbol("x0 (LU) ")],
        row[Symbol("y0 (LU) ")],
        row[Symbol("z0 (LU) ")],
        row[Symbol("vx0 (LU/TU) ")],
        row[Symbol("vy0 (LU/TU) ")],
        row[Symbol("vz0 (LU/TU) ")]
        ]

    # Numerical Integration and plotting
    local states = ni(CR3BP!, y0, (0.0, row[Symbol("Period (TU) ")]*PERIODS), mu, dt)
    plot!(plt_orbit, states[1, :], states[2, :], states[3, :], label=orbit_name, lw=1.5, color=colors[i])

    # Do Jacobian calculations and plot it
    #calculated_jacobi = [calc_jacobi(j, mu) for j in eachrow(states)] #BUG: states is a matrix, and I need calc_jacobi to take each row and do the calculation on each state vector
    JPL_J = row[Symbol("Jacobi constant (LU2/TU2) ")]
    #plot!(plt_jacobi,states.t,calculated_jacobi, lw=1.5, color=colors[i]) # How do I get the t time?
    hline!(plt_jacobi, [JPL_J], linestyle=:dash, lw=1, color=:black)
    i = i+1
end

#display(plot(plt_orbit, plt_jacobi, layout=(1,2),size=(1000,600)))