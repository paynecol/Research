using Plots, CSV, DataFrames, DifferentialEquations
plotlyjs()
include("CRTBP_dynamics.jl")
orbit_data = CSV.read("periodic_orbits.csv", DataFrame)

# CONSTANTS
const L1 = 0.836915
const L2 = 1.155682
const PERIODS = 5


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
#scatter!(plt, [ -ref_mu ], [0], [0], color=:blue, label="Earth") #Earth is pretty far from the orbits and makes the graph look worse
scatter!(plt_orbit, [1 - ref_mu], [0], [0], color=:gray, label="Moon")
scatter!(plt_orbit, [L1], [0], [0], color=:red, label="L1")
scatter!(plt_orbit, [L2], [0], [0], color=:red, label="L2")

colors = palette(:tab10)


# NUMERICAL INTEGRATION
# Define function for calculating jacobi given a state vector and mu
function calc_jacobi(y::Vector{Float64}, mu::Float64)
    r1 = sqrt((y[1]+mu)^2 + y[2]^2 + y[3]^2) # Distance to the earth 
    r2 = sqrt((y[1]+mu-1)^2 + y[2]^2 + y[3]^2) # Distance to the moon 
    return (y[1]^2) + (y[2]^2) + (2*(1-mu)/r1) + (2*mu/r2) - ((y[4]^2)+(y[5]^2)+(y[6]^2))
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

    # Numerical integration of each orbit using Vern7 (Which solver is best for this? What are good tolerances?)
    prob = ODEProblem(CR3BP!, y0, (0.0, row[Symbol("Period (TU) ")]*PERIODS), mu)
    local sol = solve(prob, Vern7(), reltol=1e-12, abstol=1e-12)
    plot!(plt_orbit, sol[1, :], sol[2, :], sol[3, :], label=orbit_name, lw=1.5, color=colors[i])

    # Do Jacobian calculations and plot it
    calculated_jacobi = [calc_jacobi(j, mu) for j in sol.u]
    JPL_J = row[Symbol("Jacobi constant (LU2/TU2) ")]
    plot!(plt_jacobi,sol.t,calculated_jacobi, lw=1.5, color=colors[i])
    hline!(plt_jacobi, [JPL_J], linestyle=:dash, lw=1, color=:black)
    i = i+1
end

display(plot(plt_orbit, plt_jacobi, layout=(1,2),size=(1000,600)))