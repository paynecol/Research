using Plots, CSV, DataFrames, DifferentialEquations, LinearAlgebra
plotlyjs()
orbit_data = CSV.read("periodic_orbits_testing.csv", DataFrame)

# CONSTANTS
const L1 = 0.836915
const L2 = 1.155682
const PERIODS = 1.5

# FUNCTIONS
function CR3BPH(x::Float64, y::Float64, z::Float64, mu::Float64)
    r1 = sqrt((x + mu)^2 + y^2 + z^2)
    r2 = sqrt((x + mu - 1)^2 + y^2 + z^2)
    
    Uxx = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*(x+mu)^2/r1^5 + 3*mu*(x+mu-1)^2/r2^5
    Uyy = 1 - (1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*y^2/r1^5      + 3*mu*y^2/r2^5
    Uzz = 0 - (1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*z^2/r1^5      + 3*mu*z^2/r2^5
    
    Uxy = 3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x+mu-1)*y/r2^5
    Uxz = 3*(1-mu)*(x+mu)*z/r1^5 + 3*mu*(x+mu-1)*z/r2^5
    Uyz = 3*(1-mu)*y*z/r1^5      + 3*mu*y*z/r2^5
    
    return [Uxx Uxy Uxz; Uxy Uyy Uyz; Uxz Uyz Uzz]
end

function CR3BP!(dy::Vector{Float64}, y::Vector{Float64}, mu::Float64, t::Float64)
    r1 = sqrt((y[1]+mu)^2 + y[2]^2 + y[3]^2)
    r2 = sqrt((y[1]+mu-1)^2 + y[2]^2 + y[3]^2)

    dy[1] = y[4]
    dy[2] = y[5]
    dy[3] = y[6]
    dy[4] = y[1]+2*y[5]-((1-mu)*(y[1]+mu))/(r1^3)-mu*(y[1]-(1-mu))/(r2^3)
    dy[5] = y[2]-2*y[4]-(1-mu)/(r1^3)*y[2]-mu*(y[2])/(r2^3)
    dy[6] = -(1-mu)/(r1^3)*y[3]-mu*(y[3])/(r2^3)

    # What is all this?
    Urr = CR3BPH(y[1], y[2], y[3], mu) # 2nd derivative of potential?
    Omega = [0.0 2.0 0.0; -2.0 0.0 0.0; 0.0 0.0 0.0] # For coriolis
    A = [zeros(3,3) I(3); Urr Omega] # Jacobian matrix of eqs of motion
    
    Phi = reshape(y[7:42], 6, 6)
    dPhi = A * Phi
    dy[7:42] .= vec(dPhi)
    
    nothing 
end

function calc_jacobi(y::AbstractVector{Float64}, mu::Float64)
    r1 = sqrt((y[1] + mu)^2 + y[2]^2 + y[3]^2)
    r2 = sqrt((y[1] + mu - 1)^2 + y[2]^2 + y[3]^2)
    return (y[1]^2 + y[2]^2) + (2*(1 - mu)/r1) + (2*mu/r2) - (y[4]^2 + y[5]^2 + y[6]^2)
end


# PLOTTING
plt_orbit = plot(
    title="Cislunar Orbits", 
    xlabel="X (LU)", 
    ylabel="Y (LU)", 
    zlabel="Z (LU)"
)
plt_jacobi = plot(
    title="Jacobi Constant Over Time", 
    xlabel="Time (TU)", 
    ylabel="Jacobi Constant (LU^2/TU^2)"
)
# Plot static points, using row one's mass ratio (mu) for reference
ref_mu = orbit_data[1, Symbol("Mass ratio")]
#scatter!(plt, [ -ref_mu ], [0], [0], color=:blue, label="Earth") #Earth is pretty far from the orbits and makes the graph look worse
scatter!(plt_orbit, [1 - ref_mu], [0], [0], color=:gray, label="Moon")
scatter!(plt_orbit, [L1], [0], [0], color=:red, label="L1")
scatter!(plt_orbit, [L2], [0], [0], color=:red, label="L2")


# NUMERICAL INTEGRATION
for row in eachrow(orbit_data)
    mu = row[Symbol("Mass ratio")]
    period = row[Symbol("Period (TU) ")]
    orbit_id = row[Symbol("Id ")]
    y0_state = [
        row[Symbol("x0 (LU) ")],
        row[Symbol("y0 (LU) ")],
        row[Symbol("z0 (LU) ")],
        row[Symbol("vx0 (LU/TU) ")],
        row[Symbol("vy0 (LU/TU) ")],
        row[Symbol("vz0 (LU/TU) ")]
    ]

    Phi0 = Matrix{Float64}(I, 6, 6)
    y0 = vcat(y0_state, vec(Phi0))

    prob = ODEProblem(CR3BP!, y0, (0.0, period * PERIODS), mu)
    sol = solve(prob, Vern7(), reltol=1e-12, abstol=1e-12)

    x_pts = [u[1] for u in sol.u]
    y_pts = [u[2] for u in sol.u]
    z_pts = [u[3] for u in sol.u]
    
    # Calculate Jacobi constant across time steps
    jacobi_vals = [calc_jacobi(u[1:6], mu) for u in sol.u]

    # Plot orbit trajectory and Jacobi verification
    plot!(plt_orbit, x_pts, y_pts, z_pts, label=orbit_id, lw=1.5)
    plot!(plt_jacobi, sol.t, jacobi_vals, label=orbit_id, lw=1.5)
end

display(plot(plt_orbit, plt_jacobi, layout=(1,2), size=(1000, 550)))