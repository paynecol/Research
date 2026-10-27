"""
This program models a rocket landing
in the vertical and horizontal dimensions
with no drag and constant gravity
using nondimensional units DU and TU.
"""

using Convex, ECOS, LinearAlgebra, Plots

# CONSTANTS

const rx0 = 1 # Initial horizontal displacement DU
const ry0 = 1 # Initial vertical displacement DU
const vx0 = -1 # Initial horizontal velocity DU/TU
const vy0 = -1 # Initial vertical velocity DU/TU

const g = -1 # Acceleration due to gravity DU/TU^2
const u_max = 2 # Acceleration due to thrust DU/TU^2

const total_time = 4 # TU
const dtau = 0.01 # Time step TU
const N = Integer(total_time/dtau) # Total time steps

# SET UP PROBLEM

# r_x_dot(t) = v_x(t)
# r_y_dot(t) = v_y(t)
# v_x_dot(t) = u_x(t)
# v_y_dot(t) = u_y(t) + g

# [r_x_dot, r_y_dot, v_x_dot, v_y_dot]T = A@[r_x, r_y, v_x, v_y]T + B@[u_x, u_y]T + C

# Linearization

A = [1.0  0.0  dtau 0.0;
     0.0  1.0  0.0  dtau;
     0.0  0.0  1.0  0.0;
     0.0  0.0  0.0  1.0]
B = [0.0  0.0;
     0.0  0.0;
     dtau 0.0;
     0.0  dtau]
C = [0.0; 0.0; 0.0; g * dtau]

x = Variable(4, N) # State variable of rx, ry, vx, vy
u = Variable(2, N-1) # Control variable of ux and uy

constraints = Constraint[
    x[:, 1] == [rx0; ry0; vx0; vy0], # Start conditions
    x[:, N] == [0.0; 0.0; 0.0; 0.0] # End conditions
]

# SOLVE PROBLEM

for k in 1:(N-1)
    push!(constraints, x[2,k] >= 0.0) # Cannot go into ground
    push!(constraints, x[:, k+1] == A * x[:, k] + B * u[:, k] + C)
    # push!(constraints, u_min <= norm(u[:, k])) # norm is already above 0
    push!(constraints, norm(u[:, k]) <= u_max)
end

solve!(minimize(sum([norm(u[:, k]) for k in 1:(N-1)]), constraints), ECOS.Optimizer; silent=false)

# PLOTTING

x_star = evaluate(x)
u_star = evaluate(u)
println("Total Thrust: ",sum(u_star))

t_x = range(0, stop=(N-1)*dtau, length=N)
t_u = range(0, stop=(N-2)*dtau, length=N-1)

p1 = plot(t_x, x_star[1, :], label="r_x", ylabel="Position (DU)")
     plot!(p1, t_x, x_star[2, :], label="r_y")
p2 = plot(t_x, x_star[3, :], label="v_x", ylabel="Velocity (DU/TU)")
     plot!(p2, t_x, x_star[4, :], label="v_y")
p3 = plot(t_u, u_star[1, :], label="u_x", ylabel="Thrust (DU/TU^2)", xlabel="Time (TU)")
     plot!(p3, t_u, u_star[2, :], label="u_y")
p_final = plot(p1, p2, p3, layout=(3, 1), size=(800, 800))
display(p_final)