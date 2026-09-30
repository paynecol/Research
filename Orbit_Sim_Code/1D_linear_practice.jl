"""
This program models a rocket landing
in only the vertical dimension
with no drag and constant gravity
using nondimensional units DU and TU.
"""

using Convex, ECOS, LinearAlgebra, Plots

# CONSTANTS

const r0 = 1 # Initial altitude DU
const v0 = -1 # Initial velocity DU/TU
const g = -1 # Acceleration due to gravity DU/TU^2
const u_min = 0 # Acceleration due to thrust DU/TU^2
const u_max = 1

const total_time = 2 # TU
const dtau = 0.01 # Time step TU
const N = Integer(total_time/dtau) # Total time steps

# SET UP PROBLEM

# h_dot(t) = v(t)
# v_dot(t) = u(t) + g

# [h_dot, v_dot]T = [[0,1],[0,0]]@[h,v]T + [0,1]T*u + [0,g]T

# x_(k+1) = x_k + (A@x_k + Bu + C)dtau
# x_(k+1) = (I + A*dtau)@x_k + B*dtau*u + C*dtau

A = [1.0  dtau; 0.0  1.0]
B = [0.0; dtau]
C = [0.0; g * dtau]

x = Variable(2, N) # State variable of h and v
u = Variable(N-1)

constraints = Constraint[
    x[:, 1] == [r0; v0], # Start conditions
    x[:, N] == [0.0; 0.0] # End conditions
]

# SOLVE PROBLEM

for k in 1:(N-1)
    push!(constraints, x[1,k] >= 0.0) # Cannot go into ground
    push!(constraints, x[:, k+1] == A * x[:, k] + B * u[k] + C)
    push!(constraints, u_min <= u[k])
    push!(constraints, u[k] <= u_max)
end

solve!(minimize(sum(u), constraints), ECOS.Optimizer; silent=true)

# PLOTTING

x_star = evaluate(x)
u_star = evaluate(u)
println("Total Thrust: ",sum(u_star))

t_x = range(0, stop=(N-1)*dtau, length=N)
t_u = range(0, stop=(N-2)*dtau, length=N-1)

p1 = plot(t_x, x_star[1, :], ylabel="Altitude (DU)", legend=false)
p2 = plot(t_x, x_star[2, :], ylabel="Velocity (DU/TU)", legend=false)
p3 = plot(t_u, u_star, ylabel="Thrust (DU/TU^2)", xlabel="Time (TU)", legend=false)

p_final = plot(p1, p2, p3, layout=(3, 1), size=(800,800))
display(p_final)