"""
This program models a rocket landing
in the vertical and horizontal dimensions
with nonlinear gravity using a = -mu/r^3
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

const max_iters = 5
const tol = 1e-4

# SET UP PROBLEM

# Initial Reference
rx_ref = range(rx0, stop=0, length=N)
ry_ref = range(ry0, stop=1, length=N)
vx_ref = range(vx0, stop=0, length=N)
vy_ref = range(vy0, stop=0, length=N)

x_ref = [rx_ref'; ry_ref'; vx_ref'; vy_ref']
x_star = copy(x_ref)
u_star = zeros(2, N-1)

for iter in 1:max_iters
     global x_star, u_star, x_ref
     local x, u, constraints, A, B, C

     x = Variable(4, N)  # State variable [rx; ry; vx; vy]
     u = Variable(2, N-1) # Control variable [ux; uy]

     constraints = Constraint[
          x[:, 1] == [rx0; ry0; vx0; vy0],
          x[:, N] == [0; 1; 0; 0]
     ]

     for k in 1:(N-1)
     # Reference positions
          rx_bar = x_ref[1, k]
          ry_bar = x_ref[2, k]
          r_bar = sqrt(rx_bar^2 + ry_bar^2)

          ag_x = g * rx_bar / (r_bar^3)
          ag_y = g * ry_bar / (r_bar^3)

          J11 = g * (r_bar^2 - 3.0 * rx_bar^2) / (r_bar^5)
          J12 = g * (-3.0 * rx_bar * ry_bar) / (r_bar^5)
          J21 = J12
          J22 = g * (r_bar^2 - 3.0 * ry_bar^2) / (r_bar^5)

          # Linearized matrices
          A = [1.0       0.0       dtau 0.0;
               0.0       1.0       0.0  dtau;
               J11*dtau  J12*dtau  1.0  0.0;
               J21*dtau  J22*dtau  0.0  1.0]

          B = [0.0   0.0;
               0.0   0.0;
               dtau  0.0;
               0.0   dtau]

          C = [0.0;
               0.0;
               (ag_x - (J11 * rx_bar + J12 * ry_bar)) * dtau;
               (ag_y - (J21 * rx_bar + J22 * ry_bar)) * dtau]

          push!(constraints, (rx_bar * x[1, k] + ry_bar * x[2, k]) / r_bar >= 1.0)
          push!(constraints, x[:, k+1] == A * x[:, k] + B * u[:, k] + C)
          push!(constraints, norm(u[:, k]) <= u_max)
     end

     solve!(minimize(sum([norm(u[:, k]) for k in 1:(N-1)]), constraints), ECOS.Optimizer; silent=true)

     x_star = copy(evaluate(x))
     u_star = copy(evaluate(u))

     change = norm(x_star - x_ref, Inf)
     println("Iteration ", iter, ": Max Change = ", change)

     x_ref = copy(x_star)

     if change < tol
          println("Converged in ", iter, " iterations.")
          break
     end
end

# PLOTTING

println("Total Thrust: ", sum([norm(u_star[:, k]) for k in 1:(N-1)]))

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