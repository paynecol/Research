"""
This program models a rocket landing
in the vertical and horizontal dimensions
with no drag and nonlinear gravity
using nondimensional units DU and TU.
"""

using Convex, ECOS, LinearAlgebra, Plots

# CONSTANTS

const r0 = 2 # Initial vertical displacement DU
const rf = 1 # Target surface (floor) radius DU
const v0 = -1 # Initial vertical velocity DU/TU

const g_coeff = -1 # Acceleration due to gravity DU/TU^2
const u_min = 0 # Acceleration due to thrust DU/TU^2
const u_max = 2 

const total_time = 4 # TU
const dtau = 0.01 # Time step TU
const N = Integer(total_time/dtau) # Total time steps

const max_iters = 10 # SCP iterations
const tol = 1e-4 # Convergence tolerance

# SCP

# r_dot(t) = v(t)
# v_dot(t) = u(t) + g/r(t)^2

r_ref = range(r0, stop=rf, length=N) # Initial guess: linear interpolation
v_ref = range(v0, stop=0, length=N)
x_ref = [r_ref'; v_ref'] # Apostrophe is transpose
x_star = copy(x_ref)
u_star = zeros(N-1)

for iter in 1:max_iters
     global x_star, u_star, x_ref
     local x, u, constraints

     x = Variable(2, N) # State variable of rx, ry, vx, vy
     u = Variable(N-1) # Control variable of ux and uy
     nu = Variable(2, N-1)

     constraints = Constraint[
     x[:, 1] == [r0; v0], # Start conditions
     x[:, N] == [rf; 0] # End conditions
     ]

     # SOLVE PROBLEM

     for k in 1:(N-1)
          local A, B, C
          r_bar = x_ref[1,k]

          # r_dot(t) = v(t)
          # v_dot(t) = u(t) + g/r(t)^2

          # Linearization by 1st-order taylor expansion of g(r) = g/r^2 around r_bar
          # g(r) = g/r^2
          # g(r_bar) = g/r_bar^2
          # g'(r_bar) = -2g/r_bar^3
          # g(r) = g(r_bar) + g'(r_bar)(r-r_bar)
          # v_dot = u + g'(r_bar)*r + (g_bar - g'r(bar)*r_bar)
          # [r_dot; v_dot] = [0 1; dg_dr 0]@[r, v] + [0; 1]u + [0; g_bar - dg_dr*r_bar]

          # Discretization
          # r_(k+1) = r_k + v_k*dtau
          # v_(k+1) = v_k + (dg_dr*r_k + u_k + g_bar-dg_dr*r_bar)*dtau

          g_bar = g_coeff/(r_bar^2)
          dg_dr = -2 * g_coeff / (r_bar^3)

          A = [1 dtau; dg_dr*dtau 1]
          B = [0; dtau]
          C = [0; (g_bar - dg_dr*r_bar)*dtau]

          push!(constraints, x[1,k] >= rf) # Cannot go into ground
          push!(constraints, x[:, k+1] == A * x[:, k] + B * u[k] + C + nu[:, k])
          push!(constraints, u_min <= u[k])
          push!(constraints, u[k] <= u_max)
     end

     solve!(minimize(sum(u) + 1e6 * norm(nu, 1), constraints), ECOS.Optimizer; silent=true)

     # CHECK CONVERGENCE
     x_star = copy(evaluate(x))
     u_star = evaluate(u)

     change = norm(x_star - x_ref, Inf)
     println("Iteration ", iter, ": Max Change = ", change)

     x_ref = copy(x_star)

     if change < tol
          println("Converged in ", iter, " iterations.")
          break
     end
end

# PLOTTING

println("Total Thrust: ",sum(u_star))


t_x = range(0, stop=(N-1)*dtau, length=N)
t_u = range(0, stop=(N-2)*dtau, length=N-1)

p1 = plot(t_x, x_star[1, :], ylabel="Radius (DU)", legend=false)
p2 = plot(t_x, x_star[2, :], ylabel="Velocity (DU/TU)", legend=false)
p3 = plot(t_u, u_star, ylabel="Thrust (DU/TU^2)", xlabel="Time (TU)", legend=false)

p_final = plot(p1, p2, p3, layout=(3, 1), size=(800, 800))
display(p_final)