"""
In this code I will be modeling a rocket landing with thrust, gravity,
and aerodynamic drag.
"""

using Convex, ECOS, LinearAlgebra, Plots


# CONSTANTS

const m = 1000.0 # rocket mass kg
const g = 9.81 # SL gravity m/s^2
const rho = 1.225 # SL density kg/m^3
const S = 2.0 # reference area m^2
const CD = 0.5 # drag coefficient
const dt = 0.5 # time set duration s
const N = 33 # number of time steps
const T_min = 2000.0 # minimum engine thrust N
const T_max = 20000.0 # maximum engine thrust N
const start_altitude = 1300.0 # m
const start_velocity = -150.0 # m/s
const max_iters = 30


# FUNCTIONS

function get_linearized_drag(v,v_bar)
    # v^2 = v_bar^2 + 2v_bar(v-v_bar) its a first order taylor series
    return 0.5 * rho * S * CD * (2.0 * v * v_bar - v_bar^2)
end

function solve_subproblem(h_bar, v_bar, T_bar, tr_v, tr_T)
    h = Variable(N)
    v = Variable(N)
    T = Variable(N-1)
    nu = Variable(N-1)
    z = Variable(N-1, :Bin)

    constraints = Constraint[]
    push!(constraints, h[1]==start_altitude)
    push!(constraints, v[1]==start_velocity)
    push!(constraints, h[N]==0.0)
    push!(constraints, v[N]==0.0)

    for k in 1:(N-1)
        push!(constraints, h[k+1] == h[k] + v[k] * dt)

        acceleration_k = (T[k] - m*g + get_linearized_drag(v[k], v_bar[k])) / m
        push!(constraints, v[k+1] == v[k] + acceleration_k * dt + nu[k])

        push!(constraints, T_min * z[k] <= T[k])
        push!(constraints, T[k] <= T_max * z[k])
        push!(constraints, v[k] - v_bar[k] <= tr_v)
        push!(constraints, v_bar[k] - v[k] <= tr_v)
        push!(constraints, T[k] - T_bar[k] <= tr_T)
        push!(constraints, T_bar[k] - T[k] <= tr_T)
    end

    objective = sum(T) + 1e6 * norm(nu, 1) # weighted penalty on nu
    problem = minimize(objective, constraints)
    solve!(problem, ECOS.Optimizer; silent=true)
    return evaluate(h), evaluate(v), evaluate(T), evaluate(nu)
end

function run_SCP()
    # Initial guess that we linearly interpolate down
    h_bar = collect(range(start_altitude,stop=0,length=N))
    v_bar = collect(range(start_velocity,stop=0,length=N))
    T_bar = fill(m*g, N-1) # Initial guess that Thrust matches weight

    tr_v = 15
    tr_T = 5000
    previous_slack = Inf

    for iter in 1:max_iters
        h_star, v_star, T_star, nu_star, z_star = solve_subproblem(h_bar, v_bar, T_bar, tr_v, tr_T)
        total_slack = norm(nu_star, 1)

        println("Iter $iter: Total  Slack: ", total_slack)
        
        if total_slack < 1e-3
            println("Trajectory Converged")
            println("Landing Velcoity: ", v_star[end], " m/s")
            println("Landing Thrusts: ", T_star[end-3:end])

            p1 = plot(0:dt:(N-1)*dt, h_star, label="Altitude (m)", color=:blue)
            p2 = plot(0:dt:(N-1)*dt, v_star, label="Velocity (m/s)", color=:red)
            p3 = plot(0:dt:(N-2)*dt, T_star, label="Thrust (N)", color=:green)
            display(plot(p1, p2, p3, layout=(3,1)))
            return
        end

        if total_slack < previous_slack
            previous_slack = total_slack
            # update references for next iteration
            h_bar = h_star
            v_bar = v_star
            T_bar = T_star

            tr_v = min(tr_v * 1.2, 30.0)
            tr_T = min(tr_T * 1.2, 10000.0)
        else
            tr_v = max(tr_v * 0.5, 0.5)
            tr_T = max(tr_T * 0.5, 500.0)
        end
    end
end

run_SCP()

