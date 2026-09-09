using Plots, CSV, DataFrames, DifferentialEquations, LinearAlgebra
plotlyjs()

orbit_data = CSV.read("periodic_orbits.csv", DataFrame)

# CONSTANTS
const L1 = 0.836915
const L2 = 1.155682
const PERIODS = 1.5

# FUNCTIONS
function CR3BPH(x::Float64, y::Float64, z::Float64, mu::Float64)
    r1_sq = (x + mu)^2 + y^2 + z^2
    r2_sq = (x + mu - 1)^2 + y^2 + z^2
    r1_5 = r1_sq^(5/2)
    r2_5 = r2_sq^(5/2)
    
    Uxx = 1 - (1-mu)/r1_sq^(3/2) - mu/r2_sq^(3/2) + 3*(1-mu)*(x+mu)^2/r1_5 + 3*mu*(x+mu-1)^2/r2_5
    Uyy = 1 - (1-mu)/r1_sq^(3/2) - mu/r2_sq^(3/2) + 3*(1-mu)*y^2/r1_5 + 3*mu*y^2/r2_5
    Uzz = -(1-mu)/r1_sq^(3/2) - mu/r2_sq^(3/2) + 3*(1-mu)*z^2/r1_5 + 3*mu*z^2/r2_5
    
    Uxy = 3*(1-mu)*(x+mu)*y/r1_5 + 3*mu*(x+mu-1)*y/r2_5
    Uxz = 3*(1-mu)*(x+mu)*z/r1_5 + 3*mu*(x+mu-1)*z/r2_5
    Uyz = 3*(1-mu)*y*z/r1_5 + 3*mu*y*z/r2_5
    
    return [Uxx Uxy Uxz;
            Uxy Uyy Uyz;
            Uxz Uyz Uzz]
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

    U_rr = cr3bp_hessian(y[1], y[2], y[3], mu)
    Omega = [0.0 2.0 0.0; -2.0 0.0 0.0; 0.0 0.0 0.0]
    A = [zeros(3,3) I(3); U_rr Omega]
    
    Phi = reshape(y[7:42], 6, 6)
    dPhi = A * Phi
    dy[7:42] .= vec(dPhi)
    
    nothing 
end