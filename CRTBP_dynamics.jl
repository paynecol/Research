export CR3BP!
"""
- dy is the derivative of the state vector y\n
- y is the state vector [x,y,z,vx,vy,vz]\n
- mu is the mass ratio moon/(earth + moon)\n
- t isn't used here?\n
"""
function CR3BP!(dy::Vector{Float64}, y::Vector{Float64}, mu::Float64, t::Float64)
    #
    r1 = sqrt((y[1]+mu)^2 + y[2]^2 + y[3]^2) # Distance to the earth 
    r2 = sqrt((y[1]+mu-1)^2 + y[2]^2 + y[3]^2) # Distance to the moon 

    dy[1] = y[4]
    dy[2] = y[5]
    dy[3] = y[6]
    dy[4] = y[1]+2*y[5]-((1-mu)*(y[1]+mu))/(r1^3)-mu*(y[1]-(1-mu))/(r2^3)
    dy[5] = y[2]-2*y[4]-(1-mu)/(r1^3)*y[2]-mu*(y[2])/(r2^3)
    dy[6] = -(1-mu)/(r1^3)*y[3]-mu*(y[3])/(r2^3)
    #
    nothing 
    #
end