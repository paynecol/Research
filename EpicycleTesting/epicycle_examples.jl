using Epicycle
using Printf

# Set the parking and final orbit radii
const R_PARK = 6878.137                   # km
const R_FINAL = 21000.0                   # km, about half geostationary
const MU_EARTH = 398600.4418              # km^3/s^2

# Configure the propagator
forces = ForceModel(PointMassGravity(earth, ()))
integ = IntegratorConfig(DP8();
                        dt = 60.0,
                        reltol = 1e-11,
                        abstol = 1e-11)
propagator = OrbitPropagator(forces, integ)

# Initialize a circular parking orbit at 28.5 degrees inclination
v_park = sqrt(MU_EARTH / R_PARK)
inclination_park = deg2rad(28.5)

sat = Spacecraft(state = CartesianState([R_PARK, 0.0, 0.0,
                                         0.0,
                                         v_park * cos(inclination_park),
                                         v_park * sin(inclination_park)]),
                 time = Time("2024-03-01T00:00:00.000", UTC(), ISOT()),
                 name = "Transfer Vehicle")

# Compute the two Hohmann burns in closed form
a_transfer = 0.5 * (R_PARK + R_FINAL)
dv_raise = v_park * (sqrt(2 * R_FINAL / (R_PARK + R_FINAL)) - 1)
dv_circularize = sqrt(MU_EARTH / R_FINAL) * (1 - sqrt(2 * R_PARK / (R_PARK + R_FINAL)))
transfer_time = π * sqrt(a_transfer^3 / MU_EARTH)

# Propagate one parking orbit
propagate!(propagator, sat,
           StopAt(sat, PropDurationSeconds(), 2π * sqrt(R_PARK^3 / MU_EARTH)))

# Raise apogee to the final radius
raise = ImpulsiveManeuver(axes = VNB(), element1 = dv_raise)
maneuver!(sat, raise)

# Coast the transfer to apogee
propagate!(propagator, sat, StopAt(sat, PropDurationSeconds(), transfer_time))

# Circularize
circularize = ImpulsiveManeuver(axes = VNB(), element1 = dv_circularize)
maneuver!(sat, circularize)

# Propagate one orbit at the final radius
propagate!(propagator, sat,
           StopAt(sat, PropDurationSeconds(), 2π * sqrt(R_FINAL^3 / MU_EARTH)))



# Draw the trajectory with about two minutes of playback
orbitview("Orbit raising", sat; speed = 300.0)

# Report the transfer against the closed-form values
@printf("apogee raise    : %.4f km/s\n", dv_raise)
@printf("circularization : %.4f km/s\n", dv_circularize)
@printf("transfer time   : %.2f hours\n", transfer_time / 3600)
@printf("final radius    : %.1f km   (target %.1f)\n",
        position_magnitude(sat), R_FINAL)
@printf("arcs in the view: %d\n", length(sat.history))