import numpy as np

G = 9.81
WATER_DENSITY = 1000
AIR_DENSITY = 1.225


def calculate_buoyancy(V, density_fluid):
    """calculates buoyancy
    Takes V for volume displaced and density_fluid for the density of the fluid"""
    return V * density_fluid * G


def will_it_float(V, mass):
    """returns if an object will float
    Takes V for volume of object and mass for mass of the object"""
    return mass / V < WATER_DENSITY


def calculate_pressure(depth):
    """calculates pressure on an object
    Takes depth for depth of object"""
    if depth <= 0:
        return 101325 - AIR_DENSITY * G * depth
    else:
        return depth * G * WATER_DENSITY + 101325


def calculate_acceleration(F, m):
    """this calculates acceleration
    Takes F for force and m for mass"""
    return F / m


def calculate_angular_acceleration(tau, I):
    """calculates angular acceleration
    Takes tau for torque and I for inertia"""
    return tau / I


def calculate_torque(F_magnitude, F_direction, r):
    """calculates torque
    Takes F_magnitude for force magnitude, F_direction for force direction, and r for length of the arm
    """
    F_direction = F_direction * (np.pi / 180)
    return F_magnitude * r * np.sin(r)


def calculate_moment_of_inertia(m, r):
    """calculates moment of inertia
    Takes m for mass and r for radius"""
    return m * (r**2)


def calculate_auv_acceleration(
    F_magnitude, F_angle, mass=100, volume=0.1, thruster_distance=0.5
):
    """calculates acceleration of the AUV
    Takes F_magnitude for force magnitude, F_angle for force direction, mass for mass of AUV, volume for volume of AUV, and thruster_distance for the thruster distance of the AUV
    """
    return F_magnitude / mass


def calculate_auv_angular_acceleration(
    F_magnitude, F_angle, inertia=1, thruster_distance=0.5
):
    """calculates angular acceleration of the AUV
    Takes F_magnitude for force magnitude, F_angle for force direction, inertia for moment of inertia, and thruster_distance for the thruster distance of the AUV
    """
    return calculate_torque(F_magnitude, F_angle, thruster_distance) / inertia


def calculate_auv2_acceleration(T, alpha, theta, mass=100):
    """calculates acceleration of the AUV
    Takes T as a matrix of thruster magnitudes, alpha as the angle of the thrusters, theta as the angle of the AUV, and mass as the mass of the AUV
    """
    T_dir = np.array(
        [
            alpha + theta,
            2 * np.pi - alpha + theta,
            np.pi + alpha + theta,
            np.pi - alpha + theta,
        ]
    )
    hor_components = np.array([0, 0, 0, 0])
    ver_components = np.array([0, 0, 0, 0])
    for i in range(4):
        hor_components[i] = T[i] * np.cos(T_dir[i])
        ver_components[i] = T[i] * np.sin(T_dir[i])

    return np.ndarray([np.sum(hor_components), np.sum(ver_components)]) / mass


def calculate_auv2_angular_acceleration(T, alpha, L, l, inertia):
    """calculates angular acceleration of the AUV
    Takes T as a matrix of thruster magnitudes, alpha as the angle of the thrusters, L and l as the longitudinal and lateral distance of the thrusters to the middle, and inertia for moment of inertia
    """
    T_dir = np.ndarray([alpha, 2 * np.pi - alpha, np.pi + alpha, 180 - alpha])
    arm_vectors = np.ndarray([[l, -L, 0], [l, L, 0], [-l, L, 0], [-l, -L, 0]])
    total_torque = 0
    for i in range(4):
        torque = np.cross(
            arm_vectors[i], [T[i] * np.cos(T_dir[i]), T[i] * np.sin(T_dir[i]), 0]
        )
        total_torque = total_torque + torque
    return total_torque / inertia


def simulate_auv2_motion(
    T, alpha, L, l, mass=100, inertia=100, dt=0.1, t_final=10, x0=0, y0=0, theta0=0
):
    """simulates the AUVs motion
    Takes T as the matrix of thruster magnitudes, alpha as the angle of the thrusters, L and l as the longitudinal and lateral distance of the thrusters to the middle, mass for the mass of AUV, inertia for moment of inertia, dt as the time step of the simulation in seconds, t_final as the final time of the simulation, x0 and y0 as the initial position, and theta0 as the initial angle of the AUV
    """
    AUV_theta = theta0
    vel = np.zeros(2)
    pos = np.zeros(2)
    # np.ndarray([x0,y0])
    ang_vel = 0
    num = int(t_final / dt)
    t = np.zeros(num)
    x = np.zeros(num)
    y = np.zeros(num)
    t = np.zeros(num)
    theta = np.zeros(num)
    v = np.zeros(num)
    omega = np.zeros(num)
    a = np.zeros(num)
    accel = 0
    time = 0

    for i in range(num):
        t[i] = time
        x[i] = pos[0]
        y[i] = pos[1]
        theta[i] = AUV_theta
        v[i] = np.sqrt(v[0] ** 2 + v[1] ** 2)
        omega[i] = ang_vel
        a[i] = accel
        accel = calculate_auv2_acceleration(T, alpha, theta, mass)
        vel = np.add(vel, accel * dt)
        pos = np.add(pos, vel * dt)
        ang_accel = calculate_auv2_angular_acceleration(T, alpha, L, l, inertia)
        ang_vel = ang_vel + ang_accel * dt
        theta = theta + ang_vel
        time += dt
    return np.ndarray([t, x, y, theta, v, omega, a])


print(
    simulate_auv2_motion(
        np.array([1, 1, 1, -1]), np.pi / 30, 2, 2, inertia=0.01, theta0=30, mass=0.01
    )
)
