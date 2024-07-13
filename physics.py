import numpy as np

G = 9.81
D_WATER = 1000


def calculate_bouyancy(V, density_fluid):
    return V * density_fluid * G


def will_it_float(V, mass):
    F_up = calculate_bouyancy(V, D_WATER)
    F_down = mass * G
    return F_up > F_down


def calculate_pressure(depth):
    return depth * G * D_WATER


def calculate_acceleration(F, m):
    return F / m


def calculate_angular_acceleration(tau, I):
    return tau / I


def calculate_torque_radian(F_magitude, F_direction, r):
    radian_direction = F_direction * (np.pi) / 180
    return r * F_magitude * np.sin(radian_direction)


def calculate_torque_degree(F_magitude, F_direction, r):
    return r * F_magitude * np.sin(F_direction)


def calculate_moment_of_inertia(m, r):
    m * (r**2)


def calculate_auv_acceleration(
    F_magnitude, F_angle, mass=100, volume=0.1, thruster_distance=0.5
):
    return F_magnitude / mass


def calculate_auv_angular_acceleration(
    F_magnitude, F_angle, intertia=1, thruster_distance=0.5
):
    assert -30 < F_angle < 30
    torque = calculate_torque_degree(F_magnitude, F_angle, thruster_distance)
    return torque / intertia


def calculate_auv2_acceleration(T_mag, alpha, theta, mass=100):
    Thruster_angle = np.array(
        [theta + alpha, theta - alpha, theta + np.pi + alpha, theta + np.pi - alpha]
    )

    # Initialize numpy arrays for components
    Thruster_hor = np.zeros((4))
    Thruster_ver = np.zeros((4))

    # Set the thruster components
    for i in range(4):
        Thruster_hor[i] = np.cos(Thruster_angle[i]) * T_mag[i]
        Thruster_ver[i] = np.sin(Thruster_angle[i]) * T_mag[i]

    # pythagorean theorem
    Thruster_force = [np.sum(Thruster_hor), np.sum(Thruster_ver)]

    return np.linalg.norm(Thruster_force) / mass


def calculate_auv2_angular_acceleration(Thruster_mag, alpha, L, l, inertia=100):
    # Thruster vectors
    Thruster_angles = np.array([alpha, -alpha, np.pi + alpha, np.pi - alpha])
    Thruster_vectors = np.zeros((4,2))
    for i in range (4):
        Thruster_vectors[i,0] = Thruster_mag[i]*np.cos(Thruster_angles[i])
        Thruster_vectors[i,1] = Thruster_mag[i]*np.sin(Thruster_angles[i])

    #Moment arm vectors
    moment_arms = np.array([[l,-L],
                            [l,L],
                            [-l,L],
                            [-l,-L]])

    #Find torques
    torque = np.zeros(4)
    for i in range (4):
        torque[i] = np.cross(moment_arms[i], Thruster_vectors[i])
    print(torque)
    
    #Find acceleration
    return np.sum(torque)/inertia

def simulate_auv2_motion():
    pass


def plot_auv2_motion():
    pass
