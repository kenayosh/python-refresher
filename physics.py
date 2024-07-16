import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

G = 9.81
D_WATER = 1000
ONE_ATM = 101325

# TODO
# Doc strings
# Unit tests


def calculate_bouyancy(
    V: float,  # Type hints (hints for docs that says V should be int), python does not enforce this.
    density_fluid: int,
):
    """This calculates bouyancy
    Takes v as volume (m^3), and density_fluid as density (kg/m^3)
    Returns bouyancy
    """
    if type(V) is not isinstance(int) or type(V) is not isinstance(float):
        raise ValueError("V must be an integer")
    return V * density_fluid * G


def will_it_float(V, mass):
    F_up = calculate_bouyancy(V, D_WATER)
    F_down = mass * G
    return F_up > F_down


def calculate_pressure(depth):
    return depth * G * D_WATER + ONE_ATM


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
    Fx = F_magnitude * np.cos(F_angle)
    Fy = F_magnitude * np.sin(F_angle)
    return (Fx / mass, Fy / mass)


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
    Thruster_force = np.array([np.sum(Thruster_hor), np.sum(Thruster_ver)])

    return Thruster_force / mass


def calculate_auv2_angular_acceleration(Thruster_mag, alpha, L, l, inertia=100):
    # Thruster vectors
    Thruster_angles = np.array([alpha, -alpha, np.pi + alpha, np.pi - alpha])
    Thruster_vectors = np.zeros((4, 3))

    for i in range(4):
        Thruster_vectors[i, 0] = Thruster_mag[i] * np.cos(Thruster_angles[i])
        Thruster_vectors[i, 1] = Thruster_mag[i] * np.sin(Thruster_angles[i])
        Thruster_vectors[i, 2] = 0

    # Moment arm vectors
    moment_arms = np.array([[l, -L, 0], [l, L, 0], [-l, L, 0], [-l, -L, 0]])

    # Find torques
    torque = np.zeros(4)
    for i in range(4):
        torque[i] = (np.cross(moment_arms[i], Thruster_vectors[i]))[2]
    # Find acceleration
    return np.sum(torque) / inertia


def simulate_auv2_motion(
    Thruster_mag,
    alpha,
    L,
    l,
    mass=100,
    inertia=100,
    dt=0.1,
    t_final=1,
    x0=0,
    y0=0,
    theta0=0,
):
    # Calculate amount of steps + initialize return arrays
    steps = int(t_final / dt)
    np_t = np.zeros(steps)
    np_x = np.zeros(steps)
    np_y = np.zeros(steps)
    np_theta = np.zeros(steps)
    np_v = np.zeros(steps)
    np_omega = np.zeros(steps)
    np_a = np.zeros(steps)

    x_old = x0
    y_old = y0
    theta_old = theta0
    t_old = 0
    vx_old = 0
    vy_old = 0
    omega_old = 0

    for i in range(steps):
        #####CALCULATE NEW VALUES
        # Calculate accelerations
        a_x, a_y = (
            calculate_auv2_acceleration(Thruster_mag, alpha, theta_old, mass)[0],
            calculate_auv2_acceleration(Thruster_mag, alpha, theta_old, mass)[1],
        )
        a_angular = calculate_auv2_angular_acceleration(
            Thruster_mag, alpha, L, l, inertia
        )

        # Calculate velocities
        vx_new = vx_old + a_x * dt
        vy_new = vy_old + a_y * dt
        omega_new = omega_old + a_angular * dt
        # Calculate positions
        x_new = x_old + vx_new * dt
        y_new = y_old + vy_new * dt
        theta_new = theta_old + omega_new * dt
        # update t
        t_new = t_old + dt
        #####ADD VALUES TO THE RETURN FUNCTION
        np_t[i] = t_new
        np_x[i] = x_new
        np_y[i] = y_new
        np_theta[i] = theta_new
        np_v[i] = np.sqrt(vx_new**2 + vy_new**2)
        np_omega[i] = omega_new
        np_a[i] = np.sqrt(a_x**2 + a_y**2)
        #####SET THE OLD VALUES TO NEW VALUES
        x_old = x_new
        y_old = y_new
        theta_old = theta_new
        t_old = t_new
        vx_old = vx_new
        vy_old = vy_new
        omega_old = omega_new

    return (np_t, np_x, np_y, np_theta, np_v, np_omega, np_a, L, l)


import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def plot_auv2_motion(
    Thruster_mag,
    alpha,
    L,
    l,
    mass=100,
    inertia=100,
    dt=0.1,
    t_final=1,
    x0=0,
    y0=0,
    theta0=0,
):

    fig = plt.figure(figsize=(10, 10))
    plt.plot()
    plt.gca().set_aspect("equal", adjustable="box")

    motion = simulate_auv2_motion(
        Thruster_mag, alpha, L, l, mass, inertia, dt, t_final, x0, y0, theta0
    )

    for i in range(int(t_final / dt)):
        plt.gca().add_patch(
            Rectangle(
                (motion[1][i], motion[2][i]),
                2 * l,
                2 * L,
                angle=motion[3][i],
                edgecolor="red",
                facecolor="none",
                lw=4,
            )
        )

    plt.savefig("plot.png")


print(
    simulate_auv2_motion(np.array([1000, 1000, 1000, 1000]), 2, 3, 3, mass=1, inertia=1)
)
plot_auv2_motion(np.array([1000, 1000, 0, 0]), 2, 3, 3, mass=1, inertia=1, t_final=1)
