import numpy as np
import matplotlib.pyplot as plt
import pickle

plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['lines.linewidth'] = 2

t = 1
hbar = 1
a = 1
beta = 10
gamma = 0.1


def E(kx, ky, kz):
    return -2*t*(np.cos(a*kx) + np.cos(a*ky) + np.cos(a*kz))


def f(E):
    """
    Fermi-Dirac

    """
    return 1/(np.exp(E*beta) + 1)


def pi(q, omega, mesh_size):
    """
    Insert q_x, omega (scalars) and a mesh size (e.g 50). Returns the value
    of the "polarization operator".

    """
    # k-values in Brioullin zone
    k = np.linspace(0, 2*np.pi/a, mesh_size)
    # Create 3d grid (momentum grid)
    kx, ky, kz = np.meshgrid(k, k, k)
    Nk = mesh_size**3

    # get a 3d matrix of energies (E(k)) in each point in the grid
    E_k = E(kx, ky, kz)

    # get a 3d matrix of occupancy probabilities (f(E(k))) at each point in the grid
    f_k = f(E_k)

    # find index in k-vector with the closest value to q
    q_index = np.argmin(np.abs(k-q))

    # shift array to find a 3d matrix of E(k-q)
    E_kq = np.roll(E_k, q_index, axis=1)

    # find 3d matrix of f(E(k-q))
    f_kq = f(E_kq)

    # 3d-matrix of the total expression
    total = (f_kq - f_k)/(E_kq - E_k + hbar*omega + 1j*hbar*gamma)
    return 1/Nk * np.sum(total)  # sum over all k in the 3d grid


def av_kin_energy(mesh_size):
    k = np.linspace(-np.pi/a, np.pi/a, mesh_size)
    kx, ky, kz = np.meshgrid(k, k, k)
    E_k = E(kx, ky, kz)
    f_k = f(E_k)
    return -np.sum(E_k*f_k)/np.sum(f_k)


def save_object(obj):
    try:
        with open("mesh.pickle", "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as ex:
        print("Error during pickling object (Possibly unsupported):", ex)


def load_object(filename):
    try:
        with open(filename, "rb") as f:
            return pickle.load(f)
    except Exception as ex:
        print("Error during unpickling object (Possibly unsupported):", ex)


def calc_mesh(mesh_size, q, omegas):
    pi_mesh = np.array([[pi(qx, omega, mesh_size) for qx in q] for omega in omegas])
    save_object(pi_mesh)


if __name__ == '__main__':
    mesh_size = 50
    q = np.linspace(0, 2*np.pi/a, mesh_size)  # q-values in color plot
    omegas = np.linspace(0, 5, mesh_size)  # omegas in color plot
    # calc_mesh(mesh_size, q, omegas)
    pi_mesh = load_object('mesh.pickle')
    V0s = np.array([1, 10, 25])
    wp = np.zeros(V0s.shape)
    for index, V0 in enumerate(V0s):
        chi = np.array([list(V0/q_i**2 * np.array(pi_mesh[:, i])) for i, q_i in enumerate(q[1::], start=1)]).T
        plot_vals = (1/(1-chi)).imag
        wp[index] = omegas[np.argmin(plot_vals[:, 0])]

    fit = np.polyfit(np.sqrt(V0s), wp, 1)
    print(fit[0])

    # Q, Omega = np.meshgrid(q, omegas)  # coordinate system for color plot
    # plt.contourf(Q, Omega, plot_vals) #colorplot
    # plt.xlabel(r'$q_x$')
    # plt.ylabel(r'$\omega$')
    # plt.title(r'$\Im \{ 1/(1-\chi(q_x,\omega))\}$')
    # plt.colorbar()

    # plt.contourf(Q, Omega, plot_vals) #colorplot
    # plt.xlabel(r'$q_x$')
    # plt.ylabel(r'$\omega$')
    # plt.title(r'$\Im \{ 1/(1-\chi(q_x,\omega))\}$')
    # plt.colorbar()
    # plt.figure(dpi=300)
    # plt.contourf(Q, Omega, -pi_mesh.imag) #colorplot
    # slope = av_kin_energy(mesh_size)
    # plt.plot(q[0:int(mesh_size/3)], slope*q[0:int(mesh_size/3)], color='white')
    # plt.xlabel(r'$q_x$')
    # plt.ylabel(r'$\omega$')
    # plt.title(r'$-\Im\{\Pi(q_x,\omega)\}$')
    # plt.colorbar()
