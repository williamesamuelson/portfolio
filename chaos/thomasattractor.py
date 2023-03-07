"""Study and calculations on Thomas' cyclically symmetric attractor"""

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

plt.rcParams["lines.linewidth"] = 1
plt.rcParams["axes.labelsize"] = 15


class Thomas:
    def __init__(self, b):
        self.b = b

    def __call__(self, t, x):
        y = np.zeros(3)
        y[0] = math.sin(x[1]) - self.b*x[0]
        y[1] = math.sin(x[2]) - self.b*x[1]
        y[2] = math.sin(x[0]) - self.b*x[2]
        return y


SIGMA = 10
RHO = 28
BETA = 8/3


def lorenz(t, x):
    """The Lorenz system"""
    y = np.zeros(3)
    y[0] = SIGMA*(x[1] - x[0])
    y[1] = x[0]*(RHO - x[2]) - x[1]
    y[2] = x[0]*x[1] - BETA*x[2]
    return y


# def thomas(t, x):
#     """The rhs of Thomas' cyclically symmetric attractor."""
#     y = np.zeros(3)
#     y[0] = math.sin(x[1]) - b*x[0]
#     y[1] = math.sin(x[2]) - b*x[1]
#     y[2] = math.sin(x[0]) - b*x[2]
#     return y


def rescale(dx_old, epsilon):
    """Rescale delta_x at each lyapunov step"""
    if np.linalg.norm(dx_old) == 0:
        return epsilon*np.array([1, 1, 1])/np.linalg.norm(np.array([1, 1, 1]))
    return epsilon*dx_old/np.linalg.norm(dx_old)


def lyapunov_step(T, x, y, thomas, h):
    """Do one lyapunov step"""
    evals = int(T/h)
    t = np.linspace(0, T, evals)
    sol_x = solve_ivp(thomas, (0, T), x, t_eval=t)
    sol_y = solve_ivp(thomas, (0, T), y, t_eval=t)
    # sol_x = solve_ivp(thomas, (0, T), x, atol=1e-15, rtol=1e-11)
    # sol_y = solve_ivp(thomas, (0, T), y, atol=1e-15, rtol=1e-11)
    delta = sol_x.y - sol_y.y
    xnew = sol_x.y[:, -1]
    return delta, xnew, t


def lyapunov(N, T, x, delta_0, epsilon, thomas, h, plot=False):
    """Calculate the largest Lyapunov exponent."""
    y = x + epsilon*delta_0/np.linalg.norm(delta_0)
    deltas = np.zeros(N)
    if plot:
        tvec = np.array([])
        deltaplot = np.array([])
    for i in range(N):
        deltavec, x, t = lyapunov_step(T, x, y, thomas, h)
        deltas[i] = np.linalg.norm(deltavec[:, -1])
        delta = rescale(deltavec[:, -1], epsilon)
        y = x - delta
        if plot:
            tvec = np.append(tvec, t+i*T)
            temp = np.linalg.norm(deltavec, axis=0)
            deltaplot = np.append(deltaplot, temp)
    if plot:
        plt.semilogy(tvec, deltaplot)
        plt.xlabel('Time [t]')
        plt.ylabel(r'$|\delta\mathbf{x}(t)|$')
        plt.grid(True)
        plt.show()
    lamda = 1/(N*T) * sum(np.log(deltas/epsilon))
    return lamda


def kaplan_yorke(lamda, b):
    return 2+lamda/abs(3*b+lamda)


def plot_lyap_kaplan(b_vec, N, T, x0, delta_0, epsilon, h):
    lyap_vec = np.zeros(len(b_vec))
    kaplan_vec = np.zeros(len(b_vec))
    thomas = Thomas(None)
    for i in range(len(b_vec)):
        thomas.b = b_vec[i]
        lamda = lyapunov(N, T, x0, delta_0, epsilon, thomas, h)
        kaplan = kaplan_yorke(lamda, thomas.b)
        lyap_vec[i] = lamda
        kaplan_vec[i] = kaplan
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(b_vec[::-1], lyap_vec[::-1])
    ax2.plot(b_vec[::-1], kaplan_vec[::-1])
    ax1.set_xlabel('b')
    ax1.set_ylabel(r'$\lambda_1$', rotation=0)
    ax2.set_xlabel('b')
    ax2.set_ylabel(r'$D_{KY}$', rotation=0)
    ax1.grid(True)
    ax2.grid(True)
    plt.show()



if __name__ == '__main__':
    N = 15
    T = 100
    x0 = np.array([1, 1.2, -0.3])
    delta_0 = np.array([1, 1, 1])
    EPSILON = 1e-14
    h = 0.05
    # lamb = lyapunov(N, T, x0, delta_0, EPSILON, Thomas(0.177), h, plot=True)
    b_vec = np.linspace(0.005, 0.21, 100)
    plot_lyap_kaplan(b_vec, N, T, x0, delta_0, EPSILON, h)

    # t_end = 5000
    # sol = solve_ivp(Thomas(0.177), (0, t_end), x0, rtol=1e-6, atol=1e-9)
    # colors = np.linspace(0, 1, len(sol.t))
    # ax = plt.figure().add_subplot(projection='3d')
    # # ax.scatter(sol.y[0, :], sol.y[1, :], sol.y[2, :], c=colors, s=1)
    # ax.plot(sol.y[0, :], sol.y[1, :], sol.y[2, :])
    # ax.scatter(sol.y[0, 0], sol.y[1, 0], sol.y[2, 0], s=10, c='black')
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_zlabel('z')
    # # #plt.savefig('sim11.png', bbox_inches='tight')
    # plt.show()
