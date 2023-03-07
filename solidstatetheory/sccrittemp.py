import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['lines.linewidth'] = 2


V0_D = 0.13
hbar_omega = 1


def create_integrand(log_delta):
    def integrand(E):
        return 1/np.sqrt((10**log_delta)**2 + E**2)
    return integrand


def create_integrand2(log_delta, kbT):
    def integrand(E):
        tanh_arg = np.sqrt((10**log_delta)**2 + E**2)/(2*kbT)
        return np.tanh(tanh_arg)/np.sqrt((10**log_delta)**2 + E**2)
    return integrand


def find_deltaT0(log_delta_plus, log_delta_minus, tol, max_steps):
    log_delta = (log_delta_plus + log_delta_minus)/2
    for i in range(max_steps):
        fun = create_integrand(log_delta)
        integral = V0_D/4*integrate.quad(fun, -hbar_omega, hbar_omega)[0]
        error = integral-1
        if np.abs(error) < tol:
            return 10**log_delta

        if error > 0:
            log_delta_minus = log_delta
        else:
            log_delta_plus = log_delta
        log_delta = (log_delta_plus + log_delta_minus)/2
    else:
        raise Exception("Max steps reached")


def calc_integral(log_delta, kbT):
    fun = create_integrand2(log_delta, kbT)
    integral = V0_D/4*integrate.quad(fun, -hbar_omega, hbar_omega)[0]
    return integral


if __name__ == '__main__':
    log_delta_plus = 0
    log_delta_minus = -15
    tol = 1e-5
    max_steps = 100
    deltaT0 = find_deltaT0(log_delta_plus, log_delta_minus, tol, max_steps)
    vec_calc_integral = np.vectorize(calc_integral)
    kbT_vals = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.56])
    kbTs = kbT_vals*deltaT0
    log_deltas = np.linspace(-10, -5, 30)
    result = np.zeros((len(log_deltas), len(kbTs)))
    for i, kbT in enumerate(kbTs):
        result[:, i] = vec_calc_integral(log_deltas, kbT)

    plt.figure(dpi=400)
    for i in range(len(kbTs)):
        plt.plot(log_deltas, result[:, i])

    leg = [r'$k_bT = $' + str(kbT) + r'$\Delta(T=0)$' for kbT in kbT_vals]
    plt.legend(leg)
    plt.title(r'$V_0D(\mu) = $' + str(V0_D))
    plt.xlabel(r'$\log_{10}\Delta$')
    plt.show()
