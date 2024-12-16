import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

from scipy.integrate import odeint
from scipy.interpolate import interp1d


t_system = 20 + 273.15  # K
p_system = 101325  # Pa
g = 9.81  # m/s2


def get_f_sigma(S):
    sigma_values = np.array([0.053, 0.051, 0.049, 0.047, 0.045, 0.043, 0.041, 0.039])
    f_sigma_values = np.array([283.15, 293.15, 303.15, 313.15, 323.15, 333.15, 343.15, 353.15])

    interp_func_sigma = interp1d(f_sigma_values, sigma_values, kind='cubic')

    return interp_func_sigma(S)


sigma = get_f_sigma(t_system)  # N/m
print(f'Межфазное натяжение Вода/Октсан: {sigma} Н/м')

rho_d = CP.PropsSI('D', 'T', t_system, 'P', p_system, 'n-Octane')  # kg/m³
rho_c = CP.PropsSI('D', 'T', t_system, 'P', p_system, 'Water')  # kg/m³
print(f'Плотность октана: {rho_d} кг/м3')
print(f'Плотность воды: {rho_c} кг/м3')

mu_d = CP.PropsSI('V', 'T', t_system, 'P', p_system, 'n-Octane')  # Pa*s
mu_c = CP.PropsSI('V', 'T', t_system, 'P', p_system, 'Water')  # Pa*s
print(f'Вязкость октана: {mu_d} Па с')
print(f'Вязкость воды: {mu_c} Па с')

# расчет размера капли
delta_rho = rho_c - rho_d

gamma_value = np.sqrt((2 * sigma) / (g * delta_rho))

d_0 = 0.002  # диаметр сопла (m)

r_value = d_0 / (2 * gamma_value)


def get_f_R(R):
    R_values = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    f_R_values = np.array([1.0, 0.74, 0.68, 0.64, 0.61, 0.60, 0.60, 0.61, 0.63, 0.64, 0.66])

    interp_func = interp1d(R_values, f_R_values, kind='cubic')

    return interp_func(R)


f_R_graph = get_f_R(r_value)

v_value = np.pi * r_value * f_R_graph

d_calc = gamma_value * ((6 * v_value) / np.pi) ** (1 / 3)
print(f'Расчетный диаметр капли {d_calc:.4f}, м')

# расчет скорости свободного осаждения капли

p_value = (rho_c ** 2 * sigma ** 3) / (g * delta_rho * mu_c ** 4)

t_value = (4 * delta_rho * g * (d_calc ** 2) * (p_value ** 0.15)) / (3 * gamma_value)

q_value = (22 * t_value) ** 0.42

Re_value = (q_value - 0.75) * (p_value ** 0.15)

w_0 = (Re_value * mu_c) / (rho_c * d_calc)

w_0 = 0.15
print(f'Скорость капли: {w_0:.3f}, м/с')

D_d = 1.4424e-8

Re_d = rho_d * w_0 * d_calc / mu_d

Sc_d = mu_d / (rho_d * D_d)

Sh_d = 2 + 0.552 * (Re_d**0.5) * (Sc_d**(1/3))

betta_d = (Sh_d * D_d) / d_calc

D_c = 1.7943e-8

Re_c = rho_c * w_0 * d_calc / mu_c

Sc_c = mu_c / (rho_c * D_c)

Sh_c = 2 + 0.552 * (Re_d**0.5) * (Sc_c**(1/3))

betta_c = (Sh_c * D_c) / d_calc

m = 24

K = (1/betta_d + m/betta_c)**-1


V_drop = np.pi * d_calc**3 / 6
F = np.pi * d_calc**2

Q_o = 0.0000001
d_column = 0.051
S = np.pi * d_column**2 / 4

e_o = Q_o / S * w_0
print(e_o)
e_w = 1 - e_o
print(e_w)

a = 6 * e_o / d_calc
print(a)

# Определяем функцию для системы уравнений
def dCdt(y, t):
    c1, c2 = y

    dc1_dt = - (K * a / e_o)  * (c1 - m * c2)
    dc2_dt = (K * a / e_w) * (c1 - m * c2)

    return [dc1_dt, dc2_dt]

# Интервал времени
t_span = np.linspace(0, 10, 20)

# Начальные условия
c0 = [12.54, 0]

# Решение системы уравнений
solution = odeint(dCdt, c0, t_span)
print(solution)

# Построение графиков
plt.plot(t_span, solution[:, 0], label="c1 (t)")
#plt.plot(t_span, solution[:, 1], label="c2 (t)")
plt.xlabel("Time")
plt.ylabel("Concentrations")
plt.legend()
plt.grid()
plt.show()
