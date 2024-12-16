import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import CoolProp.CoolProp as CP

t_system = 80 + 273.15  # K
p_system = 101325  # Pa
g = 9.81  # m/s2

d_0 = 0.002  # диаметр сопла (m)

H_max = 1.7
h = np.linspace(0, H_max, num=100)


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

C_d0 = (1.1674 * rho_d) / 100      # Начальная концентрация этантиола в капле, кг/м^3

# расчет размера капли
delta_rho = rho_c - rho_d

gamma_value = np.sqrt((2 * sigma) / (g * delta_rho))

r_value = d_0 / (2 * gamma_value)
print(r_value)


def get_f_R(R):
    R_values = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    f_R_values = np.array([1.0, 0.74, 0.68, 0.64, 0.61, 0.60, 0.60, 0.61, 0.63, 0.64, 0.66])

    interp_func = interp1d(R_values, f_R_values, kind='cubic')

    return interp_func(R)


f_R_graph = get_f_R(r_value)
print(f_R_graph)

v_value = np.pi * r_value * f_R_graph

d_calc = gamma_value * ((6 * v_value) / np.pi) ** (1 / 3)
print(f'Расчетный радиус капли {d_calc:.4f}, м')

r_calc = d_calc / 2

# расчет скорости свободного осаждения капли

p_value = (rho_c ** 2 * sigma ** 3) / (g * delta_rho * mu_c ** 4)

t_value = (4 * delta_rho * g * (d_calc ** 2) * (p_value ** 0.15)) / (3 * gamma_value)

q_value = (22 * t_value) ** 0.42

Re_value = (q_value - 0.75) * (p_value ** 0.15)
print(Re_value)

w_0 = (Re_value * mu_c) / (rho_c * d_calc)
print(f'Скорость свободного осаждения капли: {w_0:.3f}, м/с')

m = 24  # Коэффициент распределения (безразмерный)
D = 5e-7   # Коэффициент диффузии этантиола в воде, м^2/с
Sc_value = mu_d / (rho_d * D)

Sh_value = 2 + 0.552 * (Re_value**0.5) * (Sc_value**(1/3))

betta_d = (Sh_value * D) / d_calc

# Вычисление концентрации C_d(h)
C_d_h = C_d0 * np.exp(- (3 * betta_d * h) / (r_calc * m * w_0))

conc1 = (1.017 * rho_d) / 100 # % масс
h_exp1 = 0.40
print(conc1)

conc2 = (0.982 * rho_d) / 100 # % масс
h_exp2 = (0.40 + 0.40)
print(conc2)

conc3 = (0.932 * rho_d) / 100 # % масс
h_exp3 = (0.4 + 0.75)
tau_exp3 = (0.4 + 0.75) / w_0
print(conc3)

conc4 = (0.891 * rho_d) / 100 # % масс
h_exp4 = (0.4 + 0.75 + 0.4)
print(conc4)

# Построение графика зависимости концентрации от высоты
plt.figure(figsize=(8, 6))
plt.plot(h, C_d_h, label='Расчет')
plt.scatter(h_exp1, conc1, label=f'H=0.40 м, {conc1:.2} (кг/м³)', color='red')
plt.scatter(h_exp2, conc2, label=f'H=0.80 м, {conc2:.2} (кг/м³)', color='red')
plt.scatter(h_exp3, conc3, label=f'H=1.15 м, {conc3:.2} (кг/м³)', color='red')
plt.scatter(h_exp4, conc4, label=f'H=1.55 м, {conc4:.2} (кг/м³)', color='red')
plt.xlabel('Высота h (м)')
plt.ylabel('Концентрация C_d(h) (кг/м³)')
plt.title('Изменение концентрации этантиола в капле при всплытии')
plt.grid(True)
plt.legend()
plt.show()
