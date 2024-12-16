import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def odes(t, y, params):
    # Распаковка переменных
    c0R, c1R, c2Na, c3RNa = y

    # Распаковка параметров
    beta = params['beta']
    V1 = params['V1']
    V2 = params['V2']
    S = params['S']
    phi = params['phi']
    k1 = params['k1']
    Kp = params['Kp']

    # Система уравнений
    dc0R_dt = -beta * V2 * S * (c0R - phi * c1R)
    dc1R_dt = beta * V2 * S * (c0R - phi * c1R) * (V1 / V2) - k1 * c1R * c2Na + (k1 / Kp) * c3RNa
    dc2Na_dt = -k1 * c1R * c2Na + (k1 / Kp) * c3RNa
    dc3RNa_dt = k1 * c1R * c2Na - (k1 / Kp) * c3RNa

    return [dc0R_dt, dc1R_dt, dc2Na_dt, dc3RNa_dt]


# Пример задания параметров
params = {
    'beta': 0.5,  # примерные значения
    'V1': 1.0,
    'V2': 2.0,
    'S': 0.1,
    'phi': 0.2,
    'k1': 0.3,
    'Kp': 1.0
}

# Начальные условия, задайте свои
y0 = [12.584, 0.0, 0.0, 0.0]

stop = 50
start = 0

# Интервал интегрирования по времени
t_span = (start, stop)

# Решение системы
sol = solve_ivp(lambda t, y: odes(t, y, params), t_span, y0, t_eval=np.linspace(start, stop, 100))

# Визуализация решения
plt.figure(figsize=(10, 6))
plt.plot(sol.t, sol.y[0], label='c0R')
plt.plot(sol.t, sol.y[1], label='c1R')
plt.plot(sol.t, sol.y[2], label='c2Na')
plt.plot(sol.t, sol.y[3], label='c3RNa')

plt.xlabel('Время, сек.')
plt.ylabel('Концентрация, масс доли')
plt.title('Решение системы ОДУ')
plt.legend()
plt.grid(True)
plt.show()
