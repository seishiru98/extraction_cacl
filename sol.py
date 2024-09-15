import matplotlib.pyplot as plt
import numpy as np

sol_in_water = np.array([0.0684, 0.0626, 0.0569, 0.0479, 0.0387, 0.0224, 0.0115])
sol_in_organic = np.array([3.1441, 2.666, 2.4574, 2.0465, 1.6259, 0.8301, 0.4275])

# Построение и сохранение графика
plt.figure(figsize=(10, 6))
plt.plot(sol_in_organic, sol_in_water, label=r'30')
#plt.scatter(tau_exp1, conc1, label=f'H=0.40 м', color='red')
plt.ylabel(r'В водной фазе, % масс.', fontsize=12)
plt.xlabel(r'В органической фазе, % масс.', fontsize=12)
plt.title('Растворимость этантиола', fontsize=16)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.85))
plt.grid(True)
plt.show()