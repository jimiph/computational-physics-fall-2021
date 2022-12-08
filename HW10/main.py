import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss
from classes import van_der_waals_monoatomic_MD
import matplotlib.animation as animation
import scipy.constants as CONSTANTS



N = 100
max_initial_velocity_component = 1.5
h = 10**-3
length = 50
d = 2
simulation_time = 200
sampling_step = 100
sigma = 3.405 * 10 ** -10
mass = 6.63 * 10 ** -26
epsilon = 1.653 * 10 ** -21
samples_num = int(simulation_time / (h*sampling_step))
"""
RUN for 100 atoms and making trajectory. It is not necessary to consider cut_off for this number
of atoms; because run time is just about 1 min.
"""
# argon_gas = van_der_waals_monoatomic_MD(atoms_num=N, container_length=length, dimension=2,
#                                         max_initial_velocity_component=max_initial_velocity_component,
#                                         simulation_time=simulation_time, sampling_step=sampling_step,
#                                         integrating_step=h, cut_off=False)
# data = argon_gas.run_system_and_make_files()
# data.plot_energies()
# data.plot_left_side_atoms_number_fraction()
# data.plot_auto_correlation_of_velocites()
# data.plot_pressure()
# data.plot_temperature()
# data.make_animation()
# data.print_data(return_in_original_units=False)


"""
RUN for some varying initial max_initial_velocity_component and making plot P_T for 
verifecation of the Van der Wals's equation of state.
"""
# max_initial_velocity_component_s = np.linspace(1, 4, 30, endpoint=False)
# P_s = np.zeros(len(max_initial_velocity_component_s))
# P_err_s = np.zeros(len(max_initial_velocity_component_s))
# T_s = np.zeros(len(max_initial_velocity_component_s))
# T_err_s = np.zeros(len(max_initial_velocity_component_s))
# # E_s = np.zeros(len(max_initial_velocity_component_s))
# argon_gas = van_der_waals_monoatomic_MD(atoms_num=N, container_length=30, dimension=2,
#                                         max_initial_velocity_component=max_initial_velocity_component,
#                                         simulation_time=simulation_time, sampling_step=sampling_step,
#                                         integrating_step=h, cut_off=False)

# for i in range(len(max_initial_velocity_component_s)):
#     argon_gas.max_v0_comp = max_initial_velocity_component_s[i]
#     data = argon_gas.run_system_and_make_files()
#     P_s[i], P_err_s[i] = data.mean_and_error_of_pressure()
#     T_s[i], T_err_s[i] = data.mean_and_error_of_temperature()


# np.save('P_s_for_equation_of_state_plot.npy', P_s)
# np.save('P_err_s_for_equation_of_state_plot.npy', P_err_s)
# np.save('T_s_for_equation_of_state_plot.npy', T_s)
# np.save('T_err_s_for_equation_of_state_plot.npy', T_err_s)
# P_s = np.load('P_s_for_equation_of_state_plot.npy')
# P_err_s = np.load('P_err_s_for_equation_of_state_plot.npy')
# T_s = np.load('T_s_for_equation_of_state_plot.npy')
# T_err_s = np.load('T_err_s_for_equation_of_state_plot.npy')

# slope, intercept, r, pp, se = ss.linregress(T_s, P_s)
# plt.errorbar(T_s, P_s, xerr=T_err_s, yerr=P_err_s, label=f'Data', linestyle=' ', marker='o')
# plt.plot(T_s, slope*T_s+intercept, label='Fit with Eq. '+r'$P = $'+f'{round(slope, 3)}'+r'$T$'+f'+{round(intercept, 3)}')
# plt.xlabel(r'T $(\times \epsilon / k_B)$')
# plt.ylabel(r'P $(\times \epsilon / \sigma^' + f'{d}' + ')$')
# plt.title('Pressure per Temperature')
# plt.legend()
# plt.savefig('figures/P_T_plot.jpg')
# plt.show()


"""
RUN for phase transition.
"""
length = 30
max_initial_velocity_component = 2

"""
I used https://github.com/sina-moammar/2020-Fall-Computational-Physics/blob/main/HW9/main.py
to determine proper scales.
"""
# scales = 0.9 * np.ones(5)
# scales = np.concatenate((scales, 0.98 * np.ones(10)))
# scales = np.concatenate((scales, 0.95 * np.ones(10)))
# scales = np.concatenate((scales, 0.9 * np.ones(5)))
# scales = np.concatenate((scales, 0.85 * np.ones(5)))
# scales = np.concatenate((scales, 0.8 * np.ones(5)))
# scales = np.concatenate((scales, 0.75 * np.ones(5)))

# gas1 = van_der_waals_monoatomic_MD(atoms_num=N, container_length=length, dimension=2,
#                                         max_initial_velocity_component=max_initial_velocity_component,
#                                         simulation_time=simulation_time, sampling_step=sampling_step,
#                                         integrating_step=h, cut_off=False)

# data = gas1.run_system_and_make_files()

# gas2 =  van_der_waals_monoatomic_MD(atoms_num=N, container_length=length, dimension=2,
#                                         max_initial_velocity_component=max_initial_velocity_component,
#                                         simulation_time=simulation_time, sampling_step=sampling_step,
#                                         integrating_step=h, cut_off=False)

# gas2._custom_initial_positions_velocities(gas1.positions, gas1.velocities)
# gas1 = gas2
# data = gas1.run_system_and_make_files()

# E_s = np.zeros(len(scales))
# E_err_s = np.zeros(len(scales))
# T_s = np.zeros(len(scales))
# T_err_s = np.zeros(len(scales))


# for index in range(len(scales)):
#     gas1.velocities = gas1.velocities * scales[index]
#     data = gas1.run_system_and_make_files()
#     T_s[index], T_err_s[index] = data.mean_and_error_of_temperature()
#     E_s[index], E_err_s[index] = data.mean_and_error_of_energy()

# data.make_animation(phase_trans=True)

# np.save('E_s_for_phase_transition_plot.npy', E_s)
# np.save('E_err_s_for_phase_transition_plot.npy', E_err_s)
# np.save('T_s_for_phase_transition_plot.npy', T_s)
# np.save('T_err_s_for_phase_transition_plot.npy', T_err_s)
# E_s = np.load('E_s_for_phase_transition_plot.npy')
# E_err_s = np.load('E_err_s_for_phase_transition_plot.npy')
# T_s = np.load('T_s_for_phase_transition_plot.npy')
# T_err_s = np.load('T_err_s_for_phase_transition_plot.npy')
# plt.scatter(T_s, E_s)
# plt.xlabel(r'T $(\times \epsilon / k_B)$')
# plt.ylabel(r'Energy $(\times \epsilon)$')
# plt.title('Energy per Temperature in Phase Transition')
# plt.savefig(f'figures/MD_{length}_phase_transition_E_T.jpg')
# plt.show()
