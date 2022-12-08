import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.constants as CONSTANTS


class van_der_waals_monoatomic_MD:
    def __init__(self, atoms_num=100, container_length=30, dimension=2,
                 max_initial_velocity_component=1,
                 sigma=3.405 * 10 ** -10, mass=6.63 * 10 ** -26, epsilon=1.653 * 10 ** -21,
                 simulation_time=200, sampling_step=100, integrating_step=10**-3, cut_off=False):
        """
        All attributes are supposed in reduced units.
        In positons array and velocities array:
            positions[0][:] : x and positions[1][:] : y
            velocities[0][:] : v_x and velocities[1][:] : v_y
        In  accelerations array each value shows the acceleleration of an atom 
            influenced by interation with other atoms. 

        """
        self.N = atoms_num
        self.L = container_length
        self.L_half = container_length / 2
        self.d = dimension
        self.V = self.L ** self.d
        self.cut_off_L = 2*np.power(2, 1/6)
        self.max_v0_comp = max_initial_velocity_component
        self.sigma = sigma
        self.mass = mass
        self.epsilon = epsilon
        self.simulation_time = simulation_time
        self.time = 0
        self.sampling_step = sampling_step
        self.h = integrating_step
        self.h_half = integrating_step / 2
        self.samples_num = int(simulation_time / (integrating_step*sampling_step))
        self.cut_off = cut_off
        self.positions = np.zeros((self.d, self.N))
        self._initiate_positions_to_left_side_crystal()
        self._apply_periodic_boundary()
        self.velocities = np.zeros((self.d, self.N))
        self._label_initial_velocities()
        self._transform_velocities_to_cm_frame()
        self._positions_differences = np.zeros((self.d, self.N, self.N))
        self._update_positions_differences()
        self._distances_power_2 = np.zeros((self.N, self.N))
        self._distances_power_6 = np.zeros((self.N, self.N))
        self._distances_power_12 = np.zeros((self.N, self.N))
        self._update_distances_for_all_powers()
        self.accelerations = np.zeros((self.d, self.N))
        self._update_accelerations()
        self.total_energy = None
        self.temperature = None
        self.potential_energy = None
        self.kinetic_energy = None
        self.pressure = None
        self.num_of_left_half_atoms = None

    def _initiate_positions_to_left_side_crystal(self):
        """
        initial left side lattice:
            *   *   *   * ... *     |
            .   .   .   .     .     |
            .   .   .   .     .     |
            .   .   .   .     .     |
            *   *   *   * ... *     |
        
        """
        row_s_num = int(np.sqrt(self.N)+1)
        lattice_constant = self.L_half / row_s_num
        for i in range(self.N):
            self.positions[0][i] = lattice_constant * (1+(i%(row_s_num-1)))
            self.positions[1][i] = lattice_constant * (1+(i//(row_s_num-1)))

    def _apply_periodic_boundary(self):
        """
        By this function I exert the peridic boundary as follows:
            if x (or y) > L then x (or y) = x (or y) - L
            if x (or y) < 0 then x (or y) = x (or y) + L
        """
        images_indices = self.positions > self.L
        self.positions[images_indices] -= self.L
        images_indices = self.positions < 0
        self.positions[images_indices] += self.L

    def _label_initial_velocities(self):
        """
        impose to each atom its initial velocity.
        """
        self.velocities[0] = self.max_v0_comp * np.random.uniform(-1, 1, self.N)
        self.velocities[1] = self.max_v0_comp * np.random.uniform(-1, 1, self.N)

    def _custom_initial_positions_velocities(self, positions, velocities):
        self.positions[:][:] = positions
        self.velocities[:][:] = velocities
        self._update_positions_and_distances()
        self._update_accelerations()

    def _transform_velocities_to_cm_frame(self):
        """
        Transform valocities of atoms into center of mass frame.
        """
        cm_velocity = np.mean(self.velocities, axis=1)
        for axis in range(self.d):
            self.velocities[axis] -= cm_velocity[axis]

    def _update_positions_differences(self):
        """
        Update positions differences of atoms from positons of atoms; also apply periodic
        condition on the positions differences.
        """
        for axis in range(self.d):
            self._positions_differences[axis] = np.tile(self.positions[axis], (self.N, 1))\
            - np.transpose(np.tile(self.positions[axis], (self.N, 1)))
        indices_1 = self._positions_differences > self.L_half
        self._positions_differences[indices_1] -= self.L
        indices_2 = self._positions_differences < -1*self.L_half
        self._positions_differences[indices_2] += self.L

    def _update_distances_for_all_powers(self):
        """
        Update distances for r^-2, r^-6 and r^-12, using positions differences of atoms.
        """
        self._distances_power_2 = 1 / np.sum(np.square(self._positions_differences), axis=0)
        np.fill_diagonal(self._distances_power_2, 0)
        self._distances_power_6 = self._distances_power_2 * self._distances_power_2 * self._distances_power_2
        self._distances_power_12 = self._distances_power_6 * self._distances_power_6
        
    def _update_positions_and_distances(self):
        """
        Update all positions differences and distances all at once.
        """
        self._apply_periodic_boundary()
        self._update_positions_differences()
        self._update_distances_for_all_powers()

    def _find_neighbours_in_cut_off_cell(self):
        """
        Return indices matrix of neighbours of each atom.
        """
        neighbours_indices_ = np.array(np.abs(self._positions_differences) < self.cut_off_L, dtype=int)
        for axis in range(self.d):
            np.fill_diagonal(neighbours_indices_[axis], 0)
        neighbours_indices = np.zeros((self.d, self.N, self.N))
        for axis in range(self.d):
            neighbours_indices[axis] = neighbours_indices_[0] * neighbours_indices_[1]
        return neighbours_indices

    def _update_accelerations(self):
        """
        Update accelerations of atoms.
        """
        if self.cut_off:
            neighbours_indices = self._find_neighbours_in_cut_off_cell()
            self.accelerations =  24 * np.sum(
                neighbours_indices * (-2 * self._distances_power_12 + self._distances_power_6) * self._distances_power_2 *
                self._positions_differences, axis=2)
        else:
            self.accelerations =  24 * np.sum(
                 (-2 * self._distances_power_12 + self._distances_power_6) * self._distances_power_2 *
                self._positions_differences, axis=2)
    
    def _update_system_for_one_integration_step(self):
        """
        Update the trajectory of the system for one integration step.
        """
        changes = self.accelerations * self.h_half
        self.positions += self.velocities * self.h + self.h * changes
        self.velocities += changes
        self._transform_velocities_to_cm_frame()
        self._update_positions_and_distances()
        self._update_accelerations()
        self.velocities += self.h_half * self.accelerations
        self._transform_velocities_to_cm_frame()
        self.potential_energy = None
        self.kinetic_energy = None
        self.temperature = None
        self.pressure = None
        self.num_of_left_half_atoms = None

    def run_system_and_make_files(self):
        """
        Run the system for the whole interval time of evolution determined by the given
        simulation time.
        files: A dictionary including names of files and their corresponding lists.
        info: A dictionary including all important attributes of the class van_der_waals_monoatomic_MD;
        Also this dictionary is given to the class analyser as an argument by returning
        analyser(info).
        """
        x_s, y_s, vx_s, vy_s = np.array([np.zeros((self.samples_num, self.N)) for _ in range(4)])
        T_s, K_s, U_s, E_s, P_s, N_l_s = np.array([np.zeros(self.samples_num) for _ in range(6)])

        for i in range(self.samples_num):
            print(f'[Info]: Making sample {i} of {self.samples_num}')
            for _ in range(self.sampling_step):
                self._update_system_for_one_integration_step()
            x_s[i] = self.positions[0]
            y_s[i] = self.positions[1]
            vx_s[i] = self.velocities[0]
            vy_s[i] = self.velocities[1]
            T_s[i] = self.calculate_temperature()
            K_s[i] = self.caculate_kinetic_energy()
            U_s[i] = self.calculate_potential_energy()
            E_s[i] = self.calculate_energy()
            P_s[i] = self.calculate_pressure()
            N_l_s[i] = self.calculate_num_of_left_half_atoms()
        self.total_energy = round(np.mean(E_s), 2)
        self.file_name_base = f'MD_simulation_N={self.N}_dimension={self.d}_L={self.L}E={self.total_energy}'
        files = {
            f'{self.file_name_base}_(x_s).npy' : x_s,
            f'{self.file_name_base}_(y_s).npy' : y_s,
            f'{self.file_name_base}_(vx_s).npy' : vx_s,
            f'{self.file_name_base}_(vy_s).npy' : vy_s,
            f'{self.file_name_base}_(T_s).npy' : T_s,
            f'{self.file_name_base}_(K_s).npy' : K_s,
            f'{self.file_name_base}_(U_s).npy' : U_s,
            f'{self.file_name_base}_(E_s).npy' : E_s,
            f'{self.file_name_base}_(P_s).npy' : P_s,
            f'{self.file_name_base}_(N_l_s).npy' : N_l_s,
                }
        info =  {
            'sigma': self.sigma,
            'mass': self.mass,
            'epsilon': self.epsilon,
            'length': self.L,
            'number of atoms': self.N,
            'dimension': self.d,
            'max_v0_comp': self.max_v0_comp,
            'h': self.h,
            'sampling step': self.sampling_step,
            'number of samples' : self.samples_num,
            'simulation time': self.simulation_time,
            'base of file name': self.file_name_base,
            'names of files': np.array([key for key in files.keys()]),
                }
        for item in files.items():
            np.save(item[0], item[1])
        return analyser(info)
    
    def calculate_temperature(self):
        if self.temperature is None:
            if self.kinetic_energy is None:
                self._transform_velocities_to_cm_frame()
                self.temperature = np.sum(np.square(self.velocities)) / ((self.N-1)*self.d)
            else:
                self._transform_velocities_to_cm_frame()
                self.temperature = self.kinetic_energy / (.5*(self.N-1)*self.d)
        return self.temperature
    
    def caculate_kinetic_energy(self):
        if self.kinetic_energy is None:
            if self.temperature is None:
                self._transform_velocities_to_cm_frame()
                self.kinetic_energy = 0.5 * np.sum(np.square(self.velocities))
            else:
                self._transform_velocities_to_cm_frame()
                self.kinetic_energy = self.temperature * 0.5 * (self.N-1) * self.d
        return self.kinetic_energy

    def calculate_potential_energy(self):
        if self.potential_energy is None:
            self._update_distances_for_all_powers()
            if self.cut_off:
                neighbours_indices = self._find_neighbours_in_cut_off_cell()
                self.potential_energy = 0.5 + 2*np.sum(neighbours_indices*(self._distances_power_12 - self._distances_power_6))
            else:
                self.potential_energy = 2*np.sum(self._distances_power_12 - self._distances_power_6)
        return self.potential_energy
    
    def calculate_energy(self):
        return self.calculate_potential_energy() + self.caculate_kinetic_energy()
    
    def calculate_pressure(self):
        if self.pressure is None:
            self.pressure = (self.N*self.calculate_temperature() + (12/self.d) * 
            np.sum(-2*self._distances_power_12+self._distances_power_6)) / self.V
        return self.pressure
    
    def calculate_num_of_left_half_atoms(self):
        if self.num_of_left_half_atoms is None:
            self.num_of_left_half_atoms = np.count_nonzero(self.positions[0] < self.L_half)
        return self.num_of_left_half_atoms
    

class analyser:
    """
    All of data analysing tasks including ploting and making animation of the trajectory of 
    the system are implemented by this class.
    info_of_md_simulation: A dictionary made by van_der_waals_monoatomic_MD.run_system_and_make_files()
    """
    def __init__(self, info_of_md_simulation):
        self.sigma = info_of_md_simulation['sigma']
        self.mass = info_of_md_simulation['mass']
        self.epsilon = info_of_md_simulation['epsilon']
        self.tau = np.sqrt(self.mass * self.sigma * self.sigma / self.epsilon)
        self.L = info_of_md_simulation['length']
        self.N = info_of_md_simulation['number of atoms']
        self.d = info_of_md_simulation['dimension']
        self.max_v0_comp = info_of_md_simulation['max_v0_comp'] 
        self.h = info_of_md_simulation['h']
        self.sampling_step = info_of_md_simulation['sampling step']
        self.samples_num = info_of_md_simulation['number of samples']
        self.simulation_time = info_of_md_simulation['simulation time']
        self.file_name_base = info_of_md_simulation['base of file name']
        self.names_of_files = info_of_md_simulation['names of files']
        self.reduced_times_for_plotting = self.h * self.sampling_step * np.arange(self.samples_num)
        self.equilibrium_index = None
        self._relaxation_index = None
        self._temperature = None
        self._temperature_error = None
        self._potential_energy = None
        self._potential_energy_error = None
        self._kinetic_energy = None
        self._kinetic_energy_error = None
        self._energy = None
        self._energy_error = None
        self._pressure = None
        self._pressure_error = None
        self._number_of_left_side_atoms = None
        self._number_of_left_side_atoms_error = None

    def calculate_autocorrelation_of_velocities(self):
        """
        Return autocorrelation of velocities as an array.
        """
        samples = np.arange(0, self.samples_num, dtype=int)
        c_v_s = np.zeros(len(samples))
        v_x, v_y = np.load(self.names_of_files[2]), np.load(self.names_of_files[3])
        v_s = np.sqrt(v_x*v_x + v_y*v_y)
        variance_of_v_s = np.var(v_s)
        thereshold = np.exp(-1)
        for sample_num in samples:
            indices = np.arange(0, self.samples_num-sample_num)
            c_v_s[sample_num] = (np.mean(v_s[indices]*v_s[indices+sample_num]) - \
            np.mean(v_s[indices])*np.mean(v_s[indices+sample_num])) / variance_of_v_s
            if c_v_s[sample_num] <= thereshold and self._relaxation_index is None:
                self._relaxation_index = sample_num
        if self._relaxation_index is None:
            self._relaxation_index = math.ceil(-sample_num[-1] / np.log(c_v_s[-1]))
        return c_v_s
        
    def make_animation(self, delay_in_ms=50, phase_trans=False):
        """
        Make animation and save.
        """
        def update_frame(sample_num, frame, ax):
            x_s = np.load(self.names_of_files[0])[sample_num] 
            y_s = np.load(self.names_of_files[1])[sample_num]
            frame.set_data(x_s, y_s)
            temperature = np.load(self.names_of_files[4])[sample_num]
            ax.set_title(r'T='+f'{round(temperature * self.epsilon / CONSTANTS.Boltzmann,2)} K')
            return frame,
        fig = plt.figure()
        ax = plt.gca()
        frame, = plt.plot([], [], linestyle='', marker='o')
        plt.xlim(0, self.L)
        plt.ylim(0, self.L)
        anima = animation.FuncAnimation(fig, update_frame, self.samples_num, fargs=(frame, ax), 
                                        interval = delay_in_ms, blit=True)
        if phase_trans==True:
            anima.save(f'animations/{self.file_name_base}_phase_transition.gif')
        else:
            anima.save(f'animations/{self.file_name_base}.gif')


    def _calculate_equilibrium_index(self):
        """
        I determine the equilibrium time using the fraction of left side atoms; when it
        reached a neighbourhood of 0.5 (with a percent error that can determine arbitrarily)
        we call that system is equlized.
        """
        fraction_of_left_side_particles_per_sample = np.load(self.names_of_files[9]) / self.N
        percent_error = 2
        equilibrium_index = 0
        for i in range(len(fraction_of_left_side_particles_per_sample)):
            fraction = fraction_of_left_side_particles_per_sample[i]
            if (fraction < 0.5*(1+percent_error/100)) and (fraction > 0.5*(1-percent_error/100)):
                equilibrium_index = i
                break
        if self.equilibrium_index is None:
            self.equilibrium_index = equilibrium_index
                
    def _caculate_mean_and_error_of_macro_over_time(self, macro_index):
        """
        macro : {'temperature', 'potential energy', 'kinetic energy',
                 'energy', pressure', 'number of left side particles'}

        Return the time average and error of the given macroscopic quantity.
        """
        
        macro = np.load(self.names_of_files[macro_index])
        mean = np.mean(macro[self.equilibrium_index:])
        error = np.std(macro[self.equilibrium_index:]) / np.sqrt(len(macro[self.equilibrium_index:]))
        return mean, error
        
    def mean_and_error_of_temperature(self, return_in_original_units=False):
        if self._temperature is None:
            self._temperature, self._temperature_error = self._caculate_mean_and_error_of_macro_over_time(4)
        if return_in_original_units:
            return (self.epsilon/CONSTANTS.Boltzmann * self._temperature, self.epsilon/CONSTANTS.Boltzmann * self._temperature_error)
        else:
            return self._temperature, self._temperature_error
    
    def mean_and_error_of_kinetic_energy(self, return_in_original_units=False):
        if self._kinetic_energy is None:
            self._kinetic_energy, self._kinetic_energy_error = self._caculate_mean_and_error_of_macro_over_time(5)
        if return_in_original_units:
            return (self.epsilon * self._kinetic_energy, self.epsilon * self._kinetic_energy_error)
        else:
            return self._kinetic_energy, self._kinetic_energy_error
        
    def mean_and_error_of_potential_energy(self, return_in_original_units=False):
        if self._potential_energy is None:
            self._potential_energy, self._potential_energy_error = self._caculate_mean_and_error_of_macro_over_time(6)
        if return_in_original_units:
            return (self.epsilon * self._potential_energy, self.epsilon * self._potential_energy_error)
        else:
            return self._potential_energy, self._potential_energy_error
    
    def mean_and_error_of_energy(self, return_in_original_units=False):
        if self._energy is None:
            self._energy, self._energy_error = self._caculate_mean_and_error_of_macro_over_time(7)
        if return_in_original_units:
            return (self.epsilon * self._energy, self.epsilon * self._energy_error)
        else:
            return self._energy, self._energy_error
    
    def mean_and_error_of_pressure(self, return_in_original_units=False):
        if self._pressure is None:
            self._pressure, self._pressure_error = self._caculate_mean_and_error_of_macro_over_time(8)
        if return_in_original_units:
            return (self.epsilon/(np.power(self.sigma, self.d)) * self._pressure, self.epsilon/(np.power(self.sigma, self.d)) * self._pressure_error)
        else:
            return self._pressure, self._pressure_error

    def mean_and_error_of_left_side_atoms_number(self):
        if self._number_of_left_side_atoms is None:
            self._number_of_left_side_atoms, self._number_of_left_side_atoms_error = \
                self._caculate_mean_and_error_of_macro_over_time(9)
        return self._number_of_left_side_atoms, self._number_of_left_side_atoms_error
    
    def plot_left_side_atoms_number_fraction(self):
        left_side_atoms_numbers = np.load(self.names_of_files[9])
        if self.equilibrium_index is None:
            self._calculate_equilibrium_index()
        t_c = self.reduced_times_for_plotting[self.equilibrium_index]
        t_c = round(t_c, 2)
        plt.axhline(0.5, color='red', linestyle='--')
        plt.plot(self.reduced_times_for_plotting, left_side_atoms_numbers / self.N)
        plt.plot(self.reduced_times_for_plotting[self.equilibrium_index], 0.5, 'go',
                 label=f't={t_c}')
        plt.title('Fraction of Left Side Atoms Number per Time')
        plt.xlabel(r'Time $(\times \tau)$')
        plt.ylabel(r'Fraction of left side atoms')
        plt.legend()
        plt.savefig(f'figures/{self.file_name_base}_left_side_atoms_fraction_plot.jpg')
        plt.show()

    def plot_temperature(self):
        temeperatures = np.load(self.names_of_files[4])  
        plt.plot(self.reduced_times_for_plotting, temeperatures)
        plt.title('Temperature per Time')
        plt.xlabel(r'Time $(\times \tau)$')
        plt.ylabel(r'T $(\times \epsilon / k_B)$')
        plt.savefig(f'figures/{self.file_name_base}_temperatures_plot.jpg')
        plt.show()

    def plot_pressure(self):
        pressures = np.load(self.names_of_files[8])
        plt.plot(self.reduced_times_for_plotting, pressures)
        plt.title('Pressure per Time')
        plt.xlabel(r'Time $(\times \tau)$')
        plt.ylabel(r'P $(\times \epsilon / \sigma^' + f'{self.d}' + ')$')
        plt.savefig(f'figures/{self.file_name_base}_pressures_plot.jpg')
        plt.show()

    def plot_energies(self):
        kinetic_energies = np.load(self.names_of_files[5])
        potential_energies = np.load(self.names_of_files[6])
        total_energies = np.load(self.names_of_files[7])
        plt.plot(self.reduced_times_for_plotting, kinetic_energies)
        plt.plot(self.reduced_times_for_plotting, potential_energies)
        plt.plot(self.reduced_times_for_plotting, total_energies)
        plt.title('Energy per Time')
        plt.xlabel(r'Time $(\times \tau)$')
        plt.ylabel(r'Energy $(\times \epsilon)$')
        plt.legend(['Kinetic Energy', 'Potential Energy', 'Total Energy'])
        plt.savefig(f'figures/{self.file_name_base}_energies_plot.jpg')
        plt.show()

    def plot_auto_correlation_of_velocites(self):
        if self._relaxation_index is None:
            self.calculate_autocorrelation_of_velocities()
        c_v_s = self.calculate_autocorrelation_of_velocities()
        plt.plot(self.reduced_times_for_plotting, c_v_s)
        plt.plot(self.reduced_times_for_plotting[self._relaxation_index], c_v_s[self._relaxation_index],
                 'go', label=f'Relaxation threshold at ({round(self.reduced_times_for_plotting[self._relaxation_index],2)}, {round(c_v_s[self._relaxation_index],2)})')
        plt.plot(self.reduced_times_for_plotting, np.exp(-self.reduced_times_for_plotting/self.reduced_times_for_plotting[self._relaxation_index]), 'g--')
        plt.title(f'Autocorrelation of Velocities per Time')
        plt.xlabel(r'Time $(\times \tau)$')
        plt.ylabel(r'$C_v$')
        plt.legend()
        plt.savefig(f'figures/{self.file_name_base}_c_v_plot.jpg')
        plt.show()

    def print_data(self, return_in_original_units=False):
        """
        print the values corresponding to the properties of the system rounded by {r} decimal numbers.
        """
        r = 5
        if return_in_original_units==False:
            self.mean_and_error_of_energy()
            self.mean_and_error_of_kinetic_energy()
            self.mean_and_error_of_potential_energy()
            self.mean_and_error_of_left_side_atoms_number()
            self.mean_and_error_of_pressure()
            self.mean_and_error_of_temperature()
            self.calculate_autocorrelation_of_velocities()
            self._calculate_equilibrium_index()
            
            reduced_units_data =\
            {
            'Time average of total energy': f'{round(self._energy,r)}'+' ± '+f'{round(self._energy_error,r)}'+' ×ϵ',
            'Time average of kinetic energy': f'{round(self._kinetic_energy,r)}'+' ± '+f'{round(self._kinetic_energy_error,r)}'+' ×ϵ',
            'Time average of potential energy': f'{round(self._potential_energy,r)}'+' ± '+f'{round(self._potential_energy_error,r)}'+' ×ϵ',
            'Time average of pressure': f'{round(self._pressure,r)}'+' ± '+f'{round(self._pressure_error,r)}'+' ×ϵ/σ^2 ',
            'Time average of temperature': f'{round(self._temperature,r)}'+' ± '+f'{round(self._temperature_error,r)}'+' ×ϵ/k_B ',
            'Time average of fraction of left side atoms number': f'{round(self._number_of_left_side_atoms/self.N,r)}'+' ± '+f'{round(self._number_of_left_side_atoms_error/self.N,r)}',
            'Relaxation time': f'{round(self.reduced_times_for_plotting[self._relaxation_index],r)}'+' ×τ',
            'Equilibrium time': f'{round(self.reduced_times_for_plotting[self.equilibrium_index],r)}'+' ×τ',
            }
            print(f':::Properties of Gas in Reduced Units:::\n::All values are rounded by {r} decimal numbers; So some values may have been vanished.::')
            for key, value in reduced_units_data.items():
                print(key, ' = ', value)
            # with open(f'{self.file_name_base}_INFO.txt', 'w') as file: 
                # file.writelines("% s = % s\n" % data for data in reduced_units_data.items())

        else:
            E, E_err = self.mean_and_error_of_energy(return_in_original_units=True)
            K, K_err = self.mean_and_error_of_kinetic_energy(return_in_original_units=True)
            U, U_err = self.mean_and_error_of_potential_energy(return_in_original_units=True)
            N_l, N_l_err = self.mean_and_error_of_left_side_atoms_number()
            P, P_err = self.mean_and_error_of_pressure(return_in_original_units=True)
            T, T_err = self.mean_and_error_of_temperature(return_in_original_units=True)
            self.calculate_autocorrelation_of_velocities()
            self._calculate_equilibrium_index()
            
            data =\
            {
            'Time average of total energy': f'{round(E,r)}'+' ± '+f'{round(E_err,r)}'+' J',
            'Time average of kinetic energy': f'{round(K,r)}'+' ± '+f'{round(K_err,r)}'+' J',
            'Time average of potential energy': f'{round(U,r)}'+' ± '+f'{round(U_err,r)}'+' J',
            'Time average of pressure': f'{round(P,r)}'+' ± '+f'{round(P_err,r)}'+' Pa',
            'Time average of temperature': f'{round(T,r)}'+' ± '+f'{round(T_err,r)}'+' K',
            'Time average of fraction of left side atoms number': f'{round(N_l/self.N,r)}'+' ± '+f'{round(N_l_err/self.N,r)}',
            'Relaxation time': f'{round(self.tau*self.reduced_times_for_plotting[self._relaxation_index],r)}'+' s',
            'Equilibrium time': f'{round(self.tau*self.reduced_times_for_plotting[self.equilibrium_index],r)}'+' s',
            }
            print(f':::Properties of Gas in SI Units:::\n::All values are rounded by {r} decimal numbers; So some values may have been vanished.::')
            for key, value in data.items():
                print(key, ' = ', value)
            # with open(f'{self.file_name_base}_INFO.txt', 'w') as file: 
                # file.writelines("% s = % s\n" % data_ for data_ in data.items())