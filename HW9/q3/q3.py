import numpy as np
import matplotlib.pyplot as plt



class alg_instability_RC_circuit:
    """
    I make a class with some usefull methods and attributes.
    """
    def __init__(self, initial_chrage, resistance, voltage, capacity, stop_time,
                 tau_step):
        """
        Initializing.
        """
        self.q_0 = initial_chrage
        self.r = resistance
        self.v = voltage
        self.c = capacity
        self.t_stop = stop_time
        self.tau_stop = stop_time / (resistance * capacity)
        self.tau_step_num = int(self.tau_stop / tau_step)
        self.x_0 = -1 + initial_chrage / (capacity * resistance)
        self.h = tau_step

    def undo_variable_change(self, x, tau, is_it_tau):
        """
        Using the given x and tau and returns q and t which are parameters in original units.
        """
        q = self.c * self.v * (1 + x)
        if is_it_tau == True:
            t = self.r * self.c * tau
            return {
                    'charge' : q,
                    'time' : t
                   }
        else:
            t = tau
            return {
                    'charge' : q,
                    'time' : t
                    }

    def analytic_result(self, t, is_it_tau):
        """
        Using the given t or tau returns q and t.
        """
        tau = t
        x0 = self.q_0 / (self.c * self.r) - 1
        x = x0 * np.exp(-tau)
        origin_parameters = self.undo_variable_change(x, tau, is_it_tau)
        q = origin_parameters['charge']
        t = origin_parameters['time']
        return  q, t

    def analytic_solution_complete(self):
        """
        Returns q(t) and t arrays.
        """
        tau_s = np.array([i*self.h for i in range(self.tau_step_num)])
        t_s = np.zeros(self.tau_step_num)
        q_s = np.zeros(self.tau_step_num)
        for i in range(self.tau_step_num):
            q_s[i], t_s[i] = self.analytic_result(tau_s[i], True)
        return q_s, t_s

    def derivative_respect_to_tau(self, x):
        """
        Returns the value of f(x,t) in the original equation of the algorithm.
        """
        return -x
    
    def simulation_result(self):
        """
        Returns simulated q(t) and t arrays.
        """
        x_s = np.zeros(self.tau_step_num)
        x_s[0] = self.x_0
        x_s[1] = x_s[0] + self.h * self.derivative_respect_to_tau(x_s[0])
        for i in range(1, self.tau_step_num-1):
            x_s[i+1] = x_s[i-1] + 2 * self.h * self.derivative_respect_to_tau(x_s[i])
        q_s = np.zeros(self.tau_step_num)
        tau_s = np.array([i*self.h for i in range(self.tau_step_num)])
        t_s = np.zeros(self.tau_step_num)
        for i in range(self.tau_step_num):
            item = self.undo_variable_change(x_s[i], tau_s[i], True)
            q_s[i] = item['charge']
            t_s[i] = item['time']
        return q_s, t_s


if __name__ == '__main__':
    """
    Simulation for some different values of h.
    """
    R = 3 * 10 ** 3
    V = 10
    C = 10 ** (-6)
    Q0 = 0
    stop_time = 21 * 10 ** (-3)
    sample_1 = alg_instability_RC_circuit(initial_chrage=Q0, resistance=R, voltage=V,
    capacity=C, tau_step=0.02, stop_time=stop_time)
    sample_2 = alg_instability_RC_circuit(initial_chrage=Q0, resistance=R, voltage=V,
    capacity=C, tau_step=0.04, stop_time=stop_time)
    sample_3 = alg_instability_RC_circuit(initial_chrage=Q0, resistance=R, voltage=V,
    capacity=C, tau_step=0.08, stop_time=stop_time)  
    sample_analytic = alg_instability_RC_circuit(initial_chrage=Q0, resistance=R, voltage=V,
    capacity=C, tau_step=0.01, stop_time=stop_time)
    plt.figure(figsize=(10, 8))  
    plt.plot(sample_1.simulation_result()[1], sample_1.simulation_result()[0], label='h = 0.02')
    plt.plot(sample_2.simulation_result()[1], sample_2.simulation_result()[0], label='h = 0.04')
    plt.plot(sample_3.simulation_result()[1], sample_3.simulation_result()[0], label='h = 0.08')
    plt.plot(sample_analytic.analytic_solution_complete()[1], sample_analytic.analytic_solution_complete()[0], label='Analytic Solution')
    plt.legend()
    plt.xlabel('time')
    plt.ylabel('charge of the capacitor')
    plt.title(f'Instability of algorithm\n for some varying values of h.')
    plt.show()

        

            
