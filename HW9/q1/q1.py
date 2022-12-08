import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp



class solving_RC_circuit:
    """
    We define some useful methods and attributes under this class to simulate.
    """
    def __init__(self, initial_chrage, resistance, voltage, capacity, stop_time,
                 h_start=10**-1, h_stop=10**-5, h_samples_num=100, save=False, load=True):
        """
        We initialize the made object by getting initial charge, resistance, voltage, 
        capacity, stop time of simulation, and bounds of tau.
        """
        self.q_0 = initial_chrage
        self.r = resistance
        self.v = voltage
        self.c = capacity
        self.t_stop = stop_time
        self.h_s_num = h_samples_num
        self.tau_stop = stop_time / (resistance * capacity)
        self.x_0 = -1 + initial_chrage / (capacity * resistance)
        self.h_s = np.exp(np.linspace(np.log(h_start), np.log(h_stop), h_samples_num))
        self.save = save
        self.load = load
        
    def undo_variable_change(self, x, tau, is_it_tau):
        """
        By this function we can undo the variable change that we used to do the simulation.
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
        By this function we can calculate the analytical result.
        """
        if is_it_tau == True:
            tau = t
        else:
            tau = t / (self.r * self.c)
        x = self.x_0 * np.exp(-tau)
        origin_parameters = self.undo_variable_change(x, tau, is_it_tau)
        q = origin_parameters['charge']
        t = origin_parameters['time']
        return  q, t

    def f(self, x):
        """
        This is the same function f(x, t) occured in the Euler's formula.
        """
        return -x

    def euler_result_per_h_s(self):
        """
        By this fucntion we can calculate the Euler's result for each h.
        """
        euler_result_per_h_s = []
        for index, h in enumerate(self.h_s):
            print(f'[Info]: Calculating the Euler result for the h number {index+1} of {self.h_s_num}.')
            euler_result = np.zeros(int(self.tau_stop/h))
            euler_result[0] = self.x_0
            tau_s = np.array([h * i for i in range(int(self.tau_stop/h))])
            for i in range(int(self.tau_stop/h)-1):
                euler_result[i+1] = euler_result[i] + self.f(euler_result[i]) * h
            for j in range(int(self.tau_stop/h)):
                original_variables = self.undo_variable_change(euler_result[j], tau_s[j], True)
                euler_result[j] = original_variables['charge']
            euler_result_per_h_s.append(euler_result)
        return euler_result_per_h_s

    def find_error(self):
        """
        By this function we can calculate the error values due to the Euler's method.
        """
        errors = np.zeros(self.h_s_num)
        euler_result_per_h_s = self.euler_result_per_h_s()
        for i in range(self.h_s_num):
            h = self.h_s[i]
            euler_result = euler_result_per_h_s[i][int(self.tau_stop/h)-1]
            error = np.abs((self.analytic_result(self.tau_stop, True)[0]
                           - euler_result) / self.analytic_result(self.tau_stop, True)[0])
            errors[i] = error
        return errors

    def plotting(self):
        """
        Finally, we can plot all usefull data to watch our simulation results and
        compare it with anaytical result.
        """
        h_s = self.h_s
        t_s = R * C * np.array([h_s[0] * i for i in range(int(self.tau_stop/h_s[0]))])

        if self.load == True:
            euler_result_per_t = np.load('q1_Euler_results.npy')
            analytic_result_per_t = np.load('q1_analytic_results.npy')
            error_s = np.load('q1_Euler_errors.npy')
        if self.load == False:
            euler_result_per_t = self.euler_result_per_h_s()[0]
            analytic_result_per_t = np.zeros(int(self.tau_stop/self.h_s[0]))
            error_s = self.find_error()
        for index, t in enumerate(t_s):
            analytic_result_per_t[index] = self.analytic_result(t, False)[0]
        slope, intercept, r, p, se = sp.linregress(np.log10(h_s), np.log10(error_s))
        plt.scatter(np.log10(h_s), np.log10(error_s), label='Data', s=1)
        plt.xlabel('ln(h)')
        plt.ylabel('ln(error)')
        plt.title('Plot of Euler percent error ')

        plt.plot(np.log10(h_s), slope*np.log10(h_s)+intercept, color='orange',
                label='Fit with equation {:.3}*ln(h)+{:.3}'.format(slope, intercept))
        plt.legend()
        if self.save:
            plt.savefig(f'q1/q1_Euler_error_{h_start}_{h_end}_{h_samples_num}_{stop_time}.jpg')
            np.save('q1_Euler_errors.npy', error_s)
            np.save('q1_Euler_results.npy', euler_result_per_t)
            np.save('q1_analytic_results.npy', analytic_result_per_t)
        plt.show()
        plt.scatter(t_s, euler_result_per_t, s=1, label='Euler Method')
        plt.plot(t_s, analytic_result_per_t, label='Analytic Solution', color='orange')
        plt.xlabel('time')
        plt.ylabel('charge')
        plt.title('Charge of the capacitor as a function of time')
        if self.save:
            plt.savefig(f'q1/q1_Euler_result_h={h_start}_tau_stop={self.tau_stop}.jpg')
        plt.legend()
        plt.show()
        
    
if __name__ == '__main__':
    """
    We make an object from the class 'solving_RC_circuit' and initialize it.
    """
    R = 3 * 10 ** 3
    V = 10
    C = 10 ** (-6)
    Q0 = 0
    stop_time = 21 * 10 ** (-3)
    stop_tau = stop_time / (R * C)
    h_end = 10**-1
    h_start = 10**-5
    h_samples_num = 100
    RC = solving_RC_circuit(initial_chrage=Q0, resistance=R, voltage=V, capacity=C,
    stop_time=stop_time, h_start=h_start, h_stop=h_end, h_samples_num=h_samples_num,
    save=True, load=True)
    RC.plotting()

