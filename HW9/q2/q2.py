import numpy as np
import matplotlib.pyplot as plt


class SHM:
    """
    We make a class with some useful attributes and methods to do the simulation.
    """
    def __init__(self, initial_position, initial_velocity, stop_time, time_step):
        """
        First we initialize the made object. Initial position, initial velocity,
        stop time of simulation and time step are given arguments. 
        """
        self.x_0 = initial_position
        self.v_0 = initial_velocity
        self.t_stop = stop_time
        self.h = time_step
        self.step_s_num = int(stop_time / time_step)
        self.t_s = np.array([self.h * i for i in range(self.step_s_num)])
        
    
    def analytic_result(self):
        """
        Return: Two arrays x values and v values.
        """
        x_s =  self.x_0 * np.cos(self.t_s) + self.v_0 * np.sin(self.t_s)
        v_s = -self.x_0 * np.sin(self.t_s) + self.v_0 * np.cos(self.t_s)
        return x_s, v_s

    def euler_result(self):
        """
        Return: Two arrays x values and v values made by Euler's method.
        """
        x_s = np.zeros(self.step_s_num)
        v_s = np.zeros(self.step_s_num)
        x_s[0] = self.x_0
        v_s[0] = self.v_0
        for i in range(self.step_s_num-1):
            x_s[i+1] = x_s[i] + v_s[i] * self.h
            v_s[i+1] = v_s[i] + (-x_s[i]) * self.h
        return x_s, v_s

    def euler_cromer_result(self):
        """        
        Return: Two arrays x values and v values made by Euler-Cromer's method.
        """
        x_s = np.zeros(self.step_s_num)
        v_s = np.zeros(self.step_s_num)
        x_s[0] = self.x_0
        v_s[0] = self.v_0
        for i in range(self.step_s_num-1):
            v_s[i+1] = v_s[i] + (-x_s[i]) * self.h
            x_s[i+1] = x_s[i] + v_s[i+1] * self.h
        return x_s, v_s
        
    def verlat(self):
        """        
        Return: Two arrays x values and v values made by Verlat's method.
        Also we need to define extended initial conditions due to the intrinsic of the 
        algorithm.
        """
        x_s = np.zeros(self.step_s_num+2)
        v_s = np.zeros(self.step_s_num+1)
        x_s[0:2] = self.x_0, self.x_0
        v_s[0] = self.v_0
        for i in range(1, 1+self.step_s_num):
            x_s[i + 1] = 2 * x_s[i] - x_s[i - 1] + (-x_s[i]) * self.h ** 2
            v_s[i] = (x_s[i + 1] - x_s[i - 1]) / (2 * self.h)
        x_s = np.delete(x_s, [0, -1])
        v_s = np.delete(v_s, 0)
        return x_s, v_s

    def velocity_verlat_result(self):
        """        
        Return: Two arrays x values and v values made by velocity Verlat's method.
        """
        x_s = np.zeros(self.step_s_num)
        v_s = np.zeros(self.step_s_num)
        x_s[0] = self.x_0
        v_s[0] = self.v_0
        for i in range(self.step_s_num-1):
            x_s[i + 1] = x_s[i] + v_s[i] * self.h + 0.5 * self.h ** 2 * (-x_s[i])
            v_s[i + 1] = v_s[i] + 0.5 * ((-x_s[i + 1]) + (-x_s[i])) * self.h

        return x_s, v_s

    def beeman_result(self):
        """        
        Return: Two arrays x values and v values made by Beeman's method.
        Also we need to define extended initial conditions due to the intrinsic of the 
        algorithm.
        """
        x_s = np.zeros(self.step_s_num+1)
        v_s = np.zeros(self.step_s_num+1)
        x_s[0:2] = self.x_0
        for i in range(1, self.step_s_num):
            x_s[i + 1] = x_s[i] + v_s[i] * self.h + 1.0 / 6 \
                * (4 * (-x_s[i]) - (-x_s[i - 1])) * self.h ** 2

            v_s[i + 1] = v_s[i] + 1 / 6.0 \
                * ( 2 * (-x_s[i + 1]) + 5 * (-x_s[i]) - (-x_s[i - 1])) * self.h
        x_s = np.delete(x_s, 0)
        v_s = np.delete(v_s, 0)
        return x_s, v_s



if __name__ == '__main__':
    """
    Now we make an object from the class 'SHM' and do the simulation.
    """
    x0 = 1
    v0 = 0
    h = 10 ** -2
    t_stop = 10 * (2 * np.pi)
    step_s_num = int(t_stop / h)
    t_s = np.array([h * i for i in range(step_s_num)])
    shm = SHM(initial_position=x0, initial_velocity=v0, time_step=h, stop_time=t_stop)
    results_of_methods = {  'Analytic Solution' : shm.analytic_result(),
                            'Euler Solution' : shm.euler_result(),
                            'Euler Cromer Solution' : shm.euler_cromer_result(),
                            'Verlat Solution' : shm.verlat(),
                            'Velocity Verlat Solution' : shm.velocity_verlat_result(),
                            'Beeman Solution' : shm.beeman_result()
                         }

    plt.figure(figsize=(10, 5))
    plt.xlabel(r'$t$')
    plt.ylabel(r'$x(t)$')
    plt.title(f'SHM for '+r'$x_0$ '+f'={x0}'+'and '+r'$v_0$'+f'={v0} '+'and '+f'h={h}')
    for value in results_of_methods.values():
        plt.plot(t_s, value[0])
    plt.legend(results_of_methods.keys())
    plt.savefig(f'q2_x_t_{0}_{t_stop}_{h}_x0={x0}_v0={v0}.jpg')
    plt.show()
    for index, item in enumerate(results_of_methods.items()):
        plt.figure(figsize=(6, 5))
        plt.xlabel(r'$X$')
        plt.ylabel(r'$V$')
        plt.title(item[0])
        plt.plot(item[1][0], item[1][1])
        plt.axis('equal')
        plt.savefig(f'q2_v_x_{item[0]}_{0}_{t_stop}_{h}_x0={x0}_v0={v0}.jpg')
        plt.show()