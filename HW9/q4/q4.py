import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

class logistic_map:
    """
    logistic map class with some useful methods and attributes.
    """

    def __init__(self, r_start, r_end, r_num_s, x_0_num_s, n_max, x_lim_s=None, save=False, load=False):
        """
        Initializing the made object.
        """
        self.r_start = r_start
        self.r_end = r_end
        self.r_num_s = r_num_s
        self.x_0_num_s  = x_0_num_s
        self.n_max = n_max
        self.save = save
        self.load = load
        self.x_lim_s = x_lim_s

    def logistic_map(self, r, x_0_s):
        """
        Returns the logistic result for some given r and inital x_0 values.
        Number of iterations is n_max which is given by the user.
        """
        ans_s = x_0_s
        for i in range(self.n_max):
            ans_s = (4 * r) * ans_s * (1 - ans_s)
        return ans_s

    def bifurcation_plot(self):
        """
        Plotting the bifurcation plot.
        """
        x_0_s = np.random.rand(self.x_0_num_s)
        r_s = np.linspace(self.r_start, self.r_end, self.r_num_s)
        x_per_r = np.zeros((self.r_num_s, self.x_0_num_s))
        file_name = f'q4_r_start={self.r_start}r_end={self.r_end}r_num_s={self.r_num_s}x_0_num_s={self.x_0_num_s}n_max={self.n_max}.npy'

        if self.load == True:
            x_per_r = np.load(file_name)
        if self.load == False:
            for index, r in enumerate(r_s):
                print(f'[Info]: This is the r number {index+1} of the {self.r_num_s} for making data.')
                x_per_r[index][:] = self.logistic_map(r=r, x_0_s=x_0_s)
            np.save(file_name, x_per_r)

        plt.figure(figsize=(10, 5))
        for index in range(self.r_num_s):
            print(f'[Info]: This is the r number {index+1} of the {self.r_num_s} for plotting.')
            plt.scatter(r_s[index]*np.ones(self.x_0_num_s), x_per_r[index][:], color='blue')       
            plt.plot(r_s, 0.5*np.ones(self.r_num_s), color='orange')
        plt.xlabel(r'$r$')
        plt.ylabel(r'$x$')
        if self.x_lim_s != None:
            plt.ylim(self.x_lim_s)
        if self.save:
            plt.savefig(f'q4/q4_{self.r_start}_{self.r_end}_{self.r_num_s}_{self.x_0_num_s}_{self.n_max}.jpg')
        plt.show() 


if __name__ == '__main__':
    """
    Uncomment one of the code scopes below to simulate for about some r_n.
    """
    # logistic_map = logistic_map(r_start=0, r_end=1, r_num_s=10**4, x_0_num_s=100, n_max=10**3)
   
    # logistic_map = logistic_map(r_start=0.25-2*10**(-6), r_end=0.25+2*10**(-6), r_num_s=9, x_0_num_s=100, n_max=10**7)
   
    # logistic_map = logistic_map(r_start=0.75-1*10**(-4), r_end=0.75+1*10**(-4), r_num_s=20, x_0_num_s=100, n_max=10**6)
   
    # logistic_map = logistic_map(r_start=0.8623+69*10**(-6), r_end=0.8623+75*10**(-6),
                                # r_num_s=7, x_0_num_s=100, n_max=10**6, x_lim_s=(0.43900, 0.44100))
   
    # logistic_map = logistic_map(r_start=0.8860+22*10**(-6), r_end=0.8860+23*10**(-6),
                                # r_num_s=20, x_0_num_s=100, n_max=10**6, x_lim_s=(0.5232, 0.5240))
   
    # logistic_map = logistic_map(r_start=0.8911, r_end=0.8911+5*10**(-6),
                                # r_num_s=20, x_0_num_s=100, n_max=10**6, x_lim_s=(0.4895, 0.4916))
   
    # logistic_map = logistic_map(r_start=0.892189, r_end=0.892190,
                                #   r_num_s=10, x_0_num_s=100, n_max=10**6, x_lim_s=(0.5036, 0.5039))
   
    # logistic_map = logistic_map(r_start=0.892445+21*10**(-6), r_end=0.892500+24*10**(-6),
                                # r_num_s=55, x_0_num_s=100, n_max=10**6, x_lim_s=(0.3656, 0.3660))
   
    # logistic_map = logistic_map(r_start=0.8886+5*10**(-5), r_end=0.8886+7*10**(-5),
                                # r_num_s=20, x_0_num_s=100, n_max=10**5, x_lim_s=(0.5-4*10**-5,
                                                                                # 0.5+4*10**-5))
   
    # logistic_map = logistic_map(r_start=0.89166+4*10**(-6), r_end=0.89166+10*10**(-6),
    #                             r_num_s=7, x_0_num_s=100, n_max=10**6, x_lim_s=(0.5-3*10**-5,
    #                                                                              0.5+3*10**-5))    
  
    logistic_map.bifurcation_plot()


