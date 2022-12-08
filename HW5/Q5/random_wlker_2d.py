import numpy as np
import matplotlib.pyplot as plt


def random_walker_2d(iteration_num=100):
    #This function changes (x, y) of walker randomly; then returns the average 
    #distance of wlker from the origin.
    def find_avg_dist(walk):
        dist = 0
        for i in range(iteration_num):
            x = y = 0
            for i in range(walk):
                x = x + np.random.choice([1, -1])
                y = y + np.random.choice([1, -1])
            dist = dist + np.sqrt(x**2 + y**2)
        return (dist/iteration_num)

    #This funcion uses previous function to plot the log(average distances) per
    #log(time of walking).
    def analyse_and_plot():
        walk_steps = [2**x for x in range(5, 15)]
        dist_values = np.zeros(len(walk_steps))
        for i in range(len(walk_steps)):
            dist_values[i] = find_avg_dist(walk_steps[i]) 
        walk_steps_log = np.log(walk_steps)
        dist_values_log = np.log(dist_values)
        coefs = np.polyfit(walk_steps_log, dist_values_log, 1)
        slope, intercept = coefs[0], coefs[1]
        plt.scatter(walk_steps_log, dist_values_log, s=2, label='data point')
        plt.plot(walk_steps_log, slope*walk_steps_log+intercept, 
         label='fit with equation {:.3}*log(n) + {:.3}'.format(slope, intercept))
        plt.title(f'random walker in 2D. Number of iterations is {iteration_num}')
        plt.xlabel('log(time)')
        plt.ylabel('log(<d>)')
        plt.legend()
        plt.show()
    analyse_and_plot()
        

            
random_walker_2d()

        
