import numpy as np
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import size

def gaussian_disturb(sigma=1, iteration_num=10**5, isTest=False):
    """
    By this funciton you can make a gaussain distribution and plot it as a test.
    """
    def gaussian_gen():
        """
        This function return a random gaussian number. The calulations below is 
        explained in the report.
        """
        def make_random_polar():
            x1 = np.random.rand()
            x2 = np.random.rand()
            r = np.sqrt(-2*np.square(sigma)*np.log(1-x1))
            t = 2*np.pi*x2
            return r, t

        r, t = make_random_polar()
        x = r*np.cos(t)
        y = r*np.sin(t)

        return x

    def test_gaussian_gen():
        """
        This function plots a distribution of a lot of gaussian random numbers.
        The total number is equal to iteration_num; also i have used uniform steps
        equal to 0.01; but if you like you can change it.
        """
        step = 0.01
        gaussian_random_nums = np.zeros(iteration_num)
        for i in range(iteration_num):
            gaussian_random_nums[i] = gaussian_gen()

        round_gaussian_random_nums = gaussian_random_nums
        for i in range(iteration_num):
            p = round(gaussian_random_nums[i]/step)
            round_gaussian_random_nums[i] = float(step*p)
        left_bound = np.min(round_gaussian_random_nums)
        right_bound = np.max(round_gaussian_random_nums)
        disturb = np.zeros(1+int((right_bound-left_bound)/step))
        nums = [left_bound+i*step for i in range(1+int((right_bound-left_bound)/step))]
        for x in round_gaussian_random_nums:
            disturb[int((x-left_bound)/step)] += 1
        plt.bar(nums, disturb)
        plt.xlabel('Numbers')
        plt.ylabel('Distribution')
        plt.title(f'Test of my gaussian generator.\n Number of iterations is {iteration_num}\n sigma={sigma}')
        plt.show()
    if isTest == True:
        test_gaussian_gen()

    if isTest == False:
        gaussian_gen()


gaussian_disturb(isTest=True)