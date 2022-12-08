import numpy as np
import matplotlib.pyplot as plt

def clt_test(N_values = [1000], iteration_num = 10**5):
    def gen_sum_rand(N):
        """
        This function gets N; then firstly generate N random integers between 0-9 and 
        finaly return the sum of them. 
        """
        rand_sum = 0
        for i in range(N):
            rand_sum = rand_sum + np.random.randint(10)
        return rand_sum

    def find_disturb(N):
        """
        This function call the function gen_sum_rand() for iteration_num times. Then
        return the disturibution.
        """
        size = 1 + 9*N
        disturb = np.zeros(size)
        for i in range(iteration_num):
            disturb[gen_sum_rand(N)] += 1
        return disturb

    def plot_disturb_and_variance():
        """
        This function plot the histogram of the distributions. The distributions 
        are being made by calling the function find_disturb() for each values of
        N_values.
        """
        sigma_values = np.zeros(len(N_values))
        for i in range(len(N_values)):
            plt.figure(i)
            disturb = find_disturb(N_values[i])
            # sigma_values[i] = np.log(np.std(disturb))
            right_bound = len(disturb)-1
            left_bound = 0
            while disturb[right_bound]==0:
                right_bound -= 1

            while disturb[left_bound]==0:
                left_bound += 1
            
            rand_sum_values = np.arange(left_bound, right_bound)
            plt.bar(rand_sum_values, disturb[left_bound:right_bound])
            plt.title(f'CLT test. N={N_values[i]}. number of iterations = {iteration_num}\nfigure {i+1} of {len(N_values)}')
            plt.xlabel('random sum vlaues')
            plt.ylabel('plenty of each sum value')
            plt.show()

    plot_disturb_and_variance()


clt_test()