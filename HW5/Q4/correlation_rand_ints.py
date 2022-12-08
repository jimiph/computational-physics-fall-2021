import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss

#This function returns plenty of random integers 0-9 for total_iter_num iterations.
def find_disturb(total_iter_num):
    all_rands = np.random.randint(10, size=total_iter_num)
    selected_rands = []
    for i in range(total_iter_num-1):
        if all_rands[i+1] == 4:
            selected_rands.append(all_rands[i])
    plenty_of_rands = np.zeros(10)
    for i in selected_rands:
        plenty_of_rands[i] += 1
    return plenty_of_rands

#This function plot distribution and std values.
def plot_disturb_and_sigma(total_iter_nums):
    sigma_values = np.zeros(len(total_iter_nums))
    for i in range(len(total_iter_nums)):
        rands = np.arange(10)
        plt.figure(i)
        disturb = find_disturb(total_iter_nums[i])
        sigma_values[i] = np.log(np.std(disturb))
        plt.bar(rands, disturb)
        plt.title(f'unfirom distribution of 0-9 integers. total iteration is {total_iter_nums[i]}.\nfigure {i+1} of {len(total_iter_nums)}')
        plt.xticks(np.arange(10, step=1))
        plt.xlabel('random integers between 0 and 9')
        plt.ylabel('plenty of each integer')
        plt.show()
    total_iter_nums_log = np.array([np.log(x) for x in total_iter_nums])
    plt.scatter(total_iter_nums_log, sigma_values, s=2)
    slope, intercept, r, pp, se = ss.linregress(total_iter_nums_log, sigma_values)
    plt.plot(total_iter_nums_log, slope*total_iter_nums_log+intercept,
    label='fit with eqation {:.3}*log(N)+{:.3}'.format(slope, intercept), color='r')
    plt.legend()
    plt.show()


def main():
    if __name__ == '__main__':
        total_iter_nums = [2**x for x in range(11, 25)]
        plot_disturb_and_sigma(total_iter_nums)

main()