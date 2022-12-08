from matplotlib import legend
import numpy as np
import matplotlib.pyplot as plt


def metropolis_alg(num_selections=10**5, is_plot_gaussian=False, is_find_corr_lngths=False, sigma=1, mu=0):
    def gen_gaussian(sigma=sigma, mu=mu, step=1):
        def gaussian(x):
            return (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-0.5*((x-mu)/sigma)**2)
        
        x = mu
        rand_nums = np.zeros(num_selections)
        count_acception = 0

        # for i in range(num_selections):
        #     y = x + step * np.random.uniform(-1, 1)
        #     if gaussian(y) / gaussian(x) >= 1:
        #         count_acception += 1
        #         x = y
        #     else:
        #         x = np.random.choice([y, x], p=[gaussian(y)/gaussian(x), 1-gaussian(y)/gaussian(x)])
        #     rand_nums[i] = x

        for i in range(num_selections):
            y = x + step * np.random.uniform(-1, 1)
            if np.random.uniform(0, 1) < gaussian(y) / gaussian(x):
                count_acception += 1
                x = y
            rand_nums[i] = x



        return rand_nums, count_acception/num_selections

    def plot_gaussian():
        rand_nums = gen_gaussian()[0]
        plt.hist(rand_nums, density=1, bins=40, label='Histogram of distribution')
        x = np.linspace(-4*sigma, 4*sigma, num=10**3)
        y = (1/(np.sqrt(2*np.pi)*sigma))*np.exp(-0.5*((x-mu)/sigma)**2)
        plt.plot(x, y, label='Gussian fit')
        plt.xlabel('x')
        plt.ylabel('p(x)')
        plt.title(f'test histogram of gaussian generator using Metropolis algorithm\nsigma={sigma}, mean={mu}\nN={num_selections}')
        plt.legend()
        plt.show()
    

    def calculate_correlation_length(rand_nums):
        corr_lngth = np.zeros(40)
        for i in range(40):
            s = 0
            for j in range(num_selections-i):
                s += rand_nums[j] * rand_nums[i+j]
            s /= num_selections - i
            corr_lngth[i] = (s - np.mean(rand_nums) ** 2) / np.var(rand_nums)
        return len(corr_lngth[corr_lngth > np.exp(-1)])
    
    if is_plot_gaussian == True:
        plot_gaussian()
    
    if is_find_corr_lngths == True:
        # steps = [11.3, 5.65, 3.75, 2.73, 2.10, 1.56, 1.11, 0.72, 0.35]
        # steps = [23.6, 11.87, 8, 5.95, 4.75, 3.875, 3.228, 2.68, 2.205]
        steps = [0.5, 1, 1.58, 2.2, 2.95, 3.88, 5.3, 7.85, 16]
        # steps = [1]
        rand_datas = [None] * len(steps)
        acception_rates = [None] * len(steps)
        corr_lngths = [None] * len(steps)
        for i in range(len(steps)):
            rand_datas[i], acception_rates[i] = gen_gaussian(step=steps[i])
            corr_lngths[i] = calculate_correlation_length(rand_datas[i])
        
        print('accepstion rates:', acception_rates)
        np.savetxt('acception_rates_q3.csv', acception_rates)
        print('correlation lengths:', corr_lngths)
        np.savetxt('corr_lengths_q3.csv', corr_lngths)
        plt.plot(acception_rates, corr_lngths, '-o')
        plt.title('correlation length per acception rate.')
        plt.xlabel('acception rate')
        plt.ylabel('correlation length')
        plt.show()


metropolis_alg(is_find_corr_lngths=True)
# metropolis_alg(is_plot_gaussian=True)