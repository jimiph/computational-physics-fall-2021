import numpy as np
import matplotlib.pyplot as plt


def saw(num_steps, is_save=False):
    def count_saw(num_steps):
        #I use a binary garid to check that path is SAW or not.
        grid = np.zeros((2*num_steps+1, 2*num_steps+1), dtype=int)

        #This is the orging and the starting point of the walker.
        grid[num_steps][num_steps] = 1

        #This counts the total number of SAW paths for a given number of steps.
        total_saw_paths = np.zeros(1, dtype=int)

        def recursive_walk(x, y, counter):
            """
            x, y are current coordinates of the walker. The counter counts the number
            of SAW paths in each call of function.
            """
            #If counter is 0 we counted another SAW path so we should add it up by 1!
            if counter == 0:
                total_saw_paths[0] += 1
            #If counter is not 0 then we yet not find a complete SAW path. So we
            #check for other possible directions.
            else:
                if grid[x+1][y] == 0:
                    grid[x+1][y] = 1
                    recursive_walk(x+1, y, counter-1)
                if grid[x-1][y] == 0:
                    grid[x-1][y] = 1
                    recursive_walk(x-1, y, counter-1)
                if grid[x][y+1] == 0:
                    grid[x][y+1] = 1
                    recursive_walk(x, y+1, counter-1)
                if grid[x][y-1] == 0:
                    grid[x][y-1] = 1
                    recursive_walk(x, y-1, counter-1)
            #Befor returning we should erase the footprint on the current house!
            grid[x][y] = 0

        for h in range(1, num_steps):
            #We start in the downward direction.
            grid[num_steps][num_steps+h] = 1
            grid[num_steps+1][num_steps+h] = 1
            recursive_walk(num_steps+1, num_steps+h, num_steps-h-1)
        #By the symmetry the total saw paths is 8*total_saw_paths+4!
        return 8*total_saw_paths[0]+4

    def analyse_and_plot():
        total_saw_paths_per_n = np.zeros(num_steps, dtype=int)
        ratio_saw_paths_to_ordinary = np.zeros(num_steps)
        n_values = [n for n in range(num_steps)]
        for i in range(1, 1+num_steps):
            count = count_saw(i)
            total_saw_paths_per_n[i-1] = count
            ratio_saw_paths_to_ordinary[i-1] = count / (4**(i))
        print(total_saw_paths_per_n)
        print(ratio_saw_paths_to_ordinary)
        if is_save:
            np.save('total_saw_paths_per_n', total_saw_paths_per_n)
            np.save('ratio_saw_paths_to_ordinary', ratio_saw_paths_to_ordinary)
        plt.figure(1)
        plt.plot(n_values, total_saw_paths_per_n, '-o')
        plt.title('number of total SAW paths per number of steps')
        plt.xlabel('number of steps')
        plt.ylabel('number of total SAW paths')
        plt.figure(2)
        plt.plot(n_values, ratio_saw_paths_to_ordinary, '-o')
        plt.title('ratio of total SAW paths to RW paths per number of steps')
        plt.xlabel('number of steps')
        plt.ylabel('ratio of total SAW paths to RW paths')
        plt.show()
    analyse_and_plot()





saw(20, is_save=True)