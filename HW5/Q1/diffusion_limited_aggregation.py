import numpy as np
import matplotlib.pyplot as plt


def diffusion_limited_aggregation(num_sites=20, num_walkers=100, shift_up=2):

    def make_dla():
        """
        This function makes a 2D numpy array to show it.
        """
        cluster = np.zeros((num_walkers+shift_up, num_sites))
        cluster[0] = 1

        def check_neighbors(position):
            x = position[1]
            y = position[0]
            xr = (x+1) % num_sites
            xl = (x-1) % num_sites 
            yu = y+1
            yd = y-1
            check = cluster[y][xr] + cluster[y][xl] + cluster[yu][x] + cluster[yd][x]
            if check!=0:
                return False
            else:
                return True
        #Chooses a neighbor randomly and then return it.
        def rand_neighbor(position):
            x = position[1] 
            y = position[0]
            neighbors = [(y, (x-1)%num_sites), (y, (x+1)%num_sites), (y-1, x), (y+1, x)]
            neighbor = neighbors[np.random.choice([0, 1, 2, 3])]
            return neighbor
        #Finds maximum height of the cluster to initialize the initial height
        #of the next random walker.
        def find_max_height(cluster):
            height = 0
            while np.any(cluster[height]):
                height += 1
            return (height + shift_up)
        
        #Makes the final 2D array for showing by imhsow().
        for walker in range(num_walkers):
            position = (find_max_height(cluster), np.random.randint(num_sites))
            while check_neighbors(position):
                if position[0]<num_walkers:
                    position = rand_neighbor(position)
                else:
                    break
            cluster[position[0]][position[1]] = (walker / 25) + 1
        np.save(f'cluster_walkers={num_walkers}.npy', cluster)
        
        return find_max_height(cluster)

    #Shows the cluster.
    def show_dla():
        #If you want to use a new data, use the code below.
        # max_height = make_dla()
        #If you want to use the saved data, use the code below.
        max_height = 40
        fig, ax = plt.subplots(figsize=(10, 10))
        cluster = np.load(f'cluster_walkers={num_walkers}.npy')
        ax.imshow(cluster[0:max_height], cmap='jet')
        ax.invert_yaxis()
        ax.set_title(f'DLA for number of sites = {num_sites} and number of walkers={num_walkers}')
        plt.show()

    show_dla()

diffusion_limited_aggregation()























































