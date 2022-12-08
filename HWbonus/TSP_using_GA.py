import numpy as np
import matplotlib.pyplot as plt


class tsp_using_ga():
    def __init__(self, mutation_p=0.8*10**-2, cross_over_p=0.8, population_num=200, cities_num=15):
        """
        initializing. 

        threshold: I suppose a stable solution when for the last 'threshold' numbers we have 
        the same numbers.

        for probability of selections I use sorting due to the fitness parameter; which in
        our problem is the minimization of the length of the path.
        """
        
        self.mutaion_p = mutation_p
        self.cross_over_p = cross_over_p
        self.population_num = population_num
        self.cities_num = cities_num
        self.space_dimension = 2
        self.threshold = 5
        self.cross_over_cut_index = int(self.cities_num/2)
        self.min_path_s_per_generation = []
        self.min_path_length_s_per_generation = []
        self.current_generation = np.zeros((self.population_num, self.cities_num), dtype=int)
        self.position_s = np.zeros((self.space_dimension, self.cities_num), dtype=np.float16)
        self.initialize_position_s()
        self.distanc_s = np.zeros((self.cities_num, self.cities_num), dtype=np.float16)
        self.determine_distance_s()
        self.initialize_generatoin()
        self.current_path_length_s = np.zeros(self.population_num, dtype=np.float16)
        self.determine_path_length_s()
        self.update_min_length_and_path()
        self.normalized_rank_s = np.zeros(self.population_num)
        self.sort_path_s()
        self.size = 0

    def initialize_position_s(self):
        """
        assign to each city a random position.
        """
        for i in range(self.space_dimension):
            self.position_s[i] = np.random.randn(self.cities_num).astype('float16')

    def determine_distance_s(self):
        """
        determining the distances between cities as a self.cities_num * self.cities_nnum array.
        """
        for axis in range(self.space_dimension):
            self.distanc_s += np.square(np.tile(self.position_s[axis], (self.cities_num, 1))\
                - np.transpose(np.tile(self.position_s[axis], (self.cities_num, 1))))
        self.distanc_s = np.sqrt(self.distanc_s)
    
    def initialize_generatoin(self):
        """
        making an initial population using random permutation for cities.
        """
        for i in range(self.population_num):
            self.current_generation[i] = np.random.permutation(self.cities_num)
        self.initial_generatoin = self.current_generation
    
    def determine_path_length_s(self):
        """
        determining the length of each path in the generation.
        """
        for i in range(self.population_num):
            chro = self.current_generation[i]
            path_length = 0
            for j in range(self.cities_num-1):
                path_length += self.distanc_s[chro[j]][chro[j+1]]
            self.current_path_length_s[i] = path_length
     
    def update_min_length_and_path(self):
        """
        saving the minimum length of path and corresponding order of cities.
        """
        self.determine_path_length_s()
        min_path = np.min(self.current_path_length_s)
        min_path_index = np.where(self.current_path_length_s==min_path)
        self.min_path_length_s_per_generation.append(min_path)
        path = self.current_generation[min_path_index[0]][:]
        self.min_path_s_per_generation.append(path)
        self.size = len(self.min_path_length_s_per_generation)

    def sort_path_s(self):
        """
        sorting the generation and corresponding lengths for determining the probablity of selection.
        """
        self.determine_path_length_s()
        indices = np.argsort(self.current_path_length_s)
        generation_before_sort = self.current_generation
        for i in range(self.population_num):
            self.current_generation[i] = generation_before_sort[indices[i]]
        ranks = np.arange(self.population_num)
        for j in range(1, self.population_num):
            if ~np.all(self.current_generation[j]-self.current_generation[j-1]):
                ranks[j] = ranks[j-1]
            else:
                ranks[j] = ranks[j-1]+1
        self.normalized_rank_s = np.max(ranks) - ranks
        self.normalized_rank_s += 1
        self.normalized_rank_s = self.normalized_rank_s / np.sum(self.normalized_rank_s)
    
    def mutate(self, path):
        """
        mutating each offspring befor transfering it to the next generation.
        """
        p = np.random.rand()
        if p <= self.mutaion_p:
            i, j = np.random.randint(self.cities_num, size=2)
            dummy = path[i]
            path[i] = path[j]
            path[j] = dummy
        return path

    def cross_over(self, path1, path2):
        """
        making cross over. path1 and path2 are both two arrays with size of self.cities_num. 
        """
        chro1 = path1[self.cross_over_cut_index:]
        chro2 = path2[self.cross_over_cut_index:]
        new_path1 = path1
        new_path2 = path2
        def make_new_chro(chro1, path2):
            new_chro1 = np.zeros(len(chro1), dtype=int)
            i = 0
            for ch in path2:
                if ch in chro1:
                    new_chro1[i] = ch
                    i += 1
            return new_chro1

        p = np.random.rand()
        if p <= self.cross_over_p:
            new_chro1 = make_new_chro(chro1, path2)
            new_chro2 = make_new_chro(chro2, path1)
            new_path1[self.cross_over_cut_index:] = new_chro1
            new_path2[self.cross_over_cut_index:] = new_chro2 
        return new_path1, new_path2
    
    def make_two_offsprings(self, path1, path2):
        """
        using two parents we make two offsprings.
        """
        new_path1, new_path2 = self.cross_over(path1, path2)
        new_path1 = self.mutate(new_path1)
        new_path2 = self.mutate(new_path2)
        return new_path1, new_path2
        
    def make_new_generation(self):
        new_generation = self.current_generation
        selection_check = np.zeros(self.population_num, dtype=int)
        while ~np.all(selection_check):
            i, j = np.random.choice(np.arange(self.population_num), 2, p=self.normalized_rank_s)
            selection_check[i] = 1
            selection_check[j] = 1
            new_path1, new_path2 = self.make_two_offsprings(self.current_generation[i], self.current_generation[j])
            new_generation[i][:] = new_path1
            new_generation[j][:] = new_path2
        self.current_generation = new_generation
        self.sort_path_s()
        self.update_min_length_and_path()

    def make_complete_evolution(self, is_show=True):
        while True:
            for _ in range(20):
                self.make_new_generation()
                # print(self.current_generation)
                # print(self.min_path_length_s_per_generation)
            last_min_s = self.min_path_length_s_per_generation[self.size-self.threshold:self.size]
            if np.all((last_min_s-last_min_s[0])==0):
                break
            else: 
                continue
        print(self.min_path_length_s_per_generation)
        last = len(self.min_path_s_per_generation)
        some_path = self.min_path_s_per_generation[last-1][0]
        print('a path with minimum length:', some_path)
        fig, ax = plt.subplots(figsize=(10, 10))
        
        x_s = self.position_s[0]
        y_s = self.position_s[1]
        ax.plot(x_s, y_s, 'ro')
        data = []
        for i in range(self.cities_num):
            data.append(f'{i}')
        for i, num in enumerate(data):
            ax.annotate(num, (x_s[i], y_s[i]))
        x_p_s = []
        y_p_s = []    
        for i in range(self.cities_num):
            x_p_s.append(self.position_s[0][some_path[i]])
            y_p_s.append(self.position_s[1][some_path[i]])
        ax.plot(x_p_s, y_p_s)
        ax.set_title(f'TSP using GA algorithm result.\nPopulation number = {self.population_num}, Cities number = {self.cities_num}, Showen minimum path = {str(some_path)}')
        plt.show()
            
a = tsp_using_ga()
a.make_complete_evolution()