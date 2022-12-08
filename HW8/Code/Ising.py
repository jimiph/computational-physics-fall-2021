import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.stats as sp

class Ising_2d_model:
    """
    for simulating the Ising model we need a lot of attributes and methods that
    we defined them below.
    """

    def __init__(self, length: int, beta: float, num_of_samples=100):
        """
        initializing the made Ising object.
        """
        self.length = length
        self.size = length ** 2
        self.beta = beta
        self.energy_gaps = np.exp(-beta * np.array([-8, -4, 0, 4, 8]))
        self.lattice = np.random.choice([-1, 1], size=(length, length))
        self.num_of_samples = num_of_samples

    def update_beta(self, beta):
        """
        we want to use the final lattice of bet_k as an initial lattice
        for bet_(k+1), so we need to update beta value.
        """
        self.beta = beta
        self.energy_gaps = np.exp(-beta * np.array([-8, -4, 0, 4, 8]))

    def metropolis_decision(self, i, j):
        """
        using metropolis decision with periodic boundary conditions.
        """
        energy_change = 2 * self.lattice[i][j] * (
                self.lattice[(i + 1) % self.length][j] +
                self.lattice[i - 1][j] +
                self.lattice[i][(j + 1) % self.length] +
                self.lattice[i][j - 1]
        )

        energy_gap_index = int((energy_change + 8) / 4)
        return np.random.uniform(0, 1) < self.energy_gaps[energy_gap_index]

    def metropolis_evolve_whole_lattice(self):
        """
        using this function you can evolve the whole of the lattice for one time.
        """
        for n in range(self.length ** 2):
            i = np.random.randint(self.length)
            j = np.random.randint(self.length)
            if self.metropolis_decision(i, j):
                self.lattice[i][j] *= -1

    def total_energy(self):
        """
        calculating the total energy of the lattice.
        """
        energy = 0
        for i in range(self.length):
            for j in range(self.length):
                energy += self.lattice[i][j] * (self.lattice[(i + 1) % self.length][j] +
                                                self.lattice[i - 1][j] +
                                                self.lattice[i][(j + 1) % self.length] +
                                                self.lattice[i][j - 1]
                                                )
        return -1 * energy / 2.0

    def relax_lattice(self):
        """
        evolve the lattce till it reaches a relaxation state.
        """
        relaxed_energies = []
        is_relaxed = False
        limit = np.exp(-2)

        while not is_relaxed:
            for _ in range(self.num_of_samples):
                self.metropolis_evolve_whole_lattice()
                relaxed_energies.append(self.total_energy())
            relaxed_energies_array = np.array(relaxed_energies)
            j = int(len(relaxed_energies_array) / 10)
            auto_corr = ((np.dot(relaxed_energies_array, np.roll(relaxed_energies_array, j))) / len(
                relaxed_energies_array)
                         - np.mean(relaxed_energies_array) ** 2) / np.var(relaxed_energies_array)
            if np.absolute(auto_corr) < limit:
                is_relaxed = True

        return np.array(relaxed_energies)

    def magnetization(self):
        """
        calculating the magnetization of the lattice.
        """
        return np.absolute(np.sum(self.lattice) / self.size)

    def spatial_corr_length(self) -> float:
        """
        calculating the spatial correlation between spins.
        """
        corr_lengths = np.zeros(self.length)
        for i in range(self.length):
            corr_lengths[i] = find_corr_length(self.lattice[i][:])

        return np.mean(corr_lengths)


def find_corr_length(list_1d):
    """
    by giving a 1D array, this function calculates the correlation length.
    """
    n = int(len(list_1d) / 10)
    auto_corr = np.zeros(n)

    for j in range(n):
        # computing auto correlation for j
        auto_corr[j] = (np.dot(list_1d, np.roll(list_1d, j)) / len(list_1d)
                        - np.mean(list_1d) ** 2) / np.var(list_1d)

    # finaly returning the first j for which the auto_cor reaches lower than exp(-1).
    return 1 + len(auto_corr[auto_corr > np.exp(-1)])


def make_data(ising, step, num_of_samples=100):
    """
    make usefull data of energy values in a list.
    """
    energy_s = np.zeros(num_of_samples)
    magnetization_s = np.zeros(num_of_samples)
    spatial_correlation_s = np.zeros(num_of_samples)

    for i in range(num_of_samples):
        for _ in range(step):
            ising.metropolis_evolve_whole_lattice()
        energy_s[i] = ising.total_energy()

        magnetization_s[i] = ising.magnetization()
        spatial_correlation_s[i] = ising.spatial_corr_length()
    return (
        np.mean(energy_s),
        np.std(energy_s, ddof=1) / np.sqrt(num_of_samples),
        np.var(energy_s),
        np.mean(magnetization_s),
        np.std(magnetization_s, ddof=1) / np.sqrt(num_of_samples),
        np.var(magnetization_s),
        np.mean(spatial_correlation_s),
        np.std(spatial_correlation_s, ddof=1) / np.sqrt(num_of_samples)
    )


def run_ising_for_beta_s(length, beta_s, num_of_smaples=100):
    """
    now we run the simulation for varying values of beta.
    """
    ising = Ising_2d_model(length, np.min(beta_s))
    n = len(beta_s)
    mean_energy_s = np.zeros(n)
    mean_energy_err_s = np.zeros(n)
    var_energy_s = np.zeros(n)
    mean_magnetization_s = np.zeros(n)
    mean_magnetization_err_s = np.zeros(n)
    var_magnetization_s = np.zeros(n)
    spatial_correlation_s = np.zeros(n)
    spatial_correlation_err_s = np.zeros(n)

    for i in range(len(beta_s)):
        beta = beta_s[i]
        ising.update_beta(beta)
        ising.relax_lattice()
        energy_s = np.zeros(num_of_smaples)
        # giving the user some information of run situation!
        print("[Info]: Value of beta updated to", beta)

        for j in range(num_of_smaples):
            ising.metropolis_evolve_whole_lattice()
            energy_s[j] = ising.total_energy()
        mean_energy_s[i], mean_energy_err_s[i], var_energy_s[i], \
        mean_magnetization_s[i], mean_magnetization_err_s[i], var_magnetization_s[i], \
        spatial_correlation_s[i], \
        spatial_correlation_err_s[i] = make_data(ising, 10, num_of_smaples)
    chi_s = beta_s * var_magnetization_s
    heat_capacity_s = beta_s ** 2 * var_energy_s

    return mean_energy_s, mean_energy_err_s, \
           mean_magnetization_s, mean_magnetization_err_s, \
           chi_s, heat_capacity_s, spatial_correlation_s, spatial_correlation_err_s


def main(length_s, beta_s, num_of_samples=100):
    """
    by this function we make the final data such as energy, magnetization, ....
    num_of_beta_s : number of beta values between beta_min and beta_max.
    length_s : lengths values of lattices.
    """
    energy_s = np.zeros(shape=(len(length_s), len(beta_s)))
    energy_err_s = np.zeros(shape=(len(length_s), len(beta_s)))
    magnetization_s = np.zeros(shape=(len(length_s), len(beta_s)))
    magnetization_err_s = np.zeros(shape=(len(length_s), len(beta_s)))
    chi_s = np.zeros(shape=(len(length_s), len(beta_s)))
    heat_capacity_s = np.zeros(shape=(len(length_s), len(beta_s)))
    spatial_correlation_s = np.zeros(shape=(len(length_s), len(beta_s)))
    spatial_correlation_err_s = np.zeros(shape=(len(length_s), len(beta_s)))

    for i in range(len(length_s)):
        print(">>>[Info]: Length of the Ising system updated to", length_s[i])
        energy_s[i], energy_err_s[i], \
        magnetization_s[i], magnetization_err_s[i], \
        chi_s[i], heat_capacity_s[i], \
        spatial_correlation_s[i], spatial_correlation_err_s[i] = run_ising_for_beta_s(length_s[i], beta_s,
                                                                                      num_of_samples)
        # saving datas of each length separately.
        np.savetxt(f'energy L={length_s[i]}.csv', energy_s[i], delimiter=',')
        np.savetxt(f'energy_err L={length_s[i]}.csv', energy_err_s[i], delimiter=',')
        np.savetxt(f'magnetization L={length_s[i]}.csv', magnetization_s[i], delimiter=',')
        np.savetxt(f'magnetization_err L={length_s[i]}.csv', magnetization_err_s[i], delimiter=',')
        np.savetxt(f'chi L={length_s[i]}.csv', chi_s[i], delimiter=',')
        np.savetxt(f'heat_capacity L={length_s[i]}.csv', heat_capacity_s[i], delimiter=',')
        np.savetxt(f'spatial_corelation L={length_s[i]}.csv', spatial_correlation_s[i], delimiter=',')
        np.savetxt(f'spatial_correlation_err L={length_s[i]}.csv', spatial_correlation_err_s[i], delimiter=',')
    # saving datas of all lengths to CSV files.
    df_energy = pd.DataFrame(energy_s)
    df_energy.to_csv("energy.csv")
    df_energy_err = pd.DataFrame(energy_err_s)
    df_energy_err.to_csv("energy_err.csv")
    df_magnetization = pd.DataFrame(magnetization_s)
    df_magnetization.to_csv("magnetization.csv")
    df_magnetization_err = pd.DataFrame(magnetization_err_s)
    df_magnetization_err.to_csv("magnetization_err.csv")
    df_chi = pd.DataFrame(chi_s)
    df_chi.to_csv("chi.csv")
    df_heat_capacity = pd.DataFrame(heat_capacity_s)
    df_heat_capacity.to_csv("heat_capacity.csv")
    df_spatial_correlation = pd.DataFrame(spatial_correlation_s)
    df_spatial_correlation.to_csv("spatial_correlation.csv")
    df_spatial_correlation_err = pd.DataFrame(spatial_correlation_err_s)
    df_spatial_correlation_err.to_csv("spatial_correlation_err.csv")


def load_and_plot(length_s, beta_s):
    """
    finaly we load all data files, read the data and plot them as functions of beta.
    """
    energy = pd.read_csv('energy.csv').transpose()
    energy.columns = length_s[:]
    energy.drop('Unnamed: 0', inplace=True)
    energy.index = beta_s[:]

    energy_err = pd.read_csv('energy_err.csv').transpose()
    energy_err.columns = length_s[:]
    energy_err.drop('Unnamed: 0', inplace=True)
    energy_err.index = beta_s[:]

    magnetization = pd.read_csv('magnetization.csv').transpose()
    magnetization.columns = length_s[:]
    magnetization.drop('Unnamed: 0', inplace=True)
    magnetization.index = beta_s[:]

    magnetization_err = pd.read_csv('magnetization_err.csv').transpose()
    magnetization_err.columns = length_s[:]
    magnetization_err.drop('Unnamed: 0', inplace=True)
    magnetization_err.index = beta_s[:]

    chi = pd.read_csv('chi.csv').transpose()
    chi.columns = length_s[:]
    chi.drop('Unnamed: 0', inplace=True)
    chi.index = beta_s[:]

    heat_capacity = pd.read_csv('heat_capacity.csv').transpose()
    heat_capacity.columns = length_s[:]
    heat_capacity.drop('Unnamed: 0', inplace=True)
    heat_capacity.index = beta_s[:]

    spatial_correlation = pd.read_csv('spatial_correlation.csv').transpose()
    spatial_correlation.columns = length_s[:]
    spatial_correlation.drop('Unnamed: 0', inplace=True)
    spatial_correlation.index = beta_s[:]

    spatial_correlation_err = pd.read_csv('spatial_correlation_err.csv').transpose()
    spatial_correlation_err.columns = length_s[:]
    spatial_correlation_err.drop('Unnamed: 0', inplace=True)
    spatial_correlation_err.index = beta_s[:]

    sns.set_theme()

    # plt.errorbar(x=beta_s[:], y=np.array(energy[:]),
    #             yerr=np.array(energy_err[:]), marker='o', ls='-.')
    sns.lineplot(data=energy, markers=True)
    plt.tight_layout()
    plt.xlabel("beta")
    plt.ylabel("<E>")
    plt.title('plot of mean energy per spin')
    plt.savefig("energy_plot.jpg", bbox_inches='tight', dpi=300)
    plt.show()

    # plt.errorbar(x=beta_s[:], y=np.array(magnetization[:]),
    # yerr=np.array(magnetization_err[:]), marker='o', ls='-.')
    sns.lineplot(data=magnetization, markers=True)
    plt.tight_layout()
    plt.xlabel("beta")
    plt.ylabel("<|M|>")
    plt.title('plot of absolute of magnetization')
    plt.savefig("magnetization_plot.jpg", bbox_inches='tight', dpi=300)
    plt.show()

    sns.lineplot(data=chi, markers=True)
    plt.tight_layout()
    plt.xlabel("beta")
    plt.ylabel("chi")
    plt.title('plot of susceptibility')
    plt.savefig("chi_plot.jpg", bbox_inches='tight', dpi=300)
    plt.show()

    sns.lineplot(data=heat_capacity, markers=True)
    plt.tight_layout()
    plt.xlabel("beta")
    plt.ylabel("C_v")
    plt.title('plot of heat capacity per spin')
    plt.savefig("heat_capacity_plot.jpg", bbox_inches='tight', dpi=300)
    plt.show()

    # plt.errorbar(x=beta_s[:], y=spatial_correlation[:],
    # yerr=spatial_correlation_err[:], marker='o', ls='-.', label=f'{np.min(length_s)}')
    sns.lineplot(data=spatial_correlation, markers=True)
    plt.tight_layout()
    plt.legend()
    plt.xlabel("beta")
    plt.ylabel("ksi")
    plt.title('plot of spatial correlation length')
    plt.savefig("spatial_correlation_plot.jpg", bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":
    # length values of lattice we want to simulate.
    length_s = np.array([250, 184, 136, 100])
    num_of_samples, beta_s = 100, np.concatenate((np.linspace(0.2, 0.41, 10, endpoint=False, dtype=np.float16),
                                                  np.linspace(0.41, 0.46, 10, endpoint=False, dtype=np.float16),
                                                  np.linspace(0.46, 0.7, 10, endpoint=True, dtype=np.float16)))
    main(length_s=length_s, beta_s=beta_s, num_of_samples=num_of_samples)
    load_and_plot(length_s=length_s, beta_s=beta_s)


