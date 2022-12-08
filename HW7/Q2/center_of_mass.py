import numpy as np


def center_of_mass(num_samples = 10**6, num_ensembles=10):
    R = 15

    def rho(r, theta):
        return np.round((3*R + r*np.cos(theta)), 6)
    x_cm_s = np.zeros(num_ensembles)
    y_cm_s = np.zeros(num_ensembles)
    z_cm_s = np.zeros(num_ensembles)
    for i in range(num_ensembles):
        mass = 0
        zm_mass = 0
        ym_mass = 0
        xm_mass = 0
        for j in range(num_samples):
            r, theta, phi = R*np.random.rand(), np.pi*np.random.rand(), 2*np.pi*np.random.rand()
            mass += rho(r, theta) * np.square(r) * np.sin(theta)
            zm_mass += rho(r, theta) * r * np.cos(theta) * np.square(r) * np.sin(theta)
            ym_mass += rho(r, theta) * r * np.sin(theta) * np.square(r) * np.sin(theta) * np.sin(phi)
            xm_mass += rho(r, theta) * r * np.sin(theta) * np.square(r) * np.sin(theta) * np.cos(phi)
        z_cm_s[i] = np.round(zm_mass / mass, 6)
        y_cm_s[i] = np.round(ym_mass / mass, 6)
        x_cm_s[i] = np.round(xm_mass / mass, 6)
    x_cm = np.round(np.average(x_cm_s),6)
    y_cm = np.round(np.average(y_cm_s),6)
    z_cm = np.round(np.average(z_cm_s),6)
    x_cm_err = np.round(np.var(x_cm_s)/np.sqrt(num_ensembles),6)
    y_cm_err = np.round(np.var(y_cm_s)/np.sqrt(num_ensembles),6)
    z_cm_err = np.round(np.var(z_cm_s)/np.sqrt(num_ensembles),6)

    print(f'Center of mass is located at\n x={x_cm}, error={x_cm_err}\ny={y_cm}, error={y_cm_err}\nz={z_cm}, error={z_cm_err}')



center_of_mass()