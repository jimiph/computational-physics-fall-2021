from os import write
import numpy as np
import time
import csv


def exp_disturb(a=1, b=1):
    """
    By this funciton you can make an exponential distribution. g(y) = a*exp(-y/b)
    """
    x = np.random.rand()
    y = -b * np.log(1-x/(a*b))
    return y

def simple_integ(num_sample):
    """
    This function uses simple sampling method to calculate I=int(exp(-x^2)) in [0, 2].
    """
    rand_x_values = np.array([2*np.random.rand() for x in range(num_sample)])
    I_values = np.exp(-np.square(rand_x_values))
    integ_value = np.sum(I_values)/num_sample
    error = np.std(I_values)/np.sqrt(num_sample)
    return 2*integ_value, error

def important_integ(num_sample):
    """
    This function uses important sampling method to calculate I=int(exp(-x^2)) in [0, 2].
    """
    def f(x):
        return np.exp(-np.square(x))
    def g(x):
        return np.exp(-x)
    rand_x_values = 3*np.ones(num_sample)
    i = 0
    while rand_x_values[i] == 3 and i<num_sample-1:
        x = exp_disturb()
        if x>=0 and x<=2:
            rand_x_values[i] = x
            i = i+1

    I_values = np.array([f(x)/g(x) for x in rand_x_values])
    integ_value = np.sum(I_values)/num_sample
    error = np.std(I_values)/np.sqrt(num_sample)
    return (1-np.exp(-2))*integ_value, error

def make_data():
    """
    This function uses two methods and functions above to make data contains errors and 
    integral values.
    """
    num_sample_values = np.array([10**x for x in range(2, 7)], dtype=int)
    exact_value = 0.882081 * np.ones(len(num_sample_values))

    sts_errs_simple = np.zeros(len(num_sample_values))
    sts_errs_important = np.zeros(len(num_sample_values))
    integ_values_simple = np.zeros(len(num_sample_values))
    integ_values_important = np.zeros(len(num_sample_values))
    time_values_simple = np.zeros(len(num_sample_values))
    time_values_important = np.zeros(len(num_sample_values))

    for i in range(len(num_sample_values)):
        start_time = time.time()
        integ_vlaue, error = simple_integ(num_sample_values[i])
        time_values_simple[i] = time.time() - start_time
        integ_values_simple[i] = integ_vlaue
        sts_errs_simple[i] = error
    abs_errs_simple = np.abs(exact_value - integ_values_simple)
    for i in range(len(num_sample_values)):
        start_time = time.time()
        integ_vlaue, error = important_integ(num_sample_values[i])
        time_values_important[i] = time.time() - start_time
        integ_values_important[i] = integ_vlaue
        sts_errs_important[i] = error
    abs_errs_important = np.abs(exact_value - integ_values_important)
    old_data = [num_sample_values, exact_value, integ_values_simple, integ_values_important,
    sts_errs_simple, sts_errs_important, abs_errs_simple, abs_errs_important,
    time_values_simple, time_values_important]
    data_sheet = []
    for data in old_data:
        data = list(np.round(data, 6))
        data_sheet = data_sheet + [data]
        
    name_rows = ['number of samples', 'exact value of integral', 
    'value of simple method', 'value of important method',
    'statistical error of simple method', 'statistical error of important method',
    'absolute error of simple method', 'absolute error of important method',
    'runtime of simple method', 'runtime of important method']
    with open ("datas_q1.csv", 'w', newline='') as f:
        writer = csv.writer(f)
        i = 0            
        for row in data_sheet:
            row = [name_rows[i]] + row
            writer.writerow(row)
            i += 1
        




make_data()