# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 11:00:21 2023

@author: 19735
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


def compartment_model(X, t, A, E, interval_list):
    """ Function that describes the behavior of the dynamics for the compartment model

    Args:
        X (list): amount of drug (not concentration) in each compartment
        t (list): np.linspace of the time horizon of the simulation
        A (list): list of parameter values for edge weights in the compartmental model
        E (list): list of excretion weights for the system
        interval_list (list): list of ranges to check for dosing time

    Returns:
        dx1dt (float): the value of the derivative of compartment 1 in the model
        dx2dt (float): the value of the derivative of compartment 2 in the model
        dx3dt (float): the value of the derivative of compartment 3 in the model
        dx4dt (float): the value of the derivative of compartment 4 in the model
    """
    # x1, = Brain, x2 = CSF, x3 = Blood, x4 = PT
    # a12 = Brain-CSF, a13 = Brain-Blood across BBB
    # a21 = CSF-Brain, a23 = CSF-Blood across BCSFB
    # a31 = Blood-Brain across BBB, a32 = Blood-CSF across BCSFB, a34 = Blood-PT
    # a43 = PT-Blood
    x1, x2, x3, x4 = X
    a12, a13, a21, a23, a31, a32, a34, a43 = A
    e1, e2, e3, e4 = E

    D = 0
    for interval in interval_list:
        if interval[0] < t < interval[1]:
            D = 1
            break
        else:
            D = 0

    dx1dt = -(a12 + a13 + e1) * x1 + a21 * x2 + a31 * x3
    dx2dt = -(a21 + a23 + e2) * x2 + a12 * x1 + a32 * x3
    dx3dt = -(a31 + a32 + a34 + e3) * x3 + D + a13 * x1 + a23 * x2 + a43 * x4
    dx4dt = -(a43 + e4) * x4 + a34 * x3

    return dx1dt, dx2dt, dx3dt, dx4dt


def solve_compartment_model(X, A, E):
    """ Set the auxiliary variables and runs odeint for the compartmental model

    Args:
        X (list): list of initial values for each compartment
        A (list): list of edge weights between each compartment
        E (list): list of excretion weights

    Returns:
        t: time (x-value)
        drug_con: drug concentration (y-value)
        Shows plots the trajectory of each compartment
    """
    start = 0
    stop = 100
    steps = 200
    dosage_time = 24
    step_size = (stop-start) / steps
    t = np.linspace(start, stop, steps)
    current_time = start
    interval_list = []
    while current_time < stop:
        interval = (current_time - step_size * 1/2, current_time + step_size * 1/2)
        interval_list.append(interval)
        current_time += dosage_time
    drug_con = odeint(compartment_model, X, t, args=(A, E, interval_list))

    plt.plot(t, drug_con)
    plt.ylim([0, 1])
    plt.legend(labels=['Brain', 'CSF', 'Blood', 'PT'])
    plt.xlabel("Time (Days)")
    plt.ylabel("Drug Concentration (mg/L)")
    plt.show()

    return t, drug_con


if __name__ == '__main__':
    X = [0, 0, 1, 0]
    # Physiological  condition
    A = [1, 0.000001, 0.8, 0.05, 0.002, 0.0002, 0.95, 0.95]

    # High Rate condition
    A1 = [1, 0.000005, 0.8, 0.08, 0.005, 0.001, 1, 1]

    E = [0.05, 0.05, 0.15, 0.25]

    solve_compartment_model(X, A, E)
    solve_compartment_model(X, A1, E)
