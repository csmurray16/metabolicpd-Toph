# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 10:35:05 2023

@author: 19735
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


def compartment_model(x, t, interval_list):
    # x1, x2, x3, x4 = brain, CSF, blood, PT
    x1, x2, x3, x4 = x[0], x[1], x[2], x[3]

    # numbers taken from High Rate Condition
    a12 = 1
    a13 = 0.000005
    a21 = 0.8
    a23 = 0.08
    a31 = 0.005
    a32 = 0.001
    a34 = 1
    a43 = 1

    # constants(mg/L)
    E1 = 0.05
    E2 = 0.05
    E3 = 0.15
    E4 = 0.25

    # intake
    for interval in interval_list:
        print(interval)
        if t > interval[0] and t < interval[1]:
            D = 1
            break
        else:
            D = 0

    dx1dt = -(a12 + a13 + E1)*x1 + a21*x2 + a31*x3
    dx2dt = -(a21 + a23 + E2)*x2 + a12*x1 + a32*x3
    dx3dt = -(a31 + a32 + a34 + E3)*x3 + D + a13*x1 + a23*x2 + a43*x4
    dx4dt = -(a43 + E4)*x4 + a34*x3
    return dx1dt, dx2dt, dx3dt, dx4dt


x = [0, 0, 1, 0]
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
solved = odeint(compartment_model, x, t, args=(interval_list,))


plt.plot(t, solved)
plt.ylim([0, 2])
plt.legend(labels=['Brain', 'CSF', 'Blood', 'PT'])
plt.xlabel("Time (Days)")
plt.ylabel("Drug Concentration (mg/L)")
plt.show()
