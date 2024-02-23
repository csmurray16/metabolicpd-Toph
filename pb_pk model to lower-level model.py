# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:57:10 2023

@author: 19735
"""

from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import compartmentalized_pb_pk_model as pb_pk


def interpolate_brain_concentration(t, brain_concentration):
    """Perform linear interpolation for drug concentration in the Brain compartment.

    Args:
        t (numpy array): Time points at which the drug concentration is known.
        brain_concentration (numpy array): Drug concentration in the Brain compartment at corresponding time points.

    Returns:
        t_interp (numpy array): Time points for interpolated drug concentration.
        brain_concentration_interp (numpy array): Interpolated drug concentration in the Brain compartment.
    """
    f = interpolate.interp1d(t, brain_concentration, kind='linear')

    t_interp = np.arange(0, t[-1], 0.1)
    brain_concentration_interp = f(t_interp)

    return t_interp, brain_concentration_interp


if __name__ == "__main__":
    X = [0, 0, 1, 0]
    A = [1, 0.000001, 0.8, 0.05, 0.002, 0.0002, 0.95, 0.95]
    A1 = [1, 0.000005, 0.8, 0.08, 0.005, 0.001, 1, 1]
    E = [0.05, 0.05, 0.15, 0.25]
    t, pb_pk_solved = pb_pk.solve_compartment_model(X, A, E)
    t, pb_pk_solved1 = pb_pk.solve_compartment_model(X, A1, E)

    # Extracting the concentrations of each compartment from the results
    # Physiological  condition
    brain_concentration = pb_pk_solved[:, 0]
    csf_concentration = pb_pk_solved[:, 1]
    blood_concentration = pb_pk_solved[:, 2]
    pt_concentration = pb_pk_solved[:, 3]

    # High Rate condition
    brain_concentration1 = pb_pk_solved1[:, 0]
    csf_concentration1 = pb_pk_solved1[:, 1]
    blood_concentration1 = pb_pk_solved1[:, 2]
    pt_concentration1 = pb_pk_solved1[:, 3]

    plt.plot(t, brain_concentration, label='Brain (High Rate)')
    plt.plot(t, csf_concentration, label='CSF (High Rate)')
    plt.plot(t, blood_concentration, label='Blood (High Rate)')
    plt.plot(t, pt_concentration, label='PT (High Rate)')

    plt.plot(t, brain_concentration1, label='Brain (Physiological)')
    plt.plot(t, csf_concentration1, label='CSF (Physiological)')
    plt.plot(t, blood_concentration1, label='Blood (Physiological)')
    plt.plot(t, pt_concentration1, label='PT (Physiological)')

    plt.ylim([0, 1])
    plt.legend()
    plt.xlabel("Time (Days)")
    plt.ylabel("Drug Concentration (mg/L)")
    plt.show()

    # Linear interpolation for drug concentration in the Brain
    t_interp, brain_concentration_interp = interpolate_brain_concentration(t, brain_concentration)
    t_interp1, brain_concentration_interp1 = interpolate_brain_concentration(t, brain_concentration1)

    plt.plot(t, brain_concentration, 'o', label='Brain (High Rate) - Known Data')
    plt.plot(t_interp, brain_concentration_interp, '-', label='Brain (High Rate) - Interpolated')

    plt.plot(t, brain_concentration1, 'o', label='Brain (Physiological) - Known Data')
    plt.plot(t_interp1, brain_concentration_interp1, '-', label='Brain (Physiological) - Interpolated')

    plt.ylim([0, 0.006])
    plt.legend()
    plt.xlabel("Time (Days)")
    plt.ylabel("Drug Concentration (mg/L)")
    plt.show()


