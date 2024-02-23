import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


class ptr:
    def __init__(self, obj):
        self.obj = obj

    def get(self):
        return self.obj

    def set(self, obj):
        self.obj = obj

    # define the oeprators
    def __iadd__(self, other):
        self.obj = self.get() + other.get()
        return self

    def __add__(self, other):
        new_ptr = ptr(self.get() + other.get())
        return new_ptr

    def __sub__(self, other):
        new_ptr = ptr(self.get() - other.get())
        return new_ptr

    def __mul__(self, other):
        new_ptr = ptr(self.get() * other.get())
        return new_ptr

    def __truediv__(self, other):
        new_ptr = ptr(self.get() / other.get())
        return new_ptr

    # define comparisons
    def __lt__(self, other):
        return self.get() < other.get()

    def __le__(self, other):
        return self.get() <= other.get()

    def __eq__(self, other):
        return self.get() == other.get()

    def __neq__(self, other):
        return self.get() != other.get()

    def __ge__(self, other):
        return self.get() >= other.get()

    def __gt__(self, other):
        return self.get() > other.get()

    # define casting
    def __float__(self):
        return self.get()


def integrate_this(x, t, S, f):
    dxdt = np.matmul(S, f)
    return dxdt.transpose()


def min_fun(x_list):
    return min(x_list)


def euler_method(metabolites_X_init, h_step, h_span, stoich_S, flux_F):
    """Summary
    
    Parameters
    ----------
    metabolites_X_init : array of ptr
        vector that corresponds to the metabolite level in the network
    time_horizon_t : array of ptr
        holds each timepoint that we'll integrate at
    h_step : double
        how large of a step do we take
    h_span : tuple of double
        range that we integrate over
    stoich_S : 2d matrix of ptr types
        stoichiometric S matrix that encodes the network (with hyper edges)
    flux_F : arrary of ptr
        vector that correspond to the flux values of each edge
    
    Returns
    -------
    vector of ptr
        updated metabolite levels
    """

    t_current = h_span[0]
    t_ending = h_span[1]
    t_array = [t_current]
    metabolite_array = [metabolites_X_init]
    while t_current < t_ending:

        S_values = [[item.get() for item in ele] for ele in stoich_S]
        print(S_values)
        print("###################################")

        metabolites_X_init = np.matmul(stoich_S, flux_F)

        metabolite_array.append(metabolites_X_init)

        t_current += h_step
        t_array.append(t_current)

    return (t_array, metabolite_array) 


if __name__ == '__main__':

    #######################################################
    # ############### ptr example with arrays #############
    #######################################################

    num_points = 10
    x_ptr = [ptr(i) for i in range(0, num_points)]

    # test the desired behavior:
    zero_ptr = ptr(0)
    depending_array = [zero_ptr, zero_ptr, x_ptr[4], zero_ptr, zero_ptr, x_ptr[9]]

    print([ele.get() for ele in depending_array])
    x_ptr[4].set(5)
    print([ele.get() for ele in depending_array])

    ######################################################
    # ######### now with integration? #####################
    ######################################################

    num_edges = 3
    num_nodes = 2

    x_ptr = [ptr(ele) for ele in np.random.rand(num_nodes)]
    x_ptr = [ptr(1), ptr(-1)]
    zero_ptr = ptr(0)
    negative_one_ptr = ptr(-1)
    one_ptr = ptr(1)
    
    # define the stoichiometric array that encodes the network
    stoich_array_ptr = np.array([[zero_ptr, x_ptr[0], negative_one_ptr * x_ptr[1]], [x_ptr[1], zero_ptr, negative_one_ptr]])
    
    # define the fluxes for the edges, these will all be one for now
    flux_ptr = np.array([ptr(1)] * num_edges)

    # set the interval and step size for the integration
    t_current = 0
    t_step = 0.1
    t_final = 1

    # will append to t_array and values for plotting later
    t_array = [t_current]
    values = [[ele.get() for ele in x_ptr]]
    t_step_ptr = ptr(t_step)
    while t_current < t_final:
        print("######################################")
        print("######################################")
        print("######################################")
        print("######################################")
        print('current_t:', t_current)
        print('stoich_S:', stoich_array_ptr)
        print('flux_ptr:', flux_ptr)
        print('x_ptr:', x_ptr)

        # integration step: x = x + S(x) * f
        # integration_step = x_ptr + np.matmul(stoich_array_ptr, flux_ptr)
        integration_step = x_ptr + np.array([ele * t_step_ptr for ele in np.matmul(stoich_array_ptr, flux_ptr)])
        print("integration_step:", integration_step)

        # update the x_ptr object, changes will be seen in stoich_array_ptr
        for i in range(0, len(integration_step)):
            x_ptr[i].set(integration_step[i].get())

        print('x_ptr values after update:', [ele.get() for ele in x_ptr])
        print("s_values:", [[item.get() for item in ele] for ele in stoich_array_ptr])
        t_current += t_step
        t_array.append(t_current)
        values.append([ele.get() for ele in x_ptr])

    print("######################################")
    print("######################################")
    print("######################################")
    print("######################################")

    plt.plot(t_array, values)
    plt.show()

    # ############################################################
    # sol = odeint(integrate_this, x_ptr_values, t, args=(stoich_array_ptr, flux_ptr))
    # timestep = 0.1
    # time_span = (0, 1)
    # integrated_output = euler_method(x_ptr, timestep, time_span, stoich_array_ptr, flux_ptr)

    # timesteps = integrated_output[0]
    # values = integrated_output[1]

    # plt.plot(timesteps, values)
    # plt.show()