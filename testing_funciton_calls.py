from ptr_class import Pointer_like
import numpy as np


def min_F_function(x_list):
    return min(x_list)


def custom_min(x_list):
    min_val = x_list[0]
    for entry in x_list:
        if entry.get() < min_val.get():
            min_val = entry
    return min_val

def attempted_fix():
    np.random.seed(0)
    meta_array = [Pointer_like(ele) for ele in np.random.rand(4)]
    depending_array = [meta_array[0], meta_array[1], custom_min([meta_array[0], meta_array[1]])]

    print([ele.get() for ele in depending_array])
    meta_array[1].set(0)
    print([ele.get() for ele in depending_array])


def attempted_fix_2():
    np.random.seed(0)
    meta_array = [Pointer_like(ele) for ele in np.random.rand(4)]
    depending_array = [meta_array[0], meta_array[1], Pointer_like(min_F_function([meta_array[0], meta_array[1]]))]

    print([ele.get() for ele in depending_array])
    meta_array[1].set(0)
    print([ele.get() for ele in depending_array])

def highlight_issue():
    np.random.seed(0)
    meta_array = [Pointer_like(ele) for ele in np.random.rand(4)]
    depending_array = [meta_array[0], meta_array[1], min_F_function([meta_array[0], meta_array[1]])]

    print([ele.get() for ele in depending_array])
    meta_array[1].set(0)
    print([ele.get() for ele in depending_array])


def min_keeps_index():
    np.random.seed(0)
    meta_array = [Pointer_like(ele) for ele in np.random.rand(4)]
    depending_array = [meta_array[0], meta_array[1], min_F_function([meta_array[0], meta_array[1]])]

    print([ele.get() for ele in depending_array])
    meta_array[0].set(0)
    print([ele.get() for ele in depending_array])
    meta_array[1].set(-1)
    print([ele.get() for ele in depending_array])





def attempted_3():
    np.random.seed(0)
    meta_array = [Pointer_like(ele) for ele in np.random.rand(4)]
    func_eval = Pointer_like(min_F_function([meta_array[0], meta_array[1]]))
    depending_array = [meta_array[0], meta_array[1], func_eval.get()]

    print([ele.get() for ele in depending_array])
    meta_array[1].set(0)
    print([ele.get() for ele in depending_array])


if __name__ == '__main__':
    # min_keeps_index()
    # highlight_issue()
    # attempted_fix()
    # attempted_fix_2()
    attempted_3()
