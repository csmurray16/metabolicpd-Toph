from ptr_class import Pointer_like
import pandas as pd
import numpy as np
import re


def extract_metabolites(network_df):
    unique_entries = np.unique(network_df[['tail', 'head']].values)

    unique_metabolites = set()
    for entry in unique_entries:
        elements = entry.split(', ')
        for ele in elements:
            unique_metabolites.add(ele)

    return unique_metabolites


def create_meta_lookup(unique_metabolite_list):
    metabolite_lookup_dict = {}
    index_counter = 0
    for entry in unique_metabolite_list:
        metabolite_lookup_dict[entry] = index_counter
        index_counter += 1

    return metabolite_lookup_dict


def find_source_sink_index(meta_lookup_dict):

    source_sink_dict = {'source':[], 'sink':[]}
    actual_meta_list = []
    for key in meta_lookup_dict.keys():
        sink_res = re.search('^[e]+[0-9]$', key)  # starts with e, ends with ##...
        source_res = re.search('^[s]+[0-9]$', key)
        if sink_res is not None:
            source_sink_dict['sink'].append(meta_lookup_dict[key])
        elif source_res is not None:
            source_sink_dict['source'].append(meta_lookup_dict[key])
        else:
            actual_meta_list.append(key)

    return source_sink_dict, actual_meta_list


def F_min_func(x_list):
    return min(x_list)


def build_S_matrix(net_df, unique_metabolites, meta_lookup, source_sink_dict):
    # get the number of non-source/sinks nodes
    num_source_sink = 0
    for val in source_sink_dict.values():
        # count the number of indices that correspond to source/sink
        num_source_sink += len(val)

    # general information
    num_metabolites = len(unique_metabolites) - num_source_sink
    num_edges = len(net_df['tail'].values)

    # create an array of pointer likes for the metabolites
    metabolite_array = [Pointer_like(ele) for ele in np.random.rand(num_metabolites)]
    negative_one = Pointer_like(-1)  # create negative one for ease of use

    edge_columns_df = net_df[['tail', 'head']]

    S_matrix = []
    for row in edge_columns_df.itertuples():
        # get the substrates and products as a list
        substrates = row.tail.split(', ')
        products = row.head.split(', ')

        # how many subs or prods do we have? needed for hyperedge determination
        num_subs = len(substrates)
        num_prods = len(products)

        # what type of edge we have determines the format of the col in S
        current_col = [Pointer_like(0) for i in range(0, num_metabolites)]
        if (num_subs > 1 or num_prods > 1):  # hyperedge!
            index_of_subs_in_meta = []
            index_of_prods_in_meta = []
            for sub in substrates:
                index_of_subs_in_meta.append(actual_meta_lookup[sub])
            for prod in products:
                index_of_prods_in_meta.append(actual_meta_lookup[prod])

            # big question: if we call min function here, does it break
            # the pointer-like behavior? Can we make a pointer_like of a function call?

        else:  # not a hyperedge!
            if re.search('^[s]+[0-9]$', substrates[0]) is not None:  # source term
                index_of_prod_in_meta = actual_meta_lookup[products[0]]
                current_col[index_of_prod_in_meta] = Pointer_like(1)
            elif re.search('^[e]+[0-9]', products[0]) is not None:  # sink term
                index_of_sub_in_meta = actual_meta_lookup[substrates[0]]
                current_col[index_of_sub_in_meta] = negative_one * metabolite_array[index_of_sub_in_meta]
            else:  # not a source or sink
                index_of_prod_in_meta = actual_meta_lookup[products[0]]
                index_of_sub_in_meta = actual_meta_lookup[substrates[0]]

                current_col[index_of_sub_in_meta] = negative_one * metabolite_array[index_of_sub_in_meta]
                current_col[index_of_prod_in_meta] = metabolite_array[index_of_sub_in_meta]

        S_matrix.append(current_col)

    return S_matrix


def print_matrix(np_mat):
    mat_df = pd.DataFrame([[item.get() for item in ele] for ele in np_mat]).transpose()
    print(mat_df)

if __name__ == '__main__':
    file_name = 'simple_net.xlsx'
    net_df = pd.read_excel(file_name)

    uni_meta = extract_metabolites(net_df)
    meta_lookup = create_meta_lookup(uni_meta)

    source_sink_lookup, actual_meta = find_source_sink_index(meta_lookup)

    actual_meta_lookup = create_meta_lookup(actual_meta)

    s_mat = build_S_matrix(net_df, uni_meta, actual_meta_lookup, source_sink_lookup)
    s_mat_np = np.array(s_mat)

    print_matrix(s_mat_np)

    print(actual_meta_lookup)