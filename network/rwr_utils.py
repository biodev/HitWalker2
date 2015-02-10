#!/usr/bin/env python
from scipy import io, sparse
from scipy.sparse import csgraph, diags
from numpy import array, asarray, linalg, sqrt
from os import path
from string import joinfields
import itertools


def threshold_graph (initial_graph_file, edge_thresh):
    try:
        
        #first test if the desired subsetted graph already exists
        outfile_list = list(path.splitext(initial_graph_file))
        outfile_base = outfile_list.pop(0)
        outfile_list.insert(0, str(int(edge_thresh*1000)))
        outfile_suf = joinfields(outfile_list, '')
        outfile_mat = outfile_base + path.extsep + outfile_suf
        outfile_names = outfile_base + path.extsep + outfile_suf + path.extsep + "names"
        outfile_names = outfile_names.replace(path.extsep + "mtx", "")
        
        if path.isfile(outfile_mat) and path.isfile(outfile_names):
            print "Keeping old graph"
            ret_graph = io.mmread(outfile_mat)
            ret_inp = open(outfile_names, 'r')
            
            cur_line = 0
            ret_dict = {}
            
            for i in ret_inp:
                cur_i = i.strip()
                if ret_dict.has_key(cur_i):
                    raise Exception("Unexpected duplicate key found")
                else:
                    ret_dict[cur_i] = cur_line
                cur_line += 1
            
            ret_inp.close()
            
            return ret_graph, ret_dict, outfile_mat, outfile_names
        
        else:
            
            print "Making new graph"
            
            print "Reading"
            
            #otherwise create it
            initial_graph = io.mmread(initial_graph_file)
            dimnames = open(initial_graph_file.replace("mtx", "") + 'names', 'r')
            
            if initial_graph.__class__.__name__ != "coo_matrix":
                raise Exception("Matrix not a coo_matrix")
            
            dim_list = []
            
            for i in dimnames:
                dim_list.append(i)
            
            dimnames.close()
            
            print "Subsetting"
            
            keep_nodes = set()
           
            for i in xrange(len(initial_graph.data)):
                if initial_graph.data[i] > edge_thresh:
                    keep_nodes.add(initial_graph.row[i])
                else:
                    initial_graph.data[i] = 0
            
            use_graph = initial_graph.tocsr()
            use_graph.eliminate_zeros()
            
            keep_node_list = list(keep_nodes)
            
            keep_graph = use_graph[:, keep_node_list][keep_node_list,:]
            
            #reorder the row/column names
            
            ret_list=list()
            ret_dict = {}
            cur_line = 0
            
            for i in keep_node_list:
                ret_list.append(dim_list[i])
                cur_i = dim_list[i].strip()
                if ret_dict.has_key(cur_i):
                    raise Exception("Unexpected duplicate key found")
                else:
                    ret_dict[cur_i] = cur_line
                cur_line += 1
            
            #write the graph
            
            io.mmwrite(outfile_mat, keep_graph)
            
            out_p = open (outfile_names, "w")
            out_p.writelines(ret_list)
            out_p.close()
            
            return keep_graph, ret_dict, outfile_mat, outfile_names   
        
    except Exception as e:
        print e

def compute_net_prop_mat(use_graph):
    
    degree = 1/sqrt(use_graph.sum(0))
    diag_degree = diags(asarray(degree), array([0]), format="csr")
    
    return diag_degree.dot(use_graph).dot(diag_degree)

#this is the last function to finish/check
def compute_rwr(use_graph, dim_dict, c_val, conv_thresh, max_iter, seed_dict):
    
    try:
        residue = 1
        iter_num = 1
        
        #make sure all elements of seed dict are in dim_dict
        
        diff_prots = set(seed_dict).difference(set(dim_dict.keys()))
        
        for i in diff_prots:
            print "popping " + i
            seed_dict.pop(i)
        
        #then create a nx1 matrix and convert to csr format as well
        
        col_list = []
        row_list = []
        data_list = []
        
        for i in seed_dict.keys():
            data_list.append(seed_dict[i])
            row_list.append(dim_dict[i])
            col_list.append(0)
        
        data_array = array(data_list)
        
        if sum(data_array) > 1:
            data_array = data_array/float(sum(data_array))
        
        prox_vector = sparse.coo_matrix((data_array, (array(row_list), array(col_list))), shape=(len(dim_dict), 1)).tocsr()
        
        restart_vector = prox_vector.copy()
        
        while(residue > conv_thresh and iter_num < max_iter):
            
            old_prox_vector = prox_vector.copy()
            prox_vector = (1-c_val) * use_graph.dot(prox_vector) + (c_val * restart_vector)
            
            residue = linalg.norm(abs(prox_vector.todense() - old_prox_vector.todense()))
            
            iter_num += 1
        
        if iter_num == max_iter:
            print "WARNING: Did not converge in " + str(max_iter) + " iterations"
        
        return map(lambda x:x[0], prox_vector.todense().tolist())
        
    except Exception as e:
        print e
    

if __name__ == '__main__':
    
    #from rwr_utils import threshold_graph, compute_net_prop_mat, compute_rwr
   use_graph, dim_dict, use_graph_file, use_graph_names = threshold_graph("protein.links.detailed.v9.05.9606.mm.mtx", .4)
   
   use_graph_np = compute_net_prop_mat(use_graph)
   
   seed_dict = {"9606.ENSP00000292535":1, "9606.ENSP00000368349":1, "9606.ENSP00000313021":1, "9606.ENSP00000384479":1, "9606.ENSP00000401504":1}
   
   res_vec = compute_rwr(use_graph_np, dim_dict, .3, 1e-10, 100, seed_dict)
   
   io.mmwrite("test_rwr_results.mm", res_vec)
   
   graph_dict = {}
   
   new_dim_dict = {}
   
   for i in dim_dict.keys():
    if new_dim_dict.has_key(str(dim_dict[i])) == False:
        new_dim_dict[str(dim_dict[i])] = i
    else:
        raise Exception("Unexpected Duplicate")
    
   for data, row, col in itertools.izip(use_graph.data, use_graph.row, use_graph.col):
        if graph_dict.has_key(new_dim_dict[str(row)]):
            if graph_dict[new_dim_dict[str(row)]].has_key(new_dim_dict[str(col)]):
                raise Exception("Unexpected Duplicate")
            else:
                graph_dict[new_dim_dict[str(row)]][new_dim_dict[str(col)]] = data
        else:
            graph_dict[new_dim_dict[str(row)]] = {new_dim_dict[str(col)]:data}
   
    