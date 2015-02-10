from py2neo import neo4j

def match_sample(query):
    
    from config import cypher_session
    
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/')
    
    query_list = []
    
    sample_query = neo4j.CypherQuery(graph_db,'MATCH (n:CellLine) WHERE n.name =~ "'+query+'.*' +'"UNWIND [n.name] + n.alias AS name_alias WITH n, name_alias RETURN ID(n), name_alias')
    
    for i in sample_query.execute().data:
        query_list.append({'id':i.values[0], 'text':i.values[1], 'search_list':[i.values[1]]})
    
    return query, query_list



def netprop_rwr(request, seed_gene_nl, initial_graph_file, string_session_name, res_prob_session_name, conv_thresh_session_name, max_iter_session_name):
    
    from rwr_utils import threshold_graph, compute_net_prop_mat, compute_rwr
    
    #add in imports from the rwr files
    
    if (request.session.has_key(string_session_name) == False):
        raise Exception(string_session_name+" is necessary but not present in the user session")
    
    if (request.session.has_key(res_prob_session_name) == False):
        raise Exception(res_prob_session_name+" is necessary but not present in the user session")
    
    if (request.session.has_key(conv_thresh_session_name) == False):
        raise Exception(conv_thresh_session_name+" is necessary but not present in the user session")
    
    if (request.session.has_key(max_iter_session_name) == False):
        raise Exception(max_iter_session_name+" is necessary but not present in the user session")
    
    use_graph, dim_dict, use_graph_file, use_graph_names = threshold_graph(initial_graph_file, float(request.session[string_session_name]))
        
    use_graph_np = compute_net_prop_mat(use_graph)
    
    seed_dict = seed_gene_nl.todict()
    
    res_vec = compute_rwr(use_graph_np, dim_dict, float(request.session[res_prob_session_name]), float(request.session[conv_thresh_session_name]), int(request.session[max_iter_session_name]), seed_dict)
    
    prot_nl = core.NodeList()
    
    for i in dim_dict.keys():
        prot_score = [i, res_vec[dim_dict[i]]]
        temp_node = core.BasicNode(prot_score, only_child=True)
        prot_nl.add(temp_node)
    
    return core.SeedList(prot_nl)