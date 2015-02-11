from py2neo import neo4j, cypher
import core
import collections

def match_sample(query):
    
    from config import cypher_session
    
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/')
    
    query_list = []
    
    sample_query = neo4j.CypherQuery(graph_db,'MATCH (n:CellLine) WHERE n.name =~ "'+query+'.*' +'"UNWIND [n.name] + n.alias AS name_alias WITH n, name_alias RETURN ID(n), name_alias')
    
    for i in sample_query.execute().data:
        query_list.append({'id':i.values[0], 'text':i.values[1], 'search_list':[i.values[1]]})
    
    return query, query_list


def gene_seed_list_to_protein(request, seed_list):
    
    gene_to_prot_query = 'MATCH (gene:EntrezID{name:{gene_name}})-[:MAPPED_TO]->(protein) RETURN gene.name, protein.name'
    
    from config import cypher_session
    
    session = cypher.Session(cypher_session)
    tx = session.create_transaction()
    
    if isinstance(seed_list, core.SeedList):
    
        for i in seed_list.nodeList().ids():
            tx.append(gene_to_prot_query, {'gene_name':i})
    
        gene_to_prot_list = tx.commit()
        
        prot_dict = collections.defaultdict(float)
        
        for i in core.BasicResultsIterable(gene_to_prot_list):
            if len(i) > 0:
                if isinstance(i[0], tuple):
                    for j in i:
                        prot_dict[j[1]] = max(seed_list.getScores([j[0]])[0], prot_dict[j[1]])
                else:
                    prot_dict[i[1]] = max(seed_list.getScores([i[0]])[0], prot_dict[i[1]])
        
    
        #first make a BasicNode list and coerce to a SeedNode List
        
        prot_gl = core.NodeList()
        
        for i in prot_dict.items():
            prot_gl.add(core.BasicNode(i, only_child=True))
            
        prot_sl = core.SeedList(prot_gl)
        
        return prot_sl
    
    #used by protein_seed_list_to_gene
    elif isinstance(seed_list, list):
        
        for i in seed_list:
            tx.append(gene_to_prot_query, {'gene_name':i})
            
        gene_to_prot_list = tx.commit()
        
        prot_dict = collections.defaultdict(list)
        
        for i in core.BasicResultsIterable(gene_to_prot_list):
             if len(i) > 0:
                if isinstance(i[0], tuple):
                    for j in i:
                        prot_dict[j[1]].append(j[0])
                else:
                    prot_dict[i[1]].append(i[0])
        
        return prot_dict
        
    else:
        raise Exception("seed_list needs to be either of class SeedList or list")

def protein_seed_list_to_gene(request, seed_list, limit_list):
    
    if limit_list != None:
        prot_genes = gene_seed_list_to_protein(request, limit_list)
        use_prots = prot_genes.keys()
    else:
        use_prots = seed_list.nodeList().ids()
    
    prot_to_gene_query = 'MATCH (prot:StringID{name:{prot_name}})<-[:MAPPED_TO]-(gene) RETURN prot.name, gene.name'
    
    from config import cypher_session
    
    session = cypher.Session(cypher_session)
    tx = session.create_transaction()

    for i in use_prots:
        tx.append(prot_to_gene_query, {'prot_name':i})
    
    prot_to_gene_list = tx.commit()
    
    #make this into a gene SeedList similar to above
    
    gene_dict = collections.defaultdict(float)
    
    for i in core.BasicResultsIterable(prot_to_gene_list):
        if len(i) > 0:
            if isinstance(i[0], tuple):
                for j in i:
                    if seed_list.nodeList().hasNode(j[0]):
                        gene_dict[j[1]] = max(seed_list.getScores([j[0]])[0], gene_dict[j[1]])
            else:
                if seed_list.nodeList().hasNode(i[0]):
                    gene_dict[i[1]] = max(seed_list.getScores([i[0]])[0], gene_dict[i[1]])
    
    
    gene_nl = core.NodeList()
        
    for i in gene_dict.items():
        gene_nl.add(core.BasicNode(i, only_child=True))
        
    gene_sl = core.SeedList(gene_nl)
    
    return gene_sl

def make_seed_list(hit_sl):
    
    from config import cypher_session
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/') 
    
    hit_symb_query = 'MATCH(n:Gene)-[r:KNOWN_AS]-(m) WHERE r.status = "symbol" AND n.name IN {hit_genes} RETURN n.name, m.name'
    
    data=neo4j.CypherQuery(graph_db, hit_symb_query).execute(**{"hit_genes":hit_sl.nodeList().ids()})
    
    seed_list = []
        
    for i in data:
        if hit_sl.nodeList().hasNode(i.values[0]):
            seed_list.append({'gene':i.values[0], 'symbol':i.values[1], 'score':hit_sl.getScores([i.values[0]])[0]})
    
    seed_list.sort(key=lambda x: x["score"], reverse=True)
        
    return seed_list


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