from py2neo import neo4j, cypher
import core
import collections


def match_sample(query):
    
    from config import cypher_session
    
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/')
    
    query_list = []
    
    sample_query = neo4j.CypherQuery(graph_db,'MATCH (n:@SUBJECT@) WHERE n.name =~ "'+query+'.*' +'"UNWIND [n.name] + n.alias AS name_alias WITH n, name_alias RETURN ID(n), name_alias')
    
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
    
    hit_symb_query = 'MATCH(n:EntrezID)-[r:REFFERED_TO]-(m) WHERE n.name IN {hit_genes} RETURN n.name, m.name'
    
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

def get_gene_names(res_list, nodes, request):
    
    for i in core.BasicResultsIterable(res_list):
        if len(i) > 0:
            nodes.add(core.GeneNode(i))

def get_subject (res_list, nodes, request):
    for i in core.BasicResultsIterable(res_list):
        if len(i) > 0:
            nodes.add(core.SubjectNode(i))

def get_shortest_paths (request, request_post):
    
    from config import cypher_session
    session = cypher.Session(cypher_session)
    tx = session.create_transaction()
    
    if request_post.has_key('query_samples'):
        
        ##if bypassing table I think this would be necessary...
        #simply create a graph with the node(s) requested by the user 
        
        final_nodes_list = core.get_nodes(list(set(request.session['query_samples']['SampleID'].values())), 'SampleID', request)
        
        node_names = final_nodes_list.display_names()
        
        if len(node_names) == 1:
            title = 'Sample: '+ node_names[0]
        else:
            title = 'Samples: ' + string.joinfields(node_names, ',')
        
        return final_nodes_list.tolist(), [], title
    else:
        
        cur_node_set = set(request_post['var_select'] + request_post['seed_select'])
        
        cur_edge_set = core.RelationshipSet()
        
        samp_set = set()
        
        for i in request.session['query_samples']['SampleID'].values():
            if i != None:
                samp_set.add(i)
        
        final_nodes_list = core.get_nodes(list(samp_set), 'Sample', request)
        
        sp_query = string.joinfields(['MATCH (n:EntrezID{name:{var_select}})-[:MAPPED_TO]->(np) WITH np MATCH (m:EntrezID{name:{seed_select}})-[:MAPPED_TO]->(mp) WITH mp, np MATCH p=(np)-[:ASSOC*1..2]->(mp)',
                    'WHERE ALL(x IN NODES(p) WHERE (x)<-[:MAPPED_TO]-()) AND ALL(x IN RELATIONSHIPS(p) WHERE HAS(x.score) AND x.score > {score}) WITH p, REDUCE(val=0, x in RELATIONSHIPS(p)| val+x.score) as use_score,',
                    'LENGTH(RELATIONSHIPS(p)) AS min_len ORDER BY min_len, use_score DESC RETURN p limit 1'], " ")
        
        for var, seed in itertools.product(request_post['var_select'], request_post['seed_select']):
            tx.append(sp_query, {'var_select':var, 'seed_select':seed, 'score':int(request.session['string_conf']*1000)})
            
        all_path = tx.execute()
        
        print all_path
        
        nodes_to_recon = []
        
        prot_to_gene = string.joinfields([
                'MATCH (n:StringID{name:{cur_prot}})-[:MAPPED_TO]-()-[:EXTERNAL_ID]-(m) RETURN n.name,m.name'
            ], " ")
        
        for i in core.BasicResultsIterable(all_path):
            if len(i) > 0:
                for j in i[0].nodes:
                    tx.append(prot_to_gene, {"cur_prot":j["name"]})
        
        recon_prot_genes = tx.execute()
        
        prot_to_genes = collections.defaultdict(list)
        
        #the isinstance , tuple code below is to address a too many values to unpack error for 12-00192
        
        for i in core.BasicResultsIterable(recon_prot_genes):
            if len(i) > 0:
                if isinstance(i[0], tuple):
                    for j in i:
                        prot_to_genes[j[0]].append(j[1])
                else:
                    prot_to_genes[i[0]].append(i[1])
        
        for i in core.BasicResultsIterable(all_path):
            if len(i) > 0:
                for j in i[0].relationships:
                    if cur_edge_set.check(j.start_node["name"], j.end_node["name"], undirected=True, map_dict=prot_to_genes) == False:
                        cur_edge_set.add(j, prot_to_genes)
        
        for i in cur_edge_set.nodes():
            cur_node_set.add(i)
        
        #then figure out whether any of these nodes are connected to one-another
        
        dp_query = 'MATCH (n:Gene{name:{gene_1}})-[:EXTERNAL_ID]-()-[:MAPPED_TO]-()-[r:ASSOC]->()-[:MAPPED_TO]-()-[:EXTERNAL_ID]-(m:Gene{name:{gene_2}}) WHERE HAS(r.score) AND r.score > {score} WITH n,m, MAX(r.score) AS max_score RETURN n.name,m.name, max_score'
        # 
        for i,j in itertools.combinations(list(cur_node_set), 2):
            tx.append(dp_query, {'gene_1':i, 'gene_2':j, 'score':int(request.session['string_conf']*1000)})
            
        rem_rels = tx.commit()
        
        for i in  core.BasicResultsIterable(rem_rels):
            if len(i) > 0 and cur_edge_set.check(i[0], i[1], undirected=True, map_dict=None) == False:
                cur_edge_set.direct_add(i)
        
        #so now prep the nodes and edges for the view
        
        sample_dict = {}
        
        for i in request.session['query_samples']['SampleID'].items():
            if i[1] != None:
                sample_dict[i[0]] = [i[1]]
            else:
                sample_dict[i[0]] = []
        
        #print cur_node_set
        
        gene_nodes_list = core.get_nodes(list(cur_node_set), 'Gene', request, param_list=[sample_dict])
        
        #print sample_dict
        
        samp_node_ref = map(lambda x: x.id, final_nodes_list)
        
        final_edge_list = []
        
        for i_ind, i in enumerate(gene_nodes_list):
            for j in i.children():
                child_dict = j.todict()
                for k_ind, k in enumerate(child_dict["attributes"]["other_nodes"]):
                    final_edge_list.append({'source':samp_node_ref.index(k), 'target':i_ind+len(samp_node_ref),
                                            'attributes':get_link_atts(request, final_nodes_list.getNode(k), i,child_dict["attributes"]["node_type"], {},child_dict["attributes"]["meta"]["is_hit"][k_ind])})
            
            final_nodes_list.add(i)
        
        node_ref = map(lambda x: x.id, final_nodes_list)
        
        for f,t,props in cur_edge_set:
            final_edge_list.append({'source':node_ref.index(f), 'target':node_ref.index(t), 'attributes':get_link_atts(request, final_nodes_list.getNode(f), final_nodes_list.getNode(t), None, props)})
        
        return final_nodes_list.tolist(), final_edge_list, 'Initial HitWalker Result'
