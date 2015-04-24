from py2neo import neo4j, cypher
import core
import collections
import string
import itertools
import json
import copy

def match_sample(query):
    
    from config import cypher_session
    
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/')
    
    query_list = []
    
    sample_query = neo4j.CypherQuery(graph_db,'MATCH (n:@SUBJECT@)-[:DERIVED]->() WITH n, CASE n.alias WHEN null THEN [n.name] ELSE [n.name]+n.alias END AS alias_query UNWIND alias_query AS name_alias WITH n, name_alias WHERE name_alias =~ "'+query+'.*' +'" RETURN ID(n)+name_alias, name_alias, COLLECT(n.name)')
    
    for i in sample_query.execute().data:
        query_list.append({'id':i.values[0], 'text':i.values[1], 'search_list':i.values[2]})
    
    return query, query_list


def match_gene(query):
    
    from config import cypher_session
    
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/')
    query_list = []
    
    sample_query = neo4j.CypherQuery(graph_db,'MATCH (n:Symbol)<-[r:REFFERED_TO]-(m) UNWIND [n.name] + r.synonyms AS name_syns WITH n,m,name_syns WHERE name_syns =~"'+query+'.*"  RETURN m.name+name_syns, name_syns, COLLECT(m.name) ')
    
    search_res = []
    
    for i in sample_query.execute().data:
        
        query_list.append({'id':i.values[0], 'text':i.values[1], 'search_list':i.values[2]})
    
    return query, query_list

def match_pathway(query):
    
    query = json.loads(query)
    
    from config import cypher_session
    
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/')
    
    query_list = [] 
    
    if len(query) == 0:
        return query, []
    else:
        sample_query = neo4j.CypherQuery(graph_db,'MATCH (g:EntrezID)-[:PATHWAY_CONTAINS]-(path) WITH path, COUNT(g) AS g_count WHERE g_count < 200 MATCH (gene:EntrezID)-[r:PATHWAY_CONTAINS]-(path) WHERE ALL(x IN '+str(json.dumps(query))+' WHERE ANY(y in x WHERE gene.name = y)) RETURN DISTINCT ID(path), path.name + " (n=" + g_count + ")"')
        
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
            
def get_pathway(res_list, nodes, request):
    
    gene_names = []
    path_name = []
    
    for i in core.BasicResultsIterable(res_list):
        if len(i) > 0:
            if isinstance(i[0], tuple):
                use_i = i[:]
            else:
                use_i = [i[:]]
            
            for j in use_i:
                path_name.append(j[0])
                gene_names.extend(j[1])
    
    #as pathways are currently only produced from the copy_nodes interface
    import config
    
    gene_nl = core.get_nodes(gene_names, 'Gene', request, config_struct=config.edge_queries['nodes'], missing_param="skip")
    
    #a pathway in this context is a special metanode
    nodes.extend(gene_nl)

def get_link_atts(request, from_node, to_node, rel_type, props, is_hit=False):
    from_node_dict = from_node.todict()
    to_node_dict = to_node.todict()
    
    import config
    
    if from_node_dict['attributes'].has_key('node_type') and to_node_dict['attributes'].has_key('node_type'):
        
        if from_node_dict["attributes"]["node_type"] == "Subject" and to_node_dict["attributes"]["node_type"] == "Gene":
            
            is_hit = rel_type.endswith("_Hit")
            temp_rel_name = rel_type.replace("_Hit", "")
            
            if temp_rel_name == config.data_types['target']:
                
                if request.session.has_key('hitwalker_score') and request.session['hitwalker_score'].nodeList().hasNode(to_node_dict["id"]) == True and from_node_dict['id'] == request.session['query_samples']['SampleID'][config.data_types['target']]:
                    temp_rel_name = "Ranked_" + temp_rel_name
                else:
                    temp_rel_name = "Observed_" + temp_rel_name
            else:
            
                if is_hit == True:
                    temp_rel_name = "Observed_" + temp_rel_name
                else:
                    temp_rel_name = "Possible_" + temp_rel_name
            
            return {'type':temp_rel_name}
        
        elif from_node_dict["attributes"]["node_type"] == "Gene" and to_node_dict["attributes"]["node_type"] == "Gene":
            
            return {'type':'STRING'}
        
        else:
            return {'type':'Unknown'}#'rank':1000000}
    else:
        return {'type':'Unknown'}#'rank':1000000}


def no_rels(query, subj, request, config_struct_nodes, cur_graph):
    return cur_graph

def gene_to_sample (genes, sampleids, request, config_struct_nodes, cur_graph):
    
    #as the genes have a dependency on labids, retrieve the nodes again as well as the children (a little wasteful but direct for now...)
   
    #where labids should be a list of labIDs that will get directly substituted into the query
    param_list = [{'SAMPLE':list(sampleids)}]
    
    gene_nodes = core.get_nodes(list(genes), 'Gene', request, config_struct=config_struct_nodes, param_list=param_list)
    
    #for each of the genes, examine their children and record an edge corresponding the the child type and labid
    #use this list and the specified labids to determine if edges are appropriate
    
    for i in gene_nodes:
        for j in i.children():
            #add children to existing nodes
            cur_graph["nodes"].addChild(i.id, j)
            #add links
            child_dict = j.todict()
            #print j.id
            for k_ind, k in enumerate(child_dict["attributes"]["other_nodes"]):
                #print k
                cur_graph['links'].append({'source':cur_graph["nodes"].nodeIndex(k), 'target':cur_graph["nodes"].nodeIndex(i.id), 'attributes':get_link_atts(request, cur_graph["nodes"].getNode(k) , cur_graph["nodes"].getNode(i.id) , child_dict["attributes"]["node_type"], {}, child_dict["attributes"]["meta"]["is_hit"][k_ind])})
    return cur_graph

#as gene_to_gene is used for determination of pathways add in the path_conf session variable instead of the standard string_conf
def gene_to_gene(query_genes, subj_genes, request, config_struct_nodes, cur_graph):
    
    import config
    
    session = cypher.Session(config.cypher_session)
    tx = session.create_transaction()
    
    print len(query_genes), len(subj_genes)
    
    cur_edge_set = core.RelationshipSet()
    
    for i in query_genes:
        tx.append(config.gene_rels['query'], {"FROM_GENE":i, "TO_GENES":subj_genes, "string_conf":core.iterate_dict(request.session, ['path_conf'])})
    
    for i in core.BasicResultsIterable(tx.execute()):
        if len(i) > 0:
            if isinstance(i[0], tuple):
                use_i = i[:]
            else:
                use_i = [i[:]]
            for j in use_i:
                if len(j) > 0 and cur_edge_set.check(j[0], j[1]) == False:
                    cur_edge_set.direct_add(j)
                    cur_graph['links'].append({'source':cur_graph["nodes"].nodeIndex(j[0]), 'target':cur_graph["nodes"].nodeIndex(j[1]), 'attributes':get_link_atts(request, cur_graph["nodes"].getNode(j[0]), cur_graph["nodes"].getNode(j[1]), None, {'score':j[2]})})
    
    return cur_graph

def get_pathways_sample (request, request_post):
    
    import config
    print request_post
    subj_nodes = map(lambda x:{'node_type':'Subject', 'id':x}, request_post['sample_name'])
    query_nodes = map(lambda x:{'node_type':'Pathway', 'id':x}, request_post['pathway_name'])
    
    new_edge_queries = copy.deepcopy(config.edge_queries)
    new_edge_queries['edges']['Gene']['Gene']['handler'] = gene_to_gene
    
    cur_graph = core.copy_nodes(subj_nodes, query_nodes, request, new_edge_queries, never_group=True)
    
    return cur_graph['nodes'].tolist(), cur_graph['links'], query_nodes[0]['id']

def get_shortest_paths (request, request_post):
    
    from config import cypher_session
    session = cypher.Session(cypher_session)
    tx = session.create_transaction()
    
    if request_post.has_key('query_samples'):
        
        ##if bypassing table I think this would be necessary...
        #simply create a graph with the node(s) requested by the user 
        
        final_nodes_list = core.get_nodes(list(set(request.session['query_samples']['SampleID'].values())), 'Sample', request)
        
        node_names = final_nodes_list.display_names()
        
        if len(node_names) == 1:
            title = 'Subject: '+ node_names[0]
        else:
            title = 'Subject: ' + string.joinfields(node_names, ',')
        
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
                    'WHERE ALL(x IN NODES(p) WHERE (x)<-[:MAPPED_TO]->()) AND ALL(x IN RELATIONSHIPS(p) WHERE HAS(x.score) AND x.score > {score}) WITH p, REDUCE(val=0, x in RELATIONSHIPS(p)| val+x.score) as use_score,',
                    'LENGTH(RELATIONSHIPS(p)) AS min_len ORDER BY min_len, use_score DESC RETURN p limit 1'], " ")
        
        for var, seed in itertools.product(request_post['var_select'], request_post['seed_select']):
            print var, seed
            tx.append(sp_query, {'var_select':var, 'seed_select':seed, 'score':int(request.session['string_conf']*1000)})
            
        all_path = tx.execute()
        
        print all_path
        
        nodes_to_recon = []
        
        prot_to_gene = string.joinfields([
                'MATCH (n:StringID{name:{cur_prot}})-[:MAPPED_TO]-(m) RETURN n.name,m.name'
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
        
        dp_query = 'MATCH (n:EntrezID{name:{gene_1}})-[:MAPPED_TO]-()-[r:ASSOC]->()-[:MAPPED_TO]-(m:EntrezID{name:{gene_2}}) WHERE HAS(r.score) AND r.score > {score} WITH n,m, MAX(r.score) AS max_score RETURN n.name,m.name, max_score'
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
        
        print cur_node_set
        
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
