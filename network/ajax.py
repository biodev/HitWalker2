from django.utils import simplejson
from dajaxice.decorators import dajaxice_register
from network.models import user_parameters
from py2neo import neo4j , cypher
import json
import string
import itertools
import sys
import views


#need to use directed edges only as currently it is hardcoded to expect say, n, to be the first element in some of the functions below... 

def get_avail_paths_helper (name, index, indexed_name, direction):
    graph_db = neo4j.GraphDatabaseService()
    
    if direction == "outgoing":
        query_str = 'START n=node:'+index+'('+indexed_name+'="' + name + '") MATCH (n)-[r]->() RETURN DISTINCT TYPE(r)'
    elif direction == "incoming":
        query_str = 'START n=node:'+index+'('+indexed_name+'="' + name + '") MATCH (n)<-[r]-() RETURN DISTINCT TYPE(r)'
    else:
        print "ERROR: Unknown direction type"
        sys.exit()
    data, cypher_meta = cypher.execute(graph_db,query_str)
    
    ret_json = []
    #get rid of the array of arrays that is returned
    for i in data:
        ret_json.append(i[0])
    
    return ret_json

@dajaxice_register
def get_available_paths(request, name, index, indexed_name):
    
    overall_dict = set()
    ret_dict = dict()
    
    for i in ["outgoing", "incoming"]:
        ret_dict[i] = get_avail_paths_helper (name, index, indexed_name, i)
        
    #if there are duplicates, chose one (arbitrarily) for simplicity (path should be explorable either direction
    
    common_elems = set.intersection(set(ret_dict["outgoing"]), set(ret_dict["incoming"]))
    
    #remove, here assuming the elements are already distinct as was returned in the neo query
    if len(common_elems) > 0:
        list_elems = list(common_elems)
        for i in list_elems:
            ret_dict["incoming"].remove(i)
    
    return json.dumps(ret_dict)


#the problem with this function is that for highly connected nodes it takes too long to iterate over data
#can probably handle the sanity check elsewhere
def get_available_paths_dep(request, name, index, indexed_name):
    
    graph_db = neo4j.GraphDatabaseService()
    
    #query_str = 'START n=node:'+index+'('+indexed_name+'="' + name + '") MATCH (n)-[r]-() RETURN DISTINCT TYPE(r)'
    query_str = 'START n=node:'+index+'('+indexed_name+'="' + name + '") MATCH (n)-[r]-() RETURN r'
    
    data, cypher_meta = cypher.execute(graph_db,query_str)
    
    ret_list = []
    
    rel_type = dict()
    
    for i in data:
        if rel_type.has_key(i[0].type) == False:
            rel_type[i[0].type] = i[0].get_properties().keys()
        else:
            for j in rel_type[i[0].type]:
                assert i[0].get_properties().has_key(j)
    
    print rel_type
    
    ret_json  = {'num_edges':len(rel_type.keys()), 'edges':rel_type.keys()}
    
    return json.dumps(ret_json)

def get_all_permutations(iterable):
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.permutations(s, r) for r in range(len(s)+1))

#neotool cypher 'START n=node(22133), p=node(239869,26943,36621,23430,31954) MATCH (n)-[r:TRANSCRIBED]-(m)-[r2?]-(p) RETURN m,r2,p'
#neotool cypher 'START n=node(22133),p=node(239869,26943,36621,23430,31954) MATCH (n)-[r:TRANSCRIBED]-(m) WITH m skip 1 limit 1  MATCH (m)-[r2?]-(p) RETURN m,r2,p'
#neotool cypher 'START n=node(22133) MATCH (n)-[r:TRANSCRIBED]-(m) WITH distinct m skip 2 limit 1 START p=node(239869,26943,36621,23430,31954)  MATCH (m)-[r2?]-(p) RETURN m,r2,p'


def return_basic_query_dep (path_name, cur_node_id, other_node_ids, path_where, graph_db, request,keep_first, path_count):
    #return 'START n=node('+ str(cur_node_id) + '), p=node(' + string.joinfields(map(str, other_node_ids), ",") + ') MATCH (n)-[r:' + path_name + ']-(m)-[r2?]-(p) ' + path_where + ' RETURN m,r2,p'
    
    if keep_first == True:
        if path_where != '':
            path_where = 'AND (' + path_where + ')'
        return 'START n=node({cur_node_id}) MATCH (n)-[r:' + path_name + ']-(m) WHERE (NOT(ID(m) IN {other_node_ids})) ' + path_where + ' WITH DISTINCT m LIMIT 1 START p=node({other_node_ids}) MATCH (m)-[r2?]-(p) RETURN m,r2,p'
    else:
        if path_where != '':
            path_where = 'WHERE ' + path_where
        return 'START n=node({cur_node_id}), p=node({other_node_ids}) MATCH (n)-[r:' + path_name + ']-(m)-[r2?]-(p) ' + path_where + ' RETURN m,r2,p'
    
def return_basic_query (path_name, cur_node_id, other_node_ids, graph_db, request,keep_first, path_count):
    #return 'START n=node('+ str(cur_node_id) + '), p=node(' + string.joinfields(map(str, other_node_ids), ",") + ') MATCH (n)-[r:' + path_name + ']-(m)-[r2?]-(p) ' + path_where + ' RETURN m,r2,p'
    
    if keep_first == True:
        return 'START n=node({cur_node_id}) MATCH (n)-[r:' + path_name + ']-(m) WHERE (NOT(ID(m) IN {other_node_ids})) WITH DISTINCT m LIMIT 1 START p=node({other_node_ids}) MATCH (m)-[r2?]-(p) RETURN m,r2,p'
    else:
        return 'START n=node({cur_node_id}), p=node({other_node_ids}) MATCH (n)-[r:' + path_name + ']-(m)-[r2?]-(p) RETURN m,r2,p'

def retrieve_m_graph_query_dep(path_name, cur_node_id, other_node_ids, path_where, graph_db, request,keep_first, path_count):
    
    edge_thresh = request.session['edge_thresh']*1000
    
    if keep_first == True:
        if path_where != '':
            path_where = 'AND (' + path_where + ')'
        return 'START n=node({cur_node_id}) MATCH (n)-[:EXTERNAL_ID]-()-[:MAPPED_TO]-()-[r:ASSOC]->()-[:MAPPED_TO]-()-[:EXTERNAL_ID]-(m) WHERE r.score > \
                        ' + str(edge_thresh) +'AND NOT(ID(m) IN {other_node_ids}) ' + path_where + ' WITH DISTINCT m LIMIT 1 START p=node({other_node_ids}) MATCH (m)-[r2?]-(p)  RETURN m,r2,p'
    else:
        if path_where != '':
            path_where = 'WHERE ' + path_where
        return 'START n=node({cur_node_id}) MATCH (n)-[:EXTERNAL_ID]-()-[:MAPPED_TO]-()-[r:ASSOC]->()-[:MAPPED_TO]-()-[:EXTERNAL_ID]-(m) WHERE r.score > \
                        ' + str(edge_thresh) +' WITH DISTINCT m START p=node({other_node_ids}) MATCH (m)-[r2?]-(p) ' + path_where + ' RETURN m,r2,p'
                    
def retrieve_m_graph_query(path_name, cur_node_id, other_node_ids, graph_db, request,keep_first, path_count):
    
    edge_thresh = float(request.session['string_conf'])*1000
    
    if keep_first == True:
        return 'START n=node({cur_node_id}) MATCH (n)-[:EXTERNAL_ID]-()-[:MAPPED_TO]-()-[r:ASSOC]->()-[:MAPPED_TO]-()-[:EXTERNAL_ID]-(m) WHERE r.score > \
                        ' + str(edge_thresh) +'AND NOT(ID(m) IN {other_node_ids}) WITH DISTINCT m LIMIT 1 START p=node({other_node_ids}) MATCH (m)-[r2?]-(p)  RETURN m,r2,p'
    else:
        return 'START n=node({cur_node_id}) MATCH (n)-[:EXTERNAL_ID]-()-[:MAPPED_TO]-()-[r:ASSOC]->()-[:MAPPED_TO]-()-[:EXTERNAL_ID]-(m) WHERE r.score > \
                        ' + str(edge_thresh) +' WITH DISTINCT m START p=node({other_node_ids}) MATCH (m)-[r2?]-(p) RETURN m,r2,p'

def return_variant_query (path_name, cur_node_id, other_node_ids, graph_db, request,keep_first, path_count):
    
    #add the relevant variant filters
    
    where_template = request.session['where_template'] 
    necessary_vars = request.session['necessary_vars']
    
    base_query = views.add_where_input_query('Sample', 'START sample=node({cur_node_id}) MATCH (sample)-[r:DNA_DIFFERENCE]-(var)-[u:UNIT_DNA_DIFFERENCE]-(exp) WHERE (exp)<-[:GENOTYPED_USING]-(sample) WITH sample,u,r,var MATCH (var)-[r2:IMPACTS]->(o) RETURN', where_template, necessary_vars)
    
    if keep_first == True:
        return base_query.replace("RETURN", "") + "AND NOT(ID(var) IN {other_node_ids}) WITH DISTINCT var LIMIT 1 START p=node({other_node_ids}) MATCH (var)-[r2?]-(p) RETURN var,r2,p"
    else:
        return base_query.replace("RETURN", "") + 'WITH DISTINCT var START p=node({other_node_ids}) MATCH (var)-[r2?]-(p) RETURN var,r2,p'

def make_query_from_path_dep(path_name, cur_node_id, other_node_ids, path_where, graph_db, request, keep_first, path_count):
    path_funcs = {'IMPACTS_GENE':{'func':return_basic_query},
                  'IMPACTS':{'func':return_basic_query},
                  'TRANSCRIBED':{'func':return_basic_query},
                  'DNA_DIFFERENCE':{'func':return_basic_query},
                  'DNA_DIFFERENCE_GENE':{'func':return_basic_query},
                  'EXTERNAL_ID':{'func':retrieve_m_graph_query},
                  'GENOTYPED_USING':{'func':return_basic_query},
                  'UNIT_DNA_DIFFERENCE':{'func':return_basic_query}}
    if path_funcs.has_key(path_name):
        return path_funcs[path_name]['func'](path_name, cur_node_id, other_node_ids, path_where, graph_db, request, keep_first, path_count)
    else:
        raise Exception("Unknown path " + path_name)
    
def make_query_from_path(path_name, cur_node_id, other_node_ids, graph_db, request, keep_first, path_count):
    path_funcs = {'IMPACTS_GENE':{'func':return_basic_query},
                  'IMPACTS':{'func':return_basic_query},
                  'TRANSCRIBED':{'func':return_basic_query},
                  'DNA_DIFFERENCE':{'func':return_variant_query},
                  'DNA_DIFFERENCE_GENE':{'func':return_basic_query},
                  'EXTERNAL_ID':{'func':retrieve_m_graph_query},
                  'GENOTYPED_USING':{'func':return_basic_query},
                  'UNIT_DNA_DIFFERENCE':{'func':return_basic_query}}
    if path_funcs.has_key(path_name):
        return path_funcs[path_name]['func'](path_name, cur_node_id, other_node_ids, graph_db, request, keep_first, path_count)
    else:
        raise Exception("Unknown path " + path_name)

def make_group_query(group_dict, other_node_ids, query_type, path_name):
    
    print group_dict
    print path_name
    
    if query_type == "common" and len(group_dict.keys()) == 1:
        
        if group_dict.has_key("1") == False:
            raise Exception("group_dict expected to have a 1 key")
        
        query = 'START n=node({group_nodes}) MATCH (n)-[r:' + path_name + ']-(m) WITH  m,COUNT (*) as m_count WHERE m_count = {group_size} WITH DISTINCT m as dist_m limit 10 START p=node({other_node_ids}) MATCH (dist_m)-[r2?]-(p) RETURN dist_m,r2,p'
        params = {'group_nodes':group_dict['1'], 'other_node_ids':other_node_ids, 'group_size':len(group_dict['1'])}
        
    elif query_type == "group_1" and len(group_dict.keys()) == 2:
        if group_dict.has_key("1") == False or group_dict.has_key("2") == False:
            raise Exception("group_dict expected to have a 1 and 2 key")
        
        query = 'START n=node({all_group_nodes}) MATCH (n)-[:' + path_name + ']-(m) WITH COLLECT(n) AS n_col, m WHERE ALL(x IN n_col WHERE ID(x) IN {in_group_nodes} ) \
                WITH n_col, m,COUNT(n_col) AS n_count WHERE n_count = {in_group_count} WITH m AS dist_m limit 10 START p=node({other_node_ids}) MATCH (dist_m)-[r2?]-(p) RETURN dist_m,r2,p'
    
        params = {'all_group_nodes':group_dict['1'] + group_dict['2'],'in_group_nodes':group_dict['1'], 'other_node_ids':other_node_ids, 'in_group_count':len(group_dict['1'])}
        
        print params
    
    else:
        raise Exception("Unimplemented query requested")
    
    return query, params 

#test queries:

{'all_group_nodes': [617038, 617341, 617128, 617113], 'other_node_ids': [76969, 74010, 79334, 71909, 69775], 'in_group_nodes': [617038, 617341], 'in_group_count': 2}

#neotool cypher 'START n=node(617038, 617341, 617128, 617113) MATCH (n)-[:DNA_DIFFERENCE_GENE]-(m) WITH COLLECT(n) AS n_col, m WHERE ALL(x IN n_col WHERE ID(x) IN [617038,617341] ) WITH n_col, m,COUNT(n_col) AS n_count  return n_col,m, n_count limit 10'

#neotool cypher 'START n=node(239875,239872,239871,239870) MATCH p1=(n)-[:DNA_DIFFERENCE_GENE]-(m) WITH COLLECT(n) AS n_col, m WHERE ALL(x IN n_col WHERE ID(x) IN [239875,239872] ) return  count( m) '

#sanity check of result
#neotool cypher 'START n=node(33593,60751,29152,25626,21562,26003), m=node(239871,239870) MATCH (n)-[:DNA_DIFFERENCE_GENE]-(m) return n,m'
    #looks pretty good : )
    
#try a rewrite of the all query...
#old
#neotool cypher 'START n=node(239875,239872) MATCH (n)-[:DNA_DIFFERENCE_GENE]-(m) WITH  m,COUNT (*) AS m_count WHERE m_count = 2 return count(*)'
#new--not so good...
#neotool cypher 'START n=node(239875,239872), m=node:Gene("name:*") WHERE (n)-[:DNA_DIFFERENCE_GENE]-(m) return count(*)'

@dajaxice_register
def match_samples(request, query):
    
    graph_db = neo4j.GraphDatabaseService()
    query_list = []
    
    sample_query = neo4j.CypherQuery(graph_db,'MATCH (n:LabID) WHERE n.name =~ "'+query['term']+'.*' +'" RETURN ID(n), n.name')
    
    for i in sample_query.execute().data:
        query_list.append({'id':i.values[0], 'text':i.values[1]})
    
    return json.dumps({"query":query, "data":{"results":query_list}}) 

@dajaxice_register
def get_sample_rels(request, cur_id):
    
    graph_db = neo4j.GraphDatabaseService()
    
    neo_query = 'MATCH (init_id:LabID{name:"'+cur_id+'"})-[:PRODUCED]-()-[:HAS_DISEASE]-(patient) WITH patient MATCH (patient)-[:HAS_DISEASE]-(disease)-[:PRODUCED]-(specimen) WITH patient, disease, specimen '+ \
    'OPTIONAL MATCH (specimen)-[ga:ALIAS_OF]-(geno_sample)-[var:DNA_DIFF]-() WHERE ga.alias_type = "genotype" WITH patient, disease, specimen, COUNT(var) AS Variants ' + \
    'OPTIONAL MATCH (specimen)-[sr:SIRNA_RUN]-() WITH patient, disease, specimen, Variants, COUNT(sr) AS siRNA ' + \
    'OPTIONAL MATCH (specimen)-[gr:GENE_SCORE_RUN]-() RETURN patient, disease, specimen, Variants,siRNA, COUNT(gr) AS GeneScore, CASE WHEN (Variants > 0) AND (siRNA > 0 OR COUNT(gr) > 0) THEN 1 ELSE 0 END AS required_data'
    
    sample_query = neo4j.CypherQuery(graph_db,neo_query)
    
    sample_list=[]
    
    for i in sample_query.execute().data:
        print i
        temp_dict = {}
        for j_ind,j in enumerate(i.values):
            if isinstance(j, int) and i.columns[j_ind] != "required_data":
                if temp_dict.has_key('data'):
                    temp_dict['data'].append({'name':i.columns[j_ind], 'count':j})
                else:
                    temp_dict['data'] = [{'name':i.columns[j_ind], 'count':j}]
            elif isinstance(j, int) and i.columns[j_ind] == "required_data":
                temp_dict["required_data"] = j
            else:
                temp_dict[i.columns[j_ind]] = j.get_properties()
        sample_list.append(temp_dict)
        
    return json.dumps({"sample_list":sample_list, "cur_id":cur_id})

@dajaxice_register
def save_parameters(request, save_name, filter_hs, param_hs, var_logic):
    #to be consistent with the posted data processed in views, need to make the value of each element be a list
    
    for i in filter_hs.keys():
        if isinstance(filter_hs[i], list) == False:
            filter_hs[i] = [filter_hs[i]]
    
    num_hs, group_hs, necessary_vars = views.process_request_post (filter_hs)
    new_filt_dict = views.make_updated_filter_dict (num_hs, group_hs, var_logic)
    
    print new_filt_dict
    
    exist_param = user_parameters.objects.filter(user=request.user, name=save_name)
    
    if len(exist_param) > 0:
        exist_param.delete();
        
    print new_filt_dict
    
    cur_db = user_parameters(user=request.user, name=save_name, filt=json.dumps(new_filt_dict), param=json.dumps(param_hs))
    cur_db.save()
    
    return json.dumps({"save_name":save_name})

@dajaxice_register
def load_parameters(request, load_name):
    exist_param = user_parameters.objects.filter(user=request.user, name=load_name)
    
    if len(exist_param) != 1:
        raise Exception("ERROR: only expected one database entry")
    
    return json.dumps({'filters':exist_param[0].filt, 'parameters':exist_param[0].param})

@dajaxice_register
def get_group_query (request, group_dict, other_node_ids, path_names, query_type, path_filters):
    
    graph_db = neo4j.GraphDatabaseService()
    link_list = []
    new_node_list = []
    existing_node_set = set()
    
    for i in path_names:
        group_query, group_params = make_group_query(group_dict, other_node_ids,query_type, i)
        
        data, cypher_meta = cypher.execute(graph_db,group_query, group_params)
        
        m_set = set()
        r_set = set()
        
        for j in data:
            print j
            m_node = j[0]
            r_rel = j[1]
            existing_node_set.add(j[2][j[2]['indexed_name']])
            # shouldb't be equal to any of the defined p's
            m_set.add(j[2][j[2]['indexed_name']])
            
            #record all the relationships between the m nodes and the existing nodes specified through p
            #can't compare r_rel to None directly, so this works for now
            if r_rel.__class__.__name__ != "NoneType":
                if r_rel.end_node._id != m_node._id and r_rel.start_node._id == m_node._id:
                    link_list.append({'source':m_node[m_node['indexed_name']], 'target':r_rel.end_node[r_rel.end_node['indexed_name']]})
                elif r_rel.start_node._id != m_node._id and r_rel.end_node._id == m_node._id:
                    link_list.append({'source':m_node[m_node['indexed_name']], 'target':r_rel.start_node[r_rel.start_node['indexed_name']]})
            
            #keep the unique m nodes (potentially new nodes to be added)
            if (m_node[m_node['indexed_name']] in m_set) == False:
                use_node = views.get_node(m_node['index'],m_node[m_node['indexed_name']], m_node['indexed_name'], graph_db, request)
                new_node_list.append(use_node)
                m_set.add(m_node[m_node['indexed_name']])
        
        
    inp_nodes = []
        
    for i in group_dict.keys():
        inp_nodes += group_dict[i]
    
    out_nodes = []
    
    for i in new_node_list:
        out_nodes.append(i['db_id'])
    
    path_pat = string.joinfields(path_names, " | ")
    
    #need to add in support for the cases where n,r,m does not exist but n and m do...
    n_m_query = 'START n=node({inp_nodes}), m=node({out_nodes}) MATCH (n)-[r:' + path_pat + ']-(m) RETURN n,m'
    
    print inp_nodes
    print out_nodes
    
    data, cypher_meta = cypher.execute(graph_db,n_m_query, {'inp_nodes':inp_nodes, 'out_nodes':out_nodes})
    
    for i in data:
        link_list.append({'source': i[0][i[0]['indexed_name']],'target':i[1][i[1]['indexed_name']]})
        existing_node_set.add( i[0][i[0]['indexed_name']])
    
   
    
    if len(existing_node_set) != (len(inp_nodes) + len(other_node_ids)):
        
        if len(set(inp_nodes) & existing_node_set) < len(inp_nodes):
            rem_node_query = 'START n=node({other_node_ids}) RETURN n'
            query_ids = list(set(inp_nodes).difference(existing_node_set))
        else:
            raise Exception("Unexpected contents of existing_node_set")
        
        data, cypher_meta = cypher.execute(graph_db,rem_node_query, {'other_node_ids':query_ids})
        
        for i in data:
            existing_node_set.add(i[0][i[0]['indexed_name']])
    
    return json.dumps({'link_list':link_list, 'new_node_list':new_node_list, 'existing_node_list':list(existing_node_set)})


@dajaxice_register
def get_path_nodes(request, cur_node_name, cur_node_id, other_node_ids, path_names, path_filters, keep_first, path_counts):
    if views.db_type == "neo4j":
        return get_path_nodes_neo4j(request, cur_node_name, cur_node_id, other_node_ids, path_names, path_filters, keep_first, path_counts)
    elif views.db_type == "MySQL":
        return get_path_nodes_mysql(request, cur_node_name, cur_node_id, other_node_ids, path_names, path_filters, keep_first, path_counts)
    else:
        raise Exception("Unexpected database type")

def get_path_nodes_mysql (request, cur_node_name, cur_node_id, other_node_ids, path_names, path_filters, keep_first, path_counts):
    graph_db = connections['varDb'].cursor()
    
    

def get_path_nodes_neo4j(request, cur_node_name, cur_node_id, other_node_ids, path_names, path_filters, keep_first, path_counts):
    
    graph_db = neo4j.GraphDatabaseService()
    link_list = []
    new_node_list = []
    existing_node_set = set()
    
    print path_counts
    
    #for each path, get the names of the nodes to retrieve attributes on
    for i in path_names:
        #path_where = make_filter_where_stat (path_filters[i])
        
        #kind of a hack for now...
        
        #query_str = make_query_from_path(i, cur_node_id, other_node_ids, path_where, graph_db, request, keep_first, path_counts[i])
        query_str = make_query_from_path(i, cur_node_id, other_node_ids, graph_db, request, keep_first, path_counts[i])
        
        m_set = set()
        r_set = set()
        cur_params = {'cur_node_id':cur_node_id, 'other_node_ids':other_node_ids}
        
        print query_str
        print cur_params
        data, cypher_meta = cypher.execute(graph_db,query_str, cur_params)
        
        for j in data:
            print j
            m_node = j[0]
            r_rel = j[1]
            existing_node_set.add(j[2][j[2]['indexed_name']])
            # shouldb't be equal to any of the defined p's
            m_set.add(j[2][j[2]['indexed_name']])
            
            #record all the relationships between the m nodes and the existing nodes specified through p
            #can't compare r_rel to None directly, so this works for now
            if r_rel.__class__.__name__ != "NoneType":
                if r_rel.end_node._id != m_node._id and r_rel.start_node._id == m_node._id:
                    link_list.append({'source':m_node[m_node['indexed_name']], 'target':r_rel.end_node[r_rel.end_node['indexed_name']]})
                elif r_rel.start_node._id != m_node._id and r_rel.end_node._id == m_node._id:
                    link_list.append({'source':m_node[m_node['indexed_name']], 'target':r_rel.start_node[r_rel.start_node['indexed_name']]})
            #keep the unique m nodes (potentially new nodes to be added)
            if (m_node[m_node['indexed_name']] in m_set) == False:
                use_node = views.get_node(m_node['index'],m_node[m_node['indexed_name']], m_node['indexed_name'], graph_db, request)
                link_list.append({'source': cur_node_name,'target':use_node['id']})
                new_node_list.append(use_node)
                m_set.add(m_node[m_node['indexed_name']])
        
        existing_node_set.add(cur_node_name)
        
        #If no nodes are returned as part of data, then the existing_node set will be empty (other than cur_node_name), causing issues as it controls which nodes are displayed.
        
        if len(existing_node_set) == 1:
            rem_node_query = 'START n=node({other_node_ids}) RETURN n'
            
            data, cypher_meta = cypher.execute(graph_db,rem_node_query, {'other_node_ids':other_node_ids})
            
            for i in data:
                existing_node_set.add(i[0][i[0]['indexed_name']])
        
    return json.dumps({'link_list':link_list, 'new_node_list':new_node_list, 'existing_node_list':list(existing_node_set)})#'cur_node_name':cur_node_name, 'cur_path_count':path_counts})

def make_filter_where_stat (path_filters):
    if len(path_filters.keys()) > 0:
        cur_where = []
        
        for i in path_filters.keys():
            #the presence of a 'concat' key indicates categorical, otherwise assume numeric
            #if concat = "", then only do exact matching...
            print path_filters
            if path_filters[i].has_key("concat"):
                
                range_perms = []
                
                if path_filters[i]["concat"] == "":
                    for j in path_filters[i]["range"]:
                        range_perms.append("'" + j + "'")
                else:
                    #Otherwise generate every permutation of range + concat + range
                    #need to modify permutations so that it comes up with every permutations from including 1...m elements
                    for j in get_all_permutations(path_filters[i]["range"]):
                        if len(j) > 0:
                            range_perms.append("'" + string.joinfields(j, path_filters[i]["concat"]) + "'")
                        
                cur_where.append('(r.' + i + ' IN [' + string.joinfields(range_perms, ",") + '])')
            else:
                cur_where.append('(r.' + i + ' IN RANGE(' + path_filters[i]["range"][0] + ',' + path_filters[i]["range"][1] + ',1))')
    
        where_stat = string.joinfields(cur_where, " AND ")
    else:
        where_stat = ''
    
    return where_stat
