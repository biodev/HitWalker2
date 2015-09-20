from py2neo import neo4j
import re
import collections
import json
import copy
import sys

try:
    import config
    from core import get_nodes
    import custom_functions
except:
    sys.exit(0)

class TestRequest():
    
    session = {}
    
    def __init__(self):
        
        for i in config.adjust_fields.keys():
            for j in config.adjust_fields[i]['fields'].items():
                self.session[j[0]] = j[1]['default']
        
        self.session['where_vars'] = []

if __name__ == '__main__':
    
    graph_db = neo4j.GraphDatabaseService(config.cypher_session+'/db/data/') 
    
    samp_vals = collections.defaultdict(list)
    
    for i in ['gene_names', 'subject', 'sample', 'pathway']:
        cur_query = getattr(config, i)
        
        whole_query = re.sub("\{name:.+\}", "", cur_query['query'])
        
        whole_query = re.sub("WHERE.+RETURN", "RETURN", whole_query)
        
        if i == 'gene_names':
            use_query = whole_query + " LIMIT 10"
        else:
            use_query = whole_query + " LIMIT 1"
        
        query_res = neo4j.CypherQuery(graph_db, use_query).execute()
        
        for j in query_res:
            
            temp_list = []
            
            for k_ind, k in enumerate(j.columns):
                if k.endswith(".name"):   
                    temp_list.append(j.values[k_ind])
            
            #if nothing has an explicit name attribute try to get the names from the Node objects
            if len(samp_vals[i]) == 0:
                for k_ind, k in enumerate(j.values):
                    if isinstance(k, neo4j.Node):
                        temp_list.append(k['name'])
            
            if len(temp_list) == 1:
                temp_list = temp_list[0]
            samp_vals[i.replace('_names', '')].append(temp_list)
    
    #basic node retrieval

    print samp_vals

    request = TestRequest()
    
    param_key = collections.defaultdict(list)
    
    for i in config.data_list:
        param_key[i].append(samp_vals['subject'][0])
    
    for i in ['Gene', 'Sample', 'Subject']:
        
        cur_val = samp_vals[i.lower()][0]
        
        if isinstance(cur_val, list):
            cur_val = cur_val[0]
        print cur_val
        print param_key
        print get_nodes([cur_val], i, request, param_list=[param_key]).display_names()
    
    #gene, subject and pathway searches
    
    for i in ['sample','gene','pathway']:
        
        use_str = ''
        print samp_vals[i]
        #for gene, the second entry is the symbol
        if len(samp_vals[i]) > 1:
            use_str= samp_vals[i][0][1]
        else:
            use_str = samp_vals[i][0]
        
        if i == 'pathway':
            use_val = json.dumps([[samp_vals['gene'][0]]])
        else:
            use_val = use_str
        
        print getattr(custom_functions, 'match_'+i)(use_val)
    
    #template queries
    
    param_dict = copy.deepcopy(request.session)
    
    param_dict.update({'Gene':samp_vals['gene'], 'Subject':samp_vals['subject']})
    
    for i in config.node_group_content.items():
        for j in i[1]['options']:
            if j['query'] != '':
                print list(neo4j.CypherQuery(graph_db, j['query']).execute(**param_dict))
    
    #carry out the shortest path query
    
    request_post = {'var_select':map(lambda x: x[0], samp_vals['gene'][:3]), 'seed_select':map(lambda x: x[0], samp_vals['gene'][7:])}
    
    print request_post
    
    request.session['query_samples'] = {'SampleID':param_key}
    
    print custom_functions.get_shortest_paths(request, request_post)
    
    #finally carry out the network query
    
    print custom_functions.gene_to_gene(request_post['var_select'], request_post['seed_select'], request, None, {})
    