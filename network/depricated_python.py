@dajaxice_register
def match_samples(request, query):
    
    graph_db = neo4j.GraphDatabaseService()
    query_list = []
    
    sample_query = 'START n=node:Alias("name:'+ query['term']+'*' +'") WHERE n.valid_sample = true RETURN ID(n), n.name';
    data, cypher_meta = cypher.execute(graph_db,sample_query)
    
    for i in data:
        query_list.append({'id':i[0], 'text':i[1]})
    
    return json.dumps({"query":query, "data":{"results":query_list}})

##FILTER=<ID=GATKStandardFilter,Description="QUAL < 30.0 || QD < 5.0 || HRun > 5">
##FILTER=<ID=HARD_TO_VALIDATE,Description="MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=QualFilter,Description="QUAL < 10">
##FILTER=<ID=StrandBiasFilter,Description="SB >= -1.0">

#example values, should determine the real ones from the database upon updating also need a validation function
#field_dict = {
#    'DNA_DIFFERENCE':
#    {'allele_count': {'type':'numeric', 'range':[0,2], 'default':[0,2]},
#    'genotype_quality': {'type':'numeric', 'range':[0,100], 'default':[0,100]}},
#    'UNIT_DNA_DIFFERENCE':
#    {'HRun': {'type':'numeric', 'range':[0,10000], 'default':[0,10000]},
#    'MQ0':{'type':'numeric', 'range':[0,10000], 'default':[0,10000]},
#    'MQ': {'type':'numeric', 'range':[0,100], 'default':[0,100]},
#    'QD': {'type':'numeric', 'range':[0,10000], 'default':[0,10000]},
#    'SB': {'type':'numeric', 'range':[-10000,10000], 'default':[-10000,10000]},
#    'SB': {'type':'numeric', 'range':[0,10000], 'default':[0,10000]},
#    'FILTER': {'type':'character', 'range':['GATKStandardFilter', 'PASS', 'StrandBiasFilter'], 'concat':';', 'default':'PASS'}},
#    'IMPACTS':{'Cons_cat':{'type':'character', 'range':['Synonymous', 'NonSynonymous', 'Other'], 'concat':'', 'default':'NonSynonymous'}}
#}
## 
#
#default_filters = [
#            [
#                [{'filter':'Cons_cat', 'default':'NonSynon.'},    
#                {'filter':'genotype_quality', 'comparison':'>', 'default':40, 'logical':'AND'},
#                {'filter':'QD', 'comparison':'>', 'default':5, 'logical':'AND'},
#                {'filter':'HRun', 'comparison':'<', 'default':5, 'logical':'AND'},
#                {'filter':'MQ0', 'comparison':'<', 'default':4, 'logical':'AND'},
#                {'filter':'cohort_freq', 'comparison':'<', 'default':.5, 'logical':'AND'}]
#            ],
#            [
#                [{'filter':'in_1kg', 'default':'False', 'logical':'AND'},{'filter':'in_dbsnp', 'default':'False', 'logical':'AND'}],
#                [{'filter':'in_1kg', 'default':'True', 'logical':'OR'},{'filter':'freq', 'comparison':'<','default':.01, 'logical':'AND'}],
#                [{'filter':'in_dbsnp', 'default':'True', 'logical':'OR'},{'filter':'in_1kg', 'default':'False', 'logical':'AND'},{'filter':'cohort_count', 'comparison':'=', 'default':1, 'logical':'AND'}]
#            ]
#]

#default_parameters = {
#    'zscore':{'default':-2,'range':[-100, 0]},
#    'gene_score':{'default':0, 'range':[-100,100]},
#    'string_conf':{'default':.4, 'range':[0, 1]},
#    'res_prob':{'default':.3, 'range':[0,1]},
#    'max_iter':{'default':100, 'range':[0, 10000]},
#    'conv_thresh':{'default':1e-10, 'range':[0,1]}
#}


#filter_dict = {
#        'Sample':{
#            'DNA_DIFFERENCE':{'to':['in_1kg', 'in_dbsnp', 'freq'],'edge':['AD_1', 'AD_2', 'depth', 'allele_count', 'genotype_quality']}
#            },
#        'Variation':{
#            'UNIT_DNA_DIFFERENCE':[{'from':['in_1kg', 'in_dbsnp','freq'], 'edge':['HRun', 'MQ0', 'MQ', 'QD', 'SB']}],
#            'DNA_DIFFERENCE':[{'from':['in_1kg', 'in_dbsnp', 'freq'],'edge':['AD_1', 'AD_2', 'depth', 'allele_count', 'genotype_quality']}],
#            'IMPACTS':[{'from':['in_1kg', 'in_dbsnp','freq'], 'edge':['Cons_cat']}]
#        }
#            }



#graph_struct = {'Sample':{
#                    'DNA_DIFFERENCE':'Variation',
#                    'GENOTYPED_USING':'Experiment',
#                    'ALIAS_OF':'Alias'},
#                'Alias':{
#                  'ALIAS_OF':'Sample'
#                },
#                'Experiment':{
#                    'UNIT_DNA_DIFFERENCE':'Variation',
#                    'GENOTYPED_USING':'Sample'
#                },
#                'Variation':{
#                    'IMPACTS':'Transcript',
#                    'DNA_DIFFERENCE':'Sample',
#                    'UNIT_DNA_DIFFERENCE':'Experiment'
#                },
#                'Transcript':{
#                    'TRANSCRIBED':'Gene',
#                    'IMPACTS':'Variation'
#                },
#                'Gene':{
#                    'KNOWN_AS':'Symbol',
#                    'TRANSCRIBED':'Transcript'
#                },
#                'Symbol':{
#                    'KNOWN_AS':'Gene'
#                    }}

#query_dict = {
#    'Gene':{
#        
#        'self_query':'SELECT stable_id AS id, display_label AS display_name FROM homo_sapiens_core_70_37.gene JOIN homo_sapiens_core_70_37.xref ON display_xref_id = xref_id WHERE stable_id = %s',
#            
#        'path_query':{
#                'Sample':{
#                'path':'DNA_DIFFERENCE_GENE',
#                'base_query':'SELECT DISTINCT samp_ind AS Sample FROM sample_consequence WHERE Gene = %s',
#                'count_query':'SELECT COUNT(DISTINCT samp_ind) AS Sample_count FROM sample_consequence WHERE Gene = %s'
#                }
#            }
#        
#        
#    },
#    'Sample':{
#        
#        'self_query':'SELECT samp_ind AS id, alias AS display_name FROM sample_alias WHERE samp_ind=%s',
#            
#        'path_query':{
#            'Gene':{
#            'path':'DNA_DIFFERENCE_GENE',
#            'base_query': 'SELECT DISTINCT Gene AS Gene FROM sample_consequence WHERE samp_ind = %s',
#            'count_query': 'SELECT COUNT(DISTINCT Gene) AS Gene_count FROM sample_consequence WHERE samp_ind = %s'
#            } 
#        }
#        
#        
#    }  
#}

#from the django docs
def dictfetchall(cursor):
    "Returns all rows from a cursor as a dict"
    desc = cursor.description
    return [
        dict(zip([col[0] for col in desc], row))
        for row in cursor.fetchall()
    ]

def att_count_func_mysql (cur_data, index, indexed_name, name, graph_db, request):
    att_dict = {'Sample':{'name':'Sample', 'func':return_counts},
                'Transcript':{'name': 'Transcript', 'func':return_counts},
                'Variation':{'name':'Variation', 'func':return_counts},
                'Gene':{'name': 'Gene', 'func':return_counts},
                'Experiment':{'name':'Experiment', 'func':return_counts},
                'EntrezID':{'name':'Gene', 'func':get_graph_counts}}
    
    #can add a 'blacklist' for att_dict like:  'Gene':{'name': 'Gene', 'func':return_counts,'blacklist':set(['Sample'])}
   
    ret_dict = dict()
    
    for i in cur_data:
        if att_dict.has_key(i[0]) == True:
            if (att_dict[i[0]].has_key('blacklist') == False) or (att_dict[i[0]].has_key('blacklist') == True and (index in att_dict[i[0]]['blacklist']) == False):
                count = att_dict[i[0]]['func'](i, index, indexed_name, name, graph_db, request)
                ret_dict[i[0]] = {'link':i[1], 'count':count, 'name':att_dict[i[0]]['name'], 'cur_pos':0}
   
    return ret_dict

def get_path_count_atts_mysql (index, indexed_name, name, graph_db, request):
    
    data = []
    for i in query_dict[index]['path_query'].iteritems():
        print i
        graph_db.execute(i[1]['count_query'], [name])
        res = dictfetchall(graph_db)
        if len(res) != 1:
            raise Exception("Unexpected length for the database result")
        data.append([i[0], i[1]['path'], res[0][i[0]+'_count']])
    
    return att_count_func(data, index, indexed_name, name, graph_db, request)

def get_simple_node_summary_mysql (index, name, indexed_name, graph_db, request, other_dict):
    
    att_dict = get_path_count_atts_mysql(index, indexed_name, name, graph_db, request)
    print name
    graph_db.execute(query_dict[index]['self_query'], [name])
    
    data = dictfetchall(graph_db)
    
    print data
    
    if len(data) == 1:
        display_name = data[0]['display_name']
        db_id = data[0]['id']
        
        data[0].pop('display_name')
        data[0].pop('id')
        
        meta = data[0]
    else:
        display_name = name
        meta = {}
        db_id = None
    
    return {'display_name': display_name, 'id':name, 'db_id':db_id,'attributes':{'node_type':index, 'indexed_name':indexed_name, 'paths':att_dict, 'meta':meta, 'add_string':''}}

def gene_node_summary_mysql (index, name, indexed_name, graph_db, request, keydict):
    m_node_summary = get_simple_node_summary_mysql(index, name, indexed_name, graph_db, request, keydict)
    
    if request.session['hitwalker_score'].has_key(name):
        if m_node_summary.has_key('node_score') == False:
            m_node_summary['node_score'] = request.session['hitwalker_score'][name]
        else:
            raise Exception("ERROR: node_score attribute already exists in dictionary")
    
    #need to revisit siRNA and gene hits, probably do the query upfront and pass along in request.session as above...
    
    return m_node_summary

group_fore_dict = {'Gene':{'on_expand':[], 'on_follow':['ASSOC', 'DNA_DIFFERENCE_GENE']},
                    'Variant': {'on_expand':[], 'on_follow':[]},
                    'Sample':{'on_expand':[], 'on_follow':[]}}

def make_rev_dict (inp_dict):
    
    rev_dict = dict()
    
    for i in inp_dict.keys():
        for j in inp_dict[i].keys():
            if rev_dict.has_key(j):
                if rev_dict[j].has_key(inp_dict[i][j]):
                    rev_dict[j][inp_dict[i][j]].append(i)
                else:
                    rev_dict[j][inp_dict[i][j]] = [i]
            else:
                rev_dict[j] = {inp_dict[i][j]:[i]}
    
    return rev_dict

#group_rev_dict = make_rev_dict(group_fore_dict)

#print group_rev_dict

 #graph_db = neo4j.GraphDatabaseService()
 #       
 #       sample_query_str = 'START n=node({sample_alias}) MATCH (n)-[:ALIAS_OF]->(m) WITH m MATCH (m)<-[r:ALIAS_OF]-(p) WHERE r.display_type = "display" return m,p'
 #       
 #       data, cypher_meta = cypher.execute(graph_db,sample_query_str, {'sample_alias':int(request.session['sample_alias'])})
 #       
 #       ret_samp = []
 #       
 #       for i in data:
 #           ret_samp.append(i)
 #       
 #       if len(ret_samp) == 0:
 #           return HttpResponseRedirect(reverse('indexSampleNotFound'))
 #       elif len(ret_samp) > 1:
 #           return HttpResponseRedirect(reverse('indexAmbigous'))
 #       else
        
def proper_type_dep(var):
    try:
        cur_type = type_of_value(var)
        
        if cur_type == type(int()):
            return int(var)
        elif cur_type == type(float()):
            return float(var)
        elif cur_type == type(None):
            return None
        else:
            return str(var)
    except Exception:
        return str(var)
    
def return_counts (cur_data,index, indexed_name, name, graph_db, request):
    return cur_data[2]

def get_graph_counts (cur_data,index, indexed_name, name, graph_db, request):
    
    edge_thresh = float(request.session['string_conf'])*1000
    
    #note that assoc is undirected, so the edges would be counted twice with a direction constraint
    string_name_query = 'START n=node:' + index + '({indexed_name}={name}) MATCH (n)-[:EXTERNAL_ID]-()-[:MAPPED_TO]-()-[r:ASSOC]->()-[:MAPPED_TO]-()-[:EXTERNAL_ID]-(m) WHERE r.score > {score} RETURN COUNT (DISTINCT m)'
    
    data, cypher_meta = cypher.execute(graph_db,string_name_query, {'indexed_name':indexed_name, 'name':name, 'score':edge_thresh})
    print data
    if len(data) > 1:
        print "ERROR: Unexpected number of return values in get_graph_counts"
        count = -999
    else:
        count = data[0][0]
    
    return count

def att_count_func (cur_data, index, indexed_name, name, graph_db, request):
    
    att_dict = {'Sample':{'name':'Sample', 'func':return_counts, 'blacklist':set(['Experiment'])},
                'Gene':{'name': 'Gene', 'func':return_counts},
                'Experiment':{'name':'Sample', 'func':return_counts},
                'EntrezID':{'name':'Gene', 'func':get_graph_counts}
            }
    
    #can add a 'blacklist' for att_dict like:  'Gene':{'name': 'Gene', 'func':return_counts,}
   
    ret_dict = dict()
    
    for i in cur_data:
        if att_dict.has_key(i[0]) == True:
            if (att_dict[i[0]].has_key('blacklist') == False) or (att_dict[i[0]].has_key('blacklist') == True and (index in att_dict[i[0]]['blacklist']) == False):
                count = att_dict[i[0]]['func'](i, index, indexed_name, name, graph_db, request)
                ret_dict[i[0]] = {'link':i[1], 'count':count, 'name':att_dict[i[0]]['name'], 'cur_pos':0}
   
    return ret_dict

def get_path_count_atts (index, indexed_name, name, graph_db, request):
    
    query_str = 'START n=node:' + index + '({indexed_name}={name}) MATCH (n)-[r]-(m)  RETURN m.index, TYPE(r), COUNT(*)'
    data, cypher_meta = cypher.execute(graph_db,query_str, {'indexed_name':indexed_name, 'name':name})
    
    return att_count_func(data, index, indexed_name, name, graph_db, request)



def get_simple_node_summary (index, name, indexed_name, graph_db, request):
    
    att_dict = get_path_count_atts(index, indexed_name, name, graph_db, request)
    
    display_name_str = 'START n=node:' + index + '({indexed_name}={name}) RETURN n'
    data, cypher_meta = cypher.execute(graph_db,display_name_str, {'indexed_name':indexed_name, 'name':name})
    
    if len(data) == 1:
        display_name = data[0][0][data[0][0]['indexed_name']]
        db_id = data[0][0]._id
        meta = data[0][0].get_properties()
        meta.pop("index")
        meta.pop(meta['indexed_name'])
        meta.pop("indexed_name")
    else:
        display_name = name
        meta = {}
        db_id = None
    
    return {'display_name': display_name, 'id':name, 'db_id':db_id,'attributes':{'node_type':index, 'indexed_name':indexed_name, 'paths':att_dict, 'meta':meta, 'add_string':''},'children':[]}

def get_multi_node_summary (index, name, indexed_name, graph_db, request, path_name, where_statement):
    
    att_dict = get_path_count_atts(index, indexed_name, name, graph_db, request)
    
    #also add in node specific info
    
    display_name_str = 'START n=node:' + index + '({indexed_name}={name}) MATCH (n)-[r:' + path_name + ']-(m) ' + where_statement + ' RETURN n,m'
    print display_name_str
    data, cypher_meta = cypher.execute(graph_db,display_name_str,{'indexed_name':indexed_name, 'name':name})
    print data
    if len(data) == 1:
        display_name = data[0][1][data[0][1]['indexed_name']]
        
        db_id = data[0][0]._id
        meta = data[0][0].get_properties()
        meta.pop("index")
        meta.pop(meta['indexed_name'])
        meta.pop("indexed_name")
    
    else:
        raise Exception("Expected a single result for get_multi_node_summary")
    
    return {'display_name': display_name, 'id':name ,'db_id':db_id, 'attributes':{'node_type':index, 'indexed_name':indexed_name, 'paths':att_dict, 'meta':meta, 'add_string':''}, 'children':[]}
    
def get_node_summary (index, name, indexed_name, graph_db, request, keydict):
    
    if len(keydict.keys()) == 0:
        return get_simple_node_summary(index, name, indexed_name, graph_db, request)
    elif keydict.has_key('path_name') and keydict.has_key('where_statement'):
        return get_multi_node_summary(index, name, indexed_name, graph_db, request, keydict['path_name'], keydict['where_statement'])

def gene_node_summary (index, name, indexed_name, graph_db, request, keydict):
    m_node_summary = get_multi_node_summary (index, name, indexed_name, graph_db, request, keydict['path_name'], keydict['where_statement'])
    #then also determine if it is a siRNA/gene score hit and what the hitwalker score is for the current patient
    
    #if request.session['hitwalker_score'].has_key(name):
    #    if m_node_summary.has_key('node_score') == False:
    #        m_node_summary['node_score'] = request.session['hitwalker_score'][name]
    #    else:
    #        raise Exception("ERROR: node_score attribute already exists in dictionary")
    
    ###moved hitwalker score attributes to the links
    m_node_summary['node_score'] = 1
    
    #need to make this specific to the patient in question so limit m to only the current patient
    hit_query = 'START n=node: ' + index + '({indexed_name}={name})' + ', m=node({alias_node}) MATCH (n)-[r:GENE_SCORE|SIRNA_HIT]-(m) WHERE (HAS(r.zscore) AND r.zscore < {zthresh}) OR (HAS(r.gene_score) AND r.gene_score > {gene_score_thresh}) RETURN TYPE(r), r.zscore?, r.gene_score?'
    
    data, cypher_meta = cypher.execute(graph_db,hit_query, {'name':name, 'indexed_name':indexed_name, 'zthresh':float(request.session['zscore']), \
                                                            'gene_score_thresh':float(request.session['gene_score']), 'alias_node':request.session['alias_node']})
    
    if m_node_summary.has_key('children') == False:
        raise Exception ("children key expected in m_node_summary")
    
    if len(data) > 0:
        for i in data:
            if i[0] == 'GENE_SCORE':
                m_node_summary['children'].append({'attributes':{'node_type':'Gene_Score', 'score':data[0][2]}, 'display_name':m_node_summary['display_name']+'_Gene_Score', 'id':name+'_Gene_Score'})
            elif i[0] == 'SIRNA_HIT':
                m_node_summary['children'].append({'attributes':{'node_type':'siRNA_Hit','score':data[0][1]}, 'display_name':m_node_summary['display_name']+'_siRNA_Hit', 'id':name+'_siRNA_Hit'})
            else:
                raise Exception("Unexpected type of hit")
    
    #get the variants assigned to the gene for the given sample...
    
    gene_var_query = 'START n=node({gene_id}),m=node({alias_node}) MATCH (n)-[:TRANSCRIBED]-(o)-[r:IMPACTS]-(var) WITH n,o,var,r,m MATCH (m)-[:ALIAS_OF]-(sample)-[r2:DNA_DIFFERENCE]-(var)-[u:UNIT_DNA_DIFFERENCE]-(exp) WHERE (exp)<-[:GENOTYPED_USING]-(sample) RETURN var,o,r,u,exp'
    
    use_gene_var_query = add_where_input_query({'n':'Gene', 'm':'Alias'}, gene_var_query, request.session['where_template'], request.session['necessary_vars'])
    
    data, cypher_meta = cypher.execute(graph_db,use_gene_var_query, {'gene_id':m_node_summary['db_id'], 'alias_node':request.session['alias_node']})
    
    for i in data:
        temp_dict = {'attributes':{'node_type':'Variant'}, 'display_name':i[0]['name'], 'id':name+'_'+i[1]['name']+'_'+i[0]['name']}
        for j in i:
            for k in j.get_properties().items():
                temp_dict['attributes'][k[0]] = k[1]
        m_node_summary['children'].append(temp_dict)
    return m_node_summary

#old query for loucy'Sample':{'func':get_node_summary, 'argdict':{'path_name':'ALIAS_OF', 'where_statement':'WHERE r.status="display"'}}
def get_node(node_type, name, indexed_name, graph_db, request):
    
    db_type = 'neo4j'
    
    func_dict = {'neo4j':{
                'Gene':{'func':gene_node_summary, 'argdict':{'path_name':'KNOWN_AS', 'where_statement':'WHERE r.status="symbol"'}},
                'Sample':{'func':get_node_summary, 'argdict':{'path_name':'ALIAS_OF', 'where_statement':'WHERE r.display_type="display"'}},
                'Symbol':{'func':get_node_summary,'argdict':{}},
                'Variation':{'func':get_node_summary,'argdict':{}},
                'Transcript':{'func':get_node_summary,'argdict':{}},
                'Experiment':{'func':get_node_summary,'argdict':{}}}
            }
    #where db_type is a global variable...
    if func_dict.has_key(db_type) and func_dict[db_type].has_key(node_type):
        return func_dict[db_type][node_type]['func'](node_type, name, indexed_name, graph_db, request, func_dict[db_type][node_type]['argdict'])
    else:
        return None
    
def get_link_atts (request, from_node, to_node):
    if from_node['attributes'].has_key('node_type') and to_node['attributes'].has_key('node_type'):
        #at least find the name of the edge
        #still need to implement this
        #graph_db = neo4j.GraphDatabaseService()
        
        #query_str = 'START n=node:'+from_node["attributes"]["node_type"]+'(' +from_node["attributes"]["indexed_name"]  + '="' + from_node["id"] + '"' )
        
        if from_node["attributes"]["node_type"] == "Sample" and to_node["attributes"]["node_type"] == "Gene" and request.session['hitwalker_score'].has_key(to_node["id"]):
            return {'type':'HitWalker_Rank', 'rank':request.session['hitwalker_score'][to_node["id"]]+1}
        
        elif from_node["attributes"]["node_type"] == "Sample" and to_node["attributes"]["node_type"] == "Gene" and request.session['hit_dict'].has_key(to_node["id"]):
            return {'type':'Hit_Rank', 'rank':request.session['hit_dict'][to_node["id"]]}
        else:
            return {'type':'Unknown', 'rank':1000000}
    else:
        return {'type':'Unknown', 'rank':1000000}