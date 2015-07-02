import core
import custom_functions
import string

##globals

prog_type = ""

#test_methods_path="/Users/bottomly/Desktop/github/HitWalker2/populate/ccle_test_methods.R"
#hw_config_path="/Users/bottomly/Desktop/hitwalker2_paper/temp_hw_conf.RData"
#chrome_driver_path="/Users/bottomly/Desktop/chromedriver"

graph_struct_file = "/var/www/hitwalker2_inst/static/network/data/graph_struct.json"

cypher_session="http://localhost:7474"
#the maximum number of nodes in a group before it becomes a metanode
max_nodes = 1

subj_att_set = @SUBJECT_ATTRIBUTES@

#sizes for network view
network_sizes = {"w":400, "h":400, "legend_offset":200, "history_offset":10}

#sizes for pathway view
pathway_sizes = {"w":800, "h":800, "legend_offset":200, "history_offset":10}

#the seeds need to refer to all the variables that are used as 'seeds' for the ranking algorithm
#the target is the variable to be ranked (otherwise known as the query)
data_types=@DATA_TYPES@

#The data_list is simply a list of all data types whether they are used in the prioritization or not
data_list=@USE_DATA@

#fields part of the 'Adjust' dropdown on the main page that accept user input
#this input controls which genes get prioritized and how the queries are performed
#input should be in the form:
#button_name:{
#   {field_id:{type:('standard' or 'grouped'), fields:[{type:('numeric' or 'character')}, range:[numeric1, numeric2], name:string]}}    
#}
#there are other fields depending on whether these are standard input fields or should be treated as logical groups
#logical fields require more information on the graph database in order to ensure they are executed appropriately
#   The grouped type additionally needs a 'trans' field which is a function which transforms the variables prior to use in a query
#   Additionally if 'needs_has' is specified, regardless of value, a HAS statement will be added.
#   It also needs a 'var_name' field which indicates what variable should be present in the query and a 'required' field which is a dictionary {'from':'node or relationship name'}

#For the standard fields, they are stored in the session as their field names and can be retrieved as such in the below queries using the lists specified in the 'session_params' portion.
#Currently these are opt-in.

adjust_fields = {
    'General_Parameters':{'type':'standard',
                          'fields':{
                                'string_conf':{'type':'numeric', 'default':.7, 'range':[.7, 1], 'comparison':'>', 'name':'HitWalker STRING Confidence'},
                                'path_conf':{'type':'numeric', 'default':.95, 'range':[0, 1], 'comparison':'>', 'name':'Pathway STRING Confidence'},
                                'res_prob':{'type':'numeric', 'default':.3, 'range':[0,1], 'comparison':'=', 'name':'Restart Probability'},
                                'max_iter':{'type':'numeric', 'default':100, 'range':[0, 10000], 'comparison':'=', 'name':'Max Iterations'},
                                'conv_thresh':{'type':'numeric', 'default':1e-10, 'range':[0,1], 'comparison':'<', 'name':'Convergence Threshold'}@HIT_PARAMS@
                                }, 
                        }
}

#Basic info for genes:
#Needs the query to return info in the following form:
    #gene ID, gene symbol, and the collection of pathways (or [] if that info is not available).
#gene_names = {'query':'MATCH (n:Gene{name:{GENE}})-[r:KNOWN_AS]-(m) WHERE r.status="symbol" WITH n,m OPTIONAL MATCH (n)<-[:EXTERNAL_ID]-()<-[:GENESET_CONTAINS]-(p) RETURN n.name,m.name, COLLECT(DISTINCT p.name)', 'handler':custom_functions.get_gene_names, 'session_params':None}

gene_names = {'query':'MATCH (n:EntrezID{name:{GENE}})-[r:REFFERED_TO]-(m) RETURN n.name,m.name, []', 'handler':custom_functions.get_gene_names, 'session_params':None}

#used for the network view
gene_rels = {'query':'MATCH (gene:EntrezID{name:{FROM_GENE}})-[:MAPPED_TO]-(string_from)-[r:ASSOC]-(string_to)-[:MAPPED_TO]-(gene_to) WHERE gene_to.name IN {TO_GENES} AND HAS(r.score) AND r.score > ({string_conf}*1000) RETURN gene.name,gene_to.name, MAX(r.score)', 'handler':None, 'session_params':[['string_conf']]}

subject = {'query':'MATCH (n:@SUBJECT@{name:{SUBJECTID}}) RETURN n', 'handler':custom_functions.get_subject, 'session_params':None}
#The initial call after the user chooses a sample can be made in terms of sample names not necessarily subject names
sample = {'query':'MATCH (n)-[:DERIVED]->(Sample{name:{SUBJECTID}}) RETURN n', 'handler':custom_functions.get_subject, 'session_params':None}

pathway = {'query':'MATCH (path:Pathway{name:{PATHNAME}})-[:PATHWAY_CONTAINS]->(gene) RETURN path.name, COLLECT(DISTINCT gene.name)', 'handler':custom_functions.get_pathway, 'session_params':None}


#Basic node info for the seeds:

#needs to return data in the following form:
    #gene_name, sample_name, var (corresponding to datatypes), and other metadata...

@BASE_QUERIES@

@TEMPLATE_QUERIES@

#also need to add in expression


#gene_exprs_low = {'query':'MATCH(n:LabID)-[r:EXON_ARRAY_RUN]-()-[r2:ASSIGNED]-(m:Gene{name:{GENE}}) WHERE HAS(r.score) AND n.name IN {LABID} RETURN m.name AS gene, n.name AS sample, "LowExpr" AS var, MAX(r.score) AS score, ALL(x IN COLLECT(r.score) WHERE x < {LOWEXPRS}) AS is_hit',
#          'handler':custom_functions.get_exprs, 'session_params':None}
#

#sample relationships

#needs to have $$sample$$ which will be supplied by the 'text' return value of matchers['sample']
sample_rels_query = """@REL_QUERY_STR@"""
sample_rels_type = 'hierarchical'

##This specifies the functions to be used for searching on genes, pathways and samples
matchers = {
    'pathway':custom_functions.match_pathway,
    'gene':custom_functions.match_gene,
    'subject':custom_functions.match_sample,
    'sample':custom_functions.match_sample
}

##The queries used to choose the data used as seeds or otherwise in the prioritization algorithm
#core.handle_hits needs the following:

#{name} is necessary here as the sample name will come directly from getNodes
hit_session_dict = {}

for i in data_types['seeds']:
    hit_session_dict[i] = [core.customize_query(eval(i), query=lambda x:x.replace("{name:{GENE}}", "").replace("{SAMPLE}", "{name}") , handler=lambda x: core.handle_hits)]

#Here we also need to have 'query_ind', 'gene_ind' and a unique 'row_id'
query_prior_dict = dict([(data_types['target'], [core.customize_query(eval(data_types['target']), query=lambda x:x.replace("{name:{GENE}}", "").replace("{SAMPLE}", "{name}") , handler=lambda x: core.handle_query_prior)])])

score_hits=core.no_combining
convert_ids_to=custom_functions.gene_seed_list_to_protein
convert_ids_from=custom_functions.protein_seed_list_to_gene
get_seed_list=custom_functions.make_seed_list
get_query_list=custom_functions.make_seed_list
#
##Need to ensure that the query can find the appropriate parameters in the session and that the sample is refered to as '{name}'
prioritization_func={'function':custom_functions.netprop_rwr, "args":{"initial_graph_file":"/var/www/hitwalker2_inst/static/network/data/9606.protein.links.v9.1.mm.mtx",
                                                                      "string_session_name":"string_conf",
                                                                      "res_prob_session_name":"res_prob",
                                                                      "conv_thresh_session_name":"conv_thresh",
                                                                      "max_iter_session_name":"max_iter"
                                                                        }}



#['query_samples', 'LabID', 'GeneScore']
#['query_samples', 'LabID', 'siRNA']
#lambda x: [['query_samples', 'LabID', 'Variants']]

node_queries={
    'Gene':[core.customize_query(gene_names, query=lambda x: x.replace("{GENE}", "{name}"))],
    'Sample':[core.customize_query(sample, query=lambda x: x.replace("{SUBJECTID}", "{name}"))],
    'Subject':[core.customize_query(subject, query=lambda x: x.replace("{SUBJECTID}", "{name}"))]
}

#by default the user has no control over these parameters, if this was desired then these queries would need to be specified in 'adjust_fields' above and the unique id would need to be specified in session_params (e.g. core.customize_query(etc, etc, session_params=lambda x: [['gene_score']])
#assuming gene_score was the unique key in 'adjust_fields'


for i in data_list:
    node_queries['Gene'].append(core.customize_query(eval(i), query=lambda x: x.replace("{GENE}", "{name}").replace("{SAMPLE}", "{"+i+"}")))


#Whereas the node_queries are only involved in the get_shortest_paths functionality
edge_queries = {
    'nodes':{
        'Gene':[core.customize_query(gene_names, query=lambda x: x.replace("{GENE}", "{name}"))],
        'Subject':[core.customize_query(subject, query=lambda x: x.replace("{SUBJECTID}", "{name}"))],
        'Pathway':[core.customize_query(pathway, query=lambda x: x.replace("{PATHNAME}", "{name}"))]
        },
    
    'edges':{
        'Gene':{
            'Subject':{'query':None, 'handler':custom_functions.gene_to_sample, 'session_params':None},
            'Gene':{'query':None, 'handler':custom_functions.no_rels, 'session_params':None}
            },
        'Subject':{
            'Subject':{'query':None, 'handler':custom_functions.no_rels, 'session_params':None}
        }
    }

}

for i in data_list:
    edge_queries['nodes']['Gene'].append(core.customize_query(eval(i), query=lambda x: x.replace("{GENE}", "{name}")))


##This dictionary specifies how the summary modals should be displayed upon right-clicking a specified node
node_content = {
    'Subject':{'func':core.make_sample_table, 'args':{}},
    'Gene':{'func':core.make_gene_table, 'args':{'gene_link':'http://www.ncbi.nlm.nih.gov/gene/?term='}},
    'MetaNode':{'func':core.make_metanode_table, 'args':{}}
}

node_group_content={
    'Gene':{'title':'Find subjects with hits in the largest fraction of these genes using:',
            'returned_node_type':'Subject',
            'options':[]
            },
    'Subject':{'title':'Find the most frequent gene hits for these subjects using:',
            'returned_node_type':'Gene',
            'options':[]
              }
}

group_key_set = set(node_group_content.keys())
assert len(group_key_set) == 2

for i in group_key_set:
    for j in data_list:
        node_group_content[i]['options'].append(core.specify_type_query_tmp(eval(j+"_tmpl"), ret_type=list(group_key_set.difference(set([i])))[0], coll_type=i))

#needs to be a list of the form: [(Gene | Sample, value in data_list to disable)]
disable_node_group_content = []

node_abbreviations = {}

edge_abbreviations = {}

#
#plot_data_types = {
#    'siRNA':custom_functions.get_sirna_data
#}
#

graph_initializers = {
    'panel':custom_functions.get_shortest_paths,
    'image':custom_functions.get_pathways_sample
}
#
