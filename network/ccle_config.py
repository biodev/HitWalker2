import core
import custom_functions
import string

##globals

prog_type = ""

graph_struct_file = "/var/www/hitwalker_2_inst"+core.fix_prog_type(prog_type)+"/graph_struct.json"

cypher_session="http://localhost:7474"
#the maximum number of nodes in a group before it becomes a metanode
max_nodes = 1

#sizes for network view
network_sizes = {"w":400, "h":400, "legend_offset":200, "history_offset":10}

#sizes for pathway view
pathway_sizes = {"w":800, "h":800, "legend_offset":200, "history_offset":10}

#the seeds need to refer to all the variables that are used as 'seeds' for the ranking algorithm
#the target is the variable to be ranked (otherwise known as the query)
data_types={
    'seeds':['DrugScore'],
    'target':'Variants'
}

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
    'Seed_Parameters':{'type':'standard',
                       'fields':{
                                'drug_score':{'type':'numeric', 'default':10, 'range':[-100,100], 'comparison':'>', 'name':'Gene Score'}
                            }
                        },
    'General_Parameters':{'type':'standard',
                          'fields':{
                                'string_conf':{'type':'numeric', 'default':.4, 'range':[0, 1], 'comparison':'>', 'name':'HitWalker STRING Confidence'},
                                'path_conf':{'type':'numeric', 'default':.95, 'range':[0, 1], 'comparison':'>', 'name':'Pathway STRING Confidence'},
                                'res_prob':{'type':'numeric', 'default':.3, 'range':[0,1], 'comparison':'=', 'name':'Restart Probability'},
                                'max_iter':{'type':'numeric', 'default':100, 'range':[0, 10000], 'comparison':'=', 'name':'Max Iterations'},
                                'conv_thresh':{'type':'numeric', 'default':1e-10, 'range':[0,1], 'comparison':'<', 'name':'Convergence Threshold'},
                                'expression':{'type':'numeric', 'default':.75, 'range':[0,1], 'comparison':'>', 'name':'Expression (Hit) Threshold'},
                            }   
                        }
}

#Basic info for genes:
#Needs the query to return info in the following form:
    #gene ID, gene symbol, and the collection of pathways (or [] if that info is not available).
#gene_names = {'query':'MATCH (n:Gene{name:{GENE}})-[r:KNOWN_AS]-(m) WHERE r.status="symbol" WITH n,m OPTIONAL MATCH (n)<-[:EXTERNAL_ID]-()<-[:GENESET_CONTAINS]-(p) RETURN n.name,m.name, COLLECT(DISTINCT p.name)', 'handler':custom_functions.get_gene_names, 'session_params':None}

#gene_names = {'query':'MATCH (n:Gene{name:{GENE}})-[r:KNOWN_AS]-(m) WHERE r.status="symbol" RETURN n.name,m.name, []', 'handler':custom_functions.get_gene_names, 'session_params':None}


#Basic node info for the seeds:

#needs to return data in the following form:
    #gene_name, sample_name, var (corresponding to datatypes), and other metadata...

#Also note here that perhaps the weight type variable in r2 should be used...

#as our test case: 22RV1_PROSTATE

gene_score = {'query':'MATCH (n:Sample)-[r:HAS_DRUG_ASSAY]-() WITH percentileCont(r.ic50, .5) AS med_ic50 MATCH (n:Sample)-[r:HAS_DRUG_ASSAY]-(m)-[r2:ACTS_ON]-(o:EntrezID{name:{GENE}}) WHERE n.name IN {SAMPLE}  \
                        WITH n, o,  SUM(CASE WHEN r.ic50 < (med_ic50 / 5.0) THEN 1 ELSE -1 END) AS effect_score  RETURN o.name AS gene, n.name AS sample, "DrugScore" AS var, effect_score AS score, effect_score > 0 AS is_hit',
              'handler':None, 'session_params':None}

variant = {'query': 'MATCH (n:Sample)-[r:HAS_DNASEQ]-(var)-[r2:IMPACTS]-(gene) WHERE n.name IN {SAMPLE} RETURN n,r,m,r2,o' , 'handler':None, 'session_params':None}

#gene_score = {'query':'MATCH(n:LabID)-[r:GENE_SCORE_RUN]-()-[r2:SCORE_MAPPED_TO]-(m:Gene{name:{GENE}}) WHERE HAS(r.score) AND n.name IN {LABID} RETURN m.name AS gene, n.name AS sample, "GeneScore" AS var, MAX(r.score) AS score ,ANY(x IN COLLECT(r.score*r2.modifier) WHERE x > {GENESCORE}) AS is_hit',
#              'handler':custom_functions.get_gene_score, 'session_params':None}
#siRNA =  {'query':'MATCH(n:LabID)-[r:SIRNA_RUN]-()-[r2:SIRNA_MAPPED_TO]-(m:Gene{name:{GENE}}) WHERE HAS(r.zscore) AND n.name IN {LABID} RETURN m.name AS gene, n.name AS sample, "siRNA" AS var, r.run_type AS type, MIN(r.zscore) AS score, ANY(x IN COLLECT(r.zscore*r2.modifier) WHERE x < {SIRNASCORE}) AS is_hit',
#          'handler':custom_functions.get_sirna_score, 'session_params':None}
#
#gene_exprs_low = {'query':'MATCH(n:LabID)-[r:EXON_ARRAY_RUN]-()-[r2:ASSIGNED]-(m:Gene{name:{GENE}}) WHERE HAS(r.score) AND n.name IN {LABID} RETURN m.name AS gene, n.name AS sample, "LowExpr" AS var, MAX(r.score) AS score, ALL(x IN COLLECT(r.score) WHERE x < {LOWEXPRS}) AS is_hit',
#          'handler':custom_functions.get_exprs, 'session_params':None}
#
##where the distinct in the collect function addresses cases where a sample has been run multiple times...
#siRNA_node_group_tmpl = {'title':'$$ret_type$$s with siRNA hits for $$result$$','text':'siRNA Hit', 'query':'MATCH(sample:LabID)-[r:SIRNA_RUN]-()-[r2:SIRNA_MAPPED_TO]-(gene:Gene) WHERE  HAS(r.zscore) AND (r.zscore*r2.modifier) < {zscore} AND $$lower_coll_type$$.name IN {$$coll_type$$} WITH $$lower_ret_type$$.name AS ret_type, COLLECT(DISTINCT $$lower_coll_type$$.name) AS use_coll WHERE LENGTH(use_coll) = {$$coll_type$$_length} RETURN ret_type',
#                         'handler':None, 'session_params':[['zscore']]}
#
#genescore_node_group_tmpl = {'title':'$$ret_type$$s with GeneScore hits for $$result$$','text':'GeneScore Hit', 'query':'MATCH(sample:LabID)-[r:GENE_SCORE_RUN]-()-[r2:SCORE_MAPPED_TO]-(gene:Gene) WHERE HAS(r.score) AND (r.score*r2.modifier) > {gene_score} AND $$lower_coll_type$$.name IN {$$coll_type$$} WITH $$lower_ret_type$$.name AS ret_type, COLLECT(DISTINCT $$lower_coll_type$$.name) AS use_coll WHERE LENGTH(use_coll) = {$$coll_type$$_length} RETURN ret_type',
#                            'handler':None, 'session_params':[['gene_score']]}
#
#variant = {'query':'MATCH (gene:Gene{name:{GENE}})-[:TRANSCRIBED]->(trans)<-[impacts_r:IMPACTS]-(var) WITH gene,trans,var,impacts_r MATCH (labid:LabID)<-[alias_r:ALIAS_OF]-(samp)-[dna_diff_r:DNA_DIFF]->(var)<-[unit_dna_diff_r:UNIT_DNA_DIFF]-(exp) WHERE (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype" AND labid.name IN {LABID} RETURN DISTINCT gene.name, var.name, labid.name',
#           'handler':custom_functions.get_variants, 'session_params':None}
#
#gene_rels = {'query':'MATCH (gene:Gene{name:{FROM_GENE}})-[:EXTERNAL_ID]-()-[:MAPPED_TO]-(string_from)-[r:ASSOC]-(string_to)-[:MAPPED_TO]-()-[:EXTERNAL_ID]-(gene_to) WHERE gene_to.name IN {TO_GENES} AND HAS(r.score) AND r.score > ({string_conf}*1000) RETURN gene.name,gene_to.name, MAX(r.score)', 'handler':None, 'session_params':[['string_conf']]}
#
#labid = {'query':'MATCH (n:LabID{name:{LABID}})<-[:PRODUCED]-(m)<-[:HAS_DISEASE]-(k) RETURN n,m,k', 'handler':custom_functions.get_lab_id, 'session_params':None}
#
#pathway = {'query':'MATCH (path:GeneSet{name:{PATHNAME}})-[:GENESET_CONTAINS]->()-[:EXTERNAL_ID]->(gene) RETURN path.name, COLLECT(DISTINCT gene.name)', 'handler':custom_functions.get_pathway, 'session_params':None}
#

#sample relationships

#

#needs to have $$sample$$ which will be supplied by the 'text' return value of matchers['sample']
sample_rels_query = 'MATCH (cellline:CellLine)-[]->(sample) WHERE ANY(x IN [cellline.name] + cellline.alias WHERE x = "$$sample$$") WITH cellline, sample ' + \
                    'OPTIONAL MATCH (sample)-[:HAS_EXPRESSION]-(expr) WITH cellline, sample, COUNT(expr) AS Exprs ' + \
                    'OPTIONAL MATCH (sample)-[:HAS_DRUG_ASSAY]-(assay) WITH cellline, sample, Exprs, COUNT(assay) AS DrugScore ' + \
                    'OPTIONAL MATCH (sample)-[:HAS_DNASEQ]-(variant) WITH cellline, sample, Exprs, DrugScore, COUNT(variant) AS Variants ' + \
                    ' WHERE Exprs > 0 OR DrugScore > 0 OR Variants > 0 RETURN cellline.name as CellLine, cellline.histology_subtype AS Type, sample.name AS Sample, Exprs, DrugScore, Variants, CASE WHEN (DrugScore > 0) AND (Variants > 0) THEN 1 ELSE 0 END AS required_data'

sample_rels_type = 'hierarchical'

##This specifies the functions to be used for searching on genes, pathways and samples
matchers = {
    #'pathway':custom_functions.match_pathway,
    #'gene':custom_functions.match_gene,
    'sample':custom_functions.match_sample
}

##The queries used to choose the data used as seeds or otherwise in the prioritization algorithm
#core.handle_hits needs the following:

#{name} is necessary here as the sample name will come directly from getNodes
hit_session_dict = {"DrugScore":[core.customize_query(gene_score, query=lambda x:x.replace("{name:{GENE}}", "").replace("{SAMPLE}", "{name}") , handler=lambda x: core.handle_hits)]}

#Here we also need to have 'query_ind', 'gene_ind' and a unique 'row_id'
query_prior_dict = {"Variants":[
                            {'query': ('MATCH (n:Sample)-[r:HAS_DNASEQ]-(var)-[r2:IMPACTS]-(gene) WHERE n.name IN {name}' 
                              'RETURN var.name AS Variant_Position, r2.transcript AS Transcript, gene.name AS Gene,'
                              'r.ref_counts as Ref_Counts, r.alt_counts AS Alt_Counts, REPLACE(RTRIM(REDUCE(str="",n IN var.dbsnp|str+n+" ")), " ", ";") AS dbSNP,'
                              'r2.variant_classification AS variant_classification, r2.protein AS Protein_Change, 0 AS query_ind, 1 AS gene_ind, var.name + "_" + gene.name AS row_id') ,
                            'handler':core.handle_query_prior, 'session_params':None}]}

score_hits=core.no_combining
#convert_ids_to=custom_functions.gene_seed_list_to_protein
#convert_ids_from=custom_functions.protein_seed_list_to_gene
#get_seed_list=custom_functions.make_seed_list
#get_query_list = custom_functions.make_seed_list
#
##Need to ensure that the query can find the appropriate parameters in the session and that the sample is refered to as '{name}'
prioritization_func={'function':custom_functions.netprop_rwr, "args":{"initial_graph_file":"/var/www/hitwalker_2_inst/protein.links.detailed.v9.05.9606.mm.mtx",
                                                                      "string_session_name":"string_conf",
                                                                      "res_prob_session_name":"res_prob",
                                                                      "conv_thresh_session_name":"conv_thresh",
                                                                      "max_iter_session_name":"max_iter"
                                                                        }}



#['query_samples', 'LabID', 'GeneScore']
#['query_samples', 'LabID', 'siRNA']
#lambda x: [['query_samples', 'LabID', 'Variants']]

#node_queries = {
#    'Gene':[core.customize_query(gene_names, query=lambda x: x.replace("{GENE}", "{name}")),
#            core.customize_query(gene_score, query=lambda x: x.replace("{GENE}", "{name}").replace("{GENESCORE}", "{gene_score}").replace("{LABID}", "{GeneScore}"), session_params=lambda x: [['gene_score']]),
#            core.customize_query(siRNA, query=lambda x:x.replace("{LABID}", "{siRNA}").replace("{GENE}", "{name}").replace("{SIRNASCORE}", "{zscore}"), session_params=lambda x:[['zscore']]),
#            core.customize_query(variant, query=lambda x:x.replace("{GENE}", "{name}").replace("{LABID}","{Variants}"), session_params=lambda x:None),
#            #core.customize_query(gene_exprs_low, query=lambda x:x.replace("{GENE}", "{name}").replace("{LABID}","{Variants}").replace("{LOWEXPRS}", "{low_exprs}"), session_params=lambda x:[['low_exprs']]),
#            core.customize_query(gene_exprs_low, query=lambda x:x.replace("{GENE}", "{name}").replace("{LABID}","{Variants}").replace("LowExpr", "HighExpr").replace("< {LOWEXPRS}", "> {high_exprs}"), session_params=lambda x:[['high_exprs']])
#            ],
#    'LabID':[core.customize_query(labid, query=lambda x: x.replace("{LABID}", "{name}"))]
#            
#}
#
##note that Sample is used in the HitWalker2 javascript instead of LabID so it is synonymous with LabID here...
#edge_queries = {
#    
#    'nodes':{
#        
#        'Gene':[core.customize_query(gene_names, query=lambda x: x.replace("{GENE}", "{name}")),
#                core.customize_query(gene_score, query=lambda x: x.replace("{GENE}", "{name}").replace("{GENESCORE}", "{gene_score}"), session_params=lambda x: [['gene_score']]),
#                core.customize_query(siRNA, query=lambda x:x.replace("{GENE}", "{name}").replace("{SIRNASCORE}", "{zscore}"), session_params=lambda x:[['zscore']]),
#                core.customize_query(variant, query=lambda x:x.replace("{GENE}", "{name}")),
#                #core.customize_query(gene_exprs_low, query=lambda x:x.replace("{GENE}", "{name}").replace("{LOWEXPRS}", "{low_exprs}"), session_params=lambda x:[['low_exprs']]),
#                core.customize_query(gene_exprs_low, query=lambda x:x.replace("{GENE}", "{name}").replace("LowExpr", "HighExpr").replace("< {LOWEXPRS}", "> {high_exprs}"), session_params=lambda x:[['high_exprs']])],
#        'Sample':[core.customize_query(labid, query=lambda x: x.replace("{LABID}", "{name}"))],
#        'Pathway':[core.customize_query(pathway, query=lambda x: x.replace("{PATHNAME}", "{name}"))]
#        },
#    
#    'edges':{
#        'Gene':{
#            'Sample':{'query':None, 'handler':custom_functions.gene_to_lab_id, 'session_params':None},
#            'Gene':{'query':None, 'handler':custom_functions.no_rels, 'session_params':None} #custom_functions.gene_to_gene
#        },
#        'Sample':{
#            'Sample':{'query':None, 'handler':custom_functions.sample_to_sample, 'session_params':None},
#        }
#    }
#}
#
##This dictionary specifies how the summary modals should be displayed upon right-clicking a specified node
#node_content = {
#    'Sample':{'func':core.make_sample_table, 'args':{'sample_link':'https://octri.ohsu.edu/lls_scor/datamgmt/Patient/Specimen?search='}},
#    'Gene':{'func':core.make_gene_table, 'args':{'gene_link':'http://jan2013.archive.ensembl.org/Homo_sapiens/Gene/Summary?g=', 'plot_vars':{'siRNA':'make_waterfall_plot(this)'}}},
#    'MetaNode':{'func':core.make_metanode_table, 'args':{}}
#}
##{'title':'Samples with siRNA hits for $$result$$','text':'siRNA Hit', 'query':'MATCH(n:LabID)-[r:SIRNA_RUN]-()-[r2:SIRNA_MAPPED_TO]-(m:Gene) WHERE (n)<-[:PRODUCED]-()<-[:HAS_DISEASE]-() AND HAS(r.zscore) AND (r.zscore*r2.modifier) < {zscore} AND m.name IN {Gene} WITH n.name AS sample, COLLECT(m.name) AS gene_col WHERE LENGTH(gene_col) = {Gene_length}  RETURN sample', 'handler':custom_functions.return_sample_nodes, 'session_params':[['zscore']]},
#node_group_content = {
#    'Gene':{'title':'Find Samples where all Genes in this set has a(n):',
#            'returned_node_type':'LabID',
#            'options':[
#            {'title':'Samples with Variants for $$result$$', 'text':'Variant', 'query':'MATCH (gene:Gene)-[:TRANSCRIBED]->(trans)<-[impacts_r:IMPACTS]-(var) WHERE gene.name IN {Gene} WITH gene,trans,var,impacts_r \
#                        MATCH (sample:LabID)<-[alias_r:ALIAS_OF]-(samp)-[dna_diff_r:DNA_DIFF]->(var)<-[unit_dna_diff_r:UNIT_DNA_DIFF]-(exp) WHERE (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype" \
#                        WITH sample.name AS ret_type, COLLECT(DISTINCT gene.name) AS use_coll WHERE LENGTH(use_coll) = {Gene_length} RETURN ret_type', 'handler':None, 'session_params':None},
#            core.specify_type_query_tmp(siRNA_node_group_tmpl, ret_type='Sample', coll_type='Gene'),
#            core.specify_type_query_tmp(genescore_node_group_tmpl, ret_type='Sample', coll_type='Gene')]},
#    'Sample':{'title':'Find Genes where all Samples in this set has a(n):',
#        'returned_node_type':'Gene',
#        'options':[
#       
#        core.specify_type_query_tmp(siRNA_node_group_tmpl, ret_type='Gene', coll_type='Sample'),
#        core.specify_type_query_tmp(genescore_node_group_tmpl, ret_type='Gene', coll_type='Sample')]
#        }
#    }
#
## Need to work on this for Sample... {'text':'Variant', 'query':'', 'handler':'', 'session_params':None},
#
#node_abbreviations = {
#    'MYELOPROLIFERATIVE NEOPLASMS':'MPN',
#    'PRECURSOR LYMPHOID NEOPLASMS':'PLN',
#    'ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS': 'AML',
#    'MYELODYSPLASTIC SYNDROMES':'MDS',
#    'MYELODYSPLASTIC/MYELOPROLIFERATIVE NEOPLASMS':'MDS/MPN',
#    'MATURE B-CELL NEOPLASMS':'MBCN',
#    'HighExpr':'Expression',
#    'HighExpr_Hit':'HighExpression'
#}
#
#edge_abbreviations = {
#    'Possible_HighExpr':'Expression_Collected',
#    'Observed_HighExpr':'Observed_HighExpression',
#    
#}
#

#
#plot_data_types = {
#    'siRNA':custom_functions.get_sirna_data
#}
#
#graph_initializers = {
#    'panel':custom_functions.get_shortest_paths,
#    'image':custom_functions.get_pathways_sample
#}
#
