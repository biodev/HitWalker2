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
    'seeds':['siRNA', 'GeneScore'],
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
                                'zscore':{'type':'numeric', 'default':-2,'range':[-100, 0], 'comparison':'<', 'name':'siRNA Zscore'},
                                'gene_score':{'type':'numeric', 'default':10, 'range':[-100,100], 'comparison':'>', 'name':'Gene Score'}
                            }
                        },
    'General_Parameters':{'type':'standard',
                          'fields':{
                                'string_conf':{'type':'numeric', 'default':.4, 'range':[0, 1], 'comparison':'>', 'name':'HitWalker STRING Confidence'},
                                'path_conf':{'type':'numeric', 'default':.95, 'range':[0, 1], 'comparison':'>', 'name':'Pathway STRING Confidence'},
                                'res_prob':{'type':'numeric', 'default':.3, 'range':[0,1], 'comparison':'=', 'name':'Restart Probability'},
                                'max_iter':{'type':'numeric', 'default':100, 'range':[0, 10000], 'comparison':'=', 'name':'Max Iterations'},
                                'conv_thresh':{'type':'numeric', 'default':1e-10, 'range':[0,1], 'comparison':'<', 'name':'Convergence Threshold'},
                                'low_exprs':{'type':'numeric', 'default':.25, 'range':[0,1], 'comparison':'<', 'name':'Low Expression Threshold'},
                                'high_exprs':{'type':'numeric', 'default':.75, 'range':[0,1], 'comparison':'>', 'name':'High Expression Threshold'},
                            }   
                        },
    'Variant_Filters':{'type':'grouped',
                       'fields':{
                            'freq':{'type':'numeric', 'comparison':'<','default':.01, 'range':[0,1], 'name':'Global MAF', 'var_name':'gmaf','required':{'from':'Variation'}, 'trans':core.return_numeric, 'needs_has':''},
                            'cohort_freq':{'type':'numeric', 'comparison':'<', 'default':.5, 'range':[0,1], 'name':'Cohort Alt. Frequency', 'var_name':'Sample_count', 'required':{'from':'Variation'}, 'trans':custom_functions.make_freq_from_count},
                            'cohort_count':{'type':'numeric', 'comparison':'=', 'default':1, 'range':[0,500], 'name':'Cohort Count', 'var_name':'Sample_count', 'required':{'from':'Variation'}, 'trans':core.return_numeric},
                            'in_1kg':{'type':'character', 'default':'False', 'range':['True', 'False'], 'name':'In 1000 genomes','var_name':'in_1kg','required':{'from':'Variation'}, 'trans':core.return_binary},
                            'in_dbsnp':{'type':'character', 'default':'False', 'range':['True', 'False'], 'name':'In dbSNP','var_name':'in_dbsnp', 'required':{'from':'Variation'},'trans':core.return_binary},
                            'allele_count': {'type':'numeric', 'range':[0,2], 'name':'Genotype','var_name':'allele_count','required':{'from':'DNA_DIFF'},'trans':core.return_numeric},
                            'genotype_quality': {'type':'numeric', 'comparison':'>', 'default':40, 'range':[0,100], 'name':'Genotype Quality', 'var_name':'genotype_quality', 'required':{'from':'DNA_DIFF'},'trans':core.return_numeric},
                            'depth':{'type':'numeric', 'range':[0,100000], 'name':'Read Depth','var_name':'depth', 'required':{'from':'DNA_DIFF'},'trans':core.return_numeric},
                            
                            'FS': {'type':'numeric', 'range':[-10000,10000],'name':'Fisher Strand','var_name':'FS', 'required':{'from':'UNIT_DNA_DIFF'},'trans':core.return_numeric},
                            'MQ0':{'type':'numeric', 'comparison':'<', 'default':4, 'range':[0,10000], 'name':'Number Ambigous Reads','var_name':'MQ0', 'required':{'from':'UNIT_DNA_DIFF'},'trans':core.return_numeric},
                            'MQ': {'type':'numeric', 'range':[0,100], 'name':'Mapping Quality','var_name':'MQ', 'required':{'from':'UNIT_DNA_DIFF'},'trans':core.return_numeric},
                            'QD': {'type':'numeric', 'comparison':'>', 'default':5, 'range':[0,10000], 'name':'Quality / Depth','var_name':'QD', 'required':{'from':'UNIT_DNA_DIFF'},'trans':core.return_numeric},
                            'SB': {'type':'numeric', 'range':[-10000,10000], 'name':'Strand Bias','var_name':'SB', 'required':{'from':'UNIT_DNA_DIFFERENCE'},'trans':core.return_numeric},
                            
                            'Cons_cat':{'type':'character', 'default':'NonSynon.', 'range':['Synonymous', 'NonSynon.', 'Other'], 'var_name':'Cons_cat', 'name':'Consequence', 'required':{'from':'IMPACTS'}, 'trans':custom_functions.return_full_cons}
                    },
                    'default_groups':[[
                        [{'field':'Cons_cat'},
                         {'field':'genotype_quality', 'logical':'AND'},
                         {'field':'QD', 'logical':'AND'},
                         {'field':'MQ0', 'logical':'AND'},
                         {'field':'cohort_freq','logical':'AND'}]
                        ],
                        [
                            [{'field':'in_1kg', 'default':'False', 'logical':'AND'},{'field':'in_dbsnp','default':'False','logical':'AND'}],
                            [{'field':'in_1kg', 'default':'True', 'logical':'OR'}, {'field':'freq', 'logical':'AND'}],
                            [{'field':'in_dbsnp','default':'True', 'logical':'OR'},{'field':'in_1kg','default':'False', 'logical':'AND'},{'field':'cohort_count','logical':'AND'}]
                        ]]
                    
    }
}

#basic set of queries, handlers etc.

#Basic info for genes:
#Needs the query to return info in the following form:
    #gene ID, gene symbol, and the collection of pathways (or [] if that info is not available).
#gene_names = {'query':'MATCH (n:Gene{name:{GENE}})-[r:KNOWN_AS]-(m) WHERE r.status="symbol" WITH n,m OPTIONAL MATCH (n)<-[:EXTERNAL_ID]-()<-[:GENESET_CONTAINS]-(p) RETURN n.name,m.name, COLLECT(DISTINCT p.name)', 'handler':custom_functions.get_gene_names, 'session_params':None}

gene_names = {'query':'MATCH (n:Gene{name:{GENE}})-[r:KNOWN_AS]-(m) WHERE r.status="symbol" RETURN n.name,m.name, []', 'handler':custom_functions.get_gene_names, 'session_params':None}


#Basic node info for the seeds:

#needs to return data in the following form:
    #gene_name, sample_name, var (corresponding to datatypes), and other metadata...

gene_score = {'query':'MATCH(n:LabID)-[r:GENE_SCORE_RUN]-()-[r2:SCORE_MAPPED_TO]-(m:Gene{name:{GENE}}) WHERE HAS(r.score) AND n.name IN {LABID} RETURN m.name AS gene, n.name AS sample, "GeneScore" AS var, MAX(r.score) AS score ,ANY(x IN COLLECT(r.score*r2.modifier) WHERE x > {GENESCORE}) AS is_hit',
              'handler':custom_functions.get_gene_score, 'session_params':None}
siRNA =  {'query':'MATCH(n:LabID)-[r:SIRNA_RUN]-()-[r2:SIRNA_MAPPED_TO]-(m:Gene{name:{GENE}}) WHERE HAS(r.zscore) AND n.name IN {LABID} RETURN m.name AS gene, n.name AS sample, "siRNA" AS var, r.run_type AS type, MIN(r.zscore) AS score, ANY(x IN COLLECT(r.zscore*r2.modifier) WHERE x < {SIRNASCORE}) AS is_hit',
          'handler':custom_functions.get_sirna_score, 'session_params':None}

gene_exprs_low = {'query':'MATCH(n:LabID)-[r:EXON_ARRAY_RUN]-()-[r2:ASSIGNED]-(m:Gene{name:{GENE}}) WHERE HAS(r.score) AND n.name IN {LABID} RETURN m.name AS gene, n.name AS sample, "LowExpr" AS var, MAX(r.score) AS score, ALL(x IN COLLECT(r.score) WHERE x < {LOWEXPRS}) AS is_hit',
          'handler':custom_functions.get_exprs, 'session_params':None}

#where the distinct in the collect function addresses cases where a sample has been run multiple times...
siRNA_node_group_tmpl = {'title':'$$ret_type$$s with siRNA hits for $$result$$','text':'siRNA Hit', 'query':'MATCH(sample:LabID)-[r:SIRNA_RUN]-()-[r2:SIRNA_MAPPED_TO]-(gene:Gene) WHERE  HAS(r.zscore) AND (r.zscore*r2.modifier) < {zscore} AND $$lower_coll_type$$.name IN {$$coll_type$$} WITH $$lower_ret_type$$.name AS ret_type, COLLECT(DISTINCT $$lower_coll_type$$.name) AS use_coll WHERE LENGTH(use_coll) = {$$coll_type$$_length} RETURN ret_type',
                         'handler':None, 'session_params':[['zscore']]}

genescore_node_group_tmpl = {'title':'$$ret_type$$s with GeneScore hits for $$result$$','text':'GeneScore Hit', 'query':'MATCH(sample:LabID)-[r:GENE_SCORE_RUN]-()-[r2:SCORE_MAPPED_TO]-(gene:Gene) WHERE HAS(r.score) AND (r.score*r2.modifier) > {gene_score} AND $$lower_coll_type$$.name IN {$$coll_type$$} WITH $$lower_ret_type$$.name AS ret_type, COLLECT(DISTINCT $$lower_coll_type$$.name) AS use_coll WHERE LENGTH(use_coll) = {$$coll_type$$_length} RETURN ret_type',
                            'handler':None, 'session_params':[['gene_score']]}

variant = {'query':'MATCH (gene:Gene{name:{GENE}})-[:TRANSCRIBED]->(trans)<-[impacts_r:IMPACTS]-(var) WITH gene,trans,var,impacts_r MATCH (labid:LabID)<-[alias_r:ALIAS_OF]-(samp)-[dna_diff_r:DNA_DIFF]->(var)<-[unit_dna_diff_r:UNIT_DNA_DIFF]-(exp) WHERE (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype" AND labid.name IN {LABID} RETURN DISTINCT gene.name, var.name, labid.name',
           'handler':custom_functions.get_variants, 'session_params':None}

gene_rels = {'query':'MATCH (gene:Gene{name:{FROM_GENE}})-[:EXTERNAL_ID]-()-[:MAPPED_TO]-(string_from)-[r:ASSOC]-(string_to)-[:MAPPED_TO]-()-[:EXTERNAL_ID]-(gene_to) WHERE gene_to.name IN {TO_GENES} AND HAS(r.score) AND r.score > ({string_conf}*1000) RETURN gene.name,gene_to.name, MAX(r.score)', 'handler':None, 'session_params':[['string_conf']]}

labid = {'query':'MATCH (n:LabID{name:{LABID}})<-[:PRODUCED]-(m)<-[:HAS_DISEASE]-(k) RETURN n,m,k', 'handler':custom_functions.get_lab_id, 'session_params':None}

pathway = {'query':'MATCH (path:GeneSet{name:{PATHNAME}})-[:GENESET_CONTAINS]->()-[:EXTERNAL_ID]->(gene) RETURN path.name, COLLECT(DISTINCT gene.name)', 'handler':custom_functions.get_pathway, 'session_params':None}

#['query_samples', 'LabID', 'GeneScore']
#['query_samples', 'LabID', 'siRNA']
#lambda x: [['query_samples', 'LabID', 'Variants']]

node_queries = {
    'Gene':[core.customize_query(gene_names, query=lambda x: x.replace("{GENE}", "{name}")),
            core.customize_query(gene_score, query=lambda x: x.replace("{GENE}", "{name}").replace("{GENESCORE}", "{gene_score}").replace("{LABID}", "{GeneScore}"), session_params=lambda x: [['gene_score']]),
            core.customize_query(siRNA, query=lambda x:x.replace("{LABID}", "{siRNA}").replace("{GENE}", "{name}").replace("{SIRNASCORE}", "{zscore}"), session_params=lambda x:[['zscore']]),
            core.customize_query(variant, query=lambda x:x.replace("{GENE}", "{name}").replace("{LABID}","{Variants}"), session_params=lambda x:None),
            #core.customize_query(gene_exprs_low, query=lambda x:x.replace("{GENE}", "{name}").replace("{LABID}","{Variants}").replace("{LOWEXPRS}", "{low_exprs}"), session_params=lambda x:[['low_exprs']]),
            core.customize_query(gene_exprs_low, query=lambda x:x.replace("{GENE}", "{name}").replace("{LABID}","{Variants}").replace("LowExpr", "HighExpr").replace("< {LOWEXPRS}", "> {high_exprs}"), session_params=lambda x:[['high_exprs']])
            ],
    'LabID':[core.customize_query(labid, query=lambda x: x.replace("{LABID}", "{name}"))]
            
}

#note that Sample is used in the HitWalker2 javascript instead of LabID so it is synonymous with LabID here...
edge_queries = {
    
    'nodes':{
        
        'Gene':[core.customize_query(gene_names, query=lambda x: x.replace("{GENE}", "{name}")),
                core.customize_query(gene_score, query=lambda x: x.replace("{GENE}", "{name}").replace("{GENESCORE}", "{gene_score}"), session_params=lambda x: [['gene_score']]),
                core.customize_query(siRNA, query=lambda x:x.replace("{GENE}", "{name}").replace("{SIRNASCORE}", "{zscore}"), session_params=lambda x:[['zscore']]),
                core.customize_query(variant, query=lambda x:x.replace("{GENE}", "{name}")),
                #core.customize_query(gene_exprs_low, query=lambda x:x.replace("{GENE}", "{name}").replace("{LOWEXPRS}", "{low_exprs}"), session_params=lambda x:[['low_exprs']]),
                core.customize_query(gene_exprs_low, query=lambda x:x.replace("{GENE}", "{name}").replace("LowExpr", "HighExpr").replace("< {LOWEXPRS}", "> {high_exprs}"), session_params=lambda x:[['high_exprs']])],
        'Sample':[core.customize_query(labid, query=lambda x: x.replace("{LABID}", "{name}"))],
        'Pathway':[core.customize_query(pathway, query=lambda x: x.replace("{PATHNAME}", "{name}"))]
        },
    
    'edges':{
        'Gene':{
            'Sample':{'query':None, 'handler':custom_functions.gene_to_lab_id, 'session_params':None},
            'Gene':{'query':None, 'handler':custom_functions.no_rels, 'session_params':None} #custom_functions.gene_to_gene
        },
        'Sample':{
            'Sample':{'query':None, 'handler':custom_functions.sample_to_sample, 'session_params':None},
        }
    }
}

#This dictionary specifies how the summary modals should be displayed upon right-clicking a specified node
node_content = {
    'Sample':{'func':core.make_sample_table, 'args':{'sample_link':'https://octri.ohsu.edu/lls_scor/datamgmt/Patient/Specimen?search='}},
    'Gene':{'func':core.make_gene_table, 'args':{'gene_link':'http://jan2013.archive.ensembl.org/Homo_sapiens/Gene/Summary?g=', 'plot_vars':{'siRNA':'make_waterfall_plot(this)'}}},
    'MetaNode':{'func':core.make_metanode_table, 'args':{}}
}
#{'title':'Samples with siRNA hits for $$result$$','text':'siRNA Hit', 'query':'MATCH(n:LabID)-[r:SIRNA_RUN]-()-[r2:SIRNA_MAPPED_TO]-(m:Gene) WHERE (n)<-[:PRODUCED]-()<-[:HAS_DISEASE]-() AND HAS(r.zscore) AND (r.zscore*r2.modifier) < {zscore} AND m.name IN {Gene} WITH n.name AS sample, COLLECT(m.name) AS gene_col WHERE LENGTH(gene_col) = {Gene_length}  RETURN sample', 'handler':custom_functions.return_sample_nodes, 'session_params':[['zscore']]},
node_group_content = {
    'Gene':{'title':'Find Samples where all Genes in this set has a(n):',
            'returned_node_type':'LabID',
            'options':[
            {'title':'Samples with Variants for $$result$$', 'text':'Variant', 'query':'MATCH (gene:Gene)-[:TRANSCRIBED]->(trans)<-[impacts_r:IMPACTS]-(var) WHERE gene.name IN {Gene} WITH gene,trans,var,impacts_r \
                        MATCH (sample:LabID)<-[alias_r:ALIAS_OF]-(samp)-[dna_diff_r:DNA_DIFF]->(var)<-[unit_dna_diff_r:UNIT_DNA_DIFF]-(exp) WHERE (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype" \
                        WITH sample.name AS ret_type, COLLECT(DISTINCT gene.name) AS use_coll WHERE LENGTH(use_coll) = {Gene_length} RETURN ret_type', 'handler':None, 'session_params':None},
            core.specify_type_query_tmp(siRNA_node_group_tmpl, ret_type='Sample', coll_type='Gene'),
            core.specify_type_query_tmp(genescore_node_group_tmpl, ret_type='Sample', coll_type='Gene')]},
    'Sample':{'title':'Find Genes where all Samples in this set has a(n):',
        'returned_node_type':'Gene',
        'options':[
       
        core.specify_type_query_tmp(siRNA_node_group_tmpl, ret_type='Gene', coll_type='Sample'),
        core.specify_type_query_tmp(genescore_node_group_tmpl, ret_type='Gene', coll_type='Sample')]
        }
    }

# Need to work on this for Sample... {'text':'Variant', 'query':'', 'handler':'', 'session_params':None},

node_abbreviations = {
    'MYELOPROLIFERATIVE NEOPLASMS':'MPN',
    'PRECURSOR LYMPHOID NEOPLASMS':'PLN',
    'ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS': 'AML',
    'MYELODYSPLASTIC SYNDROMES':'MDS',
    'MYELODYSPLASTIC/MYELOPROLIFERATIVE NEOPLASMS':'MDS/MPN',
    'MATURE B-CELL NEOPLASMS':'MBCN',
    'HighExpr':'Expression',
    'HighExpr_Hit':'HighExpression'
}

edge_abbreviations = {
    'Possible_HighExpr':'Expression_Collected',
    'Observed_HighExpr':'Observed_HighExpression',
    
}

#This specifies the functions to be used for searching on genes, pathways and samples
matchers = {
    'pathway':custom_functions.match_pathway,
    'gene':custom_functions.match_gene,
    'sample':custom_functions.match_sample
}

plot_data_types = {
    'siRNA':custom_functions.get_sirna_data
}

graph_initializers = {
    'panel':custom_functions.get_shortest_paths,
    'image':custom_functions.get_pathways_sample
}

#sample relationships

sample_rels_query = 'MATCH (init_id:LabID{name:"$$sample$$"})-[:PRODUCED]-()-[:HAS_DISEASE]-(patient) WITH patient MATCH (patient)-[:HAS_DISEASE]-(disease)-[:PRODUCED]-(specimen) WITH patient, disease, specimen '+ \
    'OPTIONAL MATCH (specimen)-[ga:ALIAS_OF]-(geno_sample)-[var:DNA_DIFF]-() WHERE ga.alias_type = "genotype" WITH patient, disease, specimen, COUNT(var) AS Variants ' + \
    'OPTIONAL MATCH (specimen)-[sr:SIRNA_RUN]-() WITH patient, disease, specimen, Variants, COUNT(sr) AS siRNA ' + \
    'OPTIONAL MATCH (specimen)-[gr:GENE_SCORE_RUN]-() RETURN patient.name AS PatientID, disease.specific_diagnosis AS SpecificDiagnosis, specimen.name AS Sample, Variants,siRNA, COUNT(gr) AS GeneScore, CASE WHEN (Variants > 0) AND (siRNA > 0 OR COUNT(gr) > 0) THEN 1 ELSE 0 END AS required_data'
    
sample_rels_type = 'hierarchical'

#The queries used to choose the data used as seeds or otherwise in the prioritization algorithm
hit_session_dict = {"siRNA":[core.customize_query(siRNA, query=lambda x:x.replace("{name:{GENE}}", "").replace("{SIRNASCORE}", "{zscore}").replace("{LABID}", "{name}"), session_params=lambda x:[['zscore'], ], handler=lambda x: core.handle_hits)],
                    "GeneScore":[core.customize_query(gene_score, query=lambda x: x.replace("{name:{GENE}}", "").replace("{GENESCORE}", "{gene_score}").replace("{LABID}", "{name}"), session_params=lambda x: [['gene_score']], handler=lambda x: core.handle_hits)]}

##specifying a header seperately is optional, just for clarity...

variant_table = 'MATCH (labid:LabID)-[alias_r:ALIAS_OF]-(samp)-[dna_diff_r:DNA_DIFF]-(var)-[unit_dna_diff_r:UNIT_DNA_DIFF]-(exp) WHERE (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype" AND labid.name IN {name} WITH samp,dna_diff_r,unit_dna_diff_r,var MATCH (var)-[impacts_r:IMPACTS]->(trans)-[:TRANSCRIBED]-(gene)-[symb_r:KNOWN_AS]-(symbol) WHERE symb_r.status = "symbol"  RETURN '


#Note that the REPLACE(RTRIM(REDUCE part of cypher_header simply provides functionality similar to joinfields and such.
query_table_header = ['var.name AS Variant_position', 'var.var_type AS var_type', 'ROUND(dna_diff_r.allele_count) AS allele_count', 'gene.name AS Gene_name', 'symbol.name AS Symbol', 'trans.name AS Transcript_name','impacts_r.Amino_acids AS Amino_acids',
                 'impacts_r.Protein_position AS Protein_position', 'impacts_r.SIFT AS SIFT', 'impacts_r.PolyPhen AS PolyPhen', 'REPLACE(RTRIM(REDUCE(str="",n IN impacts_r.Consequence|str+n+" ")), " ", ";") AS Consequence', 'var.freq AS freq',
                 'REPLACE(RTRIM(REDUCE(str="",n IN var.Existing_variation|str+n+" ")), " ", ";") AS Existing_variation', 'REPLACE(REPLACE(STR(var.in_dbsnp), "1", "True"), "0", "False") AS in_dbsnp', 'REPLACE(REPLACE(STR(var.in_1kg), "1", "True"), "0", "False") AS in_1kg', '0 AS query_ind', '3 AS gene_ind', 'var.name + "_" + trans.name AS row_id']

query_prior_dict = {"Variants":[{'query':variant_table + string.joinfields(query_table_header, ","), 'handler':custom_functions.handle_query_prior, 'session_params':None}]}

score_hits=custom_functions.combine_sirna_gs
convert_ids_to=custom_functions.gene_seed_list_to_protein
convert_ids_from=custom_functions.protein_seed_list_to_gene
get_seed_list=custom_functions.make_seed_list
get_query_list = custom_functions.make_seed_list

#Need to ensure that the query can find the appropriate parameters in the session and that the sample is refered to as '{name}'
prioritization_func={'function':custom_functions.netprop_rwr, "args":{"initial_graph_file":"/var/www/hitwalker_2_inst/protein.links.detailed.v9.05.9606.mm.mtx",
                                                                      "string_session_name":"string_conf",
                                                                      "res_prob_session_name":"res_prob",
                                                                      "conv_thresh_session_name":"conv_thresh",
                                                                      "max_iter_session_name":"max_iter"
                                                                        }}