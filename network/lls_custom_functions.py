from py2neo import neo4j, cypher
import json
import ast
import string
import itertools
import collections
import core
import sys
import copy
from rwr_utils import threshold_graph, compute_net_prop_mat, compute_rwr


def elem_or_na(cur_dict, val):
    if cur_dict.has_key(val):
        return cur_dict[val]
    else:
        return "N/A"

def truncate_var_name(name):
    split_name = name.split('-')
    
    for i in range(1,3):
        if len(split_name[i]) > 3:
            split_name[i] = split_name[i][:3] + '...'
    
    return split_name[0]+':'+split_name[1]+'-'+split_name[2]

def get_sirna_data(request_dict):
    
    from config import cypher_session
    
    session = cypher.Session(cypher_session)
    tx = session.create_transaction()
    
    neo_query = 'MATCH(n:LabID{name:{name}})-[r:SIRNA_RUN]-(m)-[r2:SIRNA_MAPPED_TO]-(k)-[r3:KNOWN_AS]-(o) where r3.status = "symbol" return DISTINCT r.run_type, r.sample_mean,m.name, COLLECT(DISTINCT o.name)'
    plot_title = 'siRNA results for sample ' + request_dict["sample"]
    
    tx.append(neo_query, {"name":request_dict["sample"]})
    
    data_list=[]
    
    for i in core.BasicResultsIterable(tx.execute()):
        for j in i:
            data_list.append({"type":j[0], "value":j[1], "x_name":j[2], "db_name":j[3]})
    
    data_dict = {'data':data_list}
    data_dict.update(request_dict)
    
    return data_dict, plot_title


def make_freq_from_count(value):
    from config import cypher_session
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/')    
    data = neo4j.CypherQuery(graph_db, 'MATCH (n:LabID)-[r:ALIAS_OF]-() WHERE r.alias_type = "genotype" RETURN COUNT(DISTINCT n)').execute().data
    if len(data) > 1:
        raise Exception("Unexpected length of data")
    
    return (float(value)*float(data[0].values[0]))
    
    
def return_full_cons(value):
    
    if isinstance(value, list) == False:
        value = [value]
    
    ret_value = []
    
    for i in value:
        if i == "NonSynon.":
            ret_value.append('NonSynonymous')
        else:
            ret_value.append(str(value))
            
    if len(ret_value) == 1:
        return ret_value[0]
    else:
        return ret_value

def handle_query_prior(res_list, nodes, request):
    
    table_header = list(core.cypherHeader(res_list))
    
    for sample_table in core.BasicResultsIterable(res_list):
        for row in sample_table:
            nodes.add(core.RowNode(list(row), copy.copy(table_header)))
    
    table_header.remove("gene_ind")
    table_header.remove("query_ind")
    table_header.remove("row_id")
    
    nodes.attributes['header'] = table_header

def combine_sirna_gs(request, valid_gene_hits):
    
    #The siRNA hits are to be set to the max of the gene score values
    #   This explicitly deals with overlapping genes
    
    #now create a new object that contains the transformed scores along with the gene->hit metadata
    
    seed_list = core.SeedList(valid_gene_hits)
    
    #first find the max of the gene scores (if any)
    
    gs_nl = copy.deepcopy(valid_gene_hits)
    
    gs_nl.filterByChild(lambda x: x.getAttr(["attributes","node_type"]) == "GeneScore", any)
    
    if len(gs_nl) > 0:
        gene_score_maxs = max(gs_nl.summarizeChildren(lambda x: x.getAttr(["attributes","meta", "score"]), max))
    else:
        gene_score_maxs = 1
    
    #apply the max gene score to the siRNA genes
    
    si_nl = copy.deepcopy(valid_gene_hits)
    
    si_nl.filterByChild(lambda x: x.getAttr(["attributes","node_type"]) == "siRNA", any)
    
    if len(si_nl) > 0:
    
        score_dict = {}
        
        for i in si_nl.ids():
            score_dict[i] = gene_score_maxs
        
        seed_list.adjustScores(score_dict)
    
    return seed_list

def gene_seed_list_to_protein(request, seed_list):
    
    gene_to_prot_query = 'MATCH (gene:Gene{name:{gene_name}})-[:EXTERNAL_ID]-()-[:MAPPED_TO]-(string_from) RETURN gene.name, string_from.name'
    
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
    
    prot_to_gene_query = 'MATCH (prot:StringID{name:{prot_name}})-[:MAPPED_TO]-()-[:EXTERNAL_ID]-(gene) RETURN prot.name, gene.name'
    
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
    

    #then as this is still a coo_matrix, examine the data, row and col attributes to find 
def get_valid_hits_dep (requested_samps, thresh_dict, graph_db):
    
    gene_score = {}
    hit_annot = {}
    
    #print requested_samps
    
    #start with gene score
    
    if (requested_samps['GeneScore'] != None):
    
        gs_query = 'MATCH(n:LabID{name:{requested_samps}})-[r:GENE_SCORE_RUN]-()-[r2:SCORE_MAPPED_TO]-()-[:EXTERNAL_ID]-()-[:MAPPED_TO]-(m) WHERE (HAS(r.score) AND (r.score*r2.modifier) > {gene_score_thresh}) RETURN r.score,m.name'
        data = neo4j.CypherQuery(graph_db,gs_query).execute(**{'requested_samps':requested_samps['GeneScore'], 'gene_score_thresh':thresh_dict['GeneScore']}).data
        
        for i in data:
            use_name = i.values[1]
            if hit_annot.has_key(use_name) == False:
                hit_annot[use_name] = {'GeneScore':i[0]}
            else:
                hit_annot[use_name]['GeneScore'] = max([hit_annot[use_name]['GeneScore'], i[0]])
        
    #sirna query 
    
    sirna = set()
    
    if requested_samps['siRNA'] != None:
    
        hit_query = 'MATCH(n:LabID{name:{requested_samps}})-[r:SIRNA_RUN]-()-[r2:SIRNA_MAPPED_TO]-()-[:EXTERNAL_ID]-()-[:MAPPED_TO]-(m) WHERE (HAS(r.zscore) AND (r.zscore*r2.modifier) < {zthresh}) RETURN r.run_type, r.zscore ,m.name'
        data = neo4j.CypherQuery(graph_db,hit_query).execute(**{'requested_samps':requested_samps['siRNA'], 'zthresh':thresh_dict['siRNA']}).data
        
        #currently set to either tyner or anupriya will trigger a hit so type is currently ignored...
        for i in data:
            use_name = i.values[2]
            if hit_annot.has_key(use_name) == False:
                hit_annot[use_name] = {'siRNA':i[1]}
            elif hit_annot[use_name].has_key('siRNA') == False:
                hit_annot[use_name]['siRNA'] = i[1]
            else:
                hit_annot[use_name]['siRNA'] = min([hit_annot[use_name]['siRNA'], i[1]])
        
    #for hitwalker purposes, summarize the results taking the max gene score if they overlap, 1s and 0s if only siRNA and gene scores if only gene scores...
    
    gene_score = {}
    max_gene_score = 0
    
    for i in hit_annot.values():
        if i.has_key('GeneScore'):
            max_gene_score = max([max_gene_score, i['GeneScore']])
    
    for i in hit_annot.items():
        if i[1].has_key('GeneScore') and i[1].has_key('siRNA'):
            gene_score[i[0]] = max_gene_score
        elif i[1].has_key('GeneScore'):
            gene_score[i[0]] = i[1]['GeneScore']
        else:
            if max_gene_score == 0:
                gene_score[i[0]] = 1
            else:
                gene_score[i[0]] = max_gene_score
        
    return gene_score, hit_annot

def make_gene_id_dict(graph_db):

    gene_query = 'MATCH (n:Gene)-[:EXTERNAL_ID]-()-[:MAPPED_TO]-(m) RETURN n.name,m.name'
    
    data = neo4j.CypherQuery(graph_db,gene_query).execute().data

    gene_id_dict = {}

    for i in data:
        if gene_id_dict.has_key(i[1]):
            gene_id_dict[i[1]].append(i[0])
        else:
            gene_id_dict[i[1]] = [i[0]]

    return gene_id_dict

def seed_to_gene_ids(seed_dict, gene_id_dict):
    hit_dict = {}
    
    for i in seed_dict.keys():
        if gene_id_dict.has_key(i):
            for j in gene_id_dict[i]:
                if hit_dict.has_key(j):
                    hit_dict[j] = max([hit_dict[j], seed_dict[i]])
                else:
                    hit_dict[j] = seed_dict[i]
                    
    return hit_dict

def seed_annot_to_gene_ids(seed_annot_dict, gene_id_dict):
    annot_dict = {}
    for i in seed_annot_dict.keys():
        if gene_id_dict.has_key(i):
            for j in gene_id_dict[i]:
                if annot_dict.has_key(j):
                    if annot_dict[j].has_key("siRNA"):
                        annot_dict[j]["siRNA"] = min([annot_dict[j]["siRNA"], seed_annot_dict[i]["siRNA"]])
                    if annot_dict[j].has_key("GeneScore"):
                        annot_dict[j]["GeneScore"] = max([annot_dict[j]["GeneScore"], seed_annot_dict[i]["GeneScore"]])
                else:
                    annot_dict[j] = seed_annot_dict[i].copy()
    return annot_dict

def add_scores_to_genes(request, graph_db, gene_dict, property_name):
    
    #print "adding scores to genes"
    
    for i in gene_dict.keys():
        rels, metadata = cypher.execute(graph_db, "START n=node:User(name={user_name}), m=node:Gene(name={gene_name}) MATCH (n)-[r:ASSIGNED]-(m) RETURN r", {'user_name':str(request.user), 'gene_name':i})
        
        if len(rels) == 1:
            #print rels[0]
            cur_props = rels[0][0].get_properties()
            cur_props[property_name] = gene_dict[i]
            rels[0][0].set_properties(cur_props)
        else:
            raise Exception("ERROR: too many or two few relationships found")
    
    #print "done"
    
def make_result_table_dep (graph_db, use_table_query_str, samp_id, header_ord, otf_headers, hit_annot_dict, gene_assoc_dict):
    
    data = neo4j.CypherQuery(graph_db, use_table_query_str).execute(**{'samp_id':samp_id}).data
    
    #print samp_id
    #print data[0]
    
    rels = []
    
    temp_header_ord = header_ord[:]
    
    temp_header_ord.extend(otf_headers)
    
    for i in data:
        use_i = list(i.values)
        use_i.extend([None]*len(otf_headers))
        if hit_annot_dict.has_key(use_i[4]) and hit_annot_dict[use_i[4]].has_key('siRNA'):
            use_i[temp_header_ord.index("siRNA")] = hit_annot_dict[use_i[4]]['siRNA']
        elif hit_annot_dict.has_key(use_i[4]) and hit_annot_dict[use_i[4]].has_key('GeneScore'):
            use_i[temp_header_ord.index("GeneScore")] = hit_annot_dict[use_i[4]]['GeneScore']
        
        if gene_assoc_dict.has_key(use_i[4]):
            use_i[temp_header_ord.index("HitWalkerScore")] = gene_assoc_dict[use_i[4]]
        
        rels.append(use_i)
        
    rels.sort(key=lambda x: x[temp_header_ord.index("HitWalkerScore")], reverse=True)
    
    cur_rank = 1
    cur_score = str(rels[0][temp_header_ord.index("HitWalkerScore")])
    
    for i_ind, i in enumerate(rels):
        new_score = str(rels[i_ind][temp_header_ord.index("HitWalkerScore")])
        if new_score != cur_score:
            cur_rank += 1
            
        rels[i_ind][temp_header_ord.index("HitWalkerRank")] = cur_rank
        cur_score = new_score
    
    return rels, temp_header_ord


def get_hitwalker_score_from_table_dep(rels, rels_header, score_col, gene_col, rank_col):
    hitwalker_score = {}
    
    for i in rels:
        if i[rels_header.index(score_col)] != None:
            if hitwalker_score.has_key(i[rels_header.index(gene_col)]) == False:
                hitwalker_score[i[rels_header.index(gene_col)]] = i[rels_header.index(rank_col)]
    
    return hitwalker_score

def make_ranked_hit_dict(hit_dict):
    
    hit_dict_ranked_list = hit_dict.items()    
    hit_dict_ranked_list.sort(key=lambda x:x[1], reverse=True)
        
    r_hit_dict = {}
        
    for i_ind, i in enumerate(hit_dict_ranked_list):
        r_hit_dict[i[0]] = (i_ind + 1)
        
    return r_hit_dict

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

def make_seed_list_dep(graph_db, hit_dict):
        
        hit_symb_query = 'MATCH(n:Gene)-[r:KNOWN_AS]-(m) WHERE r.status = "symbol" AND n.name IN {hit_genes} RETURN n.name, m.name'
        
        data=neo4j.CypherQuery(graph_db, hit_symb_query).execute(**{"hit_genes":hit_dict.keys()})
        
        seed_list = []
        
        for i in data:
            if hit_dict.has_key(i.values[0]):
                seed_list.append({'gene':i.values[0], 'symbol':i.values[1], 'score':hit_dict[i.values[0]]})
        
        seed_list.sort(key=lambda x: x["score"], reverse=True)
        
        return seed_list

class RelationshipSet:
    def __init__(self):
        self.rel_set = {}
        self.rel_keys = []
        self.rel_key_pos = 0
    
    def add (self, rel, map_dict=None):
        if map_dict != None:
            for i in map(lambda x: x[0] + "." + x[1], itertools.product(map_dict[rel.start_node["name"]], map_dict[rel.end_node["name"]])):
                if self.rel_set.has_key(i) == False:
                    self.rel_set[i] = rel.get_properties()
                    self.rel_keys.append(i)
                else:
                    if self.rel_set[i] != rel.get_properties():
                        raise Exception("Found duplicate, discordant rels")
                
        else:
            raise Exception("Currently unimplemented")
            #self.rel_set.add(rel.start_node["name"] + "." + rel.end_node["name"])
            
    def direct_add (self,i):
        if len(i) > 0:
            use_name = str(i[0]) + "." + str(i[1])
            if self.rel_set.has_key(use_name) == False:
                self.rel_set[use_name] = {"score":i[2]}
                self.rel_keys.append(use_name)
            else:
                if  self.rel_set[use_name]["score"] != i[2]:
                     raise Exception("Found duplicate, discordant rels")
            
    
    def check (self, start_node, end_node, undirected=True, map_dict=None):
        
        if map_dict != None:
            if undirected == True:
                use_set = map(lambda x: str(x[0]) + "." + str(x[1]), itertools.product(map_dict[start_node], map_dict[end_node])) + map(lambda x: str(x[1]) + "." + str(x[0]), itertools.product(map_dict[start_node], map_dict[end_node]))
            else:
                use_set = map(lambda x: str(x[0]) + "." + str(x[1]), itertools.product(map_dict[start_node], map_dict[end_node]))
        else:
            if undirected == True:
                use_set = [start_node + "." + end_node, end_node + "." + start_node]
            else:
                raise Exception("Currently unimplemented")
            #    use_set = [[rel.start_node["name"] + "." + rel.end_node["name"]]]
        
        return any(map(lambda x: self.rel_set.has_key(x), use_set))
    
    def nodes (self):
        
        ret_set = set()
        
        for i in self.rel_set.keys():
            for j in i.split("."):
                ret_set.add(j)
                
        return list(ret_set)
    
    def __iter__(self):
        return self
    
    def next(self):
        
        if self.rel_key_pos != len(self.rel_keys):
            cur_val = self.rel_set[self.rel_keys[self.rel_key_pos]]
            split_key = self.rel_keys[self.rel_key_pos].split(".")
            split_key.append(cur_val)
            
            self.rel_key_pos += 1
            return split_key
        else:
            self.rel_key_pos = 0
            raise StopIteration

class DiagnosisNode (core.Node):
    def __init__(self,cypher_res):
        use_meta = cypher_res.get_properties()
        
        if use_meta.has_key("diagnosis") == False or use_meta["diagnosis"] == "Unknown":
            use_meta["diagnosis"] = "UnknownDiagnosis"
        
        self.node_dict = {'id':use_meta["name"]+use_meta["diagnosis"], 'display_name':use_meta["diagnosis"], 'attributes':{'node_type':use_meta["diagnosis"], 'other_nodes':[], 'indexed_name':'name', 'meta':{'node_cat':'Diagnosis'}}, 'children':core.NodeList()}
        self.id = use_meta["name"]+use_meta["diagnosis"]
        self.display_name = use_meta["diagnosis"]
    def todict (self):
        return super(DiagnosisNode, self).todict()

class GenderNode(core.Node):
    def __init__(self,cypher_node):
        cypher_res = cypher_node.get_properties()
        if cypher_res.has_key("gender") == False or cypher_res["gender"] == "Unknown":
            cypher_res["gender"] = "UnknownGender"
        self.node_dict = {'id':cypher_res["name"]+cypher_res["gender"], 'display_name':cypher_res["gender"], 'attributes':{'node_type':cypher_res["gender"], 'other_nodes':[], 'indexed_name':'name', 'meta':{'node_cat':'Gender'}}, 'children':core.NodeList()}
        self.id = cypher_res["name"]+cypher_res["gender"]
        self.display_name = cypher_res["gender"]
    def todict (self):
        return super(GenderNode, self).todict()

class RaceNode(core.Node):
    def __init__(self,cypher_node):
        cypher_res = cypher_node.get_properties()
        if cypher_res.has_key("race") == False or cypher_res["race"] == "Unknown":
            cypher_res["race"] = "UnknownRace"
        self.node_dict = {'id':cypher_res["name"]+cypher_res["race"], 'display_name':cypher_res["race"], 'attributes':{'node_type':cypher_res["race"], 'other_nodes':[], 'indexed_name':'name', 'meta':{'node_cat':'Race'}}, 'children':core.NodeList()}
        self.id = cypher_res["name"]+cypher_res["race"]
        self.display_name = cypher_res["race"]
    def todict (self):
        return super(RaceNode, self).todict()
    
class MetaNode(core.Node):
    def __init__(self, sample_nl):
        nl_types = sample_nl.types()
        type_count = collections.Counter(nl_types)
        disp_name = string.joinfields(map(lambda x: str(x[0]) + ' ('+ str(x[1]) + ')', type_count.items()), ',')
        self.node_dict = {'id':'meta_node_1', 'display_name':disp_name, 'attributes':{'node_type':'MetaNode', 'indexed_name':'name', 'meta':{}}, 'children':sample_nl}
        self.id = 'meta_node_1'
        self.display_name = disp_name
    def todict(self):
        return super(MetaNode, self).todict()
    
class PathwayNode(core.Node):
    def __init__(self, gene_nl, path_name):
        self.node_dict = {'id':path_name, 'display_name':path_name, 'attributes':{'node_type':'MetaNode', 'indexed_name':'name', 'meta':{'node_cat':'Pathway'}}, 'children':gene_nl}
        self.id = path_name
        self.display_name = path_name
    def todict(self):
        return super(PathwayNode, self).todict()

class LabIDNode(core.Node):
    def __init__(self,cypher_res):
        diag_child = DiagnosisNode(cypher_res[1])
        gend_child = GenderNode(cypher_res[2])
        race_child = RaceNode(cypher_res[2])
        
        att_meta = cypher_res[0].get_properties()
        att_meta.pop("name")
        for i in cypher_res[2].get_properties().items():
            att_meta[i[0]] = i[1]
        
        self.node_dict = {'id':cypher_res[0]["name"], 'display_name':cypher_res[0]["name"], 'attributes':{'node_type':'Sample', 'indexed_name':'name', 'meta':att_meta}, 'children':core.NodeList()}
        self.node_dict['children'].add(diag_child)
        self.node_dict['children'].add(gend_child)
        self.node_dict['children'].add(race_child)
        self.id = cypher_res[0]["name"]
        self.display_name = cypher_res[0]["name"]
        
    def todict (self):
        return super(LabIDNode, self).todict()


#need to add id and display name to below class
class SirnaChildNode(core.Node):
    def __init__(self,gene_node, samp_dict):
        self.node_dict = {'id':gene_node.id +'_siRNA', 'display_name':gene_node.display_name + '_siRNA', 'attributes':{'node_type':'siRNA', 'other_nodes':[], 'meta':{'node_cat':'Assay Result', 'type':[], 'score':[], 'is_hit':[]}}, 'children':core.NodeList()}
        for i in samp_dict.items():
            self.node_dict['attributes']['other_nodes'].append(i[0])
            self.node_dict['attributes']['meta']['type'].append(i[1][0][0])
            self.node_dict['attributes']['meta']['score'].append(i[1][0][1])
            self.node_dict['attributes']['meta']['is_hit'].append(i[1][0][2])
        if any(self.node_dict['attributes']['meta']['is_hit']):
            self.node_dict['attributes']['node_type'] += '_Hit'
        self.id = gene_node.id +'_siRNA'
        self.display_name = gene_node.display_name + '_siRNA'
    def todict (self):
        return super(SirnaChildNode, self).todict()

class BaseChildNode(core.Node):
    def __init__(self,gene_node, samp_dict, var):
        self.node_dict = {'id':gene_node.id +'_' + var, 'display_name':gene_node.display_name + '_' + var, 'attributes':{'node_type':var, 'other_nodes':[], 'meta':{'node_cat':'Assay Result', 'type':[], 'score':[], 'is_hit':[]}}, 'children':core.NodeList()}
        for i in samp_dict.items():
            self.node_dict['attributes']['other_nodes'].append(i[0])
            self.node_dict['attributes']['meta']['score'].append(i[1][0][0])
            self.node_dict['attributes']['meta']['is_hit'].append(i[1][0][1])
        if any(self.node_dict['attributes']['meta']['is_hit']):
            self.node_dict['attributes']['node_type'] += '_Hit'
        self.id = gene_node.id +'_' + var
        self.display_name = gene_node.display_name + '_' + var
    def todict (self):
        return super(BaseChildNode, self).todict()

class GeneScoreChildNode(core.Node):
    def __init__(self, gene_node, samp_dict):
        self.node_dict = {'id':gene_node.id +'_GeneScore','display_name':gene_node.display_name + '_GeneScore','attributes':{'node_type':'GeneScore', 'other_nodes':[], 'meta':{'node_cat':'Assay Result', 'score':[], 'is_hit':[]}}, 'children':core.NodeList()}
        for i in samp_dict.items():
            self.node_dict['attributes']['other_nodes'].append(i[0])
            self.node_dict['attributes']['meta']['score'].append(i[1][0][0])
            self.node_dict['attributes']['meta']['is_hit'].append(i[1][0][1])
        if any(self.node_dict['attributes']['meta']['is_hit']):
            self.node_dict['attributes']['node_type'] += '_Hit'
        self.id = gene_node.id +'_GeneScore'
        self.display_name = gene_node.display_name + '_GeneScore'
    def todict(self):
        return super(GeneScoreChildNode, self).todict()

class VariantChildNode(core.Node):
    def __init__(self,gene_node,var_name, samp_list):
        self.node_dict = {'id':gene_node.id+'_'+var_name, 'display_name':var_name, 'attributes':{'node_type':'Variants', 'other_nodes':samp_list,'meta':{'node_cat':'Assay Result', 'is_hit':[True]*len(samp_list)}}, 'children':core.NodeList()}
        self.id = gene_node.id+'_'+var_name
        self.display_name = var_name
        
    def todict(self):
        return super(VariantChildNode, self).todict()

def return_sample_nodes(res_list, nodes, request):
    dan='test'

def get_gene_names(res_list, nodes, request):
    
    for i in core.BasicResultsIterable(res_list):
        if len(i) > 0:
            nodes.add(core.GeneNode(i))
        
def get_gene_score(res_list, nodes, request):
    for i in core.BasicResultsIterable(res_list):
        #print i
        if len(i) > 0:
            if isinstance(i[0], tuple):
                use_i = i[:]
            else:
                use_i = [i[:]]
            
            gene_list = collections.defaultdict(list)
            for j in use_i:
                gene_list[j[1]].append([j[3], j[4]])
            gene_score = GeneScoreChildNode(nodes.getNode(use_i[0][0]), gene_list)
            nodes.addChild(use_i[0][0], gene_score)

def get_sirna_score(res_list, nodes, request):
    for i in core.BasicResultsIterable(res_list):
        if len(i) > 0:
            if isinstance(i[0], tuple):
                use_i = i[:]
            else:
                use_i = [i[:]]
            
            gene_list = collections.defaultdict(list)
            for j in use_i:
                gene_list[j[1]].append([j[3], j[4], j[5]])
            
            sirna = SirnaChildNode(nodes.getNode(use_i[0][0]), gene_list)
            nodes.addChild(use_i[0][0], sirna)

def get_exprs(res_list, nodes, request):
    
    for i in core.BasicResultsIterable(res_list):
        if len(i) > 0:
            if isinstance(i[0], tuple):
                use_i = i[:]
            else:
                use_i = [i[:]]
            
            #[[u'ENSG00000158258', u'07-00112', u'LowExpr', 0.199990661272727, True]]
            
            gene_list = collections.defaultdict(list)
            use_vars = set()
            
            for j in use_i:
                gene_list[j[1]].append([j[3], j[4]])
                use_vars.add(j[2])
            
            assert len(use_vars) == 1
            
            l_exp = BaseChildNode(nodes.getNode(use_i[0][0]), gene_list, list(use_vars)[0])
            nodes.addChild(use_i[0][0], l_exp)

def get_variants(res_list, nodes, request):
    
    for i in core.BasicResultsIterable(res_list):
        
        if len(i) > 0:
        
            if isinstance(i[0], tuple):
                use_i = i[:]
            else:
                use_i = [i[:]]
        
            gene_list = collections.defaultdict(list)
            for j in use_i:
                gene_list[j[1]].append(j[2])
                
            cur_gene = nodes.getNode(use_i[0][0])
            for k in gene_list.items():
                variant = VariantChildNode(cur_gene, k[0], k[1])
                nodes.addChild(cur_gene.id, variant)

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
        
        if from_node_dict["attributes"]["node_type"] == "Sample" and to_node_dict["attributes"]["node_type"] == "Gene":
            
            is_hit = rel_type.endswith("_Hit")
            temp_rel_name = rel_type.replace("_Hit", "")
            
            if temp_rel_name == config.data_types['target']:
                
                if request.session.has_key('hitwalker_score') and request.session['hitwalker_score'].nodeList().hasNode(to_node_dict["id"]) == True and from_node_dict['id'] == request.session['query_samples']['SampleID'][config.data_types['target']]:
                    temp_rel_name = "Ranked_" + temp_rel_name
                
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
    

def get_link_atts_dep (request, from_node, to_node, rel_type, props, is_hit=False):
    
    from_node_dict = from_node.todict()
    to_node_dict = to_node.todict()
    
    #print from_node_dict
    
    if from_node_dict['attributes'].has_key('node_type') and to_node_dict['attributes'].has_key('node_type'):
        
        if from_node_dict["attributes"]["node_type"] == "Sample" and to_node_dict["attributes"]["node_type"] == "Gene" and rel_type == "Variants" and request.session['hitwalker_score'].has_key(to_node_dict["id"]) and from_node_dict['id'] == request.session['query_samples']['SampleID']['Variants']:
            #print '*******'
            #print from_node_dict['id']
            #print request.session['query_samples']['LabID']['Variants']
            #print '*******'
            return {'type':'Ranked_Variants'}#'rank':request.session['hitwalker_score'][to_node_dict["id"]]+1}
        
        elif from_node_dict["attributes"]["node_type"] == "Sample" and to_node_dict["attributes"]["node_type"] == "Gene" and rel_type == "Variants":
            #print '*******'
            #print from_node_dict['id'], 'observed'
            #print '*******'
            return {'type':'Observed_Variants'}
        
        elif from_node_dict["attributes"]["node_type"] == "Sample" and to_node_dict["attributes"]["node_type"] == "Gene" and rel_type in set(["siRNA_Hit", "siRNA"]) and is_hit==True:
            return {'type':'Observed_siRNA'}#'score':samp_gene_score}
        elif from_node_dict["attributes"]["node_type"] == "Sample" and to_node_dict["attributes"]["node_type"] == "Gene" and rel_type in set(["siRNA_Hit", "siRNA"]) and is_hit==False:
            return {'type':'Possible_siRNA'}
        elif from_node_dict["attributes"]["node_type"] == "Sample" and to_node_dict["attributes"]["node_type"] == "Gene" and rel_type in set(["GeneScore_Hit", "GeneScore"]) and is_hit == True:
            return {'type':'Observed_GeneScore'}#'score':samp_gene_score}
        elif from_node_dict["attributes"]["node_type"] == "Sample" and to_node_dict["attributes"]["node_type"] == "Gene" and rel_type in set(["GeneScore_Hit", "GeneScore"]) and is_hit == False:
            return {'type':'Possible_GeneScore'}
        elif from_node_dict["attributes"]["node_type"] == "Gene" and to_node_dict["attributes"]["node_type"] == "Gene":
            return {'type':'STRING'}#'score':props['score']}
        
        else:
            return {'type':'Unknown'}#'rank':1000000}
    else:
        return {'type':'Unknown'}#'rank':1000000}

def get_lab_id (res_list, nodes, request):
    for i in core.BasicResultsIterable(res_list):
        if len(i) > 0:
            nodes.add(LabIDNode(i))

def sample_to_sample(query_sample, subj_sample, request, config_struct_nodes, cur_graph):
    return cur_graph


def no_rels(query_genes, subj_genes, request, config_struct_nodes, cur_graph):
    return cur_graph

#as gene_to_gene is used for determination of pathways add in the path_conf session variable instead of the standard string_conf
def gene_to_gene(query_genes, subj_genes, request, config_struct_nodes, cur_graph):
    
    import config
    
    session = cypher.Session(config.cypher_session)
    tx = session.create_transaction()
    
    print len(query_genes), len(subj_genes)
    
    cur_edge_set = RelationshipSet()
    
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

    #for i_ind, i in enumerate(itertools.product(query_genes, subj_genes)):
    #    if i_ind % 1000 == 0:
    #        #need to append to cur_graph['links']
    #        for i in core.BasicResultsIterable(tx.execute()):
    #            if len(i) > 0 and cur_edge_set.check(i[0], i[1]) == False:
    #                cur_edge_set.direct_add(i)
    #                cur_graph['links'].append({'source':cur_graph["nodes"].nodeIndex(i[0]), 'target':cur_graph["nodes"].nodeIndex(i[1]), 'attributes':get_link_atts(request, cur_graph["nodes"].getNode(i[0]), cur_graph["nodes"].getNode(i[1]), None, {'score':i[2]})})
    #        tx = session.create_transaction()
    #       
    #    else:
    #        tx.append(config.gene_rels['query'], {"FROM_GENE":i[0], "TO_GENES":[i[1]], "string_conf":core.iterate_dict(request.session, ['path_conf'])})
    #
    #and return cur_graph
    
    return cur_graph

def gene_to_lab_id (genes, labids, request, config_struct_nodes, cur_graph):
    
    #as the genes have a dependency on labids, retrieve the nodes again as well as the children (a little wasteful but direct for now...)
   
    #where labids should be a list of labIDs that will get directly substituted into the query
    param_list = [{'LABID':list(labids)}]
    
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

def get_pathways_sample (request, request_post):
    
    import config
    print request_post
    subj_nodes = map(lambda x:{'node_type':'Sample', 'id':x}, request_post['sample_name'])
    query_nodes = map(lambda x:{'node_type':'Pathway', 'id':x}, request_post['pathway_name'])
    
    new_edge_queries = copy.deepcopy(config.edge_queries)
    new_edge_queries['edges']['Gene']['Gene']['handler'] = gene_to_gene
    
    ###a hack until this can be better specified...
    #old_string = request.session['string_conf']
    #
    #print request.session['path_conf']
    #
    #request.session['string_conf'] = .99
    
    cur_graph = core.copy_nodes(subj_nodes, query_nodes, request, new_edge_queries, never_group=True)
    
    #request.session['string_conf'] = old_string
    
    return cur_graph['nodes'].tolist(), cur_graph['links'], query_nodes[0]['id']

def get_shortest_paths (request, request_post):
    
    from config import cypher_session
    session = cypher.Session(cypher_session)
    tx = session.create_transaction()
    
    if request_post.has_key('query_samples'):
        
        ##if bypassing table I think this would be necessary...
        ###get ride of the token...
        #request_post.pop("csrfmiddlewaretoken")
        #
        #print request.session['query_samples']
        #
        ###add this to the config file somehow...
        #for i in ['res_prob', 'zscore', 'gene_score', 'max_iter', 'conv_thresh', 'string_conf', 'sample_alias', 'query_samples']:
        #    request.session[i] = core.proper_type(request_post[i][0])
        #    request_post.pop(i)
        #where_template, necessary_vars,group_buttons, num_hs, group_hs = core.parse_filter_input(request_post)
        #
        #request.session['filter'] = core.make_updated_filter_dict (num_hs, group_hs, group_buttons)
        #
        #request.session['where_template'] = where_template
        #request.session['necessary_vars'] = necessary_vars
        #
        #then create a graph with the node(s) requested by the user 
        
        final_nodes_list = core.get_nodes(list(set(request.session['query_samples']['SampleID'].values())), 'LabID', request)
        
        node_names = final_nodes_list.display_names()
        
        if len(node_names) == 1:
            title = 'Sample: '+ node_names[0]
        else:
            title = 'Samples: ' + string.joinfields(node_names, ',')
        
        return final_nodes_list.tolist(), [], title
    else:
        
        cur_node_set = set(request_post['var_select'] + request_post['seed_select'])
        cur_edge_set = RelationshipSet()
        
        samp_set = set()
        
        for i in request.session['query_samples']['SampleID'].values():
            if i != None:
                samp_set.add(i)
        
        final_nodes_list = core.get_nodes(list(samp_set), 'LabID', request)
        
        sp_query = string.joinfields(['MATCH (n:Gene{name:{var_select}})-[:EXTERNAL_ID]-()-[:MAPPED_TO]-(np) WITH np MATCH (m:Gene{name:{seed_select}})-[:EXTERNAL_ID]-()-[:MAPPED_TO]-(mp) WITH mp, np MATCH p=(np)-[:ASSOC*1..2]->(mp)',
                    'WHERE ALL(x IN NODES(p) WHERE HAS(x.ensembl_path) AND x.ensembl_path = true) AND ALL(x IN RELATIONSHIPS(p) WHERE HAS(x.score) AND x.score > {score}) WITH p, REDUCE(val=0, x in RELATIONSHIPS(p)| val+x.score) as use_score,',
                    'LENGTH(RELATIONSHIPS(p)) AS min_len ORDER BY min_len, use_score DESC RETURN p limit 1'], " ")
        
        for var, seed in itertools.product(request_post['var_select'], request_post['seed_select']):
            tx.append(sp_query, {'var_select':var, 'seed_select':seed, 'score':int(request.session['string_conf']*1000)})
            
        all_path = tx.execute()
        
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

#was in core, now here as it really belongs part of 
#def rwr_to_gene_ids(res_vec, dim_dict, gene_id_dict):
#    #make a map of all genes to proteins
#    gene_dict = {}
#    
#    #if multiple proteins are present per gene, then the max is kept
#    for i in dim_dict.keys():
#        if gene_id_dict.has_key(i):
#            for j in gene_id_dict[i]:
#                if gene_dict.has_key(j):
#                    gene_dict[j] = max([gene_dict[j], res_vec[dim_dict[i]]])
#                else:
#                    gene_dict[j] = res_vec[dim_dict[i]]
#            
#    return gene_dict

def netprop_rwr(request, seed_gene_nl, initial_graph_file, string_session_name, res_prob_session_name, conv_thresh_session_name, max_iter_session_name):
    
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

def match_pathway(query):
    
    query = json.loads(query)
    
    from config import cypher_session
    
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/')
    
    query_list = [] 
    
    if len(query) == 0:
        return query, []
    else:
        sample_query = neo4j.CypherQuery(graph_db,'MATCH (g:Gene)-[:EXTERNAL_ID]-()-[:GENESET_CONTAINS]-(path) WITH path, COUNT(g) AS g_count WHERE g_count < 200 MATCH (gene:Gene)-[:EXTERNAL_ID]-()-[r:GENESET_CONTAINS]-(path) WHERE ALL(x IN '+str(json.dumps(query))+' WHERE ANY(y in x WHERE gene.name = y)) RETURN DISTINCT ID(path), path.name + " (n=" + g_count + ")"')
        
        for i in sample_query.execute().data:
            query_list.append({'id':i.values[0], 'text':i.values[1], 'search_list':[i.values[1]]})
        
    return query, query_list


def match_sample(query):
    
    from config import cypher_session
    
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/')
    
    query_list = []
    
    sample_query = neo4j.CypherQuery(graph_db,'MATCH (n:LabID) WHERE (n)<-[:PRODUCED]-()<-[:HAS_DISEASE]-() AND n.name =~ "'+query+'.*' +'" RETURN ID(n), n.name')
    
    for i in sample_query.execute().data:
        query_list.append({'id':i.values[0], 'text':i.values[1], 'search_list':[i.values[1]]})
    
    return query, query_list

def match_gene(query):
    
    from config import cypher_session
    
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/')
    query_list = []
    sample_query = neo4j.CypherQuery(graph_db,'MATCH (n:Symbol)<-[r:KNOWN_AS]-(m) WHERE n.name =~ "'+query+'.*' +'" RETURN DISTINCT ID(n), n.name , COLLECT(r.status), COLLECT(m.name)')
    
    #preferentially choose those labeled as symbol over synonym as the synonyms can be pretty far off...
    search_res = []
    
    for i in sample_query.execute().data:
        if any(map(lambda x: x == "symbol", i.values[2])):
            for j_ind, j in enumerate(i.values[2]):
                if j == "symbol":
                    search_res.append(i.values[3][j_ind])
        else:
            search_res = i.values[3]
        
        query_list.append({'id':i.values[0], 'text':i.values[1], 'search_list':search_res})
        
    return query, query_list
