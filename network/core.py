from py2neo import neo4j
import json
import ast
import re
import string
import copy
import random
import collections
import copy
import itertools
from py2neo import neo4j, cypher
import numpy as np
import network.models

###classes

class NodeList (object):
    def __init__(self):
        self.node_list = []
        self.node_names = []
        self.iter_pos = 0
        self.attributes = {}
    
    def getNode (self, name):
        return self.node_list[self.node_names.index(name)]
    
    def getByIndex(self, index):
        return self.node_list[index]
    
    def nodeIndex(self, name):
        return self.node_names.index(name)
    
    def hasNode(self, name):
        return name in set(self.node_names)
    
    def types(self):
        return map(lambda x: x.nodeType(), self.node_list)
    
    def display_names(self):
        return map(lambda x: x.display_name, self.node_list)
    
    def ids(self):
        return map(lambda x: x.id, self.node_list)
    
    def addChild(self, name, node):
        if isinstance(node, Node):
            cur_node = self.node_names.index(name)
            self.node_list[cur_node].addChild(node)
    
    def add (self, node):
        if isinstance(node, Node):
            self.node_list.append(node)
            self.node_names.append(node.id)
        else:
            raise Exception("node needs to be a subclass of Node")
    
    def extend(self, node_list):
        if isinstance(node_list, NodeList):
            self.node_list.extend(node_list.node_list)
            self.node_names.extend(node_list.node_names)
            
            self.attributes = merge_attributes(self.attributes, node_list.attributes)
        else:
            raise Exception("node_list needs to be a subclass of NodeList")
    
    #like extendIfNew but will add new children based on IDs if they don't already exist...
    def mergeChildren(self, node_list):
        if isinstance(node_list, NodeList):
            cur_names = set(self.node_names)
            for i in range(0, len(node_list.node_names)):
                if (node_list.node_names[i] in cur_names) == False:
                    self.node_names.append(node_list.node_names[i])
                    self.node_list.append(node_list.node_list[i])
                else:
                    self.node_list[i].node_dict['children'].extendIfNew(node_list.node_list[i].node_dict['children'])
            
            self.attributes = merge_attributes(self.attributes, node_list.attributes)
            
        else:
            raise Exception("node_list needs to be a subclass of NodeList")
    
    def extendIfNew(self, node_list):
        if isinstance(node_list, NodeList):
            cur_names = set(self.node_names)
            for i in range(0, len(node_list.node_names)):
                if (node_list.node_names[i] in cur_names) == False:
                    self.node_names.append(node_list.node_names[i])
                    self.node_list.append(node_list.node_list[i])
                    
            self.attributes = merge_attributes(self.attributes, node_list.attributes)
            
        else:
            raise Exception("node_list needs to be a subclass of NodeList")
    
    def json(self):
        dict_list = []
        for i in self.node_list:
            dict_list.append(i.todict())
        return json.dumps(dict_list)
    
    def summarizeChildren(self,att_func, sum_func):
        
        res_list = []
        
        for i in range(0, len(self.node_list)):
            
            if self.node_list[i].hasChildren():
                res_list.append(sum_func(map(att_func, self.node_list[i].children())))
            else:
                res_list.append(None)
            
        return res_list
    
    def filterByChild(self, att_func, sum_func):
        
        new_list = []
        new_list_names = []
        
        for i in range(0, len(self.node_list)):
            should_keep = sum_func(map(att_func, self.node_list[i].children()))
            if isinstance(should_keep, bool):
                if should_keep == True:
                    new_list.append(copy.deepcopy(self.node_list[i]))
                    new_list_names.append(self.node_names[i])
            else:
                raise Exception("The combination of att_func applied to each child followed by sum_func should result in a single logical value")
        
        self.node_list = copy.deepcopy(new_list)
        self.node_names = copy.deepcopy(new_list_names)
        
    
    def tolist(self):
        dict_list = []
        for i in self.node_list:
            dict_list.append(i.todict())
        return dict_list
    
    def __iter__(self):
        return self
    
    def next (self):
        if self.iter_pos == len(self.node_list):
            self.iter_pos = 0
            raise StopIteration
        else:
            cur_res = self.node_list[self.iter_pos]
            self.iter_pos += 1
            return cur_res
        
    def __len__(self):
        return len(self.node_list)

class SeedList():
    
    node_list = NodeList()
    node_scores = []
    
    def __init__(self, node_list):
        
        if isinstance(node_list, NodeList):
            self.node_list = node_list
            self.node_scores = node_list.summarizeChildren(lambda x: x.getAttr(["attributes","meta", "score"]), max)
        else:
            raise Exception("node_list needs to be a NodeList object")
    
    def subset(self, id_list):
        new_nl = NodeList()
        new_node_score = []
        
        for i in id_list:
            if self.node_list.hasNode(i):
                temp_ind = self.node_list.nodeIndex(i)
                new_nl.add(copy.deepcopy(self.node_list.getByIndex(temp_ind)))
                new_node_score.append(self.node_scores[temp_ind])
        
        self.node_list = new_nl
        self.node_scores = new_node_score
    
    def adjustScores(self, score_dict):
        if isinstance(score_dict, dict):
            for gene, score in score_dict.items():
                score_to_change = self.node_list.nodeIndex(gene)
                self.node_scores[score_to_change] = score
        else:
            raise Exception("score_dict needs to be a dictionary")
    
    def nodeList(self):
        return self.node_list
    
    def getScores(self, gene_list):
        
        ret_list = []
        for gene in gene_list:
            score_ind = self.node_list.nodeIndex(gene)
            ret_list.append(self.node_scores[score_ind])
            
        return ret_list
    
    def todict(self):
        """
            Returns A dictionary of  the form: {id:score}
        """
        ret_dict = {}
        
        node_ids = self.node_list.ids()
        
        for i in range(0, len(node_ids)):
            ret_dict[node_ids[i]] = self.node_scores[i]
            
        return ret_dict

class Node (object):
    node_dict = {'children':NodeList()}
    def todict (self):
        ret_dict = copy.deepcopy(self.node_dict)
        child_list = []
        for i in ret_dict.pop('children'):
            child_list.append(i.todict())
        ret_dict['children'] = child_list
        return ret_dict
    
    def nodeType(self):
        return self.node_dict['attributes']['node_type']
    
    def addChild (self, node):
        if isinstance(node, Node):
            self.node_dict['children'].add(node)
            
    def hasChildren(self):
        if len(self.node_dict['children']) > 0:
            return True
        else:
            return False
            
    def children (self):
        return self.node_dict['children']
    
    def getAttr(self, att_names):
        return iterate_dict(self.node_dict, att_names)

class RowNode(Node):
    
    def __init__(self, row, row_header):
        
        row_header = copy.copy(row_header)
        
        if ('query_ind' in set(row_header) and 'gene_ind' in set(row_header) and 'row_id' in set(row_header)) == False:
            raise Exception("Need to specify 'query_ind', 'row_id' and 'gene_ind' in the query result.")
        
        row_name = row[row_header.index('row_id')]
        
        self.node_dict = {'id':row_name, 'display_name':row_name, 'attributes':{'gene':str(row[row[row_header.index('gene_ind')]]), 'query':str(row[row[row_header.index('query_ind')]]), 'node_type':'row'}, 'children':NodeList()}
        
        row.pop(row_header.index('query_ind'))
        row_header.pop(row_header.index('query_ind'))
        
        row.pop(row_header.index('gene_ind'))
        row_header.pop(row_header.index('gene_ind'))
        
        row.pop(row_header.index('row_id'))
        
        self.display_name = self.node_dict['display_name']
        self.id = self.node_dict['id']
        
        self.node_dict['attributes']['row'] = row
        

class GeneNode(Node):
    
    def __init__(self,cypher_res):
        
        add_meta = {'node_cat':'Pathway_member', 'pathways':[]}
        
        for i in cypher_res[2]:
            add_meta['pathways'].append(i)
        
        self.node_dict = {'id':cypher_res[0], 'display_name':cypher_res[1], 'attributes':{'node_type':'Gene', 'indexed_name':'name', 'meta':add_meta}, 'children':NodeList()}
        self.display_name = cypher_res[1]
        self.id = cypher_res[0]
    def todict (self):
        return super(GeneNode, self).todict()
    def children(self):
        return super(GeneNode, self).children()

class MetaNode(Node):
    def __init__(self, sample_nl):
        nl_types = sample_nl.types()
        type_count = collections.Counter(nl_types)
        disp_name = string.joinfields(map(lambda x: str(x[0]) + ' ('+ str(x[1]) + ')', type_count.items()), ',')
        self.node_dict = {'id':'meta_node_1', 'display_name':disp_name, 'attributes':{'node_type':'MetaNode', 'indexed_name':'name', 'meta':{}}, 'children':sample_nl}
        self.id = 'meta_node_1'
        self.display_name = disp_name
    def todict(self):
        return super(MetaNode, self).todict()

class BasicSubjectChild(Node):
    #where cur_prop is a tuple of length 2
    def __init__(self, cur_prop):
        self.node_dict = {'id':string.joinfields(cur_prop), 'display_name':cur_prop[1], 'attributes':{'node_type':cur_prop[1], 'other_nodes':[], 'indexed_name':'name', 'meta':{'node_cat':cur_prop[0]}}, 'children':NodeList()}
        self.id = string.joinfields(cur_prop)
        self.display_name = cur_prop[1]
    def todict(self):
        return super(BasicSubjectChild, self).todict()
        
class SubjectNode(Node):
    
    def __init__(self, cypher_res):
        print cypher_res
        cypher_props = cypher_res[0].get_properties()
        
        self.node_dict = {'id':cypher_props["name"], 'display_name':cypher_props["name"], 'attributes':{'node_type':'Subject', 'indexed_name':'name', 'meta':{}}, 'children':NodeList()}
        
        for i in cypher_props.items():
            if (i[0] in set(['name', 'alias'])) == False:
                self.node_dict['children'].add(BasicSubjectChild(i))
        
        self.id = cypher_props["name"]
        self.display_name = cypher_props["name"]
                
    def todict (self):
        return super(SubjectNode, self).todict()
    

class BasicGeneChild(Node):
    def __init__(self,gene_node, samp_dict, var):
        self.node_dict = {'id':gene_node.id +'_' + var, 'display_name':gene_node.display_name + '_' + var, 'attributes':{'node_type':var, 'other_nodes':[], 'meta':{'node_cat':'Assay Result', 'type':[], 'score':[], 'is_hit':[]}}, 'children':NodeList()}
        for i in samp_dict.items():
            self.node_dict['attributes']['other_nodes'].append(i[0])
            self.node_dict['attributes']['meta']['score'].append(i[1][0][0])
            self.node_dict['attributes']['meta']['is_hit'].append(i[1][0][1])
        if any(self.node_dict['attributes']['meta']['is_hit']):
            self.node_dict['attributes']['node_type'] += '_Hit'
        self.id = gene_node.id +'_' + var
        self.display_name = gene_node.display_name + '_' + var
    def todict (self):
        return super(BasicGeneChild, self).todict()

class BasicChild(Node):
    
    def __init__(self,res_list):
        self.node_dict = {'id':res_list[0]+'_child', 'display_name':res_list[0] + '_child', 'attributes':{'meta':{'score':res_list[1]}}, 'children':NodeList()}
        
        self.id = self.node_dict['id']
        self.display_name = self.node_dict['display_name']
        
    def todict (self):
        return super(BasicChild, self).todict()
    def children(self):
        return super(BasicChild, self).children()

class BasicNode(Node):
    def __init__(self,res_list, only_child=False):
        self.node_dict = {'id':res_list[0], 'display_name':res_list[0], 'children':NodeList()}
        
        self.id = self.node_dict['id']
        self.display_name = self.node_dict['display_name']
        
        if only_child == True:
            if len(res_list) != 2:
                raise Exception("only_child == True requires res_list to be of length 2")
            self.node_dict['children'].add(BasicChild(res_list))
        
    def todict (self):
        return super(BasicNode, self).todict()
    def children(self):
        return super(BasicNode, self).children()

class SeedNode(Node):
    def __init__(self,gene_node, cypher_res, cypher_header):
        
        seed_type = cypher_res[cypher_header.index('var')]
        
        self.node_dict = {'id':gene_node.id + '_' + seed_type, 'display_name':gene_node.display_name + '_' + seed_type, 'attributes':{'node_type':seed_type, 'other_nodes':[], 'meta':{'node_cat':'Assay Result', 'type':[], 'score':[], 'is_hit':[]}}, 'children':NodeList()}
        
        #other_nodes doesn't really matter here...but technically is the other samples this hit can be observed in
        self.node_dict['attributes']['other_nodes'].append(cypher_res[1])
        
        for i in range(3, len(cypher_res)):
            self.node_dict['attributes']['meta'][cypher_header[i]] = cypher_res[i]
        
        self.id = self.node_dict['id']
        self.display_name = self.node_dict['display_name']
        
    def todict (self):
        return super(SeedNode, self).todict()
    
class BasicResultsIterable:
    """
        Produce a list of tuples corresponding to the values attribute
        of a cypher query
    """
    def __init__(self, cypher_results):
        self.cypher_results = iter(cypher_results)
    
    def __iter__(self):
        return self
    
    def next (self):
        cur_val = self.cypher_results.next()
        
        if len(cur_val) == 1:
            return list(cur_val[0].values)
        else:
            return map(lambda x: x.values, cur_val)


#need to add me to unit tests...
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

def handle_gene_hits(res_list, nodes, request):
    for i in BasicResultsIterable(res_list):
       
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
            
            if nodes.hasNode(use_i[0][0]):
                gene_score = BasicGeneChild(nodes.getNode(use_i[0][0]), gene_list, list(use_vars)[0])
                nodes.addChild(use_i[0][0], gene_score)
            else:
                print "Did not find node " + use_i[0][0] + " skipping..."
    

class TargetChildNode(Node):
    def __init__(self,gene_node,var_name, samp_list, node_type):
        self.node_dict = {'id':gene_node.id+'_'+var_name, 'display_name':var_name, 'attributes':{'node_type':node_type, 'other_nodes':samp_list,'meta':{'node_cat':'Assay Result', 'is_hit':[True]*len(samp_list)}}, 'children':NodeList()}
        self.id = gene_node.id+'_'+var_name
        self.display_name = var_name
        
    def todict(self):
        return super(TargetChildNode, self).todict()

def handle_dense_gene_targets(res_list, nodes, request):

    from config import data_types
    #print 'hello dan'
    #print res_list
    table_header = filter(lambda x: x.startswith("_")==False, res_list[0][0].__dict__.keys())
    print table_header
    
    subj_name = table_header.index('subject')
    gene_name = table_header.index('gene')
    
    print var_name

    for i in res_list:
        for j in i:
            
            j_m = map(lambda x: x[1], filter(lambda y: y[0].startswith("_")==False, j.__dict__.items()))
            j_m.append(str(j_m[ext_head[0]])+str(j_m[ext_head[1]]))
            
            print j_m       
            gene_list = collections.defaultdict(list)
            for k in use_j:
                gene_list['name'].append(k[subj_name])
            
            cur_gene = nodes.getNode(j_m[gene_name])
            
            for k in gene_list.items():
                variant = TargetChildNode(cur_gene, k[0], k[1], data_types['target'])
                nodes.addChild(cur_gene.id, variant)
                
def handle_gene_targets(res_list, nodes, request):
    
    if len(res_list) > 0:
        
        from config import data_types
        
        for i in res_list:
            
            if len(i) > 0:
                
                seed_header = cypherHeader(i)
                
                var_name = seed_header.index('query_ind')
                samp_name = seed_header.index('Sample')
                gene_name = seed_header.index('gene_ind')
                
                for j in BasicResultsIterable([i]):
                    if isinstance(j[0], tuple):
                        use_j = j[:]
                    else:
                        use_j = [j[:]]
                
                    #[(u'g.chr17:7578217G>A', u'uc002gim.2', u'7157', u'TP53', 157, 41, u'', u'Missense_Mutation', u'p.T211I', 0, 2, u'g.chr17:7578217G>A_7157', u'42MGBA_CENTRAL_NERVOUS_SYSTEM'), (u'g.chr17:7577093C>T', u'uc002gim.2', u'7157', u'TP53', 5, 21, u'', u'Missense_Mutation', u'p.R282Q', 0, 2, u'g.chr17:7577093C>T_7157', u'42MGBA_CENTRAL_NERVOUS_SYSTEM')]
                    
                    print use_j
                    
                    gene_list = collections.defaultdict(list)
                    for k in use_j:
                        gene_list[k[k[var_name]]].append(k[samp_name])
                    
                    cur_gene = nodes.getNode(use_j[0][use_j[0][gene_name]])
                    for k in gene_list.items():
                        variant = TargetChildNode(cur_gene, k[0], k[1], data_types['target'])
                        nodes.addChild(cur_gene.id, variant)

def make_sample_table (node, context, sample_link=""):
    
    if len(node) > 0:
        
        if sample_link != "":
            ret_str = '<a target="_blank" href="'+sample_link+node["id"]+'" class="list-group-item">Sample Info</a>'
        else:
            ret_str = ''
        
        ret_str += '<table class="table" style="width:300px">'
        ret_str += '<tbody>'
   
        for i in node['children']:
            ret_str += '<tr><th rowspan="1" valign="top">' + i['attributes']['meta']['node_cat'] + '</th>'
            ret_str += '<td>' + i['attributes']['node_type'] + '</td>'
            ret_str += '</tr>'
                    
        ret_str += '</tbody>'
        ret_str += '</table>'
        return ret_str
    else:
        return ''

def make_gene_table (node, context, gene_link="", plot_vars={}):
    
    import config
    
    ret_str = '<div class="container" style="width:300px">'
    
    sirna = []
    variant = []
    genescore = []
    
    ret_str += '<div class="row"><div class="media">'
    ret_str += '<div class="media-body">'
    ret_str += '<div class="media-heading"><h4>Gene Info</h4></div>'
    ret_str += '<ul class="list-group">'
    ret_str += '<a target="_blank" href="'+gene_link+node["id"]+'" class="list-group-item">Gene: '+node["id"]+'</a>'
    
    ret_str += '</ul></div></div></div>'
    
    child_dict = collections.defaultdict(list)
    
    for i in node['children']:
        child_dict[i['attributes']['node_type'].replace("_Hit", "")].append(i)
    
    #start with the query variable
    
    if child_dict.has_key(config.data_types['target']) == True:
        query_dta = child_dict.pop(config.data_types['target'])
        ret_str += make_hit_table(query_dta, config.data_types['target'], context, get_or_none(plot_vars, config.data_types['target']), True)
    
    for child_type, child_dta in child_dict.items():
        ret_str += make_hit_table(child_dta, child_type, context, get_or_none(plot_vars, child_type),False)
    
    ret_str += '</div>'
    
    #if the statement also contains a reference to gene name, insert it here...
    ret_str = ret_str.replace("$$GENE_NAME$$", node['display_name'])
    
    return ret_str

def make_hit_table(node, child_type, context, onclick_str=None, only_hits=False):
    if len(node) > 0:
        header = '<div class="row"><div class="media">'
        header += '<div class="media-body">'
        header += '<div class="media-heading"><h4>'+child_type +' Results</h4></div>'
        header += '<ul class="list-group">'
        
        if context == 'panel':
            dis_str = ""
        else:
            dis_str = "disabled"
        
        if onclick_str != None:
            use_onclick = "onclick=" + onclick_str
        else:
            use_onclick = ""
        
        for i in node:
            
            hit_samps = sum(i['attributes']['meta']['is_hit'])
            non_hit_samps = len(i['attributes']['meta']['is_hit']) - hit_samps
            
            hit_content = '<div class="media col-md-12"><div class="media-body"><ul class="list-group">'
            non_hit_content = '<div class="media col-md-12"><div class="media-body"><ul class="list-group">'
            
            if hit_samps == 0:
                hit_content = '<p class="text-danger">No samples to display</p>'
            else:
                for j_ind, j in enumerate(i['attributes']['meta']['is_hit']):
                    if j == True:
                        hit_content += '<li class="list-group-item btn" data-type="'+child_type+'" data-gene="$$GENE_NAME$$" '+use_onclick+'>' + i['attributes']['other_nodes'][j_ind] + '</li>'
                        
            if non_hit_samps == 0:
                non_hit_content = '<p class="text-danger">No samples to display</p>'
            else:
                for j_ind, j in enumerate(i['attributes']['meta']['is_hit']):
                    if j == False:
                        non_hit_content += '<li class="list-group-item btn" data-type="'+child_type+'" data-gene="$$GENE_NAME$$"  '+use_onclick+'>' + i['attributes']['other_nodes'][j_ind] + '</li>'
            
            non_hit_content += '</ul></div>'
            hit_content += '</ul></div>'
            
            if only_hits == True:
                
                header += '<li class="list-group-item">'+i['display_name']+':<span class="badge btn"  data-toggle="popover" data-title="Click to display" '+dis_str+'  data-content=\''+hit_content+'\' tabindex="0" data-trigger="focus" data-html=true>'+str(hit_samps)+'</span></li>'
                
            else:
            
                header += '<li class="list-group-item">Samples with Hits:<span class="badge btn"  data-toggle="popover" data-title="Click to display" '+dis_str+'  data-content=\''+hit_content+'\' tabindex="0" data-trigger="focus" data-html=true>'+str(hit_samps)+'</span></li>'
            
                header += '<li class="list-group-item">Samples without Hits:<span class="badge btn" data-title="Click to display", data-toggle="popover" '+dis_str+' data-content=\''+non_hit_content+'\' tabindex="0" data-trigger="focus" data-html=true>'+str(non_hit_samps)+'</span></li>'
        
        #did have extra div
        header += '</ul></div></div></div>'
        return header
    else:
        return ''

def make_metanode_table(node, context):
    import config
    #ret_str = '<div class=container>'
    
    samp_dict = {'Subject':[], 'Type':[], 'Value':[]}
    gene_dict = {}
    
    for i in node['children']:
        if i['attributes']['node_type'] in set(['Subject', 'Gene']):
            for j in i['children']:
                samp_dict['Subject'].append(i['display_name'])
                samp_dict['Type'].append(j['attributes']['meta']['node_cat'])
                
                #as it is in the form: Gene_siRNA
                if i['attributes']['node_type'] == 'Gene':
                    use_val = j['attributes']['node_type']
                else:
                    if config.node_abbreviations.has_key(j['display_name']):
                        use_val = config.node_abbreviations[j['display_name']]
                    else:
                        use_val = j['display_name']
                    
                samp_dict['Value'].append(use_val)
        else:
            raise Exception('node_type not currently implemented')

    samp_df = pd.DataFrame(samp_dict)
    samp_dict = samp_df.groupby(['Type', 'Value']).count().to_dict()
    
    table_dict = {}
    
    for i in samp_dict['Subject'].items():
        if table_dict.has_key(i[0][0]) == False:
            table_dict[i[0][0]] = [[i[0][1], i[1]]]
        else:
            table_dict[i[0][0]].append([i[0][1], i[1]])
    #ret_str += '<div class="row">'
    ret_str = '<table class="table" id="summary_table", style="width:300px">'
    ret_str += '<thead><tr><th>Type</th><th>Value</th><th>Count</th></tr></thead>'
    ret_str += '<tbody>'
   
    for i in table_dict.items():
        row_head = '<th rowspan="'+str(len(i[1]))+'" valign="top">' + i[0] + '</th>'
        for j_ind, j in enumerate(i[1]):
            ret_str += '<tr>'
            if j_ind == 0:
                ret_str += row_head
            ret_str += '<td>' + j[0].replace('_', ' ') + '</td>' + '<td><span class="badge btn" onclick="subset_meta_node(\''+i[0]+'\',\''+j[0]+'\', \''+node['id']+'\')">' + str(j[1]) + '</span></td>'
            ret_str += '</tr>'
                
    ret_str += '</tbody>'
    ret_str += '</table>'
    ret_str += '<div><a class="btn btn-default" id="nf_button" onclick="export_html_summary_csv(\'summary_table\')">Export Summary</a>'
                #<a class="btn btn-default" id="nf_button_2" onclick="send_to_new_panel()">Send to New Panel</a>'
    ret_str += '</div>'
    
    return ret_str

def get_or_none(inp_dict, val):
    if inp_dict.has_key(val):
        return inp_dict[val]
    else:
        return None

def merge_attributes(dict1, dict2):
    new_dict = dict1.copy()
    new_dict.update(dict2)
    return new_dict

def make_results_table(query_nl, seed_nl, ranking_sl):
    
    #need to add both the seed info and the ranking info to the table and return a sorted version (by ranking)
    
    table_header = copy.copy(query_nl.attributes['header'])
    
    #expand table_header by the seed types in seed_nl as well as another column for ranking_sl
    
    seed_types = set(reduce(lambda x,y: x+y, map(lambda x: x.children().types(), seed_nl)))
    
    table_header.extend(reduce(lambda x,y: x+y, [sorted(list(seed_types)), ['HitWalkerScore', 'HitWalkerRank']]))
    
    rows = []
    
    for i in query_nl:
        
        temp_i = i.getAttr(['attributes', 'row']) + [None]*(len(seed_types) + 2)
        temp_gene = i.getAttr(['attributes', 'gene'])
        
        #add in any seed info to the query table
        if seed_nl.hasNode(temp_gene):
            temp_seed_node = seed_nl.getNode(temp_gene)
            for child in temp_seed_node.children():
                row_pos = table_header.index(child.getAttr(['attributes', 'node_type']))
                temp_i[row_pos] = child.getAttr(['attributes', 'meta', 'score'])
        
        #add in any ranking info
        if ranking_sl.nodeList().hasNode(temp_gene):
            temp_i[table_header.index('HitWalkerScore')] = ranking_sl.getScores([temp_gene])[0]
        
        rows.append(temp_i)
    
    #then order the rows by 'HitWalkerScore' and record the result in HitWalkerRank
        #do this by looking for changes from the stringified hw scores
    
    cur_rank = 1
    score_pos = table_header.index("HitWalkerScore")
    rank_pos = table_header.index("HitWalkerRank")
    
    #sort first then check the following scores
    rows.sort(key=lambda x: x[score_pos], reverse=True)
    
    if rows[0][score_pos] != None:
    
        cur_score = str(rows[0][score_pos])
        
        for i_ind, i in enumerate(rows):
            if i[score_pos] != None:
                new_score = str(i[score_pos])
                if new_score != cur_score:
                    cur_rank += 1
                    
                rows[i_ind][rank_pos] = cur_rank
                cur_score = new_score
    
    return rows, table_header

def get_hitwalker_score_from_table(rels, rels_header, score_col, gene_col, rank_col):
    hitwalker_score = {}
    
    for i in rels:
        if i[rels_header.index(score_col)] != None:
            if hitwalker_score.has_key(i[rels_header.index(gene_col)]) == False:
                hitwalker_score[i[rels_header.index(gene_col)]] = i[rels_header.index(rank_col)]
    
    return hitwalker_score

def handle_hits(res_list, nodes, request):
    
    """
        A basic handler for hit-type queries, expects each element of the outer list to correspond to
        samples, with each row of the query providing the gene ID as the first element and the sample ID as the second as well as other
        hit related info especially the 'var' variable which indicates the type of hit and currently needs to be supplied 3rd in the query results.
        'res_list' needs to be the result from a Neo4j or related query compatible with BasicResultsIterable
    """
    
    seed_header = cypherHeader(res_list)
    
    for sample_result in BasicResultsIterable(res_list):
        for ind, gene_result in enumerate(sample_result):
            if nodes.hasNode(gene_result[0]) == False:
                #create the gene prior to adding if it doesn't exist 
                nodes.add(GeneNode([gene_result[0], gene_result[0], []]))
            
            seed_node = SeedNode(nodes.getNode(gene_result[0]), gene_result, seed_header)
            
            #then add the hits to the GeneNode
            nodes.addChild(gene_result[0], seed_node)
            
#recursively find the header for a CypherResult--kinda fragile
def cypherHeader(query_res):
    
    if len(query_res) > 0:
        if 'columns' in set(dir(query_res)):
            return(query_res.columns)
        else:
            return cypherHeader(query_res[0])
    else:
        return ()


def no_combining(request, valid_gene_hits):
    
    seed_list = SeedList(valid_gene_hits)
    
    #now encode hits as 1, everything else as 0 
    
    score_dict = {}
    
    for i in seed_list.nodeList().ids():
        score_dict[i] = 1
    
    seed_list.adjustScores(score_dict)
    
    return seed_list

def get_nodes_from_config(request, session_dict):
    
    """
        Carries out the queries specified in the session_dict for each sample
        contained within the session dict object similar to {'query_samples':{'SampleID':{'IDtype1':'ID', 'IDtype2':'ID'}}}.
    """
    
    sample_base_list = ['query_samples', 'SampleID']
    
    nl = NodeList()
    for var_type in session_dict:
        temp_list = copy.deepcopy(sample_base_list)
        temp_list.append(var_type)
        
        retr_samp = iterate_dict(request.session, temp_list)
        
        print retr_samp
        print var_type
        
        
        if retr_samp != None:
            nl.mergeChildren(get_nodes([[retr_samp]], var_type, request, config_struct=session_dict))
    
    return nl

def handle_dense_query(res_list, nodes, request):
    table_header = filter(lambda x: x.startswith("_")==False, res_list[0][0].__dict__.keys())
    
    nodes.attributes['header'] = table_header[:]
    
    ext_head = map(lambda x: table_header.index(x),['gene', 'name'])
    
    for i in res_list:
        for j in i:
            
            j_m = map(lambda x: x[1], filter(lambda y: y[0].startswith("_")==False, j.__dict__.items()))
            
            table_header.extend(["gene_ind", "query_ind", "row_id"])
            j_m.extend(ext_head)
            j_m.append(str(j_m[ext_head[0]])+str(j_m[ext_head[1]]))
            print j_m
            nodes.add(RowNode(j_m, copy.copy(table_header)))
    

def handle_query_prior(res_list, nodes, request):
    
    table_header = list(cypherHeader(res_list))
    
    for sample_table in BasicResultsIterable(res_list):
        for row in sample_table:
            nodes.add(RowNode(list(row), copy.copy(table_header)))
    
    table_header.remove("gene_ind")
    table_header.remove("query_ind")
    table_header.remove("row_id")
    
    nodes.attributes['header'] = table_header


def get_valid_querys(request, query_session_dict):
    
    return get_nodes_from_config(request, query_session_dict)

def get_valid_hits(request, hit_session_dict):
    
    gene_nodes = get_nodes_from_config(request, hit_session_dict)
    
    gene_nodes.filterByChild(lambda x: x.getAttr(["attributes", "meta", "is_hit"]), all)
    
    return gene_nodes

def customize_query(inp_query, **kwargs):
    new_query = {}
    
    for i in inp_query.keys():
        if kwargs.has_key(i):
            new_query[i] = kwargs[i](inp_query[i])
            
        else:
            new_query[i] = inp_query[i]
    
    diff_keys = set(kwargs.keys()).difference(set(inp_query.keys()))
    
    for i in diff_keys:
        new_query[i] = kwargs[i]
    
    return new_query

def iterate_dict(cur_dict, name_list):
    """
        Traverse the nested dict cur_dict in the order of the
        values in name_list.
    """
    for i in name_list:
        if cur_dict.has_key(i):
            cur_dict = cur_dict[i]
        else:
            return None
    return cur_dict



def specify_type_query_tmp(tmpl, ret_type, coll_type):
    
    """
        A helper function to simplify specification of a return and
        collection type for a common type of 'grouped' query template
    """
    
    new_tmpl = copy.deepcopy(tmpl)
    
    rep_list = [['$$ret_type$$', ret_type], ['$$coll_type$$', coll_type],
                ['$$lower_ret_type$$', ret_type.lower()], ['$$upper_ret_type$$', ret_type.upper()],
                ['$$lower_coll_type$$', coll_type.lower()], ['$$upper_coll_type$$', coll_type.upper()]]
    
    for i in rep_list:
        for j in ['query', 'title']:
            new_tmpl[j] = new_tmpl[j].replace(i[0], i[1])
            
    return new_tmpl

def build_dict_recursive(new_dict, key_list, value):
    
    if len(key_list) == 1:
        use_key = key_list.pop(0)
        new_dict[use_key] = value
    elif len(key_list) > 1:
        use_key = key_list.pop(0)
        if new_dict.has_key(use_key) == False:
            new_dict[use_key] = {}
        build_dict_recursive(new_dict[use_key], key_list, value)

def fix_jquery_array_keys (inp_dict):
    
    """
        Converts a dictionary with mangled keys (e.g. {'key[]':'val'}) to an unmangled version 
    """
    
    new_dict = {}
    for key, val in inp_dict.items():
        key_list = re.findall(r"\[([\w\d_]+)\]", str(key))
        main_key = str(re.sub(r"\[.*\]", "", key))
        key_list.insert(0, main_key)
        build_dict_recursive(new_dict, key_list, val)
        
    return new_dict

def copy_nodes (subj_nodes, query_nodes, request, query_dict, never_group=False, rel_types="combination", exclude_type=""):
    
    """
        Takes the nodes specified in subj_nodes, and adds in the nodes in query_nodes.
        Both subj_nodes and query_nodes need to have 'node_type' and 'id' keys.
        The nodes are grouped based on their edges using apply_grouping2 unless never_group = True.
        Returns a dict of nodes/edges where the nodes is a NodeList and the edges are a
        list of dicts.
    """
    
    all_node_dict = collections.defaultdict(set)
    
    #currently the type of subject nodes is fixed
    
    for i in subj_nodes:
        all_node_dict[i["node_type"]].add(i["id"])
    
    #first do the subj nodes
    all_nodes = NodeList()
    
    for key,val in all_node_dict.items():
        all_nodes.extendIfNew(get_nodes(val, key, request, config_struct=query_dict['nodes'], missing_param="skip"))
    
    #whereas query nodes can change type (e.g Pathway -> Gene)
    
    query_types = collections.defaultdict(list)
    
    for i in query_nodes:
        query_types[i["node_type"]].append(i["id"])
   
    query_nl = NodeList()
    
    print 'getting node info'
    
    #just want to create a base set of nodes here, will add the specific children and such below when iterating through the different rel types
    
    for key,val in query_types.items():
        temp_nl = get_nodes(val, key, request, config_struct=query_dict['nodes'], missing_param="skip")
        
        if len(temp_nl) > 0:
        
            #if they are the same type, add normally
            
            type_set = set(temp_nl.types())
            
            if key in type_set:
                for cur_id in val:
                    all_node_dict[key].add(cur_id)
            else:
                #otherwise modify all_node_dict accordingly
                if len(type_set) > 1:
                    raise Exception("Only 'Pathway' is a known case of different types")
                else:
                    node_ids = map(lambda x: x.id, temp_nl)
                    for cur_id in node_ids:
                        all_node_dict[list(type_set)[0]].add(cur_id)
                    query_nodes = filter(lambda x: x["node_type"] != key,query_nodes)
                    for i in node_ids:
                        query_nodes.append({'id':i, 'node_type':list(type_set)[0]})
            all_nodes.extendIfNew(temp_nl)
    
    print 'done'
    
    print 'getting relationship info'
    
    if (rel_types == "combination") or (len(all_node_dict.keys()) == 1):
        iter_func = itertools.combinations_with_replacement(all_node_dict.keys(), 2)
    elif (rel_types == "product"):
        iter_func = [tuple(all_node_dict.keys())]
    else:
        raise Exception("Unknown value for rel_types, should be 'combination' or 'product'")
    
    cur_graph = {'nodes':all_nodes, 'links':[]}
    
    
    #compute the edges for each compatible compbination
    for i in iter_func:
        print i
        for_query = iterate_dict(query_dict['edges'], i)
        new_i = list(i)
        new_i.reverse()
        rev_query = iterate_dict(query_dict['edges'], new_i)
        
        if for_query != None:
            cur_graph = query_dict['edges'][i[0]][i[1]]["handler"](all_node_dict[i[0]], all_node_dict[i[1]], request, query_dict['nodes'], cur_graph)
        elif rev_query != None:
            cur_graph = query_dict['edges'][i[1]][i[0]]["handler"](all_node_dict[i[1]], all_node_dict[i[0]], request, query_dict['nodes'], cur_graph)
        else:
            raise Exception("Neither for_query nor ret_query was matched, check config.py")
    
    print 'done'
    
    print 'grouping'
    
    cur_graph = apply_grouping2(cur_graph, map(lambda x: x["id"], query_nodes), never_group=never_group, exclude_type=exclude_type)
    
    print 'done'
    
    #cur_graph = core.apply_grouping(cur_graph, map(lambda x: x["id"], query_nodes))
    
    return cur_graph

def apply_grouping2(cur_graph, query_nodes, never_group=False, exclude_type=""):
    import config
    import sys
    
    new_graph = {'nodes':NodeList(), 'links':copy.deepcopy(cur_graph['links'])}
    #
    #temp_file = open('/Users/bottomly/Desktop/test_graph.json', 'w')
    #json.dump(cur_graph['links'], temp_file)
    #temp_file.close()
    
    if len(cur_graph['nodes']) > config.max_nodes:
        #divide nodes depending on connections to query nodes
        
        if len(query_nodes) == 0:
            temp_nl = NodeList()
            #TODO make metanodes by node type
            temp_nl.add(MetaNode(cur_graph['nodes']))
            return {'nodes':temp_nl, 'links':[]}
        else:
        
            #for each node, record all other nodes it has an edge with in the form: node-edge_type
            node_dict = collections.defaultdict(set)
            rev_dict = collections.defaultdict(set)
            
            #if a node only has an edge_group_bl type edge, then it is grouped with all the other nodes which belong to the same category
            temp_node_dict = collections.defaultdict(set)
            
            for i in cur_graph['links']:
                temp_node_dict[str(i['source'])].add(i['attributes']['type'])
                temp_node_dict[str(i['target'])].add(i['attributes']['type'])
            
            for i in cur_graph['links']:
                node_dict[str(i['source'])].add(str(i['target']) + '-' + i['attributes']['type'])
                node_dict[str(i['target'])].add(str(i['source']) + '-' + i['attributes']['type'])
            
            #edge_group_bl = set(config.edge_group_bl)
            #
            #for key,vals in temp_node_dict.items():
            #    if all(map(lambda x: x in edge_group_bl, vals)):
            #        node_dict[key] = set(['Blacklist'])
            
            #group into metanodes all those nodes with the same patterns
        
            for i in node_dict.keys():
                rev_dict[string.joinfields(map(str, sorted(node_dict[i])), '_')].add(i)
            
            for i in rev_dict.keys():
                #print i
                #check one node to see what type it is
                node_type = cur_graph['nodes'].getByIndex(int(list(rev_dict[i])[0])).nodeType()
                
                if (len(rev_dict[i]) > config.max_nodes) and (never_group==False) and (node_type != exclude_type):
                    temp_nl = NodeList()
                    for j in rev_dict[i]:
                        #print j
                        temp_nl.add(cur_graph['nodes'].getByIndex(int(j)))
                        for k_ind, k in enumerate(cur_graph['links']):
                            if k['source'] == int(j):
                                #print 'source'
                                new_graph['links'][k_ind]['source'] = len(new_graph['nodes'])
                            elif k['target'] == int(j):
                                #print 'target'
                                new_graph['links'][k_ind]['target'] = len(new_graph['nodes'])
                    new_graph['nodes'].add(MetaNode(temp_nl))
                else:
                    for j in rev_dict[i]:
                        #print j
                        for k_ind,k in enumerate(cur_graph['links']):
                            if k['source'] == int(j):
                                new_graph['links'][k_ind]['source'] = len(new_graph['nodes'])
                            elif k['target'] == int(j):
                                new_graph['links'][k_ind]['target'] = len(new_graph['nodes'])
                            #print new_graph['links'][k_ind]
                        new_graph['nodes'].add(cur_graph['nodes'].getByIndex(int(j)))
            
            nc_nodes = set(range(0, len(cur_graph['nodes']))).difference(set(map(int, node_dict.keys())))
            
            node_types = collections.defaultdict(list)
            
            #divide by node type and whether or not it is a query_node
            for i in nc_nodes:
                temp_node = cur_graph['nodes'].getByIndex(int(i))
                node_types[temp_node.nodeType()].append(int(i))
            
            for node_t,node_list in node_types.items():
                
                if (len(node_list) > config.max_nodes) and (never_group==False) and (node_t != exclude_type):
                    #make a metanode
                    temp_nl = NodeList()
                    for i in node_list:
                        temp_nl.add(cur_graph['nodes'].getByIndex(int(i)))
                    
                    #add back in links, if any, so that the metanode will have all the relationships of its contained nodes
                    for j in node_list:
                        for k_ind,k in enumerate(cur_graph['links']):
                            if k['source'] == j:
                                new_graph['links'][k_ind]['source'] = len(new_graph['nodes'])
                            elif k['target'] == j:
                                new_graph['links'][k_ind]['target'] = len(new_graph['nodes'])
                    new_graph['nodes'].add(MetaNode(temp_nl))
                elif len(node_list) > 0:
                    #regular set of nodes
                    for j in node_list:
                        for k_ind,k in enumerate(cur_graph['links']):
                            if k['source'] == j:
                                new_graph['links'][k_ind]['source'] = len(new_graph['nodes'])
                            elif k['target'] == j:
                                new_graph['links'][k_ind]['target'] = len(new_graph['nodes'])
                        new_graph['nodes'].add(cur_graph['nodes'].getByIndex(int(j)))
            
            #remove edges that point to the same nodes, for instance two genes who are now part of the same metanode...
            #also summarize duplicate edges
            
            del_list = []
            obs_type = set()
            
            for i_ind, i in enumerate(new_graph['links']):
                if i['source'] == i['target']:
                    del_list.append(i_ind)
                elif str(i['source']) + '.' + str(i['target']) + '.' + i['attributes']['type'] in obs_type:
                    del_list.append(i_ind)
                elif str(i['target']) + '.' + str(i['source']) + '.' + i['attributes']['type'] in obs_type:
                    del_list.append(i_ind)
                else:
                    obs_type.add(str(i['source']) + '.' + str(i['target']) + '.' + i['attributes']['type'])
            
            for i in sorted(del_list, reverse=True):
                new_graph['links'].pop(i)
            
           # print new_graph['links']
          #  print new_graph['nodes'].ids()
            
            return new_graph
    else:
        return cur_graph

def get_nodes(names, node_type, request, indexed_name="name",  config_struct = None, param_list=[], missing_param="fail", cypher_session=None):
    
    if missing_param != "fail" and missing_param != "skip":
        raise Exception("missing_param needs to be either fail or skip")
    
    if len(param_list) == 1:
        param_list = param_list*len(names)
    elif len(param_list) != 0 and len(param_list) != len(names):
        raise Exception("param_list needs to either be of length 1 or the same length as names")
    
    if config_struct == None:
        import config
        config_struct = config.node_queries
    
    if cypher_session == None:
        from config import cypher_session
    
    nodes = NodeList()
    
    if config_struct.has_key(node_type):
        for i_ind, i in enumerate(config_struct[node_type]):
            
            if (i.has_key('db_type') == False) or (i.has_key('db_type') and i['db_type'] == 'neo4j'):
                
                session = cypher.Session(cypher_session)
                tx = session.create_transaction()
                
                for j_ind, j in enumerate(names):
                    if j != None:
                        if i.has_key('session_params') and i['session_params'] != None:
                           use_param = {indexed_name:j}
                           for par in i['session_params']:
                                use_key = par[-1]
                                use_elem = iterate_dict(request.session, par)
                                if use_elem != None:
                                    use_param[use_key] = use_elem
                        else:
                            use_param = {indexed_name:j}
                            
                        if len(param_list) > 0:
                            for p_key in param_list[j_ind].keys():
                                use_param[p_key] = param_list[j_ind][p_key]
                            
                        use_query = i['query']
                        
                        for var_elem in request.session['where_vars']:
                            use_query = add_where_input_query(use_query, var_elem['where_statement'], var_elem['necessary_vars'], request.session['graph_struct'])
                        
                        missing_params = get_necessary_params(use_query).difference(set(use_param.keys()))
                        
                        if len(missing_params) == 0:
                            tx.append(use_query, use_param)
                        elif len(missing_params) > 0 and missing_param=="fail":
                            raise Exception("Cannot find parameter(s) " + str(list(missing_params)) + " and missing_param is set to 'fail'")
                        #otherwise this implies skip
                
                if len(names) > 0:
                    res_list = tx.commit()
                        
            elif i.has_key('db_type') and i['db_type'] == 'sql':
                from django.db.models import Q
                print 'db'
                print names
                print node_type
                print i['datatype']
                print param_list
                #also need to getattr Variation etc, sub something else for gene__in
                cur_mod = getattr(network.models, i['datatype'])
                res_list = []
                for j in names:
                    if isinstance(j, list) == False:
                        j = [j]
                    if i.has_key('colnames') == False:
                        comb_q = Q(**{node_type.lower()+"__in":j})
                    elif isinstance(i['colnames'], str):
                        comb_q = Q(**{i['colnames'].lower()+"__in":j})
                    elif len(i['colnames']) == 2:
                        comb_q = Q(**{i['colnames'][0].lower()+"__in":j}) | Q(**{i['colnames'][1].lower()+"__in":j})
                    
                    db_res = cur_mod.objects.using("data").filter(comb_q)
                    if len(db_res) > 0:
                        res_list.append(db_res)
                
            else:
                raise Exception("specified db_type is not defined")
            
            if len(res_list) > 0:
            
                i['handler'](res_list, nodes, request)
    else:
        raise Exception("config_struct does not have specified node_type")
    
    return nodes

def get_necessary_params(cypher_str):
    necessary_params = re.findall("\{([\w\d_]+)\}", cypher_str)
    return set(necessary_params)

def request_post_to_json(post_val):
    
    new_dict = {}
    for i in post_val.items():
        #print i
        #print new_dict
        if len(i[1]) == 1:
            cur_val = proper_type(i[1][0])
            if isinstance(cur_val, str):
                new_dict[i[0]] = cur_val
            else:
                new_dict[i[0]] = json.dumps(cur_val)
        else:
            new_dict[i[0]] = json.dumps(i[1])
    
    return new_dict

def type_of_value(var):
    try:
        return type(ast.literal_eval(var))
    except Exception:
        return str

def return_numeric(value):
    return float(value)

def return_binary(value):
    if value == "True":
        return 1
    else:
        return 0
    
def fix_prog_type (prog_type):
    use_prog_type = prog_type.strip('/')
    if use_prog_type != "":
        use_prog_type = "_" + use_prog_type
        
    return use_prog_type

#python manage.py shell
#import network.core
#network.core.compute_graph_struct (outfile="network/static/network/text_files/graph_struct.json")

def compute_graph_struct (outfile="network/static/network/text_files/graph_struct.json"):
    
    ret_dict = {}
    from config import cypher_session
    
    graph_db = neo4j.GraphDatabaseService(cypher_session+'/db/data/')
    
    #pick a random node as a starting point
    data = neo4j.CypherQuery(graph_db, 'MATCH (n) RETURN DISTINCT LABELS(n)').execute().data
    
    for i in data:
        if len(i.values[0]) > 0:
            print i
            if ret_dict.has_key(i.values[0][0]) == False:
                ret_dict[i.values[0][0]] = {}
            for j in neo4j.CypherQuery(graph_db, 'MATCH (n:'+ i.values[0][0] +')-[r]-(m) RETURN DISTINCT TYPE(r), LABELS(m)').execute().data:
                print j
                if len(j.values[1]) > 0:
                    ret_dict[i.values[0][0]][j.values[0]] = j.values[1][0]
    
    out_f = open(outfile, "w")
    json.dump(ret_dict, out_f)
    out_f.close()


def parse_parameters (param_dict, trans_funcs, request):
    
    necessary_where = []
    
    for key, val in param_dict.items():
        if val['type'] == 'standard':
            for field_key, field_val in val['fields'].items():
                #request.session[field_key] = field_val['comparison'] + ' ' + str(field_val['default'])
                request.session[field_key] = field_val['default']
        elif val['type'] == 'grouped':
            where_stat_list = []
            necessary_vars = set()
            subsubgroup_total = 0
            for group in val['default_groups']:
                where_stat_list.append("(")
                for subgroup in group:
                    where_stat_list.append("(")
                    
                    for subsubgroup in subgroup:
                        
                        val_field = val['fields'][subsubgroup['field']]
                        
                        if val_field['type'] == 'character':
                            #no comparison field for character types
                            use_comp = '='
                        else:
                            use_comp = val_field['comparison']
                        
                        if subsubgroup.has_key('default'):
                            
                            #will be the transformed value specified by the user or input as default
                            use_var =  trans_funcs[subsubgroup['field']](subsubgroup['default'])
                        else:
                            #instead will be the default specified as part of the fields
                            use_var =  trans_funcs[subsubgroup['field']](val_field['default'])
                        
                        #prepend logical statements if  conflicting '((' already in result
                        if str(val['logical_list'][subsubgroup_total]) != "null":
                            should_remove = 0
                            
                            if where_stat_list[-1] == "(":
                                for l in range(1,3):
                                    if where_stat_list[-l] == "(":
                                        should_remove += 1
                                temp_list = []
                                for l in range(0, should_remove):
                                    temp_list.append(where_stat_list.pop())
                                
                                where_stat_list.append(val['logical_list'][subsubgroup_total])
                                where_stat_list += temp_list
                            else:
                                where_stat_list.append(val['logical_list'][subsubgroup_total])
                            
                        sub_var_name = '$$' + val_field['required']['from'] + '$$'
                        
                        if val_field.has_key('needs_has'):
                            has_statement = 'HAS(' + sub_var_name + "." + val_field['var_name'] + ')'
                            basic_statement = sub_var_name + "." + val_field['var_name'] + " " + use_comp + " " + json.dumps(use_var)
                            where_stat_list.append(has_statement + ' AND ' + basic_statement)
                        else:
                            where_stat_list.append(sub_var_name + "." + val_field['var_name'] + " " + use_comp + " " + json.dumps(use_var))
                            
                        necessary_vars.add(val_field['required']['from'])
                        subsubgroup_total += 1
                    where_stat_list.append(")")
                where_stat_list.append(")")
            
            temp_where = string.joinfields(where_stat_list, " ")
            
            #make a pretty version for output
            
            pretty_where = re.sub("\$\$[\w_\d]+\$\$", "", temp_where)
            
            for temp_field in val['fields'].values():
                pretty_where = re.sub("\."+temp_field['var_name'], temp_field['name'], pretty_where)
        
            necessary_where.append({'necessary_vars':copy.deepcopy(necessary_vars), 'where_statement':temp_where, 'name':key, 'pretty_where':pretty_where})
            
        else:   
            raise Exception("Unsupported type for parameter dict")
            
    
    request.session.save()
    
    return necessary_where


def add_where_input_query(query_str, where_template, necessary_vars, graph_struct):
    
    node_type_dict, where_re = check_input_query_where(query_str, necessary_vars, graph_struct)
    
    if len(node_type_dict) == 0:
        return query_str
    else:
        for i in node_type_dict.keys():
            if i != "":
                where_template = where_template.replace("$$" + node_type_dict[i] + "$$", i)
        
        if where_re != None:
            if query_str.find(where_re) == -1:
                raise Exception("ERROR: Cannot find specified match for add_where_input_query")
            return query_str.replace(where_re, where_re + ' ' + where_template + ' AND ')
        else:
            #assume here that the RETURN statement is the best place to put it
            split_query = query_str.split("MATCH")
            
            if split_query[-1].find("WHERE") != -1:
                return query_str.replace("RETURN", " AND " + where_template + " RETURN")
            else:
                return query_str.replace("RETURN", "WHERE " + where_template + " RETURN")

#ensures that the variables defined in check_var is present in all the with statements
def check_input_query_with (query_str, check_var, on_missing="replace"):
    
    if len(check_var) != 1:
        raise Exception('check_var currently needs to be a set of length 1')
    
    withs = re.findall(r'WITH\s+(([\w,\s_\.\(\)]+?)\s+((?:MATCH)|(?:WHERE)|(?:RETURN])))', query_str)

    #should perhaps support renaming of variables, but not for now, so will just take the values on the right side of the AS statements 
    with_vars = map(lambda x: map(lambda y: re.sub(r'\s+', '', y), re.sub(r',.+AS', ",", x[1]).split(',')), withs)

    in_with = map(lambda x: len(check_var.intersection(set(x))) == len(check_var), with_vars)
    
    if sum(in_with) != len(in_with):
        
        if on_missing == "fail":
            raise Exception('Provided query does not have the necessary variables')
        elif on_missing == "replace":
            
            missing_withs = filter(lambda x: x[0]==False, itertools.izip(in_with, withs))
            
            sub_var = check_var.pop()
            
            for i in missing_withs:
                query_str = re.sub(i[1][0], sub_var+','+i[1][0], query_str)
                
            return query_str
            
        else:
            raise Exception("on_missing only accepts 'fail' or 'replace'")
            
    else:
        return query_str

def check_input_query_where (query_str, necessary_vars,graph_struct):
    
    #start match
    
    #for cypher <  neo4j 2.0
    #start_segment = re.findall(r'START\s+(.+?)\s+MATCH', query_str)
    
    #7/16/2014 adjustment to regular expression to allow querying using a label but without a constaint ie labid:LabID
    start_match = re.findall(r'\(([\w_\d]+)(:([\w_\d]+)(\{[\w_:\}\{\d]+\})*)\)', query_str)
    
    where_re = None
    
    start_node_type = {}
    
    for i in start_match:
        if start_node_type.has_key(i[0]) == False:
            start_node_type[i[0]] = i[2]
        else:
            raise Exception("Cannot have duplicated starting variables in cypher query")
    
    if len(start_node_type) == 0:
        raise Exception("Cannot find labeled starting points in cypher query")
    
    #start_node_type should be a dictionary of the form:
    #{n:variable_name}
    #use this dictionary to to start off the parsed match statements below...
    
    #the older way of doing things...
    #if len(set.intersection(set(start_vars), set(start_node_type.keys()))) != len(set.union(set(start_vars), set(start_node_type.keys()))):
    #    raise Exception("Expected start_node_type to specify all start variables in the form {start_var:variable_name}")
    
    #basic match
    #match_stats = re.findall(r'MATCH\s*([>:\[\]\(\)\w\d_<-]+)\s*[(?:WITH)(?:RETURN)]', query_str)
    match_where = re.findall(r'(MATCH\s*(\S+)\s*((?:WITH)|(?:RETURN)|(?:WHERE)))', query_str)
    
    #as tuples are not assignable...
    match_where = map(list, match_where)
    
    for i in start_match:
        for j in range(0, len(match_where)):
            match_where[j][1] = match_where[j][1].replace(i[1], "")
    
    #with_stats = re.findall(r'WITH\s+([\w,\s_]+?)\s+[(?:MATCH)(?:WHERE)(?:RETURN)]', query_str)
    with_where = re.findall(r'WITH\s+(([\w,\s_]+?)\s+((?:MATCH)|(?:WHERE)|(?:RETURN])))', query_str)

    if len(with_where) != (len(match_where) - 1):
        raise Exception("Unexpected number of with/match statements")
    
    match_res = []
    
    for i in map(lambda x: x[1], match_where):
        temp_i = map(lambda x: re.sub("[>\(\)\[\]<]", "", x), re.split("-|:", i))
        temp_dict = {'statement':i, 'graph':[]}
        temp_list = []
        for j in range(0,len(temp_i)):
            temp_list.append(temp_i[j])
            if j % 3 == 0 and j > 0:
                temp_dict['graph'].append(temp_list)
                temp_list = [temp_list[-1]]
        match_res.append(temp_dict)
    
    node_type_dict = start_node_type.copy()
    
    #Examine match_res at each step (i) the variables that are carried through are only recorded (via a WITH statement).
        #we would only be interested in the variables visible at the end of the statement which is recorded in node_type_dict.
      
    for i in range(0, len(match_res)):
        for j in range(0, len(match_res[i]['graph'])):
            cur_type = node_type_dict[match_res[i]['graph'][j][0]]
            node_type_dict[match_res[i]['graph'][j][3]] = graph_struct[cur_type][match_res[i]['graph'][j][2]]
            if match_res[i]['graph'][j][1] != "":
                node_type_dict[match_res[i]['graph'][j][1]] = match_res[i]['graph'][j][2]
        #if there is a where statement here, are all the available variables accounted for?
        if match_where[i][2] == 'WHERE':
            missing_vars = set.difference(necessary_vars, set(node_type_dict.values()))
            if len(missing_vars) == 0:
                where_re = match_where[i][0]
                break
        if len(with_where) > i:
            split_widths = set(with_where[i][1].replace(" ", "").split(","))
            for k in node_type_dict.keys():
                if k in split_widths == False:
                    node_type_dict.pop(k)
        if len(with_where) > i and with_where[i][2] == 'WHERE':
            missing_vars = set.difference(necessary_vars, set(node_type_dict.values()))
            if len(missing_vars) == 0:
                where_re = with_where[i][0]
                break
            #if there is a where statement here, are all the available variables accounted for?
    
    if node_type_dict.has_key(""):
        node_type_dict.pop("")
    
    current_vars = set(node_type_dict.values())
    
    missing_vars = set.difference(necessary_vars, current_vars)
    
    if len(missing_vars) > 0:
        #raise Exception("input query statement does not contain enough variables for requested where statement")
        return {}, where_re
    else:
        return node_type_dict, where_re

def proper_type(var):
    try:
        cur_type = type_of_value(var)
        
        if cur_type == type(None):
            return None
        else:
            if var == "":
                return ""
            else:
                return(json.loads(var))
        
    except Exception:
        return str(var)

#def add_user_input_to_dict(field_dict, request):
#    index_field_dict = copy.deepcopy(field_dict)
#    
#    #as the functions are not JSON serializable...
#    
#    for i in index_field_dict['fields'].keys():
#        index_field_dict['fields'][i].pop("trans")
#    
#    for i in index_field_dict['parameters'].keys():
#        index_field_dict['parameters'][i]['default'] = request.session[i]
#        
#    return index_field_dict

#def pop_if_empty(coll):
#    for i in coll.keys():
#        if len(coll[i]) == 0:
#            coll.pop(i)
#            
#    return coll
#
#def getNodeDictByAtt(nodeDict, attribute, value):
#    ret_vals = []
#    
#    child_list = nodeDict["children"]
#    
#    for i in child_list:
#        if i["attributes"].has_key(attribute) and i["attributes"][attribute] == value:
#            ret_vals.append(i)
#    
#    return ret_vals

#def apply_grouping(cur_graph, query_nodes, always_group=False):
#    
#    import config
#    import sys
#    from custom_functions import MetaNode
#    
#    new_graph = {'nodes':NodeList(), 'links':copy.deepcopy(cur_graph['links'])}
#    #
#    #temp_file = open('/Users/bottomly/Desktop/test_graph.json', 'w')
#    #json.dump(cur_graph['links'], temp_file)
#    #temp_file.close()
#    
#    if len(cur_graph['nodes']) > config.max_nodes:
#        #divide nodes depending on connections to query nodes
#        
#        if len(query_nodes) == 0:
#            temp_nl = NodeList()
#            #TODO make metanodes by node type
#            temp_nl.add(MetaNode(cur_graph['nodes']))
#            return {'nodes':temp_nl, 'links':[]}
#        else:
#            query_pos = set()
#            #node_dict = collections.defaultdict(lambda:collections.defaultdict(set))
#            node_dict = collections.defaultdict(set)
#            #find position of the query_nodes
#            
#            #need to determine the pattern of edges for each query_node
#            
#            for i in query_nodes:
#                query_pos.add(cur_graph['nodes'].nodeIndex(i))
#            
#            #find nodes connected to query_nodes
#            #group by connected nodes to determine query_node connection pattern
#            
#            #for i in cur_graph['links']:
#            #    if i['source'] in query_pos:
#            #        node_dict[i['target']].add(str(i['source']) + '.' + i['attributes']['type'])
#            #    elif i['target'] in query_pos:
#            #        node_dict[i['source']].add(str(i['target']) + '.' + i['attributes']['type'])
#            
#            #new since 10-17-2014 to allow seperation by node number in rev_dict
#            #for i in cur_graph['links']:
#            #    if i['source'] in query_pos:
#            #        node_dict[str(i['target'])][str(i['source'])].add(i['attributes']['type'])
#            #    elif i['target'] in query_pos:
#            #        node_dict[str(i['source'])][str(i['target'])].add(i['attributes']['type'])
#            
#            #even newer...
#            
#            for i in cur_graph['links']:
#                if i['source'] in query_pos:
#                    node_dict[i['target']].add(i['attributes']['type'])
#                elif i['target'] in query_pos:
#                    node_dict[i['source']].add(i['attributes']['type'])
#            
#            #then reverse the dictionary, sorting the lists first so that the dictionary is of edge patterns
#            
#            
#            rev_dict = collections.defaultdict(set)
#            
#            #for i in node_dict.keys():
#            #    #divide by node number first
#            #    for j in node_dict[i].keys():
#            #        rev_dict[j+string.joinfields(map(str, sorted(node_dict[i][j])), '_')].add(int(i))
#            #    #rev_dict[string.joinfields(map(str, sorted(node_dict[i])), '_')].add(i)
#            
#            for i in node_dict.keys():
#                rev_dict[string.joinfields(map(str, sorted(node_dict[i])), '_')].add(i)
#                
#            print rev_dict
#            
#            #for each of these groups, insert into metanode or leave it as is
#            
#            for i in rev_dict.keys():
#                #print i
#                if len(rev_dict[i]) > config.max_nodes or always_group==True:
#                    temp_nl = NodeList()
#                    for j in rev_dict[i]:
#                        temp_nl.add(cur_graph['nodes'].getByIndex(int(j)))
#                        for k_ind, k in enumerate(cur_graph['links']):
#                            if k['source'] == j:
#                                new_graph['links'][k_ind]['source'] = len(new_graph['nodes'])
#                            elif k['target'] == j:
#                                new_graph['links'][k_ind]['target'] = len(new_graph['nodes'])
#                    new_graph['nodes'].add(MetaNode(temp_nl))
#                else:
#                    for j in rev_dict[i]:
#                        #print j
#                        for k_ind,k in enumerate(cur_graph['links']):
#                            if k['source'] == j:
#                                new_graph['links'][k_ind]['source'] = len(new_graph['nodes'])
#                            elif k['target'] == j:
#                                new_graph['links'][k_ind]['target'] = len(new_graph['nodes'])
#                            #print new_graph['links'][k_ind]
#                        new_graph['nodes'].add(cur_graph['nodes'].getByIndex(int(j)))
#            
#            #also deal with any non-query connected nodes as well as the remaining query nodes not connected to anything else
#            
#            #need to stratify these nodes by node_type before placing into a metanode...
#            nc_nodes = set(range(0, len(cur_graph['nodes']))).difference(set(map(int, node_dict.keys())))#.union(query_pos)
#            
#            node_types = collections.defaultdict(list)
#            
#            print 'nc_nodes:', str(nc_nodes)
#            
#            #divide by node type and whether or not it is a query_node
#            for i in nc_nodes:
#                temp_node = cur_graph['nodes'].getByIndex(int(i))
#                if i in query_pos:
#                    node_types['query_'+str(i)].append(int(i))
#                else:
#                    node_types[temp_node.nodeType()].append(int(i))
#            
#            for node_t,node_list in node_types.items():
#                if len(node_list) > config.max_nodes or always_group==True:
#                    #make a metanode
#                    temp_nl = NodeList()
#                    for i in node_list:
#                        temp_nl.add(cur_graph['nodes'].getByIndex(int(i)))
#                    
#                    #add back in links, if any, so that the metanode will have all the relationships of its contained nodes
#                    for j in node_list:
#                        for k_ind,k in enumerate(cur_graph['links']):
#                            if k['source'] == j:
#                                new_graph['links'][k_ind]['source'] = len(new_graph['nodes'])
#                            elif k['target'] == j:
#                                new_graph['links'][k_ind]['target'] = len(new_graph['nodes'])
#                    new_graph['nodes'].add(MetaNode(temp_nl))
#                elif len(node_list) > 0:
#                    #regular set of nodes
#                    for j in node_list:
#                        for k_ind,k in enumerate(cur_graph['links']):
#                            if k['source'] == j:
#                                new_graph['links'][k_ind]['source'] = len(new_graph['nodes'])
#                            elif k['target'] == j:
#                                new_graph['links'][k_ind]['target'] = len(new_graph['nodes'])
#                        new_graph['nodes'].add(cur_graph['nodes'].getByIndex(int(j)))
#            
#            #remove edges that point to the same nodes, for instance two genes who are now part of the same metanode...
#            #also summarize duplicate edges
#            
#            del_list = []
#            obs_type = set()
#            
#            for i_ind, i in enumerate(new_graph['links']):
#                if i['source'] == i['target']:
#                    del_list.append(i_ind)
#                elif str(i['source']) + '.' + str(i['target']) + '.' + i['attributes']['type'] in obs_type:
#                    del_list.append(i_ind)
#                elif str(i['target']) + '.' + str(i['source']) + '.' + i['attributes']['type'] in obs_type:
#                    del_list.append(i_ind)
#                else:
#                    obs_type.add(str(i['source']) + '.' + str(i['target']) + '.' + i['attributes']['type'])
#            
#            for i in sorted(del_list, reverse=True):
#                new_graph['links'].pop(i)
#            
#            return new_graph
#        
#    else:
#        return cur_graph

#def process_request_post_dep (request_post):
#    import config
#    group_hs = {}
#    necessary_vars = set()
#    num_hs = {}
#    
#    for i in request_post.keys():
#        #print i
#        split_key = re.split("\.", i)
#        parent_group = string.joinfields(split_key[0:2], "_")
#        group = string.joinfields(split_key[2:4], "_")
#        var_type = string.joinfields(split_key[5:8], ".")
#        type_num = split_key[7]
#        
#        if num_hs.has_key(parent_group) == False:
#            num_hs[parent_group] = {}
#        if num_hs[parent_group].has_key(group) == False:
#            num_hs[parent_group][group] = {}
#        num_hs[parent_group][group][var_type] = str(type_num)
#        
#        if group_hs.has_key(parent_group) == False:
#            group_hs[parent_group] = {}
#        
#        if group_hs[parent_group].has_key(group) == False:
#            group_hs[parent_group][group] = {}
#            
#        group_hs[parent_group][group][var_type] = request_post[i][0]
#        
#        split_var = var_type.split(".")[0]
#        
#        if config.field_dict['fields'].has_key(split_var):
#            new_field_var = config.field_dict['fields'][split_var]
#            if new_field_var.has_key('required') and new_field_var['required'].has_key('from'):
#                necessary_vars.add(new_field_var['required']['from'])
#            else:
#                raise Exception("The filter field dictionary needs both a 'required' field and 'from' field")
#        else:
#            raise Exception("Need to supply the filter field dictionary with a key and values for "  + split_var)
#        
#    return num_hs, group_hs, necessary_vars

#def parse_filter_input_dep (request_post):
#    
#    import config
#    group_buttons = json.loads(request_post.pop("variant_filter_logic")[0])
#    
#    num_hs, group_hs, necessary_vars = process_request_post (request_post)
#    
#    where_stat_list = []
#    
#    for i in sorted(group_hs.keys(), key=lambda x: int(x.split("_")[1])):
#        where_stat_list.append("(")
#        
#        for j in sorted(group_hs[i].keys(), key=lambda x: int(x.split("_")[1])):
#            where_stat_list.append("(")
#            
#            for k in sorted(group_hs[i][j].keys(), key=lambda x: int(x.split(".")[2])):
#                
#                if group_buttons['comp_button'].has_key(num_hs[i][j][k]) == False:
#                    #assume input from select box
#                    if isinstance(group_hs[i][j][k], list):
#                        use_comp = 'IN'
#                    else:
#                        use_comp = '='
#                else:
#                    use_comp = group_buttons['comp_button'][num_hs[i][j][k]]
#                
#                use_var = config.field_dict['fields'][k.split(".")[0]]['trans'](group_hs[i][j][k])
#                
#                if group_buttons['logical_button'].has_key(num_hs[i][j][k]):
#                    should_remove = 0
#                    if where_stat_list[-1] == "(":
#                        for l in range(1,3):
#                            if where_stat_list[-l] == "(":
#                                should_remove += 1
#                        temp_list = []
#                        for l in range(0, should_remove):
#                            temp_list.append(where_stat_list.pop())
#                        
#                        where_stat_list.append(group_buttons['logical_button'][num_hs[i][j][k]])
#                        where_stat_list += temp_list
#                    else:
#                        where_stat_list.append(group_buttons['logical_button'][num_hs[i][j][k]])
#                
#                sub_var_name = '$$' + config.field_dict['fields'][k.split(".")[0]]['required']['from'] + '$$'
#                
#                if config.field_dict['fields'][k.split(".")[0]].has_key('needs_has'):
#                    has_statement = 'HAS(' + sub_var_name + "." + config.field_dict['fields'][k.split(".")[0]]['var_name'] + ')'
#                    basic_statement = sub_var_name + "." + config.field_dict['fields'][k.split(".")[0]]['var_name'] + " " + use_comp + " " + json.dumps(use_var)
#                    where_stat_list.append(has_statement + ' AND ' + basic_statement)
#                else:
#                    where_stat_list.append(sub_var_name + "." + config.field_dict['fields'][k.split(".")[0]]['var_name'] + " " + use_comp + " " + json.dumps(use_var))
#            
#            where_stat_list.append(")")
#        where_stat_list.append(")")
#    
#    where_stat = string.joinfields(where_stat_list, " ")
#    
#    return where_stat, necessary_vars, group_buttons, num_hs, group_hs
#

#not complete or tested...
#def simulate_graph_struct (graph_json_file, other_required_fields={}, node_range=[1,15], rel_range=[0,4]):
#    
#    graph_inp = open(graph_json_file, "r")
#    graph_json = json.load(graph_inp)
#    
#    sim_node_dict = {}
#    sim_rel_dict = {}
#    
#    for key,value in graph_json.items():
#        if sim_node_dict.has_key(key) == False:
#            num_samps = random.randint(node_range[0], node_range[1])
#            for i in range(0, num_samps):
#                if sim_node_dict.has_key(key) == False:
#                    sim_node_dict[key] = {key+"_"+str(i):{}}
#                else:
#                    sim_node_dict[key][key+"_"+str(i)] = {}
#                #also retrieve the other types of properties from the nodes first from config.field_dict, then from other_required_fields
#        #then generate relationships
#        for rel_key,rel_value in value.items():
#            #for each node
#            for i in sim_node_dict[key].keys():
#                num_rels = random.randint(rel_range[0], rel_range[1])
#                if sim_node_dict.has_key(rel_value):
#                    num_rels = max(len(sim_node_dict[rel_value].keys()), num_rels)
#                else:
#                    #add to the nodes list based on num_rels
#                    for j in range(0, num_rels):
#                        if sim_node_dict.has_key(rel_value) == False:
#                            sim_node_dict[rel_value] = {rel_value+"_"+str(j):{}}
#                        else:
#                            sim_node_dict[rel_value][rel_value+"_"+str(j)] = {}
#                perm_rels = range(0,num_rels)
#                random.shuffle(perm_rels)
#                for j in perm_rels:
#                    if sim_rel_dict.has_key(i) == False and sim_rel_dict.has_key(rel_value+"_"+str(j)) == False:
#                        if sim_rel_dict.has_key(i) == False:
#                            sim_rel_dict[i] = {rel_value+"_"+str(j):{"type":rel_key}}
#                        else:
#                            sim_rel_dict[i][rel_value+"_"+str(j)] = {"type":rel_key}
#    
#    for i in graph_json.items():
#        print i
#    
#    for i in sim_node_dict.items():
#        print i
#        
#    for i in sim_rel_dict.items():
#        print i
#    
#    
#    graph_inp.close()



#default_filters = [
#            [
#                [{'filter':'Cons_cat', 'default':'NonSynon.', 'logical':'AND'},    
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

#def make_updated_filter_dict (num_hs, group_hs, group_buttons):
#    
#    filter_list = []
#    
#    for i in sorted(group_hs.keys(), key=lambda x: int(x.split("_")[1])):
#        temp_list_1 = []
#        for j in sorted(group_hs[i].keys(), key=lambda x: int(x.split("_")[1])):
#            temp_list_2 = []
#            for k in sorted(group_hs[i][j].keys(), key=lambda x: int(x.split(".")[2])):
#                temp_dict = {}
#                temp_dict['filter'] =  k.split(".")[0]
#                temp_dict['default'] = group_hs[i][j][k]
#                if group_buttons['comp_button'].has_key(num_hs[i][j][k]):
#                    temp_dict['comparison'] = group_buttons['comp_button'][num_hs[i][j][k]]
#                if group_buttons['logical_button'].has_key(str(int(num_hs[i][j][k]))):
#                    temp_dict['logical'] = group_buttons['logical_button'][str(int(num_hs[i][j][k]))]
#                temp_list_2.append(temp_dict.copy())
#            temp_list_1.append(temp_list_2)
#        filter_list.append(temp_list_1)
#    
#    return filter_list
#            

if __name__ == "__main__":
    compute_graph_struct()