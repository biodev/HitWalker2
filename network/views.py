# Create your views here.
from django.http import HttpResponse, HttpResponseRedirect
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import SetPasswordForm
from django.views.decorators.cache import never_cache
from django.views.decorators.debug import sensitive_post_parameters
from django.contrib.auth import login, authenticate
from django.template import RequestContext, loader
from django.shortcuts import render, get_object_or_404, redirect
from django.core.urlresolvers import reverse
from django.core.files.storage import get_storage_class
from django.conf import settings
from django.contrib.staticfiles.templatetags.staticfiles import static
from django.views.decorators.csrf import ensure_csrf_cookie
from py2neo import neo4j , cypher
import numpy
import json
import sys
import re   
import string
import timeit
import copy
import ast
import itertools
import os
import config
import custom_functions
import core
import collections
import socket
import csv

#optionally need tinycss, cssselect, lxml, cairosvg

import tempfile
from network.models import user_parameters

import cssutils
import colour

use_css = ''

prog_type = config.prog_type

if prog_type != "":
    prog_type = "/" + prog_type

graph_inp = open(config.graph_struct_file, "r")
graph_struct = json.load(graph_inp)
graph_inp.close()

graph_edge_struct = {}

for i in graph_struct.keys():
    for j in graph_struct[i]:
        if graph_edge_struct.has_key(j) == False:
            graph_edge_struct[j] = {}
        graph_edge_struct[j][i]= graph_struct[i][j]
        

def get_css_name_colors(css_obj, by="color", only_default=True):
    css_dict = {}

    for rule in css_obj:
        if (only_default == True and rule.selectorText.find('.default') == 0) or (only_default == False):
            for prop in rule.style:
                if prop.name == 'fill' and prop.value != 'none':
                    if by == "color":
                        css_dict[str(colour.Color(prop.value))] = rule.selectorText
                    elif by == "name":
                        css_dict[rule.selectorText] = str(colour.Color(prop.value))
                    else:
                        raise Exception("Unknown 'by' value specified")
    return css_dict

def generate_css (user_name):
    
    #from http://stackoverflow.com/questions/20001282/django-how-to-reference-paths-of-static-files-and-can-i-use-them-in-models?rq=1
    static_storage = get_storage_class(settings.STATICFILES_STORAGE)()
    css_path =os.path.join(static_storage.location, "network/static/network/css/HitWalker2.css")
    hw_css = cssutils.parseFile(css_path)
    
    node_type_transl={'Gene':{'name':'Gene', 'class':'Gene'}, 'Sample':{'name':'Sample', 'class':'Sample'}}
    edge_type_transl={'STRING':{'name':'STRING', 'class':'STRING'}}
    #for recording the colors of the hits/targets and subtracting them from the available defaults
    important_defaults = set()
    used_defaults = set()
    
    #remove the .'s from the class name before passing up to JS
    default_css_classes = set(get_css_name_colors(hw_css, by="name", only_default = True).keys())
    
        #if the users specified a node_colors attribute for config.py
        #check the colors against the defaults in hw_css
            #if a match is found, add to node_type_transl pointing to the default value
            #otherwise keep as is and add to hw_css
    #try:
    css_dict = get_css_name_colors(hw_css, by="color", only_default = True)
    
    try:
        node_colors = getattr(config, "node_colors")
    except:
        node_colors = {}
    
    print node_colors
    
    for css_name,color in node_colors.items():
        
        if node_type_transl.has_key(css_name) == False:
            node_type_transl[css_name] = {"name":css_name, "class":css_name}
        
        use_color = str(colour.Color(color))
        if css_dict.has_key(use_color):
            node_type_transl[css_name]['class'] = css_dict[use_color]
            used_defaults.add(css_dict[use_color])
        else:
            hw_css.add('.'+css_name+'{fill:'+use_color+';}')
        
    #refresh the css dict here and create it by name
    css_name_dict = get_css_name_colors(hw_css, by="name", only_default = False)
    
    #for the entries in node_abbreviations, replace the key where appropriate
    
    try:
        node_abbreviations = getattr(config, "node_abbreviations")
    except:
        node_abbreviations = {}
    
    for node,abrev in node_abbreviations.items():
        if node_type_transl.has_key(abrev):
            node_type_transl[node] = node_type_transl.pop(abrev)
        else:
            node_type_transl[node] = {'name':abrev, 'class':None}
    
    #similarly, add in the edge_abbreviations
    try:
        edge_abbreviations = getattr(config, "edge_abbreviations")
    except:
        edge_abbreviations = {}
    
    for edge, abrev in edge_abbreviations.items():
        edge_type_transl[edge] = {'name':abrev, 'class':None}
    
    #make sure that the seeds and target specified in data_types are represented so the hits can be rendered correctly
    #use the values defined in data_list for the _Hits with the exception of the target
    
    temp_hits = set(config.data_list).difference(set([config.data_types['target']]))
    
    print temp_hits
    
    for node in map(lambda x: x+'_Hit', temp_hits) + [config.data_types['target']]:
        if node_type_transl.has_key(node) == False:
            #choose a default for it
            new_color = default_css_classes.difference(used_defaults).pop()
            node_type_transl[node] = {'name':node, 'class':new_color.replace('.', '')}
            used_defaults.add(new_color)
            
            important_defaults.add(new_color)
    
    print node_type_transl
    
    #if it is a _Hit or can be translated to one, then add in the appropriate class for the non-hit to hw_css
    #additionally, add in the edges corresponding to hits in a similar fashion
    for node_name in node_type_transl.keys():
        if node_name.find("_Hit") != -1 and node_type_transl[node_name]['class'] != None:
            node_rep = node_name.replace("_Hit", "")
            
            if node_type_transl.has_key(node_rep) == False:
                node_type_transl[node_rep] = {'name':node_rep, 'class':node_rep}
            
            hw_css.add(
                '.'+node_type_transl[node_rep]['name'] + \
                 '{' + \
                    'stroke:'+css_name_dict['.'+node_type_transl[node_name]['class']]+';\n' + \
                    'stroke-width: 2;\n' + \
                    'fill: '+css_name_dict['.Gene']+';\n' + \
                 '}'
            )
            if node_type_transl[node_name.replace("_Hit", "")]['class'] == None:
                node_type_transl[node_rep]['class'] = copy.copy(node_type_transl[node_rep]['name'])
            
            #add in edge info as well
            #'Possible_' edges
            if edge_type_transl.has_key('Possible_'+node_rep) and edge_type_transl['Possible_'+node_rep].has_key('name'):
                #set the class name to the 'name' version of this
                use_p_edge_name = edge_type_transl['Possible_'+node_rep]['name']
                edge_type_transl['Possible_'+node_rep]['class'] = use_p_edge_name
            else:
                use_p_edge_name = 'Possible_'+node_rep
                edge_type_transl['Possible_'+node_rep] = {'name':use_p_edge_name, 'class':use_p_edge_name}
                
            
            hw_css.add(
                '.'+use_p_edge_name + \
                 '{' + \
                    'stroke:'+css_name_dict['.'+node_type_transl[node_name]['class']]+';\n' + \
                    'stroke-dasharray:5,5;\n' + \
                 '}'
            )
            #'Observed_' edges
            if edge_type_transl.has_key('Observed_'+node_rep) and edge_type_transl['Observed_'+node_rep].has_key('name'):
                use_o_edge_name = edge_type_transl['Observed_'+node_rep]['name']
                edge_type_transl['Observed_'+node_rep]['class'] = use_o_edge_name
            else:
                use_o_edge_name = 'Observed_'+node_rep
                edge_type_transl['Observed_'+node_rep] = {'name':'Observed_'+node_rep, 'class':'Observed_'+node_rep}
                
            hw_css.add(
                '.'+use_o_edge_name + \
                 '{' + \
                    'stroke:'+css_name_dict['.'+node_type_transl[node_name]['class']]+';\n' + \
                 '}'
            )
    
    #Finally, fill in the target edge info:
    
    for i in ['Ranked', 'Observed']:
        key_val = i+'_'+config.data_types['target']
        if edge_type_transl.has_key(key_val) and edge_type_transl[key_val].has_key('name'):
            use_v_edge_name = edge_type_trans[key_val]['name']
            edge_type_transl[key_val]['class'] = use_v_edge_name
        else:
            use_v_edge_name = key_val
            edge_type_transl[key_val] = {'name':key_val, 'class':key_val}
        
        if i == 'Observed':
            dash_str = 'stroke-dasharray:5,1;\n'
        else:
            dash_str = ''
        
        hw_css.add(
            '.'+use_v_edge_name + \
            '{' + \
                    'stroke:'+css_name_dict['.'+node_type_transl[config.data_types['target']]['class']]+';\n' + \
                    dash_str + \
            '}'
        )
    
    #now serialize the new css appropriately in a unique file
    new_css_path = os.path.join(static_storage.location, "network/static/network/css/HitWalker2_"+user_name+".css")
    new_css_out = open(new_css_path, "w")
    new_css_out.write(hw_css.cssText)
    new_css_out.close()
    
    default_css_list = map(lambda x: x.replace(".", ""), default_css_classes.difference(important_defaults))
    
    return default_css_list, new_css_path, node_type_transl, edge_type_transl
    

@sensitive_post_parameters()
@never_cache
@login_required(login_url=prog_type + '/HitWalker2/login/')
#adapted from http://stackoverflow.com/questions/13900357/how-to-use-django-usercreationform-correctly
def password(request):
    
    if request.method =='POST':
        form = SetPasswordForm(request.user, request.POST)
        if form.is_valid():
            user = form.save()
            #should already be logged in...
            ##assume username and password are valid as they were just created....
            #auth_user = authenticate(username=request.POST['username'], password=request.POST['password1'])
            #
            #login(request, auth_user)
            
            return redirect(prog_type + '/HitWalker2/')
    else:
        form = SetPasswordForm(request.user)

    return render(request, 'network/password.html', {'form':form, 'prog_type':prog_type, 'username':request.user})

def get_default_parameters(request):
   
    index_field_dict = copy.deepcopy(config.adjust_fields)
    for i in index_field_dict.keys():
        for j in index_field_dict[i]['fields'].keys():
            if index_field_dict[i]['fields'][j].has_key('trans'):
                index_field_dict[i]['fields'][j].pop('trans')
        if index_field_dict[i]['type'] == 'grouped':
            logical_list = reduce (lambda x,y: x+y, reduce (lambda x,y: x+y,index_field_dict[i]['default_groups']))
            index_field_dict[i]['logical_list'] = map(lambda x: 'null' if x.has_key('logical')==False else x['logical'], logical_list)
    
    return HttpResponse(json.dumps({'adjust_fields':index_field_dict}),  mimetype="application/json")

@login_required(login_url=prog_type + '/HitWalker2/login/')
def index(request, retry_message):
    
    index_field_dict = copy.deepcopy(config.adjust_fields)
    
    ##as the functions are not JSON serializable...
    for i in index_field_dict.keys():
        for j in index_field_dict[i]['fields'].keys():
            if index_field_dict[i]['fields'][j].has_key('trans'):
                index_field_dict[i]['fields'][j].pop('trans')
        if index_field_dict[i]['type'] == 'grouped':
            logical_list = reduce (lambda x,y: x+y, reduce (lambda x,y: x+y,index_field_dict[i]['default_groups']))
            index_field_dict[i]['logical_list'] = map(lambda x: 'null' if x.has_key('logical')==False else x['logical'], logical_list)
            
    
    #if request.session.has_key('filter'):
    #    use_filter = request.session['filter']
    #else:
    #    use_filter = config.default_filters
    #
    #if len(set(index_field_dict['parameters'].keys()) - set(request.session.keys())) == 0:
    #    for i in index_field_dict['parameters'].keys():
    #        index_field_dict['parameters'][i]['default'] = request.session[i]
    
    param_names = []
    
    for i in user_parameters.objects.filter(user=request.user):
        param_names.append(i.name)
    
    request.session['graph_struct'] = graph_struct
    
    #
    #context = {'retry_message':retry_message, 'filter':json.dumps(index_field_dict), 'cur_filts': json.dumps(use_filter),
    #           'default_filts':json.dumps(config.default_filters), 'html_defaults':json.dumps(index_field_dict['parameters']),
    #           'param_names':json.dumps(param_names), 'prog_type':prog_type, 'username':request.user}
    
    context = {'adjust_fields':index_field_dict, 'param_names':json.dumps(param_names), 'prog_type':prog_type, 'username':request.user, 'data_types':config.data_types}
    
    return render(request, 'network/index.html', context)

@login_required(login_url=prog_type + '/HitWalker2/login/')
def qtests(request):
    return render(request, 'network/qtests.html')

@login_required(login_url=prog_type + '/HitWalker2/login/')
def table(request):
    
    if len(request.POST) == 0:
        return redirect(prog_type + '/HitWalker2/')
    elif request.POST.has_key("redirect"):
        
        seed_list = request.session['seed_list']
        query_list = request.session['query_list']
        rels = request.session['rels']
        rels_header = request.session['rels_header']
        
        where_vars = request.session['where_vars']
        inp_params = request.session['inp_params']
        
        cur_param = {}
        
        for key,val in inp_params.items():
            if val['type'] == 'standard':
                cur_param[key] = copy.deepcopy(val['fields'])
        
        cur_filts = {}
        
        for i in where_vars:
            cur_filts[i['name']] = i['pretty_where']
    
        return render(request, 'network/res_table.html', {'res_table':rels, 'res_header':rels_header, 'sample_display_name':request.session['query_samples']['SampleID'],
                                                        'cur_filts':json.dumps(cur_filts), 'cur_param':json.dumps(cur_param), 'seeds':json.dumps(seed_list[:50]), 'queries':json.dumps(query_list[:50]), 'prog_type':prog_type,
                                                        'username':request.user})
    
    #    
    #    cur_filts = re.sub("\$\$[\w_\d]+\$\$\.", "", request.session['where_template'])
    #    
    #    index_field_dict = copy.deepcopy(config.field_dict)
    #    
    #    temp_header_ord = config.header_ord[:]
    #    
    #    
    #    #probably need to add more to me...
    #    temp_header_ord.append("HitWalker_Score")
    #    
    #    #as the functions are not JSON serializable...
    #    for i in index_field_dict['parameters'].keys():
    #        index_field_dict['parameters'][i]['default'] = request.session[i]
    #    
    #    return render(request, 'network/res_table.html', {'res_table':request.session['rels'], 'res_header':temp_header_ord, 'sample_display_name':request.session['sample_info']['display_name'],
    #                                                      'sample_info':'{"node_type":"Sample", "name":"' + request.session['sample_info']['name']  + '","indexed_name":"' +  request.session['sample_info']['indexed_name'] + '"}',
    #                                                    'cur_filts':cur_filts, 'cur_param':json.dumps(index_field_dict['parameters']), 'seeds':json.dumps(request.session['seed_list']), 'prog_type':prog_type, 'username':request.user})
    else:
        
        print request.session.keys()
        
        request_post = dict(request.POST.iterlists())
        
        #get ride of the token...
        request_post.pop("csrfmiddlewaretoken")
        
        inp_params = json.loads(request_post['parameters'][0])
        trans_funcs={}
        
        #The user should only have modified the values and/or the groups of the fields, the 'trans' fields etc should not have been altered from config
        for param,param_val in config.adjust_fields.items():
          for field, field_val in param_val['fields'].items():
              if field_val.has_key('trans'):
                  trans_funcs[field] = copy.deepcopy(field_val['trans'])
        
        #the non-grouped vars are added to the session here and the grouped vars are returned as the where_vars list
        where_vars = core.parse_parameters(inp_params, trans_funcs, request)
        
        #save the where_vars list to the session
        request.session['where_vars'] = where_vars
        #though the appropriate variables are added to the session, this is also added for simplicity when a table is re-requested.
        request.session['inp_params'] = inp_params
        
        #also add in the sample info
        request.session["query_samples"] = core.proper_type(request.POST["query_samples"])
        
        valid_hit_genes = core.get_valid_hits(request, config.hit_session_dict)
        
        if len(valid_hit_genes) > 0:
        
            valid_query_res = core.get_valid_querys(request, config.query_prior_dict)
            
            comb_gene_hits = config.score_hits(request, valid_hit_genes)
            
            conv_prot_hits = config.convert_ids_to(request, comb_gene_hits)
            
            prior_prot_hits = config.prioritization_func['function'](request, conv_prot_hits, **config.prioritization_func['args'])
            
            #subset to only those genes which are seeds or queries
            
            limit_list = comb_gene_hits.nodeList().ids()
            query_only_ids = map(lambda x: x.getAttr(['attributes', 'gene']),valid_query_res)
            
            limit_list.extend(query_only_ids)
            
            query_results = config.convert_ids_from(request, prior_prot_hits, limit_list)
            
            seed_list = config.get_seed_list(comb_gene_hits)
            
            temp_query = copy.deepcopy(query_results)
            temp_query.subset(query_only_ids)
            
            query_list = config.get_query_list(temp_query)
            
            #make the results table
            
            rels, rels_header = core.make_results_table(valid_query_res, valid_hit_genes, query_results)
            
            
        else:
            
            rels = []
            rels_header = []
            seed_list = []
            query_list = []
            query_results = core.SeedList(core.NodeList())
            
        request.session['hitwalker_score'] = query_results
        request.session['seed_list'] = seed_list
        request.session['query_list'] = query_list
        request.session['rels'] = rels
        request.session['rels_header'] = rels_header
        
        #prepare cur_filts and cur_param to be ouput as part of the table
        
        cur_param = {}
        
        for key,val in inp_params.items():
            if val['type'] == 'standard':
                cur_param[key] = copy.deepcopy(val['fields'])
        
        cur_filts = {}
        
        for i in where_vars:
            cur_filts[i['name']] = i['pretty_where']
            
        return render(request, 'network/res_table.html', {'res_table':rels, 'res_header':rels_header, 'sample_display_name':request.session['query_samples']['SampleID'],
                                                    'cur_filts':json.dumps(cur_filts), 'cur_param':json.dumps(cur_param), 'seeds':json.dumps(seed_list[:50]), 'queries':json.dumps(query_list[:50]), 'prog_type':prog_type,
                                                    'username':request.user})
            
        
        
        
        

#@login_required(login_url=prog_type+'/HitWalker2/login/')
#def networkWaiting(request):
#    
#    request_post = dict(request.POST.iterlists())
#    request_post['prog_type'] = prog_type
#    ret_post = core.request_post_to_json(request_post)
#    
#    print ret_post
#    
#    #need to feed the values from this request into network_waiting which will then invisibly submit them to network and redirect to network
#    return render(request, 'network/network_waiting.html', ret_post)

#ajax related-views...

def save_parameters(request):
    #to be consistent with the posted data processed in views, need to make the value of each element be a list
    
    save_name = request.POST['save_name']
    param = request.POST['adjust_fields']
    
    print param
    
    #filter_hs = json.loads(request.POST['filter_hs'])
    #param_hs = json.loads(request.POST['param_hs'])
    #var_logic = json.loads(request.POST['var_logic'])
    #
    #for i in filter_hs.keys():
    #    if isinstance(filter_hs[i], list) == False:
    #        filter_hs[i] = [filter_hs[i]]
    #
    #num_hs, group_hs, necessary_vars = core.process_request_post (filter_hs)
    #new_filt_dict = core.make_updated_filter_dict (num_hs, group_hs, var_logic)
    #
    exist_param = user_parameters.objects.filter(user=request.user, name=save_name)
    
    if len(exist_param) > 0:
        exist_param.delete();
    
    #cur_db = user_parameters(user=request.user, name=save_name, filt=json.dumps(new_filt_dict), param=json.dumps(param_hs))
    cur_db = user_parameters(user=request.user, name=save_name, param=param)
    cur_db.save()
    
    return HttpResponse(json.dumps({"save_name":save_name}),  mimetype="application/json")

def load_parameters(request):
    
    load_name = request.POST['load_name']
    
    exist_param = user_parameters.objects.filter(user=request.user, name=load_name)
    
    if len(exist_param) != 1:
        raise Exception("ERROR: only expected one database entry")
    
    #return HttpResponse(json.dumps({'filters':exist_param[0].filt, 'parameters':exist_param[0].param}),  mimetype="application/json")
    return HttpResponse(json.dumps({'parameters':exist_param[0].param}),  mimetype="application/json")
    
def get_match(request, match_type):
    
    query, query_list = config.matchers[match_type](request.POST['query'])
    
    return HttpResponse(json.dumps({"query":query, "data":{"results":query_list}}), mimetype="application/json")


#as input, it needs the requested ID denoted 'cur_id'
#the supplied cypher query needs to be supplied this ID and return 
#if config.sample_rels_type == 'hierarchical'
    #the first n columns need to be strings with headers to display with one denoted as 'Sample' which contains the sample ids
    #the next m are integers representing counts of the number of valid records and the last column should be called required_data and needs to be
    #1 or 0 indicating whether the sample is valid to be used for prioritization or not.
#else;
#the query needs 

def get_sample_rels(request):
    
    graph_db = neo4j.GraphDatabaseService()
    
    cur_id = request.POST['cur_id']
    
    sample_query = neo4j.CypherQuery(graph_db,config.sample_rels_query.replace("$$sample$$", cur_id))
    
    sample_list=[]
    
    if config.sample_rels_type == 'hierarchical':
    
        for i in sample_query.execute().data:
            temp_dict = {}
            temp_dict['property_order'] = []
            
            for j_ind,j in enumerate(i.values):
                if isinstance(j, int) and i.columns[j_ind] != "required_data":
                    if temp_dict.has_key('data'):
                        temp_dict['data'].append({'name':i.columns[j_ind], 'count':j})
                    else:
                        temp_dict['data'] = [{'name':i.columns[j_ind], 'count':j}]
                elif isinstance(j, int) and i.columns[j_ind] == "required_data":
                    temp_dict["required_data"] = j
                else:
                    temp_dict[i.columns[j_ind]] = j#.get_properties()
                    temp_dict['property_order'].append(i.columns[j_ind])
                    
            sample_list.append(temp_dict)
    else:
        raise Exception("sample_rels_type has not been implemented")
    
    use_sample_list = sorted(sample_list, key=lambda x: x['required_data'], reverse=True)
    
    return HttpResponse(json.dumps({"sample_list":use_sample_list, "cur_id":cur_id}), mimetype="application/json")

#is overkill for this function, add in below in multi_node_query
def node_query(request):
    
    ret_nodes = json.loads(request.POST["nodes"])
    
    ret_context = request.POST["context"]
    
    unique_types = set(map(lambda x: str(x["attributes"]["node_type"]),ret_nodes))
    
    ret_content = ''
    
    if len(unique_types) == 1:
        cur_type = unique_types.pop()
        ret_content = config.node_content[cur_type]["func"](ret_nodes[0], ret_context, **config.node_content[cur_type]['args'])
    
    return HttpResponse(json.dumps({"content":ret_content}), mimetype="application/json")

def multi_node_query(request):
    ret_nodes = json.loads(request.POST["nodes"])
    
    unique_types = set(map(lambda x: str(x['attributes']['node_type']),ret_nodes))
    
    num_metas = 0
    
    if 'MetaNode' in unique_types:
        num_metas += 1
        unique_types.remove('MetaNode')
        for i in ret_nodes:
            if i['attributes']['node_type'] == "MetaNode":
                for j in i['children']:
                    unique_types.add(str(j['attributes']['node_type']))
    
    unique_type_list = sorted(unique_types, key=str.lower)
    
    cur_dict = core.iterate_dict(config.node_group_content, unique_type_list)
    
    ungroup_dis_text = "disabled"
    pathway_dis_text = "disabled"
    
    if num_metas == len(unique_types):
        ungroup_dis_text = ""
    
    type_intersect = unique_types.intersection(set(['Sample', 'Gene']))
    
    if  (len(type_intersect) == 1 and type_intersect == set(['Sample'])) or len(type_intersect) == 2:
        pathway_dis_text = ""
    
    group_buttons = '<div class="btn-group"><button type="button" class="btn btn-default" '+ungroup_dis_text+' onclick="ungroup_nodes(this)">Ungroup</button></div>'
    group_buttons += '<div class="btn-group"><button type="button" class="btn btn-default" '+pathway_dis_text+' onclick="pathway_context(this)">Pathway Context</button></div>'
    
    if cur_dict != None:
    
        question_list = map(lambda x: x['text'], cur_dict['options'])
        
        header = '<div class=container style="width:300px">'
        header += '<div class="row"><div class="list-group">'
        
        header += '<a class="list-group-item list-group-item-heading text-center" style="background-color:#f5f5f5">'+cur_dict['title']+'</a>'
        for i_ind, i in enumerate(question_list):
            header += '<a class="list-group-item text-center" onclick="post_to_fullfill(this)" style="cursor:pointer" data-value='+reduce(lambda x,y:x+'.'+y, unique_types)+'-'+str(i_ind)+'>'+i+'</a>'
        
        header += '<br>' + group_buttons + '</div></div>'
    else:
        header = '<div class=container style="width:300px"><p class="text-danger">Sorry, no queries are currently defined for type(s): '+string.joinfields(unique_type_list, ' and ')+'.  Try selecting a different type or only one type at a time.</p><br>' +group_buttons+'</div>'
        
    
    return HttpResponse(json.dumps({"content":header}), mimetype="application/json")

def download(request):
    
    use_file = request.session['tmp_file']
    
    resp_file_name = 'query_result.csv'
    
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="'+resp_file_name+'"'
    
    writer = csv.writer(response)
    
    file_inp = open(use_file, 'r')
    
    #write the parameter info
    
    inp_params = request.session['inp_params']
    where_vars = request.session['where_vars']
    
    for key,val in inp_params.items():
        if val['type'] == 'standard':
            writer.writerow([key])
            for i,j in val['fields'].items():
                writer.writerow([j['name'] + " " + j['comparison'] + " " + str(j['default'])])
    
    for i in where_vars:
        writer.writerow([i['name']])
        write.writerow([i['pretty_where']])
    
    writer.writerow([])
    writer.writerow([])
    
    #write the title
    
    writer.writerow(["Title:" + request.session['tmp_title']])
    writer.writerow([])
    
    for i in csv.reader(file_inp):
        writer.writerow(i)
    
    return response

def fullfill_node_query(request):
    
    node_req = request.POST["choice"]
    ret_node_queries = json.loads(request.POST["nodes"])
    
    #also add in the number of input genes if the queries need it
    
    node_queries = {}
    
    for i in ret_node_queries.keys():
        node_queries[i] = map(lambda x: x['id'], ret_node_queries[i])
        node_queries[i+'_length'] = len(ret_node_queries[i])
    
    #retrieve the query by parsing node_req
    
    req_table = node_req.split('-')
    
    search_list = req_table[0].split('.')
    
    query_type = core.iterate_dict(config.node_group_content, search_list)
    
    query_info = query_type['options'][int(req_table[1])]
    
    if query_info['session_params'] != None:
        for i in query_info['session_params']:
            use_key = i[-1]
            if node_queries.has_key(use_key) == False:
                node_queries[use_key] = core.iterate_dict(request.session, i)
    
    #execute the query
    
    session = cypher.Session()
    tx = session.create_transaction()
    
    use_query = query_info['query']
    
    for var_elem in request.session['where_vars']:
        use_query = core.add_where_input_query(use_query, var_elem['where_statement'], var_elem['necessary_vars'], request.session['graph_struct'])
    
    print use_query, node_queries
    
    #use_query = core.add_where_input_query(query_info['query'], request.session['where_template'], request.session['necessary_vars'], request.session['graph_struct'])
    tx.append(use_query, node_queries)
    res_list = tx.commit()
    
    temp_node_list = []
    
    for i in res_list[0]:
        temp_node_list.append(i.values[0])
    
    if len(ret_node_queries.keys()) == 1 and ret_node_queries.has_key('Gene'):
        if len(ret_node_queries['Gene']) > 3:
            use_title = query_info['title'].replace('$$result$$', ret_node_queries['Gene'][0]['display_name'] + '...' + ret_node_queries['Gene'][-1]['display_name'])
        else:
            use_title = query_info['title'].replace('$$result$$', string.joinfields(map(lambda x: x['display_name'], ret_node_queries['Gene']), ','))
    elif len(ret_node_queries.keys()) == 1 and ret_node_queries.has_key('Sample'):
        if len(ret_node_queries['Sample']) > 3:
            use_title = query_info['title'].replace('$$result$$', ret_node_queries['Sample'][0]['display_name'] + '...' + ret_node_queries['Sample'][-1]['display_name'])
        else:
            use_title = query_info['title'].replace('$$result$$', string.joinfields(map(lambda x: x['display_name'], ret_node_queries['Sample']), ','))
    else:
        use_title = 'ERROR: Unknown title...'
    
    #too many nodes to render efficiently, will allow user to download csv file...
    if len(temp_node_list) > 2000:
        
        subj_nodes = []
        query_nodes = []
        
        for i in ret_node_queries.keys():
            for j in ret_node_queries[i]:
                subj_nodes.append({'id':j['id'], 'node_type':i})
        
        for i in temp_node_list:
            query_nodes.append({'id':i, 'node_type':query_type['returned_node_type']})
        
        node_graph = core.copy_nodes(subj_nodes, query_nodes, request, config.edge_queries, never_group=True)
        
        tmp_file = tempfile.mktemp()
        
        out_p = open(tmp_file, "w")
        
        header = ['id', 'display_name'] + map(lambda x: x['id'], subj_nodes)
        
        writer = csv.writer(out_p)
        
        writer.writerow(header)
        
        for i in node_graph['nodes'].tolist():
            temp_ln = []
            if i['attributes']['node_type'] == query_type['returned_node_type']:
                temp_ln.extend([i['id'], i['display_name']] + [None]*(len(header)-2))
                for j in i['children']:
                    if j['attributes']['node_type'].replace('_Hit', '') == query_info['text']:
                        for k_ind, k in enumerate(j['attributes']['other_nodes']):
                            val_loc = header.index(k)
                            if val_loc == -1:
                                raise Exception("unexpected node found")
                            else:
                                if j['attributes']['meta'].has_key('score'):
                                    temp_ln[val_loc] = j['attributes']['meta']['score'][k_ind]
                                else:
                                    temp_ln[val_loc] = 1
                writer.writerow(temp_ln)
                            
        out_p.close()
        
        request.session['tmp_file'] = tmp_file
        request.session['tmp_title'] = use_title
        
        ret_dict = {'is_graph':False, 'graph':{}, 'title':'Sorry, your query is too large to display.  However, you may download a text version.'}
        
    else:
        
        #convert the IDs to nodes
        
        ret_nodes = core.get_nodes(temp_node_list, query_type['returned_node_type'], request,  missing_param="skip")
    
        ret_nodes = core.apply_grouping2({'nodes':ret_nodes, 'links':[]}, [])['nodes']
    
        ret_dict = {'is_graph':True, 'graph':{'nodes':ret_nodes.tolist(), 'links':[]}, 'title':use_title}
    
    return HttpResponse(json.dumps(ret_dict),mimetype="application/json")

def get_data (request):
    
    request_post = core.request_post_to_json(dict(request.POST.iterlists()))
    
    data_dict, title = config.plot_data_types[request_post["type"]](request_post)
    
    return HttpResponse(json.dumps({'graph':data_dict, 'title':title}),mimetype="application/json")

def copy_nodes(request):
    subj_nodes = json.loads(request.POST["subj"])
    query_nodes = json.loads(request.POST["query"])
    
    cur_graph = core.copy_nodes(subj_nodes, query_nodes, request, config.edge_queries)
    
    ret_dict = {'nodes':cur_graph['nodes'].tolist(), 'links':cur_graph['links']}
    
    return HttpResponse(json.dumps(ret_dict),mimetype="application/json")

def get_graph(request):
    
    request_post = core.fix_jquery_array_keys(dict(request.POST.iterlists()))
    
    nodes, links, title = core.iterate_dict(config.graph_initializers, request_post['panel_context'])(request, request_post)
    
    ret_dict = {'nodes':nodes, 'links':links, 'title':title}
    
    return HttpResponse(json.dumps(ret_dict),mimetype="application/json")

@ensure_csrf_cookie
@login_required(login_url=prog_type+'/HitWalker2/login/')
def pathway(request):
    
    if len(request.POST) == 0:
       
       if socket.gethostname() in set(["HRCC448", "HRCC255"]):
           pickle_inp = open("/var/www/hitwalker_2_inst/test_pathway_request.json", "r")
           ret_json = json.load(pickle_inp)
       return render(request, 'network/d3_test.html', ret_json)
    else:
        
        default_css_classes, new_css_path, node_type_transl, edge_type_transl = generate_css(str(request.user))
        request.session['new_css_path'] = new_css_path
        
        request_post = dict(request.POST.iterlists())
        
        #as multiple samples can be specified and they will be sent ',' delimited by the JS:
        input_vals = {'sample_name':request.POST['sample_name'].split(','), 'pathway_name':request.POST['pathway_name']}
        
        inp_params = request.session['inp_params']
        where_vars = request.session['where_vars']
        
        cur_param = {}
        
        for key,val in inp_params.items():
            if val['type'] == 'standard':
                cur_param[key] = copy.deepcopy(val['fields'])
        
        cur_filts = {}
        
        for i in where_vars:
            cur_filts[i['name']] = i['pretty_where']
        
        print input_vals
        
        ret_json = {'prog_type':str(prog_type), 'username':str(request.user), 'node_type_transl':json.dumps(node_type_transl), 'edge_type_transl':json.dumps(edge_type_transl),
                    'metanode_thresh':config.max_nodes, 'panel_context':'image', 'input_vals':json.dumps(input_vals), 'cur_param':json.dumps(cur_param), 'cur_filts':json.dumps(cur_filts),
                    'css_path':os.path.basename(new_css_path), 'default_css':json.dumps(default_css_classes)}
        
        for key, val in config.pathway_sizes.items():
            if ret_json.has_key(key) == False:
                ret_json[key] = val;
            else:
                raise Exception("Duplicate entry for " + key)
        
        if socket.gethostname() in set(["HRCC448", "HRCC255"]):
                pickle_outp = open("/var/www/hitwalker_2_inst/test_pathway_request.json", "w")
                json.dump(ret_json, pickle_outp)
        
        return render(request, 'network/d3_test.html', ret_json)

@ensure_csrf_cookie
@login_required(login_url=prog_type+'/HitWalker2/login/')
def network(request):
    
    if len(request.POST) == 0:
        
        if socket.gethostname() in set(["HRCC448", "HRCC255"]):
            pickle_inp = open("/var/www/hitwalker_2_inst/test_network_request.json", "r")
            ret_json = json.load(pickle_inp)
        return render(request, 'network/d3_test.html', ret_json)
        #return redirect(prog_type+'/HitWalker2/')
    elif request.POST.has_key("output_format") and request.POST.has_key("data"):
        
        try:
            import cairosvg
            #and all the somewhat invisible dependencies needed to use external css styling on the svgs...
            import lxml
            import tinycss
            import cssselect
        
            #from http://stackoverflow.com/questions/20001282/django-how-to-reference-paths-of-static-files-and-can-i-use-them-in-models?rq=1
            #static_storage = get_storage_class(settings.STATICFILES_STORAGE)()
            
            #cur_dir = tempfile.tempdir
            
            #css_path = os.path.relpath(os.path.join(static_storage.location, "network/static/network/css/HitWalker2.css"), cur_dir)
           
            #css_path =os.path.join(static_storage.location, "network/static/network/css/HitWalker2.css")
            
            #use_prog_type = prog_type.strip('/')
            #if use_prog_type != "":
            #    use_prog_type = "_" + use_prog_type
            #    css_path = os.path.join("/var/www/hitwalker_2_inst" + use_prog_type, "HitWalker2/network/static/network/css/HitWalker2.css")
            #else:
            #css_path =os.path.join(static_storage.location, "network/static/network/css/HitWalker2.css")
            
            if request.POST["plot_type"] == 'svg':
                
                use_svg = '<?xml-stylesheet type="text/css" href="' + request.session['new_css_path'] + '" ?>' + request.POST["data"]
                
            else:
                if request.POST["plot_type"] == 'g1':
                    width="600"
                    height="400"
                    transl = 'translate(200,0)'
                elif request.POST["plot_type"] == 'siRNA_plot':
                    width="850"
                    height="400"
                    transl='translate(50,0)'
            
                use_data = '<svg xmlns="http://www.w3.org/2000/svg" width="'+width+'" height="'+height+'" >' + request.POST["data"].replace('xmlns="http://www.w3.org/2000/svg"', "")  + '</svg>'
                
                import xml.etree.ElementTree as ET
                tree = ET.fromstring(use_data)
                
                #ensure everything will fit on screen
                #assuming only one non-legend g for now
                tree.find("./*[@class=\'"+request.POST["plot_type"]+"\']").attrib['transform'] = transl
                
                tree.find(".//*[@class='BorderRect BorderSelected']/..").remove(tree.find(".//*[@class='BorderRect BorderSelected']"))
                
                use_svg = '<?xml-stylesheet type="text/css" href="' + request.session['new_css_path'] + '" ?>' + ET.tostring(tree)
                
            temp_file = tempfile.mktemp()
            
            #temp_file = '/Users/bottomly/Desktop/test_image.xml'
            
            temp_inp = open(temp_file, "w")
            
            temp_inp.write(use_svg)
            temp_inp.close()
            
            
            if request.POST["output_format"] == "svg":
                response = HttpResponse(cairosvg.svg2svg(url=temp_file), content_type='image/svg+xml')
                response['Content-Disposition'] = 'attachment; filename="HitWalker2.svg"'
            elif request.POST["output_format"] == "pdf":
                response = HttpResponse(cairosvg.svg2pdf(url=temp_file), content_type='application/pdf')
                response['Content-Disposition'] = 'attachment; filename="HitWalker2.pdf"'
            
            return response
        except ImportError, e:
            response = HttpResponse('<?xml-stylesheet type="text/css" ?>' + request.POST["data"], content_type='image/svg+xml')
            response['Content-Disposition'] = 'attachment; filename="HitWalker2.svg"'
            return response
        
    else:
        
        request_post = dict(request.POST.iterlists())
        
        default_css_classes, new_css_path, node_type_transl, edge_type_transl = generate_css(str(request.user))
        
        request.session['new_css_path'] = new_css_path
        
        if request_post.has_key('query_samples'):
        
            #if bypassing table I think this would be necessary...
            ##get ride of the token...
            
            request_post = dict(request.POST.iterlists())
            
            #get ride of the token...
            request_post.pop("csrfmiddlewaretoken")
            
            inp_params = json.loads(request_post['parameters'][0])
            trans_funcs={}
            
            #The user should only have modified the values and/or the groups of the fields, the 'trans' fields etc should not have been altered from config
            for param,param_val in config.adjust_fields.items():
              for field, field_val in param_val['fields'].items():
                  if field_val.has_key('trans'):
                      trans_funcs[field] = copy.deepcopy(field_val['trans'])
            
            #the non-grouped vars are added to the session here and the grouped vars are returned as the where_vars list
            where_vars = core.parse_parameters(inp_params, trans_funcs, request)
            
            #save the where_vars list to the session
            request.session['where_vars'] = where_vars
            
            #also add in the sample info
            request.session["query_samples"] = core.proper_type(request.POST["query_samples"])
            request.session["inp_params"] = inp_params
            
            input_vals = {'query_samples':request.session['query_samples']}
        else:
            input_vals = {'var_select':request_post['var_select'], 'seed_select':request_post['seed_select']}
        
        
        inp_params = request.session['inp_params']
        where_vars = request.session['where_vars']
        
        cur_param = {}
        
        for key,val in inp_params.items():
            if val['type'] == 'standard':
                cur_param[key] = copy.deepcopy(val['fields'])
        
        cur_filts = {}
        
        for i in where_vars:
            cur_filts[i['name']] = i['pretty_where']
        
        
        ret_json = {'prog_type':str(prog_type), 'username':str(request.user), 'node_type_transl':json.dumps(node_type_transl), 'edge_type_transl':json.dumps(edge_type_transl),
                    'metanode_thresh':config.max_nodes, 'panel_context':'panel', 'input_vals':json.dumps(input_vals), 'cur_param':json.dumps(cur_param), 'cur_filts':json.dumps(cur_filts),
                    'css_path':os.path.basename(new_css_path),  'default_css':json.dumps(default_css_classes)}
        
        for key, val in config.network_sizes.items():
            if ret_json.has_key(key) == False:
                ret_json[key] = val;
            else:
                raise Exception("Duplicate entry for " + key)
        
        if socket.gethostname() in set(["HRCC448", "HRCC255"]):
            pickle_outp = open("/var/www/hitwalker_2_inst/test_network_request.json", "w")
            json.dump(ret_json, pickle_outp)
        
        return render(request, 'network/d3_test.html', ret_json)