# Create your views here.
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseServerError
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
import network.models as netmods
from network.models import user_parameters
from django.db.models import Q
from django.db.models.query import QuerySet

import cssutils
import colour

use_css = ''

prog_type = config.prog_type

if prog_type != "":
    prog_type = "/" + prog_type

static_base = '/var/www/hitwalker2_inst'+prog_type+'/static/'

try:
    graph_inp = open(config.graph_struct_file, "r")
    graph_struct = json.load(graph_inp)
    graph_inp.close()
except:
    graph_struct = {}

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
    
    if os.path.isfile(static_base + 'network/css/HitWalker2.css') == False:
        #from http://stackoverflow.com/questions/20001282/django-how-to-reference-paths-of-static-files-and-can-i-use-them-in-models?rq=1
        static_storage = get_storage_class(settings.STATICFILES_STORAGE)()
        css_base = static_storage.location
    else:
        css_base = static_base[:]
    
    css_path =os.path.join(css_base, "network/css/HitWalker2.css")
    hw_css = cssutils.parseFile(css_path)
    
    node_type_transl={'Gene':{'name':'Gene', 'class':'Gene'}, 'Subject':{'name':'Subject', 'class':'Subject'}}
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
    
    for node in map(lambda x: x+'_Hit', temp_hits) + [config.data_types['target']]:
        if node_type_transl.has_key(node) == False:
            #choose a default for it
            new_color = default_css_classes.difference(used_defaults).pop()
            node_type_transl[node] = {'name':node, 'class':new_color.replace('.', '')}
            used_defaults.add(new_color)
            
            important_defaults.add(new_color)
    
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
    new_css_path = os.path.join(css_base, "network/css/HitWalker2_"+user_name+".css")
    new_css_out = open(new_css_path, "w")
    new_css_out.write(hw_css.cssText)
    new_css_out.close()
    
    default_css_list = map(lambda x: x.replace(".", ""), default_css_classes.difference(important_defaults))
    
    return default_css_list, new_css_path, node_type_transl, edge_type_transl
    

@sensitive_post_parameters()
@never_cache
@login_required(login_url='/HitWalker2' + prog_type +'/login/')
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
            
            return redirect('/HitWalker2'+prog_type)
    else:
        form = SetPasswordForm(request.user)

    return render(request, 'network/password.html', {'form':form, 'prog_type':prog_type, 'username':request.user})

def get_default_parameters(request):
   
    try:
        
        index_field_dict = copy.deepcopy(config.adjust_fields)
        for i in index_field_dict.keys():
            for j in index_field_dict[i]['fields'].keys():
                if index_field_dict[i]['fields'][j].has_key('trans'):
                    index_field_dict[i]['fields'][j].pop('trans')
            if index_field_dict[i]['type'] == 'grouped':
                logical_list = reduce (lambda x,y: x+y, reduce (lambda x,y: x+y,index_field_dict[i]['default_groups']))
                index_field_dict[i]['logical_list'] = map(lambda x: 'null' if x.has_key('logical')==False else x['logical'], logical_list)
        
        return HttpResponse(json.dumps({'adjust_fields':index_field_dict}),  mimetype="application/json")
    
    except:
        
        return HttpResponseServerError()

@login_required(login_url='/HitWalker2'+prog_type+'/login/')
def index(request, retry_message):
    
    if request.session.has_key('inp_params'):
        index_field_dict = request.session['inp_params']
    else:
        
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
    
    context = {'adjust_fields':json.dumps(index_field_dict), 'param_names':json.dumps(param_names), 'prog_type':prog_type, 'username':request.user, 'data_types':config.data_types}
    
    return render(request, 'network/index.html', context)

@login_required(login_url='/HitWalker2'+prog_type+'/login/')
def table(request):
    
    if len(request.POST) == 0:
        return redirect('/HitWalker2'+prog_type)
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
    else:
        
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
            
            rels, rels_header = core.make_results_table(valid_query_res, valid_hit_genes, query_results, query_list)
            
            
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
            
#ajax related-views...

def save_parameters(request):
    #to be consistent with the posted data processed in views, need to make the value of each element be a list
    
    try:
        
        save_name = request.POST['save_name']
        param = request.POST['adjust_fields']
        
        exist_param = user_parameters.objects.filter(user=request.user, name=save_name)
        
        if len(exist_param) > 0:
            exist_param.delete();
        
        #cur_db = user_parameters(user=request.user, name=save_name, filt=json.dumps(new_filt_dict), param=json.dumps(param_hs))
        cur_db = user_parameters(user=request.user, name=save_name, param=param)
        cur_db.save()
        
        return HttpResponse(json.dumps({"save_name":save_name}),  mimetype="application/json")
    except:
        return HttpResponseServerError()

def load_parameters(request):
    
    try:
        
        load_name = request.POST['load_name']
        
        exist_param = user_parameters.objects.filter(user=request.user, name=load_name)
        
        if len(exist_param) != 1:
            raise Exception("ERROR: only expected one database entry")
        
        #return HttpResponse(json.dumps({'filters':exist_param[0].filt, 'parameters':exist_param[0].param}),  mimetype="application/json")
        return HttpResponse(json.dumps({'parameters':exist_param[0].param}),  mimetype="application/json")
    except:
        return HttpResponseServerError()
    
def get_match(request, match_type):
    
    try:
        query, query_list = config.matchers[match_type](request.POST['query'])
        
        return HttpResponse(json.dumps({"query":query, "data":{"results":query_list}}), mimetype="application/json")
    
    except:
        
        return HttpResponseServerError()


#as input, it needs the requested ID denoted 'cur_id'
#the supplied cypher query needs to be supplied this ID and return 
#if config.sample_rels_type == 'hierarchical'
    #the first n columns need to be strings with headers to display with one denoted as 'Sample' which contains the sample ids
    #the next m are integers representing counts of the number of valid records and the last column should be called required_data and needs to be
    #1 or 0 indicating whether the sample is valid to be used for prioritization or not.
#else;
#the query needs 

def get_sample_rels(request):
    
    try:
        
        graph_db = neo4j.GraphDatabaseService(config.cypher_session+'/db/data/')
        
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
    except:
        
        return HttpResponseServerError()

#is overkill for this function, add in below in multi_node_query
def node_query(request):
    
    try:
    
        ret_nodes = json.loads(request.POST["nodes"])
        
        ret_context = request.POST["context"]
        
        unique_types = set(map(lambda x: str(x["attributes"]["node_type"]),ret_nodes))
        
        ret_content = ''
        
        if len(unique_types) == 1:
            cur_type = unique_types.pop()
            ret_content = config.node_content[cur_type]["func"](ret_nodes[0], ret_context, **config.node_content[cur_type]['args'])
        
        return HttpResponse(json.dumps({"content":ret_content}), mimetype="application/json")
        
    except:
        return HttpResponseServerError()

def multi_node_query(request):
    
    try:
    
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
        
        type_intersect = unique_types.intersection(set(['Subject', 'Gene']))
        
        if  (len(type_intersect) == 1 and type_intersect == set(['Subject'])) or len(type_intersect) == 2:
            pathway_dis_text = ""
        
        #group_buttons = '<div class="btn-group"><button type="button" class="btn btn-default" '+ungroup_dis_text+' onclick="ungroup_nodes(this)">Ungroup</button></div>'
        group_buttons = '<div class="btn-group"><button type="button" class="btn btn-default" '+pathway_dis_text+' onclick="pathway_context(this)">Pathway Context</button></div>'
        
        if cur_dict != None:
        
            question_list = map(lambda x: x['text'], cur_dict['options'])
            
            header = '<div class="row"><div class="list-group">'
            header += '<a class="list-group-item list-group-item-heading text-center" style="background-color:#f5f5f5">'+cur_dict['title']+'</a>'
            
            header += '<div class="panel-group" id="pg1" role="tablist" aria-multiselectable="true">'
            
            for i_ind, i in enumerate(question_list):
                header += '''
                    <div id="ppd$id$" class="panel panel-default">
                      <div class="panel-heading" role="tab" id="heading$id$" position="relative">
                        <h4 class="panel-title">
                          <a data-toggle="false" onclick="post_to_fullfill(this)"  style="cursor:pointer" data-parent="#pg1" data-value=$data$ aria-expanded="false" aria-controls="collapse$id$">$data_title$</a>
                        </h4>
                      </div>
                    <div id="collapse$id$" class="panel-collapse collapse" role="tabpanel" aria-labelledby="heading$id$">
                        <div class="panel-body"></div>
                    </div>
                </div>
                '''.replace("$data$", reduce(lambda x,y:x+'.'+y, unique_types)+'-'+str(i_ind)).replace("$id$", str(i_ind)).replace("$data_title$", i)
                
            
            #header = '<div class=container style="width:300px">'
            #header += '<div class="row"><div class="list-group">'
            #
            #header += '<a class="list-group-item list-group-item-heading text-center" style="background-color:#f5f5f5">'+cur_dict['title']+'</a>'
            #for i_ind, i in enumerate(question_list):
            #    header += '<a class="list-group-item text-center" onclick="post_to_fullfill(this)" style="cursor:pointer" data-value='+reduce(lambda x,y:x+'.'+y, unique_types)+'-'+str(i_ind)+'>'+i+'</a>'
            #
            header += '<br>' + group_buttons + '</div></div></div>'
        else:
            header = '<div class=container style="width:300px"><p class="text-danger">Sorry, no queries are currently defined for type(s): '+string.joinfields(unique_type_list, ' and ')+'.  Try selecting a different type or only one type at a time.</p><br>' +group_buttons+'</div>'
            
        
        return HttpResponse(json.dumps({"content":header}), mimetype="application/json")
    except Exception as e:
        print e
        return HttpResponseServerError()


def download(request, title, var_type, query_set, query_info, return_set, returned_node_type):
    
    resp_file_name = 'query_result.csv'
    
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="'+resp_file_name+'"'
    
    writer = csv.writer(response)
    
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
        writer.writerow([i['pretty_where']])
    
    writer.writerow([])
    writer.writerow([])
    
    #write the title
    
    writer.writerow(["Title:" + title])
    writer.writerow([])
    
    #start writing the results
    
    #make this all only for dense datatypes
    
    #can sidestep get_nodes as we will already have the use_query QuerySet, can execute direcly and run the handler
    import network.models
    cur_mod = getattr(network.models, var_type)
    header = map(lambda x: x.name, cur_mod._meta.fields)
    
    table_header = ['Symbol']
    table_header.extend(header)
    print table_header
    writer.writerow(table_header)
    
    #get the gene symbols
    
    #further filter the QuerySet by the actual specified genes/subject (just in case)
    
    all_qs = Q()
    
    for i in query_info.items():
        all_qs = all_qs & Q(**{i[0].lower()+"__in":i[1]})
    
    query_res = query_set.filter(all_qs).values_list()
    
    gene_list = set()
    fin_table = []
    for i in query_res:
        if str(i[header.index(returned_node_type.lower())]) in return_set:
            fin_table.append(list(i))
            gene_list.add(i[header.index('gene')])
        
    temp_nl = core.NodeList()
    
    for i in gene_list:
        temp_nl.add(core.BasicNode([i]))
    
    temp_sl = core.SeedList(temp_nl)
    
    symb_list = config.get_query_list(temp_sl)
    
    print symb_list[:5]
    
    symb_dict = collections.defaultdict(list)
    
    for i in symb_list:
        symb_dict[i['gene']].append(i['symbol'])
    
    for i in fin_table:
        
        if symb_dict.has_key(i[header.index('gene')]):
            i_tab = [string.joinfields(symb_dict[i[header.index('gene')]], ',')]
        else:
            i_tab = [i[header.index('gene')]]
        
        i_tab.extend(i)
        
        writer.writerow(i_tab)
    
    return response

def fullfill_node_query(request):
    
    try:
    
        node_req = request.POST["choice"]
        ret_node_queries = json.loads(request.POST["nodes"])
        
        #also add in the number of input genes if the queries need it
        
        node_queries = {}
        
        cur_len = 0
        
        for i in ret_node_queries.keys():
            node_queries[i] = map(lambda x: x['id'], ret_node_queries[i])
            cur_len = len(ret_node_queries[i])
        
        #retrieve the query by parsing node_req
        
        req_table = node_req.split('-')
        
        search_list = req_table[0].split('.')
        
        query_type = core.iterate_dict(config.node_group_content, search_list)
        
        query_info = query_type['options'][int(req_table[1])]
        
        #execute the query
        
        temp_node_list = []
        max_vals = []
        use_query = query_info['query']
        
        if (query_info.has_key('db_type') == False) or (query_info.has_key('db_type') and query_info['db_type'] == 'neo4j'):
        
            #this can be local to this function as this info can be accessed through the session info in core.handle_dense_gene_hits
            if query_info['session_params'] != None:
                for i in query_info['session_params']:
                    use_key = i[-1]
                    if node_queries.has_key(use_key) == False:
                        node_queries[use_key] = core.iterate_dict(request.session, i)
        
            graph_db = neo4j.GraphDatabaseService(config.cypher_session+'/db/data/')
            
            #neo4j grouping filters are now depricated...
            #for var_elem in request.session['where_vars']:
            #    use_query = core.add_where_input_query(use_query, var_elem['where_statement'], var_elem['necessary_vars'], request.session['graph_struct'])
            
            #print use_query, node_queries
            
            #use_query = core.add_where_input_query(query_info['query'], request.session['where_template'], request.session['necessary_vars'], request.session['graph_struct'])
            query = neo4j.CypherQuery(graph_db, use_query)
            
            
            for i in query.stream(**node_queries):
                #print i.values
                if (len(max_vals) == 0) or ((i.values[1] < max_vals[-1]) and (len(max_vals) < 5)):
                    max_vals.append(i.values[1])
                    temp_node_list.append((str(int(round((i.values[1]/float(cur_len))*100))), i.values[0]))
                elif i.values[1] >= max_vals[-1]:
                    temp_node_list.append((str(int(round((i.values[1]/float(cur_len))*100))), i.values[0]))
                else:
                    break
                
        elif query_info.has_key('db_type') and query_info['db_type'] == 'sql':
            
            print node_queries

            query_keys = node_queries.keys()

            if len(query_keys) == 1:
                node_type = query_keys[0]
            else:
                raise Exception("Can't handle multiple keys for node queries for now...")
            
            #get database
            cur_mod = getattr(netmods, query_info['text'])
            
            header = map(lambda x: x.name, cur_mod._meta.fields)
            
            comb_q = Q()
            
            if query_info['session_params'] != None:
                #have to also get the directional info from the parameter info
                #note that query_info['session_params'] and obj_name should be the same value after unlisting of the former...
                
                for i in query_info['session_params']:
                    
                    obj_name = i[:].pop().lower()
                    if request.session.has_key(obj_name):
                            comp = request.session['inp_params']['Seed_Parameters']['fields'][obj_name]['comparison']
                            db_comp = 'eq'
                            if comp == '>':
                                db_comp = 'gt'
                            elif comp == '<':
                                db_comp = 'lt'
                            
                            comb_q = comb_q & Q(**{'score__'+db_comp:core.iterate_dict(request.session, i)})
                            
            
            comb_q = comb_q & Q(**{node_type.lower()+"__in":node_queries[node_type]})
            
            if request.session.has_key('where_vars') and len(request.session['where_vars']) > 0 and len(request.session['where_vars'][0]['necessary_vars'].difference(set(header)))==0:
                            
                db_res = cur_mod.objects.using("data").filter(comb_q).extra(where=[request.session['where_vars'][0]['where_statement'].replace("$$$$.", "")])
            
            elif request.session.has_key('where_vars') and len(request.session['where_vars']) > 1:
                raise Exception("Not currently expecting multiple grouped filters")
            else:
                db_res = cur_mod.objects.using("data").filter(comb_q)
            
            res_set_dict = collections.defaultdict(set)
            
            key_col = node_type.lower()
            append_col = query_type['returned_node_type'].lower()
            
            use_query=db_res
            
            for i in db_res:
                #get the current field as well as the other one based on the returned_node field of query_info (or parent object)
                res_set_dict[str(getattr(i, key_col))].add(str(getattr(i, append_col)))
            
            if len(res_set_dict.keys()) == 0:
                res_col = collections.Counter()
            elif len(res_set_dict.keys()) == 1:
                res_col = collections.Counter(res_set_dict.items()[0][1])
            else:
                
                temp_dict = []
                
                for i in res_set_dict.values():
                    for j in i:
                        temp_dict.append(j)
                    
                res_col = collections.Counter(temp_dict)
            
            for i in res_col.most_common():
                if (len(max_vals) == 0) or ((i[1] < max_vals[-1]) and (len(max_vals) < 5)):
                        max_vals.append(i[1])
                        temp_node_list.append((str(int(round((i[1]/float(cur_len))*100))), i[0]))
                elif i[1] >= max_vals[-1]:
                        temp_node_list.append((str(int(round((i[1]/float(cur_len))*100))), i[0]))
                else:
                        break

        else:
            raise Exception("Unknown db_type")
            
        ret_dict = collections.defaultdict(list)
        
        for k,v in temp_node_list:
            ret_dict[k].append(v)
        
        count_coll = map(lambda x: (x[0], len(x[1])), ret_dict.items())
        
        if request.session.has_key('tmp_query_results') == False:
        
            request.session['tmp_query_results'] = {str(req_table[1]):{'ret_dict':ret_dict, 'ret_node_queries':ret_node_queries,
                                                                       'query_info':query_info, 'use_query':use_query, 'node_queries':node_queries,
                                                                       }}
        else:
            request.session['tmp_query_results'][str(req_table[1])] = {'ret_dict':ret_dict, 'ret_node_queries':ret_node_queries,
                                                                       'query_info':query_info, 'use_query':use_query, 'node_queries':node_queries}
        
        #this is necessary as the keys to the session were not modified
        request.session.modified = True
        
        return HttpResponse(json.dumps({'ret_node_type':query_type['returned_node_type'], 'results':dict(count_coll)}),mimetype="application/json")
    
    except Exception as e:
        print e
        return HttpResponseServerError()


def provide_data_for_request(request):

    try:
        
        query_choice = request.POST["query_choice"]
        frequency_choice = request.POST["freq_choice"]
        returned_node_type = request.POST["ret_node_type"]
        display_type = request.POST["display_type"]
        
        session_data = request.session['tmp_query_results'][str(query_choice)]
    
        ret_node_queries = session_data['ret_node_queries']
        query_info = session_data['query_info']
        use_query = session_data['use_query']
        node_queries = session_data['node_queries']
        
        #temp_node_list will be formed from whichever frequency the user choses...
        
        temp_node_list = session_data['ret_dict'][frequency_choice]
        
        #we also need node_queries and query_info
        
        if len(ret_node_queries.keys()) == 1 and ret_node_queries.has_key('Gene'):
            if len(ret_node_queries['Gene']) > 3:
                use_title = query_info['title'].replace('$$result$$', ret_node_queries['Gene'][0]['display_name'] + '...' + ret_node_queries['Gene'][-1]['display_name'])
            else:
                use_title = query_info['title'].replace('$$result$$', string.joinfields(map(lambda x: x['display_name'], ret_node_queries['Gene']), ','))
        elif len(ret_node_queries.keys()) == 1 and ret_node_queries.has_key('Subject'):
            if len(ret_node_queries['Subject']) > 3:
                use_title = query_info['title'].replace('$$result$$', ret_node_queries['Subject'][0]['display_name'] + '...' + ret_node_queries['Subject'][-1]['display_name'])
            else:
                use_title = query_info['title'].replace('$$result$$', string.joinfields(map(lambda x: x['display_name'], ret_node_queries['Subject']), ','))
        else:
            use_title = 'ERROR: Unknown title...'
        
        #too many nodes to render efficiently, will allow user to download csv file...
        if display_type == "Download":
            
            return download(request, use_title, query_info['text'], use_query, node_queries, set(temp_node_list), returned_node_type)
            
        else:
            
            #convert the IDs to nodes
            use_nodes = core.get_nodes(temp_node_list, returned_node_type, request,  only_base=True)
            
            ret_nodes = core.apply_grouping2({'nodes':use_nodes, 'links':[]}, [])['nodes']
                
            ret_dict = {'is_graph':True, 'graph':{'nodes':ret_nodes.tolist(), 'links':[]}, 'title':use_title}
        
            return HttpResponse(json.dumps(ret_dict),mimetype="application/json")
    
    except Exception as e:
        print e
        return HttpResponseServerError()

def get_data (request):
    
    request_post = core.request_post_to_json(dict(request.POST.iterlists()))
    
    data_dict, title = config.plot_data_types[request_post["type"]](request_post)
    
    return HttpResponse(json.dumps({'graph':data_dict, 'title':title}),mimetype="application/json")

def copy_nodes(request):
    
    try:
        
        subj_nodes = json.loads(request.POST["subj"])
        query_nodes = json.loads(request.POST["query"])
        
        #print subj_nodes
        #print query_nodes
        
        temp_nodes = []
        
        for i in query_nodes:
            print i
            if i['id'].startswith('@') and i['node_type'] == 'Subject':
                node_ids = custom_functions.match_by_category(i['id'])
                
                temp_nodes.extend(map(lambda x: {'node_type':'Subject', 'id':x}, node_ids))
            
            elif i['id'].startswith('@'):
                raise Exception("Unexpected non-subject ID starting with '@'")
            else:
                temp_nodes.append(i)
        
        cur_graph = core.copy_nodes(subj_nodes, temp_nodes, request, config.edge_queries)
        
        ret_dict = {'nodes':cur_graph['nodes'].tolist(), 'links':cur_graph['links']}
        
        return HttpResponse(json.dumps(ret_dict),mimetype="application/json")
    except Exception as e:
        print e
        return HttpResponseServerError()

def get_graph(request):
    
    try:
    
        request_post = core.fix_jquery_array_keys(dict(request.POST.iterlists()))
        
        nodes, links, title = core.iterate_dict(config.graph_initializers, request_post['panel_context'])(request, request_post)
        
        ret_dict = {'nodes':nodes, 'links':links, 'title':title}
        
        return HttpResponse(json.dumps(ret_dict),mimetype="application/json")

    except:
        return HttpResponseServerError()

@ensure_csrf_cookie
@login_required(login_url='/HitWalker2'+prog_type+'/login/')
def network(request):
    
    
    if len(request.POST) == 0:
        
        return redirect('/HitWalker2'+prog_type)
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
        
        ret_json = {'prog_type':str(prog_type), 'username':str(request.user), 'node_type_transl':json.dumps(node_type_transl), 'edge_type_transl':json.dumps(edge_type_transl),
                    'metanode_thresh':config.max_nodes, 'panel_context':'image', 'input_vals':json.dumps(input_vals), 'cur_param':json.dumps(cur_param), 'cur_filts':json.dumps(cur_filts),
                    'css_path':os.path.basename(new_css_path), 'default_css':json.dumps(default_css_classes)}
        
        for key, val in config.pathway_sizes.items():
            if ret_json.has_key(key) == False:
                ret_json[key] = val;
            else:
                raise Exception("Duplicate entry for " + key)
            
        return render(request, 'network/network.html', ret_json)

@ensure_csrf_cookie
@login_required(login_url='/HitWalker2'+prog_type+'/login/')
def panel(request):
    
    if len(request.POST) == 0:
        
        return redirect('/HitWalker2'+prog_type)
    elif request.POST.has_key("output_format") and request.POST.has_key("data"):
        
        try:
            import cairosvg
            #and all the somewhat invisible dependencies needed to use external css styling on the svgs...
            import lxml
            import tinycss
            import cssselect
        
            if request.POST["plot_type"] == 'svg':
                
                use_svg = '<?xml-stylesheet type="text/css" href="' + request.session['new_css_path'] + '" ?>' + request.POST["data"]
                
            else:
                
                if request.POST['panel_context'] == 'image':
                    width = str(config.pathway_sizes['w'] + config.pathway_sizes['legend_offset'])
                    height = str(config.pathway_sizes['h'])
                    transl = 'translate('+str(config.pathway_sizes['legend_offset'])+',0)'
                elif request.POST['panel_context'] == 'panel':
                    
                    if request.POST["plot_type"] == 'g1':
                        width=str(config.network_sizes['w']+config.network_sizes['legend_offset'])
                        transl = 'translate('+str(config.pathway_sizes['legend_offset'])+',0)'
                    elif request.POST["plot_type"] == 'siRNA_plot':
                        width=str((config.network_sizes['w']*2)+50)
                        transl='translate(50,0)'
                    else:
                        raise Exception('Unknown plot_type')
                    
                    height = str(config.network_sizes['h'])
                else:
                    raise Exception('Unknown value for panel_context')
               
            
                use_data = '<svg xmlns="http://www.w3.org/2000/svg" width="'+width+'" height="'+height+'" >' + request.POST["data"].replace('xmlns="http://www.w3.org/2000/svg"', "")  + '</svg>'
                
                import xml.etree.ElementTree as ET
                tree = ET.fromstring(use_data)
                
                #ensure everything will fit on screen
                #assuming only one non-legend g for now
                tree.find("./*[@class=\'"+request.POST["plot_type"]+"\']").attrib['transform'] = transl
                
                sub_tree = tree.find(".//*[@class='BorderRect BorderSelected']/..")
                
                if sub_tree != None:
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
            raise e
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
            
            request.session.modified = True
            
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
        
        return render(request, 'network/network.html', ret_json)