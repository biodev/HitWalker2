from django.test import LiveServerTestCase

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

from hitwalker_tests import HitWalkerInteraction, SingleMetaNodeSelector

from django.utils import unittest
from django.contrib.auth.models import User
from django.conf import settings
from django.test import TestCase, RequestFactory, Client
from django.utils.importlib import import_module
from django.core.files.storage import get_storage_class

import views
import json
import custom_functions
import config
import core
import copy
import os
import string
import collections
import re
import time
import bs4
import subprocess
import sys
import pyRserve
from py2neo import neo4j, cypher
import numpy as np
import scipy.spatial

test_cypher_session = "http://localhost:7474"

class RTestSession(object):
    
    def __init__(self):
        self._conn = pyRserve.connect()
        self._conn.voidEval('options(stringsAsFactors=F)')
        self._conn.voidEval('library(hwhelper)')
        self._conn.voidEval('load("'+config.hw_config_path+'")')
        self._conn.voidEval('assign("hw2_obj", get(ls()))')
        self._conn.voidEval('source("'+config.test_methods_path+'")')
        
    def getConn(self):
        return self._conn



#@unittest.skip("Skipping selenium")
class BasicSeleniumTests(LiveServerTestCase):
    
    def setUp(self):
        
        #self.driver = webdriver.Firefox()
        self.driver = webdriver.Chrome(executable_path="/Users/bottomly/Desktop/chromedriver")
        # create user
        self.user = User.objects.create_user(username="selenium",
                                             email=None,
                                             password="test")
        #self.client.login(username="selenium", password="test") #Native django test client
        #cookie = self.client.cookies['sessionid']
        #self.driver.get(self.live_server_url + '/HitWalker2')  #selenium will set cookie domain based on current page domain
        #self.driver.add_cookie({'name': 'sessionid', 'value': cookie.value, 'secure': True, 'path': '/'})
        #self.driver.refresh() #need to update page for logged in user
        #self.driver.get(self.live_server_url + '/HitWalker2')
        
        #for some reason the above doesn't work for chrome...
        #self.driver = webdriver.Firefox()
        self.driver.implicitly_wait(1)
        self.driver.get('%s%s' % (self.live_server_url, '/HitWalker2'))
        elem = self.driver.find_element_by_id("id_username")
        elem.send_keys("selenium")
        elem = self.driver.find_element_by_id("id_password")
        elem.send_keys("test")
        #
        self.driver.find_element_by_css_selector("input[type=submit]").click()
    
    def tearDown(self):
        self.driver.quit()
    
    #for the blank text issue: String text = ((JavaScriptExecutor)driver).executeScript("return $(arguments[0]).text();", element);
    #maybe the issue is that the element technically is not displayed, so should move to it first etc...
    def get_metanode_query_tables(self):
        hw_obj = HitWalkerInteraction(self.driver, self.live_server_url)
        
        hw_obj.panel_by_query("@liver")
        
        panel_1 = hw_obj.to_panel("1")
        
        hw_obj.select_context_node(panel_1, SingleMetaNodeSelector())
        
        datatype_divs = self.driver.find_elements_by_css_selector("#pg1 > div.panel.panel-default")
        
        table_dict = collections.OrderedDict()
        
        for i in datatype_divs:
            cur_link = i.find_element_by_css_selector("div > h4 > a")
            
            self.driver.execute_script("arguments[0].scrollIntoView(true);", cur_link);
            
            link_text = cur_link.text
            cur_link.click()
            #wait until table appears
            element = WebDriverWait(self.driver, 20).until(
                EC.presence_of_element_located((By.CSS_SELECTOR,"#"+i.get_attribute("id") + " > div > div > table"))
            )
            
            self.driver.execute_script("arguments[0].scrollIntoView(true);", element);
            
            #Assuming table is of the form: Genes/Subjects, Frequency, View/Download statements
            
            cur_table = []
            tab_row = element.find_elements_by_css_selector("tbody > tr")
            
            #need to do each chunk of selecting at a single time otherwise get a stale reference exception
            for i in tab_row:
                temp_row = tuple(map(lambda x: x.text,i.find_elements_by_tag_name("td")))
                cur_table.append(temp_row[:2])
            
            table_dict[link_text] = cur_table
            
        return table_dict
    
    #def test_metanode_query_click(self):
    #    
    #    found_tables = self.get_metanode_query_tables()
    #    
    #    print found_tables
    #    
    #    hw_obj = HitWalkerInteraction(self.driver, self.live_server_url)
    #    
    #    panel_1 = hw_obj.to_panel("1")
    #    
    #    empty_portion_panel = webdriver.ActionChains(self.driver).move_to_element_with_offset(panel_1, 0, 0).click().perform()
    #    
    #    cur_panel = 1
    #    
    #    for i_ind, i in enumerate(found_tables.items()):
    #        for j_ind, j in enumerate(i[1]):
    #            
    #            print i
    #            
    #            panel_1 = hw_obj.to_panel("1")
    #            hw_obj.select_context_node(panel_1, SingleMetaNodeSelector())
    #            #need to click the link here...
    #            
    #            all_links = self.driver.find_elements_by_css_selector("#pg1 a")
    #            
    #            all_links[i_ind].click()
    #           
    #            element = WebDriverWait(self.driver, 20).until(
    #                EC.presence_of_element_located((By.CSS_SELECTOR,"#pg1 > div.panel.panel-default table"))
    #            )
    #            
    #            #this just makes sure that everything is able to be clicked (assuming they are ready when one of them is...)
    #            span_el = WebDriverWait(element, 20).until(
    #                EC.element_to_be_clickable((By.CSS_SELECTOR,"span"))
    #            )
    #            
    #            all_spans = element.find_elements_by_css_selector("span")
    #            
    #            all_spans[j_ind].click()
    #            
    #            cur_panel += 1
    #            
    #            #TODO: Need to deal with the 'Download' possibility
    #            
    #            result_panel = WebDriverWait(self.driver, 20).until(
    #                EC.element_to_be_clickable((By.CSS_SELECTOR,"g.g1[id='panel_"+str(cur_panel)+"']"))
    #            )
    #            
    #            #check the retrieved metanode against the results from the table
    #            
    #            if j[0] != '' and int(j[0]) > 1:
    #            
    #                attrs = hw_obj.get_metanode_attrs(result_panel, SingleMetaNodeSelector())
    #                
    #                print attrs
    #                print j, j_ind
    #                
    #                self.assertTrue(attrs['type'] == 'Gene')
    #                
    #                if j[0] != '':
    #                
    #                    self.assertTrue(attrs['count'] == int(j[0]))
    #                    self.assertTrue(hw_obj.count_metanode_children(cur_panel, SingleMetaNodeSelector()) == attrs['count'])
    #            elif j[0] != '':
    #                self.assertTrue(hw_obj.count_nodes(cur_panel, "Gene") == 1)
    #                
    #            #delete the new panel
    #            
    #            hw_obj.delete_panel(cur_panel)
    #    
    
    #def test_metanode_query_table(self):
    #    
    #    table_dict = self.get_metanode_query_tables()
    #    
    #    if r_obj != None:
    #        print 'r object exists! testing...'
    #        
    #        print table_dict
    #        
    #        for i in table_dict.items():
    #            
    #            print i
    #            
    #            if i[0] != '':
    #                
    #                dta = r_obj.getConn().r.getFrequency(r_obj.getConn().ref.hw2_obj, i[0], 'liver', 'Subject_Category')
    #                
    #                print dta
    #                
    #                if isinstance(dta, pyRserve.TaggedList):
    #                    
    #                    for j_ind, j in enumerate(i[1]):
    #                        if j[0] != '':
    #                            self.assertEqual(int(j[0]), dta['Genes'][j_ind])
    #                            self.assertEqual(str(j[1]), str(dta['Frequency'][j_ind]))
    #                    
    #                else:
    #                    print 'skipping test for '+i[0]+' due to invalid returned object'
    #    else:
    #        print 'r object does not exist... skipping tests'
            
            
    #def test_metanode_subsetting(self):
    #    
    #    ##NOTE: Make sure to replace liver with another subject attr
    #    
    #    hw_obj = HitWalkerInteraction(self.driver, self.live_server_url)
    #    
    #    hw_obj.panel_by_query("@liver")
    #    
    #    cur_meta = hw_obj.get_metanode("1", SingleMetaNodeSelector())
    #    
    #    webdriver.ActionChains(self.driver).move_to_element(cur_meta).context_click(cur_meta).perform()
    #    
    #    #first the subset spans, check the numbers versus the data in R
    #    
    #    trs = self.driver.find_elements_by_css_selector("#summary_table > tbody > tr")
    #    
    #    cur_th = ""
    #    
    #    ret_table = []
    #    
    #    for i in trs:
    #        try:
    #            cur_th = i.find_element_by_css_selector("th").text
    #        except Exception as e:
    #            pass
    #        
    #        cur_row = [cur_th]
    #        
    #        cur_row.extend(map(lambda x: x.text, i.find_elements_by_css_selector("td")))
    #        
    #        ret_table.append(cur_row)
    #    
    #    print ret_table
    #    
    #    if r_obj != None:
    #        print 'r object exists! testing...'
    #        
    #        dta = r_obj.getConn().r.subjectAttrs(r_obj.getConn().ref.hw2_obj, 'liver', 'Subject_Category')
    #        
    #        print dta
    #        
    #        for i in ret_table:
    #            if i[1].endswith("..."):
    #                
    #                find_str = i[1].replace(" ...", "")
    #                val_ind = np.where(np.char.find(np.char.capitalize(dta['Value']), find_str) > -1)
    #                
    #            else:
    #                val_ind = np.where(np.char.capitalize(dta['Value']) == i[1])
    #            
    #            type_ind = np.where(dta['Type'] == i[0])
    #            
    #            common_ind = set(val_ind[0]).intersection(set(type_ind[0]))
    #            
    #            if len(common_ind) > 0:
    #                
    #                self.assertTrue(int(dta['Count'][common_ind.pop()]) == int(i[2]))
    #            
    #            else:
    #                print 'common index not found: ' + str(val_ind) + str(type_ind)
    #                self.assertTrue(False)
    #            
    #        
    #    else:
    #        print 'r object does not exist, skipping r tests'
    #    
    #    #then perform several subsets and check them relative to the span text
    #    
    #    subset_spans = self.driver.find_elements_by_css_selector("#summary_table > tbody > tr > td > span")
    #    
    #    #just the first for now
    #    
    #    first_span_size = int(subset_spans[0].text)
    #    
    #    subset_spans[0].click()
    #    
    #    #check the metanode size by the value in the span element
    #    
    #    p2_meta_count = hw_obj.count_metanode_children("2", SingleMetaNodeSelector())
    #    
    #    self.assertTrue(first_span_size == p2_meta_count)
    #    
    #    #then limit the second metanode by another feature
    #    
    #    p2_meta = hw_obj.get_metanode("2", SingleMetaNodeSelector())
    #    
    #    webdriver.ActionChains(self.driver).move_to_element(p2_meta).context_click(p2_meta).perform()
    #    
    #    second_spans = self.driver.find_elements_by_css_selector("#summary_table > tbody > tr > td > span")
    #    
    #    second_span_size = int(second_spans[-1].text)
    #    
    #    second_spans[-1].click()
    #    
    #    p3_meta_count = hw_obj.count_metanode_children("3", SingleMetaNodeSelector())
    #    
    #    self.assertTrue(second_span_size == p3_meta_count)
    #    
    #    #now do the single/multiple node select box which should be the actual names back on the first panel
    #    
    #    cur_meta = hw_obj.get_metanode("1", SingleMetaNodeSelector())
    #    
    #    webdriver.ActionChains(self.driver).move_to_element(cur_meta).context_click(cur_meta).perform()
    #    
    #    self.driver.find_element_by_css_selector(".select2-input").click()
    #    all_lis = self.driver.find_elements_by_css_selector("#select2-drop > ul > li")
    #    
    #    select_text = map(lambda x: str(x.text), all_lis)
    #    
    #    subj_names = r_obj.getConn().r.subjectSubset(r_obj.getConn().ref.hw2_obj, 'liver', 'Subject_Category')
    #    
    #    self.assertListEqual(sorted(select_text), sorted(subj_names))
    #    
    #    all_lis[0].click()
    #    
    #    self.driver.find_element_by_css_selector(".select2-input").click()
    #    all_lis = self.driver.find_elements_by_css_selector("#select2-drop > ul > li")
    #    
    #    all_lis[-1].click()
    #   
    #    self.driver.find_element_by_css_selector("#subset_samp_button").click()
    #    
    #    #there should only be a single metanode on the panel with 2 children
    #    
    #    p4_meta_count = hw_obj.count_metanode_children("4", SingleMetaNodeSelector())
    #    
    #    self.assertTrue(p4_meta_count == 2)
        
    #def test_gene_addition_metanode(self):
    #    
    #    hw_obj = HitWalkerInteraction(self.driver, self.live_server_url)
    #    
    #    hw_obj.panel_by_query("@liver")
    #    
    #    panel_1 = hw_obj.to_panel("1")
    #    
    #    hw_obj.click_context_button("1", 2)
    #    
    #    self.driver.find_element_by_css_selector("ul li:first-child a").click()
    #    
    #    #select a gene that probably has hits, say TP53
    #    
    #    self.driver.find_element_by_css_selector(".select2-choice").click()
    #    self.driver.find_element_by_css_selector("input.select2-input").send_keys("TP53")
    #    
    #    input_highlight = WebDriverWait(self.driver, 20).until(
    #        EC.presence_of_element_located((By.CSS_SELECTOR,".select2-result-label"))
    #        )
    #    
    #    input_highlight.click()
    #    
    #    self.driver.find_element_by_xpath("//button[.='OK']").click()
    #    
    #    #figure out which samples have hits for TP53 via R
    #    
    #    if r_obj != None:
    #    
    #        r_obj.getConn().r.gene_hits = r_obj.getConn().r.findHits(r_obj.getConn().ref.hw2_obj, 'liver', 'TP53', 'Subject_Category')
    #    
    #        subj_groups = r_obj.getConn().r.encode_groups(r_obj.getConn().ref.gene_hits)
    #        
    #        print subj_groups
    #        
    #        hit_cat_set = collections.Counter(subj_groups['FixedDt'])
    #        
    #        print hit_cat_set
    #        
    #        metanode_list = hw_obj.get_panel_metanodes("2")
    #        
    #        trans_vals = map(lambda x: map(float, re.split("[,\)\(]", x.get_attribute("transform"))[1:3]), metanode_list)
    #        
    #        print trans_vals
    #        
    #        kd_tree = scipy.spatial.KDTree(np.array(trans_vals))
    #        
    #        metanode_count = []
    #        
    #        for i in metanode_list:
    #            metanode_count.append(hw_obj.count_metanode_children(None, i))
    #        
    #        print metanode_count
    #        
    #        links = hw_obj.to_panel("2")
    #        
    #        link_els = links.find_elements_by_css_selector("path.link")
    #        
    #        link_pos = map(lambda x: map(float, re.split("[M,L]", x.get_attribute("d"))[1:]) , link_els)
    #        
    #        #in this case all the links should point to the same position, so we are just looking at which metanode should be assigned the relationship
    #        
    #        link_type = map(lambda x: x.get_attribute("class").replace("link ", ""), link_els)
    #        
    #        link_dict = collections.defaultdict(list)
    #        
    #        for i_ind, i in enumerate(link_pos):
    #            link_match = kd_tree.query(np.array(i[:2]))
    #            link_dict[str(link_match[1])].append(link_type[i_ind])
    #            
    #        res_dict = {}
    #        
    #        for i in link_dict.items():
    #            res_dict[string.joinfields(i[1], ",")] = metanode_count[int(i[0])]
    #        
    #        self.assertDictEqual(hit_cat_set, res_dict)
    
    def compare_dicts(self, dict1, dict2):
        for i in dict1.items():
                for j in i[1].items():
                    if dict2.has_key(i[0]):
                        if dict2[i[0]].has_key(j[0]):
                            if dict2[i[0]][j[0]] != j[1]:
                                raise Exception('sets not the same:' + str(dict2[i[0]][j[0]]) + ' ' + str(j[1]))
                        else:
                            raise Exception('dict2 missing key: ' + i[0] + '->' + j[0])
                    else:
                        raise Exception('dict2 missing key: ' + i[0])
    
    #def test_hitwalker_panel(self):
    #    
    #    hw_obj = HitWalkerInteraction(self.driver, self.live_server_url)
    #    
    #    default_thresh_mat_base = config.prioritization_func['args']['initial_graph_file'].replace(".mtx", "")
    #    
    #    if r_obj != None:
    #        
    #        gene_text = hw_obj.panel_by_prioritize("HEPG2_LIVER")
    #        
    #        found_rels = hw_obj.get_node_rels("1")
    #        
    #        #needed to do it this way as it kept interpreting sub_graph as a list using the .ref approach
    #        r_obj.getConn().voidEval("sub_graph <- process_matrix_graph('"+default_thresh_mat_base+"', .4)")
    #        
    #        gene_links = r_obj.getConn().r.get_gene_connections(r_obj.getConn().ref.sub_graph, gene_text['seeds'], gene_text['targs'])
    #        
    #        r_found_dta = collections.defaultdict(lambda: collections.defaultdict(set))
    #        
    #        for i in range(0, len(gene_links['from'])):
    #            r_found_dta[gene_links['from'][i]][gene_links['to'][i]].add('STRING')
    #            r_found_dta[gene_links['to'][i]][gene_links['from'][i]].add('STRING')
    #        
    #        #make this into a unique set
    #        gene_set = set(list(gene_links['from']) + list(gene_links['to']))
    #        
    #        r_obj.getConn().r.gene_hits = r_obj.getConn().r.findHits(r_obj.getConn().ref.hw2_obj, 'HEPG2_LIVER', np.array(list(gene_set)), 'Subject', 'Gene')
    #        
    #        subj_groups = r_obj.getConn().r.encode_groups(r_obj.getConn().ref.gene_hits, True, 'HEPG2_LIVER')
    #        
    #        for i in range(0, len(subj_groups['Subject'])):
    #            split_dt = subj_groups['FixedDt'][i].split(',')
    #            for j in split_dt:
    #                r_found_dta[subj_groups['Subject'][i]][subj_groups['Gene'][i]].add(j)
    #                r_found_dta[subj_groups['Gene'][i]][subj_groups['Subject'][i]].add(j)
    #        
    #        #for i in subj_groups['Subject']:
    #        #    for j in subj_groups['Gene']:
    #        #
    #        print 'R vs Screen'
    #        self.compare_dicts(r_found_dta, found_rels)
    #        print 'Screen vs R'
    #        self.compare_dicts(found_rels, r_found_dta)
    #    
    #    self.assertEqual(r_found_dta, found_rels)
    
    #def test_gene_addition_prioritize(self):
    #    
    #    hw_obj = HitWalkerInteraction(self.driver, self.live_server_url)
    #
    #    gene_text = hw_obj.panel_by_prioritize("HEPG2_LIVER")
    #    
    #    #add a new gene
    #    hw_obj.add_gene("1", "TP53")
    #    
    #    node_rels = hw_obj.get_node_rels("2")
    #    
    #    #assume that all went well with the formation of the graph and just add in the nodes from the page...
    #    gene_list = hw_obj.get_node_rels("1").keys() + ['TP53']
    #    
    #    if r_obj != None:
    #        
    #        print gene_list
    #        
    #        r_obj.getConn().r.gene_hits = r_obj.getConn().r.findHits(r_obj.getConn().ref.hw2_obj, "HEPG2_LIVER", gene_list, 'Subject', 'Gene')
    #        
    #        r_obj.getConn().voidEval('print(deparse(gene_hits))')
    #        
    #        subj_groups = r_obj.getConn().r.encode_groups(r_obj.getConn().ref.gene_hits, True, "HEPG2_LIVER", "Gene")
    #        
    #        #as there is a single subject, then group the genes into metanodes
    #        
    #        #gene_counter = collections.Counter(subj_groups['FixedDt'])
    #        
    #        r_found_dta = collections.defaultdict(lambda: collections.defaultdict(set))
    #        
    #        for i in range(0, len(subj_groups['Subject'])):
    #            split_dt = subj_groups['FixedDt'][i].split(',')
    #            for j in split_dt:
    #                r_found_dta[subj_groups['Subject'][i]][subj_groups['Gene'][i]].add(j)
    #                r_found_dta[subj_groups['Gene'][i]][subj_groups['Subject'][i]].add(j)
    #        
    #        print 'R vs Screen'
    #        self.compare_dicts(r_found_dta, node_rels)
    #        print 'Screen vs R'
    #        self.compare_dicts(node_rels, r_found_dta)
    #        
    #        self.assertEqual(r_found_dta, node_rels)
    
    #def test_pathway_addition_prioritize(self):
    #
    #    hw_obj = HitWalkerInteraction(self.driver, self.live_server_url)
    #
    #    gene_text = hw_obj.panel_by_prioritize("HEPG2_LIVER")
    #    
    #    #add a pathway this time
    #    pathway_name = hw_obj.add_pathway("1", ["TP53"])
    #    
    #    node_rels = hw_obj.get_node_rels("2")
    #    
    #    gene_list = hw_obj.get_node_rels("1")
    #    
    #    clean_pathway_name = re.sub("\s+\(n=\d+\)\s*", "", pathway_name)
    #    
    #    if r_obj != None:
    #        
    #        r_obj.getConn().r.gene_hits1 = r_obj.getConn().r.findHits(r_obj.getConn().ref.hw2_obj, "HEPG2_LIVER", clean_pathway_name, 'Subject', 'Pathway')
    #        
    #        #also need to collect the genes that were on panel 1
    #        
    #        #r_obj.getConn().voidEval('gene_hits2 <- findHits(hw2_obj, "HEPG2_LIVER", "'+clean_pathway_name+'", "Subject", "Pathway")')
    #        
    #        r_obj.getConn().r.gene_hits2 = r_obj.getConn().r.findHits(r_obj.getConn().ref.hw2_obj, "HEPG2_LIVER", gene_list.keys(), 'Subject', 'Gene')
    #        
    #        r_obj.getConn().voidEval('gene_hits_all <- rbind(as.data.frame(gene_hits1), as.data.frame(gene_hits2))')
    #        
    #        subj_groups = r_obj.getConn().r.encode_groups(r_obj.getConn().ref.gene_hits_all, True, "HEPG2_LIVER", "Gene")
    #        
    #        #as there is a single subject, then group the genes into metanodes
    #        
    #        r_found_dta = collections.defaultdict(lambda: collections.defaultdict(set))
    #        
    #        for i in range(0, len(subj_groups['Subject'])):
    #            split_dt = subj_groups['FixedDt'][i].split(',')
    #            for j in split_dt:
    #                r_found_dta[subj_groups['Subject'][i]][subj_groups['Gene'][i]].add(j)
    #                r_found_dta[subj_groups['Gene'][i]][subj_groups['Subject'][i]].add(j)
    #        
    #        print 'R vs Screen'
    #        self.compare_dicts(r_found_dta, node_rels)
    #        print 'Screen vs R'
    #        self.compare_dicts(node_rels, r_found_dta)
    #        
    #        self.assertEqual(r_found_dta, node_rels)
    #
    
    #def test_subject_addition_prioritize(self):
    #    
    #    subjects = ["HEPG2_LIVER", "LOUCY_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE"]
    #    
    #    hw_obj = HitWalkerInteraction(self.driver, self.live_server_url)
    #
    #    gene_text = hw_obj.panel_by_prioritize(subjects[0])
    #    
    #    #add a new subject this time...
    #    
    #    hw_obj.add_subject("1", subjects[1])
    #    
    #    node_rels = hw_obj.get_node_rels("2")
    #    
    #    gene_list = hw_obj.get_node_rels("1")
    #    
    #    if r_obj != None:
    #        #it shouldn't matter that gene_list will contain subject nodes as well--they will be ignored
    #        r_obj.getConn().r.gene_hits = r_obj.getConn().r.findHits(r_obj.getConn().ref.hw2_obj, subjects, gene_list.keys(), 'Subject', 'Gene')
    #        
    #        subj_groups = r_obj.getConn().r.encode_groups(r_obj.getConn().ref.gene_hits, True, "HEPG2_LIVER", "Gene")
    #        
    #        r_found_dta = collections.defaultdict(lambda: collections.defaultdict(set))
    #        
    #        for i in range(0, len(subj_groups['Subject'])):
    #            split_dt = subj_groups['FixedDt'][i].split(',')
    #            for j in split_dt:
    #                r_found_dta[subj_groups['Subject'][i]][subj_groups['Gene'][i]].add(j)
    #                r_found_dta[subj_groups['Gene'][i]][subj_groups['Subject'][i]].add(j)
    #        
    #        print 'R vs Screen'
    #        self.compare_dicts(r_found_dta, node_rels)
    #        print 'Screen vs R'
    #        self.compare_dicts(node_rels, r_found_dta)
    #        
    #        self.assertEqual(r_found_dta, node_rels)
    
    def test_context_pathway_mode_prioritize(self):
        hw_obj = HitWalkerInteraction(self.driver, self.live_server_url)
        
        gene_text = hw_obj.panel_by_prioritize('HEPG2_LIVER')
        
        #go to the context pathway mode screen
        
        hw_obj.click_context_button("1", 3)
        
        WebDriverWait(self.driver, 20).until(
            EC.element_to_be_clickable((By.CSS_SELECTOR,".select2-search-choice-close"))
        )
        
        close_buttons = self.driver.find_elements_by_css_selector(".select2-search-choice-close")
        
        #check that five genes are selected
        
        self.assertTrue(len(close_buttons) == 5)
        
        #unselect them all and go through the subject_pathway interface
        
        for i in close_buttons:
            i.click()
        
        time.sleep(2)
        
        path_text = hw_obj.add_pathway(None, ["KRAS"])
        
        #wait until network screen is present
        
        spinner = WebDriverWait(self.driver, 20).until(
            EC.presence_of_element_located((By.CSS_SELECTOR, "div.spinner"))
        )
        
        WebDriverWait(self.driver, 120).until(
            EC.staleness_of(spinner)
        )
        
        node_rels = hw_obj.get_node_rels("1")
        
        print node_rels
        
        #if r_obj != None:
        
        #get string
        
        #get hits
        

##globally useful functions and classes

class BasicNodeWithAttr(core.BasicNode):
    def __init__(self,res_list, only_child=False):
        self.node_dict = {'id':res_list[0], 'display_name':res_list[0], 'attributes':{'node_type':res_list[2]}, 'children':core.NodeList()}
        
        self.id = self.node_dict['id']
        self.display_name = self.node_dict['display_name']
        
        if only_child == True:
            if len(res_list) < 2:
                raise Exception("only_child == True requires res_list to be at least of length 2")
            self.node_dict['children'].add(core.BasicChild(res_list))

def make_gene_list(res_list):
    nodes = core.NodeList()
    
    seed_header = ['gene', 'sample', 'var', 'score', 'is_hit']
    
    for ind, gene_result in enumerate(res_list):
        
        if nodes.hasNode(gene_result[0]) == False:
            #create the gene prior to adding if it doesn't exist 
            nodes.add(core.GeneNode([gene_result[0], gene_result[0], []]))
        
        seed_node = core.SeedNode(nodes.getNode(gene_result[0]), gene_result, seed_header)
        
        #then add the hits to the GeneNode
        nodes.addChild(gene_result[0], seed_node)
            
    return nodes
    
def simple_handler(res_list, nodes, request):
    for i in core.BasicResultsIterable(res_list):
        if len(i) > 0:
            nodes.add(core.BasicNode(i))


def simple_handler_w_child(res_list, nodes, request):
    for i in core.BasicResultsIterable(res_list):
        if len(i) > 0:
            nodes.add(BasicNodeWithAttr(i, only_child=True))

test_params = {
    'General_Parameters':{'type':'standard',
                          'fields':{
                                'string_conf':{'type':'numeric', 'default':.4, 'range':[0, 1], 'comparison':'>', 'name':'HitWalker STRING Confidence'},
                                'path_conf':{'type':'numeric', 'default':.95, 'range':[0, 1], 'comparison':'>', 'name':'Pathway STRING Confidence'},
                                'res_prob':{'type':'numeric', 'default':.3, 'range':[0,1], 'comparison':'=', 'name':'Restart Probability'},
                                'max_iter':{'type':'numeric', 'default':100, 'range':[0, 10000], 'comparison':'=', 'name':'Max Iterations'},
                                'conv_thresh':{'type':'numeric', 'default':1e-10, 'range':[0,1], 'comparison':'<', 'name':'Convergence Threshold'},
    'expression':{'type':'numeric', 'default':0.75 , 'range':[0,1], 'comparison':'>' , 'name':'Expression (Hit) Threshold'},
    'genescore':{'type':'numeric', 'default':0 , 'range':[-100,100], 'comparison':'>' , 'name':'GeneScore (Hit) Threshold'}
                                }, 
                        }
}

try:

    r_obj = RTestSession()
except:
    r_obj = None

class Test_views(TestCase):
    
    maxDiff = None
    
    def setUp(self):
        # Every test needs access to the request factory.
        self.factory = RequestFactory()
        self.user = User.objects.create_user(
            username='jacob', email='', password='top_secret')
        
        settings.SESSION_ENGINE = 'django.contrib.sessions.backends.file'
        engine = import_module(settings.SESSION_ENGINE)
        store = engine.SessionStore()
        store.save()
        
        self.store = store
    
    def get_result_list(self, content):
        html_res = bs4.BeautifulSoup(content)
        
        tab_res = html_res.find_all('tr')
        
        mut_res = []
        
        #outside of the top ten, not going to be able to make meaningful comparisions as the values tend to get pretty small and
        #into rounding errors...
        for i in tab_res[:10]:
            temp_list = list(i.stripped_strings)
            if len(temp_list) > 0:
                mut_res.append(temp_list[0])
                
        return mut_res
    
    def test_table(self):
        
        test_post_1 = {'csrfmiddlewaretoken': 'test',
            'query_samples': json.dumps({"SampleID":{"Expression":"METIS_p_NCLE_RNA1_Human_U133_Plus_2.0_A09_240866.CEL",
                                           "GeneScore":"HEPG2_LIVER",
                                           "Variants":"HEPG2_LIVER"}}),
            'parameters':json.dumps(test_params)}
        
        request = self.factory.post('/table/', test_post_1)
        
        # middleware are not supported. You can simulate a
        # logged-in user by setting request.user manually.
        request.user = self.user
        request.session = self.store
        
        #remove any pre-created matrix file
        
        default_thresh_mat = config.prioritization_func['args']['initial_graph_file'].replace(".mtx", ".400.mtx")
        
        if os.path.exists(default_thresh_mat):
            os.remove(default_thresh_mat)
            os.remove(default_thresh_mat.replace(".mtx", ".names"))
        
        response = views.table(request)
        self.assertEqual(response.status_code, 200)
        
        mut_res1 = self.get_result_list(response.content)
        
        #again, this time keeping the subsetted matrix
        
        response2 = views.table(request)
        
        mut_res2 = self.get_result_list(response2.content)
        
        self.assertListEqual(mut_res1, mut_res2)
        
        #compare this with the R equivalent
        
        def r_rwr(request, seed_gene_nl, initial_graph_file, string_session_name, res_prob_session_name, conv_thresh_session_name, max_iter_session_name):
            
            seed_dict = seed_gene_nl.todict()
            
            subprocess.call(["Rscript", "--vanilla",  os.path.join("network", "perform_rwr.R"),
                             initial_graph_file.replace(".mtx", ""),
                             str(request.session[string_session_name])] + seed_dict.keys())
            
            result_inp = open("temp_ranking_results.txt", "r")
            
            prot_nl = core.NodeList()
            
            for i in result_inp:
                split_i = i.strip('\n').split()
                prot_score = [split_i[0], float(split_i[1])]
                temp_node = core.BasicNode(prot_score, only_child=True)
                prot_nl.add(temp_node)
            
            result_inp.close()
            
            os.remove("temp_ranking_results.txt")
            
            return core.SeedList(prot_nl)
            
        config.prioritization_func['function'] = r_rwr
        
        responseR = views.table(request)
        self.assertEqual(responseR.status_code, 200)
        
        mut_res2 = self.get_result_list(responseR.content)
        
        #checks that the variants are arranged in the same order...
        self.assertListEqual(mut_res1, mut_res2)
        
        #could also check the sanity of the list to choose from for the visualizations
        

###tests for the core module

class Test_node_classes(TestCase):
    
    def setUp(self):
        
        #make a couple NodeLists
        self.sirna_nodes = make_gene_list([
                        (u'ENSG00000254087', u'12-00145', u'siRNA',-2, True),
                        (u'ENSG00000101336', u'12-00145',u'siRNA', -1, False),
                        (u'ENSG00000213341', u'12-00145', u'siRNA',-2, True)
                        ])
        
        self.gene_score = make_gene_list([
                        (u'ENSG00000254087', u'12-00145', u'GeneScore',10, True),
                        (u'ENSG00006898999', u'12-00145',u'GeneScore', 15, False),
                        (u'ENSG00000213765', u'12-00145', u'GeneScore',30, True)
                        ])
        
        #make a SeedList as well
        temp_nl = copy.deepcopy(self.sirna_nodes)
        temp_nl.mergeChildren(self.gene_score)
        
        self.test_sl = core.SeedList(temp_nl)
    
    def test_addChild(self):
       test_1 = copy.deepcopy(self.sirna_nodes)
       
       test_1_node = self.gene_score.getNode("ENSG00000254087").children().next()
       
       test_1.addChild("ENSG00000254087", test_1_node)
       
       self.assertEqual(sorted(test_1.getNode('ENSG00000254087').children().ids()), sorted(['ENSG00000254087_siRNA', 'ENSG00000254087_GeneScore']))
    
    def test_add (self):
        test_1 = copy.deepcopy(self.sirna_nodes)
       
        test_1_node = self.gene_score.getNode("ENSG00000254087")
    
        test_1.add(test_1_node)
        
        self.assertEqual(sorted(test_1.ids()), sorted(['ENSG00000254087', 'ENSG00000101336', 'ENSG00000213341','ENSG00000254087']))
        
    def test_extend(self):
        #just extend with no additional checking
        test_1 = copy.deepcopy(self.sirna_nodes)
        
        test_1.extend(self.gene_score)
        
        self.assertEqual(sorted(test_1.ids()), sorted(['ENSG00000254087', 'ENSG00000101336', 'ENSG00000213341','ENSG00000254087', 'ENSG00006898999','ENSG00000213765']))
    
    def test_mergeChildren(self):
        
        #the mergeChildren method should behave similarly to extendIfNew, but should also add in new children
        
        test_2 = copy.deepcopy(self.sirna_nodes)
        
        test_2.mergeChildren(self.gene_score)
        
        self.assertEqual(sorted(test_2.ids()), sorted(['ENSG00000254087', 'ENSG00000101336', 'ENSG00000213341','ENSG00006898999','ENSG00000213765']))
        
        #Now the overlapping gene should have both an siRNA and gene_score gene
        
        self.assertEqual(sorted(test_2.getNode('ENSG00000254087').children().ids()), sorted(['ENSG00000254087_siRNA', 'ENSG00000254087_GeneScore']))
        
    
    def test_extendIfNew(self):
        
        test_1 = copy.deepcopy(self.sirna_nodes)
        
        test_1.extendIfNew(self.gene_score)
        
        self.assertEqual(sorted(test_1.ids()), sorted(['ENSG00000254087', 'ENSG00000101336', 'ENSG00000213341','ENSG00006898999','ENSG00000213765']))
        
        #The overlapping gene 'ENSG00000254087' should only have an siRNA child
        
        self.assertEqual(test_1.getNode('ENSG00000254087').children().ids(), ['ENSG00000254087_siRNA'])
        
    def test_summarizeChildren(self):
        
        test_4 = copy.deepcopy(self.gene_score)
        
        self.assertEqual(test_4.summarizeChildren(lambda x: x.getAttr(["attributes","meta", "score"]), max), [10,15,30])
        
        #check summarizeChildren's behavior when there are no children
        test_5 = copy.deepcopy(self.gene_score)
        
        for i in range(0, len(self.gene_score)):
            test_5.node_list[i].node_dict['children'] = core.NodeList()
            
        self.assertEqual(test_5.summarizeChildren(lambda x: x.getAttr(["attributes","meta", "score"]), max), [None, None, None])
    
    def test_filterByChild(self):
        
        test_3 = copy.deepcopy(self.sirna_nodes)
        
        test_3.filterByChild(lambda x: x.getAttr(["attributes", "meta", "is_hit"]), all)
        
        self.assertEqual(sorted(test_3.ids()), sorted(['ENSG00000254087', 'ENSG00000213341']))
        
    def test_iteration(self):
        
        test_1 = copy.deepcopy(self.sirna_nodes)
        
        self.assertTrue(len(test_1) == 3)
        
        count = 0
        
        for i in test_1:
            self.assertTrue(isinstance(i, core.Node))
            count += 1
        
        self.assertTrue(count == 3)
        
        #check that it resets itself once done
        
        for i in test_1:
            count += 1
            
        self.assertTrue(count == 6)
    
    ##then onto SeedList
    
    def test_subset(self):
        
        temp_sl = copy.deepcopy(self.test_sl)
        
        temp_sl.subset(['ENSG00000101336'])
        
        self.assertTrue(len(temp_sl.nodeList()) == 1)
        
        self.assertJSONEqual(temp_sl.nodeList().json(), json.dumps([self.sirna_nodes.getNode('ENSG00000101336').todict()]))
        
        self.assertEqual(temp_sl.node_scores, self.test_sl.getScores(['ENSG00000101336']))
    
    def test_adjustScores(self):
        
        temp_sl = copy.deepcopy(self.test_sl)
        
        self.assertEqual(sorted(temp_sl.node_scores), sorted([10,15,30,-1,-2]))
        
        temp_sl.adjustScores({'ENSG00000254087':30,'ENSG00000101336':40})
        
        self.assertEqual(sorted(temp_sl.node_scores), sorted([30,15,30,40,-2]))
        
    def test_getScores(self):
        
        temp_sl = copy.deepcopy(self.test_sl)
        
        self.assertEqual(temp_sl.getScores(['ENSG00006898999', 'ENSG00000213765']), [15,30])
    
    def test_todict(self):
        
        temp_sl = copy.deepcopy(self.test_sl)
        
        sl_dict = temp_sl.todict()
        
        self.assertDictEqual(sl_dict, {'ENSG00000254087':10, 'ENSG00000101336':-1, 'ENSG00000213341':-2, 'ENSG00006898999':15, 'ENSG00000213765':30})
       
    
class Test_core_classes (TestCase):
    
    def test_BasicResultsIterable(self):
        
        #get a single result from a single transaction
        
        session = cypher.Session()
        tx = session.create_transaction()
        
        tx.append("MATCH (n) RETURN n.name limit 1")
        
        test_results = tx.commit()
        
        test_list_1 = list(core.BasicResultsIterable(test_results))
        
        #this is a simplification for the case of a single entry so the user doesn't always have to
        #iterate over an unnecessary tuple
        self.assertTrue(len(test_list_1[0]) == 1)
        self.assertTrue(isinstance(test_list_1[0][0], unicode))
        
        #get a bunch from several transactions
        
        tx = session.create_transaction()
        
        tx.append("MATCH (n) RETURN n.name limit 10")
        tx.append("MATCH (n) RETURN n.name skip 10 limit 10")
        
        test_results_2 = tx.commit()
        
        for i in core.BasicResultsIterable(test_results_2):
            self.assertTrue(isinstance(i, list))
            self.assertTrue(len(i) == 10)
            for j in i:
                self.assertTrue(isinstance(j, tuple))
                self.assertTrue(len(j) == 1)
                self.assertTrue(isinstance(j[0], unicode))
        
        #get none
        
        tx = session.create_transaction()
        
        tx.append("MATCH (n) RETURN n.name limit 0")
        
        test_results_3 = tx.commit()
        
        test_list_3 = list(core.BasicResultsIterable(test_results_3))
        
        self.assertTrue(len(test_list_3[0]) == 0)
        
    def test_proper_type (self):
        self.assertEquals(core.proper_type("1"), 1)
        self.assertEquals(core.proper_type("1.5"), 1.5)
        self.assertEquals(core.proper_type("None"), None)
        self.assertEquals(core.proper_type("[1,2,10]"), [1,2,10])
        self.assertEquals(core.proper_type('{"test":[1,2,10]}'), {"test":[1,2,10]})
    
    def test_fix_jquery_array_keys(self):
        
        test_1_res = core.fix_jquery_array_keys({'var_select[]': [u'ENSG00000010327', u'ENSG00000196132', u'ENSG00000181929'], 'panel_context': [u'panel'], 'seed_select[]': [u'ENSG00000109320', u'ENSG00000092445', u'ENSG00000182578']})
        
        self.assertEqual(test_1_res, {'var_select':[u'ENSG00000010327', u'ENSG00000196132', u'ENSG00000181929'], 'panel_context':[u'panel'], 'seed_select':[u'ENSG00000109320', u'ENSG00000092445', u'ENSG00000182578']})
        
        test_2_res = core.fix_jquery_array_keys({u'query_samples[LabID][Variants]': [u'12-00145'], u'panel_context': [u'panel'], u'query_samples[LabID][GeneScore]': [u'12-00145'], u'query_samples[LabID][siRNA]': [u'12-00145']})
        
        self.assertEqual(test_2_res, {'query_samples':{'LabID':{'Variants': [u'12-00145'], 'GeneScore':[u'12-00145'], 'siRNA':[u'12-00145']}}, u'panel_context': [u'panel']})
        
    def test_make_results_table(self):
        
        test_query_nl = core.NodeList()
        
        query_header = ['var', 'gene', 'data', 'query_ind', 'row_id', 'gene_ind']
        
        test_query_nl.add(core.RowNode(['1', 'gene_1', 10, 0, '1_gene_1', 1], query_header))
        test_query_nl.add(core.RowNode(['2', 'gene_1', 15, 0, '2_gene_1', 1], query_header))
        test_query_nl.add(core.RowNode(['2', 'gene_2', -5, 0, '2_gene_2', 1], query_header))
        test_query_nl.add(core.RowNode(['3', 'gene_3', 105, 0, '3_gene_3', 1], query_header))
        
        test_query_nl.attributes['header'] = query_header[0:3]
        
        gene_seeds = make_gene_list([
                        (u'gene_1', u'12-00145', u'GeneScore',10, True),
                        (u'gene_2', u'12-00145',u'GeneScore', 15, True)])
        
        sirna_seeds = make_gene_list([
                        (u'gene_1', u'12-00145', u'siRNA',-2, True),
                        (u'gene_3', u'12-00145',u'siRNA', -3, True),
                        (u'gene_4', u'12-00145',u'siRNA', -2.5, True)])
        
        test_seed_nl = copy.deepcopy(gene_seeds)
        
        test_seed_nl.mergeChildren(sirna_seeds)
        
        
        test_ranking_nl = core.NodeList()
        
        test_ranking_nl.add(core.BasicNode(('gene_1', .5677), only_child=True))
        test_ranking_nl.add(core.BasicNode(('gene_3', .0897), only_child=True))
        test_ranking_nl.add(core.BasicNode(('gene_4', .93455), only_child=True))
        
        test_ranking_sl = core.SeedList(test_ranking_nl)
        
        test_rows, test_header = core.make_results_table(test_query_nl, test_seed_nl, test_ranking_sl)
        
        test_res_header = query_header[0:3]
        test_res_header.extend(sorted(['GeneScore', 'siRNA']) + ['HitWalkerScore', 'HitWalkerRank'])
        
        self.assertEqual(test_header, test_res_header)
        
        test_res_rows = [
            ['1', 'gene_1', 10, 10, -2, .5677, 1],
            ['2', 'gene_1', 15, 10, -2, .5677, 1],
            ['3', 'gene_3', 105, None, -3, .0897,2],
            ['2', 'gene_2', -5, 15, None, None, None]
            ]
        
        self.assertEqual(test_rows, test_res_rows)

    def test_handle_hits(self):
        
        test_nl = core.NodeList()
        
        res_list = []
        
        test_record_obj = collections.namedtuple("Records", ("columns", "values"))
        
        test_vals = [
            (u'ENSG00000254087', u'12-00145', u'siRNA', 'test', -2, True),
            (u'ENSG00000101336', u'12-00145',u'siRNA', 'test', -1, False),
            (u'ENSG00000213341', u'12-00145', u'siRNA','test',-2, True)
        ]
        
        test_records = []
        
        for i in test_vals:
            test_records.append(test_record_obj(**{"columns":('gene', 'sample', 'var', 'type', 'score', 'is_hit'), "values":i}))
        
        core.handle_hits([test_records], test_nl, None)
        
        #It should add GeneNodes to test_nl and also add the appropriate children to the genes
        self.assertTrue(len(test_nl) == 3)
        self.assertEqual(test_nl.ids(), ['ENSG00000254087','ENSG00000101336','ENSG00000213341'])
        
        for i_ind, i in enumerate(test_nl):
            
            i_dict = i.children().json()
            test_row = test_vals[i_ind]
            
            self.assertJSONEqual(i_dict, json.dumps([{'id':test_row[0] + '_' + test_row[2], 'display_name':test_row[0] + '_' + test_row[2],
                                                      'attributes':{'node_type':test_row[2], 'other_nodes':[test_row[1]], 'meta':{'node_cat':'Assay Result', 'type':'test', 'score':test_row[4], 'is_hit':test_row[5]}},
                                                        'children':[]}]))
            
    def test_customize_query(self):
        inp_query1 = {'query':'MATCH (n:Gene{name:{GENE}})-[r:KNOWN_AS]-(m) WHERE r.status="symbol" RETURN n.name,m.name', 'handler':custom_functions.get_gene_names, 'session_params':None}
        
        inp_query2 = {'query':'MATCH(n:LabID{name:{LABID}})-[r:GENE_SCORE_RUN]-()-[r2:SCORE_MAPPED_TO]-(m:Gene{name:{GENE}}) WHERE HAS(r.score) RETURN m.name, r.score,(r.score*r2.modifier) > {GENESCORE} AS is_hit ORDER BY r.score DESC limit 1',
              'handler':lambda x: x, 'session_params':None} 
        
        query_res = core.customize_query(inp_query1, query=lambda x:x.replace("{GENE}", "{name}"))
        self.assertEqual(query_res['query'], 'MATCH (n:Gene{name:{name}})-[r:KNOWN_AS]-(m) WHERE r.status="symbol" RETURN n.name,m.name')
       
        query_res = core.customize_query(inp_query2, query=lambda x: x.replace("{GENE}", "{name}").replace("{GENESCORE}", "{gene_score}").replace("{LABID}", "{GeneScore}"), session_params=lambda x: [['query_samples', 'LabID', 'GeneScore'], ['gene_score']])
        self.assertEqual(query_res['query'], 'MATCH(n:LabID{name:{GeneScore}})-[r:GENE_SCORE_RUN]-()-[r2:SCORE_MAPPED_TO]-(m:Gene{name:{name}}) WHERE HAS(r.score) RETURN m.name, r.score,(r.score*r2.modifier) > {gene_score} AS is_hit ORDER BY r.score DESC limit 1')
        self.assertEqual(query_res['session_params'], [['query_samples', 'LabID', 'GeneScore'], ['gene_score']])

    def test_iterate_dict(self):
        test_dict = {'lev1':{'lev2':{'lev3':'done'}}}
        
        test_res1 = core.iterate_dict(test_dict, ['lev1', 'lev2', 'lev3'])
        self.assertEqual(test_res1, 'done')
        
        test_res1 = core.iterate_dict(test_dict, ['lev1', 'nonlevel2'])
        self.assertEqual(test_res1, None)
        
    def test_specify_type_query_tmp(self):
        
        test_tmpl = {'title':'$$ret_type$$s with siRNA hits for $$result$$','text':'siRNA Hit', 'query':'MATCH(sample:LabID)-[r:SIRNA_RUN]-()-[r2:SIRNA_MAPPED_TO]-(gene:Gene) WHERE  HAS(r.zscore) AND (r.zscore*r2.modifier) < {zscore} AND $$lower_coll_type$$.name IN {$$coll_type$$} WITH $$lower_ret_type$$.name AS ret_type, COLLECT(DISTINCT $$lower_coll_type$$.name) AS use_coll WHERE LENGTH(use_coll) = {$$coll_type$$_length} RETURN ret_type',
                         'handler':None, 'session_params':[['zscore']]}
        
        new_tmpl = core.specify_type_query_tmp(test_tmpl, ret_type='Sample', coll_type='Gene')
        
        self.assertDictEqual(new_tmpl, {'title':'Samples with siRNA hits for $$result$$','text':'siRNA Hit', 'query':'MATCH(sample:LabID)-[r:SIRNA_RUN]-()-[r2:SIRNA_MAPPED_TO]-(gene:Gene) WHERE  HAS(r.zscore) AND (r.zscore*r2.modifier) < {zscore} AND gene.name IN {Gene} WITH sample.name AS ret_type, COLLECT(DISTINCT gene.name) AS use_coll WHERE LENGTH(use_coll) = {Gene_length} RETURN ret_type',
                         'handler':None, 'session_params':[['zscore']]})
    
    
    
        
        
#the workaround here for session modifications is due to http://blog.mediaonfire.com/?p=36 and http://stackoverflow.com/questions/4453764/how-do-i-modify-the-session-in-the-django-test-framework
#class network_view_tests (TestCase):
#    def setUp(self):
#        self.client = Client()
#        self.request_post = {u'siRNA': [u'11-00003'],
#                            u'Variants': [u'11-00003'],
#                            u'var_select': [u'ENSG00000182578', u'ENSG00000185338', u'ENSG00000106070'],
#                            u'GeneScore': [u'None'],
#                            u'seed_select': [u'ENSG00000182578']}
#        
#        settings.SESSION_ENGINE = 'django.contrib.sessions.backends.file'
#        engine = import_module(settings.SESSION_ENGINE)
#        store = engine.SessionStore()
#        store.save()
#        self.session = store
#        self.client.cookies[settings.SESSION_COOKIE_NAME] = store.session_key
#    
#    def tearDown(self):
#        store = self.session
#        os.unlink(store._key_to_file())
#    
#    def test_get_shortest_paths(self):
#        
#        graph_inp = open("/var/www/hitwalker_2_inst/graph_struct.json", "r")
#        graph_struct = json.load(graph_inp)
#        graph_inp.close()
#        
#        session = self.session
#        session['string_conf'] = .4
#        session['query_samples'] = {'LabID':{'siRNA':'11-00003','GeneScore':None, 'Variants':'11-00003'}}
#        session['zscore'] = -2
#        session['gene_score'] = 0
#        session['necessary_vars'] = set(['UNIT_DNA_DIFF', 'Variation', 'IMPACTS', 'DNA_DIFF'])
#        session['where_template'] = '( ( $$IMPACTS$$.Cons_cat = "NonSynonymous" AND $$DNA_DIFF$$.genotype_quality > 40.0 AND $$UNIT_DNA_DIFF$$.QD > 5.0 AND $$UNIT_DNA_DIFF$$.MQ0 < 4.0 AND $$Variation$$.Sample_count < 187.5 ) ) AND ( ( $$Variation$$.in_1kg = 0 AND $$Variation$$.in_dbsnp = 0 ) OR ( $$Variation$$.in_1kg = 1 AND HAS($$Variation$$.gmaf) AND $$Variation$$.gmaf < 0.01 ) OR ( $$Variation$$.in_dbsnp = 1 AND $$Variation$$.in_1kg = 0 AND $$Variation$$.Sample_count = 1.0 ) )'
#        session['graph_struct'] = graph_struct
#        session['hitwalker_score'] = {}
#        session['hit_dict'] = {}
#        
#        session.save()
#        
#        custom_functions.get_shortest_paths (self.client, self.request_post)

class ajax_and_core_utils_tests(TestCase):
    def setUp(self):
        
        self.client = Client()
        settings.SESSION_ENGINE = 'django.contrib.sessions.backends.file'
        engine = import_module(settings.SESSION_ENGINE)
        store = engine.SessionStore()
        store.save()
        self.session = store
        self.client.cookies[settings.SESSION_COOKIE_NAME] = store.session_key
        
        self.graph_db = neo4j.GraphDatabaseService(test_cypher_session+'/db/data/') 
        
        self.session.save()
    
    def tearDown(self):
        store = self.session
        os.unlink(store._key_to_file())
        
     #ajax related methods
    
    def test_get_nodes(self):
        
        #pull out a sample of nodes
        
        #should look like: [{'necessary_vars':set(), 'where_statement':'', 'name':'', 'pretty_where':''}]
        self.session['where_vars'] = []
        self.session.save()
        
        data = neo4j.CypherQuery(self.graph_db, 'MATCH (n) RETURN n.name, LABELS(n) limit 10').execute().data
        
        inp_dict = collections.defaultdict(list)
        
        for i in data:
            inp_dict[i[1][0]].append(i[0])
        
        #make a test config_struct
        
        use_label = inp_dict.keys()[1]
        not_used_label = inp_dict.keys()[0]
        
        test_struct = {
            use_label:[
                {'query':'MATCH (n:'+use_label+'{name:{name}}) RETURN n.name', 'handler':simple_handler, 'session_params':None}
            ]
        }
        
        node_res = core.get_nodes(names=inp_dict[use_label], node_type=use_label, request=self.client, indexed_name="name",  config_struct = test_struct, param_list=[], missing_param="fail", cypher_session = test_cypher_session)
        
        self.assertTrue(isinstance(node_res, core.NodeList))
        self.assertTrue(len(node_res) == len(inp_dict[use_label]))
        self.assertEqual(node_res.ids(), inp_dict[use_label])
        
        ##assert that it will fail if a non referenced label is queried on
        
        self.assertRaises(Exception, core.get_nodes, names=inp_dict[not_used_label], node_type=not_used_label, request=self.client, indexed_name="name",  config_struct = test_struct, param_list=[], missing_param="fail", cypher_session = test_cypher_session)
        
        #assert that it fails upon missing a specified parameter in the query
        
        test_struct_2 = copy.deepcopy(test_struct)
        test_struct_2[use_label][0]['query'] = 'MATCH (n:'+use_label+'{name:{name}}) WHERE HAS(no_var.missing) AND no_var.missing = {MISSING} RETURN n.name'
        
        self.assertRaises(Exception, core.get_nodes, names=inp_dict[not_used_label], node_type=not_used_label, request=self.client, indexed_name="name",  config_struct = test_struct_2, param_list=[], missing_param="fail", cypher_session = test_cypher_session)
        
        ##assert that it can use parameters saved in a session
        
        self.session['param_score'] = 10
        self.session.save()
        
        ##Note that using only_child = True will automatically name i[1] to score.  Will make the cypher query match this
        def simple_handler_w_child(res_list, nodes, request):
            for i in core.BasicResultsIterable(res_list):
                if len(i) > 0:
                    nodes.add(core.BasicNode(i, only_child=True))
                    
        test_struct_3 = copy.deepcopy(test_struct)
        test_struct_3[use_label][0]['query'] = 'MATCH (n:'+use_label+'{name:{name}}) RETURN n.name, {param_score} AS score'
        test_struct_3[use_label][0]['handler'] = simple_handler_w_child
        test_struct_3[use_label][0]['session_params'] = [['param_score']]
                    
        node_res_3 = core.get_nodes(names=inp_dict[use_label], node_type=use_label, request=self.client, indexed_name="name",  config_struct = test_struct_3, param_list=[], missing_param="fail", cypher_session = test_cypher_session)
        
        self.assertTrue(isinstance(node_res_3, core.NodeList))
        self.assertTrue(len(node_res_3) == len(inp_dict[use_label]))
        self.assertEqual(node_res_3.ids(), inp_dict[use_label])
        
        self.assertEqual(node_res_3.summarizeChildren(lambda x: x.getAttr(["attributes","meta", "score"]), len), [1]*len(inp_dict[use_label]))
        self.assertEqual(node_res_3.summarizeChildren(lambda x: x.getAttr(["attributes","meta", "score"]), max), [10]*len(inp_dict[use_label]))
        
        ##assert that parameters can be supplied via param_list
        
        test_struct_3[use_label][0]['session_params'] = None
        
        node_res_4 = core.get_nodes(names=inp_dict[use_label], node_type=use_label, request=self.client, indexed_name="name",  config_struct = test_struct_3, param_list=[{'param_score':10}], missing_param="fail", cypher_session = test_cypher_session)
        
        self.assertJSONEqual(node_res_3.json(), node_res_4.json())
        
        ##assert that it can inject parameters into the queries in a basic sense using 'where_vars'
        
        ####also need to add in graph_struct, make a fake one???
        self.session['graph_struct'] = {}
        self.session['where_vars'] = [{'necessary_vars':set([use_label]), 'where_statement':'$$'+use_label+'$$.name = "'+inp_dict[use_label][0]+'"'}]
        self.session.save()
        
        node_res_5 = core.get_nodes(names=inp_dict[use_label], node_type=use_label, request=self.client, indexed_name="name",  config_struct = test_struct, param_list=[], missing_param="fail", cypher_session = test_cypher_session)
        
        self.assertTrue(isinstance(node_res_5, core.NodeList))
        self.assertTrue(len(node_res_5) == 1)
        self.assertEqual(node_res_5.ids(), [inp_dict[use_label][0]])
    
    def test_copy_nodes(self):
        
        self.session['where_vars'] = []
        self.session.save()
        
        data = neo4j.CypherQuery(self.graph_db, 'MATCH (n) RETURN DISTINCT LABELS(n) limit 2').execute().data
        
        use_labels = map(lambda x: x.values[0][0],data)
        
        #pull some (somewhat) randomly from both labels
        
        inp_dict = collections.defaultdict(list)
        
        for i in use_labels:
            for j in neo4j.CypherQuery(self.graph_db, 'MATCH (n:'+i+') RETURN n.name LIMIT 100').execute().data:
                inp_dict[i].append(j.values[0])
        
        test_struct = {'nodes':{}, 'edges':{}}
        
        for i in use_labels:
            test_struct['nodes'][i] = [{'query':'MATCH (n:'+i+'{name:{name}}) RETURN n.name, 10,"'+i+'"', 'handler':simple_handler_w_child, 'session_params':None}]
        
        #need to specify handlers for the edges
        
        def make_fixed_edges(query_genes, subj_genes, request, config_struct_nodes, cur_graph):
            
            cur_graph['links'] = [
                {'source':10, 'target':8, 'attributes':{'type':'test'}},
                {'source':6, 'target':5, 'attributes':{'type':'test'}}
            ]
            
            return cur_graph
        
        def no_rels(query_genes, subj_genes, request, config_struct_nodes, cur_graph):
            return cur_graph
        
        test_struct['edges'] = {
            use_labels[0]:{
                    use_labels[0]:{'query':None, 'handler':no_rels, 'session_params':None},
                    use_labels[1]:{'query':None, 'handler':make_fixed_edges, 'session_params':None}
                },
            use_labels[1]:{
                use_labels[1]:{'query':None, 'handler':no_rels, 'session_params':None}
            }}
        
        test_1_subj = map(lambda x: {'id':x, 'node_type':use_labels[0]}, inp_dict[use_labels[0]][:10])
        test_1_query = map(lambda x: {'id':x, 'node_type':use_labels[1]}, inp_dict[use_labels[1]][:10])
        
        ##Note we will only run copy_nodes with never_group = True, as the apply_grouping tests will cover that one as it can be run on its own...
        
        res_g_1 = core.copy_nodes(test_1_subj, test_1_query, request=self.client, query_dict=test_struct, never_group=True)
        
        self.assertTrue(len(res_g_1['nodes']) == 20)
        self.assertTrue(len(res_g_1['links']) == 2)
        
        #now try again but with a few of the same nodes as both the subject and query
        
        print test_1_subj + copy.deepcopy(test_1_query[5:7])
        
        res_g_2 = core.copy_nodes(test_1_subj + copy.deepcopy(test_1_query[5:7]), test_1_query, request=self.client, query_dict=test_struct, never_group=True)
        
        self.assertTrue(len(res_g_2['nodes']) == 20)
        self.assertTrue(len(res_g_2['links']) == 2)
        
        #need to finish this and also add in a check to ensure the children are appropriate
        
        #try with the exact same subjects and queries
        
        
        
        #try with subject empty
        
        #Finally try with query empty
        
    
    def check_node_ids(self, gene_node_res, spec_nodes, verbose=False):
        
        if isinstance(gene_node_res['nodes'], core.NodeList):
            use_nodes = gene_node_res['nodes'].tolist()
        else:
            use_nodes = gene_node_res['nodes']
        
        found_nodes = set()
        for i in use_nodes:
            found_nodes.add(i["id"])
            
        if verbose == True:
            print "**************"
            print spec_nodes, found_nodes
            print "**************"
        self.assertTrue(len(spec_nodes.intersection(found_nodes)) == len(spec_nodes.union(found_nodes)))
    
    def check_links(self, gene_node_res, true_links, verbose=False):
        
        if isinstance(gene_node_res['nodes'], core.NodeList):
            use_nodes = gene_node_res['nodes'].tolist()
        else:
            use_nodes = gene_node_res['nodes']
        
        found_links = 0
        
        if verbose == True:
            print "**************"
            print map(lambda x: x["id"], use_nodes)
        
        for i in gene_node_res['links']:
                did_find = False
                if (use_nodes[i["source"]]["id"] + '.' + use_nodes[i["target"]]["id"] + '.' + i['attributes']['type'] in true_links) or (use_nodes[i["target"]]["id"] + '.' + use_nodes[i["source"]]["id"] + '.' + i['attributes']['type'] in true_links):
                    found_links += 1
                    did_find = True
                if verbose ==True:
                    print i, did_find
        if verbose == True:
            print "**************"
       
        self.assertTrue(found_links == len(true_links) and len(gene_node_res['links']) == found_links)
    
    #core utils methods
    
    def test_apply_grouping(self):
        
        self.skipTest('Fix grouping and such related bugs/ambiguities')
        
        samp_resp = self.client.post('/HitWalker2/match_sample/', {'query':'12-'})
        
        all_samps = map(lambda x: x['text'], json.loads(samp_resp.content)['data']['results'])
        
        #usage for fullfill_node_query
        
        #get < 20 sample nodes and no query_nodes
        
        test_1_nl = core.get_nodes(all_samps[0:15], 'LabID', self.client)
        
        test_res_1 = core.apply_grouping({'nodes':test_1_nl, 'links':[]}, [])
        
        self.check_node_ids(test_res_1, set(all_samps[0:15]))
        self.assertTrue(len(test_res_1['links']) == 0)
        
        #get >20 samples nodes and no query_nodes
        
        test_2_nl = core.get_nodes(all_samps[0:30], 'LabID', self.client)
        test_res_2 = core.apply_grouping({'nodes':test_2_nl, 'links':[]}, [])
        
        self.check_node_ids(test_res_2, set(['meta_node_1']))
        self.assertTrue(len(test_res_2['links']) == 0)
        
        #usage for copy_nodes
        
        #< 20 samples with 1 query
        
        test_3_start = self.client.post('/HitWalker2/fullfill_node_query/', {'choice':'Gene-1', 'nodes':'{"Gene": [{"display_name": "TYRO3", "id": "ENSG00000092445"}]}'})
        
        test_3_start_dict = json.loads(test_3_start.content)
        
        #this is a metanode with 26 children
        
        #first add the gene back to it and make sure that the metanode is connected to it via a single edge
        
        get_children = map(lambda x:x['id'],  test_3_start_dict['graph']['nodes'][0]['children'])
        #
        test_3_nl = core.get_nodes(get_children, 'LabID', self.client)
        test_3_nl.extend(core.get_nodes(['ENSG00000092445'], 'Gene', self.client, missing_param="skip"))
        
        test_3_graph = custom_functions.gene_to_lab_id (['ENSG00000092445'], get_children, self.client, config.edge_queries['nodes'], {'nodes':test_3_nl, 'links':[]})
        
        #first check this with only siRNA type relationships...
        
        test_3_graph['links'] = filter(lambda x:x['attributes']['type'] == "Observed_siRNA", test_3_graph['links'] )
        
        test_res_3 = core.apply_grouping(test_3_graph, ['ENSG00000092445'])
        
        self.check_node_ids(test_res_3, set(['meta_node_1', 'ENSG00000092445']))
        self.check_links(test_res_3, set(['ENSG00000092445.meta_node_1.Observed_siRNA']), verbose=True)
        #
        #then try with siRNA and gene score
        
        test_4_nl = core.get_nodes(get_children, 'LabID', self.client)
        test_4_nl.extend(core.get_nodes(['ENSG00000092445'], 'Gene', self.client, missing_param="skip"))
        test_4_graph = custom_functions.gene_to_lab_id (['ENSG00000092445'], get_children, self.client, config.edge_queries['nodes'], {'nodes':test_4_nl, 'links':[]})
        
        test_4_graph['links'] = filter(lambda x:x['attributes']['type'] in set(["Observed_siRNA", "Observed_GeneScore"]), test_4_graph['links'] )
        
        
        #siRNA should connect to everyone, tally which nodes also have gene scores
        gs_nodes = set()
        
        q_node = test_4_graph['nodes'].nodeIndex('ENSG00000092445')
        
        for i in test_4_graph['links']:
            if i['attributes']['type'] == 'Observed_GeneScore':
                if i['source'] ==  q_node:
                    gs_nodes.add(test_4_graph['nodes'].getByIndex(i['target']).id)
                else:
                    gs_nodes.add(test_4_graph['nodes'].getByIndex(i['source']).id)
        
        test_res_4 = core.apply_grouping(test_4_graph, ['ENSG00000092445'])
        
        self.check_node_ids(test_res_4, set(['meta_node_1', 'ENSG00000092445']).union(gs_nodes))
        
        base_set = set(['ENSG00000092445.meta_node_1.Observed_siRNA'])
        base_set = base_set.union(set(map(lambda x: 'ENSG00000092445.'+x+'.Observed_GeneScore', gs_nodes)))
        base_set = base_set.union(set(map(lambda x: 'ENSG00000092445.'+x+'.Observed_siRNA', gs_nodes)))
        
        self.check_links(test_res_4, base_set)
        
        #now add in variants to the mix
        
        test_5_nl = core.get_nodes(get_children, 'LabID', self.client)
        test_5_nl.extend(core.get_nodes(['ENSG00000092445'], 'Gene', self.client, missing_param="skip"))
        test_5_graph = custom_functions.gene_to_lab_id (['ENSG00000092445'], get_children, self.client, config.edge_queries['nodes'], {'nodes':test_5_nl, 'links':[]})
        
        #figure out which samples have variants
        
        var_nodes = set()
        
        q_node = test_5_graph['nodes'].nodeIndex('ENSG00000092445')
        
        for i in test_5_graph['links']:
            if i['attributes']['type'] == 'Observed_Variants':
                if i['source'] ==  q_node:
                    var_nodes.add(test_5_graph['nodes'].getByIndex(i['target']).id)
                else:
                    var_nodes.add(test_5_graph['nodes'].getByIndex(i['source']).id)
        
        self.assertTrue(len(var_nodes.intersection(gs_nodes)) == len(var_nodes))
        
        test_res_5 = core.apply_grouping(test_5_graph, ['ENSG00000092445'])
        
        #should be the same as 4
        self.check_node_ids(test_res_5, set(['meta_node_1', 'ENSG00000092445']).union(gs_nodes))
        
        #links should be the same as 4 plus the variant nodes
        self.check_links(test_res_5, base_set.union(set(map(lambda x: 'ENSG00000092445.'+x+'.Observed_Variants', var_nodes))))
        
        #test 6 is the same as test 5 but with several query genes instead of 1.
        
        test_6_nl = core.get_nodes(get_children, 'LabID', self.client)
        test_6_nl.extend(core.get_nodes(['ENSG00000092445'], 'Gene', self.client, missing_param="skip"))
        test_6_nl.extend(core.get_nodes(['ENSG00000010327'], 'Gene', self.client, missing_param="skip"))
        test_6_graph = custom_functions.gene_to_lab_id (['ENSG00000092445', 'ENSG00000010327'], get_children, self.client, config.edge_queries['nodes'], {'nodes':test_6_nl, 'links':[]})
        
        var_nodes2_map = collections.defaultdict(set)
        
        q_node = test_6_graph['nodes'].nodeIndex('ENSG00000010327')
        
        #here, need to figure out what the different types of relationships are for the new node
        for i in test_6_graph['links']:
            if i['source'] ==  q_node:
                var_nodes2_map[i['attributes']['type']].add(test_6_graph['nodes'].getByIndex(i['target']).id)
            elif i['target'] ==  q_node:
                var_nodes2_map[i['attributes']['type']].add(test_6_graph['nodes'].getByIndex(i['source']).id)
        
        #from test 5
        six_set = base_set.union(set(map(lambda x: 'ENSG00000092445.'+x+'.Observed_Variants', var_nodes)))
        
        for key,val in var_nodes2_map.items():
            six_set = six_set.union(set(map(lambda x: 'ENSG00000010327.'+x+'.'+key, val)))
        
        #also add in any relationships with ENSG00000092445
        
        dist_nodes_6 = reduce(set.union, var_nodes2_map.values())
        
        q_node = test_6_graph['nodes'].nodeIndex('ENSG00000092445')
        
        dist_nodes_6_pos = set()
        for i in dist_nodes_6:
            dist_nodes_6_pos.add(test_6_graph['nodes'].nodeIndex(i))
        
        #here, need to figure out what the different types of relationships are for the new node
        for i in test_6_graph['links']:
            if i['source'] ==  q_node and i['target'] in dist_nodes_6_pos:
                six_set.add('ENSG00000092445.'+test_6_graph['nodes'].getByIndex(i['target']).id+'.'+i['attributes']['type'])
                print i
            elif i['target'] ==  q_node and i['source'] in dist_nodes_6_pos:
                six_set.add('ENSG00000092445.'+test_6_graph['nodes'].getByIndex(i['source']).id+'.'+i['attributes']['type'])
                print i
                
        for i in six_set:
            print i
        
        test_res_6 = core.apply_grouping(test_6_graph, ['ENSG00000092445', 'ENSG00000010327'])
        
        self.check_node_ids(test_res_6, set(['meta_node_1', 'ENSG00000092445', 'ENSG00000010327']).union(gs_nodes).union(dist_nodes_6))
        
        self.check_links(test_res_6, six_set)
        
    def get_exp_nodes_edges(self, use_graph, query_ids):
        
        edge_set = set()
        node_set = set(query_ids)
        
        link_nodes = set()
        
        for i in query_ids:
            q_node = use_graph['nodes'].nodeIndex(i)
            edge_col = collections.defaultdict(set)
            for j in use_graph['links']:
                if j['source'] == q_node:
                    edge_col[j['attributes']['type']].add(use_graph['nodes'].getByIndex(j['target']).id)
                elif j['target'] == q_node:
                    edge_col[j['attributes']['type']].add(use_graph['nodes'].getByIndex(j['source']).id)
                else:
                    node_set.add(use_graph['nodes'].getByIndex(j['source']).id)
                    node_set.add(use_graph['nodes'].getByIndex(j['target']).id)
            
            for key,vals in edge_col.items():
                link_nodes = set.union(link_nodes, vals)
                if len(vals) > config.max_nodes:
                    edge_set.add(i + '.meta_node_1.'+ key)
                    node_set.add('meta_node_1')
                else:
                    for k in vals:
                        edge_set.add(i + '.' + k + '.' + key)
                        node_set.add(k)
        
        lo_nodes = set(map(lambda x: x['id'], use_graph['nodes'].tolist()))
        
        if len(lo_nodes) > config.max_nodes:
            node_set.add('meta_node_1')
        else:
            node_set = set.union(node_set, lo_nodes)
        
        return node_set, edge_set
#//*[@cx ='120']
        
        


class query_parser_tests (TestCase):
    def setUp(self):
        
        graph_inp = open("/var/www/hitwalker_2_inst/graph_struct.json", "r")
        graph_struct = json.load(graph_inp)
        graph_inp.close()
        
        self.graph_struct = graph_struct
        self.necessary_vars = set(['UNIT_DNA_DIFF', 'Variation', 'IMPACTS', 'DNA_DIFF'])
        self.where_template = '( ( $$IMPACTS$$.Cons_cat = "NonSynonymous" AND $$DNA_DIFF$$.genotype_quality > 40.0 AND $$UNIT_DNA_DIFF$$.QD > 5.0 AND $$UNIT_DNA_DIFF$$.MQ0 < 4.0 AND $$Variation$$.Sample_count < 187.5 ) ) AND ( ( $$Variation$$.in_1kg = 0 AND $$Variation$$.in_dbsnp = 0 ) OR ( $$Variation$$.in_1kg = 1 AND HAS($$Variation$$.gmaf) AND $$Variation$$.gmaf < 0.01 ) OR ( $$Variation$$.in_dbsnp = 1 AND $$Variation$$.in_1kg = 0 AND $$Variation$$.Sample_count = 1.0 ) )'
    
    def leave_same(self, x):
        return x
    
    def param_test_1(self):
        return {'Variant_Filters':{'type':'grouped',
                       'fields':{
                            'freq':{'type':'numeric', 'comparison':'<','default':.01, 'range':[0,1], 'name':'Global MAF', 'var_name':'gmaf','required':{'from':'Variation'}, 'trans':self.leave_same, 'needs_has':''},
                            'cohort_freq':{'type':'numeric', 'comparison':'<', 'default':.5, 'range':[0,1], 'name':'Cohort Alt. Frequency', 'var_name':'Sample_count', 'required':{'from':'Variation'}, 'trans':self.leave_same},
                            'cohort_count':{'type':'numeric', 'comparison':'=', 'default':1, 'range':[0,500], 'name':'Cohort Count', 'var_name':'Sample_count', 'required':{'from':'Variation'}, 'trans':self.leave_same},
                            'in_1kg':{'type':'character', 'default':'False', 'range':['True', 'False'], 'name':'In 1000 genomes','var_name':'in_1kg','required':{'from':'Variation'}, 'trans':self.leave_same},
                            'in_dbsnp':{'type':'character', 'default':'False', 'range':['True', 'False'], 'name':'In dbSNP','var_name':'in_dbsnp', 'required':{'from':'Variation'},'trans':self.leave_same},
                            'allele_count': {'type':'numeric', 'range':[0,2], 'name':'Genotype','var_name':'allele_count','required':{'from':'DNA_DIFF'},'trans':self.leave_same},
                            'genotype_quality': {'type':'numeric', 'comparison':'>', 'default':40, 'range':[0,100], 'name':'Genotype Quality', 'var_name':'genotype_quality', 'required':{'from':'DNA_DIFF'},'trans':self.leave_same},
                            'depth':{'type':'numeric', 'range':[0,100000], 'name':'Read Depth','var_name':'depth', 'required':{'from':'DNA_DIFF'},'trans':self.leave_same},
                            
                            'FS': {'type':'numeric', 'range':[-10000,10000],'name':'Fisher Strand','var_name':'FS', 'required':{'from':'UNIT_DNA_DIFF'},'trans':self.leave_same},
                            'MQ0':{'type':'numeric', 'comparison':'<', 'default':4, 'range':[0,10000], 'name':'Number Ambigous Reads','var_name':'MQ0', 'required':{'from':'UNIT_DNA_DIFF'},'trans':self.leave_same},
                            'MQ': {'type':'numeric', 'range':[0,100], 'name':'Mapping Quality','var_name':'MQ', 'required':{'from':'UNIT_DNA_DIFF'},'trans':self.leave_same},
                            'QD': {'type':'numeric', 'comparison':'>', 'default':5, 'range':[0,10000], 'name':'Quality / Depth','var_name':'QD', 'required':{'from':'UNIT_DNA_DIFF'},'trans':self.leave_same},
                            'SB': {'type':'numeric', 'range':[-10000,10000], 'name':'Strand Bias','var_name':'SB', 'required':{'from':'UNIT_DNA_DIFFERENCE'},'trans':self.leave_same},
                            
                            'Cons_cat':{'type':'character', 'default':'NonSynon.', 'range':['Synonymous', 'NonSynon.', 'Other'], 'var_name':'Cons_cat', 'name':'Consequence', 'required':{'from':'IMPACTS'}, 'trans':self.leave_same}
                    },
                    'default_groups':[[
                        [{'field':'Cons_cat'},
                         {'field':'genotype_quality', 'logical':'AND'}]
                        ],
                        [
                            [{'field':'in_1kg', 'default':'False', 'logical':'AND'},{'field':'in_dbsnp','default':'False','logical':'AND'}],
                            [{'field':'in_1kg', 'default':'True', 'logical':'OR'}, {'field':'freq', 'logical':'AND'}]
                        ]]}}, '(($$IMPACTS$$.Cons_cat = "NonSynon." AND $$DNA_DIFF$$.genotype_quality > 40)) AND (($$Variation$$.in_1kg = "False" AND $$Variation$$.in_dbsnp = "False") OR ($$Variation$$.in_1kg = "True" AND HAS($$Variation$$.gmaf) AND $$Variation$$.gmaf < 0.01))'
    
    def make_logical_list(self, inp_params):
        trans_funcs= {}
        
        for i in inp_params.keys():
            for j in inp_params[i]['fields'].keys():
                if inp_params[i]['fields'][j].has_key('trans'):
                    trans_funcs[j] = inp_params[i]['fields'][j].pop('trans') 
            if inp_params[i]['type'] == 'grouped':
                logical_list = reduce (lambda x,y: x+y, reduce (lambda x,y: x+y,inp_params[i]['default_groups']))
                inp_params[i]['logical_list'] = map(lambda x: 'null' if x.has_key('logical')==False else x['logical'], logical_list)
                
        return trans_funcs, inp_params
    
    def test_parse_parameters(self):
        
        #set up a session
        
        self.client = Client()
        settings.SESSION_ENGINE = 'django.contrib.sessions.backends.file'
        engine = import_module(settings.SESSION_ENGINE)
        store = engine.SessionStore()
        store.save()
        self.session = store
        self.client.cookies[settings.SESSION_COOKIE_NAME] = store.session_key
        
        param_test1_param, param_test_res = self.param_test_1()
        
        
        params_1_tf, params_1 = self.make_logical_list(copy.deepcopy(param_test1_param))
        
        #need to adjust the logical values as in index/views.py
        
        nec_where_1 = core.parse_parameters(params_1, params_1_tf, self)
        
        #the standard parameters should just be added to the session
        #but the below is for this case...
        #not using this test for the first one...
        #for i,j in params_1.items():
        #    if j['type'] == 'standard':
        #        for k,l in j['fields'].items():
        #            self.assertEqual(self.client.session[k], l['comparison'] + ' ' + str(l['default']))
        
        self.assertEqual(nec_where_1[0]['where_statement'].replace(" ", ""), param_test_res.replace(" ", ""))
        
        nec_vars_1 = set(list(re.findall(r"\$\$(\w+)\$\$", param_test_res)))
        
        self.assertEqual(nec_where_1[0]['necessary_vars'], nec_vars_1)
        
        store = self.session
        os.unlink(store._key_to_file())
    
    def test_check_input_query_with(self):
        temp_str= '''MATCH (meth)-[r:METHYL_MAPPED_TO]-(gene:EntrezID) WHERE HAS(r.score) AND r.score < .05 WITH gene MATCH (subject:MergeID)-[d:DERIVED]-(s) WHERE subject.name IN ["103051"]
        WITH gene, COLLECT(DISTINCT subject.name) AS use_coll WHERE LENGTH(use_coll) = 1  MATCH (n:MergeID)-[d:DERIVED]-(samp) WHERE n.name IN ["103051"] WITH n, gene
        MATCH (meth)-[r:METHYL_MAPPED_TO]-(gene) WHERE HAS(r.score) WITH n, gene, MIN(r.score) AS methyl_pvalue RETURN gene.name AS gene, n.name AS sample, "Methylation" AS var,
        methyl_pvalue AS score, methyl_pvalue < .05 AS is_hit'''
        
        check_vars = set(['gene'])
        
        res_str = core.check_input_query_with(temp_str, check_vars)
        
        self.assertEqual(temp_str, res_str)
        
        temp_str2= '''MATCH (meth)-[r:METHYL_MAPPED_TO]-(gene:EntrezID) WHERE HAS(r.score) AND r.score < .05 WITH gene MATCH (subject:MergeID)-[d:DERIVED]-(s) WHERE subject.name IN ["103051"]
        WITH gene, COLLECT(DISTINCT subject.name) AS use_coll WHERE LENGTH(use_coll) = 1  MATCH (n:MergeID)-[d:DERIVED]-(samp) WHERE n.name IN ["103051"] WITH n
        MATCH (meth)-[r:METHYL_MAPPED_TO]-(gene) WHERE HAS(r.score) WITH n, gene, MIN(r.score) AS methyl_pvalue RETURN gene.name AS gene, n.name AS sample, "Methylation" AS var,
        methyl_pvalue AS score, methyl_pvalue < .05 AS is_hit'''
        
        res_str2 = core.check_input_query_with(temp_str2, check_vars)
        
        goal_str = '''MATCH (meth)-[r:METHYL_MAPPED_TO]-(gene:EntrezID) WHERE HAS(r.score) AND r.score < .05 WITH gene MATCH (subject:MergeID)-[d:DERIVED]-(s) WHERE subject.name IN ["103051"]
        WITH gene, COLLECT(DISTINCT subject.name) AS use_coll WHERE LENGTH(use_coll) = 1  MATCH (n:MergeID)-[d:DERIVED]-(samp) WHERE n.name IN ["103051"] WITH gene,n
        MATCH (meth)-[r:METHYL_MAPPED_TO]-(gene) WHERE HAS(r.score) WITH n, gene, MIN(r.score) AS methyl_pvalue RETURN gene.name AS gene, n.name AS sample, "Methylation" AS var,
        methyl_pvalue AS score, methyl_pvalue < .05 AS is_hit'''
        
        self.assertEqual(res_str2, goal_str)
    
    #checks both add_where_input_query and check_input_query_where which is called internally
    def test_add_where_input_query(self):
        
        self.skipTest('Fix grouping and such related bugs/ambiguities')
        
        cypher_header = ['var.name', 'var.var_type', 'ROUND(r.allele_count)', 'q.name', 'p.name', 'o.name','r2.Amino_acids', 'r2.Protein_position', 'r2.SIFT', 'r2.PolyPhen', 'REPLACE(RTRIM(REDUCE(str="",n IN r2.Consequence|str+n+" ")), " ", ";")', 'var.freq', 'REPLACE(RTRIM(REDUCE(str="",n IN var.Existing_variation|str+n+" ")), " ", ";")', 'REPLACE(REPLACE(STR(var.in_dbsnp), "1", "True"), "0", "False")', 'REPLACE(REPLACE(STR(var.in_1kg), "1", "True"), "0", "False")']
        query_str_1 = 'MATCH (labid:LabID{name:{samp_id}})-[:ALIAS_OF]-(sample)-[r:DNA_DIFF]-(var)-[u:UNIT_DNA_DIFF]-(exp) WHERE (exp)<-[:GENOTYPED_USING]-(sample) WITH sample,u,r,var MATCH (var)-[r2:IMPACTS]->(o)-[:TRANSCRIBED]-(p)-[r4:KNOWN_AS]-(q) WHERE r4.status = "symbol"  RETURN ' + string.joinfields(cypher_header, ",")

        cor_where_template = self.where_template[:]
        
        for i,j in [["$$IMPACTS$$", "r2"],["$$DNA_DIFF$$", "r"], ["$$UNIT_DNA_DIFF$$", "u"], ["$$Variation$$", "var"]]:
           cor_where_template = cor_where_template.replace(i, j)
            
        cor_query_str_1 = query_str_1.replace('WHERE r4.status = "symbol"', 'WHERE' + cor_where_template +' AND r4.status = "symbol" ')
        
        rep_query = core.add_where_input_query (query_str_1, self.where_template, self.necessary_vars,self.graph_struct)
        
        self.assertEqual(rep_query.replace(" ", ""), cor_query_str_1.replace(" ", ""))
        
        #then check multi-starting point query
        
        query_str_2 = 'MATCH (gene:Gene{name:{name}})-[:TRANSCRIBED]->(trans)<-[impacts_r:IMPACTS]-(var) WITH gene,trans,var,impacts \
                      MATCH (labid:LabID{name:{Variants}})<-[alias_r:ALIAS_OF]-(samp)-[dna_diff_r:DNA_DIFF]->(var)<-[unit_dna_diff_r:UNIT_DNA_DIFF]-(exp) \
                      WHERE (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype" RETURN var, trans,impacts_r,dna_diff_r,unit_dna_diff_r,exp'
        
        cor_where_template = self.where_template[:]
        
        for i,j in [["$$IMPACTS$$", "impacts_r"],["$$DNA_DIFF$$", "dna_diff_r"], ["$$UNIT_DNA_DIFF$$", "unit_dna_diff_r"], ["$$Variation$$", "var"]]:
           cor_where_template = cor_where_template.replace(i, j)
        
        cor_query_str_2 = query_str_2.replace('WHERE (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype"', 'WHERE ' + cor_where_template + ' AND (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype"')
        
        rep_query = core.add_where_input_query(query_str_2, self.where_template, self.necessary_vars,self.graph_struct)
        
        self.assertEqual(rep_query.replace(" ", ""), cor_query_str_2.replace(" ", ""))
        
        #check query that does not utilize all of the necessary_vars--should simply pass the query_str through
        
        query_str_3 = 'MATCH (gene:Gene{name:{name}})-[:TRANSCRIBED]->(trans)<-[impacts_r:IMPACTS]-(var) RETURN gene,trans,var'
        
        rep_query = core.add_where_input_query(query_str_3, self.where_template, self.necessary_vars,self.graph_struct)
        
        self.assertEqual(rep_query, query_str_3)
        
        query_str_4 = 'MATCH (gene:Gene{name:{name}})-[:TRANSCRIBED]->(trans)<-[impacts_r:IMPACTS]-(var) WITH gene,trans,var,impacts \
                      MATCH (labid:LabID)<-[alias_r:ALIAS_OF]-(samp)-[dna_diff_r:DNA_DIFF]->(var)<-[unit_dna_diff_r:UNIT_DNA_DIFF]-(exp) \
                      WHERE (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype" RETURN var, trans,impacts_r,dna_diff_r,unit_dna_diff_r,exp'
        
        cor_where_template = self.where_template[:]
        
        for i,j in [["$$IMPACTS$$", "impacts_r"],["$$DNA_DIFF$$", "dna_diff_r"], ["$$UNIT_DNA_DIFF$$", "unit_dna_diff_r"], ["$$Variation$$", "var"]]:
           cor_where_template = cor_where_template.replace(i, j)
        
        cor_query_str_4 = query_str_4.replace('WHERE (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype"', 'WHERE '+cor_where_template+ ' AND (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype" ')
        
        rep_query = core.add_where_input_query(query_str_4, self.where_template, self.necessary_vars,self.graph_struct)
        
        self.assertEqual(rep_query.replace(" ", ""), cor_query_str_4.replace(" ", ""))
        
        query_str_5 = 'MATCH (gene:Gene)-[:TRANSCRIBED]->(trans)<-[impacts_r:IMPACTS]-(var) WHERE gene.name IN ["ENSG00000196132","ENSG00000010327"] WITH gene,trans,var,impacts_r \
                        MATCH (sample:LabID)<-[alias_r:ALIAS_OF]-(samp)-[dna_diff_r:DNA_DIFF]->(var)<-[unit_dna_diff_r:UNIT_DNA_DIFF]-(exp) WHERE (exp)<-[:GENOTYPED_USING]-(samp) AND alias_r.alias_type="genotype" \
                        WITH sample.name AS ret_type, COLLECT(DISTINCT gene.name) AS use_coll WHERE LENGTH(use_coll) = 2 RETURN ret_type'
        
        cor_where_template = self.where_template[:]
        
        for i,j in [["$$IMPACTS$$", "impacts_r"],["$$DNA_DIFF$$", "dna_diff_r"], ["$$UNIT_DNA_DIFF$$", "unit_dna_diff_r"], ["$$Variation$$", "var"]]:
           cor_where_template = cor_where_template.replace(i, j)
        
        cor_query_str_5 = query_str_5.replace("WHERE (exp)<-[:GENOTYPED_USING]-(samp)", "WHERE "+ cor_where_template + "AND (exp)<-[:GENOTYPED_USING]-(samp)")
        
        rep_query = core.add_where_input_query (query_str_5, self.where_template, self.necessary_vars,self.graph_struct)
        
        self.assertEqual(rep_query.replace(" ", ""), cor_query_str_5.replace(" ", ""))
        
        query_str_6 = "MATCH (n:LabID{name:{name}})<-[:PRODUCED]-(m)<-[:HAS_DISEASE]-(k) RETURN n,m,k"
        
        rep_query = core.add_where_input_query (query_str_6,self.where_template, self.necessary_vars,self.graph_struct)
        
        self.assertEqual(rep_query, query_str_6)
        
        query_str_7 = 'MATCH (n:Gene{name:{name}})-[r:KNOWN_AS]-(m) WHERE r.status="symbol" WITH n,m MATCH (n)<-[:EXTERNAL_ID]-()<-[:PATHWAY_CONTAINS]-(p) RETURN n.name,m.name, COLLECT(DISTINCT p.name)'
        
        rep_query = core.add_where_input_query (query_str_7,self.where_template, self.necessary_vars,self.graph_struct)
        
        self.assertEqual(rep_query, query_str_7)

#as these have been added to core, can be used there...
class CustomFunctionsTests(TestCase):
    
    def setUp(self):
        
        #make a test session
        
        self.client = Client()
        settings.SESSION_ENGINE = 'django.contrib.sessions.backends.file'
        engine = import_module(settings.SESSION_ENGINE)
        store = engine.SessionStore()
        store.save()
        self.session = store
        self.client.cookies[settings.SESSION_COOKIE_NAME] = store.session_key
        
        graph_inp = open("/var/www/hitwalker_2_inst/graph_struct.json", "r")
        graph_struct = json.load(graph_inp)
        graph_inp.close()
        
        #self.session['graph_struct'] = graph_struct
        #self.session['zscore'] = -2
        #self.session['gene_score'] = 10
        #self.session['string_conf'] = .4
        #self.session['hitwalker_score'] = {}
        
        self.session.save()
        
        #Gene 1 -> many
        #   neotool cypher 'MATCH (gene:Gene)-[:EXTERNAL_ID]-()-[:MAPPED_TO]-(string_from) WITH gene.name AS g, COUNT(string_from.name) AS g_c WHERE g_c > 1 RETURN g  limit 10'
        #   ENSG00000240764 | 9606.ENSP00000252085 
        #   ENSG00000240764 | 9606.ENSP00000312070 
        #Protein 1 -> many
        #   neotool cypher 'MATCH (gene:Gene)-[:EXTERNAL_ID]-()-[:MAPPED_TO]-(string_from) WITH string_from.name AS p, COUNT(gene.name) AS p_c WHERE p_c > 1 RETURN p  limit 10'
        #   ENSG00000263328 | 9606.ENSP00000279783 
        #   ENSG00000150261 | 9606.ENSP00000279783
        
        sirna_seeds = make_gene_list([
                        (u'ENSG00000254087', u'12-00145', u'siRNA',-2, True),
                        (u'ENSG00000101336', u'12-00145',u'siRNA', -1, True),
                        (u'ENSG00000213341', u'12-00145', u'siRNA',-2, True),
                        (u'ENSG00000240764', u'12-00145', u'siRNA', -3, True)])
        
        gene_seeds = make_gene_list([
                        (u'ENSG00000254087', u'12-00145', u'GeneScore',10, True),
                        (u'ENSG00006898999', u'12-00145',u'GeneScore', 15, True),
                        (u'ENSG00000213765', u'12-00145', u'GeneScore',30, True),
                        (u'ENSG00000263328',  u'12-00145', u'GeneScore',40, True),
                        (u'ENSG00000150261',  u'12-00145', u'GeneScore',50, True)])
        
        sirna_seeds.mergeChildren(gene_seeds)
        
        self.test_seed_list_1 = sirna_seeds
        
    def tearDown(self):
        store = self.session
        os.unlink(store._key_to_file())
    
    #def make_gene_list(self,res_list):
    #    nodes = core.NodeList()
    #    
    #    seed_header = ['gene', 'sample', 'var', 'score', 'is_hit']
    #    
    #    for ind, gene_result in enumerate(res_list):
    #        
    #        if nodes.hasNode(gene_result[0]) == False:
    #            #create the gene prior to adding if it doesn't exist 
    #            nodes.add(core.GeneNode([gene_result[0], gene_result[0], []]))
    #        
    #        seed_node = core.SeedNode(nodes.getNode(gene_result[0]), gene_result, seed_header)
    #        
    #        #then add the hits to the GeneNode
    #        nodes.addChild(gene_result[0], seed_node)
    #            
    #    return nodes
    
    @unittest.skip("Working on graph db")
    def test_combine_sirna_gs(self):
        
        #should return a SeedList
        test_1 = custom_functions.combine_sirna_gs(self, self.test_seed_list_1)
        
        self.assertTrue(isinstance(test_1, core.SeedList))
        
        self.assertEqual(sorted(test_1.nodeList().ids()), sorted(list(set(self.test_seed_list_1.ids()))))
        
        gene_score = ['ENSG00000254087', 'ENSG00000101336', 'ENSG00000213341', 'ENSG00000240764', 'ENSG00006898999', 'ENSG00000213765', 'ENSG00000263328', 'ENSG00000150261']
        
        gene_score_vals = [50, 50, 50, 50, 15, 30, 40, 50]
        
        self.assertEqual(test_1.getScores(gene_score), gene_score_vals)
    
    @unittest.skip("Working on graph db")
    def test_gene_seed_list_to_protein(self):
        
        test_1 = custom_functions.combine_sirna_gs(self, self.test_seed_list_1)
        
        test_1_prot = custom_functions.gene_seed_list_to_protein(self, test_1)
        
        graph_db = neo4j.GraphDatabaseService()  
        
        data = neo4j.CypherQuery(graph_db,'MATCH (gene:Gene)-[:EXTERNAL_ID]-()-[:MAPPED_TO]-(string_from) WHERE gene.name IN {query_genes} RETURN gene.name, string_from.name').execute(**{'query_genes':test_1.nodeList().ids()}).data
        
        test_prot_map = collections.defaultdict(list)
        
        for i in data:
            test_prot_map[i.values[1]].extend(test_1.getScores([i.values[0]]))
            
        for i in test_prot_map.keys():
            test_prot_map[i] = max(test_prot_map[i])
        
        self.assertEqual(sorted(test_1_prot.nodeList().ids()), sorted(test_prot_map.keys()))
        
        self.assertEqual(test_1_prot.getScores(test_prot_map.keys()), test_prot_map.values())
        
#class Test_get_valid_hits(TestCase):
#    
#    def setUp(self):
#        self.graph_db = neo4j.GraphDatabaseService() 
#        self.sirna_only = neo4j.CypherQuery(self.graph_db,'MATCH (n:LabID) WHERE (n)-[:SIRNA_RUN]->() AND NOT (n)-[:GENE_SCORE_RUN]->() RETURN n').execute_one()["name"]
#        #currently doesn't exist
#        #self.gs_only = neo4j.CypherQuery(graph_db,'MATCH (n:LabID) WHERE (n)-[:GENE_SCORE_RUN]->() AND NOT (n)-[:SIRNA_RUN]->() RETURN n').execute_one()["name"]
#        self.neither = neo4j.CypherQuery(self.graph_db,'MATCH (n:LabID) WHERE NOT (n)-[:GENE_SCORE_RUN]->() AND NOT (n)-[:SIRNA_RUN]->() RETURN n').execute_one()["name"]
#        self.both = neo4j.CypherQuery(self.graph_db,'MATCH (n:LabID) WHERE (n)-[:GENE_SCORE_RUN]->() AND (n)-[:SIRNA_RUN]->() RETURN n').execute_one()["name"]
#        self.thresh_dict = {'siRNA':-2, 'GeneScore':0}
#        #self.conn = sqlite3.connect('/Users/bottomly/tyner_results/hitwalker_sqlite_test/HitWalker2.db')
#        
#    def test_sirna_only(self):
#        gene_score, hit_annot = custom_functions.get_valid_hits ({'siRNA':self.sirna_only, 'GeneScore':self.sirna_only}, self.thresh_dict, self.graph_db)
#        
#        self.assertEqual(set(gene_score.keys()), set(hit_annot.keys()))
#        self.assertTrue(all(map(lambda x: x.has_key('siRNA'), hit_annot.values())))
#        self.assertFalse(any(map(lambda x: x.has_key('GeneScore'), hit_annot.values())))
#        self.assertTrue(max(gene_score.values()) == 1)
#        
#    def test_both(self):
#        gene_score, hit_annot = custom_functions.get_valid_hits ({'siRNA':self.both, 'GeneScore':self.both}, self.thresh_dict, self.graph_db)
#        
#        self.assertEqual(set(gene_score.keys()), set(hit_annot.keys()))
#        
#        max_gene_score = 0
#        
#        for i in hit_annot.values():
#            if i.has_key('GeneScore'):
#                max_gene_score = max([max_gene_score, i['GeneScore']])
#        
#        for i in hit_annot.items():
#            if i[1].has_key('siRNA'):
#                self.assertEqual(gene_score[i[0]], max_gene_score)
#    
#    def test_neither(self):
#        gene_score, hit_annot = custom_functions.get_valid_hits ({'siRNA':self.neither, 'GeneScore':self.neither}, self.thresh_dict, self.graph_db)
#        self.assertEqual(gene_score, {})
#        self.assertEqual(hit_annot, {})
#
#class Test_make_gene_id_dict(TestCase):
#    
#    def setUp(self):
#        self.graph_db = neo4j.GraphDatabaseService()
#        
#    def test_gene_list(self):
#        id_dict = custom_functions.make_gene_id_dict(self.graph_db)
#        
#        self.assertTrue(isinstance(id_dict, dict))
#        self.assertTrue(all(map(lambda x:isinstance(x, list), id_dict.values())))
#
#class Test_gene_id_helpers(TestCase):
#    
#    def setUp(self):
#        self.seed_dict = {'gene1':1, 'gene2':15, 'gene3':10, 'gene4':20}
#        self.hit_annot = {'gene1':{'GeneScore':10}, 'gene2':{'GeneScore':15}, 'gene3':{'GeneScore':1}, 'gene4':{'siRNA':-5, 'GeneScore':5}}
#        self.gene_id_dict = {'gene1':['mg1'], 'gene2':['mg2'], 'gene3':['mg1','mg3'], 'gene4':['mg4']}
#    
#    def test_seed_to_gene_ids(self):
#        hit_dict = custom_functions.seed_to_gene_ids(self.seed_dict, self.gene_id_dict)
#        
#        self.assertEquals(hit_dict, {'mg1':10, 'mg2':15, 'mg3':10, 'mg4':20})
#    
#    def test_seed_annot_to_gene_ids(self):
#        hit_annot_dict = custom_functions.seed_annot_to_gene_ids(self.hit_annot, self.gene_id_dict)
#        
#        self.assertEquals(hit_annot_dict, {'mg1':{'GeneScore':10}, 'mg2':{'GeneScore':15}, 'mg3':{'GeneScore':1}, 'mg4':{'GeneScore':5, 'siRNA':-5}})