from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.support.ui import Select

import time
import re
import scipy.spatial
import numpy as np
import collections

class SelectedNodes(object):
    
    _display_names = []
    _css_classes = []
    _child_classes = []
    
    def __init__(self, sel_obj):
        #display names would come from the text fields
        #css classes would come from the node obj itself
        #child classes would be the subcircles (if any)
        print 'in progress'

class NodeSelectors(object):
    
    _selector=""
    _text_selector=""
    
    def selector(self):
        return self._selector
    
    def text_selector(self):
        return self._text_selector;

class SingleSubjectSelector(NodeSelectors):
    
    def __init__(self):
        self._selector = "g.node > circle.Subject"

class SingleMetaNodeSelector(NodeSelectors):
    
    def __init__(self):
        self._selector = "g.node > circle.MetaNode"
        self._text_selector = "g.labelNode > text"

class HitWalkerInteraction(object):
    
    def __init__(self, driver, live_server_url):
        self.driver = driver
        self.live_server_url = live_server_url
        
    def panel_by_query(self, selection_name):
        self.driver.get('%s%s' % (self.live_server_url, '/HitWalker2'))
        self.driver.find_element_by_css_selector(".select2-choice").click()
        self.driver.find_element_by_css_selector("#select2-drop input.select2-input").send_keys(selection_name)
        self.driver.find_element_by_css_selector(".select2-result-label").click()
        
        element = WebDriverWait(self.driver, 20).until(
            EC.element_to_be_clickable((By.ID, "query"))
        )
        
        element.click()
        
    def panel_by_prioritize(self, selection_name, collect_genes=True):
        self.driver.get('%s%s' % (self.live_server_url, '/HitWalker2'))
        self.driver.find_element_by_css_selector(".select2-choice").click()
        self.driver.find_element_by_css_selector("#select2-drop input.select2-input").send_keys(selection_name)
        self.driver.find_element_by_css_selector(".select2-result-label").click()
        
        p_button = WebDriverWait(self.driver, 20).until(
            EC.element_to_be_clickable((By.ID, "prioritize"))
        )
        
        p_button.click()
        
        modal_button = WebDriverWait(self.driver, 20).until(
            EC.element_to_be_clickable((By.CSS_SELECTOR, "button[data-target='#select_modal']"))
        )
        
        modal_button.click()
        
        if collect_genes:
            WebDriverWait(self.driver, 20).until(
                EC.element_to_be_clickable((By.CSS_SELECTOR, "option[selected='true']"))
            )
            
            targs = self.driver.find_elements_by_css_selector("#var_select > option[selected='true']")
            
            targ_text = map(lambda x: x.text, targs)
            
            seeds = self.driver.find_elements_by_css_selector("#seed_select > option[selected='true']")
            
            seed_text = map(lambda x: x.text, seeds)
            
            gene_text = {'targs':targ_text, 'seeds':seed_text}
            
        else:
            
            gene_text = None
        
        submit_button = self.driver.find_element_by_css_selector("button[type='submit']")
        submit_button.click()
        
        WebDriverWait(self.driver, 100).until(
                 EC.element_to_be_clickable((By.ID,"panel_1"))
        )
        
        return gene_text
    
    def click_panel(self, panel_num):
        use_panel = self.to_panel(panel_num)
        
        to_panel = webdriver.ActionChains(self.driver).move_to_element_with_offset(use_panel, 1, 100)
        
        to_panel.click().perform()
    
    def click_context_button(self, panel_num, button_num):
        use_panel = self.to_panel(panel_num)
        
        to_panel = webdriver.ActionChains(self.driver).move_to_element_with_offset(use_panel, 1, 100)
        
        to_panel.context_click().perform()
        
        use_button = WebDriverWait(self.driver, 20).until(
            EC.element_to_be_clickable((By.CSS_SELECTOR, "div.popover-content > div > div > div > div:nth-of-type("+str(button_num)+") > button"))
        )
        
        use_button.click()
    
    def delete_panel(self, panel_num):
        self.click_context_button(panel_num, 1)   
    
    def get_metanode(self, panel_num, node_sel):
        
        cur_panel = self.to_panel(panel_num)
        
        return WebDriverWait(self.driver, 20).until(
                EC.element_to_be_clickable((By.CSS_SELECTOR,node_sel.selector()))
            )
    
    def get_panel_nodes(self, panel_num):
        cur_panel = self.to_panel(panel_num)
        
        all_nodes = cur_panel.find_elements_by_css_selector("g.node")
        
        return all_nodes
    
    def check_if_metanode(self, node_sel):
        try:
            
            node_sel.find_element_by_css_selector("circle.MetaNode")
            
            return True
            
        except NoSuchElementException:
            
            return False
            
        
    def get_panel_metanodes(self, panel_num):
        
        panel_nodes = self.get_panel_nodes(panel_num)
        
        return filter(self.check_if_metanode, panel_nodes)
    
    def to_panel(self, panel_num):
    
        return WebDriverWait(self.driver, 20).until(
                EC.element_to_be_clickable((By.ID,"panel_"+str(panel_num)))
            )
    
    #to find all panels
                #panels = self.driver.find_elements_by_css_selector("g.g1[id^='panel_']")
                #
                #panel_nums = map(lambda x:x.get_attribute("id").replace("panel_",""),panels)
    
    def select_context_node(self, panel, node_sel):
        
        cur_node = panel.find_element_by_css_selector(node_sel.selector())
        
        node_move = webdriver.ActionChains(self.driver).move_to_element(cur_node)
        
        node_move.click()
        
        node_move.context_click(cur_node).perform()

    def get_metanode_attrs(self, panel, node_sel):
        
        cur_node = panel.find_element_by_css_selector(node_sel.selector())
        
        text_els = panel.find_elements_by_css_selector(node_sel.text_selector())
        
        #should be the second one
        
        text = text_els[1].text
        
        #for a metanode text should be of the form:
        #type (count)
        
        split_text = re.findall("(\w+)\\s+\((\d+)\)", text)
        
        return {'type':split_text[0][0], 'count':int(split_text[0][1])}
    
    def count_nodes(self, panel, node_type):
        count_panel = self.to_panel(panel)
        
        nodes = count_panel.find_elements_by_css_selector("g > circle."+node_type)
        
        return len(nodes)
    
    def count_metanode_children(self, panel, node_sel, child_type="Subject"):
        
        if isinstance(node_sel, NodeSelectors):
        
            count_panel = self.to_panel(panel)
            
            children = count_panel.find_elements_by_css_selector(node_sel.selector() + " ~ g > circle")
            
        else:
            #assumes its a WebElement...
            #just want the outer g's...
            children = node_sel.find_elements_by_css_selector("g > circle." + child_type)
        
        return len(children)
        
    
    def get_node_rels(self, panel_num):
        
        node_list = self.get_panel_nodes(panel_num)
        
        #should settle enough after a couple of seconds
        
        self.click_panel(panel_num)
        
        #time.sleep(2)
        
        trans_vals = map(lambda x: map(float, re.split("[,\)\(]", x.get_attribute("transform"))[1:3]), node_list)
        
        node_labels = map(lambda x: x.find_element_by_css_selector("circle").get_attribute("name"), node_list)
        
        kd_tree = scipy.spatial.KDTree(np.array(trans_vals))
        
        links = self.to_panel(panel_num)
            
        link_els = links.find_elements_by_css_selector("path.link")
        
        link_pos = map(lambda x: map(float, re.split("[M,L]", x.get_attribute("d"))[1:]) , link_els)
        
        #in this case all the links should point to the same position, so we are just looking at which metanode should be assigned the relationship
        
        link_type = map(lambda x: x.get_attribute("class").replace("link ", ""), link_els)
        
        link_dict = collections.defaultdict(lambda:collections.defaultdict(set))
        
        for i_ind, i in enumerate(link_pos):
            left_match = kd_tree.query(np.array(i[:2]))
            right_match = kd_tree.query(np.array(i[2:]))
            link_dict[node_labels[left_match[1]]][node_labels[right_match[1]]].add(link_type[i_ind])
            link_dict[node_labels[right_match[1]]][node_labels[left_match[1]]].add(link_type[i_ind])
        
        return link_dict
    
    def click_by_text(self, selector, text):
        
        WebDriverWait(self.driver, 20).until(
            EC.presence_of_element_located((By.CSS_SELECTOR,selector))
            )
        
        clickers = self.driver.find_elements_by_css_selector(selector)
        
        pathway_clicker = filter(lambda x: x.text == text, clickers)
        
        assert len(pathway_clicker) == 1
        
        pathway_clicker[0].click()
    
    def use_select2_box(self, sel_text):
        self.driver.find_element_by_css_selector(".select2-choice").click()
        self.driver.find_element_by_css_selector("input.select2-input").send_keys(sel_text)
        
        input_highlight = WebDriverWait(self.driver, 20).until(
            EC.element_to_be_clickable((By.CSS_SELECTOR,".select2-result-label"))
            )
        
        labels = self.driver.find_elements_by_css_selector(".select2-result-label")
        
        labels[0].click()
    
    def add_pathway(self, panel_num, gene_names, path_num=None):
        
        if panel_num != None:
        
            self.click_context_button(panel_num, 2)
            
            self.click_by_text("ul.dropdown-menu li a", "Pathway")
        
        for i in gene_names:
            
            input_text = WebDriverWait(self.driver, 20).until(
                EC.element_to_be_clickable((By.CSS_SELECTOR,"input.select2-input"))
                )
            
            input_text.send_keys(i)
            
            input_highlight = WebDriverWait(self.driver, 20).until(
                EC.element_to_be_clickable((By.CSS_SELECTOR,".select2-result-label"))
            )
        
            labels = self.driver.find_elements_by_css_selector(".select2-result-label")
            
            labels[0].click()
        
        time.sleep(2)
        
        if path_num == None:
            path_num = 0
        
        print path_num
        
        select = Select(self.driver.find_element_by_css_selector("#pathway_select"))
        
        select.deselect_all()
        
        select.select_by_index(path_num)
        
        pathway_sel = select.all_selected_options[0].text
        
        buttons = self.driver.find_elements_by_css_selector("button")
        
        #print buttons
        
        ok_button = filter(lambda x: x.text == "OK", buttons)
        
        #print ok_button
        
        assert len(ok_button) == 1
        
        ok_button[0].click()
        
        return pathway_sel
    
    
    def add_subject(self, panel_num, subj_name):
        
        self.click_context_button(panel_num, 2)
        
        self.click_by_text("ul.dropdown-menu li a", "Subject")
        
        self.use_select2_box(subj_name)
        
        self.click_by_text("button", "OK")
        
    
    def add_gene(self, panel_num, gene_name):
    
        self.click_context_button(panel_num, 2)
        
        self.driver.find_element_by_css_selector("ul li:first-child a").click()
        
        #select a gene that probably has hits, say TP53
        
        self.driver.find_element_by_css_selector(".select2-choice").click()
        self.driver.find_element_by_css_selector("input.select2-input").send_keys(gene_name)
        
        input_highlight = WebDriverWait(self.driver, 20).until(
            EC.element_to_be_clickable((By.CSS_SELECTOR,".select2-result-label"))
            )
        
        labels = self.driver.find_elements_by_css_selector(".select2-result-label")
        
        labels[0].click()
        
        buttons = self.driver.find_elements_by_css_selector("button")
        
        #print buttons
        
        ok_button = filter(lambda x: x.text == "OK", buttons)
        
        #print ok_button
        
        assert len(ok_button) == 1
        
        ok_button[0].click()
        
        #time.sleep(5)
        
        #self.driver.find_element_by_xpath("//button[.='OK']").click()