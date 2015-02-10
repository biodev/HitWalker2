from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.common.exceptions import NoSuchElementException
from selenium.common.exceptions import NoAlertPresentException
import unittest, time, re
import sys
import os
sys.path.append(os.path.abspath("../network/"))
import config
import random

class seleniumTests (unittest.TestCase):
    def setUp(self):
        self.driver = webdriver.Firefox()
        self.driver.implicitly_wait(1)
        self.driver.get("http://127.0.0.1:8000/HitWalker2")
        elem = self.driver.find_element_by_id("id_username")
        elem.send_keys("selenium")
        elem = self.driver.find_element_by_id("id_password")
        elem.send_keys("test")
        
        self.driver.find_element_by_css_selector("input[type=submit]").click()
        
        #webdriver.ActionChains(self.driver).move_to_element(button_elem).click(button_elem).perform()
    
    def is_element_present(self, select):
        try: self.driver.find_element_by_css_selector(select)
        except NoSuchElementException, e: return False
        return True
    
    def enter_vals_modals(self, inp_vals_func, should_fail=False):
        self.driver.get("http://127.0.0.1:8000/HitWalker2")
        for key, value in config.adjust_fields.items():
            self.driver.find_element_by_css_selector("button.dropdown-toggle").click()
            self.driver.find_element_by_css_selector("a[data-target$="+key+"]").click()
            if value["type"] in set(["standard", "grouped"]):
                print value
                for field_name, field_val in value["fields"].items():
                    cur_elem = self.driver.find_element_by_id(field_name)
                    if field_val["type"] == "numeric":
                        inp_vals = inp_vals_func(field_val)
                        cur_elem.clear()
                        time.sleep(5)
                        cur_elem.send_keys(inp_vals)
                    else:
                        raise Exception("Tests not currently implemented for type character")
            else:
                raise Exception("Unimplemented type: " + value["type"])
            #click OK
            #note that we are not checking the x, as there is some issue with selenium here...
            self.driver.find_element_by_css_selector("#"+key+" > div.modal-dialog.modal-lg > div.modal-content > div.modal-footer > button").click()
            
            time.sleep(5)
            #body class should be "" as opposed to "modal-open"
            if should_fail == False:
                self.assertFalse(self.is_element_present("body.modal-open"))
            else:
                #modal is still open, go back to main page to clear
                self.assertTrue(self.is_element_present("body.modal-open"))
                self.driver.get("http://127.0.0.1:8000/HitWalker2")
    @unittest.skip("Not to this part yet...")
    def test_submit(self):
        self.driver.find_element_by_css_selector(".select2-choice").click()
        self.driver.find_element_by_css_selector("#select2-drop input.select2-input").send_keys("12-00145")
        self.driver.find_element_by_css_selector(".select2-result-label").click()
        time.sleep(1)
        self.driver.find_element_by_css_selector("#prioritize").click()
        
        #also submit the table data
        self.driver.find_element_by_css_selector("button[data-target='#select_modal']").click()
        self.driver.find_element_by_css_selector("button[type='submit']").click()
        time.sleep(5)
    
    @unittest.skip("Need to work on me")
    def test_index_modals(self):
        
        #self.enter_vals_modals(lambda field_val: str(field_val["default"]), False)
        #test other random values within defined ranges do not result in errors
        #self.enter_vals_modals(lambda field_val: str(round(random.uniform(field_val["range"][0], field_val["range"][1]), 2)), False)
        #test out of range on the low side
        self.enter_vals_modals(lambda field_val: str(round(random.uniform(field_val["range"][0]-50, field_val["range"][0]), 2)), True)
        #test out of range on the high side
        #self.enter_vals_modals(lambda field_val: str(round(random.uniform(field_val["range"][1], field_val["range"][1]+50), 2)), True)
        #enter in invalid character values
        #self.enter_vals_modals(lambda field_val: "test", True)
    
    @unittest.skip("Not to this part yet...")
    def test_sample_selection(self):
        #for the LLS version can use
        #12-00145 for single sample
        #12-00156 for multi sample
        #something else for invalid sample...
        print 'hello'
    
    @unittest.skip("Not to this part yet...")
    def test_save_and_load(self):
        print 'hello'
    
    @unittest.skip("Not to this part yet...")
    def test_reset_to_defaults(self):
        print 'hello'
    
    def test_gene_addition(self):
    
        self.driver.get("http://127.0.0.1:8000/HitWalker2")
        self.driver.find_element_by_css_selector(".select2-choice").click()
        self.driver.find_element_by_css_selector("#select2-drop input.select2-input").send_keys("07-00112")
        self.driver.find_element_by_css_selector(".select2-result-label").click()
        time.sleep(1)
        self.driver.find_element_by_css_selector("#query").click()
        
        #right click on a panel
        
        use_panel = self.driver.find_element_by_css_selector("rect.BorderRect")
        
        webdriver.ActionChains(self.driver).move_to_element(use_panel).context_click(use_panel).perform()
        
        self.driver.find_element_by_css_selector("div.btn-group-vertical div.btn-group:nth-child(2) button").click()
        
        self.driver.find_element_by_css_selector("ul li:first-child a").click()
        
        #enter the gene name
        self.driver.find_element_by_css_selector(".select2-choice").click()
        #self.driver.find_element_by_css_selector("input.select2-input").send_keys("CLSTN2")
        self.driver.find_element_by_css_selector("input.select2-input").send_keys("ROR1")
        self.driver.find_element_by_css_selector(".select2-result-label").click()
        time.sleep(1)
        
        self.driver.find_element_by_xpath("//button[.='OK']").click()
        
        
        time.sleep(10)
    
    @unittest.skip("Not to this part yet...")
    def test_query_path(self):
        
        self.driver.get("http://127.0.0.1:8000/HitWalker2")
        self.driver.find_element_by_css_selector(".select2-choice").click()
        self.driver.find_element_by_css_selector("#select2-drop input.select2-input").send_keys("12-00145")
        self.driver.find_element_by_css_selector(".select2-result-label").click()
        time.sleep(1)
        self.driver.find_element_by_css_selector("#query").click()
        
        #a couple of tests to make sure it basically works
        
        #check click events
        
        use_panel = self.driver.find_element_by_css_selector(".Sample~circle")
        
        
        #to go into the pathway view
        webdriver.ActionChains(self.driver).move_to_element(use_panel).click(use_panel).context_click(use_panel).perform()
        
        self.driver.find_element_by_xpath("(//button[@type='button'])[2]").click()
        self.driver.find_element_by_css_selector(".select2-choices").click()
        self.driver.find_element_by_css_selector("input.select2-input").send_keys("IP2")
        self.driver.find_element_by_css_selector(".select2-result-label").click()
        
        path_select = Select(self.driver.find_element_by_id("pathway_select"))
        path_select.deselect_all()
        path_select.select_by_visible_text("Innate Immune System")
        
        self.driver.find_element_by_css_selector("div.row > div.btn-group > button.btn.btn-default").click()
        
        #bug here, need to ensure that request.session['hitwalker_score'] is initialized or the edge finding code is more robust
        
        time.sleep(10)
        
    def tearDown(self):
       self.driver.quit()

#class Test(unittest.TestCase):
#    def setUp(self):
#        self.driver = webdriver.Firefox()
#        self.driver.implicitly_wait(30)
#        self.base_url = "http://127.0.0.1:8000/"
#        self.verificationErrors = []
#        self.accept_next_alert = True
#    
#    def test_(self):
#        driver = self.driver
#        driver.get(self.base_url + "/HitWalker2/")
#        # ERROR: Caught exception [ERROR: Unsupported command [clickAt | css=.select2-choice | ]]
#        driver.find_element_by_css_selector("#select2-drop input.select2-input").send_keys("12-00145")
#        # ERROR: Caught exception [ERROR: Unsupported command [clickAt | css=.select2-result-label:contains('12-00145') | ]]
#        driver.find_element_by_id("query").click()
#    

#    
#    def is_alert_present(self):
#        try: self.driver.switch_to_alert()
#        except NoAlertPresentException, e: return False
#        return True
#    
#    def close_alert_and_get_its_text(self):
#        try:
#            alert = self.driver.switch_to_alert()
#            alert_text = alert.text
#            if self.accept_next_alert:
#                alert.accept()
#            else:
#                alert.dismiss()
#            return alert_text
#        finally: self.accept_next_alert = True
#    
    #def tearDown(self):
    #   self.driver.quit()
#        self.assertEqual([], self.verificationErrors)

if __name__ == "__main__":
    unittest.main()

