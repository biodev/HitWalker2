from unittest import TestCase, main
from py2neo import neo4j
import config
import custom_functions
import sqlite3


class Test_get_valid_hits(TestCase):
    
    def setUp(self):
        self.graph_db = neo4j.GraphDatabaseService() 
        self.sirna_only = neo4j.CypherQuery(self.graph_db,'MATCH (n:LabID) WHERE (n)-[:SIRNA_RUN]->() AND NOT (n)-[:GENE_SCORE_RUN]->() RETURN n').execute_one()["name"]
        #currently doesn't exist
        #self.gs_only = neo4j.CypherQuery(graph_db,'MATCH (n:LabID) WHERE (n)-[:GENE_SCORE_RUN]->() AND NOT (n)-[:SIRNA_RUN]->() RETURN n').execute_one()["name"]
        self.neither = neo4j.CypherQuery(self.graph_db,'MATCH (n:LabID) WHERE NOT (n)-[:GENE_SCORE_RUN]->() AND NOT (n)-[:SIRNA_RUN]->() RETURN n').execute_one()["name"]
        self.both = neo4j.CypherQuery(self.graph_db,'MATCH (n:LabID) WHERE (n)-[:GENE_SCORE_RUN]->() AND (n)-[:SIRNA_RUN]->() RETURN n').execute_one()["name"]
        self.thresh_dict = {'siRNA':-2, 'GeneScore':0}
        #self.conn = sqlite3.connect('/Users/bottomly/tyner_results/hitwalker_sqlite_test/HitWalker2.db')
        
    def test_sirna_only(self):
        gene_score, hit_annot = custom_functions.get_valid_hits ({'siRNA':self.sirna_only, 'GeneScore':self.sirna_only}, self.thresh_dict, self.graph_db)
        
        self.assertEqual(set(gene_score.keys()), set(hit_annot.keys()))
        self.assertTrue(all(map(lambda x: x.has_key('siRNA'), hit_annot.values())))
        self.assertFalse(any(map(lambda x: x.has_key('GeneScore'), hit_annot.values())))
        self.assertTrue(max(gene_score.values()) == 1)
        
    def test_both(self):
        gene_score, hit_annot = custom_functions.get_valid_hits ({'siRNA':self.both, 'GeneScore':self.both}, self.thresh_dict, self.graph_db)
        
        self.assertEqual(set(gene_score.keys()), set(hit_annot.keys()))
        
        max_gene_score = 0
        
        for i in hit_annot.values():
            if i.has_key('GeneScore'):
                max_gene_score = max([max_gene_score, i['GeneScore']])
        
        for i in hit_annot.items():
            if i[1].has_key('siRNA'):
                self.assertEqual(gene_score[i[0]], max_gene_score)
    
    def test_neither(self):
        gene_score, hit_annot = custom_functions.get_valid_hits ({'siRNA':self.neither, 'GeneScore':self.neither}, self.thresh_dict, self.graph_db)
        self.assertEqual(gene_score, {})
        self.assertEqual(hit_annot, {})

class Test_make_gene_id_dict(TestCase):
    
    def setUp(self):
        self.graph_db = neo4j.GraphDatabaseService()
        
    def test_gene_list(self):
        id_dict = custom_functions.make_gene_id_dict(self.graph_db)
        
        self.assertTrue(isinstance(id_dict, dict))
        self.assertTrue(all(map(lambda x:isinstance(x, list), id_dict.values())))

class Test_gene_id_helpers(TestCase):
    
    def setUp(self):
        self.seed_dict = {'gene1':1, 'gene2':15, 'gene3':10, 'gene4':20}
        self.hit_annot = {'gene1':{'GeneScore':10}, 'gene2':{'GeneScore':15}, 'gene3':{'GeneScore':1}, 'gene4':{'siRNA':-5, 'GeneScore':5}}
        self.gene_id_dict = {'gene1':['mg1'], 'gene2':['mg2'], 'gene3':['mg1','mg3'], 'gene4':['mg4']}
    
    def test_seed_to_gene_ids(self):
        hit_dict = custom_functions.seed_to_gene_ids(self.seed_dict, self.gene_id_dict)
        
        self.assertEquals(hit_dict, {'mg1':10, 'mg2':15, 'mg3':10, 'mg4':20})
    
    def test_seed_annot_to_gene_ids(self):
        hit_annot_dict = custom_functions.seed_annot_to_gene_ids(self.hit_annot, self.gene_id_dict)
        
        self.assertEquals(hit_annot_dict, {'mg1':{'GeneScore':10}, 'mg2':{'GeneScore':15}, 'mg3':{'GeneScore':1}, 'mg4':{'GeneScore':5, 'siRNA':-5}})

if __name__ == "__main__":
    main()