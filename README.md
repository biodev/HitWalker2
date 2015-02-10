# HitWalker2
=========

Installation
----------

Requires: 

Python 2.7.2
Django 1.5.1
scipy
numpy
py2neo 1.6.4
neo4j-community 2.1.6

Suggested (for exporting PDFs):
tinycss
cairosvg
lxml
cssselect

Accessing the neo4j database
----------

ssh -N -L 7474:localhost:7474 username@10.132.57.20


Running the development instance
----------

From HitWalker2 directory: 

./change_hw2_instance.sh lls
python manage.py runserver
