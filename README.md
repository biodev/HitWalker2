# HitWalker2

Installation
----------

**Requires:** 

* Python 2.7.2

* Django 1.5.1

* scipy

* numpy

* pandas

* gunicorn

* py2neo 1.6.4

* neo4j-community 2.1.6

**Suggested (for exporting PDFs):**

* tinycss

* cairosvg

* lxml

* cssselect

Accessing the neo4j database
----------

ssh -N -L 7474:localhost:7474 username@noble.ohsu.edu


Running the development instance
----------

From HitWalker2 directory, create the appropriate config.py and custom_functions.py files: 

./change_hw2_instance.sh lls

Make sure to specify the graph structure file.  By default it is located at:

/var/www/hitwalker_2_inst/graph_struct.json

If desired, this can be changed in the config file by modifying the line:

graph_struct_file = "/path/to/file"

The current file for the lls database is at HitWalker2/lls_graph_struct.json

Finally run:

python manage.py runserver


**Additional notes**
* On the Mac, if Neo4j returns an error like "Unable to find any JVMs matching version 1.7," check which JVMs you have by running:

/usr/libexec/java_home -V

If you do not see Java SE 7, install this:
http://www.oracle.com/technetwork/java/javase/downloads/jdk7-downloads-1880260.html

* If Django reports "unable to open database file," check the permissions on the directory (and parent directories) where your database lives. By default this is /var/www/hitwalker_2_inst/  

SQLite needs to be able to write to this directory. 

* If Django reports "DatabaseError: no such table: django_site", create necessary tables by running: 

python manage.py syncdb




