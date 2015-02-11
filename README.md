# HitWalker2

Installation
----------

**Requires:** 
* neo4j-community 2.1.6

* Python 2.7.2

  * Django 1.5.1

  * scipy

  * numpy

  * py2neo 1.6.4
  
  * gunicorn
  
  * pandas

**Suggested Python packages (for exporting PDFs):**

* tinycss

* cairosvg

* lxml

* cssselect

For Linux installation instructions see this [page](INSTALL.md)

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
