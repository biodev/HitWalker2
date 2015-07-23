

from /Users/bottomly/Desktop/hitwalker2_paper/base_vagrant

Adapted from the blog: https://scotch.io/tutorials/how-to-create-a-vagrant-base-box-from-an-existing-one

```
    vagrant init ubuntu/trusty64
    vagrant up
    vagrant ssh
    
    Changed memory to 8GB
```

```

    sudo apt-get update -y
     
    sudo apt-get install -y nginx python-pip python-dev build-essential  python-numpy python-scipy python-cairosvg python-lxml wget openjdk-7-jdk 
    sudo pip install --upgrade pip 

    sudo pip install Django==1.5.1
    sudo pip install py2neo==1.6.4
    
    sudo pip install tinycss
    sudo pip install colour
    sudo pip install selenium
    sudo pip install pyRserve
    
    sudo apt-get install -y gunicorn python-cssselect python-cssutils python-bs4
   
    sudo pip install eventlet==0.17.1
   
    sudo bash -c "echo 'vagrant   soft    nofile  40000' >> /etc/security/limits.conf"
    sudo bash -c "echo 'vagrant   hard    nofile  40000' >> /etc/security/limits.conf"
      
    sudo bash -c "echo 'session    required   pam_limits.so' >> /etc/pam.d/su"
    
    sudo mkdir -p /var/www/hitwalker2_inst
    
    sudo bash -c "echo 'deb http://ftp.osuosl.org/pub/cran/bin/linux/ubuntu trusty/' >> /etc/apt/sources.list"
    
    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
    
    sudo apt-get update -y
    
    sudo apt-get install -y r-base-dev libxml2-dev libcurl4-openssl-dev
    
    sudo Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite(c("Rserve", "igraph", "reshape2", "Biobase", "rjson", "affy","hgu133plus2.db", "tm", "devtools", "sva", "genefilter"))'
    
    sudo Rscript -e 'library(devtools)' -e 'install_github("nicolewhite/RNeo4j")'
    
    sudo cp /vagrant/neo4j-community-2.1.8-unix.tar.gz /opt/
  
  cd /opt/
  
  sudo tar -xzvf /opt/neo4j-community-2.1.8-unix.tar.gz
      
      
    echo "
113c113 
<   read -p 'Press any key to continue'
---
> #  read -p 'Press any key to continue'
123c123
< HEADLESS=false
---
> HEADLESS=true
124c124
< DEFAULT_USER='neo4j'
---
> DEFAULT_USER='vagrant'
" > /vagrant/temp.diff

  cd /opt/neo4j-community-2.1.8/bin/

  sudo patch neo4j-installer /vagrant/temp.diff
  
  sudo rm /vagrant/temp.diff
  
    sudo ./neo4j-installer install
    
    sudo ln -s /opt/neo4j-community-2.1.8/bin/neo4j-shell /usr/local/bin/neo4j-shell
  
    sudo service neo4j-service stop
    
    echo "
25,26c25,26
< #wrapper.java.initmemory=512
< #wrapper.java.maxmemory=512
---
> wrapper.java.initmemory=512
> wrapper.java.maxmemory=10240
" > /vagrant/temp_2.diff

  sudo patch /opt/neo4j-community-2.1.8/conf/neo4j-wrapper.conf /vagrant/temp_2.diff
  
  sudo rm /vagrant/temp_2.diff
  
  
  
echo "244a245
>       initctl emit neo4j-started" > /vagrant/temp_3.diff
  
    sudo patch /opt/neo4j-community-2.1.8/conf/neo4j /vagrant/temp_3.diff
  
    sudo rm /vagrant/temp_3.diff
  
    sudo chown -R vagrant:vagrant /opt/neo4j-community-2.1.8/
    
    sudo service neo4j-service start
    
    
    sudo apt-get clean
    
    cat /dev/null > ~/.bash_history && history -c 
  
    vagrant package --output HitWalker2_base.box







    ###set up upstart for unicorn -- Part of HW2 vagrant file
    
 sudo rm -rf /opt/neo4j-community-2.1.8/data
    sudo cp -r /vagrant/hitwalker2_base_data /opt/neo4j-community-2.1.8/data
    
    echo '
description "HitWalker2"
start on neo4j-started
stop on runlevel [016]
respawn
setuid vagrant
setgid vagrant

pre-start script

python /home/vagrant/HitWalker2/network/warm_up.py

end script

chdir /home/vagrant/HitWalker2/

exec gunicorn -k eventlet HitWalker2.wsgi:application

' > HitWalker2.conf
    
  sudo cp HitWalker2.conf /etc/init/
  
  rm HitWalker2.conf
  
  #set up nginx
  
  #this is only for a non-ssl version
  sudo cp /vagrant/HitWalker2/hw2-nginx /etc/nginx/sites-available/
  sudo ln -sf /etc/nginx/sites-available/hw2-nginx /etc/nginx/sites-enabled/default

  cp -r /vagrant/HitWalker2 /home/vagrant/
  
  sudo chown -R vagrant:vagrant /home/vagrant/HitWalker2
  
  cd /vagrant/HitWalker2/populate
  
  sudo ./roxygen_build.sh install
  
  sudo chown -R vagrant:vagrant /var/www/
  
  python /home/vagrant/HitWalker2/manage.py collectstatic --noinput
  
  sudo chown -R vagrant:vagrant /var/www/
  
  
    

  
  
  


```