# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.
Vagrant.configure(2) do |config|
  # The most common configuration options are documented and commented below.
  # For a complete reference, please see the online documentation at
  # https://docs.vagrantup.com.

  # Every Vagrant development environment requires a box. You can search for
  # boxes at https://atlas.hashicorp.com/search.
  config.vm.box = "ubuntu/trusty64"
  # Disable automatic box update checking. If you disable this, then
  # boxes will only be checked for updates when the user runs
  # `vagrant box outdated`. This is not recommended.
   config.vm.box_check_update = false

  # Create a forwarded port mapping which allows access to a specific port
  # within the machine from a port on the host machine. In the example below,
  # accessing "localhost:8080" will access port 80 on the guest machine.
  #config.vm.network "forwarded_port", guest: 80, host: 7080

  # Create a private network, which allows host-only access to the machine
  # using a specific IP.
   config.vm.network "private_network", type: "dhcp"

  # Create a public network, which generally matched to bridged network.
  # Bridged networks make the machine appear as another physical device on
  # your network.
  # config.vm.network "public_network"

  # Share an additional folder to the guest VM. The first argument is
  # the path on the host to the actual folder. The second argument is
  # the path on the guest to mount the folder. And the optional third
  # argument is a set of non-required options.
  #config.vm.synced_folder "data/", "/vagrant_data",owner: "neo4j", group: "neo4j"
  
  config.ssh.insert_key = false
  
  # Provider-specific configuration so you can fine-tune various
  # backing providers for Vagrant. These expose provider-specific options.
  # Example for VirtualBox:
  #
  config.vm.provider "virtualbox" do |vb|
  #   # Display the VirtualBox GUI when booting the machine
  #vb.gui = true
  vb.name="ccle"
  #   # Customize the amount of memory on the VM:
  vb.memory = "16384"
  vb.cpus = 4   
end
  #
  # View the documentation for the provider you are using for more
  # information on available options.

  # Define a Vagrant Push strategy for pushing to Atlas. Other push strategies
  # such as FTP and Heroku are also available. See the documentation at
  # https://docs.vagrantup.com/v2/push/atlas.html for more information.
  # config.push.define "atlas" do |push|
  #   push.app = "YOUR_ATLAS_USERNAME/YOUR_APPLICATION_NAME"
  # end

  # Enable provisioning with a shell script. Additional provisioners such as
  # Puppet, Chef, Ansible, Salt, and Docker are also available. Please see the
  # documentation for more information about their specific syntax and use.
  
  config.vm.provision "basic", type:"shell", inline: <<-SHELL
      sudo apt-get update
     
      sudo apt-get install -y nginx python-pip python-dev build-essential  python-numpy python-scipy python-cairosvg python-lxml wget openjdk-7-jdk 
      sudo pip install --upgrade pip 

      sudo pip install Django==1.5.1
      sudo pip install py2neo==1.6.4
      
      sudo pip install tinycss
      sudo pip install colour
      
      sudo apt-get install -y gunicorn python-cssselect python-cssutils
     
      sudo pip install eventlet==0.17.1
     
    sudo bash -c "echo 'vagrant   soft    nofile  40000' >> /etc/security/limits.conf"
    sudo bash -c "echo 'vagrant   hard    nofile  40000' >> /etc/security/limits.conf"
      
    sudo bash -c "echo 'session    required   pam_limits.so' >> /etc/pam.d/su"
      
      sudo mkdir -p /var/www/hitwalker2_inst
      
      sudo bash -c "echo 'deb http://ftp.osuosl.org/pub/cran/bin/linux/ubuntu trusty/' >> /etc/apt/sources.list"
      
      sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
      
      sudo apt-get update -y
      
      sudo apt-get install -y r-base-dev libxml2-dev libcurl4-openssl-dev
      
      sudo Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite(c("igraph", "reshape2", "Biobase", "rjson", "affy","hgu133plus2.db", "tm", "devtools", "sva", "genefilter"))'
      
      sudo Rscript -e 'library(devtools)' -e 'install_github("mlbernauer/Entrez")'
      
      #set up upstart for unicorn
    
echo '
description "HitWalker2" 
start on (filesystem)
stop on runlevel [016]
respawn
setuid vagrant
setgid vagrant
chdir /home/vagrant/HitWalker2
    
exec gunicorn -k 'eventlet' HitWalker2.wsgi:application
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
  
  sudo R CMD INSTALL hwhelper
  
  sudo chown -R vagrant:vagrant /var/www/
  
  python /home/vagrant/HitWalker2/manage.py collectstatic --noinput
  
  sudo chown -R vagrant:vagrant /var/www/
  
   SHELL
   
  config.vm.provision "neo4j", type:"shell", inline: <<-SHELL
   
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
  
    sudo rm -rf /opt/neo4j-community-2.1.8/data
    sudo cp -r /vagrant/hitwalker2_base_data /opt/neo4j-community-2.1.8/data
    
    sudo chown -R vagrant:vagrant /opt/neo4j-community-2.1.8/
    
    sudo service neo4j-service start
      
      
  SHELL
  
  config.vm.provision "sweave", type:"shell", inline: <<-SHELL
  
  sudo apt-get install -y texlive-latex-recommended
  sudo apt-get install -y texinfo
  sudo apt-get install -y texlive-latex-extra
  
  SHELL
  
end
