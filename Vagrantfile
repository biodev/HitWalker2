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
  config.vm.box = "bottomly/HitWalker2_base"
   config.vm.define :hw2default do |t|
        end

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
  config.ssh.forward_x11 = true
  # Provider-specific configuration so you can fine-tune various
  # backing providers for Vagrant. These expose provider-specific options.
  # Example for VirtualBox:
  #
  config.vm.provider "virtualbox" do |vb|
  #   # Display the VirtualBox GUI when booting the machine
  #vb.gui = true
  vb.name="hw2default"
  #   # Customize the amount of memory on the VM:
  vb.memory = "10240"
  vb.cpus = 2   
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
  

  config.vm.provision "HitWalker2", type:"shell", inline: <<-SHELL
  
  #this is only for a non-ssl version
  sudo cp /vagrant/HitWalker2/hw2-nginx /etc/nginx/sites-available/
  sudo ln -sf /etc/nginx/sites-available/hw2-nginx /etc/nginx/sites-enabled/default

  sudo Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite(c("VariantAnnotation", "SeqVarTools", "GenomicRanges", "jsonlite"))' 
 
  cp -r /vagrant/HitWalker2 /home/vagrant/
  
  sudo chown -R vagrant:vagrant /home/vagrant/HitWalker2
  
  cd /vagrant/HitWalker2/populate
  
  sudo ./roxygen_build.sh install
  
  sudo chown -R vagrant:vagrant /var/www/
  
  python /home/vagrant/HitWalker2/manage.py collectstatic --noinput
  
  sudo chown -R vagrant:vagrant /var/www/
  
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
  
  ##now neo4j
 
  sudo service neo4j-service stop
 
  sudo rm -rf /opt/neo4j-community-2.1.8/data
  sudo tar -xvzf /vagrant/hitwalker2_base_data*
  sudo mv hitwalker2_base_data /opt/neo4j-community-2.1.8/data
  
  sudo chown -R vagrant:vagrant /opt/neo4j-community-2.1.8/
  
  sudo service neo4j-service start
    
      
  SHELL
  
  config.vm.provision "markdown", type:"shell", inline: <<-SHELL
  
  sudo Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("rmarkdown")'
  
  wget https://github.com/jgm/pandoc/releases/download/1.15/pandoc-1.15-1-amd64.deb
  sudo dpkg -i pandoc-1.15-1-amd64.deb 
  
  rm pandoc-1.15-1-amd64.deb

  SHELL

  config.vm.provision "pdf", type:"shell", inline: <<-SHELL

  sudo apt-get install -y texlive-latex-recommended
  sudo apt-get install -y texinfo
  sudo apt-get install -y texlive-latex-extra

  SHELL
  
  config.vm.provision "patch", type:"shell", inline: <<-SHELL
  
  cd /home/vagrant
  
  rm -rf HitWalker2
  
  cp -r /vagrant/HitWalker2 .
  
  cd HitWalker2/populate
  
  sudo ./roxygen_build.sh install
  
  cd ..
  
  python manage.py collectstatic --noinput

  if [ -f /vagrant/hw2_config.RData ]
  then 
 
  	Rscript -e 'library(hwhelper)' -e 'load("/vagrant/hw2_config.RData")' -e 'configure(hw2.conf)'
  fi

  cd /vagrant
  
  sudo restart HitWalker2
  
  sudo nginx -s reload
  
  SHELL
  
end
