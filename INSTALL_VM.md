#get VirtualBox



Following the instructions:

https://help.ubuntu.com/community/Ubuntu_as_Guest_OS

Download the following:

http://www.ubuntu.com/download/server/thank-you?country=US&version=14.04.2&architecture=amd64

Open VirtualBox and create a new one e.g.:

Name: HitWalker2_test
Type: Linux
Version Ubuntu (64bit)

Set the recommended memory, will set to 12GB (12288 MB) for now

Create a virtual hard drive now--Fixed size in VDI format, 8GB called HitWalker2_test
Try: Create HitWalker2_base dynamically allocated 50GB

Start the machine and click on the little CD icon in the lower right side

Select the iso image, restart if necessary.

Install. Making sure to select minimal virtual machine installation press F4 or fn->F4 on Macs

Select English

setup keyboard etc...

select a name, by default it is ubuntu

full and user name was set to hw_user

password: hw_user

No encryption

Use guided partitioning with entire disk

No proxy and no automatic updates

Basic ubuntu server--for now...

sudo apt-get install git

#get HitWalker2

cd ~/hw_user

git clone https://github.com/HitWalker2.git

#install necessary dependencies for python etc

sudo apt-get install emacs24-nox

sudo apt-get install git nginx python-pip python-dev build-essential  python-numpy python-scipy python-cairosvg python-lxml libxml2-devel.x86_64 libxslt-devel.x86_64 wget openjdk-7-jdk vim

sudo pip install --upgrade pip 

sudo pip install Django==1.5.1
sudo pip install py2neo==1.6.4

sudo pip install tinycss cssselect cssutils colour

sudo pip install pandas gunicorn eventlet


sudo mkdir -p /var/www/hitwalker2_inst

##install neo4j, as of now v2.1.7

wget -O temp.key http://debian.neo4j.org/neotechnology.gpg.key

sudo apt-key add temp.key

rm temp.key

sudo vim /etc/apt/sources.list.d/neo4j.list

#added:

deb http://debian.neo4j.org/repo stable/

 apt-get update
 
 sudo apt-get install neo4j
 
 #may need:
  service neo4j-service start
  
  #remove limit on number of open files:
  
  #via: http://neo4j.com/docs/2.1.7/linux-performance-guide.html
  
  sudo su -
 
 vim /etc/security/limits.conf
 
 #add:
 
#neo4j   soft    nofile  40000
#neo4j   hard    nofile  40000
 
 
 vim /etc/pam.d/su
 
 #add
 
 #session    required   pam_limits.so
 
 
#Add a shared folder--macs, right click on VM box->settings->Shared Folders->e.g. /Users/bottomly/Desktop/hitwalker2_paper
 
 #First install Guest Extensions to support shared folders
 
 sudo apt-get install dkms build-essential linux-headers-generic
 
 #click on the VM window->devices->insert guest additions cd
 
 sudo mount /dev/sr0 /media/cdrom
 
 cd /media/cdrom
 
 sudo sh ./VBoxLinuxAdditions.run
 
cd ..;sudo umoun cdrom
 
 sudo mkdir /mnt/share
 
 sudo mount -t vboxsf hitwalker2_paper /mnt/share
 
 #get the base data
 
 cd /mnt/share
 
 cp -r hitwalker2_base_data /home/hw_user/
 

#set the data directory for neo4j to this location

sudo service neo4j-service stop

sudo mv /var/lib/neo4j/data /var/lib/neo4j/data-old


sudo ln -s /hom	e/hw_user/hitwalker2_base_data /var/lib/neo4j/data

cd /home/hw_user

#from /var/lib/neo4j
sudo chown -R neo4j:adm data

#from /home/hw_user
sudo chown -R neo4j:adm hitwalker2_base_data

sudo service neo4j-service start


#does it work:
neo4j-shell -c 'match (n) return n limit 10'


#networking

#shutdown the server

#go to settings->networking and select 'bridged adaptor'

#restart

#Can access machine through the ifconfig ip of the host computer...


#Install R in /home as that is where most of the space resides...

wget http://cran.us.r-project.org/src/base/R-3/R-3.1.3.tar.gz

tar -xzvf R-3.1.3.tar.gz

cd R-3.1.3

sudo apt-get install gfortran libreadline-dev

./configure --with-x=no

make

vim ~/.bash_profile
export PATH=/home/hw_user/R-3.1.3/bin:$PATH

R

source("http://bioconductor.org/biocLite.R")

biocLite(c("igraph", "reshape2", "Biobase", "rjson"))

#install HitWalker2

git clone https://github.com/biodev/HitWalker2.git HitWalker2

#cloned at this point--full clone

#started hitwalker_2_test_image

#also need to install a few additional packages for ccle

sudo apt-get install libxml2-dev libcurl4-dev

source("http://bioconductor.org/biocLite.R")

biocLite(c("SCAN.UPC","hgu133plus2.db" ))

#back to bash

git clone https://github.com/mlbernauer/Entrez.git

install.packages("tm")

R CMD INSTALL Entrez


#set up upstart for unicorn

in file: /etc/init/HitWalker2.conf

description "HitWalker2"

start on (filesystem)
stop on runlevel [016]

respawn
setuid hw_user
setgid hw_user
chdir /home/hw_user/HItWalker2

exec gunicorn -k 'eventlet' HitWalker2.wsgi:application

from /home/hw_user/HitWalker2

#probably not the best approach...
sudo chmod -R 777 /var/www

python manage.py collectstatic

sudo mkdir /etc/nginx/ssl

cd /etc/nginx/ssl

sudo openssl genrsa -des3 -out server.key 1024

passphrase hw_user

openssl req -new -key server.key -out server.csr

#this removes the need for the password on startup

cp server.key server.key.org
openssl rsa -in server.key.org -out server.key

openssl x509 -req -days 365 -in server.csr -signkey server.key -out server.crt

#from /home/hw_user/HitWalker2

sudo cp default-nginx /etc/nginx/sites-available/

sudo ln -s /etc/nginx/sites-available/default-nginx /etc/nginx/sites-enabled/hw2


##mounting shared drive

sudo mount -t vboxsf -o uid=1000,gid=1000 ccle_data /mnt/share 

#From /home/hw_user/HitWalker2/populate/

 cp -r ccle_example.Rnw /mnt/share/
 
 #from /mnt/share
 
 R CMD Sweave --pdf ccle_example.Rnw