#install HitWalker2

#CentOS 6 version

From http://nginx.org/en/linux_packages.html#mainline:

Add: 

[nginx]
name=nginx repo
baseurl=http://nginx.org/packages/centos/6/$basearch/
gpgcheck=0
enabled=1

to /etc/yum.repos.d/nginx.repo (with sudo privliges)

then: sudo yum install nginx.x86_64

Install python 2.7, following this blog post: http://toomuchdata.com/2014/02/16/how-to-install-python-on-centos/

sudo yum groupinstall "Development tools"
sudo yum install zlib-devel bzip2-devel openssl-devel ncurses-devel sqlite-devel readline-devel tk-devel gdbm-devel db4-devel libpcap-devel xz-devel

wget http://python.org/ftp/python/2.7.6/Python-2.7.6.tar.xz
tar xf Python-2.7.6.tar.xz
cd Python-2.7.6
./configure --prefix=/usr/local --enable-unicode=ucs4 --enable-shared LDFLAGS="-Wl,-rpath /usr/local/lib"
make
sudo make altinstall

#install pip
sudo /usr/local/bin/python2.7 get-pip.py

sudo /usr/local/bin/pip2.7 install virtualenv

Install requirements in a virtualenv

#from /var/www

sudo virtualenv HitWalker2_sandbox

from /var/www/HitWalker2_sandbox/bin

source ./activate

sudo chown -R bottomly:bottomly HitWalker2_Sandbox/

pip install Django==1.5.1
pip install py2neo==1.6.4

pip install numpy
#may need to install lapack first...
pip install scipy

sudo yum install libffi-devel.x86_64
pip install cairosvg

sudo yum install libxml2-devel.x86_64
sudo yum install libxslt-devel.x86_64

pip install lxml

pip install tinycss

pip install cssselect

pip install pandas

pip install gunicorn

deactivate

#setting up ancilary files:

sudo mkdir /var/www/hitwalker_2_inst
#from /var/www/
sudo chmod -R 777 /var/www/hitwalker_2_inst

scp protein.links.detailed.v9.05.9606.mm.[mn]* noble.ohsu.edu:/var/www/hitwalker_2_inst/
scp graph_struct.json noble.ohsu.edu:/var/www/hitwalker_2_inst/

Setting up a Nginx proxy

cp HitWalker2/default-nginx /etc/nginx/conf.d/default-nginx.conf

Change server name and ssl info

##ssl info for noble..
From http://www.cyberciti.biz/faq/nginx-self-signed-certificate-tutorial-on-centos-redhat-linux/

sudo mkdir -p /etc/nginx/ssl/hitwalker2
cd /etc/nginx/ssl/hitwalker2

openssl genrsa -des3 -out self-ssl.key 2048
openssl req -new -key self-ssl.key -out self-ssl.csr

openssl x509 -req -days 365 -in self-ssl.csr -signkey self-ssl.key -out self-ssl.crt

Make sure that 443 and 80 are open via:

sudo system-config-firewall
sudo service iptables restart

Starting Djanog/gunicorn

source /var/www/HitWalker2_Sandbox/bin/activate
cd /path/to/HitWalker2

#for example
./change_hw2_instance.sh lls

python manage.py syncdb

In addition need an appropriate *mm.mtx and *mm.names files as well as graph_struct.json file specified in /var/www/hitwalker2_inst/
or wherever it is specified in network/config.py

nohup gunicorn HitWalker2.wsgi:application &

Adding users:

python manage.py shell

from django.contrib.auth.models import User

user = User.objects.create_user('user', 'user@place.com', 'password')
user.save()
