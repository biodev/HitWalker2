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

Select the iso image and it will run.

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

sudo apt-get install nginx
sudo apt-get install python-pip python-dev build-essential
sudo pip install --upgrade pip 

sudo pip install Django==1.5.1
sudo pip install py2neo==1.6.4

sudo apt-get install python-numpy python-scipy

sudo apt-get install python-cairosvg python-lxml


sudo yum install libxml2-devel.x86_64
sudo yum install libxslt-devel.x86_64

sudo pip install tinycss cssselect cssutils colour

sudo pip install pandas gunicorn eventlet
