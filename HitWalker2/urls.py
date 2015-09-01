from django.conf.urls import patterns, include, url
from django.conf import settings
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from network.config import prog_type

#from django.contrib import admin
#admin.autodiscover()

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

if prog_type != "":
    urlpatterns=patterns('',
    url(r'^HitWalker2/'+prog_type+'$', include('network.urls')),
    url(r'^HitWalker2/'+prog_type+'/', include('network.urls')))
else:
    urlpatterns=patterns('',
    url(r'^HitWalker2$', include('network.urls')),
    url(r'^HitWalker2/', include('network.urls')))

urlpatterns += staticfiles_urlpatterns()
