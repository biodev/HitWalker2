from django.conf.urls import patterns, include, url
from django.conf import settings
from django.contrib.staticfiles.urls import staticfiles_urlpatterns

from django.contrib import admin
admin.autodiscover()

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    url(r'^HitWalker2$', include('network.urls')),
    url(r'^HitWalker2/', include('network.urls'))
    # Examples:
    # url(r'^$', 'HitWalker2.views.home', name='home'),
    # url(r'^HitWalker2/', include('HitWalker2.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    #url(r'^admin/', include(admin.site.urls)),
)

urlpatterns += staticfiles_urlpatterns()
