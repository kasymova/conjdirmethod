from django.conf.urls import url
from .views import index, surface

urlpatterns = [
    url(r'^conjdirmethod/$', index, name='index'),
    url(r'^surface/$', surface, name='surface')
]
