from django.shortcuts import redirect
from django.urls import path, re_path

from retroapp import views

urlpatterns = [
    path('', views.index, name='index'),
    # path('', lambda request: redirect('home', permanent=False), name='index'),
    path('home', views.home, name='home'),
    path('search', views.search_lite, name='search'),
    path('pks', views.pks_lite, name='pks'),
    path('about', views.about, name='about'),
    re_path(r'(\S+)/(\d{1,4})$', views.retrotide_usage, name='smilestr'),
]