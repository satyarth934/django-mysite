from django.urls import path, re_path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('home', views.home, name='home'),
    path('search', views.search_form, name='search'),
    re_path(r'(\S+)/(\d{1,4})$', views.retrotide_usage, name='smilestr'),
]