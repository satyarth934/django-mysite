from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('home', views.home, name='home'),
    path('search', views.search, name='search'),
    path('showtable/', views.showtable, name='showtable'),
    path('showtablejson/', views.showtablejson, name='showtablejson'),
]