from django.shortcuts import redirect
from django.urls import path, re_path, include
from django.contrib.auth.views import LogoutView

from retroapp import views

urlpatterns = [
    path('', views.index, name='index'),
    # path('', lambda request: redirect('home', permanent=False), name='index'),
    
    path('home', views.home, name='home'),
    path('search', views.search, name='search'),
    
    # TODO: Use this to get results of previous queries from the database.
    # path('pks', views.pks_search_result, name='pks_search_result'),
    
    # path('history', views.history, name='history'),
    path('history', views.QueryHistoryView.as_view(), name='history'),
    path('oldresult', views.QueryHistoryResultView.as_view(), name='history_results'),
    path('about', views.about, name='about'),
    re_path(r'(\S+)/(\d{1,4})$', views.retrotide_usage, name='smilestr'),
    
    path('accounts', include('allauth.urls'), name='oauth_accounts'),
    path('logout', LogoutView.as_view(), name='oauth_logout'),
]