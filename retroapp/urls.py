from django.urls import path, re_path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    re_path(r'(\S+)/(\d{1,4})$', views.retrotide_usage, name='smilestr'),
]