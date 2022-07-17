from collections import defaultdict
from django.shortcuts import render, redirect
from django.urls import path, re_path

from django.http import Http404, HttpResponse, JsonResponse
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS, AllChem

# from django.utils.http import urlunquote
from urllib.parse import unquote as urlunquote

from django.views.decorators.cache import cache_page
import re
import xml.etree.ElementTree as ET
from io import BytesIO

import pandas as pd
from django.template import loader
from retroapp.constants import ASCENDING, NO_SORT

import retrotide
from pprint import pprint
from mysite import utils
import warnings
warnings.formatwarning = utils.warning_format    # Defining a specific format to print warnings on console

import numpy as np

import retroapp.urls as retroapp_urls
from retroapp.forms import SearchForm

######################
# Webpage views here #
######################

def index(request):
    # url_patterns = [str(u.pattern) for u in retroapp_urls.urlpatterns]
    # url_pattern_names = [u.name for u in retroapp_urls.urlpatterns]
    # print(url_pattern_names)
    # return HttpResponse(f"You're at retroapp. Append smile string in the URL to invoke the retrotide API.<br><br>Patterns:<br>{url_patterns}<br><br>Extension calls:<br>{url_pattern_names}")

    return redirect('home', permanent=False)


def home(request):
    template = loader.get_template("retroapp/home.html")
    context = {}
    rendered_str = template.render(context, request)
    return HttpResponse(rendered_str)


def search(request):
    # If this is a GET (or any other method) create the default form.
    if request.method != "POST":
        form = SearchForm(
            initial={
                'smiles_string': "",
                'molecular_property': "Cetane Number",
                'molecular_property_min': "",
                'molecular_property_max': "",
                'sorting_mode': "Higher is better",
            }
        )
        context = {
            'form': form,
        }
        return render(request, "retroapp/search.html", context=context)
    
    else:
        form = SearchForm(request.POST)
        if form.is_valid():
            request.session['_search_query'] = form.cleaned_data
            return redirect('pks', permanent=False)


def pks(request):
    if "_search_query" in request.session:
        request.session.set_expiry(value=0)    # user’s session cookie will expire when the user’s web browser is closed
        search_query = request.session.get('_search_query')

        # Calling retrotide API
        retro_df = retrotide_call(smiles=search_query["smiles_string"], properties=search_query["molecular_property"])

        # Filtering and sorting based on search query
        min_val = search_query["molecular_property_min"]
        max_val = search_query["molecular_property_max"]
        target_val = search_query["molecular_property_target"]
        mol_property = search_query["molecular_property"]
        retro_df = retro_df[(retro_df[mol_property]>=min_val ) & (retro_df[mol_property]<=max_val)]

        sort_optn = search_query["sorting_mode"]
        if sort_optn != NO_SORT:
            if target_val is None:
                retro_df = retro_df.sort_values(by=[mol_property], ascending=(sort_optn == ASCENDING))
            else:
                retro_df = retro_df.sort_values(by=[mol_property], ascending=(sort_optn == ASCENDING), key=lambda x: np.abs(target_val-x))

        # Render the filtered dataframe to a webpage
        table_rendered_str = showtable(
            request, 
            query_dict=search_query, 
            retro_df=retro_df,
        )
        return HttpResponse(table_rendered_str)

    else:
        warnings.warn("404: No search query!")
        return HttpResponse("404: No search query!")
        

def about(request):
    template = loader.get_template("retroapp/about.html")
    context = {}
    rendered_str = template.render(context, request)
    return HttpResponse(rendered_str)


def retrotide_usage(request, smiles, width=243):
    # retrotide API call
    # retro_df = retrotideAPI_dummy(request, smiles)
    retro_df = retrotide_call(smiles=smiles)
    print(retro_df)         # DELETE: only for debugging purposes

    # show table for the return data frame from the API call
    table_rendered_str = showtable(
        request, 
        query_dict={"smiles_string": smiles}, 
        retro_df=retro_df, 
        width=width
    )
    
    return HttpResponse(table_rendered_str)

########################
# Dependency functions #
########################

def retrotide_call(smiles, properties=None):
    designs = retrotide.designPKS(Chem.MolFromSmiles(smiles))

    if properties is not None:    # Convert to a list if properties is not None
        properties = [properties] if isinstance(properties, str) else properties
        warnings.warn(f"Populating {properties} with RANDOM values.")    # DELETE
    else:
        properties = []
        warnings.warn(f"No properties mentioned.")

    df_dict = defaultdict(list)

    for i in range(len(designs[-1])):
        df_dict["SMILES"].append(Chem.MolToSmiles(designs[-1][i][2]))
        df_dict["Retrotide_Similarity_SCORE"].append(designs[-1][i][1])
        df_dict["DESIGN"].append(designs[-1][i][0].modules)
        for property in properties:
            df_dict[property].append(np.random.randint(0,100))

    return pd.DataFrame(df_dict)


def showtable(request, query_dict, retro_df, width=243):
    template = loader.get_template("retroapp/showtable.html")
    context = {
        "query_dict": query_dict,
        "keys": ["Rendered Molecule", *retro_df.keys()],
        "df": retro_df,
        "width": width,
        }
    rendered_str = template.render(context, request)
    return rendered_str











###########
# Archive #
###########

def retrotideAPI_dummy(request, smiles):
    example_df = pd.read_csv("retroapp/example_df.tab", delimiter="\t")
    return example_df



from django.http import JsonResponse
from retroapp.forms import TestForm


def form_test_view(request):
    """View function for form testing."""
    print("Testing the form view.")

    # If this is a POST request then process the Form data
    if request.method == 'POST':

        # Create a form instance and populate it with data from the request (binding):
        form = TestForm(request.POST)

        # Check if the form is valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            query = dict()
            query['data'] = form.cleaned_data['test_form_data']

            # redirect to a new URL:
            return JsonResponse(query)

    # If this is a GET (or any other method) create the default form.
    else:
        form = TestForm(initial={'test_form_data': ">>> Some default  value <<<"})

    context = {
        'form': form,
        'delete': "delete stuff",
    }

    return render(request, 'retroapp/form_test_view.html', context)
