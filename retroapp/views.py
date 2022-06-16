from collections import defaultdict
from django.shortcuts import render

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

import retrotide

import retroapp.urls as retroapp_urls
from retroapp.forms import SearchForm


# Create your views here.
def index(request):
    url_patterns = [str(u.pattern) for u in retroapp_urls.urlpatterns]
    url_pattern_names = [u.name for u in retroapp_urls.urlpatterns]
    print(url_pattern_names)
    return HttpResponse(f"You're at retroapp. Append smile string in the URL to invoke the retrotide API.<br><br>Patterns:<br>{url_patterns}<br><br>Extension calls:<br>{url_pattern_names}")


def home(request):
    template = loader.get_template("retroapp/home.html")
    context = {}
    rendered_str = template.render(context, request)
    return HttpResponse(rendered_str)


# def search(request):
#     template = loader.get_template("retroapp/search.html")
#     context = {
#             "properties": [
#                 "Cetane Number",
#                 "Research Octane Number",
#                 "Melting Point",
#                 "Flash Point",
#                 "Yield Sooting Index",
#                 "H1 Receptor pKd",
#                 "M2 Receptor pKd",
#             ]
#         }
#     rendered_str = template.render(context, request)

#     return HttpResponse(rendered_str)


def search_form(request):
    if request.method == 'POST':
        form = SearchForm(request.POST)

        if form.is_valid():
            # redirect to a new URL:
            return JsonResponse(form.cleaned_data)
    
    # If this is a GET (or any other method) create the default form.
    else:
        form = SearchForm(
            initial={
                'smiles_string': "",
                'molecular_property': "Cetane Number",
                'molecular_property_min': "",
                'molecular_property_max': "",
                # 'sorting_mode': "Select Sorting",
            }
        )

    template = loader.get_template("retroapp/search_form.html")
    context = {
        'form': form,
    }
    rendered_str = template.render(context, request)

    return HttpResponse(rendered_str)


def retrotideAPI_dummy(request, smiles):
    example_df = pd.read_csv("retroapp/example_df.tab", delimiter="\t")
    return example_df


def retrotide_call(smiles):
    designs = retrotide.designPKS(Chem.MolFromSmiles(smiles))

    df_dict = defaultdict(list)

    for i in range(len(designs[-1])):
        df_dict["SMILES"].append(Chem.MolToSmiles(designs[-1][i][2]))
        df_dict["Retrotide_Similarity_SCORE"].append(designs[-1][i][1])
        df_dict["RON"].append("###")
        df_dict["DESIGN"].append(designs[-1][i][0].modules)

    return pd.DataFrame(df_dict)


def showtable(request, query, retro_df, width=243):
    template = loader.get_template("retroapp/showtable.html")
    context = {
        "query": query,
        "keys": ["Rendered Molecule", *retro_df.keys()],
        "df": retro_df,
        "width": width,
        }
    rendered_str = template.render(context, request)
    return rendered_str


def retrotide_usage(request, smiles, width=243):
    # http_page_str = f"smiles: {smiles}<br>"
    # http_page_str += f"width: {width}<br>"

    # retrotide API call
    # retro_df = retrotideAPI_dummy(request, smiles)
    retro_df = retrotide_call(smiles)
    print(retro_df)         # DELETE: only for debugging purposes

    # show table for the return data frame from the API call
    table_rendered_str = showtable(request, smiles, retro_df, width)
    
    return HttpResponse(table_rendered_str)











# vvvvvvv
# DELETE 
# vvvvvvv


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
