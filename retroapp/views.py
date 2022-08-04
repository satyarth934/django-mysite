from collections import defaultdict
from django.shortcuts import render, redirect
# from django.urls import path, re_path

# from django.http import Http404, HttpResponse, JsonResponse
from django.http import HttpResponse, JsonResponse
from rdkit import Chem
# from rdkit.Chem import Draw, rdFMCS, AllChem

# from django.utils.http import urlunquote
# from urllib.parse import unquote as urlunquote

# from django.views.decorators.cache import cache_page
# import re
# import xml.etree.ElementTree as ET
# from io import BytesIO

import pandas as pd
from django.template import loader
from django.forms import formset_factory
from retroapp import constants
from retroapp.constants import ASCENDING, DESCENDING, NO_SORT, MOLECULE_PROPERTIES, SORTING_OPTIONS

import retrotide
from pprint import pprint
from mysite import utils
import warnings
warnings.formatwarning = utils.warning_format    # Defining a specific format to print warnings on console

import numpy as np

# import retroapp.urls as retroapp_urls
from retroapp.forms import SMILESForm, PropertyForm

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


def search_lite(request):
    smiles_form = SMILESForm(
        initial={
            'smiles_string': None,
            'notes': None,
        },
    )
    PropertyFormSet = formset_factory(PropertyForm)

    # If this is a GET (or any other method) create the default form.
    if request.method != "POST":
        pfset = PropertyFormSet()
        context = {
            "smiles_form": smiles_form,
            "property_formset": pfset,
        }
        return render(request, "retroapp/search.html", context)
    
    else:
        smiles_form = SMILESForm(request.POST)
        property_formset = PropertyFormSet(request.POST)

        if smiles_form.is_valid():
            print("[VALID] smiles_form!!")
            print(smiles_form.cleaned_data)
        else:
            print("Invalid smiles_form!!")
            pprint(smiles_form.__dict__)
        
        if property_formset.is_valid():
            print("[VALID] property_formset")
            for i, form in enumerate(property_formset):
                print(f"---> FORM {i}")
                print(form.cleaned_data)
        else:
            print("Invalid property_Formset!!")
            print(property_formset.errors)
        
        # return render(request, "formsetapp/index.html", {})
        if smiles_form.is_valid() and property_formset.is_valid():
            request.session['_search_query_smiles'] = smiles_form.cleaned_data
            request.session['_search_query_properties'] = [form.cleaned_data for form in property_formset]
            return redirect('pks', permanent=False)


def pks_lite(request):
    if "_search_query" in request.session:
        request.session.set_expiry(value=0)    # user’s session cookie will expire when the user’s web browser is closed
        # search_query = request.session.get('_search_query')
        search_query_smiles = request.session['_search_query_smiles']
        search_query_properties = request.session['_search_query_properties']
        print(f"{search_query_smiles = }")
        print(f"{search_query_properties = }")
        search_query = dict()
        search_query["smiles_string"] = search_query_smiles["smiles_string"]
        search_query["notes"] = search_query_smiles["notes"]
        search_query["properties"] = search_query_properties
        print(f"{search_query = }")

        # Calling retrotide API
        mol_properties = [d["molecular_property"] for d in search_query["properties"]]
        retro_df = retrotide_call(smiles=search_query["smiles_string"], properties=mol_properties)

        # Update the dataframe based on min-max ranges of all properties.
        for search_query_property in search_query["properties"]:
            mol_property = search_query_property["molecular_property"]
            min_val = constants.MOLECULE_PROPERTIES[mol_property]['min']
            max_val = constants.MOLECULE_PROPERTIES[mol_property]['max']
            search_query_property["target_val"] = None

            property_val_range = search_query_property["property_value_range"]
            if property_val_range:
                if isrange(property_val_range):    # Range is specified
                    min_val, max_val = [float(range_val) for range_val in property_val_range.split(" - ")]
                else:    # target value is specified
                    search_query_property["target_val"] = float(property_val_range)                

            retro_df = retro_df[(retro_df[mol_property]>=min_val ) & (retro_df[mol_property]<=max_val)]
        
        # Get all the properties that have a sorting_mode specified or have a target value.
        sorting_cols = [
            sqprop["molecular_property"] for sqprop in search_query["properties"] if sqprop["sorting_mode"] in [ASCENDING, DESCENDING] or not isrange(sqprop["property_value_range"])
        ]

        # Bool list to specify property sorting in ascending order.
        sorting_cols_asc = list()
        for sqprop in search_query["properties"]:
            if not isrange(sqprop["property_value_range"]):
                sorting_cols_asc.append(True)
            else:
                sorting_cols_asc.append(sqprop["sorting_mode"] == ASCENDING)

        # Target values corresponding to each property
        prop_targets = [sqprop["target_val"] for sqprop in search_query["properties"]]

        # function defining the key for sorting
        def sorting_key(pd_series):
            idx = sorting_cols.index(pd_series.name)
            target = prop_targets[idx]
            if target is None:
                return pd_series
            
            target = [target] * len(pd_series)
            return np.abs(np.array(target) - pd_series)

        retro_df = retro_df.sort_values(
            by=sorting_cols, 
            ascending=sorting_cols_asc, 
            key=sorting_key
        )

        # To reset the index value to start from 0 again
        retro_df = retro_df.reset_index(drop=True)

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

def predict_property_api(property):
    warnings.warn(f"Populating {property} with RANDOM values.")    # DELETE
    return np.random.randint(0,100)


def retrotide_call(smiles, properties=None):
    designs = retrotide.designPKS(Chem.MolFromSmiles(smiles))

    if properties is not None:    # Convert to a list if properties is not None
        properties = [properties] if isinstance(properties, str) else properties
        
    else:
        properties = []
        warnings.warn(f"No properties mentioned.")

    df_dict = defaultdict(list)

    for i in range(len(designs[-1])):
        df_dict["SMILES"].append(Chem.MolToSmiles(designs[-1][i][2]))
        df_dict["Retrotide_Similarity_SCORE"].append(designs[-1][i][1])
        df_dict["DESIGN"].append(designs[-1][i][0].modules)
        for property in properties:
            df_dict[property].append(predict_property_api(property=property))

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


def isrange(property_value_str):
    return " - " in property_value_str