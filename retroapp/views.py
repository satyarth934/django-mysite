from collections import defaultdict
from email import message
from django.shortcuts import render, redirect

from django.http import HttpResponse, JsonResponse
from django.contrib import messages

from rdkit import Chem

import pandas as pd
from django.template import loader
from django.forms import formset_factory, modelformset_factory
from retroapp import constants
from retroapp.constants import ASCENDING, DESCENDING, NO_SORT, MOLECULE_PROPERTIES, SORTING_OPTIONS

import retrotide
from pprint import pprint
from mysite import utils
import warnings
warnings.formatwarning = utils.warning_format    # Defining a specific format to print warnings on console

from datetime import datetime
import numpy as np
import logging
logger = logging.getLogger(f"{__name__}")

from retroapp.models import QueryDB, QueryPropertyDB, QueryResultsDB
from retroapp.forms_2 import FormQuery, FormQueryProperty
from retroapp import utils

######################
# Webpage views here #
######################

# View function for index (blank) URL
@utils.log_function
def index(request):
    return redirect('home', permanent=False)


# View function for home URL
@utils.log_function
def home(request):

    # DELETE: only for debugging purposes
    if request.user.is_authenticated:
        print("request.user:", request.user)
        print("request.user.is_authenticated:", request.user.is_authenticated)
        print("request.user.socialaccount_set.all()[0].get_avatar_url():", request.user.socialaccount_set.all()[0].get_avatar_url())

    template = loader.get_template("retroapp/home.html")
    context = {}
    rendered_str = template.render(context, request)
    return HttpResponse(rendered_str)


# View function for search URL
@utils.log_function
def search(request):
    #####################################################################
    # When user is NOT logged in. 
    #####################################################################
    if not request.user.is_authenticated:
        context = {
            "message": "Please login to use the Search functionality!!",
            "message_tag": "ERROR",
            # "remove_tabs": True,
        }

        return render(request, "retroapp/guest_landing_page.html", context)
    
    #####################################################################
    # Open the Search Page only when the user is logged in. 
    #####################################################################
    # smiles_form = SMILESForm(
    smiles_form = FormQuery(
        initial={
            'Q_smiles': None,
            'Q_notes': None,
        },
    )
    PropertyFormSet = formset_factory(FormQueryProperty)

    # If this is a GET (or any other method) create the default form.
    if request.method != "POST":
        pfset = PropertyFormSet()
        context = {
            "smiles_form": smiles_form,
            "property_formset": pfset,
        }
        return render(request, "retroapp/search.html", context)
    
    else:
        smiles_form = FormQuery(request.POST)
        property_formset = PropertyFormSet(request.POST)
        property_formset.form_kwargs["empty_permitted"] = False

        if smiles_form.is_valid():
            logger.info("[VALID] smiles_form!!")
            logger.info(smiles_form.cleaned_data)

            smiles_form_to_db = smiles_form.save(commit=False)
            smiles_form_to_db.Username = request.user
            smiles_form_to_db.Full_name = request.user.socialaccount_set.all()[0].extra_data['name']
            smiles_form_to_db.Email = request.user.socialaccount_set.all()[0].extra_data['email']
            
            smiles_form_to_db.save()

        else:
            logger.error("[INVALID] smiles_form!!")
            logger.info(smiles_form.__dict__)
        
        if property_formset.is_valid():
            logger.info("[VALID] property_formset")

            Q_uuid_fk = QueryDB.objects.last()

            for i, form in enumerate(property_formset):
                logger.info(f"---> FORM {i}")
                logger.info(form.cleaned_data)

                form_to_db = form.save(commit=False)
                if isrange(form.cleaned_data["property_value_range"]):
                    min_val, max_val = form.cleaned_data['property_value_range'].split(" - ")
                    form_to_db.Min_value = float(min_val)
                    form_to_db.Max_value = float(max_val)
                    
                elif (form.cleaned_data["property_value_range"] != "") and (form.cleaned_data["property_value_range"] is not None):
                    form_to_db.Target_value = float(form.cleaned_data["property_value_range"])
                else:
                    pass
                form_to_db.Q_uuid = Q_uuid_fk
                form_to_db.save()            

        else:
            logger.error("[INVALID] property_Formset!!")
            logger.info(property_formset.errors)
        
        # return render(request, "formsetapp/index.html", {})
        if smiles_form.is_valid() and property_formset.is_valid():
            request.session['_search_query_smiles'] = smiles_form.cleaned_data
            request.session['_search_query_properties'] = [form.cleaned_data for form in property_formset]
            # return redirect('pks', permanent=False)
            return pks_search_result(request)


# View function for search results in the same search page
@utils.log_function
def pks_search_result(request):
    if ("_search_query_smiles" in request.session) and ("_search_query_properties" in request.session):
        request.session.set_expiry(value=0)    # user’s session cookie will expire when the user’s web browser is closed
        # search_query = request.session.get('_search_query')
        search_query_smiles = request.session['_search_query_smiles']
        search_query_properties = request.session['_search_query_properties']
        logger.info(f"search_query_smiles = {search_query_smiles}")
        logger.info(f"search_query_properties = {search_query_properties}")
        search_query = dict()
        search_query["smiles_string"] = search_query_smiles["Q_smiles"]
        search_query["notes"] = search_query_smiles["Q_notes"]
        search_query["properties"] = search_query_properties
        logger.info(f"search_query = {search_query}")

        # Calling retrotide API
        # mol_properties = [d["molecular_property"] for d in search_query["properties"]]
        # mol_properties = [d["Property_name"] for d in search_query["properties"]]
        mol_properties = [
            list(constants.MOLECULE_PROPERTIES.keys())[d["Property_name"]] 
            for d in search_query["properties"]
        ]
        retro_df = retrotide_call(smiles=search_query["smiles_string"], properties=mol_properties)
        
        # Update the dataframe based on min-max ranges of all properties.
        for search_query_property in search_query["properties"]:
            # mol_property = search_query_property["molecular_property"]
            # mol_property = search_query_property["Property_name"]
            mol_property = list(constants.MOLECULE_PROPERTIES.keys())[search_query_property["Property_name"]]
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
        
        # Get all properties and their corresponding sorting modes (bool value)
        sorting_cols = list()
        sorting_cols_asc = list()
        for sqprop in search_query["properties"]:
            property_name_str = list(constants.MOLECULE_PROPERTIES.keys())[sqprop["Property_name"]]
            
            # Target value is defined (always ascending)
            if (sqprop["property_value_range"] != "") and (not isrange(sqprop["property_value_range"])):
                sorting_cols.append(property_name_str)
                sorting_cols_asc.append(True)
            
            else:    # Either blank or range but not target value
                # Sorting mode is defined (ascending | descending)
                if (sqprop["Sorting_mode"] in [ASCENDING, DESCENDING]):
                    sorting_cols.append(property_name_str)
                    sorting_cols_asc.append(sqprop["Sorting_mode"] == ASCENDING)

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

        # TODO: map the sulfur replacement function here for all the molecules.
        # TODO: Find an alternate more efficient implementation. This will take around 5-7 minutes for a million rows.
        # with utils.ExecTimeCM("pd map") as map_et:
        retro_df['SMILES'] = retro_df['SMILES'].map(clean_sulfur_from_smiles_string)
        
        # Write results to DB
        print("--- retro_df ---")
        print(retro_df.columns)
        print("--- xxxxxxxx ---")
        retro_df_dict = retro_df.to_dict('records')

        # TODO: this is currently breaking!!!
        model_instances = [QueryResultsDB(
            Q_uuid_id=QueryDB.objects.last(),
            SMILES=record['SMILES'],
            Retrotide_Similarity_SCORE=record['Retrotide_Similarity_SCORE'],
            DESIGN=record['DESIGN'],
            Cetane_number=None if 'Cetane_number' not in record else record['Cetane_number'],
            Research_octane_number=None if 'Research_octane_number' not in record else record['Research_octane_number'],
            Melting_point=None if 'Melting_point' not in record else record['Melting_point'],
            Flash_point=None if 'Flash_point' not in record else record['Flash_point'],
            Yield_sooting_index=None if 'Yield_sooting_index' not in record else record['Yield_sooting_index'],
            H1_receptor_pKd=None if 'H1_receptor_pKd' not in record else record['H1_receptor_pKd'],
            M2_receptor_pKd=None if 'M2_receptor_pKd' not in record else record['M2_receptor_pKd'],
        ) for record in retro_df_dict]

        QueryResultsDB.objects.bulk_create(model_instances)

        context = {
            "query_dict": search_query,
            "keys": ["Rendered Molecule", *retro_df.keys()],
            "df": retro_df,
            "width": 243,
        }
        return render(request, "retroapp/show_results.html", context)

    else:
        logger.warn("404: No search query!")
        return HttpResponse("404: No search query!")


# Placeholder for the history page (contains list of user's all old query runs)
@utils.log_function
def history(request):
    # TODO: Implement this view function
    # messages.info(request, "Checking django msg.")

    if request.user.is_authenticated:
        # messages.info(request, "Your past queries...")
        context = {
            "message": "Your past queries...",
            "message_tag": "INFO",
        }
        ...
    
    else:
        # messages.error(request, "Please login to view your past queries!!")
        context = {
            "message": "Please login to view your past queries!!",
            "message_tag": "ERROR",
            # "remove_tabs": True,
        }

    return render(request, "retroapp/guest_landing_page.html", context)


# View function for about URL
@utils.log_function
def about(request):
    template = loader.get_template("retroapp/about.html")
    context = {}
    rendered_str = template.render(context, request)
    return HttpResponse(rendered_str)


# View function for smilesstr URL
@utils.log_function
def retrotide_usage(request, smiles, width=243):
    # retrotide API call
    # retro_df = retrotideAPI_dummy(request, smiles)
    retro_df = retrotide_call(smiles=smiles)
    print(retro_df)         # DELETE: only for debugging purposes

    context = {
        "query_dict": {"smiles_string": smiles},
        "keys": ["Rendered Molecule", *retro_df.keys()],
        "df": retro_df,
        "width": width,
    }
    return render(request, "retroapp/showtable.html", context)

########################
# Dependency functions #
########################

def predict_property_api(property):
    warnings.warn(f"Populating {property} with RANDOM values.")    # DELETE
    return np.random.randint(0,100)


@utils.log_function
def retrotide_call(smiles, properties=None):
    designs = retrotide.designPKS(Chem.MolFromSmiles(smiles))

    if properties is not None:    # Convert to a list if properties is not None
        properties = [properties] if isinstance(properties, str) else properties
        
    else:
        properties = []
        logger.warn(f"No properties mentioned.")

    df_dict = defaultdict(list)

    for i in range(len(designs[-1])):
        df_dict["SMILES"].append(Chem.MolToSmiles(designs[-1][i][2]))
        df_dict["Retrotide_Similarity_SCORE"].append(designs[-1][i][1])
        df_dict["DESIGN"].append(designs[-1][i][0].modules)
        for property in properties:
            df_dict[property].append(predict_property_api(property=property))

    # DELETE
    # print("--- properties ---")
    # print(properties)
    # print("--- df_dict ---")
    # for k in df_dict:
    #     print(k)
    #     print(len(df_dict[k]))
    # print("--- xxxxxxx ---")

    try:
        output_df = pd.DataFrame(df_dict)
    except ValueError as ve:
        raise ValueError(f"{ve} - Check if multiple constraints are defined for the same Molecular Property.")

    return pd.DataFrame(df_dict)


def isrange(property_value_str):
    return " - " in property_value_str


def clean_sulfur_from_smiles_string(input_smiles, idx=0):
    try:
        input_mol = Chem.MolFromSmiles(input_smiles)
        Chem.SanitizeMol(input_mol)
        rxn = Chem.AllChem.ReactionFromSmarts('[C:1](=[O:2])[S:3]>>[C:1](=[O:2])[O].[S:3]')
        product = rxn.RunReactants((input_mol,))[idx][0]
        Chem.SanitizeMol(product)
    except Exception as e:
        print("input_smiles:", input_smiles)
        raise e

    return Chem.MolToSmiles(product)