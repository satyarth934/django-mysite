###########
# Imports #
###########
import warnings
from mysite import utils
warnings.formatwarning = utils.warning_format    # Defining a specific format to print warnings on console

import logging
logger = logging.getLogger(f"{__name__}")

import numpy as np
import pandas as pd
from pprint import pprint
from collections import defaultdict

from django.shortcuts import render, redirect
from django.http import HttpResponse
from django.template import loader
from django.forms import formset_factory
from django.views.generic.base import TemplateView

from retroapp import constants
from retroapp.constants import ASCENDING, DESCENDING, NO_SORT, MOLECULE_PROPERTIES, SORTING_OPTIONS
from retroapp.models import QueryDB, QueryPropertyDB, QueryResultsDB
from retroapp.forms_2 import FormQuery, FormQueryProperty
from retroapp import utils

import retrotide
from rdkit import Chem

######################
# Webpage views here #
######################

@utils.log_function
def index(request):
    """View function for index (blank) URL.
    """
    return redirect('home', permanent=False)


@utils.log_function
def home(request):
    """View function for home URL.
    """
    template = loader.get_template("retroapp/home.html")
    context = {}
    rendered_str = template.render(context, request)
    return HttpResponse(rendered_str)


@utils.log_function
def search(request):
    """View function for search URL.
    Renders ERROR page if the user is not logged in.
    Renders an empty form if the search page is opened.
    When the submit button is clicked:
    - Saves the search query in the database.
    - Renders the result below the form. (To be modified)
    - Render an information page below the form informing that the search job is submitted for the specified query. (Replacement for the above feature)
    """
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
    # Initialize the query form
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
    
    else:    # If this is a POST request, update database and invoke search.
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

            Q_uuid_fk = QueryDB.objects.latest('Timestamp')

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
            request.session['_form_state'] = request.POST
            request.session['_search_query_smiles'] = smiles_form.cleaned_data
            request.session['_search_query_properties'] = [form.cleaned_data for form in property_formset]
            # return redirect('pks', permanent=False)
            return pks_search_result(request)


# View function for search results in the same search page
@utils.log_function
def pks_search_result(request):
    """[Deprecated]
    This function is called when the search query is submitted.
    Tasks performed by this function:
    1. Getting retrotide results. (Dummy values for all the properties. Only RON value is not dummy)
    2. Apply filtering constraints on the results.
    3. Apply specified sorting on the results.
    4. Add resulting data to QueryResultsDB.
    5. Render the results as HTML page to show below the query form.

    [To be modified]
    Above mentioned logic will be modified to incorporate SF API. 
    The results will no longer directly be available.
    The user will be notified that the query request has been submitted. 
    The query status and the results can be checked on the `History` webpage.

    # TODO: Implement the above.
    """
    
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
        mol_properties = [
            list(constants.MOLECULE_PROPERTIES.keys())[d["Property_name"]] 
            for d in search_query["properties"]
        ]
        retro_df = retrotide_call(smiles=search_query["smiles_string"], properties=mol_properties)
        
        # Update the dataframe based on min-max ranges of all properties.
        for search_query_property in search_query["properties"]:
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
            key=sorting_key,
        )

        # To reset the index value to start from 0 again
        retro_df = retro_df.reset_index(drop=True)

        # TODO: map the sulfur replacement function here for all the molecules.
        # TODO: Find an alternate more efficient implementation. This will take around 5-7 minutes for a million rows.
        retro_df['SMILES'] = retro_df['SMILES'].map(clean_sulfur_from_smiles_string)
        
        # Write results to QueryResultsDB
        retro_df_dict = retro_df.to_dict('records')

        Q_uuid_fk = QueryDB.objects.latest('Timestamp')

        model_instances = [
            QueryResultsDB(
                Q_uuid_id=Q_uuid_fk.Q_uuid,
                SMILES=record['SMILES'],
                Retrotide_Similarity_SCORE=record['Retrotide_Similarity_SCORE'],
                DESIGN=record['DESIGN'],
                Cetane_number=record.get('Cetane Number', None),
                Research_octane_number=record.get('Research Octane Number', None),
                Melting_point=record.get('Melting Point', None),
                Flash_point=record.get('Flash Point', None),
                Yield_sooting_index=record.get('Yield Sooting Index', None),
                H1_receptor_pKd=record.get('H1 Receptor pKd', None),
                M2_receptor_pKd=record.get('M2 Receptor pKd', None),
            ) for record in retro_df_dict
        ]

        QueryResultsDB.objects.bulk_create(model_instances)

        # Render return webpage
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


class QueryHistoryView(TemplateView):
    """View for history page.
    The submitted queries and their results can be accessed from here.
    """

    template_name = "retroapp/history.html"

    @utils.log_function
    def get_context_data(self, **kwargs):
        # User IS NOT authenticated
        if not self.request.user.is_authenticated:
            context = {
                "message": "Please login to view your past queries!!",
                "message_tag": "ERROR",
                # "remove_tabs": True,
            }
            return context

        # User IS authenticated
        context = super(QueryHistoryView, self).get_context_data(**kwargs)
        context['MOLECULE_PROPERTIES'] = list(constants.MOLECULE_PROPERTIES.keys())
        context['SORTING_OPTIONS'] = [s_optn[1] for s_optn in constants.SORTING_OPTIONS]

        queries = QueryDB.objects.filter(Username=self.request.user).order_by('-Timestamp')

        context['query_object_list'] = list()
        for q in queries:
            q_properties = QueryPropertyDB.objects.filter(Q_uuid=q.Q_uuid)
            context['query_object_list'].append((q, q_properties))

        return context


# # Placeholder for the history page (contains list of user's all old query runs)
# @utils.log_function
# def history(request):
#     if request.user.is_authenticated:
#         # messages.info(request, "Your past queries...")
#         context = {
#             "message": "Your past queries...",
#             "message_tag": "INFO",
#         }
#         ...
#         return render(request, "retroapp/history.html", context)

#         # return QueryHistoryView.as_view()
    
#     else:
#         # messages.error(request, "Please login to view your past queries!!")
#         context = {
#             "message": "Please login to view your past queries!!",
#             "message_tag": "ERROR",
#             # "remove_tabs": True,
#         }
#         return render(request, "retroapp/guest_landing_page.html", context)


class QueryHistoryResultView(TemplateView):
    """View for history results page.
    The results for the submitted queries can be accessed from here.
    Each submitted SMILES string has a hyperlink to this results page for that query request.
    """

    template_name = "retroapp/history_result.html"

    @utils.log_function
    def get_context_data(self, **kwargs):
        # User IS NOT authenticated
        if not self.request.user.is_authenticated:
            context = {
                "message": "Please login to view your past queries!!",
                "message_tag": "ERROR",
                # "remove_tabs": True,
            }
            return context

        # User IS authenticated
        context = super(QueryHistoryResultView, self).get_context_data(**kwargs)

        query_res_object = QueryResultsDB.objects.filter(Q_uuid_id=self.request.GET['name'])
        context["query_res_object"] = query_res_object

        return context


# View function for about URL
@utils.log_function
def about(request):
    """View function for about URL.
    """

    template = loader.get_template("retroapp/about.html")
    context = {}
    rendered_str = template.render(context, request)
    return HttpResponse(rendered_str)


# View function for smilesstr URL
@utils.log_function
def retrotide_usage(request, smiles, width=243):
    """Django wrapper for using retrotide API.
    """

    # retrotide API call
    # retro_df = retrotideAPI_dummy(request, smiles)
    retro_df = retrotide_call(smiles=smiles)

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
    """Dummy function call to replicate the molecule property prediction.
    Returns random values for the query properties requested.
    """

    warnings.warn(f"Populating {property} with RANDOM values.")    # DELETE
    return np.random.randint(0,100)


@utils.log_function
def retrotide_call(smiles, properties=None):
    """Retrotide API call.
    """

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

    try:
        output_df = pd.DataFrame(df_dict)
    except ValueError as ve:
        raise ValueError(f"{ve} - Check if multiple constraints are defined for the same Molecular Property.")

    return pd.DataFrame(df_dict)


def isrange(property_value_str):
    """Checks whether the specified property is a range or not.

    Args:
        property_value_str (str): Input property value string from the form.

    Returns:
        bool: True if the property is a range, False otherwise.
    """

    return " - " in property_value_str


def clean_sulfur_from_smiles_string(input_smiles, idx=0):
    """Cleans sulfur molecule from the SMILES string.
    """

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