###########
# Imports #
###########
from typing import List, Tuple, Optional, Union
import warnings
from mysite import utils
warnings.formatwarning = utils.warning_format    # Defining a specific format to print warnings on console

import logging
logger = logging.getLogger(f"{__name__}")

import copy
import numpy as np
import pandas as pd
from pprint import pprint
from collections import defaultdict

import uuid
from django.shortcuts import render, redirect
from django.http import HttpResponse
from django.template import loader
from django.forms import formset_factory
from django.views.generic.base import TemplateView
from django.db import models

from retroapp import constants
from retroapp.constants import MOLECULE_PROPERTIES, MOLECULE_PROPERTY_CONSTRAINTS, SORT_OPTIONS
from retroapp.models import QueryDB, QueryPropertyDB, QueryResultsDB
from retroapp.forms_2 import FormQuery, FormQueryProperty
from retroapp import views_utils as vutils
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

        # Inserting values in QueryDB
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
        
        # Inserting values in QueryPropertyDB
        if property_formset.is_valid():
            logger.info("[VALID] property_formset")

            Q_uuid_fk = QueryDB.objects.latest('Timestamp')

            for i, form in enumerate(property_formset):
                logger.info(f"---> FORM {i}")
                logger.info(form.cleaned_data)

                form_to_db = form.save(commit=False)
                if vutils.isrange(form.cleaned_data["property_value_range"]):
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
            request.session['_q_uuid'] = str(Q_uuid_fk.Q_uuid)    # Will be used to update the Q_Job_id and Q_Status.
            request.session['_search_query_smiles'] = smiles_form.cleaned_data
            request.session['_search_query_properties'] = [form.cleaned_data for form in property_formset]
            # return redirect('pks', permanent=False)
            # return pks_search_result(request)
            return pks_search_result_sfapi(request)


@utils.log_function
def pks_search_result_sfapi(request):
    """Newer version of `pks_search_result`.
    """

    if ("_search_query_smiles" in request.session) and ("_search_query_properties" in request.session):
        request.session.set_expiry(value=0)    # user’s session cookie will expire when the user’s web browser is closed
        # search_query = request.session.get('_search_query')
        search_query_smiles = request.session['_search_query_smiles']
        search_query_properties = request.session['_search_query_properties']
        logger.info(f"search_query_smiles = {search_query_smiles}")
        logger.info(f"search_query_properties = {search_query_properties}")
        search_query = dict()
        search_query["q_uuid"] = request.session['_q_uuid']
        search_query["smiles_string"] = search_query_smiles["Q_smiles"]
        search_query["notes"] = search_query_smiles["Q_notes"]
        search_query["properties"] = search_query_properties
        logger.info(f"search_query = {search_query}")

        # Calling retrotide API
        # mol_properties = [
        #     list(constants.MOLECULE_PROPERTIES.keys())[d["Property_name"]] 
        #     for d in search_query["properties"]
        # ]
        mol_properties = [
            d["Property_name"] for d in search_query["properties"]
        ]

        # Get PKS results from retrotide and submit job using SFAPI for property prediction.
        retro_output_df = vutils.retrotide_call_new(
            smiles=search_query["smiles_string"],
            properties=mol_properties,
            clean_sulfur=True,
            insert_into_db=True,
        )
        job_id, status = vutils.sfapi_call(
            smiles_list=retro_output_df["SMILES"],
            update_query_uuid=uuid.UUID(search_query["q_uuid"]),
        )
        logger.info("Submitted the query job!")

        # Landing page
        context = {
            "message": "Query submitted!"
        }
        return render(request, "retroapp/query_submitted.html", context)

    else:
        logger.warn("404: No search query!")
        return HttpResponse("404: No search query!")


# @utils.log_function
# def pks_search_result(request):
#     """[Deprecated]
#     This function is called when the search query is submitted.
#     Tasks performed by this function:
#     1. Getting retrotide results. (Dummy values for all the properties. Only RON value is not dummy)
#     2. Apply filtering constraints on the results.
#     3. Apply specified sorting on the results.
#     4. Add resulting data to QueryResultsDB.
#     5. Render the results as HTML page to show below the query form.

#     [To be modified]
#     Above mentioned logic will be modified to incorporate SF API. 
#     The results will no longer directly be available.
#     The user will be notified that the query request has been submitted. 
#     The query status and the results can be checked on the `History` webpage.
#     """
    
#     if ("_search_query_smiles" in request.session) and ("_search_query_properties" in request.session):
#         request.session.set_expiry(value=0)    # user’s session cookie will expire when the user’s web browser is closed
#         # search_query = request.session.get('_search_query')
#         search_query_smiles = request.session['_search_query_smiles']
#         search_query_properties = request.session['_search_query_properties']
#         logger.info(f"search_query_smiles = {search_query_smiles}")
#         logger.info(f"search_query_properties = {search_query_properties}")
#         search_query = dict()
#         search_query["smiles_string"] = search_query_smiles["Q_smiles"]
#         search_query["notes"] = search_query_smiles["Q_notes"]
#         search_query["properties"] = search_query_properties
#         logger.info(f"search_query = {search_query}")

#         # Calling retrotide API
#         mol_properties = [
#             list(constants.MOLECULE_PROPERTIES.keys())[d["Property_name"]] 
#             for d in search_query["properties"]
#         ]
#         retro_df = vutils.retrotide_call(smiles=search_query["smiles_string"], properties=mol_properties)
        
#         # Update the dataframe based on min-max ranges of all properties.
#         for search_query_property in search_query["properties"]:
#             mol_property = list(constants.MOLECULE_PROPERTIES.keys())[search_query_property["Property_name"]]
#             min_val = constants.MOLECULE_PROPERTIES[mol_property]['min']
#             max_val = constants.MOLECULE_PROPERTIES[mol_property]['max']
#             search_query_property["target_val"] = None

#             property_val_range = search_query_property["property_value_range"]
#             if property_val_range:
#                 if vutils.isrange(property_val_range):    # Range is specified
#                     min_val, max_val = [float(range_val) for range_val in property_val_range.split(" - ")]
#                 else:    # target value is specified
#                     search_query_property["target_val"] = float(property_val_range)                

#             retro_df = retro_df[(retro_df[mol_property]>=min_val ) & (retro_df[mol_property]<=max_val)]
        
#         # Get all properties and their corresponding sorting modes (bool value)
#         sorting_cols = list()
#         sorting_cols_asc = list()
#         for sqprop in search_query["properties"]:
#             property_name_str = list(constants.MOLECULE_PROPERTIES.keys())[sqprop["Property_name"]]
            
#             # Target value is defined (always ascending)
#             if (sqprop["property_value_range"] != "") and (not vutils.isrange(sqprop["property_value_range"])):
#                 sorting_cols.append(property_name_str)
#                 sorting_cols_asc.append(True)
            
#             else:    # Either blank or range but not target value
#                 # Sorting mode is defined (ascending | descending)
#                 if (sqprop["Sorting_mode"] in [ASCENDING, DESCENDING]):
#                     sorting_cols.append(property_name_str)
#                     sorting_cols_asc.append(sqprop["Sorting_mode"] == ASCENDING)

#         # Target values corresponding to each property
#         prop_targets = [sqprop["target_val"] for sqprop in search_query["properties"]]

#         # function defining the key for sorting
#         def sorting_key(pd_series):
#             idx = sorting_cols.index(pd_series.name)
#             target = prop_targets[idx]
#             if target is None:
#                 return pd_series
            
#             target = [target] * len(pd_series)
#             return np.abs(np.array(target) - pd_series)

#         retro_df = retro_df.sort_values(
#             by=sorting_cols, 
#             ascending=sorting_cols_asc, 
#             key=sorting_key,
#         )

#         # To reset the index value to start from 0 again
#         retro_df = retro_df.reset_index(drop=True)

#         # TODO: map the sulfur replacement function here for all the molecules.
#         # TODO: Find an alternate more efficient implementation. This will take around 5-7 minutes for a million rows.
#         retro_df['SMILES'] = retro_df['SMILES'].map(clean_sulfur_from_smiles_string)
#         # retro_df.to_csv("Demo.csv")    # DELETE
        
#         # Write results to QueryResultsDB
#         retro_df_dict = retro_df.to_dict('records')

#         Q_uuid_fk = QueryDB.objects.latest('Timestamp')

#         model_instances = [
#             QueryResultsDB(
#                 Q_uuid_id=Q_uuid_fk.Q_uuid,
#                 SMILES=record['SMILES'],
#                 Retrotide_Similarity_SCORE=record['Retrotide_Similarity_SCORE'],
#                 DESIGN=record['DESIGN'],
#                 Cetane_number=record.get('Cetane Number', None),
#                 Research_octane_number=record.get('Research Octane Number', None),
#                 Melting_point=record.get('Melting Point', None),
#                 Flash_point=record.get('Flash Point', None),
#                 Yield_sooting_index=record.get('Yield Sooting Index', None),
#                 H1_receptor_pKd=record.get('H1 Receptor pKd', None),
#                 M2_receptor_pKd=record.get('M2 Receptor pKd', None),
#             ) for record in retro_df_dict
#         ]

#         QueryResultsDB.objects.bulk_create(model_instances)

#         # Render return webpage
#         context = {
#             "query_dict": search_query,
#             "keys": ["Rendered Molecule", *retro_df.keys()],
#             "df": retro_df,
#             "width": 243,
#         }
#         return render(request, "retroapp/show_results.html", context)

#     else:
#         logger.warn("404: No search query!")
#         return HttpResponse("404: No search query!")


class QueryHistoryView(TemplateView):
    """View for history page.
    The submitted queries and their results can be accessed from here.
    """

    template_name = "retroapp/history.html"
    

    def _results_in_database(self, query_entry: models.QuerySet) -> bool:
        """Private method to check if the results for the specified query_uuid already exist in the results database.

        Args:
            query_entry (models.QuerySet): Query for which the database results are to be checked.
        
        Returns:
            bool: `True` if distinct values of the first query property contains None, `False` otherwise.
        """

        # Fetching the first query property for the specified query
        qprop = QueryPropertyDB.objects.filter(Q_uuid_id=query_entry.Q_uuid)
        if qprop is not None:
            # qprop_name = list(constants.MOLECULE_PROPERTIES.keys())[qprop[0].Property_name]
            qprop_name = constants.MOLECULE_PROPERTIES[qprop[0].Property_name]
        
        # Checking the result table values of the selected query property
        qres = QueryResultsDB.objects.filter(Q_uuid_id=query_entry.Q_uuid)
        qres_vals = qres.values_list(constants.DB_COLUMN_NAMES[qprop_name])
        distinct_qres_vals = qres_vals.distinct()

        # Exception when no results are stored against a query. Logically this shouldn't be happening. Hence returning True to avoid any changes because of this.
        if distinct_qres_vals.count() <= 0:
            logging.warning("This should not be happening!!! Ignoring this for now...")
            return True

        # Result does not exist in the database if property values contain None
        if distinct_qres_vals[0][0] is None:
            return False
        
        return True


    @utils.log_function
    def _update_completed_queries_into_db(
        self, 
        property_predictor_obj=None,
    ) -> None:
        """Private method to update the results table (QueryResultsDB) in the database for all the completed query submissions.

        Args:
            property_predictor_obj (PropertyPredictor, optional): Object to use the SFAPI. Defaults to None.
        """

        # Get all the completed queries
        query_db_completed = QueryDB.objects.filter(Q_Status="COMPLETED")
        print("query_db_completed.count() = ", query_db_completed.count())

        if property_predictor_obj is None:
            property_predictor_obj = vutils.get_PropertyPredictor_obj()
        property_predictor_obj.open_session()

        for cquery_i in query_db_completed:
            # Continue if the query already has results in the database
            if self._results_in_database(cquery_i):
                continue
            
            # Get results from the SFAPI
            # TODO: Test on Spin NERSC.
            pp_qresults = property_predictor_obj.get_query_results(cquery_i.Q_Job_id)
            # pp_qresults = {
            #     "CN": [vutils.predict_property_api("CN") for _ in range(100)],
            #     "MP": [vutils.predict_property_api("MP") for _ in range(100)],
            # }    # DELETE

            # Update database
            existing_qres_db = QueryResultsDB.objects.filter(Q_uuid_id=cquery_i.Q_uuid)

            update_list = []
            for row_i, model_obj in enumerate(existing_qres_db):
                for pkey, pvalues in pp_qresults.items():
                    property_col_name = constants.DB_COLUMN_NAMES.get(pkey, None)

                    if property_col_name:
                        try:
                            model_obj.__dict__[property_col_name] = pvalues[row_i]
                        except Exception as e:
                            raise f"[ERROR] Unable to update values in the dictionary. {e}"

                update_list.append(model_obj)
                
            QueryResultsDB.objects.bulk_update(
                objs=update_list, 
                fields=[constants.DB_COLUMN_NAMES[pkey_i] for pkey_i in pp_qresults.keys() if pkey_i in constants.DB_COLUMN_NAMES.keys()],
                batch_size=100,
            )


    @utils.log_function
    def setup(self, request, *args, **kwargs):
        super().setup(request, *args, **kwargs)

        # Update Databases only if the user is authenticated.
        if self.request.user.is_authenticated:
            # Update QueryDB status for each query entry.
            # Exclude all the rows that do not have a Q_Job_id or have 'COMPLETED' or 'FAILED' Q_Status.
            query_db_unfinished = QueryDB.objects.exclude(Q_Status__in=["COMPLETED", "FAILED"]).exclude(Q_Job_id="")

            # # DELETE
            # def get_dummy_status(del_jobid):
            #     if del_jobid=="-1":
            #         return "NEG"
            #     elif del_jobid=="":
            #         return "EMPTY"
            #     else:
            #         return "COMPLETED"
            
            # Open PropertyPredictor session to fetch status
            pp = vutils.get_PropertyPredictor_obj()
            pp.open_session()
            # pp = None    # DELETE
            for query_i in query_db_unfinished:
                query_i.Q_Status = pp.job_status(query_i.Q_Job_id)
                # query_i.Q_Status = get_dummy_status(query_i.Q_Job_id)    # DELETE
                query_i.save()     # Writing updated value to DB
            
            # Update QueryResultsDB if status is COMPLETED
            # TODO: Test on Spin NERSC
            self._update_completed_queries_into_db(property_predictor_obj=pp)


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

        # context['MOLECULE_PROPERTIES'] = list(constants.MOLECULE_PROPERTIES.keys())
        context['MOLECULE_PROPERTIES'] = constants.MOLECULE_PROPERTIES
        # context['SORTING_OPTIONS'] = [s_optn[1] for s_optn in constants.SORTING_OPTIONS]
        context['SORTING_OPTIONS'] = constants.SORT_OPTIONS

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

    def _query_user_constraints(self, db_object) -> pd.DataFrame:
        """Queries all the user defined constraints on the results received from the property prediction.

        Args:
            db_object (QuerySet): The results from the property prediction being fetched directly from the database.

        Returns:
            pd.DataFrame: Resulting DataFrame after applying all the filters and sorting specified by the user in the query.
        """

        # get all the constraints
        query_uuid = self.request.GET['name']    # will have to fetch constraints using this UUID
        
        query_property_constaints = QueryPropertyDB.objects.filter(Q_uuid_id=query_uuid)

        db_object_temp = copy.copy(db_object)

        # Apply constraints to the database before viewing in the webpage
        for constraint_i in query_property_constaints:
            mol_property = constants.MOLECULE_PROPERTIES[constraint_i.Property_name].label
            mol_property = constants.DB_COLUMN_NAMES[mol_property]

            if constraint_i.Min_value:
                db_object_temp = db_object_temp.filter(**{
                    f"{mol_property}__gte": constraint_i.Min_value
                })

            if constraint_i.Max_value:
                db_object_temp = db_object_temp.filter(**{
                    f"{mol_property}__lte": constraint_i.Max_value
                })
        
        db_obj_df = pd.DataFrame(db_object_temp.values())

        # Get all properties and their corresponding sorting modes (bool value)
        sorting_cols = list()
        sorting_cols_asc = list()
        for constraint_i in query_property_constaints:
            mol_property = constants.MOLECULE_PROPERTIES[constraint_i.Property_name].label
            mol_property = constants.DB_COLUMN_NAMES[mol_property]

            if constraint_i.Target_value:
                # Target value is defined (always ascending)
                sorting_cols.append(mol_property)
                sorting_cols_asc.append(True)
            
            else:    
                # Either blank or range but not target value
                # Sorting mode is defined (ascending | descending)
                if (constraint_i.Sorting_mode in [
                    constants.SORT_OPTIONS.ASCENDING, 
                    constants.SORT_OPTIONS.DESCENDING
                ]):
                    sorting_cols.append(mol_property)
                    sorting_cols_asc.append(constraint_i.Sorting_mode == constants.SORT_OPTIONS.ASCENDING)
        
        # Get Target values for each property
        prop_targets = [constraint_i.Target_value for constraint_i in query_property_constaints]

        # function defining the key for sorting
        def sorting_key(pd_series):
            idx = sorting_cols.index(pd_series.name)
            target = prop_targets[idx]
            if target is None:
                return pd_series
            
            target = [target] * len(pd_series)
            return np.abs(np.array(target) - pd_series)

        db_obj_df = db_obj_df.sort_values(
            by=sorting_cols, 
            ascending=sorting_cols_asc, 
            key=sorting_key,
        )

        # To reset the index value to start from 0 again
        db_obj_df = db_obj_df.reset_index(drop=True)

        return db_obj_df


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

        query_res_df = self._query_user_constraints(query_res_object)

        context["query_res_object"] = query_res_df

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
    retro_df = vutils.retrotide_call(smiles=smiles)

    context = {
        "query_dict": {"smiles_string": smiles},
        "keys": ["Rendered Molecule", *retro_df.keys()],
        "df": retro_df,
        "width": width,
    }
    return render(request, "retroapp/showtable.html", context)
