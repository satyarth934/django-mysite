from django.shortcuts import render

from django.http import Http404, HttpResponse
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


# Create your views here.
def index(request):
    return HttpResponse("You're at retroapp. Append smile string in the URL to invoke the retrotide API.")


def retrotideAPI_dummy(request, smiles):
    # smiles_strs = ['CCCO','C=CCC(CC)CCC','COCC=O']
    # sim_metric_vals = [12.34, 23.45, 34.56]
    # big_ugly_strs = ['lol1', 'lol2', 'lol3']

    # retro_df = {
    #     "rendered_mol": smiles_strs,
    #     "smiles_str": smiles_strs,
    #     "sim_metric": sim_metric_vals,
    #     "big_ugly_str": big_ugly_strs,
    # }
    # return pd.DataFrame(retro_df)

    example_df = pd.read_csv("retroapp/example_df.tab", delimiter="\t")

    return example_df


def showtable(request, retro_df, width=243):
    template = loader.get_template("retroapp/showtable.html")
    context = {
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
    retro_df = retrotideAPI_dummy(request, smiles)
    print(retro_df)         # DELETE: only for debugging purposes

    # show table for the return data frame from the API call
    table_rendered_str = showtable(request, retro_df, width)
    
    return HttpResponse(table_rendered_str)