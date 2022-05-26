from collections import defaultdict
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

import retrotide


# Create your views here.
def index(request):
    return HttpResponse("You're at retroapp. Append smile string in the URL to invoke the retrotide API.")


def retrotideAPI_dummy(request, smiles):
    example_df = pd.read_csv("retroapp/example_df.tab", delimiter="\t")
    return example_df


def retrotide_call(smiles):
    designs = retrotide.designPKS(Chem.MolFromSmiles(smiles))

    df_dict = defaultdict(list)

    for i in range(len(designs[-1])):
        df_dict["SNo"].append(i+1)
        df_dict["SMILES"].append(Chem.MolToSmiles(designs[-1][i][2]))
        df_dict["SCORE"].append(designs[-1][i][1])
        df_dict["DESIGN"].append(designs[-1][i][0].modules)

    return pd.DataFrame(df_dict)


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
    # retro_df = retrotideAPI_dummy(request, smiles)
    retro_df = retrotide_call(smiles)
    print(retro_df)         # DELETE: only for debugging purposes

    # show table for the return data frame from the API call
    table_rendered_str = showtable(request, retro_df, width)
    
    return HttpResponse(table_rendered_str)