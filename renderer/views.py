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


# Create your views here.
def index(request):
    return HttpResponse("You're at renderer. Append smile string in the URL to render the smile string as a molecule.")


def renderer_dummy(request, smiles, check, width=243):
    return HttpResponse(f"Running the render function on the input => {smiles}.\nChecking 2nd argument => {check}")

@cache_page(60 * 60 * 24 * 7) # cache for one week
def renderer(request, smiles, width=243):
    # this function renders an input SMILES chemical structure as an SVG
    # and snaps it to a template if possible, for stable orientation

    # parse smiles input
    smiles = urlunquote(smiles)
    smiles = re.match(r'^(\S{1,10000})', str(smiles)).group(1)
    mol = Chem.MolFromSmiles(smiles)
    Chem.SanitizeMol(mol)
    width=int(width)
    assert 0 < width < 2000

    template = Chem.MolFromSmiles('C(=O)[S]')
    highlightAtomLists = None
    legend = ''

    # lock molecule to template
    AllChem.Compute2DCoords(template)
    height = int(width * ((30.0 + mol.GetNumAtoms() * 7.0) / 243.0))

    align = True 
    if 'align' in request.GET:
        if request.GET['align'] == 'False':
            align = False

    if align:
        try:
            AllChem.GenerateDepictionMatching2DStructure(mol, template, acceptFailure=False)
        except ValueError:
            template = Chem.MolFromSmiles('C(=O)O')
            AllChem.Compute2DCoords(template)
            try:
                AllChem.GenerateDepictionMatching2DStructure(mol, template, acceptFailure=False)
            except ValueError:
                height = width 
    else:
        height = width

    # draw SVG
    svg = Draw.MolsToGridImage([mol], highlightAtomLists=[highlightAtomLists],
                 molsPerRow=1,
                 subImgSize=(width, height), useSVG=True)

    # fix XML namespace in RDKit export
    xmlTree = ET.fromstring(svg) 
    xmlTree.attrib['xmlns'] = 'http://www.w3.org/2000/svg'
    xmlTree.attrib['xmlns:svg'] = 'http://www.w3.org/2000/svg'

    # convert XML to string
    fixedSVG = BytesIO()
    ET.ElementTree(xmlTree).write(fixedSVG, encoding='utf-8', xml_declaration=True)

    return HttpResponse(fixedSVG.getvalue(), content_type='image/svg+xml')