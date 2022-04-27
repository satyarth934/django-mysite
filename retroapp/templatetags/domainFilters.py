from django import template
# from django.utils.http import urlquote
from urllib.parse import quote as urlquote
from rdkit import Chem
import re


register = template.Library()


@register.filter
def classname(obj):
    return obj.__class__.__name__


@register.filter
def smiles(mol):
    if mol:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        return urlquote(smiles)
    else:
        return False 


@register.filter
def unquotedSmiles(mol):
    return Chem.MolToSmiles(mol, isomericSmiles=True)


@register.filter
def stripTrailingVersion(accession):
    return re.sub("\.\d+$", "", accession)


@register.filter
def urlq(str):
    return urlquote(str)


@register.filter
def sigfig(inputFloat):
    return '%s' % float('%.3g' % inputFloat)


@register.filter
def replace(value, arg):
    """
    Replacing filter
    Use `{{ "aaa"|replace:"a|b" }}`
    """
    if len(arg.split('|')) != 2:
        return value

    what, to = arg.split('|')
    return value.replace(what, to)