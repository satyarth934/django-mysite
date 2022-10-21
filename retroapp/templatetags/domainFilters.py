from django import template
# from django.utils.http import urlquote
from urllib.parse import quote as urlquote

from rdkit import Chem
import re
import numbers


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
    if isinstance(inputFloat, float):
        return '%s' % float('%.3g' % inputFloat)
    return inputFloat


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


@register.filter
def isnumeric(value):
    return isinstance(value, numbers.Number)


@register.filter
def isfloat(value):
    return isinstance(value, float)


@register.filter
def islist(value):
    return isinstance(value, list)

@register.filter
def index(indexable, i):
    """Get value at index 'i' in 'indexable' iterable.

    Args:
        indexable (Iterable): Iterable to get the inde value from.
        i (int): Index position.

    Returns:
        Object: Value at index 'i'.

    Example usage:
        {{ iterable|index:i }}
    """
    if i >= len(indexable):
        return f"ERROR: Index out of range. index={i}; len(list)={len(indexable)}"
    return indexable[i]


@register.filter
def fetch_choice_label(indexable, key):
    if key in indexable:
        return indexable[key].label
    return f"ERROR: {key} does not exist in {indexable}."