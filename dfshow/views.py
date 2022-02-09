from multiprocessing import context
from re import template
from django.shortcuts import render

# Create your views here.
from django.http import HttpResponse, JsonResponse
import pandas as pd
from django.template import loader


table = {
        "pk": [1,2,3],
        "col1": ["val1", "val2", "val3"],
        "col2": ["v1", "v2", "v3"],
        "col3": ["value1", "value2", "value3"],
    }
df = pd.DataFrame(table)


def index(request):
    return HttpResponse("dfshow -> view_index")


# http://127.0.0.1:8000/dfshow/showtable/
def showtable(request):
    template = loader.get_template("dfshow/showtable.html")
    context = {
        "keys": df.keys(),
        "table_rows": [row[1] for row in df.iterrows()],
        }
    rendered_str = template.render(context, request)
    return HttpResponse(rendered_str)


# http://127.0.0.1:8000/dfshow/showtablejson/
def showtablejson(request):
    return JsonResponse(df.to_json(), safe=False)
