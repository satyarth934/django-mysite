from contextlib import contextmanager
from tempfile import tempdir
from tkinter.messagebox import QUESTION
from django.shortcuts import render

# Create your views here.
from django.http import HttpResponse, Http404
from .models import Question
from django.template import loader

def index(request):
    latest_question_list = Question.objects.order_by('-pub_date')[:5]
    template = loader.get_template("polls/index.html")
    context = {
        "latest_question_list": latest_question_list,
    }
    # return HttpResponse("Hello, world. You're at the polls index.")
    return HttpResponse(template.render(context, request))


def lolpoint(request):
    return HttpResponse("LolPoint!")


def detail(request, question_id):
    try:
        question = Question.objects.get(pk=question_id)
    except Question.DoesNotExist:
        raise Http404("Question does not exist.")

    template = loader.get_template("polls/detail.html")
    context = {
        "question": question,
    }
    # return HttpResponse("You're looking at question %s." % question_id)
    return HttpResponse(template.render(context, request))


def results(request, question_id):
    response = "You're looking at the results of question %s."
    return HttpResponse(response % question_id)


def vote(request, question_id):
    return HttpResponse("You're voting on question %s." % question_id)