from django.shortcuts import render
from web.winter_project.enrichment_app.enrichment_backup import enrichment
from django.http import JsonResponse
# Create your views here.

def enrichment_page(request):
    return render(request, 'enrichment.html')

def enrichment_ajax(request):
    seq_data = request.POST["seq"]
    correction = request.POST["Correction"]
    p_limit = float(request.POST["p_limit"])
    enrichment_result = enrichment(seq_data, correction, p_limit)["result"]
    return JsonResponse({"enrichment_result":enrichment_result})