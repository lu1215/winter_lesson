from django.urls import path
from enrichment_app import views

urlpatterns = [
    path("",views.enrichment_page),
    path("enrichment_ajax/",views.enrichment_ajax),
]
