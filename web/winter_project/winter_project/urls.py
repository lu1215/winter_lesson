"""
URL configuration for winter_project project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include

urlpatterns = [
    path('admin/', admin.site.urls),
    path('gene_search/', include('gene_search_app.urls')),
    path('survival_analysis/', include('survival_analysis_app.urls')),
    path('DE_analysis_app/', include('DE_analysis_app.urls')),
    path('enrichment_app/', include('enrichment_app.urls')),
    path('miRNA_screener_app/', include('miRNA_screener_app.urls')),
]
