from django.db import models

# Create your models here.

class user_parameters(models.Model):
    user=models.CharField(max_length=200)
    name=models.CharField(max_length=200)
    filt=models.CharField(max_length=10000)
    param=models.CharField(max_length=1000)