from email.policy import default
import uuid
from retroapp import constants
from django.db import models
from django.utils.translation import gettext_lazy as _


# Create your models here.
class QueryDB(models.Model):
    Q_uuid = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    Timestamp = models.DateTimeField(auto_now_add=True)
    Username = models.CharField(max_length=100, blank=True, null=True)
    Full_name = models.CharField(max_length=100, blank=True, null=True)
    Email = models.CharField(max_length=100, blank=True, null=True)
    Q_smiles = models.CharField(max_length=2000, blank=True, null=True)
    Q_notes = models.CharField(max_length=2000, blank=True, null=True)
    Q_Task_id = models.CharField(max_length=200, blank=True, null=True)
    Q_Job_id = models.CharField(max_length=200, blank=True, null=True)
    Q_Status = models.CharField(max_length=50, blank=True, null=True)    # Choices: SUCCESS | FAILED | SUBMITTED - to be pulled from slurm status API.



class QueryPropertyDB(models.Model):
    id = models.BigAutoField(primary_key=True, editable=False)
    # id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    Q_uuid = models.ForeignKey(
        QueryDB,
        on_delete=models.CASCADE,
    )
    # Property_name = models.IntegerField(
    #     # choices=enumerate(constants.MOLECULE_PROPERTIES.keys()),
    #     # default=list(constants.MOLECULE_PROPERTIES.keys()).index('Cetane Number'),
    #     choices=tuple(constants.PROPERTY_CODES.items()),
    #     default="CN",
    #     blank=True,
    #     null=True,
    # )

    Property_name = models.CharField(
        max_length=3,
        choices=constants.MOLECULE_PROPERTIES.choices,
        default=constants.MOLECULE_PROPERTIES.CN,
    )
    Target_value = models.FloatField(blank=True, null=True)
    Min_value = models.FloatField(blank=True, null=True)
    Max_value = models.FloatField(blank=True, null=True)
    
    # Sorting_mode = models.IntegerField(
    #     choices=constants.SORTING_OPTIONS,
    #     default=constants.NO_SORT,
    #     blank=True,
    #     null=True,
    # )

    Sorting_mode = models.CharField(
        max_length=10,
        choices=constants.SORT_OPTIONS.choices,
        default=constants.SORT_OPTIONS.DESCENDING,
        blank=True,
        null=True,
    )


class QueryResultsDB(models.Model):
    id = models.BigAutoField(primary_key=True, editable=False)
    # id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    Q_uuid = models.ForeignKey(
        QueryDB,
        on_delete=models.CASCADE,
    )
    # ResultBlob    # Expand the table with multiple Q_uuid
    SMILES = models.CharField(max_length=2000, blank=True, null=True)
    Retrotide_Similarity_SCORE = models.FloatField(blank=True, null=True)
    DESIGN = models.CharField(max_length=2000, blank=True, null=True)
    Cetane_number = models.FloatField(blank=True, null=True)
    Research_octane_number = models.FloatField(blank=True, null=True)
    Melting_point = models.FloatField(blank=True, null=True)
    Flash_point = models.FloatField(blank=True, null=True)
    Yield_sooting_index = models.FloatField(blank=True, null=True)
    H1_receptor_pKd = models.FloatField(blank=True, null=True)
    M2_receptor_pKd = models.FloatField(blank=True, null=True)
