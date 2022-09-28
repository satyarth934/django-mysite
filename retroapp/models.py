import uuid
from retroapp import constants
from django.db import models


# Create your models here.
class QueryDB(models.Model):
    Q_uuid = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    Timestamp = models.DateTimeField(auto_now_add=True)
    Username = models.CharField(max_length=100)
    Full_name = models.CharField(max_length=100)
    Email = models.CharField(max_length=100)
    Q_smiles = models.CharField(max_length=2000)
    Q_notes = models.CharField(max_length=2000, blank=True)
    Q_Job_id = models.CharField(max_length=200)
    Q_Status = models.CharField(max_length=50, blank=True)    # Choices: SUCCESS | FAILED | SUBMITTED - to be pulled form slurm status API.



class QueryPropertyDB(models.Model):
    Q_uuid = models.ForeignKey(
        'QueryDB',
        on_delete=models.CASCADE,
    )
    Property_name = models.IntegerField(
        choices=enumerate(constants.MOLECULE_PROPERTIES.keys()),
        default=list(constants.MOLECULE_PROPERTIES.keys()).index('Cetane Number'),
    )
    Target_value = models.FloatField(blank=True, null=True)
    Min_value = models.FloatField(blank=True, null=True)
    Max_value = models.FloatField(blank=True, null=True)
    Sorting_mode = models.IntegerField(
        choices=constants.SORTING_OPTIONS,
        default=constants.NO_SORT,
        blank=True,
        null=True,
    )


class QueryResultsDB(models.Model):
    Q_uuid = models.ForeignKey(
        'QueryDB',
        on_delete=models.CASCADE,
    )
    # ResultBlob    # Expand the table with multiple Q_uuid
    SMILES = models.CharField(max_length=2000)
    Retrotide_Similarity_SCORE = models.FloatField(blank=True, null=True)
    DESIGN = models.CharField(max_length=2000, blank=True, null=True)
    Cetane_number = models.FloatField(blank=True, null=True)
    Research_octane_number = models.FloatField(blank=True, null=True)
    Melting_point = models.FloatField(blank=True, null=True)
    Flash_point = models.FloatField(blank=True, null=True)
    Yield_sooting_index = models.FloatField(blank=True, null=True)
    H1_receptor_pKd = models.FloatField(blank=True, null=True)
    M2_receptor_pKd = models.FloatField(blank=True, null=True)
