import uuid
from retroapp import constants
from django.db import models


# Create your models here.
class QueryDB(models.Model):
    Q_uuid = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    Username = models.CharField(max_length=100)
    Full_name = models.CharField(max_length=100)
    Email = models.CharField(max_length=100)
    # Q_smiles = models.TextField(max_length=2000)
    Q_smiles = models.CharField(max_length=2000)
    # Q_notes = models.TextField()
    Q_notes = models.CharField(max_length=2000, blank=True)
    Job_id = models.CharField(max_length=200)


class QueryPropertyDB(models.Model):
    Q_uuid = models.ForeignKey(
        'QueryDB',
        on_delete=models.CASCADE,
    )
    Property_name = models.CharField(
        max_length=200, 
        choices=enumerate(constants.MOLECULE_PROPERTIES.keys()),
        default=list(constants.MOLECULE_PROPERTIES.keys()).index('Cetane Number'),
    )
    Target_value = models.FloatField()
    Min_val = models.FloatField()
    Max_val = models.FloatField()
    Sorting_mode = models.CharField(
        max_length=200, 
        choices=enumerate(constants.SORTING_OPTIONS),
        default=constants.SORTING_OPTIONS.index(constants.DESCENDING),
    )

# TODO
# class QueryResultsDB(models.Model):
#     Q_uuid = models.ForeignKey(
#         'QueryDB',
#         on_delete=models.CASCADE,
#     )
#     # ResultBlob    # Expand the table with multiple Q_uuid
#     SMILES_string = models.CharField(max_length=2000)
#     Retrotide_Similarity_SCORE = models.FloatField()
#     Design = models.CharField(max_length=2000)
#     Cetane_number = models.FloatField()
#     Research_octane_number = models.FloatField()
#     Melting_point = models.FloatField()
#     Flash_point = models.FloatField()
#     Yield_sooting_index = models.FloatField()
#     H1_receptor_pKd = models.FloatField()
#     M2_receptor_pKd = models.FloatField()

    