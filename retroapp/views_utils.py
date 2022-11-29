###########
# Imports #
###########
from multiprocessing.sharedctypes import Value
from typing import List, Tuple, Union
import warnings
from mysite import utils
warnings.formatwarning = utils.warning_format    # Defining a specific format to print warnings on console

import logging
logger = logging.getLogger(f"{__name__}")

import numpy as np
import pandas as pd
from collections import defaultdict

import uuid
from retroapp.models import QueryDB, QueryPropertyDB, QueryResultsDB
from retroapp import utils
from django.conf import settings

from backend import prop_sfapi as prop

import retrotide
from rdkit import Chem


@utils.log_function
def is_user_authenticated(request):
    """Checks if the user is authenticated using the google oauth API and also checks if the user is from the WHITELISTED domains.

    Args:
        request (django.http.request): request object

    Returns:
        bool: True if the user is authenticated and from the accepted domain, False otherwise.
    """

    num_tries = 5
    for i in range(num_tries):
        try:
            login_authentication = request.user.is_authenticated

            if not login_authentication:
                return login_authentication

            user_email = request.user.socialaccount_set.all()[0].extra_data['email']
            domain_authentication = user_email.split("@")[-1] in settings.SOCIAL_AUTH_GOOGLE_OAUTH2_WHITELISTED_DOMAINS

            return (login_authentication and domain_authentication)
        except django.db.utils.OperationalError as dj_op:
            logger.warn(f"Attempting again ({i}/{num_tries}) because of the following error: {dj_oe}")

    logger.error(f"Failed even after {num_tries} attempts. {dj_oe}")


@utils.log_function
def insert_data_into_db(data_df: pd.DataFrame) -> None:
    """This function is used to insert the retrotide output to the results table (QueryResultsDB). The property values are predicted after this insertion.

    Args:
        data_df (pd.DataFrame): Input DataFrame that is to be inserted into the database.
    """
    # Write results to QueryResultsDB
    data_df_dict = data_df.to_dict('records')

    Q_uuid_fk = QueryDB.objects.latest('Timestamp')

    model_instances = [
        QueryResultsDB(
            Q_uuid_id=Q_uuid_fk.Q_uuid,
            SMILES=record['SMILES'],
            Retrotide_Similarity_SCORE=record['Retrotide_Similarity_SCORE'],
            DESIGN=record['DESIGN'],
            Cetane_number=record.get('Cetane Number', None),
            Research_octane_number=record.get('Research Octane Number', None),
            Melting_point=record.get('Melting Point', None),
            Flash_point=record.get('Flash Point', None),
            Yield_sooting_index=record.get('Yield Sooting Index', None),
            H1_receptor_pKd=record.get('H1 Receptor pKd', None),
            M2_receptor_pKd=record.get('M2 Receptor pKd', None),
        ) for record in data_df_dict
    ]

    QueryResultsDB.objects.bulk_create(model_instances)


@utils.log_function
def retrotide_call_new(
    smiles: str, 
    properties: Union[str, List]=None,
    clean_sulfur: bool=True,
    insert_into_db: bool=True,
) -> pd.DataFrame:
    """Newer version of `retrotide_call`.
    The PKS Design is computed using the retrotide API.
    The sulfur molecules are removed from the PKS Design SMILES strings if the `clean_sulfur` argument is True.
    The data is stored into the database if the `insert_into_db` argument is True.

    Args:
        smiles (str): Input SMILES string.
        properties (Union[str, List], optional): List of properties which are to be estimated using GraphDot. Defaults to None.
        clean_sulfur (bool): Boolean flag to remove the sulfur molecules from the PKS Design SMILES strings.
        insert_into_db (bool): Boolean flag to insert the PKS Design output into the database.

    Returns:
        pd.DataFrame: Output result from retrotide and GraphDot returned as a DataFrame.
    """

    designs = retrotide.designPKS(Chem.MolFromSmiles(smiles))
    if properties is not None:    # Convert to a list if properties is not None
        properties = [properties] if isinstance(properties, str) else properties
        
    else:
        properties = []
        logger.warn(f"No properties mentioned.")

    # Store retrotide results in a dict.
    df_dict = defaultdict(list)
    for i in range(len(designs[-1])):
        df_dict["SMILES"].append(Chem.MolToSmiles(designs[-1][i][2]))
        df_dict["Retrotide_Similarity_SCORE"].append(designs[-1][i][1])
        df_dict["DESIGN"].append(designs[-1][i][0].modules)
    
    try:
        output_df = pd.DataFrame(df_dict)
    except ValueError as ve:
        raise ValueError(f"{ve} - Check if multiple constraints are defined for the same Molecular Property.")

    if clean_sulfur:
        output_df['SMILES'] = output_df['SMILES'].map(clean_sulfur_from_smiles_string)
    
    if insert_into_db:
        insert_data_into_db(data_df=output_df)
        logger.info("[DB UPDATE] Added retrotide PKS Designs into the results database!")

    return output_df


@utils.log_function
def get_PropertyPredictor_obj():
    debug = 0 # produces minimal output
    # debug = 1 # produces more output
    # debug = 2 # produces a lot of output

    # client_sys = "spin"
    target_job = "perl_test"
    system = "perlmutter"

    # path to where batch job runs
    # path = "/global/cfs/cdirs/m3513/molinv/rev5"
    path = "/global/cfs/cdirs/m3513/molinv/prod"
    
    client_id = "BIOARC_SPIN_SFAPI_CLIENT" # this environment variable = SFAPI client id
    client_pem = "BIOARC_SPIN_SFAPI_PEM_KEY" # this environment variable = private PEM key
    client_env = True

    pp = prop.PropertyPredictor(
        path, 
        client_id, 
        client_pem, 
        client_env, 
        target_job, 
        debug, 
    )

    # Check if PropertyPredictor was initialized correctly.
    pp.open_session()
    status = pp.check_status()
    logger.info(msg=f"System: {system}, status: {status}")
    if status=='unavailable':
        raise Exception(f"Target system: {system} status='unavailable', exiting")

    return pp


@utils.log_function
def sfapi_call(
    smiles_list: List[str],
    properties: Union[str, List]=None,
    update_query_uuid: Union[str, uuid.UUID]=None,
) -> Tuple:
    """The properties are currently being estimated using GraphDot. This might change in the future.

    Args:
        smiles_list (List[str]): List of SMILES strings.
        properties (Union[str, List], optional): List of properties which are to be estimated using GraphDot. Defaults to None.
        update_query_uuid (bool): Query UUID against which the results are to be stored. Defaults to None.

    Returns:
        Tuple: A tuple of job id and the submitted job status.
    """

    if not isinstance(smiles_list, list):
        smiles_list = list(smiles_list)

    if properties is not None:    # Convert to a list if properties is not None
        properties = [properties] if isinstance(properties, str) else properties
        
    else:
        properties = []
        logger.warn(f"No properties mentioned.")

    # TODO: Make SFAPI call here. 
    pp = get_PropertyPredictor_obj()
    pp.open_session()
    props = {property: None for property in properties}
    job_id = pp.submit_query(smiles_list, props)
    if job_id == "ERROR":
        status = "ERROR"
    else:
        status = pp.job_status(job_id)

    # # TODO: DELETE these dummy values and use the above code.
    # job_id = 12345
    # status = "PENDING"

    if update_query_uuid is not None:
        if isinstance(update_query_uuid, str):
            update_query_uuid = uuid.UUID(update_query_uuid)
        elif isinstance(update_query_uuid, uuid.UUID):
            pass
        else:
            raise ValueError("Invalid `update_query_uuid`. Must be a string or a uuid.UUID.")

        # Update the Job ID and Status for the submitted query in the database.
        qdb = QueryDB.objects.get(Q_uuid=update_query_uuid)
        qdb.Q_Job_id = job_id  # change Q_Job_id
        qdb.Q_Status = status  # change Q_Status
        qdb.save() # this will update the row in the QueryDB database
        logger.info("[DB UPDATE] Updated job_id and status in the database!")

    return job_id, status


def isrange(property_value_str):
    """Checks whether the specified property is a range or not.

    Args:
        property_value_str (str): Input property value string from the form.

    Returns:
        bool: True if the property is a range, False otherwise.
    """

    return " - " in property_value_str


def clean_sulfur_from_smiles_string(input_smiles, idx=0):
    """Cleans sulfur molecule from the SMILES string.
    """

    try:
        input_mol = Chem.MolFromSmiles(input_smiles)
        Chem.SanitizeMol(input_mol)
        rxn = Chem.AllChem.ReactionFromSmarts('[C:1](=[O:2])[S:3]>>[C:1](=[O:2])[O].[S:3]')
        product = rxn.RunReactants((input_mol,))[idx][0]
        Chem.SanitizeMol(product)
    except Exception as e:
        print("input_smiles:", input_smiles)
        raise e

    return Chem.MolToSmiles(product)











##########################
# DELETE THESE FUNCTIONS #
##########################

def predict_property_api(property):
    """[Dummy] # DELETE
    Dummy function call to replicate the molecule property prediction.
    Returns random values for the query properties requested.
    """

    warnings.warn(f"Populating {property} with RANDOM values.")    # DELETE
    return np.random.randint(0,100)


@utils.log_function
def retrotide_call(smiles, properties=None):
    """[Deprecated]
    Retrotide API call.

    Newer version is implemented in `retrotide_call_new`.
    """

    designs = retrotide.designPKS(Chem.MolFromSmiles(smiles))

    if properties is not None:    # Convert to a list if properties is not None
        properties = [properties] if isinstance(properties, str) else properties
        
    else:
        properties = []
        logger.warn(f"No properties mentioned.")

    df_dict = defaultdict(list)

    for i in range(len(designs[-1])):
        df_dict["SMILES"].append(Chem.MolToSmiles(designs[-1][i][2]))
        df_dict["Retrotide_Similarity_SCORE"].append(designs[-1][i][1])
        df_dict["DESIGN"].append(designs[-1][i][0].modules)
        for property in properties:
            df_dict[property].append(predict_property_api(property=property))

    try:
        output_df = pd.DataFrame(df_dict)
    except ValueError as ve:
        raise ValueError(f"{ve} - Check if multiple constraints are defined for the same Molecular Property.")

    return pd.DataFrame(df_dict)
