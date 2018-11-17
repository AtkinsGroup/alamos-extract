import re
from collections import OrderedDict
import warnings

import bs4
import requests
import pandas as pd

__author__ = "Stephen Gaffney"
__copyright__ = "Stephen Gaffney"
__license__ = "gpl3"


def load_cluster(cluster_id):
    """Obtains cluster name, description, patient names+ids, and accession names+ids.

    The IDs are required for URLs.

    Args:
        cluster_id (int): a cluster ID (as in URL id)

    Returns:
        data (dict): holds cluster_str, desc, patients, accessions.
    """
    cluster_path = 'https://www.hiv.lanl.gov/components/sequence/HIV/search/cluster.comp?clu_id={}'

    # @TODO: find a way to get around "certificate verify failed" error without verify=False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c_page = requests.get(cluster_path.format(cluster_id), verify=False).content

    cluster_name_str = 'Cluster Name'
    patient_url_format = r"patient.comp\?pat_id=(\d+)"  # includes backslash escape
    genbank_url_format = r"query_one.comp\?se_id=(\d+)"

    col_dict = {'cluster_str': 'Cluster Name',
                'desc': 'Cluster Description',
                'patients': 'Patient(s)',
                'accessions': 'Accession(s)',
                }

    soup = bs4.BeautifulSoup(c_page, features="lxml", from_encoding='utf8')
    tables = soup('table')

    # Get main table
    main_table = [i for i in tables if cluster_name_str in i.text]
    assert(len(main_table) == 1), 'Failed to find main table: multiple have cluster name field'
    main_table = main_table[0]

    df = pd.read_html(str(main_table))[0]
    df = df.set_index([0]).transpose()
    cluster_str = df[col_dict['cluster_str']].iloc[0]
    desc = df[col_dict['desc']].iloc[0].strip(""" '" """)

    # GET PATIENT IDS
    patients = OrderedDict()
    patient_links = soup.find_all(href=re.compile(patient_url_format))
    for i in patient_links:
        patient_name = i.text.strip()
        # href is of form 'patient.comp?pat_id=38516'
        href = i['href']
        patient_id = int(re.findall(patient_url_format, href)[0])
        patients[patient_name] = patient_id

    # GET SEQUENCE IDS
    seqs = OrderedDict()
    seq_links = soup.find_all(href=re.compile(genbank_url_format))
    for i in seq_links:
        seq_name = i.text.strip()
        # href is of form 'patient.comp?pat_id=38516'
        href = i['href']
        seq_id = int(re.findall(genbank_url_format, href)[0])
        seqs[seq_name] = seq_id

    data = {
        'cluster_name': cluster_str,
        'description': desc,
        'patients': patients,
        'accessions': seqs
    }

    return data
