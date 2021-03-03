#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is a skeleton file that can serve as a starting point for a Python
console script. To run this script uncomment the following lines in the
[options.entry_points] section in setup.cfg:

    console_scripts =
         fibonacci = alamos_extract.skeleton:run

Then run `python setup.py install` which will install the command `fibonacci`
inside your current environment.
Besides console scripts, the header (i.e. until _logger...) of this file can
also be used as template for Python modules.

Note: This skeleton file can be safely removed if not needed!
"""

import sys
import argparse
import logging

from alamos_extract import __version__
from alamos_extract.load_data import load_cluster, search_db, Cluster
from alamos_extract.form_dicts import region_dict, subtype_dict, virus_dict


__author__ = "Stephen Gaffney"
__copyright__ = "Stephen Gaffney"
__license__ = "gpl3"

_logger = logging.getLogger(__name__)

REGION_CHOICES = tuple(region_dict)
SUBTYPE_CHOICES = tuple(subtype_dict)
VIRUS_CHOICES = tuple(virus_dict)


def parse_args(args):
    """Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Los Alamos HIV Database cluster/patient extractor.")
    parser.add_argument(
        '--version',
        action='version',
        version='alamos-extract {ver}'.format(ver=__version__))
    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser')

    parser_c = subparsers.add_parser('cluster', help='Cluster search help')
    parser_c.add_argument(
        'cluster_id',
        metavar='cluster_id',
        help="Cluster ID",
        type=int,
        default=701
    )

    # max_rec=100, virus='HIV-1', subtype='A1*', region='GENOME'
    parser_s = subparsers.add_parser('cluster_name', help='Cluster name search')
    parser_s.add_argument('cluster_name', default=None, help='Cluster name (not integer ID)')
    parser_s.add_argument('-t', '--virus', nargs='?', default='HIV-1', help='Virus', choices=VIRUS_CHOICES)
    parser_s.add_argument('-s', '--subtype', nargs='?', default='any', help='Subtype', choices=SUBTYPE_CHOICES)
    parser_s.add_argument('-r', '--region', nargs='?', default='any', help='Region', choices=REGION_CHOICES)
    parser_s.add_argument('-m', '--maxrows', nargs='?', default=100, type=int, help='Max row count')

    parser.add_argument(
        '-v',
        '--verbose',
        dest="loglevel",
        help="set loglevel to INFO",
        action='store_const',
        const=logging.INFO)
    parser.add_argument(
        '-vv',
        '--very-verbose',
        dest="loglevel",
        help="set loglevel to DEBUG",
        action='store_const',
        const=logging.DEBUG)
    return parser.parse_args(args)


def setup_logging(loglevel):
    """Setup basic logging

    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    logformat = "[%(asctime)s] %(message)s"
    logging.basicConfig(level=loglevel, stream=sys.stdout,
                        format=logformat, datefmt="%Y-%m-%d %H:%M:%S")


def main(args):
    """Main entry point allowing external calls

    Args:
      args ([str]): command line parameter list
    """
    args = parse_args(args)
    setup_logging(args.loglevel)

    if args.subparser == 'cluster':
        cluster_id = args.cluster_id
        _logger.debug("Parsing cluster ID {}".format(cluster_id))
        c = Cluster(cluster_id)
        path_accession = 'cluster_{}_accessions.tsv'.format(cluster_id)
        path_clinical = 'cluster_{}_clinical.tsv'.format(cluster_id)
        c.acc_df.to_csv(path_accession, sep='\t', index=False)
        c.desc_df.to_csv(path_clinical, sep='\t', index=True)

        patient_names = ', '.join(c.comb_patients.keys())
        n_patients = len(c.patient_dict)
        n_accessions = len(c.acc_df)

        # data = load_cluster(cluster_id)
        print('Cluster: {}'.format(c.cluster_name))
        print('Description: {}'.format(c.description))
        print('{} patients: {}'.format(n_patients, patient_names))
        print('{} accessions.'.format(n_accessions))
        print('Clinical data written to {}'.format(path_clinical))
        print('Accession data written to {}'.format(path_accession))

    elif args.subparser == 'cluster_name':
        cluster_name = args.cluster_name
        virus = virus_dict[args.virus]
        subtype = subtype_dict[args.subtype]
        df = search_db(max_rec=args.maxrows, virus=virus,
                       subtype=subtype, cluster_name=args.cluster_name,
                       region=None)
        df.drop(['blast', 'blast2'], axis=1, inplace=True)
        clusters = list(df['cluster_comb'].unique())
        n_clusters = len(clusters)
        if n_clusters == 1:
            clusters = clusters[0]
        _logger.info("%d clusters identified: %s", n_clusters, clusters)
        tsv_path = f'cluster_{cluster_name}_info.tsv'
        df.to_csv(tsv_path, sep='\t', index=False)
        _logger.info(f"Cluster sequence metadata saved to {tsv_path}.")
    _logger.info("Script complete.")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
