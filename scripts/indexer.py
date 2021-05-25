#!/usr/bin/env python3

import os
import sys
import json
import argparse
from elasticsearch import Elasticsearch

base_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

es = Elasticsearch()
es_mapping_file = os.path.join(base_dir, 'mappings', 'song-analysis.mapping.json')
with open(es_mapping_file, 'r') as m:
    es_mapping = json.load(m)


def main(index_name, analysis_dump):
    if not es.indices.exists(index_name):
        es.indices.create(index_name, es_mapping)

    for f in analysis_dump:
        with open(f, 'r') as j:
            analysis_obj = json.load(j)
            print(f'Start indexing {len(analysis_obj)} SONG analysis objects ...')
            for a in analysis_obj:
                doc_id = a['analysisId']
                if 'library_strategy' in a['experiment'] and 'experimental_strategy' not in a['experiment']:
                    a['experiment']['experimental_strategy'] = a['experiment']['library_strategy']

                es.index(index=index_name, id=doc_id, body=a)

            print(f'Indexed {len(analysis_obj)} SONG analysis objects ...')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--index-name',
        required=True,
        help='Name of the ES index'
    )
    parser.add_argument(
        '-d',
        '--analysis-dump',
        nargs='+',
        required=True,
        help='SONG JSON dump created from SONG analysis endpoint output'
    )

    args = parser.parse_args()

    main(args.index_name, args.analysis_dump)
