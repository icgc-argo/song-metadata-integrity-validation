#!/usr/bin/env python3

import os
import sys
import json
import argparse
from elasticsearch import Elasticsearch


BASE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
ES_AGGS = os.path.join(BASE_DIR, 'es-queries', 'es-aggs.json')


def get_donors_by_program(program_id, es, index_name):
    ret = es.search(body={
        "query": {
            "terms": {
                "studyId": [program_id]
            }
        },
        "aggs": {
            "donor_id": {
                "terms": {
                    "field": "samples.donor.donorId",
                    "size": 100000
                }
            }
        }
    }, index=index_name, _source=False)

    donor_ids = []
    for d in ret['aggregations']['donor_id']['buckets']:
        donor_ids.append(d['key'])

    return donor_ids


def validate_sample_and_seq_exp(program_id, donor_id, analysis_objs):
    """
    # sample_info dict populated from using all analysis objects
    {
        "{sample_id}": {
            "analysis_ids": [],
            "sample_id": "",
            "submitter_sample_id": {"{value}": []},  # keys are analysis_ids
            "sample_type": [],
            "matched_normal_sample_id": [],
            "submitter_matched_normal_sample_id": [],
            "specimen_id": [],
            "specimen_type": [],
            "tumour_normal_designation": [],
            "donor_id": [],
            "study_id": [],
            "sequencing_experiment": {
                "WGS": [],
                "WXS": [],
                "RNA-Seq": []
            }
        }
    }
    """
    sample_fields = {
        'submitterSampleId': 'samples.submitterSampleId',
        'sampleType': 'samples.sampleType',
        'matchedNormalSubmitterSampleId': 'samples.matchedNormalSubmitterSampleId',
        'specimenId': 'samples.specimen.specimenId',
        'specimenType': 'samples.specimen.specimenType',
        'tumourNormalDesignation': 'samples.specimen.tumourNormalDesignation',
        'specimenTissueSource': 'samples.specimen.specimenTissueSource',
        'donorId': 'samples.donor.donorId',
        'gender': 'samples.donor.gender',
        'studyId': 'studyId'
    }

    # gather sample information from all analysis objects
    sample_info = dict()
    for a in analysis_objs:
        analysisId = a['analysisId']
        sampleId = a['samples'][0]['sampleId']
        if sampleId not in sample_info:
            sample_info[sampleId] = {
                'sampleId': sampleId,
                'analysisId': [],
                'sequencing_experiment': {}
            }

        sample_info[sampleId]['analysisId'].append(analysisId)

        for field in sample_fields:
            if field not in sample_info[sampleId]:
                sample_info[sampleId][field] = dict()

            source_path = sample_fields[field].split('.')
            if len(source_path) == 1:
                source_value_str = str(a[source_path[0]])
            elif source_path[0] == 'samples':
                if len(source_path) == 2:
                    source_value_str = str(a['samples'][0][source_path[1]])
                elif len(source_path) == 3:
                    source_value_str = str(a['samples'][0][source_path[1]][source_path[2]])
                else:
                    assert False  # not supposed to reach here
            else:
                assert False  # not supposed to reach here

            if source_value_str not in sample_info[sampleId][field]:
                sample_info[sampleId][field][source_value_str] = []

            sample_info[sampleId][field][source_value_str].append(analysisId)

        # add sequencing_experiment
        if a.get('analysisType', {}).get('name') in ('sequencing_experiment', 'rna_sequencing_experiment'):
            strategy = a['experiment']['experimental_strategy']
            matchedNormalSubmitterSampleId = a['samples'][0]['matchedNormalSubmitterSampleId']

            if strategy not in sample_info[sampleId]['sequencing_experiment']:
                sample_info[sampleId]['sequencing_experiment'][strategy] = {
                        'sequencing_experiment_analysis_id': [analysisId],
                        'matchedNormalSubmitterSampleId': [matchedNormalSubmitterSampleId]
                    }
            else:
                sample_info[sampleId]['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id'].append(analysisId)
                sample_info[sampleId]['sequencing_experiment'][strategy]['matchedNormalSubmitterSampleId'].append(matchedNormalSubmitterSampleId)

    # print(json.dumps({'donor_id': donor_id}))
    # print(json.dumps(sample_info))

    samples, issues = value_discrepancy_check(sample_info, sample_fields)
    # print(json.dumps(samples))
    # print(json.dumps(issues))

    submitterSampleId2SampleId = mapping_sumbitter_sample_id_to_sample_id(samples)
    # print(json.dumps(submitterSampleId2SampleId))

    # figure out tumour-normal pairs
    tumour_normal_pairs, tumour_not_paired, normal_not_paired, normal_paired, issues = resolve_tumour_normal_pairs(samples, submitterSampleId2SampleId)
    # print(json.dumps(tumour_normal_pairs))
    # print(json.dumps(normal_paired))
    # print(json.dumps(tumour_not_paired))
    # print(json.dumps(normal_not_paired))
    # print(json.dumps(issues))


def resolve_tumour_normal_pairs(samples, submitterSampleId2SampleId):
    tumour_normal_pairs = dict()
    normal_paired = dict()
    tumour_not_paired = dict()
    normal_not_paired = dict()
    issues = dict()

    for s in samples:
        sample = samples[s]
        if sample.get('tumourNormalDesignation') != 'Tumour' or \
                not sample.get('sequencing_experiment'):
            continue

        # go through sequencing experiment with different experimental strategies
        for strategy in sample['sequencing_experiment']:
            matched_normal_sample_id = submitterSampleId2SampleId.get(
                    sample['sequencing_experiment'][strategy]['matchedNormalSubmitterSampleId']
                )

            if not matched_normal_sample_id:
                if sample['sampleId'] not in issues:
                    issues[sample['sampleId']] = []

                issues[sample['sampleId']].append(
                        f"Failed to get matched normal sample, perhaps no sequencing data has been submitted for the normal sample. Strategy: {strategy}, "
                        f"sequencing_experiment: {sample['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id']}, "
                        f"matchedNormalSubmitterSampleId: {sample['sequencing_experiment'][strategy]['matchedNormalSubmitterSampleId']}"
                    )
            else:
                if samples.get(matched_normal_sample_id, {}).get('sequencing_experiment', {}).get(strategy):
                    matched_normal_sequencing_experiment_analysis_id = samples[matched_normal_sample_id]['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id']
                    matchedNormalSubmitterSampleId = sample['sequencing_experiment'][strategy]['matchedNormalSubmitterSampleId']

                    if sample['sampleId'] not in tumour_normal_pairs:
                        tumour_normal_pairs[sample['sampleId']] = dict()

                    tumour_normal_pairs[sample['sampleId']][strategy] = {
                            'sequencing_experiment_analysis_id': sample['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id'],
                            'matched_normal_sequencing_experiment_analysis_id': matched_normal_sequencing_experiment_analysis_id,
                            'matchedNormalSubmitterSampleId': matchedNormalSubmitterSampleId,
                            'matchedNormalSampleId': matched_normal_sample_id
                        }

                    if matched_normal_sample_id not in normal_paired:
                        normal_paired[matched_normal_sample_id] = dict()

                    if strategy not in normal_paired[matched_normal_sample_id]:
                        normal_paired[matched_normal_sample_id][strategy] = {
                            'sequencing_experiment_analysis_id': matched_normal_sequencing_experiment_analysis_id,
                            'matched_tumour_sequencing_experiment_analysis_ids': [sample['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id']]
                        }
                    else:
                        normal_paired[matched_normal_sample_id][strategy]['matched_tumour_sequencing_experiment_analysis_ids'].append(sample['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id'])

                else:  # not able to find seq exp with the same strategy for matched normal sample
                    if sample['sampleId'] not in tumour_not_paired:
                        tumour_not_paired[sample['sampleId']] = dict()

                    tumour_not_paired[sample['sampleId']][strategy] = {
                            'sequencing_experiment_analysis_id': sample['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id'],
                            'matched_normal_sequencing_experiment_analysis_id': None
                        }

    # figure out normals that are not paired
    for s in samples:
        sample = samples[s]
        if sample.get('tumourNormalDesignation') != 'Normal' or \
                not sample.get('sequencing_experiment'):
            continue

        if sample['sampleId'] in normal_paired:
            continue

        for strategy in sample['sequencing_experiment']:
            if sample['sampleId'] in normal_paired and strategy in normal_paired[sample['sampleId']]:
                continue

            if sample['sampleId'] not in normal_not_paired:
                normal_not_paired[sample['sampleId']] = {
                        strategy: {
                            'sequencing_experiment_analysis_id': sample['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id']
                        }
                    }

    return tumour_normal_pairs, tumour_not_paired, normal_not_paired, normal_paired, issues


def mapping_sumbitter_sample_id_to_sample_id(samples):
    submitterSampleId2SampleId = dict()
    for s_id in samples:
        if 'submitterSampleId' in samples[s_id]:
            submitterSampleId2SampleId[samples[s_id]['submitterSampleId']] = samples[s_id]['sampleId']

    return submitterSampleId2SampleId


def value_discrepancy_check(sample_info, sample_fields):
    samples = dict()
    issues = dict()

    for sample_id in sample_info:
        sample = {
            'sampleId': sample_id
        }
        for field in sample_fields:
            if len(sample_info[sample_id][field]) > 1:
                if sample_id not in issues:
                    issues[sample_id] = {field: sample_info[sample_id][field]}
                else:
                    issues[sample_id].update({field: sample_info[sample_id][field]})
            else:
                sample[field] = list(sample_info[sample_id][field].keys())[0]
                if sample[field] == 'None':
                    sample[field] = None

        # sequencing_experiment
        if 'sequencing_experiment' in sample_info[sample_id]:
            for strategy in sample_info[sample_id]['sequencing_experiment']:
                if len(sample_info[sample_id]['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id']) > 1:
                    if sample_id not in issues:
                        issues[sample_id] = {
                                'sequencing_experiment': {
                                    strategy: sample_info[sample_id]['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id']
                                }
                            }
                    elif 'sequencing_experiment' not in issues[sample_id]:
                        issues[sample_id]['sequencing_experiment'] = {
                                strategy: sample_info[sample_id]['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id']
                            }
                    else:
                        issues[sample_id]['sequencing_experiment'][strategy] = \
                            sample_info[sample_id]['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id']

                else:
                    if 'sequencing_experiment' not in sample:
                        sample['sequencing_experiment'] = dict()

                    sample['sequencing_experiment'][strategy] = {
                                'sequencing_experiment_analysis_id': sample_info[sample_id]['sequencing_experiment'][strategy]['sequencing_experiment_analysis_id'][0],
                                'matchedNormalSubmitterSampleId': sample_info[sample_id]['sequencing_experiment'][strategy]['matchedNormalSubmitterSampleId'][0]
                            }

        samples[sample_id] = sample

    return (samples, issues)


def validate_donor(program_id, donor_id, analysis_objs, aggs, output_dir):
    """
    # for debugging
    with open(os.path.join(output_dir, f"{donor_id}.analyses.json"), "w") as d:
        d.write(json.dumps(analysis_objs, indent=2))

    with open(os.path.join(output_dir, f"{donor_id}.aggs.json"), "w") as d:
        d.write(json.dumps(aggs, indent=2))
    """
    samples = validate_sample_and_seq_exp(program_id, donor_id, analysis_objs)


def get_data_from_es(es, donor_id, index_name):
    aggs_query = json.load(open(ES_AGGS, 'r'))
    aggs_query.update({
            "query": {
                "terms": {
                    "samples.donor.donorId": [donor_id]
                }
            }
        })

    ret = es.search(
        body=aggs_query,
        index=index_name,
        size=10000  # assume a donor does not have over 10,000 analysis objects
    )

    hits = ret['hits'].get('hits', [])
    aggs = ret['aggregations']
    analysis_objs = []
    for hit in hits:
        analysis_objs.append(hit['_source'])

    return analysis_objs, aggs


def main(
            es,
            index_name,
            program_id,
            donor_ids,
            output_dir
        ):

    if not donor_ids:
        donor_ids = get_donors_by_program(program_id, es, index_name)

    for donor_id in donor_ids:
        analysis_objs, aggs = get_data_from_es(es, donor_id, index_name)
        validate_donor(program_id, donor_id, analysis_objs, aggs, output_dir)


if __name__ == '__main__':
    es = Elasticsearch()

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='index_name', required=True, help='Index name')
    parser.add_argument('-p', dest='program_id', required=True, help='Program ID')
    parser.add_argument('-d', dest='donor_ids', nargs='+', help='Donor ID(s)')
    parser.add_argument('-f', dest='donor_file', help='File contains donor IDs')
    parser.add_argument('-w', dest='workflow', choices=['dna-alignment', 'sanger-wgs', 'sanger-wxs', 'mutect2'],
                        help='Program ID')
    parser.add_argument('-o', dest='output_dir', default='.', help='Path for output directory')
    args = parser.parse_args()

    donor_ids = []
    if args.donor_ids and args.donor_file:
        sys.exit("Can not specify both 'donor_ids' and 'donor_file'")
    elif args.donor_file:
        with open(args.donor_file, 'r') as d:
            for row in d:
                donor_ids.append(row.strip())
    elif args.donor_ids:
        donor_ids = args.donor_ids

    main(
        es,
        args.index_name,
        args.program_id,
        donor_ids,
        args.output_dir
    )
