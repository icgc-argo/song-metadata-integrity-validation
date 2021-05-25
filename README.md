# ICGC ARGO RDPC SONG Metadata Integrity Validation

This repository contain scripts to perform integrity validation on metadata recorded in
SONG analysis objects. Examples of integrity violation: duplicated alignment results,
missing qc_metrics analysis, WXS tumour sequence paired with WGS normal sequence etc.


## Create Elasticsearch index for development

Some SONG metadata have been dumped as JSON documents and stored under `song-analysis-dump`,
the following example command load SONG analysis objects from `PACA-CA` study (also known as program) into an
Elasticsearch index called `song-2021-05-25`, run:

```
./scripts/indexer.py -i song-2021-05-25 -d song-analysis-dump/PACA-CA.20210507.json song-analysis-dump/PACA-CA.inter.20210512.json
```

Note: please make sure an Elasticsearch server run at localhost:9200

## Run validation

The follow example command generates validation reports for all SONG objects in `PACA-CA` study.
```
./scripts/validator.py -i song-2021-05-25 -p PACA-CA
```

You can report on a particular donor, for example use `-d DO230463` to report only for this donor.

This still a **WIP**, no actual reports yet.
