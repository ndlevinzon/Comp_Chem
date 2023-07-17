# Working with ChEMBL
ChEMBL database is the primary source of data for doing any large-scale analysis in early Drug Discovery. Accessing data from ChEMBL is often seamless, which I have been using for the last three years. 
But sometimes, getting data from ChEMBL API is slow, and when you are still exploring what kind of data you need and how the query should look, it can take forever. Local installation of the ChEMBL database helps alleviate the issue and makes data curation fast and more workable. 
In this tutorial, I am using the local ChEMBL SQL database and the psycopg2 python package to curate a list of chemical descriptors for all active compounds against single-protein targets.
What do you need to follow this blog?

To follow this tutorial, you need to have a local installation of the ChEMBL SQL cartridge. Follow the steps or the link here: https://iwatobipen.wordpress.com/2020/06/02/communicate-chembl27-with-rdkit-postgres-cartridge-and-sqlalchemy-rdkit-chembl-postgres-razi/comment-page-1/
```
$ conda install -c conda-forge postgresql
$ conda install -c rdkit rdkit-postgresql
$ initdb -D ~/postgresdata
$ postgers -D ~/postgresdaa
$ cd ~/Downloads
$ wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_33_postgresql.tar.gz # Or whatever the latest DB version is (https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/)
$ tar -vxf  chembl_33_postgresql.tar.gz
$ cd chembl_33_postgresql/chembl_33/chemb_33_postgresql
$ psql -U postgres
pgdb=# create database chembl_33;
pgdb=# \q;
$ pg_restore --no-owner -h localhost -p 5432 -U iwatobipen -d chembl_27 chembl_27_postgresql.dmp
```
Now we have installed chembl33 to our postgresql DB. From here, we can make rdk schema. Details are described in the rdkit original document (https://www.rdkit.org/docs/Cartridge.html)
```
chembl_33=# create extension if not exists rdkit;
chembl_33=# create schema rdk;
chembl_33=# select * into rdk.mols from (select molregno,mol_from_ctab(molfile::cstring) m  from compound_structures) tmp where m is not null;
chembl_33=# create index molidx on rdk.mols using gist(m);
CREATE INDEX
chembl_33=# alter table rdk.mols add primary key (molregno);
ALTER TABLE
chembl_33=# select molregno,torsionbv_fp(m) as torsionbv,morganbv_fp(m) as mfp2,featmorganbv_fp(m) as ffp2 into rdk.fps from rdk.mols;
chembl_33=# create index fps_ttbv_idx on rdk.fps using gist(torsionbv);
CREATE INDEX
chembl_33=# create index fps_mfp2_idx on rdk.fps using gist(mfp2);
CREATE INDEX
chembl_33=# create index fps_ffp2_idx on rdk.fps using gist(ffp2);
CREATE INDEX
chembl_33=# alter table rdk.fps add primary key (molregno);
ALTER TABLE
```
Now we can make rdk.mols and rdk.fps table in our chembl_33 DB. This command below starts a server that we can use to interact with the ChEMBL33 SQL database.
```
postgres -D ~/postgresdata
```
Using local ChEMBL SQL database

Now we can use SQL language to interact with the database. For this demo, let us try to extract all the approved drugs. Here is how that SQL query looks like
```
SELECT
    act.standard_type,
    act.standard_relation,
    act.standard_value,
    act.standard_units,
    act.activity_comment,
    cmp.canonical_smiles,
    cmp.chembl_id,
    tgt.pref_name AS target_name,
    tgt.chembl_id AS target_chembl_id
FROM
    activities act
JOIN
    assays ass ON act.assay_id = ass.assay_id
JOIN
    compound_structures cmp ON act.molregno = cmp.molregno
JOIN
    target_dictionary tgt ON ass.tid = tgt.tid
WHERE
    act.standard_type IS NOT NULL
    AND act.standard_relation IS NOT NULL
    AND act.standard_value IS NOT NULL
    AND act.standard_units IS NOT NULL
    AND act.potential_duplicate = 0
    AND act.pchembl_value >= 5
    AND act.standard_relation IN ('=', '<')
    AND act.data_validity_comment IS NULL
    AND act.standard_flag = 1
    AND ass.confidence_score >= 8
    AND ass.relationship_type = 'D'
    AND tgt.target_type = 'SINGLE PROTEIN'
```
When you run the above query, it generates a long list of ~4000 molecules that are in Phase 4 (or are approved)

ChEMBL documentation has an excellent tutorial on crafting SQL queries and the schema and was very helpful to get started. I am a SQL novice and whatever I learned is from the Intro to SQL course from DataCamp. I highly recommend it for folks who want to learn the basics of SQL.

I mainly prefer using Python as I do all my exploratory data analysis, so I wanted to figure out a way to interact with the SQL database with Python. Packages such as pychembldb, razi are summarized in Iwatobipen’s blogs and are great starting points. I found the psycopg2 package to be useful as it can keep the feel of using SQL queries intact, and there is a ton of support for this package in StackOverflow.

# Search_zinc22.py
search_zinc22.py is a script for looking up zinc ids on the zinc22 system. The operation is simple- provide a file containing a list of zincids and a destination file to write to. The script will give you a progress bar as it searches the system. If a database is down, the script will let you know and continue gathering the results it can.
## Usage: 
```
search_zinc22.py [-h] [--vendor-search] [--get-vendors] [--configuration-server-url CONFIGURATION_SERVER_URL] input_file results_out
```

search for smiles by zinc22 id

positional arguments:
```
  input_file            file containing list of zinc ids or vendor codes to look up
  results_out           destination file for output
```
optional arguments:
```
  -h, --help            show this help message and exit
  --vendor-search       look up molecules by vendor code instead of zinc id
  --get-vendors         get vendor supplier codes associated with zinc id
  --configuration-server-url CONFIGURATION_SERVER_URL
                        database containing configuration for zinc22 system
```
The output format is as follows:
```
SMILES ZINC_ID TRANCHE_NAME
```
With --get-vendors or --vendor-search the output format looks like this:
```
SMILES ZINC_ID VENDOR_ID TRANCHE_NAME CATALOG</nowiki>
```
Meaning the script will find all vendor information and smiles associated with the provided zinc ids or vendor codes.

### Usage w/ Bash
```
source /nfs/soft/zinc22/search_zinc/env/bin/activate
python /nfs/soft/zinc22/search_zinc/search_zinc22.py input_zinc_ids.txt output_zinc_ids.txt
python /nfs/soft/zinc22/search_zinc/search_zinc22.py --get-vendors input_zinc_ids.txt output_vendor_ids.txt</nowiki>
```
### Usage w/ Csh
```
source /nfs/soft/zinc22/search_zinc/env/bin/activate.csh
python /nfs/soft/zinc22/search_zinc/search_zinc22.py input_zinc_ids.txt output_zinc_ids.txt
python /nfs/soft/zinc22/search_zinc/search_zinc22.py --get-vendors input_zinc_ids.txt output_vendor_ids.txt</nowiki>
```
## Dealing with NULL
Sometimes a ZINC ID will fail to look up. This could be because a server is down (the script will notify you if this is the case), or because the ID is missing from the system for some reason. In this case, it may be helpful to separate the molecules that didn't look up from the molecules that did. You may want to save them for later when the servers come back online, or to run a deeper search with comb_legacy_files.py (more on this below).

## Example
```
[env]$ cat legitimate_ids.txt > input.txt
[env]$ echo ZINCzz00ZZZZZZZZ >> input.txt
[env]$ echo ZINCyy00AAAAAAAA >> input.txt
[env]$ echo ZINCxx00BBBBBBBB >> input.txt
[env]$ python search_zinc.py input.txt output.txt
Searching Zinc22:  |XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX| 100.0%  0.00s 23/23 complete!
[env]$ grep "_null_" output.txt
_null_ ZINCzz00ZZZZZZZZ H33P270
_null_ ZINCyy00AAAAAAAA H34P280
_null_ ZINCxx00BBBBBBBB H35P290
```

search_zinc.py will not omit IDs that don't look up from the output, instead it will return the zinc id with "_null_" in every other field. Therefore we can use grep to filter our results.

```
[env]$ grep "_null_" output.txt > missing.txt
[env]$ grep -v "_null_" output.txt > found.txt
```

