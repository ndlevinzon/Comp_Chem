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
SELECT DISTINCT m.chembl_id AS compound_chembl_id,
s.canonical_smiles,
r.compound_key
FROM compound_structures s
RIGHT JOIN molecule_dictionary m ON s.molregno = m.molregno
JOIN compound_records r ON m.molregno = r.molregno
AND m.max_phase      = 4;
```
When you run the above query, it generates a long list of ~4000 molecules that are in Phase 4 (or are approved)

ChEMBL documentation has an excellent tutorial on crafting SQL queries and the schema and was very helpful to get started. I am a SQL novice and whatever I learned is from the Intro to SQL course from DataCamp. I highly recommend it for folks who want to learn the basics of SQL.

I mainly prefer using Python as I do all my exploratory data analysis, so I wanted to figure out a way to interact with the SQL database with Python. Packages such as pychembldb, razi are summarized in Iwatobipen’s blogs and are great starting points. I found the psycopg2 package to be useful as it can keep the feel of using SQL queries intact, and there is a ton of support for this package in StackOverflow.
