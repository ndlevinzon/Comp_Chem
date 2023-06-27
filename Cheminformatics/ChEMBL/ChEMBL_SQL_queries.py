import psycopg2
import pandas as pd
import pandas.io.sql as sqlio

con = psycopg2.connect(database="chembl_29", user="leelasd", password="", host="127.0.0.1", port="5432")

print("Database opened successfully")

sql="""SELECT
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
"""
dat = sqlio.read_sql_query(sql, con)
dat.to_csv('single_protein_actives.csv',index=False)
print("Operation done successfully")
con.close()

