select distinct
    G.stable_id,
    X.`display_label`
 from
   transcript as G,
   object_xref as OX,
   xref as X ,
   external_db as D
 where
   D.external_db_id=X.external_db_id and
   OX.xref_id=X.xref_id and
   OX.ensembl_object_type="Gene" and
   G.gene_id=OX.ensembl_id and
   D.db_name like 'HGNC%';
