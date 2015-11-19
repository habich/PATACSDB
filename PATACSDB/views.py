#!/usr/bin/python
from PATACSDB import app
from flask import render_template
from flask import url_for
import schema
import json
from schema import db
from sqlalchemy import func
from sqlalchemy import distinct
import math
genomes = None
list_of_organisms = None
the_logest_AAA = {}
margin = 7

def ensembl_link(database):
    if database:
        return "http://"+database+".ensembl.org/"
    else:
        return "http://www.ensembl.org/"

def gene_link(cdna,species):
    if species.biomart_database:
        return "http://central.biomart.org/martreport/?report=report&mart=gene_ensembl_report&link_ensembl_gene_id="+cdna.gene_id+"&datasets="+species.biomart_database
    else:
        return None

def transcript_link(cdna,species):
    return ensembl_link(species.database) + cdna.organism_name+"/Transcript/Summary?db=core;g="+cdna.gene_id+";t="+cdna.transcript_id

def sequence_link(cdna,gtf,species,polyA):
    return ensembl_link(species.database) +cdna.organism_name+"/Location/View?db=core;g="+cdna.gene_id+";r="+gtf.chromosome_name+":"+\
        str(polyA.AAA_global_start)+"-"+str(polyA.AAA_global_stop)+";t="+cdna.transcript_id

def location(cdna,polyA):
    return str("{0:.2f}".format((float(polyA.AAA_start)*100)/len(cdna.nucleotide_sequence)))

def protein_name(gtf,pep):
    if gtf.protein_name:
        return gtf.protein_name
    else:
        return pep.protein_id

def make_sequence(cdna,polyA,species):
    global the_logest_AAA
    if species not in the_logest_AAA:
        the_logest_AAA[species] = db.session.query(func.max(schema.PolyA.AAA_stop-schema.PolyA.AAA_start)).join(schema.Cdna).filter(schema.Cdna.organism_name == species).all()[0][0]
    sequence = ""
    start_character = "[START]"
    stop_charaterer = "[STOP]"
    start_margin = (the_logest_AAA[species]+margin)*"-"+start_character
    stop_margin = stop_charaterer+(the_logest_AAA[species]+margin)*"-"
    AAA_len = polyA.AAA_stop-polyA.AAA_start
    left_len = int(math.floor(float(the_logest_AAA[species]+2*margin-AAA_len)/2))
    right_len = int(math.ceil(float(the_logest_AAA[species]+2*margin-AAA_len)/2))
    if polyA.AAA_start - left_len < 0:
        sequence+= start_margin[polyA.AAA_start - left_len:]+cdna.nucleotide_sequence[:polyA.AAA_start]
    else:
        sequence+=cdna.nucleotide_sequence[polyA.AAA_start - left_len:polyA.AAA_start]
    sequence += cdna.nucleotide_sequence[polyA.AAA_start:polyA.AAA_stop]
    if polyA.AAA_stop + right_len > len(cdna.nucleotide_sequence):
        sequence+= cdna.nucleotide_sequence[polyA.AAA_stop:]+stop_margin[:polyA.AAA_stop + right_len - len(cdna.nucleotide_sequence)]
    else:
        sequence+= cdna.nucleotide_sequence[polyA.AAA_stop:polyA.AAA_stop+right_len]
    return (sequence,left_len,left_len+AAA_len)

def make_dict(cdna,gtf,pep,polyA,species):
    sequence_info = make_sequence(cdna,polyA,species.id)
    return {
        "protein_name" : protein_name(gtf,pep),
        "gene_id" : cdna.gene_id,
        "gene_link" : gene_link(cdna,species),
        "transcript_id" : cdna.transcript_id,
        "transcript_link" : transcript_link(cdna,species),
        "location" : location(cdna,polyA),
        "sequence_info" : { "sequence": sequence_info[0],
                           "AAA_start":sequence_info[1],
                           "AAA_stop":sequence_info[2],
                           "sequence_link" : sequence_link(cdna,gtf,species,polyA)}
     }

@app.route("/")
def index():
    mainMetadata = db.session.query(schema.MainMetadata).first()
    return render_template("index.html", genomes = mainMetadata.genome, genes = mainMetadata.genes, polya = mainMetadata.polyA)

@app.route("/database.json")
def database_json():
    list_of_organisms = [{"organism_name":organismMetadata.organism_name,
                    "organism_link":url_for('database',organism=organismMetadata.organism_id),
                    "protein_number":organismMetadata.protein_number,
                    "transcript_number":organismMetadata.transcript_number,
                    "transcript_with_polya":organismMetadata.trans_pol } for organismMetadata in db.session.query(schema.OrganismMetadata)]
    return json.dumps(list_of_organisms)

@app.route("/<organism_id>.json")
def show_organism(organism_id):
    database = []
    for cdna,gtf,pep,polyA,species in db.session.query(schema.Cdna, schema.Gtf , schema.Pep, schema.PolyA, schema.Species ).join(schema.PolyA).join(schema.Pep).join(schema.Gtf).join(schema.Species).filter(schema.Cdna.organism_name==organism_id):
        database.append(make_dict(cdna,gtf,pep,polyA,species))
    return json.dumps(database)

@app.route("/database")
@app.route("/database/<organism>")
def database(organism=None):
    if organism:
        return render_template("organism.html",organism_url=url_for("show_organism",organism_id=organism))
    else:
        return render_template("database.html")
