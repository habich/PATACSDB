#!/usr/bin/python
from PATACSDB import app
from flask import render_template
from flask import url_for
import schema
import json
from schema import db

def make_link(polyA):
    return "http://www.ensembl.org/"+polyA.organism_name+"/Location/View?db=core;g="+polyA.gene_id+";r="+polyA.chromosome_name+":"+\
        str(polyA.chromosome_location_start)+"-"+str(polyA.chromosome_location_stop)+";t="+polyA.transcript_id

def make_dict(polyA):
    return {
     "gene_id":polyA.gene_id,
     "transcript_id":polyA.transcript_id,
     "protein_id":polyA.protein_id,
     "protein_name":polyA.protein_name,
     "AAA_sequence":polyA.AAA_sequence,
     "AAA_hem_sequence":polyA.AAA_hem_sequence,
     "ensembl_link": make_link(polyA),
     "GO_ids": ", ".join([x.GO_id for x in polyA.GOs])}   

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/database.json")
def database_json():
    list_of_organisms = [{"organism_name":x.species_name,"organism_link":url_for('database',organism=x.id)} for x in db.session.query(schema.Species)]
    return json.dumps(list_of_organisms)

@app.route("/<organism>.json")
def show_organism(organism):
    database = []
    for i in db.session.query(schema.PolyA).filter(schema.PolyA.organism_name == organism).all():
            database.append(make_dict(i))
            
    return json.dumps(database)
    
@app.route("/database")    
@app.route("/database/<organism>")
def database(organism=None):
    if organism:
        return render_template("organism.html",organism_url=url_for("show_organism",organism=organism))
    else:
        return render_template("database.html")

