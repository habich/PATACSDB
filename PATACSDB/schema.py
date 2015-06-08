#!/usr/bin/python
from PATACSDB import app
from flask.ext.sqlalchemy import SQLAlchemy
from sqlalchemy import Column, String, Integer, Text, ForeignKey
from sqlalchemy.orm import relationship

db = SQLAlchemy(app)

class Species(db.Model):
    id = Column(String(100),primary_key=True)
    species_name = Column(String(150),nullable = False)

class PolyA(db.Model):
    id = Column(Integer, primary_key = True)
    organism_name = Column(String(100), ForeignKey('species.id'))
    species = relationship('Species')
    protein_id = Column(String(30), nullable = False)
    transcript_id = Column(String(30),nullable=False)
    gene_id = Column(String(30),nullable = False)
    cdna_start = Column(Integer,nullable = False)
    cdna_stop = Column(Integer, nullable = False)
    AAA_start = Column(Integer, nullable = False)
    AAA_len = Column(Integer, nullable = False)
    AAA_stop = Column(Integer, nullable = False)
    AAA_sequence = Column(String(250),nullable = False)
    AAA_hem_start = Column(Integer, nullable = False)
    AAA_hem_stop = Column(Integer, nullable = False)
    AAA_hem_sequence = Column(String(300), nullable = False)
    protein_sequence = Column(Text, nullable = False)
    nucleotide_sequence = Column(Text, nullable = False)
    chromosome_name = Column(String(100), nullable = True)
    chromosome_location_start = Column(Integer, nullable = True)
    chromosome_location_stop = Column(Integer, nullable = True)
    strand = Column(String(50),nullable = True)
    protein_name = Column(String(150),nullable = True)
    GOs = relationship('GO')
   
class GO(db.Model):
    id = Column(Integer, primary_key = True)
    polyA_id = Column(Integer, ForeignKey('polyA.id'))
    GO_id = Column(String(30),nullable = False )
    
