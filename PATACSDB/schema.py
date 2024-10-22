#!/usr/bin/python
from PATACSDB import app
from flask.ext.sqlalchemy import SQLAlchemy
from sqlalchemy import Column, String, Integer, Text, ForeignKey
from sqlalchemy.orm import relationship

db = SQLAlchemy(app)

class Species(db.Model):
    id = Column(String(100),primary_key=True)
    database =  Column(String(20),nullable = False)
    species_name = Column(String(150),nullable = False)
    transcript = relationship('Cdna', backref='species')
    biomart_database = Column(String(150),nullable=True)

class PolyA(db.Model):
    id = Column(Integer, primary_key = True) 
    transcript_id = Column(String(30), ForeignKey('cdna.transcript_id'), nullable=False, index=True)
    AAA_start = Column(Integer, nullable = False)
    AAA_stop = Column(Integer, nullable = False)
    AAA_global_start = Column(Integer, nullable = False)
    AAA_global_stop = Column(Integer, nullable = False)
    
class Cdna(db.Model):
    transcript_id = Column(String(30),nullable=False, primary_key = True)
    gene_id = Column(String(30),nullable = False)
    nucleotide_sequence = Column(Text, nullable = False)
    protein = relationship('Pep', backref='cdna')    
    cdna_start = Column(Integer,nullable = False)
    cdna_stop = Column(Integer, nullable = False)
    organism_name = Column(String(100), ForeignKey('species.id'), index=True)
    polyA = relationship('PolyA', backref='cdna')

class Pep(db.Model):
    protein_id = Column(String(30), nullable = False, primary_key = True)
    transcript_id = Column(String(30), ForeignKey('cdna.transcript_id'), nullable=False, index=True)
    protein_sequence = Column(Text, nullable = False)

class Gtf(db.Model):
    transcript_id = Column(String(30),ForeignKey('cdna.transcript_id'),nullable=False, primary_key = True)
    chromosome_name = Column(String(100), nullable = True)
    chromosome_location_start = Column(Integer, nullable = True)
    chromosome_location_stop = Column(Integer, nullable = True)
    strand = Column(String(50),nullable = True)
    protein_name = Column(String(150),nullable = True)
    
class MainMetadata(db.Model):
    id = Column(Integer, primary_key = True)
    genome = Column(Integer,nullable = False)
    genes = Column(Integer, nullable = False)
    polyA = Column(Integer, nullable = False)
    
class OrganismMetadata(db.Model):
    organism_id = Column(String(100), ForeignKey('species.id'), nullable = False, primary_key = True)
    organism_name = Column(String(100), ForeignKey('species.species_name'), nullable = False)
    protein_number = Column(Integer, nullable = False)
    transcript_number = Column(Integer, nullable = False)
    trans_pol = Column(Integer, nullable = False)
