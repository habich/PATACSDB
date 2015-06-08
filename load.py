#!/usr/bin/python

import glob
import sys
from itertools import izip_longest
import PATACSDB.schema as schema
from PATACSDB.schema import db
    
if __name__ == "__main__":    
    if len(sys.argv) != 2:
        print "Usage: "+sys.argv[0]+" <directory with AAA data>"
        exit(1)    

    path = sys.argv[1]
    files = glob.glob(path+"/*_pairs.data")

    db.drop_all()
    db.create_all()

    header = ["protein_id","transcript_id", "gene_id", "cdna_start", "cdna_stop", "AAA_start", "AAA_len", "AAA_stop",
              "AAA_sequence", "AAA_hem_start", "AAA_hem_stop", "AAA_hem_sequence", "protein_sequence", "nucleotide_sequence",
              "chromosome_name", "chromosome_location_start", "chromosome_location_stop", "strand", "protein_name"]#, "GOs" 
    for f in files:
        species_name = "_".join(f.split("/")[-1].split("_")[:-1])
        species = schema.Species(id = species_name, species_name = species_name[0].upper() + " ".join(species_name[1:].split("_")))
        db.session.add(species)
        with open(f) as csv:
            for l in csv:
                line = [x.rstrip("\n") for x in l.split(",")]
                values = dict(izip_longest(header,line,fillvalue=None))
                values['organism_name'] = species.id
                polyA = schema.PolyA(**values)
                db.session.add(polyA)
        db.session.commit()
                