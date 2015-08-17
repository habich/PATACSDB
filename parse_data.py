#!/usr/bin/python
import glob
import sys
import os
import time
from collections import defaultdict
import PATACSDB.schema as schema
from PATACSDB.schema import db
from Bio import SeqIO
import math
from Bio.SeqUtils import seq3
from Bio.SeqUtils import seq1
from Bio.Blast.Applications import NcbiblastpCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests
ensembl_request = "http://rest.ensemblgenomes.org/lookup/id/"
AAA_id = 0
GO_id = 0

def organisms_files_dictionary(raw_data_directory):
    dic_func = lambda :{"pep":"","gtf":"","cdna":"","database":""}
    org_dict = defaultdict(dic_func)
    for directory,_,files in os.walk(raw_data_directory):
        for f in files:
            name = f.split(".")
            if "abinito" in name:
                continue
            organism = f.split("/")[-1].split(".")[0]
            if organism:
                if name[-3] == "pep":
                    org_dict[organism]["pep"]+= directory+"/"+f
                elif name[-3] == "cdna":
                    org_dict[organism]["cdna"]+= directory+"/"+f
                elif name[-1] == "gtf":
                    org_dict[organism]["gtf"]+= directory+"/"+f
                database_looking = directory.split("/")
                if "metazoa" in database_looking:
                    org_dict[organism]["database"]="metazoa"
                elif "protists" in database_looking:
                    org_dict[organism]["database"]="protists"
    return org_dict

def findPolyA(seq):
    """Finds poly A sequence in |seq|.
    
    >>> findPolyA("A")
    []
    >>> findPolyA("AAAAAAAAAAAA")
    [[0, 12]]
    >>> findPolyA("AAAAAAAAAAAAAAAA")
    [[0, 16]]
    >>> findPolyA("AGAAAAAAAAAAAAAGA")
    [[0, 17]]
    >>> findPolyA("AGAAAAAAAAAAAAAGAGAAAAAAAA\n")
    [[0, 17], [16, 12]]
    >>> findPolyA("GAAAAAAAAAAAAAAAAAAAAAAG")
    [[0, 24]]
    >>> findPolyA("AAAAAAATAAAAAT")
    [[0, 13]]
    """
    polyALength = 12    
    if len(seq) < polyALength:
        return []    
    polyARanges = []
    start = -1
    length = -1
    errors = len([x for x in seq[:polyALength-1] if x != 'A'])
    for i in range(len(seq) - polyALength + 1):
        errors += seq[i + polyALength - 1] != 'A'
        if errors <= 1:
            if start == -1:
                start = i
                length = polyALength
            else:
                length += 1
        elif start != -1:
            polyARanges.append([start, length])
            start = -1
            length = -1
        errors -= seq[i] != 'A'
    if start != -1:
        polyARanges.append([start, length])
        start = -1
        length = -1
    return polyARanges

def database_append(AAA):
        global AAA_id
        AAA_id += 1
        polyA = schema.PolyA(id = AAA_id, transcript_id = cdna_transcript_id, AAA_start = cdna_start + AAA[0], AAA_stop = cdna_start+AAA[0]+AAA[1])
        db.session.add(polyA)

def out_write(organism_out, cdna_transcript_id,cdna_start,AAA,proper_seq):
    organism_out.write("\t".join([cdna_transcript_id,str(cdna_start + AAA[0]),str(cdna_start+AAA[0]+AAA[1]), proper_seq[AAA[0]:AAA[0]+AAA[1]],"\n"]))    

def grab_AAA_information(AAA_list,organism_out, cdna_transcript_id,cdna_start,proper_seq):
    for AAA in AAA_list:
        database_append(AAA)
        out_write(organism_out, cdna_transcript_id,cdna_start,AAA,proper_seq)


def name_import(transcript_id):
    r = requests.get(ensembl_request+transcript_id+"?", headers={ "Content-Type" : "application/json", "species":organism_id.lower()})
    req_count = 0
    while req_count < 5 and not r.ok:
        r = requests.get(ensembl_request+transcript_id+"?", headers={ "Content-Type" : "application/json"})
        req_count += 1
    if not r.ok:
        print "Ensembl request failed"+transcript_id
        protein_name = None
    else:
        decoded = r.json()
        if "display_name" in decoded.keys():
            protein_name = decoded["display_name"]
        else:
            protein_name = None
    return protein_name       

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage: "+sys.argv[0]+" <directory with ensembl data>"
        exit(1)    
    #raw_data input
    raw_data_directory = sys.argv[1]
    db.drop_all()
    db.create_all()
    
    org_dict = organisms_files_dictionary(raw_data_directory)
     
    for organism in org_dict:
        print "Organism: ", organism
        
        organism_out = open(organism+"_pairs.data",'w')
        failed = open(organism+"_failed",'w')
        
        organism_id = organism
        
        species = schema.Species(id = organism_id,database=org_dict[organism]["database"], species_name = " ".join(organism.split("_")))
        db.session.add(species)
        
        pep = SeqIO.parse(open(org_dict[organism]["pep"],'r'),"fasta")
        pep_dict = defaultdict(list)
        
        gtf_file = open(org_dict[organism]["gtf"],'r')
        print "GTF reading in progress"
        
        for gtf in gtf_file:
            t0 = time.time()
            if gtf[0]=="#":
                continue
            else:
                chrom, source, feature, start, stop, score,strand,frame,attribute = gtf.split("\t")
                info = dict([x.split() for x in attribute.split(";") if len(x.split()) == 2])
                if feature == "transcript":
                    transcript_id = info["transcript_id"][1:-1]
                    protein_name = None
                    #protein_name = name_import(transcript_id)
                    qtf_feature = schema.Gtf(transcript_id = transcript_id, chromosome_name = chrom, chromosome_location_start = start, chromosome_location_stop = stop, strand = strand, protein_name = protein_name)
                    db.session.add(qtf_feature)
                else:
                    continue

                    
        cdna = list(SeqIO.parse(open(org_dict[organism]["cdna"],'r'),"fasta"))
        cdna_size = len(cdna)
        cdna_counter = 0
        next_step = 0
        print "Protein reading in progress"
        for q in pep:
            description = dict([z.split(":",1) for z in q.description.split()[1:] if ":" in z ])
            transcript = description["transcript"]
            gene = description["gene"]
            protein = schema.Pep(protein_id = q.id, transcript_id = transcript, protein_sequence = str(q.seq))
            db.session.add(protein)
            
            pep_dict[(gene,transcript)].append(q)
        print "Making database in progress"
        for cd in cdna:
            cdna_counter += 1
            if  (float(cdna_counter*100)/cdna_size) >= next_step:
                next_step += 10
                print str(float(cdna_counter*100)/cdna_size)+"%"
            description =  dict([z.split(":",1) for z in cd.description.split()[1:] if ":" in z])
            cdna_transcript_id = cd.id
            cdna_gene = description["gene"]
            
            gt = (cdna_gene,cdna_transcript_id)
    
            if gt in pep_dict:
                
                # TODO GO onthology
                
                for p in pep_dict[gt]:
                    
                    protein_id = p.id
                    gene_id = str(gt[0])
                    pep_sequence = str(seq3(p.seq))
                    cdna_sequence = str(cd.seq)
                    cdna_translated_list = []
                    
                    cdna_start = 0
                    cdna_stop = 0
                    
                    for x in range(3):
                        try:
                            cdna_translated_list.append(seq3(str(Seq(cdna_sequence[x:]).translate())))
                        except:
                            cdna_translated_list.append(seq3(str(Seq(cdna[x:(len(cdna_sequence)-x)-(len(cdna_sequence)-x)%3]).translate())))
                    cut_found = [v for v in range(len(cdna_translated_list)) if pep_sequence in cdna_translated_list[v]]

                    #easy
                    if pep_sequence == cdna_translated_list[0]:
                        cdna_start = 0
                        cdna_stop = len(cdna_sequence)
                        proper_seq  = cdna_sequence
                        AAA_list = findPolyA(proper_seq)
                        grab_AAA_information(AAA_list,organism_out, cdna_transcript_id,cdna_start,proper_seq)
                            
                    #cutting
                    
                    elif cut_found:
                        for c in cut_found:
                            idx = c+ cdna_translated_list[c].find(pep_sequence)
                            cdna_start = idx
                            cdna_stop = idx+len(pep_sequence)
                            proper_seq = cdna_sequence[idx:(idx+len(pep_sequence))]
                            AAA_list = findPolyA(proper_seq)
                            grab_AAA_information(AAA_list,organism_out, cdna_transcript_id,cdna_start,proper_seq)

                #alignment
                    else:
                        prot_seq = SeqRecord(Seq(seq1(pep_sequence)),id = "prot_seq")
                        y = open(organism+"prot.fasta",'w')
                        SeqIO.write(prot_seq, y, "fasta")
                        y.close()
                        best = []
                        for i in range(3):
                            w = cdna_sequence        
                            cdna_seq = SeqRecord(Seq(cdna_sequence[i:len(cdna_sequence)-((len(cdna_sequence)-i)%3)]).translate(stop_symbol="W"),id="cdna_seq") 
                            k = open(organism+"cdna.fasta",'w')
                            SeqIO.write(cdna_seq, k , "fasta")
                            k.close()
                            output = NcbiblastpCommandline(query=organism+"prot.fasta", subject=organism+"cdna.fasta", outfmt=5)()[0]
                            blast_result_records = list(NCBIXML.parse(StringIO(output)))
                            for bl_res in blast_result_records:
                                for z in bl_res.alignments:
                                    for h in z.hsps:
                                        best.append((h.query,h.sbjct,i,h.sbjct_start,h.score))
                        if best:
                            l = sorted(best,key=lambda x:x[-1])[-1]
                            proper_seq = cdna_sequence[l[2]+(int(l[3])-1)*3:l[2]+((int(l[3])-1)+len(l[1]))*3]
                            AAA_list = findPolyA(proper_seq)
                            cdna_start = l[2]+(int(l[3])-1)*3
                            cdna_stop = l[2]+((int(l[3])-1)+len(l[1]))*3
                            grab_AAA_information(AAA_list,organism_out, cdna_transcript_id,cdna_start,proper_seq)
                        else:
                            failed.write(",".join([protein_id,cdna_transcript_id,gene_id,pep_sequence,cdna_sequence])+"\n")
                        os.remove(organism+"cdna.fasta")
                        os.remove(organism+"prot.fasta")
            
                    cdna = schema.Cdna(transcript_id = cdna_transcript_id, gene_id = cdna_gene, nucleotide_sequence=str(cd.seq),organism_name =organism_id, cdna_start = cdna_start, cdna_stop =cdna_stop )
                    db.session.add(cdna)
                    
                    
                    
        db.session.commit()
            