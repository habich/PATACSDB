#!/usr/bin/python
import glob
import sys
import os
import time
import logging
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
import urllib2
import xmltodict
import subprocess

#ensembl_request = "http://rest.ensemblgenomes.org/lookup/id/"
ensembl_request = "http://rest.ensembl.org/lookup/id/"
biomart_database = "http://central.biomart.org/martservice/datasets?config=gene_ensembl_report"
AAA_id = 0
GO_id = 0
logger = logging.getLogger(__name__)
FAST_NAMES= True


def configure_logging(level = logging.INFO):
    logger.setLevel(level)
    ch = logging.StreamHandler()
    logger.addHandler(ch)

# Finds available organisms in Biomart database 
def get_biomart_database():
    biomart_dict = {}
    response = urllib2.urlopen(biomart_database)
    xml = response.read()
    obj = xmltodict.parse(xml)
    for dataset in obj['datasets']['dataset']:
        organism_name = dataset['@displayName'].split(' genes')[0]
        database = dataset['@name']
        biomart_dict[organism_name] = database
    return biomart_dict

def get_ensembl_databases():
    dic = {}
    i = 0
    raw_list = ""
    databases = ["mysql-eg-publicsql.ebi.ac.uk --port 4157","ensembldb.ensembl.org --port 3306"]
    for database in databases:
        while i < 5:
            try:
                sp = subprocess.check_output("mysql -h "+database+" --user=anonymous < get_databases.sql",
                                    shell = True)
                raw_list = sp
                break
            except Exception as e:
                i+= 1
                print e
                print sys.exc_info()[0]
                print "database has problem: ",database
                time.sleep(1)
        if raw_list:
            dic_list = defaultdict(list)
            for l in raw_list.split():
                if "_core_" in l:
                    (name,number) = l.split("_core_")
                    dic_list[name].append((number.split("_")[0],number))
            for x in dic_list.keys():
                dic_list[x] = x+"_core_"+max(dic_list[x])[1]
            dic[database]= dic_list
    return dic

def organisms_files_dictionary(raw_data_directory):
    dic_func = lambda :{"pep":"","gtf":"","cdna":"","database":""}
    org_dict = defaultdict(dic_func)
    for directory,_,files in os.walk(raw_data_directory):
        for f in files:
            name = f.split(".")
            if "abinitio" in name:
                continue
            organism = f.split("/")[-1].split(".")[0]
            if organism:
                if name[-3] == "pep":
                    org_dict[organism]["pep"]= directory+"/"+f
                elif name[-3] == "cdna":
                    org_dict[organism]["cdna"]= directory+"/"+f
                elif name[-1] == "gtf":
                    org_dict[organism]["gtf"]= directory+"/"+f
                database_looking = directory.split("/")
                if "metazoa" in database_looking:
                    org_dict[organism]["database"]="metazoa"
                elif "protists" in database_looking:
                    org_dict[organism]["database"]="protists"
    return org_dict

def load_gtf(ensembl_names,paths,organism):
    exon_places = defaultdict(list)
    prot_names = get_protein_names(organism.lower(),ensembl_names)
    gtf_file = open(paths["gtf"],'r')
    for gtf in gtf_file:
        if gtf[0]=="#":
            continue
        else:
            chrom, source, feature, start, stop, score,strand,frame,attribute = gtf.split("\t")
            info = dict([x.split() for x in attribute.split(";") if len(x.split()) == 2])
            if feature == "transcript":
                transcript_id = info["transcript_id"][1:-1]
                if transcript_id in prot_names:
                    protein_name = prot_names[transcript_id]
                elif FAST_NAMES == True:
                    protein_name = None
                else:
                    protein_name = name_import(transcript_id,organism)

                qtf_feature = schema.Gtf(transcript_id = transcript_id, chromosome_name = chrom, chromosome_location_start = start, chromosome_location_stop = stop, strand = strand, protein_name = protein_name)
                db.session.add(qtf_feature)
            elif feature == "exon":
                transcript_id = info["transcript_id"][1:-1]
                exon_places[transcript_id].append((int(start),int(stop)+1,strand))
            else:
                continue
    for tid in exon_places:
        exon_places[tid] = translate_exon_locations(exon_places[tid])
    return exon_places

def translate_exon_locations(exon_list):
    l = []
    if exon_list[0][2] == "-":
        l.append("-")
    else:
        l.append("+")
    current_end = 0
    current_start = 0
    for p in sorted(exon_list):
        current_end = current_start+ p[1]-p[0]
        local_position = (current_start,current_end,p[0])
        l.append(local_position)
        current_start = current_end    
    return l

def name_import(transcript_id,organism):
    logger.debug("Making request " + ensembl_request+transcript_id+"?")
    r = requests.get(ensembl_request+transcript_id+"?", headers={ "Content-Type" : "application/json", "species":organism.lower()})
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

def get_protein_names(organism,dic):
    for database in dic:
        i = 0
        if organism in dic[database]:
            while i < 5:
                try:
                    g = subprocess.check_output("mysql -h "+database+" --user=anonymous "+dic[database][organism]+" < get_names.sql",
                                            shell = True)
                    break
                except Exception as e:
                    i+= 1
                    print e
                    print sys.exc_info()[0]
                    print "getting protein names attempt: ",i,"/5"
            if i >= 5:
                print "Getting protein names failed"
                return {}
            return  dict([(x.split("\t")[0],x.split("\t")[1]) for x in g.split("\n") if x])
    return {}

def load_pep_dict(paths,organism):
    pep = SeqIO.parse(open(paths["pep"],'r'),"fasta")
    pep_dict = defaultdict(list)
    for q in pep:
        description = dict([z.split(":",1) for z in q.description.split()[1:] if ":" in z ])
        transcript = description["transcript"]
        gene = description["gene"]
        pep_dict[(gene,transcript)].append(q)
        
        protein = schema.Pep(protein_id = q.id, transcript_id = transcript, protein_sequence = str(q.seq))
        db.session.add(protein)
    return pep_dict

def load_species(biomart_database, paths,organism):
    species_name = " ".join(organism.split("_"))
    if species_name in biomart_database:
        biomart = biomart_database[species_name]
    else:
        biomart = None
    
    species = schema.Species(id = organism,database=paths["database"], species_name = species_name, biomart_database=biomart)
    db.session.add(species)

def load_cdna_and_polyA(paths,organism,pep_dict,exon_info):
    organism_out = open(organism+"_polyA.data",'w')
    failed = open(organism+"_failed",'w')
    cdna = list(SeqIO.parse(open(paths["cdna"],'r'),"fasta"))
    cdna_size = len(cdna)
    cdna_counter = 0
    next_step = 0
    for cd in cdna:
        cdna_counter += 1
        if  (float(cdna_counter*100)/cdna_size) >= next_step:
            next_step += 10
            print str(int(float(cdna_counter*100)/cdna_size))+"%"
        description =  dict([z.split(":",1) for z in cd.description.split()[1:] if ":" in z])
        cdna_transcript_id = cd.id
        cdna_gene = description["gene"]
        
        gt = (cdna_gene,cdna_transcript_id)
        
        if gt in pep_dict:
            
            for p in pep_dict[gt]:
                
                protein_id = p.id
                gene_id = str(gt[0])
                pep_sequence = str(seq3(p.seq))
                cdna_sequence = str(cd.seq)
                cdna_translated_list = []
                
                cdna_start = 0
                cdna_stop = 0
                
                for x in range(3):
                    cdna_translated_list.append(seq3(str(Seq(cdna_sequence[x:]+"N"*(3-len(cdna_sequence[x:])%3)).translate())))
                cut_found = [v for v in range(len(cdna_translated_list)) if pep_sequence in cdna_translated_list[v]]
                
                #easy
                if pep_sequence == cdna_translated_list[0]:
                    cdna_start = 0
                    cdna_stop = len(cdna_sequence)
                    proper_seq  = cdna_sequence
                    AAA_list = findPolyA(proper_seq)
                    grab_AAA_information(AAA_list,organism_out, cdna_transcript_id,cdna_start,proper_seq,exon_info[cdna_transcript_id])
                        
                #cutting
                elif cut_found:
                    for c in cut_found:
                        idx = c+ cdna_translated_list[c].find(pep_sequence)
                        cdna_start = idx
                        cdna_stop = idx+len(pep_sequence)
                        proper_seq = cdna_sequence[idx:(idx+len(pep_sequence))]
                        AAA_list = findPolyA(proper_seq)
                        grab_AAA_information(AAA_list,organism_out, cdna_transcript_id,cdna_start,proper_seq,exon_info[cdna_transcript_id])
                        
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
                        grab_AAA_information(AAA_list,organism_out, cdna_transcript_id,cdna_start,proper_seq,exon_info[cdna_transcript_id])
                    else:
                        failed.write(",".join([protein_id,cdna_transcript_id,gene_id,pep_sequence,cdna_sequence])+"\n")
                    os.remove(organism+"cdna.fasta")
                    os.remove(organism+"prot.fasta")
        
                cdna = schema.Cdna(transcript_id = cdna_transcript_id, gene_id = cdna_gene, nucleotide_sequence=str(cd.seq),organism_name =organism, cdna_start = cdna_start, cdna_stop =cdna_stop)
                db.session.add(cdna)

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

def grab_AAA_information(AAA_list,organism_out, cdna_transcript_id,cdna_start,proper_seq,exon_info):
    if exon_info[0] == "+":
        exon_sum = 0
    else:
        exon_sum = sum([y-x for x,y,z in exon_info[1:]])
    for AAA in AAA_list:
        AAA_start = abs(exon_sum - (cdna_start + AAA[0]))
        AAA_stop = abs(exon_sum - (cdna_start+AAA[0]+AAA[1]))
        if AAA_start > AAA_stop:
            AAA_start,AAA_stop = AAA_stop, AAA_start
        AAA_global_start = -1
        AAA_global_stop = -1
        assert exon_info
        for p in exon_info[1:]:
            if AAA_start >= p[0] and AAA_start <= p[1]:
                AAA_global_start = p[2]+(AAA_start-p[0])
            if AAA_stop >= p[0] and AAA_stop <= p[1]:
                AAA_global_stop = p[2]+(AAA_stop-p[0])
            if AAA_global_start > -1 and AAA_global_stop > -1:
                break
        AAA_global_stop -= 1
        AAA_start = cdna_start + AAA[0]
        AAA_stop = cdna_start+AAA[0]+AAA[1]
        database_append(AAA_start ,AAA_stop ,cdna_transcript_id,AAA_global_start,AAA_global_stop)
        out_write(organism_out, cdna_transcript_id,AAA_start,AAA_stop,AAA_global_start,AAA_global_stop,proper_seq,AAA)

def database_append(AAA_start,AAA_stop,cdna_transcript_id,AAA_global_start,AAA_global_stop):
        global AAA_id
        AAA_id += 1
        polyA = schema.PolyA(id = AAA_id, transcript_id = cdna_transcript_id, AAA_start = AAA_start, AAA_stop = AAA_stop,AAA_global_start=AAA_global_start, AAA_global_stop=AAA_global_stop)
        db.session.add(polyA)

def out_write(organism_out, cdna_transcript_id,AAA_start,AAA_stop,AAA_global_start,AAA_global_stop,proper_seq,AAA):
    organism_out.write("\t".join([cdna_transcript_id,str(AAA_start),str(AAA_stop),str(AAA_global_start),str(AAA_global_stop), proper_seq[AAA[0]:AAA[0]+AAA[1]],"\n"]))

def main():
    configure_logging()
    print "Hello!"
    if len(sys.argv) != 2:
        print "Usage: "+sys.argv[0]+" <directory with ensembl data>"
        exit(1)
    
    raw_data_directory = sys.argv[1]
    db.drop_all()
    db.create_all()
    
    org_dict = organisms_files_dictionary(raw_data_directory)
    biomart_database = get_biomart_database()
    ensembl_names = get_ensembl_databases()
    
    for organism,paths in org_dict.items():
        print "Organism: ", organism
        load_species(biomart_database,paths,organism)

	print "Paths: ", paths
	        
        print "GTF reading in progress"
        exon_info = load_gtf(ensembl_names,paths,organism)

        print "Protein reading in progress"
        pep_dict = load_pep_dict(paths,organism)
        
        print "Making database in progress"
        load_cdna_and_polyA(paths,organism,pep_dict,exon_info)
                    
        db.session.commit()

if __name__ == "__main__":
    main()
