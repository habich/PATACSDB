#!/bin/bash
rmdir -rf raw_ensembl_data
mkdir raw_ensembl_data
cd raw_ensembl_data
mkdir standard protists metazoa

# standard ensembl
cd standard
wget --mirror -A pep.all.fa.gz ftp://ftp.ensembl.org/pub/current_fasta/
wget --mirror -A cdna.all.fa.gz ftp://ftp.ensembl.org/pub/current_fasta/
wget --mirror -A gtf.gz ftp://ftp.ensembl.org/pub/current_gtf/

# protists ensembl
cd ../protists
wget --mirror -A pep.all.fa.gz ftp://ftp.ensemblgenomes.org/pub/protists/current/fasta/
wget --mirror -A cdna.all.fa.gz ftp://ftp.ensemblgenomes.org/pub/protists/current/fasta/
wget --mirror -A gtf.gz ftp://ftp.ensemblgenomes.org/pub/protists/current/gtf/

# metazoa ensembl
cd ../metazoa
wget --mirror -A pep.all.fa.gz ftp://ftp.ensemblgenomes.org/pub/metazoa/current/fasta/
wget --mirror -A cdna.all.fa.gz ftp://ftp.ensemblgenomes.org/pub/metazoa/current/fasta/
wget --mirror -A gtf.gz ftp://ftp.ensemblgenomes.org/pub/metazoa/current/gtf/

cd ..
find -name *gz |xargs gunzip
