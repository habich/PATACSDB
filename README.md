Install: SQLAlchemy, Flask, Flask-SQLAlchemy , xmltodict, mysql, biopython, ncbi-blast+

Run download.sh to download data from ensemble database. It will destroy existing data (previous versions) and download current version (even if there is no change). 
Run parse_data.py to analyse existing data about organisms, create files with informations about founded polyA (not required for database) and to initialize database.
! Warning: parse_data.py destroys database content.
Run runserver.py to run the server.

Details:

parse_database.py:
Input - Directory with raw ensembl data (pep, cdna, gtf)
Output 
Database with polyA data
Files with polyA data
Information about failed operations


 get_biomart_database

We use links to biomart database to provide comprehensive information about genes. If Biomart does not work or organism is not in biomart database, the link is not produced.
libraries: xmltodict, urllib2
used functions: None
input: None
output: dictionary <organism_name> : name of dataset in biomart

get_ensembl_databases

We use ensembl database server to get proteins names. There are two sources of information: ensembldb.ensembl.org and mysql-eg-publicsql.ebi.ac.uk . To connect with server we use local mysql script (get_databases.sql) and analyze response. If interaction did not succeeded, exception is raised. After five futile attempts the effort is aborted. From available datasets we choose the newest core information.
libraries: mysql, time, sys, collections
used functions: None
input: *get_databases.sql
output: dictionary <organism_name> : name of ensembl dataset

organisms_files_dictionary

After downloading raw data from ensembl server, datafiles are structured in a specific way. Using information about this structure we extract access paths to datafiles. There are three categories of datafiles: pep, gtf and cdna. As the names suggest files contain information about peptide sequences, standard gtf files, and cdna sequences. Additionally information about general dataset is remembered so we can distinguish between Metazoa, Protists, and standard ensembl database. We ignore abinito files.
libraries: collections
Used functions: None
input: raw_data_directory -> place where raw data is keeped
output: dictionary <organism_name> : { "pep" : <path to pep file> , "gtf" : <path to gtf file> , "cdna" : <path to cdna file> , "database": <type : "metazoa" | "protists" | "" (standard> }

load_gtf

From gtf files we extract information about transcripts and exact exon locations. Information about transcripts are stored in database in class Gtf. Remembered information is: transcript id, chromosome name, location in chromosome, strand and protein name. Protein name is extracted from ensembl server (-> get_ensembl_database) if possible or from ensembl rest service if constant FAST_NAMES is set on False (default: True). Warning: rest service is very slow. 
For every transcript information about exons location and strand is returned.
libraries: collections
used functions: get_protein_names, translate_exon_locations
input: ensembl_names (get_biomart_database output), paths (organisms_files_dictionary output values), organism (organisms_files_dictionary output keys)
output: **schema.Gtf data , dictionary: <transcript id> : [list of: (exon start, stop, global chromosome location start)] -> more information in translate_exon_locations

translate_exon_locations

Locations extracted from gtf files are global i.e. those are location in chromosome. To find polyA we also need local positions - where exons are counted from the beginning of cdna.
libraries: None
used functions: None
input: exon list -> metadata from load_gtf : [list of: (global exon start, stop, strand) ]
output: [list of: "+"|"-" (strand), (local exon start, stop, global start)]

name_import

If ensemble server does not provide information about protein names and the "FAST_NAMES" constant is set to False (FAST_NAMES=True at the beginning of parse_data.py; used in load_gtf), this function uses ensembl rest service to obtain protein name. Warning: It is very slow and unreliable.
libraries: requests
used functions: None
input: transcript_id, organism (organism name), *ensembl_request -> global constant
output: protein name

get_protein_names

This function uses ensembl databases found by get_ensembl_databases to make inquiry to ensemble to get information about protein names. Using mysql and local script (get_names.sql) program makes five attempts to get information. If all attempts fail method returns empty dictionary. When successful, produced output contains known popular names for proteins.
Note: if process of obtaining name failed, a protein will still have value in "protein name" column. The value will be the standard ensembl code for this protein for example: ENSACAP00000000394 .
libraries: subprocess, mysql, sys
input: organism (organism name) , dic (dictionary of ensembl databases; output of get_ensembl_databases)
output: dictionary: <transcript id> : protein name

load_pep_dict

In this part of code we extract information about protein sequences from fasta files. Information is stored in database in class named Pep and returned from the function. If gtf file did not contain information about specific transcript associated with this protein, information about this protein is ignored due to impossibility of further analysis. This event is reported in "failed" file. In database we store information about protein id, associated transcript id and protein sequence.
libraries: Bio, collections
used functions: None
input:  paths (organisms_files_dictionary output values), organism (organisms_files_dictionary output keys), exon_info (output of load_gtf) , failed (file with failed examples)
output: **schema.Pep, pep_dict -> dictionary: <(gene_id, transcript_id)> : Bio fasta object

load_species

For proper display information about species on web service, some information about species have to be stored in database. In class Species we keep organism id which is organism name using just small letters and with spaces replaced by "_", species name - organism name in the form to display and information about availability of biomart database.
libraries: None
input : bimart_database (get_biomart_database output), paths (organisms_files_dictionary output values), organism (organisms_files_dictionary output keys)
output : **schema.Species

load_cdna_and_polyA

The most time-consuming process during preparation of database takes place in here. In general corresponding sequences of cdna and pep are paired according to their ids and actual sequences and then polyA fragments are found and stored.
In the first part of code cdna sequences are read from cdna file. For each sequence an attempt is made to pair it to stored pep data. Pairing is made by comparison of gene and transcript ids in fasta header. If pairing fails the attempt is reported in "failed" file. Otherwise we try to pair sequences to obtain corresponding sequences. 
Note: Some sequences in ensembl database are not their exact translations. Sometimes there are frameshifts or additional nucleotides at the end of sequences. However there are also cases where trivial attempts to match fail and in those cases we use blast.
In the next part of code we pair sequences. We use three methods with increasing time consumption.
The first attempt is to try if translated cdna sequence is exact copy of pep sequence. If it is true, the information is stored.
Second attempt is to check if in any frame we can find the protein sequence somewhere in translated cdna sequence. All three frames are checked. Information about changed frame and place where protein sequence is, is stored.
The last attempt is using blast. For all three frames we store blast best matches, and then we choose the best from the best. In this case both protein and cdna sequence can be shortened, and there is a chance of mismatches. Using Blast requires making temporary files named by us "cdna.fasta" and "prot.fasta". After every protein analysis we remove those files. 
After every attempt, if it succeeded, we search for polyA fragments and store information about that. You can find more about that process in: findPolyA, grab_AAA_information.
In the end in database we store information in class Cdna - transcript id, gene id, cdna sequence, organism name, start and stop of cdna sequence (necessary id cdna was cut during pairing process)
libraries: Bio, blast, os
used functions: findPolyA, grab_AAA_information
input:  paths (organisms_files_dictionary output values), organism (organisms_files_dictionary output keys), pep_dict (load_pep_dict output),  exon_info (output of load_gtf) , failed (file with failed examples)
output: **schema.Cdna

findPolyA

If you want to change pattern that you are looking for, this is the place to change. FindPolyA as name suggests finds polyA fragments. Length of polyA fragment can be changed in the first line of this function. Here are some examples of typical problems and answers:
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
As a result function gives ranges between which you can find polyA sequence.
libraries: None
used functions: None
input : sequence (in fact any string)
output: list of: ranges where you can find polyA [start,stop]

grab_AAA_information

Depending on the strand global locations of exons are increasing or falling. In order to find global location for found polyA fragments we need to calculate relative position of polyA fragment on exons and from that on chromosome. After that information about location is transferred to output files and database_append function.
libraries: None
used functions: database_append, out_write
input: AAA_list (findPolyA output) ,organism_out (output file) , cdna_transcript_id (transcript id) ,cdna_start (place where cdna start -> more in load_cdna_and_polyA) ,proper_seq (cutted cdna sequence -> more in load_cdna_and_polyA) , exon_info (load_gtf output)
output: None

database_append

This function uses information obtained by grab_AAA_information and adds data to database, class PolyA. It also use global variable "AAA_id" which counts found polyA fragments.
libraries: None
used functions: None
input: AAA_start (local polyA start) ,AAA_stop (local polyA stop) ,cdna_transcript_id (transcript id) ,AAA_global_start (global polyA start) ,AAA_global_stop (global polyA stop)
output: ** schema.PolyA

out_write

Write output file for polyA fragments. PolyA file is not used in any other place. Those files (in Data directory) are just for user application.
libraries: None
used functions: None
input: organism_out (output file) , cdna_transcript_id (transcript id) ,AAA_start (polyA local start) ,AAA_stop (polyA local stop) ,AAA_global_start (polyA global start) ,AAA_global_stop (polyA global stop) ,proper_seq (cdna sequence used to find polyA) ,AAA (polyA sequence)
output: *file

main

Coordination and time optimizations.
In the first part of code you can find interaction with the user. Then database is created and all universal datasets are obtained. The directory with raw ensembl data is necessary to start the program. List of organisms is made based on directory content.
In the next part, for every organism on the list, if organism was not analyzed before, following steps are made: load information about species (load_species), read gtf file (load_gtf), read protein file (load_pep_dict), analyse organism data (load_cdna_and_polyA). After every organism we do a database commit.
Some information exposed in the web service would require time-consuming queries to database. Therefore after all organisms are in database we store some metadata to make our service faster. 
