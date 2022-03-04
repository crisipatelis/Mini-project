# Mini-project

This project was developed for COMP 483. It was created using a Virtual Machine with an Ubuntu OS.  


**Packages and Prerequisites**
------------------------------
A number of tools must be installed to run this project. The packages and their documentation are listed below:


1. SRA Tools, which can be found here: https://github.com/ncbi/sra-tools.

2. SPAdes, which can be found here: https://cab.spbu.ru/software/spades/.

3. BioPython*, which can be found here: https://www.tutorialspoint.com/biopython/biopython_installation.htm.
      
      *Note: in the development of this project, an older version of pip was used to accomodate Python 3.5 on virtual machines. 

4. GeneMarkS-2, which can be found here: http://topaz.gatech.edu/Genemark/license_download.cgi.

5. TopHat and Cufflinks, both of which can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3334321/.

The script runs under the assumption that packages are installed to the user's home directory. Otherwise, paths must be specified. Additionally, if packages are downloaded under different operating systems (other than Ubuntu appropriate versions), The paths to access each packages programs must be specified. For example, 

**Test Data**
-------------
In its development and in order to run it, this project used two sets of test data. The script is designed to retrieve these reads. 
- In the first portion of the project (#1-7), Illumina sequence reads with SRA ID SRR8185310 are used for resequencing the E. coli K-12 genome. In step 1 of the script, the ```sra_id``` variable is assigned to this SRA ID. Its SRA Database page can be found here: https://www.ncbi.nlm.nih.gov/sra/SRX5005282
- In the second portion of the project (#8-9), Illumina RNA sequences with SRA ID SRR1411276 are used as transcriptome data. In step 8 of the script, the ```sra_id2``` variable is assigned to this SRA ID. Its SRA Database page can be found here: https://www.ncbi.nlm.nih.gov/sra/SRX604287
- In step 6, a given multi-fasta file was used to create an E. coli K-12 database as subject sequences. It was obtained at /home/aene/Ecoli.fasta. 


**```import``` Information**
----------------------------
At the top of the script, several packages are imported for use. They are explained below.

```import subprocess``` : Used to create processes that can be run in the command line.

```import os``` : Used to call arguments to the command line.

```from Bio import SeqIO``` : Biopython's Sequence Input/Output interface.

```from Bio.Blast import NCBIWWW``` : Biopython's access to WWW version of NCBI's BLAST.

```from Bio.Blast.Applications import NcbiblastpCommandline``` : Biopython's interface to run blastp from command line.

**Parameter and Important Step Information**
------------------------------------------
#1.

```admin = os.getcwd()``` : this is a variable used many times throughout the script to begin paths. Have the script in the user's home directory so that the os call will get that current directory.  

```sra_id``` : test data information described in **Test Data** section above. 

#2. 

```spades_command```:
- ```-k 55, 77, 99, 127``` : Run SPAdes for _de novo_ assembly.
- ```-t 2``` : use two processors.
- ```--only-assembler``` : only do the assembly.
- ```-s " + sra_id + ".fastq```: one input file, just the sra_id.fastq file (for single end reads).
- ```-o " +admin+"/results/SRR8185310_assembly/```: the directory the spades command will output. 

#5. 

```protein_seqs```:
- ```--SEQ " + admin + "/results/1000_contigs.fasta```: genome sequences for contigs > 1000 as input.
- ```--genome-type bacteria```: is a bacterial genome.
- ```--output " + admin + "/results/gms2_OUT```: path to output the output file.
- ```--faa " + admin + "/results/predicted_genes"```: path to output the file that will contain predicted gene sequences.

#6. 

```make_db_command```: 
- ```makeblastdb```: invoke makeblastdb program.
- ```-in " + admin + "/Ecoli.fasta```: database input path and file.
- ```-out Ecoli```: database output name (outputs to home directory).
- ```-title Ecoli```: Ecoli is the title of this database.
- ```-dbtype prot```: protein database type.

```blastp_cline```: 
- ```query=admin+"/results/1000_contigs.fasta"```: the query sequence to enter into blastp.
- ```db="Ecoli"```: name of database. 
- ```max_target_seqs=1```: includes only the best hits.
- ```outfmt='"10 qseqid sseqid pident qcovs"'```: format of the output file (query sequence id, subject sequence id, percent identity, percent query coverage).
- ```out=admin+'/results/predicted_functionalities.csv'```: path to the output file (outputs in results directory).

```stout, stderr```: runs ```blastp_cline``` command while outputting errors (if reference to them is needed).

#8. 

```build_bowtie_command```:
- ```bowtie-build```: invokes the bowtie-build program.
- ```admin+"/results/" + sra_id + "_assembly/contigs.fasta```: input file to build bowtie index.
- ```EcoliK12```: name of the index.

```tophat_command```: 
- ```admin+"/tophat-2.1.1.Linux_x86_64/tophat2```: invokes tophat2 program.
- ```-o tophat_out```: name and location of tophat output (outputs to home directory).
- ```--no-novel-juncs```: does not include novel splice site discovery.
- ```EcoliK12```: bowtie index.
- ```sra_id2 + ".fastq"```: transcriptome fastq file.

```cufflinks_command```:
- ```admin+"/cufflinks-2.2.1.Linux_x86_64/cufflinks```: invoke cufflinks program
- ```-o "+admin+"/results/cufflinks_out "```: name and location of cufflinks output (outputs to results directory)
- ```+ admin +"/tophat_out/accepted_hits.bam"```: bam file from tophat output

 For more information on parameter options for each step/package, search the documentation at the links provided in the **Packages and Prerequisites** section.

**Results**
-----------

After running the script (with the test data), the following files can be found in the results directory:

```1000_contigs.fasta``` ```cufflinks_out``` ```gms2_OUT``` ```predicted_functionalities.csv``` ```predicted_genes``` ```SRR8185310_assembly``` ```transcriptome_data.fpkm``` ```miniproject.log```
