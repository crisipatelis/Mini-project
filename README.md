# Mini-project

This project was developed for COMP 483.


**Packages and Prerequisites**
------------------------------
A number of tools must be installed to run this project. The packages and their documentation are listed below:


1. SRA Tools, which can be found here: https://github.com/ncbi/sra-tools.

2. SPAdes, which can be found here: https://cab.spbu.ru/software/spades/.

3. BioPython*, which can be found here: https://www.tutorialspoint.com/biopython/biopython_installation.htm.
      
      *Note: in the development of this project, an older version of pip was used to accomodate Python 3.5 on virtual machines. 

4. GeneMarkS-2, which can be found here: http://topaz.gatech.edu/Genemark/license_download.cgi.

5. TopHat and Cufflinks, both of which can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3334321/.


**Test Data**
-------------
In its development and in order to run it, this project used two sets of test data. 
- In the first portion of the project (#1-7), Illumina sequence reads with SRA ID SRR8185310 are used for resequencing the E. coli K-12 genome. In step 1 of the script, the ```sra_id``` variable is assigned to this SRA ID. Its SRA Database page can be found here: https://www.ncbi.nlm.nih.gov/sra/SRX5005282
- In the second portion of the project (#8-9), Illumina RNA sequences with SRA ID SRR1411276 are used as transcriptome data. In step 8 of the script, the ```sra_id2``` variable is assigned to this SRA ID. Its SRA Database page can be found here: https://www.ncbi.nlm.nih.gov/sra/SRX604287


**```import``` Information**
----------------------------
At the top of the script, several packages are imported for use. They are explained below.

```import subprocess``` : Used to create processes that can be run in the command line.

```import os``` : Used to call arguments to the command line.

```from Bio import SeqIO``` : Biopython's Sequence Input/Output interface.

```from Bio.Blast import NCBIWWW``` : Biopython's access to WWW version of NCBI's BLAST.

```from Bio.Blast.Applications import NcbiblastpCommandline``` : Biopython's interface to run blastp from command line.

**Parameter and Step Information**
----------------------------------
