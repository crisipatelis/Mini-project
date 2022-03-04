import subprocess
import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline


# Create log file and results directory##############
log = open('miniproject.log', 'w')

# Create a variable that will get user's working directory
admin = os.getcwd()

# Create a variable that will contain the name of the new results directory
results_directory = 'results'

# Create a path variable that will join the cwd to the new results directory
path = os.path.join(admin, results_directory)

# If the results directory does not exist, create it (done to avoid dir exists error)
if not os.path.isdir(path):
	os.mkdir(path)

# 1. ###############################################


# Create a string containing the SRA ID of the sequences
sra_id = "SRR8185310"

# Create a string prefetch that will be the concatenation of the prefetch command in SRA tools and the SRA ID
prefetch = admin + "/sratoolkit.2.11.2-ubuntu64/bin/prefetch " + sra_id

# Use the call function to run the prefetch command in the command line 
subprocess.call(prefetch, shell=True)

# Create a string fasterq_dumb that will be the concatenation of the fasterq-dump command in SRA tools and the SRA ID
fasterq_dump = admin + "/sratoolkit.2.11.2-ubuntu64/bin/fasterq-dump " + sra_id

# Use the call function to run the fasterq-dump command in the command line
subprocess.call(fasterq_dump, shell=True)


# 2. #################################################

# Create a spades_command variable that will contain a string to run SPAdes
spades_command = "python3 " +  admin + "/SPAdes-3.15.4-Linux/bin/spades.py -k 55, 77, 99, 127 -t 2 --only-assembler -s " + sra_id + ".fastq -o " + admin + "/" + "results/" + sra_id + "_assembly/"

# Write the spades command to the log file
log.write(spades_command + '\n')

# Use an os.system call to run the spades command through the command line
os.system(spades_command)


# 3. ##################################################

# Create an empty list that will contain the lengths of contigs > 1000 
bp_lengths = []

# Open the contigs fasta file that was output by the spades assembly
input_fasta = open(admin+"/results/" + sra_id + "_assembly/contigs.fasta", 'r')

# Create a new fasta file that will contain the sequences of the contigs > 1000
new_fasta = open(admin+'/results/1000_contigs.fasta', 'w')

# Create a for loop that will parse through the input fasta 
# If a sequence in the fasta is greater than 1000, write the sequence to the new fasta and append the length to the list
for record in SeqIO.parse(input_fasta, "fasta"):
	length = len(record.seq)
	if length > 1000:
		SeqIO.write(record, new_fasta, "fasta")
		bp_lengths.append(length)

# Close the new fasta file
new_fasta.close()

# Create a new variable that will determine how many contigs are > 1000
number_contigs = len(bp_lengths)

# Write the number of contigs > 1000 to the log file
log.write('There are ' +  str(number_contigs) + ' contigs > 1000 in the assembly.' + '\n')


# 4. ##################################################

# Calculate the sum of the bps of contigs > 1000
number_bp = sum(bp_lengths)

# Write the number of bps to the log file
log.write('There are ' + str(number_bp) + ' bp in the assembly.' + '\n')

# 5. ##################################################

# Create a protein_seqs variable that will contain a string to run GeneMarkS-2
protein_seqs = "perl " + admin + "/" + "gms2_linux_64/gms2.pl --SEQ " + admin + "/results/1000_contigs.fasta --genome-type bacteria --output " + admin + "/results/gms2_OUT --faa " + admin + "/results/predicted_genes"

# Write the GeneMarkS-2 protein sequence command to the outfile
log.write(protein_seqs + '\n')

# Run the GeneMarkS-2 protein sequence command in the terminal
os.system(protein_seqs)

# 6. ##################################################

# Create a string for a make_db_command that will create a database for the multi-fasta Ecoli.fasta file
make_db_command = "makeblastdb -in " + admin + "/Ecoli.fasta -out Ecoli -title Ecoli -dbtype prot"

# Write the make_db_command to the log file
log.write(make_db_command + '\n')

# Use an os.system call to run the make_db_command
os.system(make_db_command)

# Create a variable blastp_cline that will run blastp
blastp_cline = NcbiblastpCommandline(query=admin+"/results/1000_contigs.fasta", db="Ecoli", max_target_seqs=1, outfmt='"10 qseqid sseqid pident qcovs"', out=admin+'/results/predicted_functionalities.csv')

# Write the blastp_cline command to the log file
log.write(str(blastp_cline) + '\n')

# Run the blastp_cline command
stout, stderr = blastp_cline()

# 7. ##################################################

# Open the predicted_genes file output from GeneMarkS-2
infile = open(admin+'/results/predicted_genes', 'r')

# Read the file in as a string
instring = infile.read()

# Separate the instring by each new line character 
split = instring.split('\n')

# Create a count variable for the number of CDS in the predicted_genes file
count = 0

# Create a for loop that will calculate how many sequences are in the predicted_genes file
for string in split:
	if string.startswith('>'):
		count = count + 1

# Create a variable that represents the number of CDS in RefSeq for E. coli K-12 (referenced in miniproject requirements)
refseq_cds = 4140

# Calculate the difference between the GeneMarkS-2 count and the RefSeq count
difference = count - refseq_cds

# Write a string of the differences needed to be written in the log file
difference_string = "GeneMarkS found " + str(difference) + " additional CDS than the RefSeq." 

# Write the difference to the log file
log.write(difference_string + '\n')

# 8. ##################################################

# Create a string containing the SRA ID of the RNA sequences
sra_id2 = "SRR1411276"

# Create a string prefetch that will be the concatenation of the prefetch command in SRA tools and the SRA ID
prefetch = admin + "/sratoolkit.2.11.2-ubuntu64/bin/prefetch " + sra_id2

# Use the call function to run the prefetch command in the command line 
subprocess.call(prefetch, shell=True)

# Create a string fasterq_dumb that will be the concatenation of the fasterq-dump command in SRA tools and the SRA ID
fasterq_dump = admin + "/sratoolkit.2.11.2-ubuntu64/bin/fasterq-dump " + sra_id2

# Use the call function to run the fasterq-dump command in the command line
subprocess.call(fasterq_dump, shell=True)

# Create a string build_bowtie_command that will build the bowtie index as the EcoliK12 index
build_bowtie_command = "bowtie-build " + admin+"/results/" + sra_id + "_assembly/contigs.fasta EcoliK12"

# Write the build_bowtie_command to the log file
log.write(build_bowtie_command + '\n')

# Use an os.system call to run the build_bowtie_command
os.system(build_bowtie_command)

# Create a string that will run tophat
tophat_command = admin+"/tophat-2.1.1.Linux_x86_64/tophat2 -o tophat_out --no-novel-juncs EcoliK12 " + sra_id2 + ".fastq"

# Write the tophat_command to the logfile
log.write(tophat_command + '\n')

# Use an os.system call to run the tophat_command
os.system(tophat_command)

# Create a string that will run cufflinks 
cufflinks_command = admin+"/cufflinks-2.2.1.Linux_x86_64/cufflinks -o "+admin+"/results/cufflinks_out " + admin +"/tophat_out/accepted_hits.bam" 

# Write the cufflinks_command to the log file
log.write(cufflinks_command + '\n')

# Use an os.system call to run the cufflinks_command
os.system(cufflinks_command)

# 9. ##################################################

# Open the transcripts.gtf cufflinks output file
transcript_infile = open(admin+'/results/cufflinks_out/transcripts.gtf', 'r')

# Read the file in as a string
transcript_instring = transcript_infile.read()

# Remove any new line characters in the instring and replace them with tabs
transcript_instring = transcript_instring.replace('\n', '\t')

# Split the instring at each tab
transcript_split = transcript_instring.split('\t')

# Create a new list for seqnames, selecting every ninth element from the split starting at 0
seqname_list = transcript_split[0::9]

# Create a new list for start, selecting every ninth element from the split starting at 3
start_list = transcript_split[3::9]

# Create a new list for end, selecting every ninth element from the split starting at 4
end_list = transcript_split[4::9]

# Create a new list for strand, selecting every ninth element from the split starting at 6
strand_list = transcript_split[6::9]

# Create a new list for the attributes, selecting every ninth element from the split starting at 8
fpkm_attributes = transcript_split[8::9]

# Join all of the elements in the fpkm_attributes list into a string
fpkm_attr_join = ''.join(map(str,fpkm_attributes))

# Split the new string at each semicolon
fpkm_split = fpkm_attr_join.split(';')

# Create a new list to put the FPKM and value into
fpkms = []

# For each element in the strip, if it starts with FPKM (with a space in front), add it to the fpkms list
for i in fpkm_split:
	if i.startswith(" FPKM"):
		fpkms.append(i)

# Create a new list for FPKM values
fpkm_list = []

# For values in fpkms, redefine the value starting at index 7 up until the last index (this includes just the numeric value)
for val in fpkms:
	val = val[7:-1]
	fpkm_list.append(val)

# Zip each of the newly created lists
transcript_zip = list(zip(seqname_list, start_list, end_list, strand_list, fpkm_list))

# Create a new file transcriptome_data.fpkm that will output to the results directory
fpkm_csv = open(admin+'/results/transcriptome_data.fpkm', 'w')

# For each tuple in the zipped list, write each tuple element, separated by a comma, ending with a new line character to the file
for tuples in transcript_zip:
	fpkm_csv.write(tuples[0] + ',' + tuples[1] + ',' + tuples[2] + ',' + tuples[3] + ',' + tuples[4] + '\n')

# Close the file
fpkm_csv.close()

#######################################################
# Close log file
log.close()

# Move the log file into the results directory
move_file = "mv " + admin + "/miniproject.log " + admin + "/results"

# Use an os.system call to run move_file command
os.system(move_file)


