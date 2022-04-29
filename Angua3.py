import os
import sys
import subprocess
import multiprocessing
import argparse
from Bio import SeqIO

def main():

	# Angua3 script version
	angua3_version = 3

	### Input arguments
	options = parseArguments()

### Create project folder 
	if(os.path.isdir(options.output)):
		print("Project folder already created.")
	else:
		print("Creating project folder.")
		os.mkdir(options.output)

### Intial directory setup
	if(options.create_dirs == "Y"):
		if(options.trimmer == "sickle"):
			os.makedirs(f"{options.output}/sickle/unpaired")
		if(options.trimmer == "bbduk"):
			os.mkdir(f"{options.output}/bbduk")
		os.mkdir(f"{options.output}/trinity")
		os.mkdir(f"{options.output}/blastn")
		os.mkdir(f"{options.output}/blastx")
		os.mkdir(f"{options.output}/mmseq2")
		os.makedirs(f"{options.output}/contigs/200")
		os.makedirs(f"{options.output}/contigs/1000")
		os.makedirs(f"{options.output}/megan/blastn")
		os.makedirs(f"{options.output}/megan/blastx")
		os.mkdir(f"{options.output}/unmapped_reads")

	else:
		print("Directory structure creation skipped.")

### Run trimming
	if(options.trimmer != "N"):
		print(f"Beginning trimming with {options.trimmer}.")
		
		# Single end
		if(options.single_end == "Y"):
			for file in os.listdir(f"{options.input}"):
				sample_name = file.split(".")[0]

				if(os.path.isfile(f"{options.input}/{file}")):
					trimmer_input = f"{options.input}/{file}"
					trimmer_output = f"{options.output}/{options.trimmer}/{file.replace('_R1_001', '_R1')}"

					if(options.trimmer == "bbduk"):
						subprocess.call(f"bbduk.sh in={trimmer_input} out={trimmer_output} minlen={options.bbduk_minl} ktrim=r k=23 mink=11 hdist=1 ref={options.bbduk_adapters} qtrim=r trimq={options.bbduk_q}", shell = True)
					if(options.trimmer == "sickle"):
						subprocess.call(f"sickle se -f {trimmer_input} -t sanger -q {options.sickle_q} -l {options.sickle_minl} -o {trimmer_output} -x -g", shell = True)

		# Paired end
		if(options.single_end == "N"):
			for file in os.listdir(f"{options.input}"):
				if("_R1" in file):
					file_R2 = file.replace("_R1", "_R2")
					sample_name_R1 = file.split(".")[0]
					sample_name_R2 = sample_name_R1.replace("_R1", "_R2")

					trimmer_input_R1 = f"{options.input}/{file}"
					trimmer_input_R2 = f"{options.input}/{file_R2}"
					trimmer_output_R1 = f"{options.output}/{options.trimmer}/{file.replace('_L001_R1_001', '_R1')}"
					trimmer_output_R2 =f"{options.output}/{options.trimmer}/{file_R2.replace('_L001_R2_001', '_R2')}"
					trimmer_output_single =f"{options.output}/{options.trimmer}/unpaired/{file}"

					if(options.trimmer == "bbduk"):
						subprocess.call(f"bbduk.sh in1={trimmer_input_R1} in2={trimmer_input_R2} out1={trimmer_output_R1} out2={trimmer_output_R2} minlen={options.bbduk_minl} ktrim=r k=23 mink=11 hdist=1 ref={options.bbduk_adapters} qtrim=r trimq={options.bbduk_q}", shell = True)
					if(options.trimmer == "sickle"):
						subprocess.call(f"sickle pe -f {trimmer_input_R1} -r {trimmer_input_R2} -t sanger -q {options.sickle_q} -l {options.sickle_minl} -o {trimmer_output_R1} -p {trimmer_output_R2} -s {trimmer_output_single} -x -g", shell = True)
	else:
		print("Trimming skipped.")

### Run assembly

	if(options.trinity == "Y"):
		print(f"Beginning assembly with Trinity.")

		# Single end
		if(options.single_end == "Y"):
			for file in os.listdir(f"{options.output}/{options.trimmer}"):
				sample_name = file.split (".")[0]

				if(os.path.isfile(f"{options.input}/{file}")):
					trinity_input = f"{options.input}/{file}"
					trinity_output = f"{options.output}/trinity/{sample_name.replace('_R1', '_trinity')}"
					trinity_log = f"{options.output}/trinity/{sample_name.replace('_R1', '.log')}"

					subprocess.call(f"Trinity --seqType fq --max_memory {options.trinity_mem} --single {trinity_input} --CPU {options.trinity_cpu} --full_cleanup --output {trinity_output} > {trinity_log}", shell = True)

		# Paired end
		if(options.single_end == "N"):
			for file in os.listdir(f"{options.output}/{options.trimmer}"):
				if("_R1" in file):
					file_R2 = file.replace("_R1", "_R2")
					sample_name_R1 = file.split (".")[0]
					sample_name_R2 = sample_name_R1.replace("_R1", "_R2")

					trinity_input_R1 = f"{options.output}/{options.trimmer}/{file}"
					trinity_input_R2 = f"{options.output}/{options.trimmer}/{file_R2}"
					trinity_output = f"{options.output}/trinity/{sample_name_R1.replace('_R1', '_trinity')}"
					trinity_log = f"{options.output}/trinity/{sample_name_R1.replace('_R1', '.log')}"

					subprocess.call(f"Trinity --seqType fq --max_memory {options.trinity_mem} --left {trinity_input_R1} --right {trinity_input_R2} --CPU {options.trinity_cpu} --full_cleanup --output {trinity_output} > {trinity_log}", shell = True)

	else:
		print("Assembly skipped.")

### Sort and rename contigs

	if(options.sort == "Y"):
		print(f"Beginning sorting contigs into bins of 200 and 1000.")
		for file in os.listdir(f"{options.output}/trinity"):
			if(file.endswith(".fasta")):
				sample_name = file.split("_trinity")[0]
				with open(f"{options.output}/contigs/200/sorted_200_{sample_name}.fasta", "w") as c200, open(f"{options.output}/contigs/1000/sorted_1000_{sample_name}.fasta", "w") as c1000:
					for seq_record in SeqIO.parse(open(f"{options.output}/trinity/{file}", mode = 'r'), 'fasta'):
						seq_record.id = f"{options.output.split('/')[0]}_{sample_name.split('_')[-1]}_{seq_record.id}"
						seq_record.description = f"{options.output.split('/')[0]}_{sample_name.split('_')[-1]}_{seq_record.description}"
						if(len(seq_record.seq) >= 200):
							SeqIO.write(seq_record, c200, 'fasta')
							if(len(seq_record.seq) >= 1000):
								SeqIO.write(seq_record, c1000, 'fasta')

	else:
		print("Sorting skipped.")

### Cluster contigs for blastx

	if(options.cluster == "Y"):
		print(f"Beginning clustering of contigs >= 1000.")
		for file in os.listdir(f"{options.output}/contigs/1000/"):
			sample_name = file.split(".fasta")[0]
			subprocess.call(f"mmseqs easy-cluster -c {options.cluster_perc} --threads {options.cluster_threads} -v 0 {options.output}/contigs/1000/{file} {options.output}/mmseq2/{sample_name}.fasta {options.output}/mmseq2/tmp", shell = True)
			os.remove(f"{options.output}/mmseq2/{file}_all_seqs.fasta")
			os.rename(f"{options.output}/mmseq2/{file}_rep_seq.fasta", f"{options.output}/mmseq2/{sample_name}_rep_seq.fasta")
			os.rename(f"{options.output}/mmseq2/{file}_cluster.tsv", f"{options.output}/mmseq2/{sample_name}_cluster.tsv")

	else:
		print("Clustering skipped.")

### Run Blastn

	# Blastn seems to work better with multiple processes instead of multiple threads.
	# As such, a pool was create that can hold 8 processes at once. Once a process has been completed, a new process will be added until everything queued has been completed.  
	# args correspond to the runBlast function

	if(options.blastn == "Y"):
		print("Beginning Blastn.")
		pool = multiprocessing.Pool(processes = int(options.blastn_pool))
		results = [pool.apply_async(runBlast, args = ("blastn", "megablast", options.nt_db, f"{options.output}/contigs/200/{file}", options.blastn_threads, options.blastn_descriptions, options.blastn_alignments, f"{options.output}/blastn/{file.split('.')[0]}.megablast.blastn")) for file in sorted(os.listdir(f"{options.output}/contigs/200/"))]

		for p in results:
			p.get()
	else:
		print("Blastn skipped.")

### Run Megan

	# A pool was created that can hold 2 processes at once. Once a process has been completed, a new process will be added until everything quened has been completed. 
	# args correspond to the runMegan function

	if(options.megan_blastn == "Y"):
		print("Beginning Megan for Blastn files.")
		pool = multiprocessing.Pool(processes = int(options.megan_processes))
		results = [pool.apply_async(runMegan, args = (f"{options.output}/blastn/{file}", "BlastN", f"{options.output}/contigs/200/{file.split('.')[0]}.fasta", options.megan_na2t, f"{options.output}/megan/blastn/")) for file in sorted(os.listdir(f"{options.output}/blastn/"))]

		for p in results:
			p.get()
	else:
		print("Megan for blastn skipped.")

### Find Unmapped Reads

	if(options.unmapped_reads == "Y"):
		print("Beginning identifying reads which do not map to a contig.")
		for file in os.listdir(f"{options.output}/contigs/200/"):
			trimmed_input = f"{options.output}/bbduk/{file.split('.')[0].replace('sorted_200_', '')}_R1.fastq.gz"
			unmapped_output = f"{options.output}/unmapped_reads/{file.split('.')[0].replace('sorted_200_', '')}.fastq.gz"
			ref = f"{options.output}/contigs/200/{file}"
			subprocess.call(f"bbmap.sh in1={trimmed_input} in2={trimmed_input.replace('_R1', '_R2')} outu={unmapped_output} ref={ref} nodisk", shell = True)

	else:
		print(f"Unmapped reads generated.")
		
### Run Blastx and Megan

	# Blastx does not seem to benefit much from multiple processes, so this is run sequentially with megan.
	# For every file in contigs/1000, call blastx and then megan
	# Output to blastx directory

	if(options.blastx == "Y"):
		print("Beginning Blastx.")

		# Assign contig path
		contig_path = ""
		if(options.cluster == "Y"):
			contig_path = f"{options.output}/mmseq2/"
		else:
			contig_path = f"{options.output}/contigs/1000/"

		# Run blastx
		for file in os.listdir(contig_path):
			if(file.endswith(".fasta")):	
				print(f"Beginning Blastx for sample: {file.split('.')[0]}")

				runBlast("blastx", "blastx", options.nr_db, f"{contig_path}{file}", options.blastx_threads, options.blastn_descriptions, options.blastn_alignments, f"{options.output}/blastx/{file.split('.')[0]}.blastx")

				if(options.megan_blastx == "Y"):
					print(f"Beginning Megan for sample: {file.split('.')[0]}")
					pool = multiprocessing.Pool(processes = 1)
					runMegan(f"{options.output}/blastx/{file.split('.')[0]}.blastx", "BlastX", f"{contig_path}{file}", options.megan_pa2t, f"{options.output}/megan/blastx/")
				else:
					print("Megan for blastx skipped.")
	else:
		print("Blastx skipped.")

### Program Details

	print("Printing Angua3 version information")
	with open(f"{options.output}/Angua3_env.txt", "w") as angua_out:
		angua_out.write(f"Angua3 Version: {angua3_version}\nMegan nucleotide abin file: {options.megan_na2t}\nMegan protein abin file: {options.megan_pa2t}\nblastn nt database:  {options.nt_db}\nblastx nr database:  {options.nr_db}\nTrimmer:  {options.trimmer}\nCluster:  {options.cluster}\n")

	subprocess.call(f"conda list >> {options.output}/Angua3_env.txt", shell = True)

	print("Angua pipeline completed!")

######################################## Blast Function ########################################
# Can be used either for blastn or blastx.

def runBlast(blast_type, blast_task, database_path, query, threads, descriptions, alignments, output):
	blast = f"{blast_type} -task {blast_task} -db {database_path} -query {query} -num_threads {threads} -num_descriptions {descriptions} -num_alignments {alignments} -out {output}"

	blastchild = subprocess.Popen(str(blast), stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines = True, shell = (sys.platform != "win32"))
	blastOutput, blastError = blastchild.communicate()
	print(f"Blast Complete for query: {query}")

################################################################################################

######################################## Megan Function ######################################## 

def runMegan(megan_input, mode, reads, a2t, output):
	megan = f"blast2rma -i {megan_input} -f BlastText -bm {mode} -r {reads} -ms 75 -sup 1 -a2t {a2t} -o {output}"

	meganchild = subprocess.Popen(str(megan), stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines = True, shell = (sys.platform != "win32"))
	meganOutput, meganError = meganchild.communicate()
	print(f"Megan complete for sample: {megan_input.split('.')[0]}")

################################################################################################

#################################### Get Arguments Function ####################################
def parseArguments():
	parser = argparse.ArgumentParser(description = "Runs the Angua3 pipeline.")

	# Main arguments
	parser.add_argument("--input", help = "This is the location of the raw data directory.", required = True)
	parser.add_argument("--output", help = "This where the output data will be generated.", required = True)
	
	parser.add_argument("--nr_db", help = "This is the path to the nr database.", default = "/biostore/bigbio_00/smcgreig/Blast_databases/nr/nr_db_12072021/nr")
	parser.add_argument("--nt_db", help = "This is the path to the nt database.", default = "/biostore/bigbio_00/smcgreig/Blast_databases/nt/nt_db_29062021/nt")
	parser.add_argument("--megan_na2t", help = "This is the path to the megan nucl_acc2tax file.", default = "/home/smcgreig/miniconda2/envs/Angua3/opt/megan-6.12.3/megan-nucl-map-Jul2020.db")
	parser.add_argument("--megan_pa2t", help = "This is the path to the megan prot_acc2tax file.", default = "/home/smcgreig/miniconda2/envs/Angua3/opt/megan-6.12.3/megan-map-Jul2020-2.db")

	# Extra arguments, useful for if a specific job has failed and you don't want to start from scratch

	parser.add_argument("--create_dirs", help = "Creates the directory structure including folders for sickle/bbduk, trinity, blastx, blastn, contigs and megan, including sub directories. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--trimmer", help = "Run trimming. Default bbduk.", choices = ["sickle", "bbduk", "N"], default = "bbduk")
	parser.add_argument("--trinity", help = "Run trinity. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--sort", help = "Sort contigs, based on length, into >=200 and >=1000. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--cluster", help = "Clusters contigs >= 1000. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--blastn", help = "Run blastn. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--blastx", help = "Run blastx. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--megan_blastn", help = "Run megan for blastn. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--megan_blastx", help = "Run megan for blastx. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--unmapped_reads", help = "Map trimmed reads back to contigs to identify unmapped reads. Default N.", choices = ["Y", "N"], default = "N")

	# Paired end data
	parser.add_argument("--single_end", help = "Activate paired end mode. Expects file format *_R1_001*. Default expected format for paired end is *_L001_R1_001*,*_L001_R2_001*. Default N.", choices = ["Y", "N"], default = "N")

	# Tool specific parameters
	# Sickle
	parser.add_argument("--sickle_q", help = "Sickle phred quality trim parameter. Default 10", default = 10)
	parser.add_argument("--sickle_minl", help = "Sickle minimum length. Default 75", default = 75)

	# Bbduk
	parser.add_argument("--bbduk_adapters", help = "Bbduk adapter references.", default = "/home/smcgreig/Scripts/Angua3/RNA_adapters.fasta")
	parser.add_argument("--bbduk_q", help = "Bbduk phred quality trim parameter. Default 10", default = 10)
	parser.add_argument("--bbduk_minl", help = "Bbduk minimum length. Default 75", default = 75)

	# Trinity
	parser.add_argument("--trinity_cpu", help = "Trinity CPU parameter. Default 60.", default = 60)
	parser.add_argument("--trinity_mem", help = "Trinity max memory parameter. Default 200G", default = "200G")

	# MMseq2
	parser.add_argument("--cluster_perc", help = "What percentage identity to cluster at. Default 0.95.", default = 0.95)
	parser.add_argument("--cluster_threads", help = "Number of threads to run mmseq2 with. Default 60.", default = 60)

	# Blastn
	parser.add_argument("--blastn_pool", help = "This is the maximum number of blastn processes allowed in the pool at any one time. Default 8.", default = 8)
	parser.add_argument("--blastn_threads", help = "This is the number of threads used for each blastn process. Default 16.", default = 16)
	parser.add_argument("--blastn_descriptions", help = "This is the number of descriptions shown. Default 25.", default = 25)
	parser.add_argument("--blastn_alignments", help = "This is the number of alignments shown. Default 25.", default = 25)

	# Blastx
	parser.add_argument("--blastx_threads", help = "This is the number of threads used for running blastx. Default 130.", default = 130)
	parser.add_argument("--blastx_descriptions", help = "This is the number of descriptions shown. Default 25.", default = 25)
	parser.add_argument("--blastx_alignments", help = "This is the number of alignments shown. Default 25.", default = 25)

	# Megan
	parser.add_argument("--megan_processes", help = "This is the maximum number of megan processes allowed in the pool at any one time. Default 2.", default = 2)

	return parser.parse_args()

################################################################################################

############################################################################################################################## Functions End ###############################################################################################################################

if __name__ == '__main__':
	main()
