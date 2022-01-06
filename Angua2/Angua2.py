import os
import sys
import subprocess
import multiprocessing
import argparse
import re
import random
import in_place

def main():

	# Angua2 script version
	angua2_version = 1.5

	### Input arguments
	options = parseArguments()

	### Work out project folder by selecting the directory immediately before the raw data directory. If there is a trailing slash then remove this character.
	if (options.project_data_dir[-1] == "/"):
		project_path = options.project_data_dir[0:-1]
	project_path = project_path.split("/")
	project_directory = "/".join(project_path[0:-1]) + "/"
	raw_data_directory = project_path[-1] + "/"

### Intial directory setup
	if(options.create_dirs == "Y"):
		if(options.trimming == "sickle"):
			subprocess.call("mkdir -p " + project_directory + "sickle/unpaired", shell=True)
		if(options.trimming == "bbduk"):
			subprocess.call("mkdir -p " + project_directory + "bbduk", shell=True)
		subprocess.call("mkdir " + project_directory + "trinity", shell=True)
		subprocess.call("mkdir " + project_directory + "blastx", shell=True)
		subprocess.call("mkdir " + project_directory + "blastn", shell=True)
		subprocess.call("mkdir -p " + project_directory + "contigs/200", shell=True)
		subprocess.call("mkdir -p " + project_directory + "contigs/1000", shell=True)
		subprocess.call("mkdir -p " + project_directory + "megan/blastn", shell=True)
		subprocess.call("mkdir -p " + project_directory + "megan/blastx", shell=True)

		if(options.viraminer == "Y"):
			subprocess.call("mkdir -p " + project_directory + "viraminer/input", shell=True)
			subprocess.call("mkdir -p " + project_directory + "viraminer/priority_contigs", shell=True)
			subprocess.call("mkdir -p " + project_directory + "viraminer/priority_blastx", shell=True)
			subprocess.call("mkdir -p " + project_directory + "viraminer/priority_megan", shell=True)
	else:
		print("Directory structure creation skipped.")

### Run bbduk
	
	if(options.trimming == "bbduk"):
		print("Beginning bbduk")
		if(options.single_end == "Y"):
			for file in os.listdir(project_directory + raw_data_directory):
				samplename = file.split(".")
				if os.path.isfile(project_directory + raw_data_directory + file):
					output = samplename[0]
					subprocess.call(f"bbduk.sh in={project_directory}{raw_data_directory}{file} out={project_directory}bbduk/{output}.fastq minlen=100 ktrim=r k=23 mink=11 hdist=1 tbo ref={options.bbduk_adapters} qtrim=r trimq={options.bbduk_q}", shell=True)
					print("Bbduk: " + file + " completed trimming.")
		else:

			for fileR1 in os.listdir(project_directory + raw_data_directory):
				samplenameR1 = fileR1.split(".")
				if("R1" in fileR1):
					fileR2 = fileR1.replace("R1", "R2")
					if os.path.isfile(project_directory + raw_data_directory + fileR2):
						samplenameR2 = fileR2.split(".")
						output1 = samplenameR1[0]
						output2 = samplenameR2[0]
						subprocess.call(f"bbduk.sh in1={project_directory}{raw_data_directory}{fileR1} in2={project_directory}{raw_data_directory}{fileR2} out1={project_directory}bbduk/{output1}.fastq out2={project_directory}bbduk/{output2}.fastq minlen=100 ktrim=r k=23 mink=11 hdist=1 tbo ref={options.bbduk_adapters} qtrim=r trimq={options.bbduk_q}", shell=True)
						print("Bbduk: " + output1 + " completed trimming.")

		### Simplify Sample Names
		if(options.single_end == "Y"):
			subprocess.call("rename 's/_R1_001/_R1/' " + project_directory + "bbduk/*", shell=True)
		else:
			subprocess.call("rename 's/_L001_R1_001/_R1/' " + project_directory + "bbduk/*", shell=True)
			subprocess.call("rename 's/_L001_R2_001/_R2/' " + project_directory + "bbduk/*", shell=True)
		print("Files renamed.")
	else:
		print("Bbduk skipped.")

### Run Sickle

	# For each file in the given directory, get the name of the sample before the .fastq
	# If R1 is in the filename, then create the R2 filename variable
	# If the R2 file also exists, then create the output names for each file
	# Call sickle command

	if(options.trimming == "sickle"):
		print("Beginning Sickle")
		if(options.single_end == "Y"):
			for file in os.listdir(project_directory + raw_data_directory):
				samplename = file.split(".")
				if os.path.isfile(project_directory + raw_data_directory + file):
					output = samplename[0]
					subprocess.call("sickle se -f " + project_directory + raw_data_directory + file + " -t sanger -q " + str(options.sickle_q) + " -l 20 -o " + project_directory + "sickle/" + output + ".fastq -x", shell=True)
					print("Sickle: " + file + " completed trimming.")
		else:
			for fileR1 in os.listdir(project_directory + raw_data_directory):
				samplenameR1 = fileR1.split(".")
				if("R1" in fileR1):
					fileR2 = fileR1.replace("R1", "R2")
					if os.path.isfile(project_directory + raw_data_directory + fileR2):
						samplenameR2 = fileR2.split(".")
						output1 = samplenameR1[0]
						output2 = samplenameR2[0]
						subprocess.call("sickle pe -f " + project_directory + raw_data_directory + fileR1 + " -r" + project_directory + raw_data_directory + fileR2 + " -t sanger -q " + str(options.sickle_q) + " -l 100 -o " + project_directory + "sickle/" + output1 + ".fastq -p " + project_directory + "sickle/" + output2 + ".fastq -s " + project_directory + "sickle/" +  "unpaired/" + output1 + ".fastq -x", shell=True)
						print("Sickle: " + fileR1 + " completed trimming.")

		### Simplify Sample Names
		if(options.single_end == "Y"):
			subprocess.call("rename 's/_R1_001/_R1/' " + project_directory + "sickle/*", shell=True)
		else:
			subprocess.call("rename 's/_L001_R1_001/_R1/' " + project_directory + "sickle/*", shell=True)
			subprocess.call("rename 's/_L001_R2_001/_R2/' " + project_directory + "sickle/*", shell=True)
		print("Files renamed.")
	else:
		print("Sickle skipped.")

### Run Trinity

	# For each file in the given directory, get the name of the sample before the .fastq
	# If R1 is in the filename, then create the R2 filename variable
	# If the R2 file also exists, then create the output names for each file
	# Call Trinity command

	if(options.trinity == "Y"):
		trimmer = raw_data_directory

		if(options.trimming == "sickle"):
			trimmer = "sickle/"
		elif(options.trimming == "bbduk"):
			trimmer = "bbduk/"
		else:
			trimmer = raw_data_directory
		
		print("Beginning Trinity")
		if(options.single_end == "Y"):
			for file in os.listdir(project_directory + trimmer):
				samplename = file.split(".")
				if os.path.isfile(project_directory + trimmer + file):
					output = samplename[0]
					subprocess.call("Trinity --seqType fq --max_memory " + str(options.trinity_mem) + " --single " + project_directory + trimmer + file + " --CPU " + str(options.trinity_cpu) + " --full_cleanup --output " + project_directory + "trinity/" + output + "_trinity > " + project_directory + "trinity/" + output + ".log", shell=True)
					print("Trinity: " + file + " completed assembly.")
		else:
			for fileR1 in os.listdir(project_directory + trimmer):
				samplenameR1 = fileR1.split(".")
				if("R1" in fileR1):
					fileR2 = fileR1.replace("R1", "R2")
					if os.path.isfile(project_directory + trimmer + fileR2):
						samplenameR2 = fileR2.split(".")
						output1 = samplenameR1[0]
						output2 = samplenameR2[0]
						subprocess.call("Trinity --seqType fq --max_memory " + str(options.trinity_mem) + " --left " + project_directory + trimmer + fileR1 + " --right " + project_directory + trimmer + fileR2 + " --CPU " + str(options.trinity_cpu) + " --full_cleanup --output " + project_directory + "trinity/" + output1 + "_trinity > " + project_directory + "trinity/" + output1 + ".log", shell=True)
						print("Trinity: " + fileR1 + " completed assembly.")

		### Simplify Sample Names
		subprocess.call("rename 's/_R1_trinity.Trinity//' " + project_directory + "trinity/*", shell=True)
		print("Files renamed.")

		### Adjust Trinity names
		for file in os.listdir(project_directory + "trinity/"):
			if(file.endswith(".fasta")):
				file_split = file.split("_")
				sample_number = file_split[-1][:-6]

				with in_place.InPlace(project_directory + "trinity/" + file) as trinity_file:
					for line in trinity_file:
						if(line.startswith(">")):
							line = line.replace(">", ">" + project_directory[:-1] + "_" + sample_number + "_")
							trinity_file.write(line)
						else:
							trinity_file.write(line)
		print("Run and sample ID prepended to fasta file IDs.")

	else:
		print("Trinity skipped.")

### Run contig filter

	# For each _trinity file in the trinity directory, call filter script
	# Repeat for 200 and 1000 bases, and place them in the desired output directories

	if(options.filter == "Y"):
		print("Beginning contig size filter")
		for file in os.listdir(project_directory + "trinity/"):
			if(file.endswith(".fasta")):
				subprocess.call("python Scripts/contigSizeFilter2.py " + project_directory + "trinity/" + file + " " + str(200) + " " + project_directory + "contigs/200/", shell=True)
				subprocess.call("python Scripts/contigSizeFilter2.py " + project_directory + "trinity/" + file + " " + str(1000) + " " + project_directory + "contigs/1000/", shell=True)
				print("Contig filter: " + file + " completed sorting.")
	else:
		print("Contig filter skipped.")

### Run Blastn

	# Blastn seems to work better with multiple processes instead of multiple threads.
	# As such, a pool was create that can hold 8 processes at once. Once a process has been completed, a new process will be added until everything quened has been completed.  
	# args correspond to the runBlast function

	if(options.blastn == "Y"):
		print("Beginning Blastn")
		pool = multiprocessing.Pool(processes = options.blastn_pool)
		results = [pool.apply_async(runBlast, args = ("blastn", options.nt_db, project_directory, "contigs/200/", file, "blastn/", options.blastn_threads)) for file in sorted(os.listdir(project_directory + "contigs/200/"))]
		for p in results:
			p.get()
	else:
		print("Blastn skipped.")

### Run Megan

	# A pool was created that can hold 2 processes at once. Once a process has been completed, a new process will be added until everything quened has been completed. 
	# args correspond to the runMegan function

	if(options.megan == "Y"):
		print("Beginning Megan for Blastn files")
		pool = multiprocessing.Pool(processes = 2)
		results = [pool.apply_async(runMegan, args = (project_directory, "blastn/", file, "contigs/200/", "megan/blastn/", options.megan_na2t)) for file in sorted(os.listdir(project_directory + "blastn/"))]
		for p in results:
			p.get()
	else:
		print("megan blastn skipped.")

### Run ViraMiner

	# This section of code runs the viraminer priority blast experimental setup.
	# A rough outline is as follows:
	# Select all contigs in the contigs 200 folder. Split these sequences into length of 300 and then discard anything shorter than that length. This length is required for the default viraminer model.
	# The sequences of length 300 are then formatted ready for the viraminer script. This includes removing any underscores and foramtting the line to be ID, Seq, Group.
	# These files are ran with viraminer, which uses a different environment as it is a python 2.7 script (this is python 3). The predictions are combined with the input files and returned to the priority contigs directory.
	# The returned file is then filtered to identify contigs which viraminer predicts, with 80% certainty, are viral. The IDs of these contigs are stored in a file and the contigs that they represent are extracted.
	# These contigs are then subject to priority blastx + megan, which happens before the main blastx step in the pipeline.
	
	if(options.viraminer == "Y"):
		print("Beginning ViraMiner for contigs 200 files")
		for file in os.listdir(project_directory + "contigs/200/"):
			if(os.stat(project_directory + "contigs/200/" + file).st_size != 0):
				subprocess.call("splitter -sequence " + project_directory + "contigs/200/" + file + " -size 300 -outseq " + project_directory + "viraminer/input/split_" + file, shell = True) 
				subprocess.call("reformat.sh in=" + project_directory + "viraminer/input/split_" + file + " out=" + project_directory + "viraminer/input/split_300_" + file + " minlength=300 fastawrap=300", shell = True)

				# ### Simplify Sample Names
				# subprocess.call("rename 's/_R1_trinity.Trinity//' " + project_directory + "trinity/*", shell=True)
				# print("Files renamed.")

				viraminer_list = []
				with open(project_directory + "viraminer/input/split_300_" + file, "r") as trinity:
					fasta = ""
					for line in trinity:
						line = line.rstrip()
						if("len" in line):
							line = re.sub(" len.*", "", line)
							line = re.sub("_", "-", line)
							fasta = line
						else:
							fasta = fasta + "," + line + "," + str(random.randint(0, 1))
							viraminer_list.append(fasta)

				with open(project_directory + "viraminer/input/viraminer" + file, "w") as vir_out:
					for seq in viraminer_list:
						vir_out.write("%s\n" % seq)

		for file in os.listdir(project_directory + "viraminer/input/"):
			if("viraminersorted" in file):
				subprocess.call("/home/smcgreig/miniconda2/envs/MachineLearning/bin/python Models/ViraMiner_code/predict_only.py --input_file " + project_directory + "viraminer/input/" + file + " --model_path Models/ViraMiner/ViraMiner_model_afterFT.hdf5", shell = True)
				with open("Models/ViraMiner/ViraMiner_model_afterFT_TEST_predictions.txt") as preds:
					with open(project_directory + "viraminer/input/" + file) as seqs:
						with open(project_directory + "viraminer/priority_contigs/predictions_" + file + ".txt", "w") as preds_out:
							predline = preds.readlines()
							seqline = seqs.readlines()
							for line1, line2 in zip(predline, seqline):
								preds_out.write("{},{}\n".format(line2.rstrip(), line1.rstrip()))

		for file in os.listdir(project_directory + "viraminer/priority_contigs/"):
			with open(project_directory + "viraminer/priority_contigs/" + file) as preds:
				output = {}
				for line in preds:
					line = line.rstrip()
					line = line.split(",")
					if(float(line[3]) >= 0.8):
						contig_name = line[0].split("-")
						contig_name = "_".join(contig_name[0:4])
						output[contig_name[1:]] = line[3]

				with open(project_directory + "viraminer/priority_contigs/virus_" + file, "w") as virus_contigs:
					for seq in output:
						print(seq, file = virus_contigs)

		for vpreds in os.listdir(project_directory + "viraminer/priority_contigs/"):
			vpreds_ids = project_directory + "viraminer/priority_contigs/" + vpreds
			if("virus_predictions" in vpreds):
				for contigs in os.listdir(project_directory + "contigs/200/"):
					contigfile = project_directory + "contigs/200/" + contigs
					if(contigs in vpreds):
						subprocess.call("filterbyname.sh in=" + contigfile + " out=" + project_directory + "viraminer/priority_contigs/contigs_" + vpreds + ".fasta include=t substring=t prefix=t names=" + vpreds_ids, shell = True)

		for file in os.listdir(project_directory + "viraminer/priority_contigs/"):
			if("contigs_virus_predictions" in file):
				subprocess.call("blastx -db " + options.nr_db + " -query " + project_directory + "viraminer/priority_contigs/" + file + " -num_threads " + str(options.blastx_threads) + " -num_alignments 25 -num_descriptions 25 -out " + project_directory + "viraminer/priority_blastx/" + file + ".blastx", shell=True)

		pool = multiprocessing.Pool(processes = 2)
		results = [pool.apply_async(runMegan, args = (project_directory, "viraminer/priority_blastx/", file, "viraminer/priority_contigs/", "viraminer/priority_megan/", options.megan_pa2t)) for file in sorted(os.listdir(project_directory + "viraminer/priority_blastx/"))]
		for p in results:
			p.get()
		
### Run Blastx and Megan

	# Blastx does not seem to benefit much from multiple processes, so this is run sequentially with megan.
	# For every file in contigs/1000, call blastx and then megan
	# Output to blastx directory

	if(options.blastx == "Y"):
		print("Beginning Blastx and Megan")
		for file in os.listdir(project_directory + "contigs/1000/"):
			print("Beginning Blastx for file: " + file)
			if(os.stat(project_directory + "contigs/1000/" + file).st_size != 0):
				subprocess.call("blastx -db " + options.nr_db + " -query " + project_directory + "contigs/1000/" + file + " -num_threads " + str(options.blastx_threads) + " -num_alignments 25 -num_descriptions 25 -out " + project_directory + "blastx/" + file + ".blastx", shell=True)
				print(file + " blastx completed")

			if(options.megan == "Y"):
				print("Beginning Megan for Blastx file: " + file)
				pool = multiprocessing.Pool(processes = 1)
				runMegan(project_directory, "blastx/", file + ".blastx", "contigs/1000/", "megan/blastx/", options.megan_pa2t)
			else:
				print("megan blastx skipped.")
	else:
		print("blastx skipped.")

### Program Details

	print("Printing Angua2 version information")
	subprocess.call("echo \"Angua2 Version: " + str(angua2_version) + "\" > " + project_directory + "Angua2_env.txt", shell = True)
	subprocess.call("echo \"Megan abin file: " + options.megan_na2t + "\" >> " + project_directory + "Angua2_env.txt", shell = True)
	subprocess.call("echo \"Megan abin file: " + options.megan_pa2t + "\" >> " + project_directory + "Angua2_env.txt", shell = True)
	subprocess.call("echo \"blastn file: " + options.nt_db + "\" >> " + project_directory + "Angua2_env.txt", shell = True)
	subprocess.call("echo \"blastx file: " + options.nr_db + "\" >> " + project_directory + "Angua2_env.txt", shell = True)
	subprocess.call("conda list >> " + project_directory + "Angua2_env.txt", shell = True)

	print("Angua pipeline completed!")


################################################################################################################################ Functions ################################################################################################################################

################################ Blast Multiprocessing Function ################################
# Takes database path, project directory, contig directory, blastn/x file, blastn/x output directory and the number of threads as arguments.
# Can be used either for blastn or blastx.

def runBlast(blast_type, databasepath, project_directory, contig_dir, blast_file, blast_out_dir, threads):
	blast = blast_type + " -db " + databasepath + " -query " + project_directory + contig_dir + blast_file + " -num_threads " + str(threads) + " -num_alignments 25 -num_descriptions 25 -out " + project_directory + "blastn/" + blast_file + ".blastn"
	print(blast)
	blastchild = subprocess.Popen(str(blast),
									stdout=subprocess.PIPE,
									stderr=subprocess.PIPE,
									universal_newlines=True,
									shell=(sys.platform != "win32"))
	blastOutput, blastError = blastchild.communicate()
	print("Blast Complete for file: " + blast_file)

################################################################################################

################################ Megan Multiprocessing Function ################################
# Takes project directory, blast directory, blast file, contig directory and megan output directory as inputs.
# Can be used for either blastn or blastx.

def runMegan(project_directory, blast_directory, blast_file, contig_directory, megan_out_dir, a2t_file):

	samplename = ""
	if(".blastn" in blast_file):
		samplenamesplit = blast_file.split(".blastn")
		samplename = samplenamesplit[0]
		megan = "blast2rma -i " + project_directory + blast_directory + blast_file + " -f BlastText -bm BlastN -r " + project_directory  + contig_directory + samplename + " -ms 75 -sup 1 -a2t " + a2t_file + " -o " + project_directory + megan_out_dir
		meganchild = subprocess.Popen(str(megan),
										stdout=subprocess.PIPE,
										stderr=subprocess.PIPE,
										universal_newlines=True,
										shell=(sys.platform != "win32"))
		meganOutput, meganError = meganchild.communicate()
		print("Megan Complete for file: " + samplename)


	if(".blastx" in blast_file):
		samplenamesplit = blast_file.split(".blastx")
		samplename = samplenamesplit[0]
		megan = "blast2rma -i " + project_directory + blast_directory + blast_file + " -f BlastText -bm BlastX -r " + project_directory  + contig_directory + samplename + " -ms 75 -sup 1 -a2t " + a2t_file + " -o " + project_directory + megan_out_dir
		meganchild = subprocess.Popen(str(megan),
										stdout=subprocess.PIPE,
										stderr=subprocess.PIPE,
										universal_newlines=True,
										shell=(sys.platform != "win32"))
		meganOutput, meganError = meganchild.communicate()
		print("Megan Complete for file: " + samplename)


################################################################################################

#################################### Get Arguments Function ####################################
def parseArguments():
	parser = argparse.ArgumentParser(description = "Runs the Angua2 pipeline. This pipeline consists of sickle, trinity, blastn, blastx and megan.")

	# Main arguments
	parser.add_argument("--project_data_dir", help = "This is the location of the raw data directory, which should be in the directory directly beneath the project directory.", required = True)
	parser.add_argument("--nr_db", help = "This is the path to the nr database.", default = "/biostore/bigbio_00/smcgreig/Blast_databases/nr/nr_db_29052020/nr")
	parser.add_argument("--nt_db", help = "This is the path to the nt database.", default = "/biostore/bigbio_00/smcgreig/Blast_databases/nt/nt_db_11082020/nt")
	parser.add_argument("--megan_na2t", help = "This is the path to the megan nucl_acc2tax file.", default = "/home/smcgreig/miniconda2/envs/Angua2/opt/megan-6.12.3/megan-nucl-map-Jul2020.db")
	parser.add_argument("--megan_pa2t", help = "This is the path to the megan prot_acc2tax file.", default = "/home/smcgreig/miniconda2/envs/Angua2/opt/megan-6.12.3/megan-map-Jul2020-2.db")

	# Extra arguments, useful for if a specific job has failed and you don't want to start from scratch

	parser.add_argument("--create_dirs", help = "Creates the directory structure including folders for sickle/bbduk, trinity, blastx, blastn, contigs and megan, including sub directories. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--trimming", help = "Run trimming. Default sickle.", choices = ["sickle", "bbduk", "N"], default = "sickle")
	parser.add_argument("--trinity", help = "Run trinity. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--filter", help = "Run contig filter. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--blastn", help = "Run blastn. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--blastx", help = "Run blastx. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--megan", help = "Run megan. Default Y.", choices = ["Y", "N"], default = "Y")
	parser.add_argument("--viraminer", help = "Run viraminer. Default N.", choices = ["Y", "N"], default = "N")

	# Paired end data
	parser.add_argument("--single_end", help = "Activate paired end mode. Expects file format *_R1_001*. Default expected format for paired end is *_L001_R1_001*,*_L001_R2_001*. Default N.", choices = ["Y", "N"], default = "N")

	# Tool specific parameters
	# Sickle
	parser.add_argument("--sickle_q", help = "Sickle phred quality trim parameter. Default 20", default = 20)

	# Bbduk
	parser.add_argument("--bbduk_adapters", help = "Bbduk adapter references.", default = "Scripts/RNA_adapters.fasta")
	parser.add_argument("--bbduk_q", help = "Bbduk phred quality trim parameter. Default 10", default = 10)

	# Trinity
	parser.add_argument("--trinity_cpu", help = "Trinity CPU parameter. Default 60.", default = 60)
	parser.add_argument("--trinity_mem", help = "Trinity max memory parameter. Default 200G", default = "200G")

	# Blastn
	parser.add_argument("--blastn_pool", help = "This is the maximum number of blastn processes allowed in the pool at any one time. Default 8.", default = 8)
	parser.add_argument("--blastn_threads", help = "This is the number of threads used for each blastn process. Default 16.", default = 16)

	# Blastx
	parser.add_argument("--blastx_threads", help = "This is the number of threads used for running blastx. Default 130.", default = 130)

	# Megan
	parser.add_argument("--megan_processes", help = "This is the maximum number of megan processes allowed in the pool at any one time. Default 2.", default = 2)

	return parser.parse_args()

################################################################################################

############################################################################################################################## Functions End ###############################################################################################################################

if __name__ == '__main__':
	main()
