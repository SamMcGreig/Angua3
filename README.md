# Angua3

**Angua3 Virus Pipeline**

Angua3 is a pipeline for the discovery of viruses from Illumina RNA sequencing. The pipeline follows the steps outlined below. The pipeline is based off of the older Angua2 pipeline. It has been updated to run with Python3 and with a tidier code base.

1.  Directory setup
2.  Quality control with Sickle or bbduk
3.  Assembly with Trinity
4.  Contig filtering
5.  Blastn on contigs >=200
6.  Blastx on contigs >=1000. Optionally cluster contigs beforehand to reduce blastx time.
7.  Produce Megan outputs for blastn and blastx results

**Setting up a Conda Environment**

The Angua3.yml file provided allows for easily setting up a conda environment, which will contain all of the tools needed to run the script.
Ensure that conda is installed and active, and then execute `conda env create -f Angua3.yml` to create and install all tools into the environment.
Please activate the conda environment before running the scripts.

**Installation**

Create a directory named `Angua`. Inside this, create a directory called `Scripts`. Copy the `Angua3.py` script into this directory. You can also place the `RNA_adapters.fasta` file in this location, if you like.

Megan mapping files can be download from the MEGAN6 download page (https://software-ab.informatik.uni-tuebingen.de/download/megan6/welcome.html). 
The files required will look something like:
-   megan-map-Jul2020-2.db.zip `# needed for --megan_pa2t` 
-   megan-nucl-Jul2020.db.zip  `# needed for --megan_na2t` 

Make sure to unzip these files and then pass the file path of these files to the respective  `--megan_pa2t` and `--megan_na2t` parameters.

Blast databases can be downloaded with the `update_blastdb.pl` command, which will be downloaded along with the blast software upon creation of the Angua3 conda environment. The required databases can be downloaded using the commands:

`update_blastdb.pl --decompress --blastdb_version 5 nt`

`update_blastdb.pl --decompress --blastdb_version 5 nr`

It is often a good idea to save these databases in a folder with the date of download.

**Getting Started**

Create a directory named i.e `Project` in the `Angua` directory. Inside this new directory, create a `raw_data` directory and place your raw data there.

Usage:

`python Scripts/Angua3.py --input Project/raw_data/ --output Project --nr_db <path to nrdb> --nt_db <path to ntdb> --megan_na2t <path to megan nucl_acc2tax file> --megan_pa2t <path to megan prot_acc2tax>`

Execute the above command from the `Angua` directory.

**Future Features**
-   Machine learning features (based off of tools such as ViraMiner and DeepVirFinder) to detect contigs that are most likely to be of viral origin, and inspecting these first with a more rigorous blastx search.
-   Back mapping trimmed reads to filtered contigs, and then checking for viruses in any reads which did not map.
-   Priorty sample queueing i.e the ability to specify high priority samples to be analysed before lower priority samples.
-   Contig polishing with Pilon.

**Citation**

If you are using Angua, please consider citing:

Fowkes, A.R.; McGreig, S.; Pufal, H.; Duffy, S.; Howard, B.; Adams, I.P.; Macarthur, R.; Weekes, R.; Fox, A. Integrating High throughput Sequencing into Survey Design Reveals Turnip Yellows Virus and Soybean Dwarf Virus in Pea (Pisum Sativum) in the United Kingdom. Viruses 2021, 13, 2530. https://doi.org/10.3390/v13122530 

https://www.mdpi.com/1999-4915/13/12/2530

**Parameters**


```
usage: Angua3.py [-h] --input INPUT --output OUTPUT [--nr_db NR_DB]
                 [--nt_db NT_DB] [--megan_na2t MEGAN_NA2T]
                 [--megan_pa2t MEGAN_PA2T] [--create_dirs {Y,N}]
                 [--trimmer {sickle,bbduk,N}] [--trinity {Y,N}] [--sort {Y,N}]
                 [--cluster {Y,N}] [--blastn {Y,N}] [--blastx {Y,N}]
                 [--megan_blastn {Y,N}] [--megan_blastx {Y,N}]
                 [--unmapped_reads {Y,N}] [--single_end {Y,N}]
                 [--sickle_q SICKLE_Q] [--sickle_minl SICKLE_MINL]
                 [--bbduk_adapters BBDUK_ADAPTERS] [--bbduk_q BBDUK_Q]
                 [--bbduk_minl BBDUK_MINL] [--trinity_cpu TRINITY_CPU]
                 [--trinity_mem TRINITY_MEM] [--cluster_perc CLUSTER_PERC]
                 [--cluster_threads CLUSTER_THREADS]
                 [--blastn_pool BLASTN_POOL] [--blastn_threads BLASTN_THREADS]
                 [--blastn_descriptions BLASTN_DESCRIPTIONS]
                 [--blastn_alignments BLASTN_ALIGNMENTS]
                 [--blastx_threads BLASTX_THREADS]
                 [--blastx_descriptions BLASTX_DESCRIPTIONS]
                 [--blastx_alignments BLASTX_ALIGNMENTS]
                 [--megan_processes MEGAN_PROCESSES]

Runs the Angua3 pipeline.

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         This is the location of the raw data directory.
  --output OUTPUT       This where the output data will be generated.
  --nr_db NR_DB         This is the path to the nr database.
  --nt_db NT_DB         This is the path to the nt database.
  --megan_na2t MEGAN_NA2T
                        This is the path to the megan nucl_acc2tax file.
  --megan_pa2t MEGAN_PA2T
                        This is the path to the megan prot_acc2tax file.
  --create_dirs {Y,N}   Creates the directory structure including folders for
                        sickle/bbduk, trinity, blastx, blastn, contigs and
                        megan, including sub directories. Default Y.
  --trimmer {sickle,bbduk,N}
                        Run trimming. Default bbduk.
  --trinity {Y,N}       Run trinity. Default Y.
  --sort {Y,N}          Sort contigs, based on length, into >=200 and >=1000.
                        Default Y.
  --cluster {Y,N}       Clusters contigs >= 1000. Default Y.
  --blastn {Y,N}        Run blastn. Default Y.
  --blastx {Y,N}        Run blastx. Default Y.
  --megan_blastn {Y,N}  Run megan for blastn. Default Y.
  --megan_blastx {Y,N}  Run megan for blastx. Default Y.
  --unmapped_reads {Y,N}
                        Map trimmed reads back to contigs to identify unmapped
                        reads. Default N.
  --single_end {Y,N}    Activate paired end mode. Expects file format
                        *_R1_001*. Default expected format for paired end is
                        *_L001_R1_001*,*_L001_R2_001*. Default N.
  --sickle_q SICKLE_Q   Sickle phred quality trim parameter. Default 10
  --sickle_minl SICKLE_MINL
                        Sickle minimum length. Default 75
  --bbduk_adapters BBDUK_ADAPTERS
                        Bbduk adapter references.
  --bbduk_q BBDUK_Q     Bbduk phred quality trim parameter. Default 10
  --bbduk_minl BBDUK_MINL
                        Bbduk minimum length. Default 75
  --trinity_cpu TRINITY_CPU
                        Trinity CPU parameter. Default 60.
  --trinity_mem TRINITY_MEM
                        Trinity max memory parameter. Default 200G
  --cluster_perc CLUSTER_PERC
                        What percentage identity to cluster at. Default 0.95.
  --cluster_threads CLUSTER_THREADS
                        Number of threads to run mmseq2 with. Default 60.
  --blastn_pool BLASTN_POOL
                        This is the maximum number of blastn processes allowed
                        in the pool at any one time. Default 8.
  --blastn_threads BLASTN_THREADS
                        This is the number of threads used for each blastn
                        process. Default 16.
  --blastn_descriptions BLASTN_DESCRIPTIONS
                        This is the number of descriptions shown. Default 25.
  --blastn_alignments BLASTN_ALIGNMENTS
                        This is the number of alignments shown. Default 25.
  --blastx_threads BLASTX_THREADS
                        This is the number of threads used for running blastx.
                        Default 130.
  --blastx_descriptions BLASTX_DESCRIPTIONS
                        This is the number of descriptions shown. Default 25.
  --blastx_alignments BLASTX_ALIGNMENTS
                        This is the number of alignments shown. Default 25.
  --megan_processes MEGAN_PROCESSES
                        This is the maximum number of megan processes allowed
                        in the pool at any one time. Default 2.

```
