**Angua Virus Pipeline**

**Note: This repository is no longering being suppoted. Please visit the Angua3 repo for the latest supported version of Angua https://fred.fera.co.uk/smcgreig/angua3)**

Angua is a pipeline for the discovery of viruses from short read RNA data. The pipeline follows the steps outlined below.

1.  Directory setup
2.  Quality control with Sickle or bbduk
3.  Assembly with Trinity
4.  Contig filtering in Python
5.  Blastn on contigs >=200
6.  Blastx on contigs >=1000
7.  Produce Megan outputs for blastn and blastx results  
8.  (Experimental) Identifying viral contigs with a viraminer model, then conduct a priority blastx on these contigs 

**Setting up a Conda Environment**

The Angua2.yml file provided allows for easily setting up a conda environment, which will contain all of the tools needed to run the script.
Ensure that conda is installed and active, and then execute `conda env create -f Angua2.yml` to create and install all tools into the environment.
Please activate the conda environment before running the scripts.

**Installation**

Create a directory named `Angua`. Inside this, create a directory called `Scripts`. Copy the `Angua2.py` and `contigSizeFilter2.py` scripts into this directory.

Megan mapping files can be download from the MEGAN6 download page (http://ab.inf.uni-tuebingen.de/data/software/megan6/download/welcome.html). 
The files required will look something like:
-   megan-map-Jul2020-2.db.zip `# needed for --megan_pa2t` 
-   megan-nucl-Jul2020.db.zip  `# needed for --megan_na2t` 

Make sure to unzip these files and then pass the file path of these files to the respective  `--megan_pa2t` and `--megan_na2t` parameters.

Blast databases can be downloaded with the `update_blastdb.pl` command. The required databases can be download using the command:

`update_blastdb.pl --decompress --blastdb_version 5 nt`
`update_blastdb.pl --decompress --blastdb_version 5 nr`

It is often a good idea to save these databases in a folder with the date of download.

**Getting Started**

Create a directory named `Project_<Project number>` in the `Angua` directory. Inside this new directory, create a `raw_data` directory and place your raw data there.

Usage:

`python Scripts/Angua2.py --project_raw_data_dir Project_<Project number>/raw_data/ --nr_db <path to nrdb> --nt_db <path to ntdb> --megan_na2t <path to megan nucl_acc2tax file> --megan_pa2t <path to megan prot_acc2tax>`

Execute the above command from the `Angua` directory.

**Parameters**


```
usage: Angua2.py [-h] --project_data_dir PROJECT_DATA_DIR [--nr_db NR_DB]
                 [--nt_db NT_DB] [--megan_na2t MEGAN_NA2T]
                 [--megan_pa2t MEGAN_PA2T] [--create_dirs {Y,N}]
                 [--trimming {sickle,bbduk,N}] [--trinity {Y,N}]
                 [--filter {Y,N}] [--blastn {Y,N}] [--blastx {Y,N}]
                 [--megan {Y,N}] [--viraminer {Y,N}] [--single_end {Y,N}]
                 [--sickle_q SICKLE_Q] [--bbduk_adapters BBDUK_ADAPTERS]
                 [--bbduk_q BBDUK_Q] [--trinity_cpu TRINITY_CPU]
                 [--trinity_mem TRINITY_MEM] [--blastn_pool BLASTN_POOL]
                 [--blastn_threads BLASTN_THREADS]
                 [--blastx_threads BLASTX_THREADS]
                 [--megan_processes MEGAN_PROCESSES]

Runs the Angua2 pipeline. This pipeline consists of sickle, trinity, blastn,
blastx and megan.

optional arguments:
  -h, --help            show this help message and exit
  --project_data_dir PROJECT_DATA_DIR
                        This is the location of the raw data directory, which
                        should be in the directory directly beneath the
                        project directory.
  --nr_db NR_DB         This is the path to the nr database.
  --nt_db NT_DB         This is the path to the nt database.
  --megan_na2t MEGAN_NA2T
                        This is the path to the megan nucl_acc2tax file.
  --megan_pa2t MEGAN_PA2T
                        This is the path to the megan prot_acc2tax file.
  --create_dirs {Y,N}   Creates the directory structure including folders for
                        sickle/bbduk, trinity, blastx, blastn, contigs and
                        megan, including sub directories. Default Y.
  --trimming {sickle,bbduk,N}
                        Run trimming. Default sickle.
  --trinity {Y,N}       Run trinity. Default Y.
  --filter {Y,N}        Run contig filter. Default Y.
  --blastn {Y,N}        Run blastn. Default Y.
  --blastx {Y,N}        Run blastx. Default Y.
  --megan {Y,N}         Run megan. Default Y.
  --viraminer {Y,N}     Run viraminer. Default N.
  --single_end {Y,N}    Activate paired end mode. Expects file format
                        *_R1_001*. Default expected format for paired end is
                        *_L001_R1_001*,*_L001_R2_001*. Default N.
  --sickle_q SICKLE_Q   Sickle phred quality trim parameter. Default 20
  --bbduk_adapters BBDUK_ADAPTERS
                        Bbduk adapter references.
  --bbduk_q BBDUK_Q     Bbduk phred quality trim parameter. Default 10
  --trinity_cpu TRINITY_CPU
                        Trinity CPU parameter. Default 60.
  --trinity_mem TRINITY_MEM
                        Trinity max memory parameter. Default 200G
  --blastn_pool BLASTN_POOL
                        This is the maximum number of blastn processes allowed
                        in the pool at any one time. Default 8.
  --blastn_threads BLASTN_THREADS
                        This is the number of threads used for each blastn
                        process. Default 16.
  --blastx_threads BLASTX_THREADS
                        This is the number of threads used for running blastx.
                        Default 130.
  --megan_processes MEGAN_PROCESSES
                        This is the maximum number of megan processes allowed
                        in the pool at any one time. Default 2.

```
