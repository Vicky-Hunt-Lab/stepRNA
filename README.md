# stepRNA

stepRNA is a RNA-seq read processor, based on bowtie2, that will align small RNA-seq query reads (passenger sequences) to reference sequences (siRNAs) and output information about the type and length of overhangs uncovered. It was originally developed for uncovering Dicer processing signatures but is not limited to just this.

#### Installation

In order for stepRNA to run you must have:
- Bowtie2 >= vX.X.X (see BOWTIE2_URL for more details)
- Bio >= vX.X.X (automatically installed with pip)
- numpy >= vX.X.X (autoamtically installed with pip)
- pysam >= vX.X.X (autoamtically installed with pip)

**To install with pip:**

```pip3 install stepRNA```

**To install from source:**

1) Download stepRNA-X.X.X.tar.gz from PyPI_URL:
2) Unzip it:

```tar -xvzf stepRNA-X.X.X.tar.gz```

3) Move into stepRNA-X.X.X repository:

```cd stepRNA-X.X.X```

4) Run the installation:

```pip3 install .```


**Alternative:**

Download the GitHub repository then make the script executable with:
```chmod +x bin/stepRNA```

Then either:
- Make a symbolic link for stepRNA.py to your bin directory
```ln FULL_PATH_TO_stepRNA.py FULL_PATH_TO_BIN```
- Direct PATH to the directory with stepRNA in it
```export PATH=$PATH:FULL_PATH_TO_GitHubRepo```

#### Use:
See the documentation for a detailed description on how to use stepRNA.
The quickest way to use stepRNA:

```stepRNA --reference REFERENCE --reads QUERY```
 
This will align the reads to the reference sequences and output into the current diretory using the READS filename as the prefix. The following should be met:
- **REFERENCE** and **QUERY** must have unique FASTA headers (if not use ```-u```)
- If **REFERENCE** have been extracted from **QUERY** then it is advised to remove the exact matches from the QUERY (use ```-e```) - ***This takes a long time!***

Helpful options:
- ```--name``` can be used to customise the prefix name
- ```--directory``` can be used to specify an output directory
