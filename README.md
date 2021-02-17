# stepRNA

stepRNA is a RNA-seq read processor, based on bowtie2, that will align query reads to reference sequences and output information about the type and length of overhangs uncovered. It was originally developed for processing DICER RNAs but is not limited to just this.

#### Installation

To install with pip:

```pip3 install stepRNA```

To install from source:
Make sure you have dist and build direcories for stepRNA (download from GitHub repository)
Then run:
```pip3 install -e .```

Alternative:
Download the GitHub repository then make the file executable with:
```chmod +x scripts/stepRNA```

Then either:
- Make a symbolic link for stepRNA.py to your bin directory
```ln FULL_PATH_TO_stepRNA.py FULL_PATH_TO_BIN```
- Copy the stepRNA into your bin directory
```cp stepRNA.py FULL_PATH_TO_BIN```
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
