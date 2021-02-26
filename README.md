# stepRNA

## Overview

stepRNA is a RNA-seq read processor, based on bowtie2, that will align small RNA-seq query reads (passenger sequences) to reference sequences (siRNAs) and output information about the type and length of overhangs uncovered. It was originally developed for uncovering Dicer processing signatures but is not limited to just this.

**Table of Contents:**
- Installation
- Use
- Example
- News
- Licence
- Additional Information

## Installation

In order for stepRNA to run you must have:
- Bowtie2 >= v2.3.4 (see [BOWTIE2 website](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for more details)
- Bio >= v0.3.0 (automatically installed with pip)
- numpy >= v1.20.1 (autoamtically installed with pip)
- pysam >= v0.16.0.1 (autoamtically installed with pip)

**To install with pip (recommended):**

```pip install stepRNA```

If this doesn't work due to non-root access issues see **Installing into a Virtual Environment**

**To install from source:**

1) Download stepRNA-X.X.X.tar.gz from PyPI_URL:
2) Unzip it:

```tar -xvzf stepRNA-X.X.X.tar.gz```

3) Move into stepRNA-X.X.X repository:

```cd stepRNA-X.X.X```

4) Run the installation:

```pip install .```

**Installing into a Virtual Environment**

This can be useful if root privalages are not available to the user.

1) Create a virutal environment:

We recommend using virtualenv
```pip install virtualenv```

2) Create and activate the environment:

```virtualenv stepRNA_env```

```source /stepRNA_env/bin/activate```

3) Install with pip:

```pip install stepRNA```

**Alternative:**

Download the GitHub repository then make the script executable with:
```chmod +x bin/stepRNA```

Then either:
- Make a symbolic link for stepRNA.py to your bin directory
```ln FULL_PATH_TO_stepRNA.py FULL_PATH_TO_BIN```
- Direct PATH to the directory with stepRNA in it
```export PATH=$PATH:FULL_PATH_TO_GitHubRepo```

## Use:

See the documentation for a detailed description on how to use stepRNA ([MANUAL](URL_LINK))

The quickest way to use stepRNA:

```stepRNA --reference REFERENCE --reads QUERY```
 
This will align the reads to the reference sequences and output into the current diretory using the READS filename as the prefix. **REFERENCE** and **QUERY** must have unique FASTA headers (if not use ```-u```)

Helpful options:
- ```--name``` can be used to customise the prefix name
- ```--directory``` can be used to specify an output directory

## Example:

Using the reads from stepRNA/example_data:

```stepRNA --reference stepRNA/example_data/FILENAME -- reads stepRNA/example_data/FILENAME --directory stepRNA_example```

This will create a new direxotry called *stepRNA_example* that contains:
- AligmentFile/; a directory containing BAM files for each of the overhang lengths
- CSVs; containing count information (see [MANUAL](URL_LINK) for more information)

## News

Latest release notes:

**Version 1.0.0 - 25 Feb, 2020**

- First public release (v1.0.0)
- Caveats: Manual is still incomplete, untested on Windows
- Features such as log file outputs still to be finalised


See [NEWS](URL_LINK) for historical updates of release notes

## Licence

stepRNA is licensed under the MIT license.  See [LICENSE](URL_LINK) file for details.

## Additional Information

For more information:
- Go to [FAQs](URL_LINK) to see commonly asked quesitons
- Loo at the [USER MANUAL](URL_LINK) to see detailed instructions and all of the available options

If you use stepRNA in your work please cite:

[stepRNA: identification of Dicer-processing signatures and passenger strand lengths in small RNA (Murcott, Pawluk & Hunt, 2021)](URL_LINK)
