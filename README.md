# stepRNA
[![PyPI version](https://badge.fury.io/py/stepRNA.svg)](https://badge.fury.io/py/stepRNA)
## Overview

stepRNA is a RNA-seq read processor, based on bowtie2, that will align small RNA-seq query reads (passenger sequences) to reference sequences (siRNAs) and output information about the type and length of overhangs uncovered. It was originally developed for uncovering Dicer processing signatures but is not limited to just this.

**Table of Contents:**
- [Reporting Issues](#Reporting-Issues)
- [Installation](#Installation)
- [Use](#Use)
- [Example](#Example)
- [News](#News)
- [Licence](#Licence)
- [Additional Information](#Additional-Information)

## Reporting Issues

Please report any issues to the stepRNA GitHub page or via email:
- https://github.com/bmm514/stepRNA/issues
- benmurcott96@gmail.com (Ben Murcott)
- v.l.hunt@bath.ac.uk (Vicky Hunt)

## Installation

In order for stepRNA to run you must have:
- Bowtie2 >= v2.3.4 (see [BOWTIE2 website](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for more details)
- biopython >= v0.3.0 (automatically installed with pip)
- numpy >= v1.19.0 (autoamtically installed with pip)
- pysam >= v0.16.0.1 (autoamtically installed with pip)

**To install with pip (recommended):**

```pip install stepRNA```

If this doesn't work due to non-root access issues see **Installing into a Virtual Environment**

**To install from source:**

1) Download [stepRNA-X.X.X.tar.gz](https://pypi.org/project/stepRNA/#files)
2) Unzip it
3) Move into stepRNA-X.X.X repository
4) Run the installation

```
tar -xvzf stepRNA-X.X.X.tar.gz

cd stepRNA-X.X.X

python3 setup.py install
```

**Installing into a Virtual Environment**

This can be useful if root privalages are not available to the user.

1) Create a virutal environment (recommend using virtualenv):

```pip install virtualenv```

2) Create and activate the environment:
3) Run the installation 

```
virtualenv stepRNA_env

source /stepRNA_env/bin/activate

pip install stepRNA
```

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

```stepRNA --reference stepRNA/example_data/26G_embryo.fa --reads stepRNA/example_data/LF_embryo.fa --directory stepRNA_example```

This will create a new direcotry called *stepRNA_example* that contains:
- LF_embryo_AligmentFile/; a directory containing BAM files for each of the overhang lengths
- CSVs; containing count information

See the [MANUAL](URL_LINK) for more information

## News

Latest release notes:

**Version 1.0.0 - 25 Feb, 2020**

- First public release (v1.0.0)
- Caveats: Manual is still incomplete, untested on Windows
- Features such as log file outputs still to be finalised


See [NEWS](https://github.com/bmm514/stepRNA/blob/master/NEWS.md) for historical updates of release notes

## Licence

stepRNA is licensed under the MIT license.  See [LICENSE](https://github.com/bmm514/stepRNA/blob/master/LICENSE) file for details.

## Additional Information

For more information:
- Go to [FAQs](https://github.com/bmm514/stepRNA/blob/master/FAQs.md) to see commonly asked quesitons
- Loo at the [USER MANUAL](URL_LINK) to see detailed instructions and all of the available options

If you use stepRNA in your work please cite:

[stepRNA: identification of Dicer-processing signatures and passenger strand lengths in small RNA (Murcott, Pawluk & Hunt, 2021)](URL_LINK)
