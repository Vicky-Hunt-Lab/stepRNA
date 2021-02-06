# stepRNA

stepRNA is a RNA-seq read processor, based on bowtie2, that will align query reads to reference sequences and output information about the type and length of overhangs uncovered. It was originally developed for processing DICER RNAs but is not limited to just this.

#### Installation

(INCLUDE ONCE I HAVE SORTED)

#### Use:
See the documentation for a detailed description on how to use stepRNA.
The quickest way to use stepRNA:

```stepRNA --reference REFERENCE --reads QUERY```
 
This will align the reads to the reference sequences and output into the current diretory using the READS filename as the prefix. The following should be met:
- **REFERENCE** and **QUERY** must have unique FASTA headers (if not use ```-u```)
- If **REFERENCE** have been extracted from **QUERY** then it is advised to remove the exact matches from the QUERY (use ```-e```) - ***THIS TAKES A LONG TIME!***

Helpful options:
- ```--name``` can be used to customise the prefix name
- ```--directory``` can be used to specify an output directory
