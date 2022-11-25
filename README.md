# PlasIDome
A tool to detect plasmids and contamination in bacterial and archaeal genome assemblies

## Introduction
In bacterial genome assembly, it is important to account for all of the contigs. PlasIDome finds all contigs in an assembly file that are less than a specified length (default: 200,000bp) and uses sequences homology to categorize each as chromosomal, plasmid, or contamination. With this tool a key step in assessing assembly quality has been automated and can be integrated into existing workflows.


## Getting Started
### Requirements
* python 3.7+
* blast+ 2.10.1+
* bioseq

### Installation
PlasIDome is available on PYPI and can be installed using pip
``` pip install plasIDome ```


## Usage
As input, PlasIDome takes a genome assembly file in fasta format (.fasta, .fa, .fna).

Example:

``` plasidome -b path/to/blastn -f genome.fasta ```


Examine contigs up to 50,000 bp in length:

``` plasidome -b path/to/blastn -f genome.fasta -l 50000 ```


Extensive Usage:

``` plasidome -b path/to/blastn -f genome.fasta -p path/to/output -o output_name ```



## Output Files
PlasIDome generates three outputs:
* report.tsv
* alignment_results.tsv
* directory of contigs

The report.tsv file contains a summary of the alignment results in table format. From left to right, the table includes the contig name, its classification, and its contamination status followed by the number of matches that were to chromosomal, plasmid, undetermined, or human sequences in the database.

A file called alignment_results.tsv is created that contains the raw, uninterpreted alignment results in table format. From left to right, the table includes the contig name, subject name, subject taxonomic ID, the percent sequence identity, query coverage, qcovhsp, alignment length, and e-value. The raw data allows the investigator to parse the alignment results manually. The alignment results written to this file are limited to homology matches where the subject and query shared at least 95% sequence identity and the subject covered at least 95% of the length of the length of the query sequence. 

FASTA files with each contig sequence are saved to a subdirectory called single_contigs so the user can re-align any sequence without having to find and isolate the contig themselves.

At this time, PlasIDome is not able to determine if a plasmid is complete or if multiple contigs are fragments of the same plasmid. PlasIDome also does not comment on the taxonomic identification of the species. If the user wants to see the taxonomic makeup of the sample, it is best to review the raw alignment results in the alignment_results.tsv file.
