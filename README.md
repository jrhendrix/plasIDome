# plasIDome
A tool to detect plasmids and contamination in bacterial and archaeal genome assemblies

## Getting Started
### Requirements
* python 3.7+
* blast+ 2.10.1+

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
