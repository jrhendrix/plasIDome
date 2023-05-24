
'''
FILE:	plasIDome
AUTHOR:	J.R. Hendrix
URL: 	http://stronglab.org
		https://github.com/jrhendrix/plasIDome
DESC:	This script extracts contigs less than a certain length 
		and aligns the sequence using blastn
		and reports hits to chromosomes and plasmids
'''

# IMPORT FROM PYTHON STANDARD LIBRARY
import argparse
import os
import subprocess
import sys

from Bio import SeqIO	# Source: https://biopython.org/wiki/SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class File:
	""" Base class for all file-types """

	def __init__(self, path, file_type=None):
		self._path = None
		self.path = path
		self.file_type = file_type

	@property
	def path(self):
		return self._path
	
	@path.setter
	def path(self, value):
		if not os.path.isabs(value):
			value = os.path.join(os.getcwd(), value)
		if os.path.isfile(value):
			self._path = value
		else:
			raise FileNotFoundError(value)

	@property
	def filename(self):
		return os.path.basename(self.path)

	@property
	def file_prefix(self):
		return self.filename.split(".")[0]

	@property
	def extension(self):
		return self.filename.split(".")[-1]


class Dir:
	""" Base class for system directories """

	def __init__(self, path):
		self._path = None
		self.path = path

	@property
	def path(self):
		return self._path
	
	@path.setter
	def path(self, value):
		if not os.path.isabs(value):
			value = os.path.join(os.getcwd(), value)
		if os.path.isdir(value):
			self._path = value
		else:
			raise NotADirectoryError(value)

	@property
	def dirname(self):
		return self.path.strip("/").split("/")[-1]

	@property
	def children(self):
		children = [Dir(os.path.join(self.path, subdir)) 
			for subdir in os.listdir(self.path) 
			if os.path.isdir(os.path.join(self.path, subdir))]
		if len(children) > 0:
			return children
		else:
			return None

	@property
	def files(self):
		files = [File(os.path.join(self.path, file))
			for file in os.listdir(self.path)
			if os.path.isfile(os.path.join(self.path, file))]
		if len(files) > 0:
			return files
		else:
			return None

	def join(self, *args):
		return os.path.join(self.path, *args)

	def make_subdir(self, *args):
		""" Makes recursive subdirectories from 'os.path.join' like arguments """
		subdir = self.join(*args)
		return self.make(subdir)

	@classmethod
	def make(cls, path):
		try:
			os.makedirs(path)
			return cls(path)
		except FileExistsError:
			return cls(path)

	def __repr__(self):
		return self.path
	


def check_input(args):
	''' Check input file for correct extension '''

	# CHECK IF FILE EXISTS
	try:
		f = File(args.fasta)
	except IOError:
		print("ERROR: Could not find file")
		return 1

	# CHECK THAT FILE IS FASTA FORMAT
	extensions = ('fasta', 'fa', 'fna', 'faa')
	if f.extension not in extensions:
		print("ERROR: Input file was not in FASTA format.")
		return 2

	return 0


def get_contigs(args, BASEDIR):
	# OUTPUT CONTIGS OF INTEREST
	ccount = 0
	conDir = BASEDIR.make_subdir('single_contigs')
	for record in SeqIO.parse(args.fasta_input, "fasta"):
		if len(record.seq) <= args.length:

				out = 'contig_', record.id, '.fasta'
				outName = ''.join(out)
				outPath = '/'.join((conDir.path, outName))
				SeqIO.write(record, outPath, "fasta")
				ccount = ccount + 1

	if ccount < 1:
		print('No contigs met the length requirement. Done.')
		exit()

	return conDir, ccount


def main(program):
	''' The main worker function to dictate processes '''

	cwd = os.getcwd()

	# PARSER : ROOT
	parent_parser = argparse.ArgumentParser(prog='plasidome')
	parent_parser.add_argument('-b', '--blastn_path', default="/Strong/proj/shared_code/ncbi-blast-2.10.1+/bin/blastn", help='Path to blastn')
	parent_parser.add_argument('-f', '--fasta_input', help='Path to fasta file')
	parent_parser.add_argument('-l', '--length', default=200000, help='Contigs shorter than this length will be tested', type=int)
	parent_parser.add_argument('-o', '--out_directory', default="contig_assignments", help='Prefix of output directory', type=str)
	parent_parser.add_argument('-p', '--path_to_output',  default=cwd, help='Path to output', type=str)
	parent_parser.add_argument('-r', '--report_file', default='report', help='Name of output summary report')

	args = parent_parser.parse_args()

	
	# CREATE OUTPUT STRUCTURE
	try:
		TOP_DIR = Dir(args.path_to_output)
		BASEDIR = TOP_DIR.make_subdir(args.out_directory)
	except IOError:
		print('ERROR: could not establish output directory.')

	# GET CONTIGS THAT MET LENGTH THRESHOLD
	conDir, ccount = get_contigs(args, BASEDIR)

	# ESTABLISH RAW ALIGNMENT OUTPUT TABLE
	outTab = '/'.join((BASEDIR.path, 'alignment_results.tsv'))
	f1 = open(outTab, 'w')
	header = 'contig', 'staxids', 'title', 'percent_ident', 'query_coverage', 'qcovhsp', 'length', 'e_value'
	head = '\t'.join(header) + '\n'
	f1.write(head)
	
	notFound = 'no alignments met required parameters'

	# RUN BLASTN and REPORT RESULTS
	print(f'There are {str(ccount)} contigs to align')
	path = args.blastn_path
	outfmt = '6 qseqid staxids stitle pident qcovs qcovhsp length evalue'
	count = 0
	for file in conDir.files:
		count = count + 1

		seqName = file.file_prefix
		print(f'File #{str(count)}: {seqName}')
		seqID = seqName.split('_')[0]

		# call BLAST remotely
		command = [path, '-remote', '-db', 'nr', '-query', file.path, '-outfmt', outfmt, '-perc_identity', '95', '-qcov_hsp_perc', '95', '-max_target_seqs', '5']
		process = subprocess.run(command, capture_output=True)
		result = process.stdout.decode("utf-8")
		result = ''.join(('\n', result))
		#print(len(result))
		#print(result)

		if len(result) < 2:
			result = '\t'.join((seqName, notFound)) + '\n\n'

		f1.write(result)

	f1.close()
	

	# EVALUATE BLAST RESULT
	print('Evaluating BLAST results')
	f2 = open(outTab, 'r')
	f2.readline() # skip header
	sumDic = {}
	chrom = 0
	plasmid = 0
	plastid = 0
	und = 0
	for line in f2:
		line = line.strip()
		#print('\n')
		#print('line: ', line)
		if len(line) == 0:
			continue

		l = line.split('\t')
		

		key = l[0]

		if key not in sumDic:
			sumDic[key] = {}
			sumDic[key]['plas'] = 0
			sumDic[key]['chrom'] = 0
			sumDic[key]['und'] = 0
			sumDic[key]['contam'] = 0
			sumDic[key]['notfound'] = 0

		#print(l)
		taxid = str(l[1])
		#print('taxid: ', taxid, notFound)
		if taxid == notFound:
			#print('\tmatch')
			sumDic[key]['notfound'] += 1
			continue


		title = l[2]

		# Counts matches to human samples
		if taxid == "9606":				# Detect companination (limited to homo sapien)
			sumDic[key]['contam'] += 1
			continue

		# Counts matches to plasmid samples
		if "plasmid" in title:			# Count number of plasmid annotations
			sumDic[key]['plas'] += 1

		# Counts matches to chromosomal samples
		if "chromosome" in title:		# Count number of chromosome annotations
			sumDic[key]['chrom'] += 1
		elif "plasmid" not in title:		# Count undetermined
			sumDic[key]['und'] += 1
	f2.close()

	# REPORT SUMMARY OF RESULTS
	print('Reporting summary of results')
	rep = BASEDIR.path, '/', args.report_file, '.tsv'

	repFile = ''.join(rep)
	f3 = open(repFile, 'w')
	header = 'contig', 'classification', 'contaminated', 'chromosome', 'plasmid', 'undetermined', 'human'
	head = ''.join(('\t'.join(header), '\n'))
	f3.write(head)
	for key in sumDic:
		con = sumDic[key]['contam']
		c = sumDic[key]['chrom']
		p = sumDic[key]['plas']
		u = sumDic[key]['und']
		n = sumDic[key]['notfound']

		#r = c/p    # Note to check for denominator = 0
		# Failure -> undetermined
		if n > c+p+u:
			ment = 'no sig. hits'
		#elif c == p & p == u & c == 0: # Novel plasmid - no hits
		#	ment = 'novel'
		elif c > (p + u):
			ment = 'chromosome'
		elif p > (c + u):
			ment = 'plasmid'
		else:
			ment = 'undetermined'

		#  Override if  contamination
		if con > (c+u+p):
			ment = 'contamination'
		

		if sumDic[key]['contam'] > 0:
			contamination = True
		else:
			contamination = False

		entry = str(key), ment, str(contamination), str(c), str(p), str(u), str(con)
		record = ''.join(('\t'.join(entry), '\n'))
		f3.write(record)

	f3.close()

	print("Done.")



if __name__== "__main__":
	main(sys.argv[1])









