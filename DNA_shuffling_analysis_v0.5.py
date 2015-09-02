import subprocess as sp
import svgwrite
from operator import itemgetter
import sys
import os

# color lookup table to color code by parental sequence 
color_palette = {
				'aldh1':'rgb(245,138,94)',
				'aldh2':'rgb(250,172,97)',
				'aldh3':'rgb(255,239,134)',
				'aldh4':'rgb(248,211,169)',
				'aldh5':'rgb(177,255,103)',
				'aldh6':'rgb(117,198,169)',
				'aldh7':'rgb(183,230,215)',
				'aldh8':'rgb(133,218,233)',
				'aldh9':'rgb(132,176,220)',
				'aldh10':'rgb(158,175,210)',
				'aldh11':'rgb(199,176,227)',
				'aldh12':'rgb(255,156,205)',
				'aldh13':'rgb(214,178,149)',
				'aldh14':'rgb(213,150,135)',
				'aldh15':'rgb(180,171,172)',
				'aldh16':'rgb(198,201,209)'}


def fasta_length_parser(fastafilename):
	#current_gene = ""   # Start with an empty string, just in case
	genes = { }         # Make an empty dictionary of genes
	try:
		fh = open(fastafilename, 'r')
	except IOError:
		print 'Could not find file with filename %s' % (fastafilename)
		result = 'Please verify that your filename is correct and try again.'
		return result
	for lineInd, line in enumerate(fh.readlines()):
		if lineInd == 0:
			if not line[0] == '>':
				print 'File does not conform to FASTA format.'
				result = 'Please try again with FASTA formatted file.'
				fh.close( )
				return result
			else:
				pass
		else:
			pass
		line = line.strip()  # Clear out leading/trailing whitespace
		line = line.upper()  # Deals with whatever case the
							 # sequence is by making it all upper case
		if len(line) > 0 and line[0] == ">":   # This one is a new gene
			seq_name = line[1:]
			#genes[current_gene] = ""
		else:                # Add onto the current gene
			seq_length = len(line)
	fh.close()
	
	seq_info = seq_name, seq_length
	return seq_info


def shuffled_blocks_analysis(fasta_file, seq_name, seq_length):
	# some arguments for running BLAST
	# aldh.fsa database comprised of aldh1-16
	# ungapped and mismatch penalty of -15 give desired blast alignments so far
	program = 'blastn'
	queryseq = fasta_file
	database = 'aldh.fsa'
	gap_mode = '-ungapped'
	penalty = '-15'
	outfmt = '6'
	 
	# run the blast search as a process and capture the output 
	proc = sp.Popen([program, '-query', queryseq, '-db', database, gap_mode, '-penalty', penalty,
	'-outfmt', outfmt], stdout=sp.PIPE)
	output = proc.communicate()

	# split the blast output by newlines and remove the empty final line
	outlist = output[0].split('\n')[:-1]

	# empty list to hold the aligned sequence blocks
	seq_blocks = []

	# read the blast output line by line, split into a list by tabs, and capture the parental sequence,
	# length, and start/end position of each block
	for line in range(len(outlist)):
		out_line = outlist[line].split('\t')
		#seq_name = out_line[0]
		seq_blocks.append([out_line[1], int(out_line[3]), int(out_line[6]), int(out_line[7])])

	# sort the seq_blocks by size, small to large
	seq_blocks = sorted(seq_blocks, key = itemgetter(1))

	print 'sorted by size'
	for x in range(len(seq_blocks)):
		print seq_blocks[x]
	print '\n'

	blocks_to_filter = []

	# save the start and end of a block (block A), starting with the largest
	# then loop through the blocks a second time (block B), if block A encompasses block B,
	# B is added to a list to be filtered
	for block in range(len(seq_blocks)):
		start = seq_blocks[block][2]
		end = seq_blocks[block][3]
		for block in range(len(seq_blocks)):
			if start == seq_blocks[block][2] and end == seq_blocks[block][3]:
				continue
			elif start <= seq_blocks[block][2] and end >= seq_blocks[block][3]:
				blocks_to_filter.append(seq_blocks[block])
			else:
				continue

	# create a filtered blocks list by subtracting out blocks from the filtered list
	seq_blocks_filtered = [x for x in seq_blocks if x not in blocks_to_filter]

	# sort blocks from 5' to 3' to prepare for resolving overlaps
	seq_blocks_filtered = sorted(seq_blocks_filtered, key = itemgetter(2))

	print 'sorted by position'
	for x in range(len(seq_blocks_filtered)):
		print seq_blocks_filtered[x]
	print '\n'

	# loop through the blocks (making sure to stop at the last block) 
	# if block A ends after the start of block B, update the start/end of each block to be the average of the overlap
	# also update the length of the block (used for making the figure)
	for block in range(len(seq_blocks_filtered)):
		#end = seq_blocks_filtered[block][3]
		if block < len(seq_blocks_filtered) - 1:
			if seq_blocks_filtered[block][3] > seq_blocks_filtered[block + 1][2]:
				junction_position = (seq_blocks_filtered[block][3] + seq_blocks_filtered[block + 1][2]) / 2
				seq_blocks_filtered[block][3] = junction_position - 1
				seq_blocks_filtered[block][1] = seq_blocks_filtered[block][3] - seq_blocks_filtered[block][2] + 1
				seq_blocks_filtered[block + 1][2] = junction_position
				seq_blocks_filtered[block + 1][1] = seq_blocks_filtered[block + 1][3] - seq_blocks_filtered[block + 1][2] + 1

	# initialize the svg file with a filename and resolution
	svg_document = svgwrite.Drawing(filename = seq_name + "_v0.5.svg",
									size = (str(seq_length + 100) + "px", "48px"))

	# add a black bar to represent the full length gene
	svg_document.add(svg_document.rect(insert = (100, 15),
										size = (str(seq_length) + "px", "16px"),
										fill = "black"))

	# draw a rectangle of the correct size and shape for each sequence block
	offset = 0
	for block in range(len(seq_blocks_filtered)):
		svg_document.add(svg_document.rect(insert = (seq_blocks_filtered[block][2] + 100, offset),
											size = (str(seq_blocks_filtered[block][1]) + 'px', "48px"),
											stroke_width = "1",
											stroke = "black",
											fill = color_palette[seq_blocks_filtered[block][0]]))
		#offset += 12

	svg_document.add(svg_document.text(seq_name, insert = (5, 30)))

	svg_document.save()
	print '%s done' % fasta_file


cwdfiles = os.listdir('.')

for cwdfile in cwdfiles:
	if cwdfile.endswith('.fasta'):
		#print cwdfile
		seq_info = fasta_length_parser(cwdfile)
		shuffled_blocks_analysis(cwdfile, seq_info[0], seq_info[1])
	
