#!/usr/bin/python
# Programmer : Yang Zhang
# Contact: yzhan116@illinois.edu
# Last-modified: 01 Feb 2023 04:21:51 PM (YZ)
# Last-modified: 23 Feb 2023 15:08:45 (GJB)

import os,sys,argparse
from math import log
import numpy as np
import pysam
import re

def parse_argument():
    ''' This Function Parse the Argument '''
    p = argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog = 'Library dependency : pysam')
    p.add_argument('-v','--version',action = 'version', version = '%(prog)s 0.5')
    p.add_argument('-r','--res',dest = "res", type = int, default = 10, help = "wig file resolution (bp)")
    p.add_argument('-w','--win',dest = "win", type = int, default = 1000,  help = "Sliding window step size (bp)")
    p.add_argument('-e','--exp',dest = "bamexp", type = str, help = "pulldown group bam file")
    p.add_argument('-c','--con',dest = "bamcon", type = str, help = "control group bam file")
    p.add_argument('-g','--genome',dest = "genome", type = str, help = "genome chromosome size file (used for wigToBigWig)")
    p.add_argument('-o','--output',dest = "output", type = str, help = "the prefix of output wig files without .wig, eg. ") 
    p.add_argument('--wig2bw',dest = "wig2bw", type = str, help = "program location (full path eg. /home/zocean/wigToBigWig of wigToBigWig program, 'sys' means wigToBigWig is the system path")
    if len(sys.argv)<2:
        sys.exit(p.print_help())
    return p.parse_args()

##########
# Utility
##########

def logging(text):
    print >>sys.stderr, "Logging: " + text

def warning(text):
    print >>sys.stderr, "Warning: " + text

def error(text):
    print >>sys.stderr, "Error: " + text

############
# Functions
############

def filterPick(list,filter):
    '''
    return the items match the filter pattern
    '''
    return [ ( l ) for l in list for m in (filter(l),) if m]

def delete_key(dic, key_list):
    for key in key_list:
        del dic[key]
    return dic

def BamToBin(bamfile,binsize):
    '''
    This function counts the number of reads in each bin of each chromosome. Bin number is determined by resolution (-r).
    '''
    # Create a hash dictionary using the chromosome names in the header of the bam file.
    bin_hash = {}
    # Read the indexed bam file using pysam
    samfile = pysam.Samfile(bamfile,"rb")
    # Create the chr_table dictionary. Store the chromosome name and length.
    chr_table = {}
    for nn in range(len(samfile.references)):
        chr_table[samfile.references[nn]] = samfile.lengths[nn]
    # Remove user-defined chromosomes
    chr_removed_list = filterPick(chr_table.keys(), chr_filter)
    chr_table = delete_key(chr_table, chr_removed_list)
    # Calculate the total number of mapped reads
    total_mapped = 0
    idx_stats = samfile.get_index_statistics()
    for idx in idx_stats:
        if idx.contig in chr_table.keys():
            total_mapped += idx.mapped
    logging("Total mapped reads is: %d" % (total_mapped))
    # Initialize the read count table.
    # If the remainder of chromosome length / bin size is not 0, add one bin to account for the extra
    for chrom, length in chr_table.items():
        bin_num = length/binsize
        if length%binsize != 0:
            bin_num += 1
        # Initiate the count in each bin with 0's
        bin_hash[chrom] = [0.000000 for row in range(bin_num)]
    # look through every read in bamfile
    n = 0
    last_percent = 0
    for alignment in samfile:
        # Skip unmapped reads
        if alignment.tid < 0:
            continue
        # Define position as the alignment position + the aligned fragment length.
        # Testing reference_start and reference_length instead of the deprecated alen
        # pos = alignment.pos + alignment.alen
        pos = alignment.reference_start + alignment.reference_length
        # Get the chromosome name for the corresponding alignment id.
        chr_name = samfile.getrname(alignment.tid)
        if chr_table.get(chr_name, None) is None:
            continue
        try:
            # Increase the read count in bin number pos/binsize for chromosome chr_name
            # pos/binsize is the index. Python truncates floats as indeces towards zero (floor operation)
            bin_hash[chr_name][pos/binsize] += 1.0
        except KeyError:
            warning("read position %d is larger than chromosome size %d" % (pos, chr_table[chr_name]))
            continue
        # Add processing information/progress bar.
        n += 1
        percent = float(n) / float(total_mapped)*100
        if int(percent) - last_percent == 2:
            print >>sys.stderr, "="*int(int(percent) / 100.0 * 80.0) + '>' + '{:.2f}'.format(percent) + '%' + '\r',
            last_percent = int(percent)
    print >>sys.stderr, ""
    samfile.close()
    return bin_hash

def BinToWin(bin_hash,step):
    '''
    This function calculates window counts by summing over bin counts. The "step" argument tells the function how many bins to combine.
    '''
    # Get the chromosome name from bin_hash
    chr_list = bin_hash.keys()
    # Create a hash dictionary for the counts in a given window.
    win_hash = {}
    # If step == 1, window counts are just bin counts
    if step == 1:
        for chrom in chr_list:
            win_hash[chrom] = list(bin_hash[chrom])
    # If step > 1, step bins are combined to make a window. Bin counts are summed. 
    else:
        for chrom in chr_list:
            # Create window hash dictionary and sum step bin counts together to get window counts.
            win_no = len(bin_hash[chrom]) - step + 1
            # Initiation win_hash with 0's.
            win_hash[chrom] = [0 for n in range(win_no)]
            count_list = bin_hash[chrom]
            # Sum the counts in bins [0:step).
            buffer = sum(count_list[0:step])
            # Assign index 0 in win_hash the value buffer
            win_hash[chrom][0] = buffer
            # Repeat the calculation over all other win_no indeces.
            for j in range(1,win_no):
                win_hash[chrom][j] = buffer + count_list[j - 1 + step] - count_list[j - 1]
                buffer = win_hash[chrom][j] # Update buffer
    return win_hash

# Calculate the average value of an array.
def CalAve(table):
    value_list = table.values()
    array = [ aa for dd in value_list for aa in dd if aa is not None ]
    return np.nanmean(array)

# Normalize window counts.
def WinNormalize(exp_win, con_win):
    '''
    Normalize pulldown to input.
    '''
    chrs_A = exp_win.keys()
    chrs_B = con_win.keys()
    # Verify that the chromosome names are the same in both conditions.
    if sorted(chrs_A) != sorted(chrs_B):
        print "The chromosome lists for the two bam files are different. Continue? (Yes/No)"
        go_on = raw_input("> ")
        if go_on == "No":
            sys.exit("User interrupt.")
        elif go_on == "Yes":
            pass
        else:
            sys.exit("Unknown answer.")
    # Calculate the average number of reads across all control windows.
    # This is used in count normalization.
    con_ave_win = CalAve(con_win)
    # Create list of the chromosome names in common between A and B.
    chrom_list = list(set(chrs_A) & set(chrs_B))
    # Create norm_win dictionary
    norm_win = {}
    # Populate norm_win, first with 0's, then with normalized values.
    for chrom in chrom_list:
        norm_win[chrom] = [0 for nn in range(len(exp_win[chrom]))]
        out_list = norm_win[chrom]
        # Experimental values
        exp_list = exp_win[chrom]
        # Control values
        con_list = con_win[chrom]
        for nn in range(len(out_list)):
            if con_list[nn] > 1e-6 and exp_list[nn] > 1e-6:
                # This is the actual normalization step.
                # experimental value * ( avg control value across all bins / control value )
                out_list[nn] = (exp_list[nn])*(con_ave_win/(con_list[nn]))
            else:
                # No reads in either sample.
                out_list[nn] = None
    return norm_win

def GetRatio(norm_win):
    '''
    Calculate the ratio of normalized counts to average normalized counts.
    '''
    table = {}
    ave_win = float(CalAve(norm_win))
    for chrom in norm_win.keys():
        table[chrom] = []
        for value in norm_win[chrom]:
            if value is not None:
                table[chrom].append(float(value) / ave_win)
            else:
                # If None, assign a value of 1. This becomes 0 on a log-scale.
                table[chrom].append(1.0)
    return table

def WriteWig(output_file, hash_table, genome_table, resolution, step):
    '''
    Write normalized counts to wig file.
    '''
    win_half = int(resolution * step / 2)
    chrs = hash_table.keys()
    chrs.sort()
    fo = open(output_file + ".wig", "w")

    for ii in range(len(chrs)):
        chr_list = hash_table[chrs[ii]]
        is_end = False
        print >>fo, "variableStep chrom=" + chrs[ii] + " span=" + str(resolution)
        for jj in range(len(chr_list)):
            if is_end:
                break
            if jj * resolution + win_half + resolution > genome_table[chrs[ii]]: # Last window
                new_span = genome_table[chrs[ii]] - jj * resolution - win_half
                print >>fo, "variableStep chrom=" + chrs[ii] + " span=" + str(new_span)
                is_end = True
            # Conversion to log2-scale happens here!
            print >>fo, "%d\t%.6f" % (jj * resolution + win_half, log(chr_list[jj],2))
    fo.flush()
    fo.close()

def WriteBedGraph(output_file, hash_table, genome_table, resolution, step):
    '''
    Write bedgraph file.
    '''
    win_half = int(resolution * step / 2)
    chrs = hash_table.keys()
    chrs.sort()
    fo = open(output_file + ".bedgraph", "w")
    for ii in range(len(chrs)):
        chr_list = hash_table[chrs[ii]]
        for jj in range(len(chr_list)):
            if jj * resolution + win_half > genome_table[chrs[ii]]:
                warning("window position %d is larger than chromosome size %d for chromosome %s" % (jj * resolution + win_half, genome_table[chrs[ii]], chrs[ii]))
                continue
            if jj * resolution + win_half + resolution > genome_table[chrs[ii]]: # Last window
                new_span = jj * resolution + win_half + resolution - genome_table[chrs[ii]]
                print >>fo, "%s\t%s\t%s\t%.6f" % (chrs[ii], jj * resolution + win_half, jj * resolution + win_half + new_span, chr_list[jj])
            else:
                print >>fo, "%s\t%s\t%s\t%.6f" % (chrs[ii], jj * resolution + win_half, jj * resolution + win_half + resolution, chr_list[jj])
    fo.flush()
    fo.close()

def LoadGenome(filename):
    genome = {}
    try:
        fin = open(filename, "r")
    except IOError:
        error("Can't open file: %s" % (filename))
        exit()
    for line in fin:
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        chrom = row[0]
        size = int(row[1])
        genome[chrom] = size
    return genome

def ReportOptions():
    text = "# TSA-seq_normalize.py version: 0.2\n"
    text += "# experimental bam file: %s\n" % (args.bamexp)
    text += "# control bam file: %s\n" % (args.bamcon)
    text += "# size of sliding window: %d\n" % (args.win)
    text += "# resolution of output wig files: %d\n" % (args.res)
    text += "# genome chromosome size file: %s\n" % (args.genome)
    text += "# output file name prefix: %s" % (args.output)
    print >>sys.stderr, text

def Main():
    global args
    global chr_filter 
    chr_filter = re.compile('(random|chrM|hap|Un)').search
    args = parse_argument()
    # Report parameters
    ReportOptions()
    # Check that given parameters are within the allowed range(s)
    if (args.res > (args.win / 2)) and args.res != args.win:
        sys.exit("Resolution cannot be larger than half the size of the sliding window.")
    elif (args.res < 1):
        sys.exit("Resolution cannot be less than 1.")
    else:
        if args.win % args.res == 0:
            step = args.win / args.res
        else:
            error("Window size (dividend) must be evenly divisible by resolution (divisor).")
            exit(1)
    # Load genome file
    genome_table = LoadGenome(args.genome)
    # Begin pipeline
    logging("*************** Begin Analysis ***************") 
    # Step 1: Allocate control reads into bins and aggregate bins into windows.
    logging("Step 1: Count the number of reads in control file.")
    logging("Step 1.1 Allocate control reads into bins.")
    Con_bins = BamToBin(args.bamcon, args.res)
    logging("Step 1.2: Aggregate bins into windows.")
    Con_win = BinToWin(Con_bins, step)
    # Step 2: Allocate pulldown reads into bins and aggregate bins into windows.
    logging("Step 2: Count the number of reads in pulldown file.")
    logging("Step 2.1: Allocate pulldown reads into bins.")
    Exp_bins = BamToBin(args.bamexp, args.res)
    logging("Step 2.2 Aggregate bins into windows.")
    Exp_win = BinToWin(Exp_bins, step)
    # Release memory
    del Exp_bins
    del Con_bins
    # Step 3: Normalize pulldown counts to input.
    logging("Step 3: Normalize pulldown counts using input control.")
    norm_win = WinNormalize(Exp_win, Con_win)
    # Release memory
    del Exp_win
    del Con_win
    # Step 4: Calculate ratio
    logging("Step 4: Calculate enrichment ratio in pulldown relative to control.")
    ratio_win = GetRatio(norm_win)
    # Step 5: Write wig and bigwig files.
    logging("Step 5: Write log2(ratio) normalized read counts to wig file.")
    WriteWig(args.output, ratio_win, genome_table, args.res, step)
    if args.wig2bw is not None:
        logging("Step 6: Convert wig to bigwig.")
        if args.wig2bw == 'sys':
            os.system("wigToBigWig %s %s %s" % (args.output+".wig", args.genome, args.output+".bw"))
        else:
            os.system(args.wig2bw + " %s %s %s" % (args.output+".wig", args.genome, args.output+".bw")) 
    logging("*************** Done! ***************")

if __name__=="__main__":
    Main()
