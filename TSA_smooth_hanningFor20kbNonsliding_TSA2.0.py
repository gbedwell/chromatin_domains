#!/usr/bin/python
# Programmer : Liguo Zhang, modified from Yang Zhang's code
# Date: 
# Last-modified: 1 Aug 2018
# Used for 20kb non-sliding window, take wig file to smooth and generate smoothed wig, take average for a bigger window and generate wig

import os,sys,argparse
import math 
import numpy as np
from scipy import exp2
from TSA_utility import *

def ParseArg():
    ''' This function parse the argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('--wig',type=str,dest="wig",help="wig file")
    p.add_argument('-g','--genome',type=str,dest="genome",help="human genome file (two column, first col is chromosome name, second chromosome is chromosome size")
    p.add_argument('-w','--window',type=int,dest="window",help="window size of input wig")
    p.add_argument('-aggwin','--aggwin',type=int,dest="aggwin",help="window size used to aggregate data")
    p.add_argument('--smooth',action="store_true",help="if set, will smooth data")
    p.add_argument('-n1','--name1',type=str,dest="name1",help="output smoothed file prefix, wig file will be generated")
    p.add_argument('-n2','--name2',type=str,dest="name2",help="output smoothed and then averaged file prefix, wig file will be generated")
    if len(sys.argv) < 2:
        print p.print_help()
        exit(1)
    return p.parse_args()

def Smooth(x,window_len=21):
    ''' This function defines 'hanning' smoothing method'''
    s = np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    w=eval('np.hanning(window_len)')
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]

def Main():
    global args
    args=ParseArg()
    wout = WriteToFile(args.name1 + '.wig')
    agg_wout = WriteToFile(args.name2 + '.wig')
    genome = LoadGenome(args.genome)
    hash={}
    array_temp=[]
    num=1
    for line in ReadFromFile(args.wig):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'track':
            continue
        if num % 10000 == 0:
            logging("Process: %d line processed in input wig" % (num))
        num += 1
        if row[0] == 'variableStep': # new track chrom begin
            if len(array_temp)!=0:
                if chrom != row[1].split('=')[1]:
                    hash[chrom]=array_temp
                    array_temp=[]
            chrom = row[1].split('=')[1]
        else:
            array_temp.append(float(row[1]))
            hash[chrom]=array_temp
    if args.smooth:
        logging("Options: turn on smooth mode")
    for chr in SortGenome(genome):
        chrom_size = genome[chr]
        logging("Process: %s\t%d" % (chr,chrom_size))
        array=hash[chr]
        array=np.array(array)
        if args.smooth:
            smooth_array = Smooth(array)
        else:
            smooth_array = array
        print >>wout, "variableStep chrom=%s span=%d" % (chr,args.window)
        for nn,value in enumerate(smooth_array):
            if nn == 0: 
                print >>wout, "%d\t%.6f" % (nn+10000,value) 
            elif nn == len(smooth_array) - 1:
                print >>wout, "variableStep chrom=%s span=%d" % (chr,chrom_size-((nn)*args.window)-10000)
                print >>wout, "%d\t%.6f" % (nn*args.window+10000,float(value))
            else:
                print >>wout, "%d\t%.6f" % (nn*args.window+10000,float(value))
        agg_array = []
        start = 0
        stop = 0 + args.aggwin/args.window
        for m in range(len(smooth_array)):
            if stop >= len(smooth_array):
                stop = len(smooth_array)
                agg_array.append(np.mean(smooth_array[start:stop]))
                break
            agg_array.append(np.mean(smooth_array[start:stop]))
            start += args.aggwin/args.window
            stop += args.aggwin/args.window
        agg_array = np.array(agg_array)
        print >>agg_wout, "variableStep chrom=%s span=%d" % (chr,args.aggwin)
        for mm,value in enumerate(agg_array):
            if mm == 0: 
                print >>agg_wout, "%d\t%.6f" % (mm+10000,value) 
            elif mm == len(agg_array) - 1:
                print >>agg_wout, "variableStep chrom=%s span=%d" % (chr,chrom_size-((mm)*args.aggwin)-10000)
                print >>agg_wout, "%d\t%.6f" % (mm*args.aggwin+10000,float(value))
            else:
                print >>agg_wout, "%d\t%.6f" % (mm*args.aggwin+10000,float(value))    
    wout.flush()
    agg_wout.flush()
    wig2bw1 = "/mnt/storage/apps/UCSC-tools/wigToBigWig -clip %s %s %s" % (args.name1 + '.wig', args.genome, args.name1 + '.bw')
    os.system(wig2bw1)
    wig2bw2 = "/mnt/storage/apps/UCSC-tools/wigToBigWig -clip %s %s %s" % (args.name2 + '.wig', args.genome, args.name2 + '.bw')
    os.system(wig2bw2)
    logging("Finish: TSA_smooth DONE!!!")

if __name__=="__main__":
    Main()
