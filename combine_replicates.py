#!/usr/bin/python

import os, sys, argparse
from collections import OrderedDict
from TSA_utility import *

def ParseArg():
    ''' This function parse the argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency:')
    p.add_argument('--wig1', type=str,dest="wig1",help="First wig file")
    p.add_argument('--wig2',type=str,dest="wig2",help="Second wig file")
    p.add_argument('-n','--name',type=str,dest="name",help="name prefix of the combined file")
    if len(sys.argv) < 3:
        print p.print_help()
        exit(1)
    return p.parse_args()

def parse_wig(filename):
    dd = OrderedDict()
    for line in ReadFromFile(filename):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'variableStep': # new track chrom begin
            chrom = row[1].split('=')[1]
            span = int(row[2].split('=')[1])
            cs = (chrom + "_" + str(span))
            if cs not in dd:
                dd[cs] = []
        else:
            coordinate = int(row[0])
            score = float(row[1])
            dd[cs].append((coordinate, score))
            
    return(dd)

def average_data(dd1, dd2):
    avg_dd = {}
    
    for key in dd1.keys():
        if key in dd2:
            dat1 = dd1[key]
            dat2 = dd2[key]
            
            if(len(dat1) == len(dat2)):
                data_len = len(dat1)
                avg_dat = []
                    
                for i in range(data_len):
                    coord1, score1 = dat1[i]
                    coord2, score2 = dat2[i]              
                    average_score = (score1 + score2) / 2.0
                    avg_dat.append((coord1, average_score))           
                    avg_dd[key] = avg_dat
            else:
                sys.exit("Data for " + key + " not the same lengths.")
        else:
            sys.exit("Key " + key + " not in both wig files.")
            
    return( avg_dd )

def Main():
    global args
    args=ParseArg()
    wout = WriteToFile(args.name + '.wig')
    dict1 = parse_wig(args.wig1)
    dict2 = parse_wig(args.wig2)
    avg_dict = average_data(dict1, dict2)

    ordered_keys = dict1.keys()
    for key in ordered_keys:
        if key in avg_dict:
            chrom, span = key.split("_")
            print >> wout, "variableStep chrom=%s span=%d" % (chrom, int(span))
            values = avg_dict[key]
            for coord, score in values:
                print >> wout, "%d\t%f" % (coord, score)

if __name__ == "__main__":
    Main()
