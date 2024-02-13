#!/usr/bin/python
# Programmer : zocean
# Date: 
# Last-modified: 01 Nov 2014 04:57:44 PM

import os,sys
import tempfile
import tabix

def logging(text):
    print >>sys.stderr, text

def warning(text):
    print >>sys.stderr, "Warning: " + text

def error(text,code=None):
    if code is not None:
        print >>sys.stderr, "Error Num:%d %s" % (code,text)
    else:
        print >>sys.stderr, "Error: %s" % (text)

def ReadFromFile(fin_name):
    try:
        fin = open(fin_name,"r")
    except IOError:
        error("Can't open file:%s to read." % (fin_name))
        exit(1)
    return fin

def WriteToFile(fout_name,is_cover=True):
    if fout_name == 'stdout' or not fout_name:
        fout = sys.stdout
        warning("use stdout to report result")
    else:
        if CheckFileExist(fout_name):
            if is_cover:
                warning("file:%s already exist. Old file will be covered by new file." % (fout_name))
                try:
                    fout = open(fout_name,"w")
                except IOError:
                    warning("Can't open file:%s to write. Use stdout instead." % (fout_name))
                    fout = sys.stdout
            else:
                error("file:%s already exist" % (fout_name))
                exit(1)
        else:
            try:
                fout = open(fout_name,"w")
            except IOError:
                warning("Can't open file:%s to write. Use stdout instead." % (fout_name))
    return fout
                    
def CheckFolderExist(folder,is_cover=True):
    """Check the existence of folder, it not exist will create it"""
    if not os.path.exists(folder):
        os.mkdir(folder)
    else:
        if is_cover:
            warning("output folder %s already exists, may overwrite the content inside it." % (folder))
        else:
            error("output folder %s already exists, will exit." % (folder))
            exit(1)

def CheckFileExist(filename,is_cover=True):
    """Check the existence of folder"""
    if not os.path.isfile(filename):
        pass
    else:
        if is_cover:
            warning("output file %s already exsits, will overwrite the content inside it." % (filename))
        else:
            error("output file %s already exists, will exit." % (filename))
            exit(1)

def LoadGenome(filename):
    genome_table = {}
    for line in ReadFromFile(filename):
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        genome_table[row[0]] = int(row[1])
    return genome_table

def SortGenome(genome_table):
    """sort genome dict in number order, return sorted chrom name"""
    chrlist = genome_table.keys()
    sort_table = {}
    for chrom in chrlist:
        if chrom[0:3].lower() == 'chr':
            try:
                sort_table[int(chrom[3:])] = chrom
            except ValueError:
                sort_table[chrom[3:]] = chrom
        else:
            sort_table[chrom] = chrom
    sorted_chrom = sorted(sort_table.keys())
    return [sort_table[name] for name in sorted_chrom]

def SortBed(filename):
    """sort bed file by chrom then by start"""
    tmpout = tempfile.NamedTemporaryFile()
    genome = {}
    for line in ReadFromFile(filename):
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        chrom = row[0]
        genome[chrom] = 0
    sorted_chrom = SortGenome(genome)
    table = {}
    for line in ReadFromFile(filename):
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        chrom = row[0]
        start = int(row[1])
        stop = int(row[2])
        other = "\t".join(item for item in row[3:])
        try:
            table[chrom].append([chrom,start,stop,other])
        except KeyError:
            table[chrom] = [[chrom,start,stop,other]]
    sorted_table = {}
    for chrom in table.keys():
        sorted_table[chrom] = sorted(table[chrom],key=lambda bed:bed[1])
    for chrom in sorted_chrom:
        for bed in sorted_table[chrom]:
            print >>tmpout, bed[0] + '\t' + str(bed[1]) + '\t' + str(bed[2]) + '\t' + bed[3]
    tmpout.flush()
    os.system("cp %s %s" % (tmpout.name,filename))

def CreateTabix(filename,zippedname,mode):
    os.system("bgzip -c %s >%s" % (filename,zippedname))
    os.system("tabix -p %s %s" % (mode,zippedname))
