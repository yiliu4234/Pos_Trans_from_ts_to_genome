#!/yjdata/hxzk/utilities/miniconda3/bin/python3

import sys
import importlib
import os
import getopt
import re
import numpy as np

bin_path = os.path.dirname(__file__)
module_path = sys.path.append("%s/../module" % bin_path)
sys.path.append(bin_path) 
sys.path.append(module_path)
import read_arg


def USAGE(script):
    helpMess = '''
     ProgramName:\t%s
         Version:\tV1.0
         Contact:\tLiaoyw yiliu4234@163.com
    Program Date:\t2024.4.25
          Modify:\t-
     Description:\tThis is a software transforming enst position to genome position 
           Usage:\tpython3 %s -g <gtf file> -b <bed file>  -o <outputfile>
                 \t~/script/pipeline_v3.0/bin/TS_to_genome.py -g /yjdata/database/star/gtf/human/Homo_sapiens.GRCh38.99.gtf -b RNAIP_RNA_TestVsInput_peaks.xls -o peak_genome.txt
         Options:
            -h --help              Show this help message.
            -g --gtf   gtf file
            -b --bed  bed file
            -o --output outputfile
    ''' %(script,script)
    return helpMess

def index_of_min(lst):
    '''
    index of min value
    '''
    return min(range(len(lst)), key=lst.__getitem__)

def read_args(argv):
    '''
    read cmmand arg    
    '''
    argDict = {}
    argDict["status"] = False
    try:
        opts, args = getopt.getopt(argv,"hg:b:o:",["help","gtf_path=","bed_path=","output="])
    except getopt.GetoptError:
        print(USAGE(sys.argv[0]))
        sys.exit(2)
    if len(opts) == 0:
        print(USAGE(sys.argv[0]))
        sys.exit()
    for opt, arg in opts:
        if opt in ('-h','--help'):
            print(USAGE(sys.argv[0]))
            sys.exit()
        elif opt in ("-g", "--gtf"):
            argDict["gtf_path"] = os.path.abspath(arg)
        elif opt in ("-b", "--bed"):
            argDict["bed_path"] = os.path.abspath(arg)
        elif opt in ("-o", "--output"):
            argDict["output"] = os.path.abspath(arg)
    return(argDict)

class TS(object):
    '''
    Transcript ID class, include all peak information of in that TS
    '''

    def __init__(self, TS_name):
        self.name = TS_name
        self.Peaks = []
        self.chr = "" 
        self.exon = []
        self.exon_len = []
    def add_peak(self, peak):
        self.Peaks.append(peak)
    def add_info(self, chr, start, end, strand, symbol):
        self.chr = chr
        self.strand = strand
        self.symbol = symbol
        self.exon.append((int(start), int(end)))
    def get_pos(self):
        self.exon_len = [i[1] - i[0] +1 for i in self.exon]
        self.exon_cumsum = np.cumsum(self.exon_len)

        for peak in self.Peaks:
            peak.chr = self.chr
            peak.symbol = self.symbol
            peak.strand = self.strand
            for i in range(len(self.exon_cumsum)):
                if self.exon_cumsum[i] >= peak.ts_start:
                    if self.strand == "+":
                        peak.start = self.exon[i][1] - (self.exon_cumsum[i] - peak.ts_start)
                        exon_num_start = i
                    else:
                        peak.end = self.exon[i][0] + (self.exon_cumsum[i] - peak.ts_start)
                        exon_num_start = i
                    break 
            for i in range(exon_num_start, len(self.exon_cumsum)):
                if self.exon_cumsum[i] >= peak.ts_end:
                    if self.strand == "+":
                        peak.end = self.exon[i][1] - (self.exon_cumsum[i] - peak.ts_end)
                    else:
                        peak.start = self.exon[i][0] + (self.exon_cumsum[i] - peak.ts_end)
                    break 
            peak.position = "Trans exon" if i != exon_num_start else "Same exon"
    def print_pos(self,f):
        for peak in self.Peaks:
            peak.print_pos(f)
         
class peak(object):
    '''
    Peak class, include all star and end  of peek in transcript, and corresponding chromsome, star and end in genome 
    '''
    def __init__(self, ts_start, ts_end, TS_name, info):
        self.ts_start = int(ts_start)
        self.ts_end = int(ts_end)
        self.TS_name = TS_name
        self.symbol = ""
        self.chr = "" 
        self.start = 1
        self.end = 1  
        self.info = "\t".join(info)
        self.position = "" # Transexon or Same exon 
    def print_pos(self,f):
        '''
        print position of of genome
        '''
        f.writelines("\t".join([self.chr, str(self.start), str(self.end), self.position, self.TS_name, self.symbol, self.strand, self.info]) + "\n")

def write_out_bed(TSs, output, raw_header):
    '''
    Write out the bed result of al peaks postion in genome
    '''
    raw_header = raw_header.split("\t")
    header = "\t".join(["chr", "start", "end", "position", "ENST", "symbol", "strand"]) + "\t" + "\t".join(raw_header[3:len(raw_header)])
    with open(output, 'w+') as f:
         f.writelines(header+"\n")
         for ts_name in TSs:
             TSs[ts_name].print_pos(f)


if __name__ == '__main__':
    # read args
    argv = sys.argv[1:]
    argDict = read_args(argv)    
    gtf = argDict["gtf_path"]
    bed = argDict["bed_path"]
    output = argDict["output"]

    TSs = {}
    raw_header = ""
    with open(bed, 'r') as f:
        for eachline in f.readlines():
            if eachline.startswith("#"):
                continue
            if eachline.startswith("chr"):
                raw_header = eachline.strip()
                continue
            if eachline.strip()=="":
                continue
            info = eachline.strip().split("\t")
            ts_name = info[0]
            start = info[1]
            end = info[2]
            otherinfo = info[3:len(info)]
            if not ts_name in TSs:
                TSs.update({ts_name:TS(ts_name)})
            TSs[ts_name].add_peak(peak(start, end, ts_name, otherinfo))

    with open(gtf, 'r') as f:
        for eachline in f.readlines():
            if eachline.startswith("#"):
                continue
            info = eachline.strip().split("\t")
            if info[2] != "exon":
                continue
            chr = info[0]
            start = info[3]
            end = info[4]
            strand = info[6] 
            otherinfo = info[8]
            ts_name = re.findall(r"transcript_id \"(.+?)\"",otherinfo)[0]
            symbol = re.findall(r"gene_name \"(.+?)\"",otherinfo)[0]
            if ts_name not in TSs:
                continue
            TSs[ts_name].add_info(chr, start, end, strand,symbol)
    for ts_name in TSs:
       TSs[ts_name].get_pos()
    write_out_bed(TSs, output, raw_header)  

