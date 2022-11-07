import sys
import re
import os.path
import time
import datetime
import csv
import shutil

INFILE      = ""
FAILFOLDER  = ""
OUTFILE     = ""
CORRECTFILE = ""
dict        = {}

def exitparser():
    # try to print the contents to an output file
    global dict, INFILE, FAILFOLDER
    counter = 0
    while True:
        counter += 1
        failname = "{}/{}_{}_{}_{}.out".format(FAILFOLDER,dict["model"],dict["alg"],dict["workers"],counter)
        if not (os.path.isfile(failname)):
            break
    shutil.copy2(INFILE, failname)
   

def checkfile(file):
    if not (os.path.isfile(file)):
        print("ERROR: cannot find file: {}".format(file))
        exitparser()


def parsevar(varname, line, regex):
    pattern = re.compile(regex)
    m = pattern.search(line)
    if (m):
        global dict
        if (dict.get(varname)):
            print("ERROR: multiple matches for {}".format(varname))
            exitparser()
        else:
            dict[varname] = m.group(1)
        if varname == "assert": print(line) #show assert info

 #COLUMN="model,N,M,workers,alg,init-mstime,insert-num,remove-num,mstime,V*,V+,S,tag,assert,date"
def parseline(line):
    #assert 
    parsevar("assert", line, r"ASSERT([\S]+)")
    parsevar("assert", line, r"([\S]+)Assertion[\S]+")

    parsevar("N",   line, r"N = ([\S]+), M = [\S]+")
    parsevar("M",   line, r"N = [\S]+, M = ([\S]+)")

    parsevar("init-mstime",      line, r"initialization costs\(ms\): ([\S]+)")
    parsevar("insert-num",      line, r"# of edges to insert: ([\S]+)")
    parsevar("remove-num",      line, r"# of edges to delete: ([\S]+)")
    parsevar("mstime",      line, r"core IorR costs\(ms\): ([\S]+)")

    parsevar("V*",      line, r"V\* size: ([\S]+)")
    parsevar("V+",      line, r"V\+ size: ([\S]+)")
    parsevar("S",      line, r"S size: ([\S]+)")
    parsevar("batch-repeat", line, r"batch insertion repeat: ([\S]+)")
    
    parsevar("snum", line, r"edges size: ([\S]+)")
    
    parsevar("tag",      line, r"our adjust tag: ([\S]+)")

def parsefile(file):
    f = open(file, 'r')
    for line in f:
        parseline(line)
    f.close()


def trytoprint(varname):
    global dict
    if (dict.get(varname)):
        return dict.get(varname)
    else:
        return "-1";


def printtofile(outfile):
    # First line of OUTFILE should contain comma-separated info on column names
    global dict
    f = open(outfile, 'r+')
    s = f.readline().strip()
    names = s.split(",")
    output  = ""
    for name in names:
        output += trytoprint(name) + ","
    output = output[:-1] # remove last ","
    f.read() # go to the last line
    f.write(output+"\n") # write the new line
    f.close()


def printtostdout():
    # find longest key name (for formatting)
    global dict
    maxlen = 1
    for varname in dict:
        maxlen = max(maxlen, len(varname))
    for varname in dict:
        print(varname + ":" +  " " * (maxlen+1-len(varname)) + dict[varname])


def addtimestamp():
    global dict
    ts = time.time()
    dict["date"] = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d.%H:%M:%S')


def main():
    global dict, INFILE, OUTFILE, FAILFOLDER, CORRECTFILE
    N_ARG = len(sys.argv)
    if (N_ARG == 6):
        dict["model"]   = str(sys.argv[1])
        dict["alg"]     = str(sys.argv[2])
        dict["workers"] = str(sys.argv[3])
        FAILFOLDER      = str(sys.argv[4])
        INFILE          = str(sys.argv[5])
        checkfile(INFILE)
        parsefile(INFILE)
        addtimestamp()
        printtostdout()
    elif (N_ARG == 7):
        dict["model"]   = str(sys.argv[1])
        dict["alg"]     = str(sys.argv[2])
        dict["workers"] = str(sys.argv[3])
        FAILFOLDER      = str(sys.argv[4])
        INFILE          = str(sys.argv[5])
        OUTFILE         = str(sys.argv[6])
        checkfile(INFILE)
        parsefile(INFILE)
        addtimestamp()
        #print(dict)
        printtofile(OUTFILE)
    else:
        print("ERROR: invalid command")
        print("Usage:")
        print(" - python parse_output.py  MODEL  ALG  WORKERS  FAILFOLDER  INFILE           # writes all output to the stdout")
        print(" - python parse_output.py  MODEL  ALG  WORKERS  FAILFOLDER  INFILE  OUTFILE  # appends the data to the OUTFILE in the same format used")


if __name__ == "__main__":
    main()

