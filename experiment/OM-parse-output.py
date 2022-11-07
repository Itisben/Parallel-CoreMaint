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
    # counter = 0
    # while True:
    #     counter += 1
    #     failname = "{}/{}_{}_{}_{}.out".format(FAILFOLDER,dict["model"],dict["alg"],dict["workers"],counter)
    #     if not (os.path.isfile(failname)):
    #         break
    # shutil.copy2(INFILE, failname)
   

def checkfile(file):
    if not (os.path.isfile(file)):
        print("ERROR: cannot find file: {}".format(file))
        exitparser()


def parsevar(varname, line, regex):
    pattern = re.compile(regex)
    #m = pattern.search(line)
    m = pattern.match(line)
    if (m):
        global dict
        if (dict.get(varname)):
            print("ERROR: multiple matches for {}".format(varname))
            exitparser()
        else:
            dict[varname] = m.group(1)
        if varname == "assert": print(line) #show assert info

def parseline(line):
    #parsevar("I",           line, r"insert size ([\S]+)")
    parsevar("l",           line, r"Lock Type: ([\S]+)")
    parsevar("w",           line, r"set number of worker: ([\S]+)")


    parsevar("Itime",      line, r"TestInsert costs\(ms\): ([\S]+)")
    parsevar("Otime",      line, r"TestOrder costs\(ms\): ([\S]+)")
    parsevar("Dtime",      line, r"TestDelete costs\(ms\): ([\S]+)")
    parsevar("Mtime",      line, r"TestMixed costs\(ms\): ([\S]+)")
    
    parsevar("tag#",      line, r"tag #: ([\S]+)")
    parsevar("subtag#",   line, r"subtag #: ([\S]+)")
    parsevar("orepeat#",      line, r"order repeat #: ([\S]+)")
    parsevar("relabel#",      line, r"relabel #: ([\S]+)")

    parsevar("maxtag#",      line, r"rebalance max tag #: ([\S]+)")
    parsevar("morepeat#",      line, r"Mixed order repeat #: ([\S]+)")

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
        return "-11";


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
    if (N_ARG == 7):
        dict["I"] = str(sys.argv[1])
        dict["T"] = str(sys.argv[2])
        dict["t"] = str(sys.argv[3])
        FAILFOLDER      = str(sys.argv[4])
        INFILE          = str(sys.argv[5])
        OUTFILE         = str(sys.argv[6])
        checkfile(INFILE)
        parsefile(INFILE)
        addtimestamp()
        printtofile(OUTFILE)
    else:
        print("ERROR: invalid command")
        print("Usage:")
        print(" - python parse_output.py  MODEL  ALG  WORKERS  FAILFOLDER  INFILE           # writes all output to the stdout")
        print(" - python parse_output.py  MODEL  ALG  WORKERS  FAILFOLDER  INFILE  OUTFILE  # appends the data to the OUTFILE in the same format used")


if __name__ == "__main__":
    main()

