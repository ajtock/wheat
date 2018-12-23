#!/usr/bin/env python
import csv
import sys

# ajt: added line for generating Unix (rather than Windows) text file
csv.register_dialect("unixEOL", delimiter="\t", lineterminator="\n")

f = csv.reader(sys.stdin, dialect="unixEOL")
of = csv.writer(sys.stdout, dialect="unixEOL")
last_read = None
for line in f :
    #take care of the header
    if(line[0][0] == "@") :
        of.writerow(line)
        continue

    if(last_read == None) :
        last_read = line
    else :
        if(last_read[0] == line[0]) :
            of.writerow(last_read)
            of.writerow(line)
            last_read = None
        else :
            last_read = line
