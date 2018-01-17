### #! path to your python
import csv
total = 0
import pandas as pd
#Raw inputs
ftw = raw_input("What file to use: ")
taj = raw_input("Search for higher TAJIMA D values than: ")
theta = raw_input("Seach for higher THETA values than: ")

files = ftw + '.txt.thetasWindow.gz.pestPG'
ofiles = ftw + '.txt'
f = open(files, 'r')
reader = csv.reader(f, delimiter = '\t')
headers = reader.next()
# lines = r.readlines()
of = open(ofiles, 'w')

#To use headers in file
column = {}
converters = [str.strip] + [float] * (len(headers)-1)
for h in headers:
    column[h] = []
for row in reader:
    for h, v, conv in zip(headers, row, converters):
        column[h].append(conv(v))

#Empty lists for line numbers of given true values for Tajimas D and LinesT
#This is for column reading. I think..

linesD = []
linesT = []
for i,h in enumerate(column['Tajima'],0):
    if float(h) >= float(taj):
        num = i
        linesD.append(num)
for i,h in enumerate(column['tW'],0):
    if float(h) >= float(theta):
        num = i
        linesT.append(num)
# f.close()
#Find if both are true
sD = set(linesD)
sT = set(linesT)
tl = sD.intersection(sT)

#Need to add 1 to each element in intersection because the above function did it without the header
tll = [x+1 for x in tl]
f.close()

# This is to print out the actual lines
r = open(files)
g = "#(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)	Chr	WinCenter	tW	tP	tF	tH	tL	Tajima	fuf	fud	fayh	zeng	nSites \n"
with open(files) as r:
    n = 0
    of.write(g)
    for line in r:
        if n in sorted(tll):
            print n
            of.write(line)
        n+=1
r.close()
