import sys
#import MySQLdb
import sqlite3
from datetime import datetime

# Binary search
def recurse(a,n1,n2):
    if n2-n1==1 or n2==n1:
        return n1
    else:
        if a<text[int((n2-n1)/2)+n1][2]:
            return recurse(a,n1,int((n2-n1)/2)+n1)
        else:
            return recurse(a,int((n2-n1)/2)+n1,n2)

def helpme():
    print ''
    print 'Program: getPairCount.py (Tools for counting the mapped exons)'
    print ''
    print 'Usage: python getPairCount.py <inputprefix> <outputfile> <sqlitefile>'
    print ''
    print 'Command: inputprefix\tIf the input files are human1.chr1, ..., human1.chrY, then inputprefix is human1'
    print '         outputfile\tthe count will be outputted to this file'
    print '         sqlitefile\tThe sqlite of txdb'
    print ''
    print 'Example: python getPairCount.py human1 outputCount.txt hg18.sqlite'
    print ''

if len(sys.argv)<3 or sys.argv[1]=='--help':
        helpme()
else:
    #f=open("deleteme21.txt", "a")
    f=open(sys.argv[2], "a")
    #fileinput="reg"
    fileinput=sys.argv[1]

    # Connect to DB
    db = sqlite3.connect(sys.argv[3])
    #db = MySQLdb.connect("D0069MEB.meb.ki.se", "everyone", "password", "dhany")
    cursor=db.cursor()
    alone=0
    multiple=0
    for i in xrange(1,24):
        exon={}
        print >> sys.stderr, "Processing Chromosome #"+str(i)
        print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        cursor.execute("SELECT region_id,exon_name,start,cigar,end FROM hg18 where Chromosome='\"chr"+str(i)+"\"' order by start")
        text=cursor.fetchall()
        header=''
        flag_counter=0
        for line in file(fileinput+".chr"+str(i), 'rU'):
            pos=int(line.split('\t')[3])
            idx=recurse(pos,0,len(text))
            if idx<len(text)-1:
                newidx=idx
                while text[newidx+1][2]==text[newidx][2]:
                    newidx+=1
                if int(line.split('\t')[3])+5>text[newidx+1][2] and int(line.split('\t')[3])+5<text[newidx][2] and int(line.split('\t')[3])+5<text[newidx][4]