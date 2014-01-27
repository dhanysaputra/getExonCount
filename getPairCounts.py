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
    print 'Command: inputprefix\tIf the input files are tumor.chr1, ..., tumor.chrY, then inputprefix is tumor'
    print '         outputfile\tthe count will be outputted to this file'
    print '         sqlitefile\tThe sqlite of txdb'
    print ''
    print 'Example: python getPairCount.py tumor outputCount.txt hg18.sqlite'
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
                if int(line.split('\t')[3])+5>text[newidx+1][2] and int(line.split('\t')[3])+5<text[newidx][2] and int(line.split('\t')[3])+5<text[newidx][4]:
                    idx=newidx+1
                    newidx=-1
                if int(line.split('\t')[3])+5<text[newidx+1][2] and int(line.split('\t')[3])+5>text[newidx][4]:
                    newidx=-1
                    continue
            if header!=line.split('\t')[0].split('.')[1]:
                header=line.split('\t')[0].split('.')[1]
                # doesnt have couple
                exon_txt=""
                gene_txt=""
                oldpos=0
                if flag_counter==1:
                    alone+=1
                # a pair consist of more than 2 reads, or multiple aligned
                elif flag_counter>2:
                    multiple+=1
                if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0]):
                    exon_txt=text[idx][1]
                    gene_txt=text[idx][0]
                    oldpos=pos
                flag_counter=1
            else:
                if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5] and gene_txt==text[idx][0]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0] and gene_txt==text[idx][0]):
                    if oldpos<pos:
                        try:
                            exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]+=1
                        except KeyError:
                            exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]=1
                        #if exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)=='ex_5-ex_6.ex_6__15255':
                        #    f.write(line)
                    else:
                        try:
                            exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]+=1
                        except KeyError:
                            exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]=1
                        #if text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)=='ex_5-ex_6.ex_6__15255':
                        #    f.write(line)
        for k in exon.keys():
            f.write(str(k)+"\t"+str(exon[k])+"\n")
    exon={}
    print >> sys.stderr, "Processing Chromosome X"
    print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    cursor.execute("SELECT region_id,exon_name,start,cigar,end FROM hg18 where Chromosome='X' order by start")
    text=cursor.fetchall()
    header=''
    flag_counter=0
    for line in file(fileinput+".chrX", 'rU'):
        pos=int(line.split('\t')[3])
        idx=recurse(pos,0,len(text))
        if idx<len(text)-1:
            newidx=idx
            while text[newidx+1][2]==text[newidx][2]:
                newidx+=1
            if int(line.split('\t')[3])+5>text[newidx+1][2] and int(line.split('\t')[3])+5<text[newidx][2] and int(line.split('\t')[3])+5<text[newidx][4]:
                idx=newidx+1
                newidx=-1
            if int(line.split('\t')[3])+5<text[newidx+1][2] and int(line.split('\t')[3])+5>text[newidx][4]:
                newidx=-1
                continue
        if header!=line.split('\t')[0].split('.')[1]:
            header=line.split('\t')[0].split('.')[1]
            # doesnt have couple
            exon_txt=""
            gene_txt=""
            oldpos=0
            if flag_counter==1:
                alone+=1
            # a pair consist of more than 2 reads, or multiple aligned
            elif flag_counter>2:
                multiple+=1
            if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0]):
                exon_txt=text[idx][1]
                gene_txt=text[idx][0]
                oldpos=pos
            flag_counter=1
        else:
            if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5] and gene_txt==text[idx][0]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0] and gene_txt==text[idx][0]):
                if oldpos<pos:
                    try:
                        exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]+=1
                    except KeyError:
                        exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]=1
                    #if exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)=='ex_5-ex_6.ex_6__15255':
                    #    f.write(line)
                else:
                    try:
                        exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]+=1
                    except KeyError:
                        exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]=1
                    #if text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)=='ex_5-ex_6.ex_6__15255':
                    #    f.write(line)
    for k in exon.keys():
        f.write(str(k)+"\t"+str(exon[k])+"\n")
    exon={}
    print >> sys.stderr, "Processing Chromosome Y"
    print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    cursor.execute("SELECT region_id,exon_name,start,cigar,end FROM hg18 where Chromosome='Y' order by start")
    text=cursor.fetchall()
    header=''
    flag_counter=0
    for line in file(fileinput+".chrY", 'rU'):
        pos=int(line.split('\t')[3])
        idx=recurse(pos,0,len(text))
        if idx<len(text)-1:
            newidx=idx
            while text[newidx+1][2]==text[newidx][2]:
                newidx+=1
            if int(line.split('\t')[3])+5>text[newidx+1][2] and int(line.split('\t')[3])+5<text[newidx][2] and int(line.split('\t')[3])+5<text[newidx][4]:
                idx=newidx+1
                newidx=-1
            if int(line.split('\t')[3])+5<text[newidx+1][2] and int(line.split('\t')[3])+5>text[newidx][4]:
                newidx=-1
                continue
        if header!=line.split('\t')[0].split('.')[1]:
            header=line.split('\t')[0].split('.')[1]
            # doesnt have couple
            exon_txt=""
            gene_txt=""
            oldpos=0
            if flag_counter==1:
                alone+=1
            # a pair consist of more than 2 reads, or multiple aligned
            elif flag_counter>2:
                multiple+=1
            if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0]):
                exon_txt=text[idx][1]
                gene_txt=text[idx][0]
                oldpos=pos
            flag_counter=1
        else:
            if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5] and gene_txt==text[idx][0]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0] and gene_txt==text[idx][0]):
                if oldpos<pos:
                    try:
                        exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]+=1
                    except KeyError:
                        exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]=1
                    #if exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)=='ex_5-ex_6.ex_6__15255':
                    #    f.write(line)
                else:
                    try:
                        exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]+=1
                    except KeyError:
                        exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]=1
                    #if text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)=='ex_5-ex_6.ex_6__15255':
                    #    f.write(line)
    for k in exon.keys():
        f.write(str(k)+"\t"+str(exon[k])+"\n")
        print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        conn.close()
        f.close()