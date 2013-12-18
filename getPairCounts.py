import sys
import MySQLdb
from datetime import datetime

# Binary search to match exon start-end coordinates & read's start-end coordinates
def recurse(a,n1,n2):
    if n2-n1==1 or n2==n1:
        return n1
    else:
        if a<text[int((n2-n1)/2)+n1][2]:
            return recurse(a,n1,int((n2-n1)/2)+n1)
        else:
            return recurse(a,int((n2-n1)/2)+n1,n2)
        
def helpme():
    print '\n--------------------------------------------------------------'
    print 'Usage: python getPairCount.py <inputfile> <outputfile> <mysql-hostname> <mysqluser> <mysqlpass> <mysqlDBname>'
    print 'This program will count the mapped exons'
    print 'First you must run the preprocessing bash script (for pair sorting/coupling, filtering MAPQ>=30, removing reads having chr_of_pair1 != chr_of_pair2, and separating BAM reads into chromosomes)'
    print 'Example: python getPairCount.py human1 outputCount.txt hostname.meb.ki.se user pass dbname'
    print '--------------------------------------------------------------\n'

if len(sys.argv)==2 and sys.argv[1]=='--help':
	helpme()
else:
	f=open(sys.argv[2], "a")
	fileinput=sys.argv[1]

	# Connect to DB
    db = MySQLdb.connect(sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    cursor=db.cursor()
    
    for i in xrange(1,23):
        exon={}
        print >> sys.stderr, "Processing Chromosome #"+str(i)
        print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        cursor.execute("SELECT SQL_CACHE region_id,exon_name,start,cigar,end FROM human where Chromosome='\"chr"+str(i)+"\"' order by start")
        text=cursor.fetchall()
        header=''
        for line in file(fileinput+".chr"+str(i), 'rU'):
            pos=int(line.split('\t')[3])
            idx=recurse(pos,0,len(text))
            # find possibility of other exon coordinates
            if idx<len(text)-1:
                newidx=idx
		# to evade different seqnames but same coordinates
                while text[newidx+1][2]==text[newidx][2]:
                    newidx+=1
                # to allow 5bp window
		if int(line.split('\t')[3])+5>text[newidx+1][2] and int(line.split('\t')[3])+5<text[newidx][2] and int(line.split('\t')[3])+5<text[newidx][4]:
                    idx=newidx+1
                    newidx=-1
                # to skip reads that's outside the 5bp window
		if int(line.split('\t')[3])+5<text[newidx+1][2] and int(line.split('\t')[3])+5>text[newidx][4]:
                    newidx=-1
                    continue
            # IF it's the new pair, treat as first pair, count, save for next loop
            if header!=line.split('\t')[0].split('.')[1]:
                header=line.split('\t')[0].split('.')[1]
                exon_txt=""
                gene_txt=""
                oldpos=0
		# get exon & region name if the cigar match
                if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0]):
                    exon_txt=text[idx][1]
                    gene_txt=text[idx][0]
                    oldpos=pos
            # ELSE IF the first pair has been counted
            else:
                if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5] and gene_txt==text[idx][0]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0] and gene_txt==text[idx][0]):
                    if oldpos<pos:
                        try:
                            exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]+=1
                        except KeyError:
                            exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]=1
                        # for debugging, uncomment these
			# if exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)=='ex_9.ex_9__17850':
                        #    f.write(line)
                    else:
                        try:
                            exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]+=1
                        except KeyError:
                            exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]=1
                        # for debugging, uncomment these
                        # if text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)=='ex_9.ex_9__17850':
                        #    f.write(line)
        # print exon-pair name and the counts
	for k in exon.keys():
            f.write(str(k)+"\t"+str(exon[k])+"\n")

    # copy paste from above for chrX and chrY, because after 22 is 23, not X and Y :P
    exon={}
    print >> sys.stderr, "Processing Chromosome X"
    print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    cursor.execute("SELECT SQL_CACHE region_id,exon_name,start,cigar,end FROM human where Chromosome='\"chrX\"' order by start")
    text=cursor.fetchall()
    header=''
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
            exon_txt=""
            gene_txt=""
            oldpos=0
            if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0]):
                exon_txt=text[idx][1]
                gene_txt=text[idx][0]
                oldpos=pos
        else:
            if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5] and gene_txt==text[idx][0]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0] and gene_txt==text[idx][0]):
                if oldpos<pos:
                    try:
                        exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]+=1
                    except KeyError:
                        exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]=1
                    # if exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)=='ex_3.ex_3__19479':
                    #    f.write(line)
                else:
                    try:
                        exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]+=1
                    except KeyError:
                        exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]=1
                    # if text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)=='ex_3.ex_3__19479':
                    #    f.write(line)
    for k in exon.keys():
        f.write(str(k)+"\t"+str(exon[k])+"\n")

    exon={}
    print >> sys.stderr, "Processing Chromosome Y"
    print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    cursor.execute("SELECT SQL_CACHE region_id,exon_name,start,cigar,end FROM human where Chromosome='\"chrY\"' order by start")
    text=cursor.fetchall()
    header=''
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
            exon_txt=""
            gene_txt=""
            oldpos=0
            if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0]):
                exon_txt=text[idx][1]
                gene_txt=text[idx][0]
                oldpos=pos
        else:
            if ('N' not in text[idx][3] and 'N' not in line.split('\t')[5] and gene_txt==text[idx][0]) or ('N' in text[idx][3] and 'N' in line.split('\t')[5] and text[idx][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[1].split('N')[0] and gene_txt==text[idx][0]):
                if oldpos<pos:
                    try:
                        exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]+=1
                    except KeyError:
                        exon[exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)]=1
                    #if exon_txt.replace("\"", "")+'.'+text[idx][1].replace("\"", "")+'__'+str(gene_txt)=='ex_3.ex_3__19479':
                    #    f.write(line)
                else:
                    try:
                        exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]+=1
                    except KeyError:
                        exon[text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)]=1
                    #if text[idx][1].replace("\"", "")+'.'+exon_txt.replace("\"", "")+'__'+str(gene_txt)=='ex_3.ex_3__19479':
                    #    f.write(line)
    for k in exon.keys():
        f.write(str(k)+"\t"+str(exon[k])+"\n")
print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
f.close()
