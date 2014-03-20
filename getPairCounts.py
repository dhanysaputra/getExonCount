import sys
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
	f=open(sys.argv[2], "a")
	fileinput=sys.argv[1]

    # Connect to DB
	db = sqlite3.connect(sys.argv[3])
	cursor=db.cursor()
	oldpos=-1
	for i in xrange(1,23):
		exon={}
		print >> sys.stderr, "Processing Chromosome #"+str(i)
		print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
		cursor.execute("SELECT region_id,exon_name,start,cigar,end FROM hg19 where Chromosome='\"chr"+str(i)+"\"' order by start")
		text=cursor.fetchall()
		flag_finish=1
		temptext=''
		for line in file(fileinput+".chr"+str(i), 'rU'):
			if flag_finish==1:
				temptext=line
				flag_finish=0
			else:
				# Skip counting if there's many N's in cigars
				if len(line.split('\t')[5].split('N'))<2 and len(temptext.split('\t')[5].split('N'))<2:
					# Fast search algo for finding exon (Binary Tree)
					pos1=int(temptext.split('\t')[3])
					pos2=int(line.split('\t')[3])
					idx1=recurse(pos1,0,len(text))
					idx2=recurse(pos2,0,len(text))
	
					# Treating exceptional cases: INDELS in cigars!
					indels1=0
					indels2=0
					if 'I' in temptext.split('\t')[5]:
						for i in xrange(0,len(temptext.split('\t')[5].split('I'))-1):
							indels1+=int(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N')[len(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N'))-1])
					if 'I' in line.split('\t')[5]:
						for i in xrange(0,len(line.split('\t')[5].split('I'))-1):
							indels2+=int(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N')[len(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N'))-1])
					if 'D' in temptext.split('\t')[5]:
						for i in xrange(0,len(temptext.split('\t')[5].split('D'))-1):
							indels1-=int(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N')[len(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N'))-1])
					if 'D' in line.split('\t')[5]:
						for i in xrange(0,len(line.split('\t')[5].split('D'))-1):
							indels2-=int(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N')[len(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N'))-1])
	
					# Another exceptional cases: this coordinate is in exon and junction ranges, or weird rows (same exons/junctions, same cigard)!
					if idx1<len(text)-1:
						newidx=idx1
						while text[newidx+1][2]==text[newidx][2]:
							if ('N' in text[idx1][3] and 'N' in temptext.split('\t')[5]):
								if (text[idx1][3].split('M')[1].split('N')[0]==temptext.split('\t')[5].split('M')[len(temptext.split('\t')[5].split('N')[0].split('M'))].split('N')[0]):
									newidx+=1
									break
							newidx+=1
						if int(temptext.split('\t')[3])+0>=text[newidx+1][2] and int(temptext.split('\t')[3])+49-indels1<=text[newidx+1][4]:
							idx1=newidx+1
						newidx=-1
					if idx2<len(text)-1:
						newidx=idx2
						while text[newidx+1][2]==text[newidx][2]:
							if ('N' in text[idx2][3] and 'N' in line.split('\t')[5]):
								if (text[idx2][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[len(line.split('\t')[5].split('N')[0].split('M'))].split('N')[0]):
									newidx+=1
									break
							newidx+=1
						if int(line.split('\t')[3])+0>=text[newidx+1][2] and int(line.split('\t')[3])+49-indels2<=text[newidx+1][4]:
							idx2=newidx+1
						newidx=-1
					if (('N' not in text[idx1][3] and 'N' in temptext.split('\t')[5]) or ('N' in text[idx1][3] and 'N' not in temptext.split('\t')[5])):
						newidx=idx1-1
						while (('N' not in text[newidx][3] and 'N' in temptext.split('\t')[5]) or ('N' in text[newidx][3] and 'N' not in temptext.split('\t')[5])):
							newidx=newidx-1
						if int(temptext.split('\t')[3])+0>=text[newidx][2] and int(temptext.split('\t')[3])+49-indels1<=text[newidx][4]:
							idx1=newidx
					if (('N' not in text[idx2][3] and 'N' in line.split('\t')[5]) or ('N' in text[idx2][3] and 'N' not in line.split('\t')[5])):
						newidx=idx2-1
						while (('N' not in text[newidx][3] and 'N' in line.split('\t')[5]) or ('N' in text[newidx][3] and 'N' not in line.split('\t')[5])):
							newidx=newidx-1
						if int(line.split('\t')[3])+0>=text[newidx][2] and int(line.split('\t')[3])+49-indels2<=text[newidx][4]:
							idx2=newidx
	
					# Increment count if both reads fulfil the condition
					if ('N' not in text[idx1][3] and 'N' not in temptext.split('\t')[5] and int(temptext.split('\t')[3])+0>=text[idx1][2] and int(temptext.split('\t')[3])+49-indels1<=text[idx1][4]) or ('N' in text[idx1][3] and 'N' in temptext.split('\t')[5] and text[idx1][3].split('M')[1].split('N')[0]==temptext.split('\t')[5].split('M')[len(temptext.split('\t')[5].split('N')[0].split('M'))].split('N')[0] and int(text[idx1][3].split('M')[0].split('\"')[1])+int(text[idx1][2])==pos1+int(temptext.split('\t')[5].split('M')[0]) and int(temptext.split('\t')[3])+0>=text[idx1][2] and int(temptext.split('\t')[3])+49-indels1<=text[idx1][4]):
						if ('N' not in text[idx2][3] and 'N' not in line.split('\t')[5] and text[idx1][0]==text[idx2][0] and int(line.split('\t')[3])+0>=text[idx2][2] and int(line.split('\t')[3])+49-indels2<=text[idx2][4]) or ('N' in text[idx2][3] and 'N' in line.split('\t')[5] and text[idx2][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[len(line.split('\t')[5].split('N')[0].split('M'))].split('N')[0] and int(text[idx2][3].split('M')[0].split('\"')[1])+int(text[idx2][2])==pos2+int(line.split('\t')[5].split('M')[0]) and text[idx1][0]==text[idx2][0] and int(line.split('\t')[3])+0>=text[idx2][2] and int(line.split('\t')[3])+49-indels2<=text[idx2][4]):
	                    				if pos1<pos2:
								try:
									exon[text[idx1][1].replace("\"", "")+'.'+text[idx2][1].replace("\"", "")+'__'+str(text[idx1][0])]+=1
								except KeyError:
									exon[text[idx1][1].replace("\"", "")+'.'+text[idx2][1].replace("\"", "")+'__'+str(text[idx1][0])]=1
							else:
								try:
									exon[text[idx2][1].replace("\"", "")+'.'+text[idx1][1].replace("\"", "")+'__'+str(text[idx1][0])]+=1
								except KeyError:
									exon[text[idx2][1].replace("\"", "")+'.'+text[idx1][1].replace("\"", "")+'__'+str(text[idx1][0])]=1
	
							# UNCOMMENT 4 lines below for debugging only:
							# debug='ex_3.ex_3__11550'
							# if text[idx1][1].replace("\"", "")+'.'+text[idx2][1].replace("\"", "")+'__'+str(text[idx1][0])==debug or text[idx2][1].replace("\"", "")+'.'+text[idx1][1].replace("\"", "")+'__'+str(text[idx1][0])==debug:
								# f.write(line)
								# f.write(temptext)
				temptext=''
				flag_finish=1
		# COMMENT both lines below for debugging only:
		for k in exon.keys():
			f.write(str(k)+"\t"+str(exon[k])+"\n")

# Chromosome X
	exon={}
	print >> sys.stderr, "Processing Chromosome X"
	print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	cursor.execute("SELECT region_id,exon_name,start,cigar,end FROM hg19 where Chromosome='\"chrX\"' order by start")
	text=cursor.fetchall()
	flag_finish=1
	temptext=''
	for line in file(fileinput+".chrX", 'rU'):
		if flag_finish==1:
			temptext=line
			flag_finish=0
		else:
			if len(line.split('\t')[5].split('N'))>1 or len(temptext.split('\t')[5].split('N'))>1:
				flag_finish=1
				continue
			pos1=int(temptext.split('\t')[3])
			pos2=int(line.split('\t')[3])
			idx1=recurse(pos1,0,len(text))
			idx2=recurse(pos2,0,len(text))
			indels1=0
			indels2=0
			if 'I' in temptext.split('\t')[5]:
				for i in xrange(0,len(temptext.split('\t')[5].split('I'))-1):
					indels1+=int(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N')[len(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N'))-1])
			if 'I' in line.split('\t')[5]:
				for i in xrange(0,len(line.split('\t')[5].split('I'))-1):
					indels2+=int(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N')[len(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N'))-1])
			if 'D' in temptext.split('\t')[5]:
				for i in xrange(0,len(temptext.split('\t')[5].split('D'))-1):
					indels1-=int(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N')[len(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N'))-1])
			if 'D' in line.split('\t')[5]:
				for i in xrange(0,len(line.split('\t')[5].split('D'))-1):
					indels2-=int(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N')[len(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N'))-1])
			if idx1<len(text)-1:
				newidx=idx1
				while text[newidx+1][2]==text[newidx][2]:
					if ('N' in text[idx1][3] and 'N' in temptext.split('\t')[5]):
						if (text[idx1][3].split('M')[1].split('N')[0]==temptext.split('\t')[5].split('M')[len(temptext.split('\t')[5].split('N')[0].split('M'))].split('N')[0]):
							newidx+=1
							break
					newidx+=1
				if int(temptext.split('\t')[3])+0>=text[newidx+1][2] and int(temptext.split('\t')[3])+49-indels1<=text[newidx+1][4]:
					idx1=newidx+1
				newidx=-1
			if idx2<len(text)-1:
				newidx=idx2
				while text[newidx+1][2]==text[newidx][2]:
					if ('N' in text[idx2][3] and 'N' in line.split('\t')[5]):
						if (text[idx2][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[len(line.split('\t')[5].split('N')[0].split('M'))].split('N')[0]):
							newidx+=1
							break
					newidx+=1
				if int(line.split('\t')[3])+0>=text[newidx+1][2] and int(line.split('\t')[3])+49-indels2<=text[newidx+1][4]:
					idx2=newidx+1
				newidx=-1
			if (('N' not in text[idx1][3] and 'N' in temptext.split('\t')[5]) or ('N' in text[idx1][3] and 'N' not in temptext.split('\t')[5])):
				newidx=idx1-1
				while (('N' not in text[newidx][3] and 'N' in temptext.split('\t')[5]) or ('N' in text[newidx][3] and 'N' not in temptext.split('\t')[5])):
					newidx=newidx-1
				if int(temptext.split('\t')[3])+0>=text[newidx][2] and int(temptext.split('\t')[3])+49-indels1<=text[newidx][4]:
					idx1=newidx
			if (('N' not in text[idx2][3] and 'N' in line.split('\t')[5]) or ('N' in text[idx2][3] and 'N' not in line.split('\t')[5])):
				newidx=idx2-1
				while (('N' not in text[newidx][3] and 'N' in line.split('\t')[5]) or ('N' in text[newidx][3] and 'N' not in line.split('\t')[5])):
					newidx=newidx-1
				if int(line.split('\t')[3])+0>=text[newidx][2] and int(line.split('\t')[3])+49-indels2<=text[newidx][4]:
					idx2=newidx
			if ('N' not in text[idx1][3] and 'N' not in temptext.split('\t')[5] and int(temptext.split('\t')[3])+0>=text[idx1][2] and int(temptext.split('\t')[3])+49-indels1<=text[idx1][4]) or ('N' in text[idx1][3] and 'N' in temptext.split('\t')[5] and text[idx1][3].split('M')[1].split('N')[0]==temptext.split('\t')[5].split('M')[len(temptext.split('\t')[5].split('N')[0].split('M'))].split('N')[0] and int(text[idx1][3].split('M')[0].split('\"')[1])+int(text[idx1][2])==pos1+int(temptext.split('\t')[5].split('M')[0]) and int(temptext.split('\t')[3])+0>=text[idx1][2] and int(temptext.split('\t')[3])+49-indels1<=text[idx1][4]):
				if ('N' not in text[idx2][3] and 'N' not in line.split('\t')[5] and text[idx1][0]==text[idx2][0] and int(line.split('\t')[3])+0>=text[idx2][2] and int(line.split('\t')[3])+49-indels2<=text[idx2][4]) or ('N' in text[idx2][3] and 'N' in line.split('\t')[5] and text[idx2][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[len(line.split('\t')[5].split('N')[0].split('M'))].split('N')[0] and int(text[idx2][3].split('M')[0].split('\"')[1])+int(text[idx2][2])==pos2+int(line.split('\t')[5].split('M')[0]) and text[idx1][0]==text[idx2][0] and int(line.split('\t')[3])+0>=text[idx2][2] and int(line.split('\t')[3])+49-indels2<=text[idx2][4]):
					if pos1<pos2:
						try:
							exon[text[idx1][1].replace("\"", "")+'.'+text[idx2][1].replace("\"", "")+'__'+str(text[idx1][0])]+=1
						except KeyError:
							exon[text[idx1][1].replace("\"", "")+'.'+text[idx2][1].replace("\"", "")+'__'+str(text[idx1][0])]=1
					else:
						try:
							exon[text[idx2][1].replace("\"", "")+'.'+text[idx1][1].replace("\"", "")+'__'+str(text[idx1][0])]+=1
						except KeyError:
							exon[text[idx2][1].replace("\"", "")+'.'+text[idx1][1].replace("\"", "")+'__'+str(text[idx1][0])]=1
					# UNCOMMENT 4 lines below for debugging only:
					# debug='ex_3.ex_3__11550'
					# if text[idx1][1].replace("\"", "")+'.'+text[idx2][1].replace("\"", "")+'__'+str(text[idx1][0])==debug or text[idx2][1].replace("\"", "")+'.'+text[idx1][1].replace("\"", "")+'__'+str(text[idx1][0])==debug:
						# f.write(line)
						# f.write(temptext)
			temptext=''
			flag_finish=1

	# COMMENT both lines below for debugging only:
	for k in exon.keys():
		f.write(str(k)+"\t"+str(exon[k])+"\n")
	exon={}
	print >> sys.stderr, "Processing Chromosome Y"
	print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	cursor.execute("SELECT region_id,exon_name,start,cigar,end FROM hg19 where Chromosome='\"chrY\"' order by start")
	text=cursor.fetchall()
	flag_finish=1
	temptext=''
        for line in file(fileinput+".chrX", 'rU'):
                if flag_finish==1:
                        temptext=line
                        flag_finish=0
                else:
                        if len(line.split('\t')[5].split('N'))>1 or len(temptext.split('\t')[5].split('N'))>1:
                                flag_finish=1
                                continue
                        pos1=int(temptext.split('\t')[3])
                        pos2=int(line.split('\t')[3])
                        idx1=recurse(pos1,0,len(text))
                        idx2=recurse(pos2,0,len(text))
                        indels1=0
                        indels2=0
                        if 'I' in temptext.split('\t')[5]:
                                for i in xrange(0,len(temptext.split('\t')[5].split('I'))-1):
                                        indels1+=int(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N')[len(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(temptext.split('\t')[5].split('I')[i].split('M')[len(temptext.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N'))-1])
                        if 'I' in line.split('\t')[5]:
                                for i in xrange(0,len(line.split('\t')[5].split('I'))-1):
                                        indels2+=int(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N')[len(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D')[len(line.split('\t')[5].split('I')[i].split('M')[len(line.split('\t')[5].split('I')[i].split('M'))-1].split('D'))-1].split('N'))-1])
                        if 'D' in temptext.split('\t')[5]:
                                for i in xrange(0,len(temptext.split('\t')[5].split('D'))-1):
                                        indels1-=int(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N')[len(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(temptext.split('\t')[5].split('D')[i].split('M')[len(temptext.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N'))-1])
                        if 'D' in line.split('\t')[5]:
                                for i in xrange(0,len(line.split('\t')[5].split('D'))-1):
                                        indels2-=int(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N')[len(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I')[len(line.split('\t')[5].split('D')[i].split('M')[len(line.split('\t')[5].split('D')[i].split('M'))-1].split('I'))-1].split('N'))-1])
			if idx1<len(text)-1:
                                newidx=idx1
                                while text[newidx+1][2]==text[newidx][2]:
                                        if ('N' in text[idx1][3] and 'N' in temptext.split('\t')[5]):
                                                if (text[idx1][3].split('M')[1].split('N')[0]==temptext.split('\t')[5].split('M')[len(temptext.split('\t')[5].split('N')[0].split('M'))].split('N')[0]):
                                                        newidx+=1
                                                        break
                                        newidx+=1
                                if int(temptext.split('\t')[3])+0>=text[newidx+1][2] and int(temptext.split('\t')[3])+49-indels1<=text[newidx+1][4]:
                                        idx1=newidx+1
                                newidx=-1
                        if idx2<len(text)-1:
                                newidx=idx2
                                while text[newidx+1][2]==text[newidx][2]:
                                        if ('N' in text[idx2][3] and 'N' in line.split('\t')[5]):
                                                if (text[idx2][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[len(line.split('\t')[5].split('N')[0].split('M'))].split('N')[0]):
                                                        newidx+=1
                                                        break
                                        newidx+=1
                                if int(line.split('\t')[3])+0>=text[newidx+1][2] and int(line.split('\t')[3])+49-indels2<=text[newidx+1][4]:
                                        idx2=newidx+1
                                newidx=-1
                        if (('N' not in text[idx1][3] and 'N' in temptext.split('\t')[5]) or ('N' in text[idx1][3] and 'N' not in temptext.split('\t')[5])):
                                newidx=idx1-1
                                while (('N' not in text[newidx][3] and 'N' in temptext.split('\t')[5]) or ('N' in text[newidx][3] and 'N' not in temptext.split('\t')[5])):
                                        newidx=newidx-1
                                if int(temptext.split('\t')[3])+0>=text[newidx][2] and int(temptext.split('\t')[3])+49-indels1<=text[newidx][4]:
                                        idx1=newidx
                        if (('N' not in text[idx2][3] and 'N' in line.split('\t')[5]) or ('N' in text[idx2][3] and 'N' not in line.split('\t')[5])):
                                newidx=idx2-1
                                while (('N' not in text[newidx][3] and 'N' in line.split('\t')[5]) or ('N' in text[newidx][3] and 'N' not in line.split('\t')[5])):
                                        newidx=newidx-1
                                if int(line.split('\t')[3])+0>=text[newidx][2] and int(line.split('\t')[3])+49-indels2<=text[newidx][4]:
                                        idx2=newidx
                        if ('N' not in text[idx1][3] and 'N' not in temptext.split('\t')[5] and int(temptext.split('\t')[3])+0>=text[idx1][2] and int(temptext.split('\t')[3])+49-indels1<=text[idx1][4]) or ('N' in text[idx1][3] and 'N' in temptext.split('\t')[5] and text[idx1][3].split('M')[1].split('N')[0]==temptext.split('\t')[5].split('M')[len(temptext.split('\t')[5].split('N')[0].split('M'))].split('N')[0] and int(text[idx1][3].split('M')[0].split('\"')[1])+int(text[idx1][2])==pos1+int(temptext.split('\t')[5].split('M')[0]) and int(temptext.split('\t')[3])+0>=text[idx1][2] and int(temptext.split('\t')[3])+49-indels1<=text[idx1][4]):
                                if ('N' not in text[idx2][3] and 'N' not in line.split('\t')[5] and text[idx1][0]==text[idx2][0] and int(line.split('\t')[3])+0>=text[idx2][2] and int(line.split('\t')[3])+49-indels2<=text[idx2][4]) or ('N' in text[idx2][3] and 'N' in line.split('\t')[5] and text[idx2][3].split('M')[1].split('N')[0]==line.split('\t')[5].split('M')[len(line.split('\t')[5].split('N')[0].split('M'))].split('N')[0] and int(text[idx2][3].split('M')[0].split('\"')[1])+int(text[idx2][2])==pos2+int(line.split('\t')[5].split('M')[0]) and text[idx1][0]==text[idx2][0] and int(line.split('\t')[3])+0>=text[idx2][2] and int(line.split('\t')[3])+49-indels2<=text[idx2][4]):
                                        if pos1<pos2:
                                                try:
                                                        exon[text[idx1][1].replace("\"", "")+'.'+text[idx2][1].replace("\"", "")+'__'+str(text[idx1][0])]+=1
                                                except KeyError:
                                                        exon[text[idx1][1].replace("\"", "")+'.'+text[idx2][1].replace("\"", "")+'__'+str(text[idx1][0])]=1
                                        else:
                                                try:
                                                        exon[text[idx2][1].replace("\"", "")+'.'+text[idx1][1].replace("\"", "")+'__'+str(text[idx1][0])]+=1
                                                except KeyError:
                                                        exon[text[idx2][1].replace("\"", "")+'.'+text[idx1][1].replace("\"", "")+'__'+str(text[idx1][0])]=1
					# UNCOMMENT 4 lines below for debugging only:
					# debug='ex_3.ex_3__11550'
					# if text[idx1][1].replace("\"", "")+'.'+text[idx2][1].replace("\"", "")+'__'+str(text[idx1][0])==debug or text[idx2][1].replace("\"", "")+'.'+text[idx1][1].replace("\"", "")+'__'+str(text[idx1][0])==debug:
						# f.write(line)
						# f.write(temptext)
			temptext=''
			flag_finish=1

	# COMMENT both lines below for debugging only:
        for k in exon.keys():
                f.write(str(k)+"\t"+str(exon[k])+"\n")
	print >> sys.stderr, datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	db.close()
	f.close()
