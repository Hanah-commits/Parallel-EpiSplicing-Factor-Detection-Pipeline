import os
#for i in range(0,23):
#    os.system("rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr" + str(i)+ ".fa.gz")

#for x in ['X', 'Y']:
#    os.system("rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr" + x + ".fa.gz")


for i in range(1,23):
	os.system("wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr"+ str(i) + ".fa.gz' -O chr" +str(i) + ".fa.gz")


for x in ['X', 'Y']:
	os.system("wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr"+ x  + ".fa.gz' -O chr" + x + ".fa.gz")

os.system("gunzip -r *.fa.gz")
