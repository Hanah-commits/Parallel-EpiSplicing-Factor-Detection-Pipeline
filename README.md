# EpiSplicing Factor Detection Pipeline

## Usage

```
user@ubuntu:~/EpiSplicing-Factor-Detection-Pipeline/Scripts$ python master.py -w
```

The weight flag (optional) `-w` is used when the expression levels of the splicing factors are used as weights for the binding scores.

## Installations

_Note: If using a conda environemnt, check the R version for which the required R packages exist before creating the conda environment._

1. Install the following tools:

- [MAJIQ](https://bitbucket.org/biociphers/majiq_academic/src/main/)

```
$ pip install git+https://bitbucket.org/biociphers/majiq_academic.git/voila
```

- MANorm
- samtools
- bedtools
- FeatureCounts (if using weights for splicing factor binding scores)
  #TODO add versions

2. Install pipeline requirements.

- Install the python packages from requirements.txt

```
$ pip install -r requirements.txt
```

- Install R packages
```
GenomicFeatures
```

- If using the weight flag `-w` , install the following R packages:

```
$ edger
$ limma
$ install.packages("rjson")
```

3.  Install [RBPMap](http://rbpmap.technion.ac.il/download.html#requirements)

    _Note: [Helper scripts](#Helper-Scripts) for the installation can be found in_

    `~/EpiSplicing-Factor-Detection-Pipeline/Scripts/HelperScripts`

- Download the required files from http://rbpmap.technion.ac.il/download.html#requirements and follow the installation instructions.
- Copy the perl script `~/EpiSplicing-Factor-Detection-Pipeline/Scripts/HelperScripts/RBPmap_EpiSplicing.pl` to your main RBPmap directory. This will be used as the main RBPmap script.
- In `RBPmap_Episplicing.pl`, please edit the paths in the following variables: $scripts_dir, $results_dir to your local *RBPmap* and *EpiSplicing-Factor-Detection-Pipeline* paths, respectively.

## Required Files

1. Reference Genome: GFF3 file, fasta

- Download the gff3 file from gencode. It fits the [requirements of MAJIQ.](https://biociphers.bitbucket.io/majiq/quick.html)
- Download the corresponding .fa file

2.  RNASeq files : .bam and .bam.bai

- Name all the bam and indexed bam files using the following convention:

```
<tissue type>_<identifier>.bam
endodermalcell_ENCFF489LAR.bam
```

3. Create config file for MAJIQ
   Example:

```
[info]
readlen=76
bamdirs=/data/MGP/ERP000591/bam[,/data/MGP2/ERP000591/bam]
genome=mm10
[experiments]
Hippocampus=Hippocampus1,Hippocampus2
Liver=Liver1,Liver2
[optional]
Hippocampus1=strandness:None,
Liver2=strandness:reverse,
```

4 ChIPSeq Files: .bed

- Name all the bed files using the following convention:

```
<histone modification>_<tissue type>_alignment.bed
H3K27ac_ectodermalcell_alignment.bed

<histone modification>_<tissue type>_peak.bed
H3K27ac_ectodermalcell_peak.bed
```

5. `~/EpiSplicing-Factor-Detection-Pipeline/Scripts/paths.json`

- Fill the fields in the .json file. Example:

```json
{
  "tissue1": "ectodermalcell",
  "tissue2": "endodermalcell",
  "Histone modifications": ["H3K27ac", "H3K9me3"],
  "RNASeq files": "/home/user/bam_files",
  "ChIPSeq files": "/home/user/bed_files",
  "Reference genome": "/home/user/gencode.gtf",
  "threads": "8",
  "MAJIQ config": "/home/user/majiq/sample.config",
  "RBPmap directory": "/home/user/RBPmap"
}
```

## Helper Scripts

### RBPMap installation helper scripts

*Note: Please replace the given paths with the appropriate paths in your machine*

1. To download .fa files of required chromosomes.
   ftp.py

```python
import os

for i in range(1,23):
	os.system("wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr"+ str(i) + ".fa.gz' -O chr" +str(i) + ".fa.gz")


for x in ['X', 'Y']:
	os.system("wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr"+ x  + ".fa.gz' -O chr" + x + ".fa.gz")

os.system("gunzip -r *.fa.gz")
```

2. To convert all.fa files to .nib
   allFaToNib.sh

```shell
for i in *.fa; do
  j="/home/ubuntu/RBPmap_1.2/UCSC/hg38/${i%.*}.nib";
  "../faToNib" "/home/ubuntu/RBPmap_1.2/UCSC/hg38/${i}" "${j}";
  done
```
