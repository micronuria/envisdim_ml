
# Sequences pre-processing

Sequences pre-processing prior to ASVs calculation.


## Remove primers

Decompress fastq.gz files of Envision project, Dimension is not compressed.

```{bash}
$ cd analysis/intermediate/raw_16_fastq
$ gzip -d *.gz
```

### Check primers presence in raw sequences

Use *TagCleaner* to check primers in sequences.

Create two files containing all R1 and R2 fastq files respectively:

```{bash}
$ cat *R1*.fastq > envdimR1.fastq
$ cat *R2*.fastq > envdimR2.fastq
```

Move concatenated fastq files to folder containing `tagcleaner.py`
and change to that folder:

```{bash}
$ mv envdimR*.fastq ../../tagcleaner/.
$ cd ../../tagcleaner/
```

In both Envision and Dimension, primers are the same:

* 515F-Y: 5'-GTGYCAGCMGCCGCGGTAA-3'
* 926R: 5'-CCGYCAATTYMTTTRAGTTT-3' 

Predict tag sequences:

```{bash}
$ perl tagcleaner.pl -64 -predict -fastq envdimR1.fastq 

Param  Tag_Sequence    Tag_Length      Percentage_Explained
tag5    GTGNNNNNNNNNNNNNNN      18      98.21
tag3    NNNNNNNNNNNNNNNACG      18      62.37

$ perl tagcleaner.pl -64 -predict -fastq envdimR2.fastq 
Param  Tag_Sequence    Tag_Length      Percentage_Explained
tag5    CCGTCAATTNNNNNNNNNNNNNNN        24      97.86
tag3    TGAGCCGCCTTCGCCACTGGTGTTCCTCCAAATATCTACGAATTTCACCTCTACACTTGGAATT        64      10.50
```

Once primers presence is confirmed, remove concatenated files:

```{bash}
$ rm -r envdimR*.fastq
```

### Trim primers

To trim primers from fastq sequences, use Mothur. 

To use Mothur, sequences must be transformed from fastq file to fasta+qual files:

```{bash}
$ cd ../intermediate/raw_16_fastq/
$ mothur

# First, make file that identifies: sample - read1 - read2
mothur > make.file(inputdir=., type=fastq, prefix=trimming)

# Next, transform fastq files into fasta + qual files 
mothur > fastq.info(file=trimming.files)

mothur > quit()
```

Move files, putting fasta and qual files of each type of read into separated folders.

```{bash}
$ cd ..
$ mv raw_16_fastq/trimming.files .

$ mkdir fasta_qual_R1
$ mkdir fasta_qual_R2

$ mv raw_16_fastq/*R1.fasta fasta_qual_R1/.
$ mv raw_16_fastq/*R1.qual fasta_qual_R1/.
$ mv raw_16_fastq/*R2.fasta fasta_qual_R2/.
$ mv raw_16_fastq/*R2.qual fasta_qual_R2/.
```

Create a text file with the primers information:

`primer_fwd.txt` that contains: forward GTGYCAGCMGCCGCGGTAA

`primer_rev.txt` that contains: forward CCGYCAATTYMTTTRAGTTT

After creating the needed directories, run a bash loop to trim forward primer, and then repeat with the reverse. 

Forward primer:

```{bash}
#!/bin/bash

# Bash script to remove forward primers from R1 seqs 

FQFILES='/data/mcm/nfernandez/envisdim/analysis/intermediate/fasta_qual_R1'
FPRIMER='/data/mcm/nfernandez/envisdim/analysis/intermediate/primer_fwd.txt'
RPRIMER='/data/mcm/nfernandez/envisdim/analysis/intermediate/primer_rev.txt'

TRIMDIR='/data/mcm/nfernandez/envisdim/analysis/intermediate/test'
SCRAPDIR='/data/mcm/nfernandez/envisdim/analysis/intermediate/test'

# create a list of files
cd $FQFILES ; ls *R1.fasta > list_R1files.txt
sed -i 's/\.fasta//g' list_R1files.txt

# Trim forward primers
while read -r file ; do
   mothur "#trim.seqs(fasta=$FQFILES/$file.fasta, qfile=$FQFILES/$file.qual, oligos=$FPRIMER, pdiffs=3, processors=1)"
done < list_R1files.txt

# once it is done, move files
mv $FQFILES/*.trim* $TRIMDIR/.
mv $FQFILES/*.scrap* $SCRAPDIR/.

# clean up
rm $FQFILES/trim.count_table
rm $FQFILES/scrap.count_table

mv $FQFILES/mothur.*.logfile /data/mcm/nfernandez/envisdim/analysis/intermediate/mothur_logfiles/. 
```

Reverse primer:

```{bash}
#!/bin/bash

# Bash script to remove primers from R2 reads 

FQFILES='/data/mcm/nfernandez/envisdim/analysis/intermediate/fasta_qual'
FPRIMER='/data/mcm/nfernandez/envisdim/analysis/intermediate/mothur_analysis/primer_fwd.txt'
RPRIMER='/data/mcm/nfernandez/envisdim/analysis/intermediate/mothur_analysis/primer_rev.txt'

TRIMDIR='/data/mcm/nfernandez/envisdim/analysis/intermediate/fasta_qual_trim'
SCRAPDIR='/data/mcm/nfernandez/envisdim/analysis/intermediate/fasta_qual_scrap'


# Create a list of files
cd $FQFILES ; ls *R2.fasta > list_R2files.txt
sed -i 's/\.fasta//g' list_R2files.txt

# Trim forward primers
while read -r file ; do
   mothur "#trim.seqs(fasta=$FQFILES/$file.fasta, qfile=$FQFILES/$file.qual, oligos=$RPRIMER, pdiffs=3, processors=1)"
done < list_R2files.txt

# once it is done, move files
mv $FQFILES/*.trim* $TRIMDIR/.
mv $FQFILES/*.scrap* $SCRAPDIR/. 

# clean up
rm $FQFILES/trim.count_table
rm $FQFILES/scrap.count_table

mv $FQFILES/mothur.*.logfile /data/mcm/nfernandez/envisdim/analysis/intermediate/mothur_logfiles/. 
```

Check that primers have been removed with TagCleaner:

* 515F-Y: 5'-GTGYCAGCMGCCGCGGTAA-3'
* 926R: 5'-CCGYCAATTYMTTTRAGTTT-3' 

```{bash}
$ cd /data/mcm/nfernandez/envisdim/analysis/intermediate/fasta_qual_trim
$ cat *_R1.trim.fasta > R1trimall.fasta
$ cat *_R2.trim.fasta > R2trimall.fasta
$ mv R*all.fasta ../../tagcleaner/.
$ cd ../../tagcleaner/
$ perl tagcleaner.pl -predict -fasta R1trimall.fasta 

#Param  Tag_Sequence    Tag_Length      Percentage_Explained
tag5    TACGGAGGGNNNNNNNNNNNNNNN        24      81.50
tag3    AGGATTAGATACCCTGGTAGTCCACGCCGTAA        32      31.96

$ perl tagcleaner.pl -predict -fasta R2trimall.fasta 

#Param  Tag_Sequence    Tag_Length      Percentage_Explained
tag5    TAANNNNNNNNNNNNNNN      18      28.43
tag3    NNNNNNNNNNNNNNN 15      36.00
```


## QC filtering and merge

To remove sequencing errors and merge paired reads we are using Moira 

### Prepare files to work with MOIRA 

#### Create an accnoss file

We need to create, for each pair of fasta files, a list with the name of the paired sequences in which primers have been successfully removed in both R1 and R2.

For that, we first create a list of samples:
```{bash}
cd /data/mcm/nfernandez/envisdim/analysis/intermediate/fasta_qual_trim
# make a sample list
ls *R1.trim.fasta > list_samples.txt
sed -i 's/_R1\.trim\.fasta//g' list_samples.txt
```

And then use a python script to create those lists, which are called accnos files:
The script `make_accnoss.py` is based in code created by o, Fernando Puente-Sánchez. 

```{python}
import os
os.getcwd()


lista = open("list_samples.txt").readlines()

for x, sample in enumerate(lista): 
    read1 = "_".join([sample.rstrip(), 'R1.trim.fasta'])
    read2 = "_".join([sample.rstrip(), 'R2.trim.fasta'])
    accnos = ".".join([sample.rstrip(), 'trim.common.accnos'])
    
    print(x)
    
    hf=set()
    hr=set()
    
    for line in open(read1):
        if line.startswith('>'):
            hf.add( line.strip().lstrip('>').split('\t')[0] )
    
    for line in open(read2):
        if line.startswith('>'):
            hr.add( line.strip().lstrip('>').split('\t')[0] )
    
    with open(accnos, 'w') as outfile:
        for h in hf & hr:
            outfile.write('{}\n'.format(h))

```

An example of the accnos file is:
```{bash}
$ head DIMprok19.trim.common.accnos 
M02696_99_000000000-ALTJF_1_1106_6475_16135
M02696_99_000000000-ALTJF_1_1118_4634_12531
M02696_99_000000000-ALTJF_1_2101_18214_13581
M02696_99_000000000-ALTJF_1_1110_19836_21701
M02696_99_000000000-ALTJF_1_2108_20371_13110
M02696_99_000000000-ALTJF_1_2104_12595_7335
M02696_99_000000000-ALTJF_1_2104_18825_14625
M02696_99_000000000-ALTJF_1_1103_3145_7852
M02696_99_000000000-ALTJF_1_1112_18980_12464
M02696_99_000000000-ALTJF_1_2112_21703_18395
```

#### Pick good sequences

Next, we need to select the trimmed sequence pairs included in the accnos files.
For that, we can use the Mothur `get.seqs` function in a Bash script to loop over all files.
After creating the needed directories, run the script:

```{bash}
#!/bin/bash

# pick sequences included in the accnos files (good pair of reads after trimming both primers)

# define paths:
# directory with input files
FQFILES='/data/mcm/nfernandez/envisdim/analysis/intermediate/fasta_qual_trim'
# text file with sample names
SAMPLES='/data/mcm/nfernandez/envisdim/analysis/intermediate/fasta_qual_trim/list_samples.txt'
# Output directory for results
OUTDIR='/data/mcm/nfernandez/envisdim/analysis/intermediate/fasta_pick'
# Directory for accnos files
ACNNDIR='/data/mcm/nfernandez/envisdim/analysis/intermediate/accnos_files'
# Directories to move the mothur logfiles (from -> to)
LOGS='/data/mcm/nfernandez/envisdim/analysis/scripts'
LOGSDIR='/data/mcm/nfernandez/envisdim/analysis/intermediate/mothur_logfiles'

while read -r file ; do
    mothur "#get.seqs(fasta=${FQFILES}/${file}_R1.trim.fasta, qfile=${FQFILES}/${file}_R1.trim.qual, accnos=${FQFILES}/${file}.trim.common.accnos)"
    mothur "#get.seqs(fasta=${FQFILES}/${file}_R2.trim.fasta, qfile=${FQFILES}/${file}_R2.trim.qual, accnos=${FQFILES}/${file}.trim.common.accnos)"

done < $SAMPLES

# Clean up, 
mv $FQFILES/*pick* $OUTDIR/.
mv $FQFILES/*.accnos $ACCNDIR/.
mv $LOGS/*.logfile $LOGSDIR/.
```


#### Translate sequences format from fasta to fastq

Although Moira can use both fasta and fastq formats, from now on, we are going to use fastq format. 
We use Mothur function `make.fastq` function to translate the selected sequences back to fastq.

```{bash}
#!/bin/bash

# Translate trimmed sequences from fasta to fastq 

# define paths
# Input fasta and qual files directory
PFILES='/data/mcm/nfernandez/envisdim/analysis/intermediate/fasta_pick'
# list of samples
SAMPLES='/data/mcm/nfernandez/envisdim/analysis/intermediate/fasta_qual_trim/list_samples.txt'
# Directory for results
OUTDIR='/data/mcm/nfernandez/envisdim/analysis/intermediate/fastq_pick'
# Directories for logfiles
LOGS='/data/mcm/nfernandez/envisdim/analysis/scripts'
LOGSDIR='/data/mcm/nfernandez/envisdim/analysis/intermediate/mothur_logfiles'

while read -r file ; do
    mothur "#make.fastq(fasta=${PFILES}/${file}_R1.trim.pick.fasta , qfile=${PFILES}/${file}_R1.trim.pick.qual)"
    mothur "#make.fastq(fasta=${PFILES}/${file}_R2.trim.pick.fasta , qfile=${PFILES}/${file}_R2.trim.pick.qual)"
done < $SAMPLES

# Clean up
mv $PFILES/*.fastq $OUTDIR/.
mv $LOGS/*.logfile $LOGSDIR/.
```


### Reads QC filter and merge into contigs

The next step is to quality-filter reads with Moira, that also builts contigs from the paired-end reads.
For it Moira applies a Poisson binomial distribution to estimate sequencing errors.

To run Moira, we have to joing all selected sequences in a single file depending their sense.

```{bash}
$ cat *_R1*.fastq > all_pick_R1.fastq
$ cat *_R2*.fastq > all_pick_R2.fastq
```

With those files, we can launch Moira after creating the corresponding folders.

```{bash}
nohup moira.py --forward_fastq=all_pick_R1.fastq --reverse_fastq=all_pick_R2.fastq --paired --output_format fastq --consensus_qscore posterior --qscore_cap 0 --collapse TRUE 
 --ambigs disallow --processors 12 --output_prefix moira_primerfree_ > moira.primerfree.logfile &
```

Moira output consist on unique sequences, therefore to properly calculate ASVs and their count table, we have to make fastq files with the corresponding abundances of each unique sequence.
To achieve it, first we make a file of groups of identical sequences using the information we already have in the accnos files:
`make_grups.py`
```{python}

import os
import csv

os.getcwd()


lista = open("list_samples.txt").readlines()

my_dic = {}
for x, sample in enumerate(lista): 
    accnos = "".join([sample.rstrip(), '.trim.common.accnos'])
    for line in open(accnos):
        my_dic[line.strip()] = sample.strip()

with open('envisdim.groups', 'w') as f:
    for key in my_dic.keys():
        f.write("%s\t%s\n"%(key,my_dic[key]))
```

The file created contains two columns, one with the read names and a second column with the sample they belong to: 
```{bash}
$ head envisdim.groups
M02696_99_000000000-ALTJF_1_2115_10748_5304	DIMprok10
M02696_99_000000000-ALTJF_1_2110_3905_8114	DIMprok10
M02696_99_000000000-ALTJF_1_1117_10018_12632	DIMprok10
M02696_99_000000000-ALTJF_1_1116_16928_20985	DIMprok10
```


Now, we can use that information to build the required fastq.
Script based on Diego Jimenez and Fernando Puente-Sanchez code.

```{python}
# Import Operative Systems tools
import os

# Get the Current Working Directory
os.getcwd()

# Make avaible info in moira output.fastq file
dict_moira = {}
with open('../fastq_moira_output/moira_primerfree_.qc.good.fastq','r') as infile:
    while True:
        holoseq = infile.readline().lstrip("@").rstrip("\n")
        if not holoseq:
            break
        ntseq = infile.readline()
        plus = infile.readline()
        qual = infile.readline().rstrip("\n")
        seq_info = ntseq + plus + qual
        dict_moira[holoseq] = seq_info


# Make avaible info in moira.names file
dict_names = {}
with open('../fastq_moira_output/moira_primerfree_.qc.good.names','r') as infile:
    for line in infile:
        holoseq, names = line.rstrip("\n").split("\t")
        list_names = names.split(",")
        for name_seq in list_names:
            dict_names[name_seq] = holoseq

# Make avaible the info in groups file
dict_groups = {}
with open('../accnos_files/envisdim.groups','r') as infile:
    for line in infile:
        name_seq, sample_id = line.rstrip("\n").split("\t")
        dict_groups[name_seq] = sample_id


# Open and Create the new files once only using a set to constraint these values in a dictionary where we keep the new files 'route' to entry them
output_files = {}
for sample_id in set(dict_groups.values()):
	output_files[sample_id] = open('moira.resample.primerfree_{}.preR.fastq'.format(sample_id),'w')


# Write the files according to the request task
for name_seq, holoseq in dict_names.items():
	sample_id = dict_groups[name_seq]
	output_files[sample_id].write('@{}\n{}\n'.format(name_seq, dict_moira[holoseq]))

for sample_id in output_files:
	output_files[sample_id].close()

```

