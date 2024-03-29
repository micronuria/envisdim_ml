
# DATA COLLECTION 

## Sequences of Envision project

ENVISION sequence data can be found in the The European Nucleotide Archive (ENA) or The Sequence Read Archive (SRA) databases under bioproject ID PRJEB36188.

Using Aspera and SRA-explorer (<https://sra-explorer.info/#>) to download fastq-files.

```{bash}
# Install aspera
wget https://ak-delivery04-mul.dhe.ibm.com/sar/CMA/OSA/0adrj/0/ibm-aspera-connect_4.1.3.93_linux.tar.gz
tar zxvf ibm-aspera-connect_4.1.3.93_linux.tar.gz
bash ibm-aspera-connect_4.1.3.93_linux.sh

# Temporary add the aspera executable to $path
PATH="$HOME/.aspera/connect/bin:$PATH"
```

Generate Aspera commands for downloading all fastq files of the project with SRA-explorer 
that are saved in a script: `sra_explorer_download.sh`

```{bash}
$ head sra_explorer_download.sh 
#!/usr/bin/env bash
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR381/001/ERR3813281/ERR3813281_1.fastq.gz .
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR381/001/ERR3813281/ERR3813281_2.fastq.gz .
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR381/000/ERR3813320/ERR3813320_1.fastq.gz .
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR381/000/ERR3813320/ERR3813320_2.fastq.gz .
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR381/007/ERR3813287/ERR3813287_1.fastq.gz .
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR381/007/ERR3813287/ERR3813287_2.fastq.gz .
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR381/009/ERR3813319/ERR3813319_1.fastq.gz .
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR381/009/ERR3813319/ERR3813319_2.fastq.gz .
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/ERR381/008/ERR3813288/ERR3813288_1.fastq.gz .
```

The data table with the bioproject information was directly downloaded from the SRA website.

For that, go to SRA and search for the project or go to <https://www.ncbi.nlm.nih.gov/sra/?term=PRJEB36188>

Then, click on **Send results to Run selector** and in the **Select box** click on **Download Total-Metadata** (we will need the full table, not only the accession list).
Save table as csv file: `SRARunTable_metadata.csv`.

After visualizing the metadata table in Excel to check all variables included, remove those not needed keeping: study_accession,	sample_accession,	experiment_accession,	run_accession,	tax_id,	scientific_name,	library_name,	read_coun,	fastq_ftp,	submitted_ftp,	sra_ft,	sample_alias,	sample_titleand. Save the file as a tab-delimited file: `filereport_read_run_PRJEB36188.tsv`.

