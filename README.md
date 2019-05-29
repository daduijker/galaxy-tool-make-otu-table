# galaxy-tool-make-otu-table
A pipeline for clustering and making otu tables in galaxy. This repo can be used for the new (03-04-2019) galaxy 19.01 Naturalis server. The old galaxy 16.04 server is not supported anymore with this pipeline.
## Getting Started
### Prerequisites

**USEARCH**<br />
(user: **galaxy**)  
```
mkdir /home/galaxy/Tools/usearch 
cd /home/galaxy/Tools/usearch
wget [your usearch licence]
mv [your usearch licence] usearch11
chmod 777 /home/galaxy/Tools/usearch/usearch11
```
(user: **ubuntu**)  
```
sudo ln -s /home/galaxy/Tools/usearch/usearch11 /usr/local/bin/usearch11
```
The following tools/packages are needed but included in the conda environment (make_otu_table_environment.yml)
* **VSEARCH**
* **R**
* **DADA2**
* **python2**
* **biopython**

### Installing  
(user: **galaxy**)  
Installing the tool for use in Galaxy
```
cd /home/galaxy/Tools
```
```
git clone https://github.com/naturalis/galaxy-tool-make-otu-table
```
```
chmod 777 galaxy-tool-make-otu-table/make_otu_table.py
```
Create the conda environment
```
conda env create -f make_otu_table_environment.yml
```
Add the following line to /home/galaxy/galaxy/config/tool_conf.xml
```
<tool file="/home/galaxy/Tools/galaxy-tool-make-otu-table/make_otu_table.xml" />
```
If you need to create the conda environment manally:
```
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults
conda create -n dada2env python=3 anaconda
conda activate dada2env
conda install -c bioconda bioconductor-dada2=1.10.0
conda install python=2.7.16
conda install biopython
conda install vsearch
conda deactivate
```
## Source
Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584

Edgar, R.C. (2016), UNOISE2: Improved error-correction for Illumina 16S and ITS amplicon reads.http://dx.doi.org/10.1101/081257

Edgar, R.C. (2013) UPARSE: Highly accurate OTU sequences from microbial amplicon reads, Nature Methods [Pubmed:23955772,  dx.doi.org/10.1038/nmeth.2604].
