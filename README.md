# LoReMINE: Long Read-based Microbial genome mining pipeline

LoreMiNE (Long Read-based Microbial genome mining pipeline) is an end-to-end pipeline for microbial natural product discovery directly from long-read sequencing data. LoreMiNE integrates multiple modules into a unified pipeline: (i) de-novo genome assembly using multiple long-read assemblers with automated selection of the highest-quality assembly, (ii) taxonomic classification, (iii) BGC detection, and (iv) clustering of identified BGCs to facilitate comparative analysis across different databases and taxa. By combining these steps into a unified pipeline, LoreMiNE facilitates comprehensive exploration of natural product diversity which has the potential to yield novel drug candidates.


## Installation

LoReMINE can be used by installing it via [conda](https://anaconda.org/kalininalab/loremine) 

To install it into the new environment:

`````shell
conda create -n loremine -c conda-forge -c kalininalab -c bioconda loremine
conda activate loremine
`````

Or to install it into the already existing environment (for e.g: loremine):

`````shell
conda create -n loremine
conda activate
conda install -c conda-forge -c kalininalab -c bioconda loremine
`````
We still need to install [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE) before using the LoReMINE pipeline. You can install it anywhere, but we prefer to install in the same conda environment directory (````` ~/anaconda3/envs/loremine/ `````) so that all tools are in one place. To install [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE), follow the instructions given below

`````shell
cd ~/anaconda3/envs/loremine/
git clone https://github.com/medema-group/BiG-SCAPE
cd BiG-SCAPE
pip install .
`````

## Setting up local databases
