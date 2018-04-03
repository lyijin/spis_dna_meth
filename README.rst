===============================================================================================
Epigenome-associated phenotypic acclimatization to ocean acidification in a reef-building coral
===============================================================================================

bioRxiv link to manuscript: https://www.biorxiv.org/content/early/2017/09/13/188227.full.pdf

Methylation pipeline
--------------------
A fuller description of the pipeline used (+ scripts, + theoretical considerations) is detailed at https://github.com/lyijin/working_with_dna_meth.

Bismark produces ``*.cov`` files, which are annotated by that pipeline (specifically, ``annotate_bismark_cov.py``) to produce ``*.annot.cov`` files. These files are the key files used in many analyses, and can be found in this folder, one for each sample. To replicate analyses in this repo, these files should be decompressed (i.e. from ``*.cov.gz`` to ``*.cov``).

Brief description of folder contents
------------------------------------
I have uploaded scripts that carry out the key analyses--but it's likely that I've missed a few. Email me if you'd like to know more about how I carried out certain things, and I'll upload them here.

Here, I provide brief descriptions of the folders I have made available. A longer explanation can be found within the folders themselves.

``bias_density_medians/``: calculates the CpG bias, methylation density, and median methylation level on a per-gene basis.

``GLM/``: how we defined differentially methylated genes.

``jnk_mapk_meth/``: code to illustrate the hypermethylation of negative regulators of JNK and MAPK pathways, and the general hypomethylation of positive regulators of the same pathways.

``primer_design/``: scripts I wrote to help (i.e. brute force) possibly optimal primer pairs from regions of interest.

``spurious_transcription/``: checks whether methylation reduced spurious transcription (it did).
