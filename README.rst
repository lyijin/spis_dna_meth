===============================================================================================
Epigenome-associated phenotypic acclimatization to ocean acidification in a reef-building coral
===============================================================================================

bioRxiv link to manuscript: https://www.biorxiv.org/content/early/2017/09/13/188227.full.pdf

Methylation pipeline
--------------------
A fuller description of the pipeline used (+ scripts, + theoretical considerations) is detailed at https://github.com/lyijin/working_with_dna_meth.

Bismark produces ``*.cov`` files, which are annotated by that pipeline (specifically, ``annotate_bismark_cov.py``) to produce ``*.annot.cov`` files. These files are the key files used in many analyses, and can be found in the ``cov_files`` folder, one for each sample.

Brief description of folder contents
------------------------------------
I have uploaded scripts that carry out the key analyses--but I'm fairly sure I've missed a couple. Email me if you'd like to know more about how I carried out certain things, and I'll upload them here.

Here, I provide brief descriptions of the folders I have made available. A longer explanation can be found within the folders themselves.

``bias_density_medians``: calculates the CpG bias, methylation density, and median methylation level on a per-gene basis.

``GLM``: how we defined differentially methylated genes.

``jnk_mapk_meth``: code to illustrate the hypermethylation of negative regulators of JNK and MAPK pathways, and the general hypomethylation of positive regulators of the same pathways.

``spurious_transcription``: checks whether methylation reduced spurious transcription (it did).

Stuff used in this project, but uploaded elsewhere
--------------------------------------------------
I designed qPCR primers using ``design_generic_primers.py``, hosted at https://github.com/lyijin/common. This script depends on other scripts i.e. ``natural_sort.py``, ``parse_fasta.py`` and ``primer3_api.py``, all of which are in the same folder. Note that primer3 has to be installed on your system for ``primer3_api.py`` to work (and there's a few hardcoded stuff in that script, which you need to change for it to work).

Similarly, I used ``design_miseq_amplicon_primers.py`` to design MiSeq amplicons for the amplicon-specific bisulphite sequencing. As bisulphite conversion causes the top and bottom strand to no longer be reverse-complements of each other, it was a pain to manually design primers to target a specific strand. This complication was what pushed me to write stuff to automated primer designs: other scripts with similar names are variants of this original script!
