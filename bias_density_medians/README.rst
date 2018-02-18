=========================================================
Calculating CpG bias, meth density and median meth levels
=========================================================

``calc_bias_density_medians.py``, when fed the *S. pistillata* genome and GFF3 (available on http://spis.reefgenomics.org), and the cov files (in ../cov_files), produces files of pattern ``spis_pH_?.??.bias_density_medians.tsv``.

Some command-line magic produced ``compiled_median_meth.tsv`` from the other ``*.tsv`` files. I think I used ``tabulate_tsvs.py *.tsv -c 6``, then manually modified the header line. That script can be found in my ``common`` GitHub repository.
