===============================
Spurious transcription analysis
===============================

We started off by using HISAT2 to map RNA-seq reads against the S. pistillata genome.

These raw files (the important bits are scaffold / position / coverage) were too large to be uploaded to GitHub. Here's how it kinda looks like::

  liewy@kw14764:~/kaust/styl/dna_meth/spurious_transcription/raw_data$ zcat spis_pH7.2A-exon_cov-norm.bed.gz | head
  GeneName        Type    Rank    Strand  Chr     Position        GeneStartPos    GeneLength      NormalizedPos   Coverage
  Spis22565       exon    1       +       Spis.scaffold1000|size75405     748     748     575     0.00173913043478261     3
  Spis22565       exon    1       +       Spis.scaffold1000|size75405     749     748     575     0.00347826086956522     3
  Spis22565       exon    1       +       Spis.scaffold1000|size75405     750     748     575     0.00521739130434783     3
  Spis22565       exon    1       +       Spis.scaffold1000|size75405     751     748     575     0.00695652173913044     3
  Spis22565       exon    1       +       Spis.scaffold1000|size75405     752     748     575     0.00869565217391304     3
  Spis22565       exon    1       +       Spis.scaffold1000|size75405     753     748     575     0.0104347826086957      3
  Spis22565       exon    1       +       Spis.scaffold1000|size75405     754     748     575     0.0121739130434783      3
  Spis22565       exon    1       +       Spis.scaffold1000|size75405     755     748     575     0.0139130434782609      3
  Spis22565       exon    1       +       Spis.scaffold1000|size75405     756     748     575     0.0156521739130435      3

These files were the input for ``parse_mean_coverages.py``. The output of this script are the ``mean_cov.pH_?.??.tsv`` files in this folder.

These tsv files were then fed into ``compile_overall_coverage.py``, to produce ``mean_cov.all.compiled.tsv``.

The plotting script ``plot_dotplot.spurious.meth_vs_unmeth.py`` uses the compiled tsv to produce a plot (which then was prettified with Illustrator).
