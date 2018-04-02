=============================
Python-assisted primer design
=============================

Dependencies
------------
The ``design_*.py`` files are the ones doing the heavy lifting, but both of them require three scripts (provided here. in the same folder): ``natural_sort.py``, ``parse_fasta.py`` and ``primer3_api.py``.

Note that ``primer3`` has to be installed on your system for ``primer3_api.py`` to work (and there's a few hardcoded stuff in that script, which you need to change for it to work).

Usage
-----
Both scripts are very similar to each other--I wrote the MiSeq script (``design_miseq_amplicon_primers.py``) first, as designing appropriate primers gets fairly complicated when one has to factor in the C-to-T conversion due to bisulphite treatment.

The file ``miseq_amplicon.in`` shows an example of an input file that I feed into the script, and ``miseq_amplicon.designed_primers.xlsx`` is the tidied output file, which includes the appropriately labelled primers I ordered from Sigma.

Later on, I wondered whether I could extend this laziness to designing normal primers. I felt that the time spent writing the code would be easily saved by being able to automate all future primer designs (and if my colleagues use them, it'll save their time as well). This resulted in ``design_generic_primers.py``.

``qpcr.in`` was the input file that I fed into the script to design qPCR primers; ``qpcr.w10.out`` is one of the output files produced when I set the "wobble" parameter to 10. I tend to run the script concurrently with three wobble values: 5, 10 and 20, and choose primer pairs that look the most appropriate. Yes, it's a bit voodoo-ish, but that's what primer design is usually anyway! :p