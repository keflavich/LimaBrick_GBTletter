Reduction Scripts
-----------------

The `do_all.py` script will run all stages of data reduction required to create
the optical depth cubes described in the paper.

In order for it to work properly, you must edit the `paths.py` file with the
correct input and output paths.  The inputs are the paths to the two GBT data
files, which can be acquired from the NRAO Archive.  The files are 736 MB
(AGBT12B_221, C-band) and 3.3 GB (AGBT14A_110, Ku-band) plus an additional few MB for the index files.
