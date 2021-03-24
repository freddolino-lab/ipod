Here is the list os current recommendations for standardizing the IPOD v 2 code.

- Make test use dataset

# General #

- Retab everything
- Update to py3
- Add support for numba.jit(nopython=True)
    - Especially helpful for bootstrapping and summary funcs
- add support for multiple chromosomed organisms

# TO DO #

## Main driver ##

- run_all_driver.py

## Alignment ##

- run_all_alignments.py 
- run_preprocessing_alignment_dna_onedir.py
    - switch str formatting
    - upgrade python version:

## Bootstrapping ##

- bootstrap_sam_file.py
    - modify sys.path.append('mbwolfe') stuff at top of file
    - switch to numba.jit(nopython=True) for map step
    - add multi-chromosome support
    - add strandedness
        - idea here, add dimension to array for strand,
          encode strand in 
