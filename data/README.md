This directory is to save intermediate files to input to GAP.

`Au5Ag_AFIR.g` and `AuCu4_AFIR.g` are cached preprocessing inputs for the bundled
`Metal/Au5Ag` and `Metal/AuCu4` samples, carried over from the
[development fork](https://github.com/mochimochi-brian/reproduce_rrm/pull/8).
The sequential regression suite uses them so it can run with GAP and Python's
standard library, without repeating the pymatgen preprocessing step.
