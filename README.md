# reproduce_rrm
This project is to construct a RRM in the shape space from an output of GRRM program. For the detail, refer to our paper, "Reproducing RRM on the Shape Space from its Quotient by CNPI group" by Hiroshi Teramoto, Takuya Saito, Masamitsu Aoki, Burai Murayama, Masato Kobayashi, Takenobu Nakamura, and Tetsuya Taketsugu. In what follows, we denote this paper as "the paper". In Appendix of this paper, there exists a detailed description of this code.

## To run the code, you need to install
* GAP - Groups, Algorithms, Programming (https://www.gap-system.org/)
* python code (rrm_reconstruction_v12.py) depends on 
  - pymatgen
  - numpy
  - networkx
* graphviz (https://graphviz.org/)

## To run the code, download zip file and in the file, type
./reproduce_rrm_demo.sh

This command computes the RRM of pentane in the shape space from the sample output files in the directory (pentane_Restruct) and output rrm_Re_pentane_AFIR.dot and rrm_Re_pentane_AFIR.png. For the detail of the dot format, refer to graphviz (https://graphviz.org/).

## Restrictions
* As is written in the paper, the code does not support RRMs with DC states and saddle connections. 
* If the input molecule is too big, GAP program (generate_rrm_v9.g) may stop with error. In this case, consider to increase the available memory for GAP. For the detail, see the instruction of GAP (https://www.gap-system.org/)

## Important Parameters
* tolerance - Distance tolerance to consider sites as symmetrically equivalent in rrm_reconstruction_v12.py
  - default value is set to 0.3 (the same as the default value in pymatgen)
  - it may happen that the assigned point group may depend on tol. The current algorithm output warning if the assigned point group is different from that assigned by GRRM program. 
