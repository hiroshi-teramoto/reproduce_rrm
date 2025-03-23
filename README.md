# reproduce_rrm
This project is to construct a RRM in the shape space from an output of GRRM program. For the detail, refer to our paper, [Reproducing RRM on the Shape Space from its Quotient by CNPI group](https://pubs.acs.org/doi/full/10.1021/acs.jctc.3c00500) by Hiroshi Teramoto, Takuya Saito, Masamitsu Aoki, Burai Murayama, Masato Kobayashi, Takenobu Nakamura, and Tetsuya Taketsugu. In what follows, we denote this paper as "the paper". In Appendix of this paper, there exists a detailed description of this code.

## To run the code, you need to install
* GAP - Groups, Algorithms, Programming (https://www.gap-system.org/)
* python code (rrm_reconstruction_v18.py) depends on 
  - pymatgen
  - numpy
  - networkx
  - some features (*-operator) for python 3.5 or later.
* graphviz (https://graphviz.org/)

## To run the code, download zip file and in the directory, type
./reproduce_rrm_demo.sh

This command computes the RRM of Au5Ag in the shape space from the sample output files in the directory Metal/Au5Ag and output rrm_Au5Ag_AFIR.dot (graphviz dot file) and rrm_Re_pentane_AFIR.png (graph figure). For the detail of the dot format, refer to graphviz (https://graphviz.org/). 

## Options
* vlabel = true or false, if it is set to true, the vertex labels are included in the file rrm_Re_pentane_AFIR.dot. Each vertex label comprises the corresponding EQ number n (EQn in the input file \*EQ_list.log) or n\* if it is an inversion isomer of EQn, and the permutation occuring from the reference structure (EQn or EQn*). 
* elabel = true or false, if it is set to true, the edge labels are included in the file rrm_Re_pentane_AFIR.dot. Each edge label comprises the corresponding TS number n (TSn in the input file \*TS_list.log) or n\* if it is an inversion isomer of TSn, and the permutation occursing from the reference structures (TSn or TSn*).

## Restrictions
* Sample data of GRRM output is in the directory (pentane_Restruct). The files required are ***EQ_list.log, ***TS_list.log, and ***TSn.log
* As is written in the paper, the code does not support RRMs with DC (Dissociation Channel) states and saddle connections. 
* If the input molecule is too big, GAP program (generate_rrm_v9.g) may stop with error. In this case, consider to increase the available memory for GAP. For the detail, see the instruction of GAP (https://www.gap-system.org/)
* If the resulting RRM in shape space is too big, it may take a while for graphviz to visualize the graph. In that case, you might want to consider changing the options of the visualization or using other software to visualize.
* If the output of GRRM contains reaction paths violating [Pechukas theorem](https://pubs.aip.org/aip/jcp/article-abstract/64/4/1516/786979/On-simple-saddle-points-of-a-potential-surface-the) and [its extention (Hiroshi Teramoto, Pontential Energy Function, symmetry and its consequences, Japanese ](https://www.jstc.org/frontier15/), the code outputs assertion errors. For example, in case of AuCu4 (in the sample directory), its output assertion error 
'''
Violation of Pechukus theorem:
TS8:
the reactant and product are permutation isomers
EQ2
in what follows, minus 1 to convert to the atom labels
U(r) (for reactant):
Group( [ () ] )
U(r) (for product):
Group( [ () ] )
stabEQs:
Group( [ (2,4) ] )
U(r) (for transition state):
Group( [ (2,4)(3,5), () ] )
Violation of Pechukus theorem:
TS16:
reactant:
EQ3
product:
EQ5
in what follows, minus 1 to convert to the atom labels
U(r) (for reactant):
Group( [ (2,5)(3,4), () ] )
U(r) (for product):
Group( [ (2,3)(4,5), (2,4)(3,5), (2,5)(3,4), () ] )
U(r) (for transition state):
Group( [ (2,3)(4,5), (2,4)(3,5), (2,5)(3,4), () ] )
GAP computation done.
'''

## Important Parameters
* tolerance - Distance tolerance to consider sites as symmetrically equivalent in rrm_reconstruction_v18.py
  - default value is set to 0.3 (the same as the default value in pymatgen)
  - it may happen that the assigned point group may depend on tol. The current algorithm outputs warning if the assigned point group is different from that assigned by GRRM program. 
