# reproduce_rrm
This project is to construct a Reaction Route Map (RRM) in the shape space from an output of Global Reaction Route Mapping (GRRM) program. GRRM finds reaction pathways but identifies CNPI isomers as one; this tool reconstructs the full reaction network in the shape space accounting for all permutations of identical atoms. This code also can be used to verify your GRRM outputs. For more details, refer to our paper (Hiroshi Teramoto et al., J. Chem. Theory Comput. 2023, 19, 17, 5886–5896), [Reproducing RRM on the Shape Space from its Quotient by CNPI group](https://pubs.acs.org/doi/full/10.1021/acs.jctc.3c00500) by Hiroshi Teramoto, Takuya Saito, Masamitsu Aoki, Burai Murayama, Masato Kobayashi, Takenobu Nakamura, and Tetsuya Taketsugu. Hereafter, we refer to this publication as ‘the paper.’ The appendix of the paper provides a detailed description of this code.

## To run the code, you need to install (Tested with Ubuntu 22.04.5 LTS)
* GAP - Groups, Algorithms, Programming (https://www.gap-system.org/) (Tested with GAP 4.11 (or 4.13))
* The Python script `rrm_reconstruction_v18.py` depends on the following packages:
  - pymatgen (Tested with 2024.8.9)
  - numpy (Tested with 1.26.4)
  - networkx (Tested with 3.3)
  - some features (*-operator) for Python 3.5 or later (Tested with 3.10.12).
* Graphviz (https://graphviz.org/) (Tested with 2.43.0 (0))

## To run the code, download zip file and in the directory
* Edit the script `reproduce_rrm_demo.sh`
  - The Python command supposed to be `python3`. If that is not the case, modify the command appropriately.
  - The path for GAP is supposed to be `/usr/local/gap-4.13.0/gap`. Modify the path to match your environment path. 
* Type `./reproduce_rrm_demo.sh` and press the enter key.

This command computes the RRM of Au5Ag in the shape space from the sample output files in the directory `Metal/Au5Ag` and outputs `rrm_Au5Ag_AFIR.dot` (a Graphviz DOT file) and `rrm_Au5Ag_AFIR.png` (the rendered graph image). For the detail of the dot format, refer to graphviz (https://graphviz.org/). If the code ran correctly, you should see the output messages like
```
Metal/Au5Ag//Au5Ag_AFIR_EQ_list.log
Metal/Au5Ag//Au5Ag_AFIR_TS_list.log
Metal/Au5Ag//Au5Ag_AFIR_TS
241 lines read
---with several warnings indicating discripancy between point-group identifications of GRRM and pymatgen)
ur computation finished!!!
urt computation finished!!!
connections computation finished!!!
output GRRM graph as a dot file!!
symc computation finished!!!
the index of symc in sym : 1
GAP computation done.
All vertices with the same number have consistent degrees.
DONE
```
the figure like: ![RRM of Au5Ag cluster](./rrm_Au5Ag_AFIR.png) Sometimes it may be too complicated to visualize the resulting RRM in shape space. In that case, you should edit the corresponding dot file or extract some of the characteristics to quantify some of the graph properties. The line
```
the index of symc in sym : 1
```
indicates that the resulting RRM in shape space is connected. In general, if the output contains the line,
```
the index of symc in sym : n
```
it indicates that the resulting RRM in shape space has `n` connected components (that are mutually graph-isomorphic with each other. See the paper and its supplimentary material for details.). If the latter is the case, this codes output one of the connected components.

## Constitution of the code
* Core Python script `rrm_reconstruction_v18.py` that handles parsing GRRM output and preparing data
* GAP script `generate_rrm_v11.g` that performs group-theoretic computations
* Helper Python script for validation `check_number_of_edges_v3.py`
* Shell script to tie it all together `reproduce_rrm_demo.sh`

## Options
* `vlabel = true or false`, if it is set to true, the vertex labels are included in the file `rrm_Au5Ag_AFIR.dot`. Each vertex label comprises the corresponding EQ number n (EQn in the input file \*EQ_list.log) or n\* if it is an inversion isomer of EQn, and the permutation from the reference structure (EQn or EQn*). 
* `elabel = true or false`, if it is set to true, the edge labels are included in the file `rrm_Au5Ag_AFIR.dot`. Each edge label comprises the corresponding TS number n (TSn in the input file \*TS_list.log) or n\* if it is an inversion isomer of TSn, and the permutation from the reference structures (TSn or TSn*).

## Restrictions
* Sample data of GRRM output is in the directory Metal. The files required are `***EQ_list.log`, `***TS_list.log`, and `***TSn.log` (`n` is the indices of the transition states.).
* As mentioned in the paper, the code does not support RRMs that include DC (dissociation channel) states or saddle connections.
* The current version of the code only accept the connected RRMs as inputs (otherwise the assertion error `assert nx.is_connected(G)` occurs in rrm_reconstruction_v18.py.
* The GAP program may stop with an error if the input molecule is too large. In this case, consider increasing the memory available to GAP. For details, see the GAP documentation. (https://www.gap-system.org/)
* If the resulting RRM in shape space is too big, it may take a while for Graphviz to visualize the graph. In that case, you might want to consider changing the options of the visualization or using other software to visualize.
* If the GRRM output contains reaction paths that violate [Pechukas's theorem](https://pubs.aip.org/aip/jcp/article-abstract/64/4/1516/786979/On-simple-saddle-points-of-a-potential-surface-the) (and [its extension (Hiroshi Teramoto, Pontential Energy Function, symmetry and its consequences, in Japanese)](https://www.jstc.org/frontier15/)), the code will throw an assertion error. For example, in the provided AuCu4 sample, TS8 and TS16 produce such errors​:

```
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
```
If this occurs, the results are not guaranteed to be correct – you should carefully examine your GRRM output in such cases.

## Important Parameters
* tolerance - Distance tolerance to consider sites as symmetrically equivalent in rrm_reconstruction_v18.py
  - default value is set to 0.3 (the same as the default value in pymatgen)
  - it may happen that the assigned point group may depend on tol. The current algorithm outputs warning if the assigned point group is different from that assigned by GRRM program. 

## Using your own GRRM outputs
* Create a directory with the name of your molecule (e.g. Au5Ag). In what follows, we suppose it is ${MOL}.
* Put all the output files of GRRM, `${MOL}_AFIR_EQ_list.log`, `${MOL}_AFIR_TS_list.log`, `${MOL}_AFIR_TSn.log` (`n` is supposed to be the indices of the transition states.) under the directory.
* Modify `MOL=${MOL}` in the `reproduce_rrm_demo.sh`.
* Run `./reproduce_rrm_demo.sh`.
* Watch out warnings and errors. If Assertion error occurred, it indicates there is a bug in this code (in that case, kindly report the bug to us!) or there is a problem in your GRRM output (like the case AuCu4 mentioned above, we observed the violation of Pechukas theorem occurred in case if Vallay-Ridge transitions occur in the middle of a reaction path or other possibly more primitive error.). This code can also used to verify your GRRM output.
* If the code ran successfuly, it will output `rrm_${MOL}_AFIR.dot` and `rrm_${MOL}_AFIR.png` (and `data/${MOL}_AFIR.g` for an intermediate file). If the png figure is too complicated to show, consider extracting some features of the graph from the dot file. For example, we use persistent homology to extract some features of output graphs. 
