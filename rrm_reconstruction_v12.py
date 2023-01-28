import pymatgen.core
from pymatgen.core.structure import Molecule
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.analysis.molecule_matcher import GeneticOrderMatcher
from pymatgen.analysis.molecule_matcher import HungarianOrderMatcher
import re
import math
import numpy as np
import sys
import networkx as nx
import warnings

# tolerance - Distance tolerance to consider sites as symmetrically equivalent. Defaults to 0.3 Angstrom. 
tol = 0.3

# compute the Euclidean distance between cor1 and cor2
def dist(cor1, cor2):
    assert len(cor1) == len(cor2)
    s = 0
    for i in range(len(cor1)):
        s += (cor1[i]-cor2[i])**2
    return math.sqrt(s)

# compute the closest atom in list2 for each atom in list1
# for example, plist = [2,3,1,4,0] means that the closest atom in list2 to the 0th atom in list1 is 2, that to the 1st atom in list1 is 3 etc.
def perm(list1,list2):
    assert len(list1) == len(list2)
    plist = []
    rlist = list(range(len(list2)))
    for i in range(len(list1)):
        dmin = 1.0
        for j in rlist:
            if dmin > dist(list1[i],list2[j]):
                jmin = j
                dmin = dist(list1[i],list2[j])

        assert dmin <= tol 
        plist.append(jmin)
        rlist.remove(jmin)

    return plist

# judge whether the set of SymmOps sop contains a volume-inversion operation.
# returns true if the molecule is chiral
def wchiral(sop):
    for sp in sop:
       if np.linalg.det(sp.rotation_matrix) < 0:
           return False
    return True

# apply inversion to the molecule mol
def invmol(mol):
    molsym = PointGroupAnalyzer(mol,tolerance=tol)
    imol = mol.copy()
    imol.apply_operation(molsym.inversion_op)
    return imol 

# finding all the fundamental cycles of a given multigraph
def cycles(G,root):
    tG = nx.MultiGraph(nx.minimum_spanning_edges(G,data=True))
    Gc = G.copy()
    Gc.remove_edges_from(tG.edges())
    fedges=Gc.edges(data=True)

    cycle = []
    for e in fedges:
        cyclee = []
        path = nx.shortest_path(tG,root,e[0])
        for i in range(len(path)-1):
            cyclee.append((path[i],path[i+1],*list(tG[path[i]][path[i+1]].values())))
        cyclee.append(e)
        path = nx.shortest_path(tG,e[1],root)
        for i in range(len(path)-1):
            cyclee.append((path[i],path[i+1],*list(tG[path[i]][path[i+1]].values())))

        cycle.append(cyclee)
    return cycle

args = sys.argv

# loading the equilibrium structures
feqlist = open(args[1],'r')
lines = feqlist.readlines()
print(len(lines),'lines read')
moleq = []
wstructure = False
for line in lines:
    if re.search(r'Geometry',line):
        wstructure = True
        geo = line.split()[3] + line.split()[4]
        sym_grrm = line.split()[7]
        atoms = []
        coords = []
    elif wstructure:
        if re.search(r'Energy',line):
            wstructure = False 
            mol = Molecule(atoms,coords).get_centered_molecule()

            molsym = PointGroupAnalyzer(mol,tolerance=tol)
            if str(molsym.get_pointgroup()) != sym_grrm:
                warnings.warn('{0} sym (grrm) = {1}, sym (pymatgen) = {2}'.format(geo,sym_grrm,molsym.get_pointgroup()))

            moleq.append(mol)
        else:
            m = line.split()
            atoms.append(m[0])
            coords.append([float(m[1]),float(m[2]),float(m[3])])

# loading the transition structures
ftslist = open(args[2],'r')
lines = ftslist.readlines()
print(len(lines),'lines read')
molts = []
gconns = []
wstructure = False
for line in lines:
    if re.search(r'Geometry',line):
        wstructure = True
        geo = line.split()[3] + line.split()[4]
        sym_grrm = line.split()[7]
        atoms = []
        coords = []
    elif re.search(r'CONNECTION',line):
        m = line.split()
        gconns.append([int(m[2]),int(m[4])])
    elif wstructure:
        if re.search(r'Energy',line):
            wstructure = False 
            mol = Molecule(atoms,coords).get_centered_molecule()

            molsym = PointGroupAnalyzer(mol,tolerance=tol)
            if str(molsym.get_pointgroup()) != sym_grrm:
                warnings.warn('{0} sym (grrm) = {1}, sym (pymatgen) = {2}'.format(geo,sym_grrm,molsym.get_pointgroup()))

            molts.append(mol)
        else:
            m = line.split()
            atoms.append(m[0])
            coords.append([float(m[1]),float(m[2]),float(m[3])])

# loading irc
ftsdata = args[3] 
molrp = []
for i in range(len(molts)):
    with open(ftsdata+str(i)+'.log',mode='r') as f:
        lines = f.readlines()
        wstructure = False
        mols = []
        for line in lines:
            if re.search(r'Optimized structure',line):
                wstructure = True
                atoms = []
                coords = []
            elif wstructure:
                if re.search(r'ENERGY',line):
                    wstructure = False
                    mols.append(Molecule(atoms,coords).get_centered_molecule())
                else:
                    m = line.split()
                    atoms.append(m[0])
                    coords.append([float(m[1]),float(m[2]),float(m[3])]) 
        molrp.append(mols)

fglist = args[4]

G = nx.MultiGraph()

# computing ur
neq = len(moleq)
inv = [] 
org_eq = list(range(neq)) 
for ind, mol in enumerate(moleq):
    molsym = PointGroupAnalyzer(mol,tolerance=tol)
    sop = molsym.get_symmetry_operations()
    if ind < neq and wchiral(sop):
        moleq.append(invmol(mol))
        inv.append(len(moleq)-1)
        org_eq.append(ind)
    else:
        inv.append(org_eq[ind])
        
    gp = 'Group(['
    for sp in sop:
        if np.linalg.det(sp.rotation_matrix) > 0:
            gp = gp + 'PermList('+str([i+1 for i in perm(mol.cart_coords,sp.operate_multi(mol.cart_coords))])+'),'
    gp = gp + '])'
    G.add_node(ind,group=gp,org_eq=org_eq[ind])

print("ur computation finished!!!")

# computing urt
urt = []
nts = len(molts)
org_ts = list(range(nts)) 
for ind, mol in enumerate(molts):
    molsym = PointGroupAnalyzer(mol,tolerance=tol)
    sop = molsym.get_symmetry_operations()
    if ind < nts and wchiral(sop):
        molts.append(invmol(mol))
        molrp.append([invmol(molrp[ind][0]),invmol(molrp[ind][1])])
        gconns.append(gconns[ind])
        org_ts.append(ind)

    gp = 'Group(['
    for sp in sop:
        if np.linalg.det(sp.rotation_matrix) > 0:
            gp = gp + 'PermList('+str([i+1 for i in perm(mol.cart_coords,sp.operate_multi(mol.cart_coords))])+'),'
    gp = gp + '])'
    urt.append(gp)

print("urt computation finished!!!")

# computing connections
for i in range(len(gconns)):
    s = []
    conni = []
    for j in [0,1]:
        hom = HungarianOrderMatcher(molrp[i][j])
        tran1 = hom.match(moleq[gconns[i][j]])
        tran2 = hom.match(moleq[inv[gconns[i][j]]])

        if tran1[3] > tol and tran2[3] > tol:
            gom = GeneticOrderMatcher(molrp[i][j],tol)
            tran1 = gom.match(moleq[gconns[i][j]])
            tran2 = gom.match(moleq[inv[gconns[i][j]]])

            if len(tran1) > len(tran2):
                conni.append(inv[gconns[i][j]])
                tran = tran1[0]
            else:
                conni.append(gconns[i][j])
                tran = tran2[0]
        else:
            if tran1[3] > tran2[3]:
                conni.append(inv[gconns[i][j]])
                tran = tran2
            else:
                conni.append(gconns[i][j])
                tran = tran1

        s.append([x+1 for x in tran[0]])

    G.add_edge(conni[0],conni[1],group=urt[i],org=org_ts[i],perms=[[conni[0],s[0]],[conni[1],s[1]]])

print("connections computation finished!!!")

# computing symc (in this implementation, we take the reference to 0) 
cycles = cycles(G,0)
symc = []
for c in cycles:
    gp = '()'
    for e in c:
        if e[2]['perms'][0][0] == e[0] and e[2]['perms'][1][0] == e[1]:
            gp = gp + '*Inverse(PermList('+str(e[2]['perms'][0][1])+'))'
            gp = gp + '*PermList('+str(e[2]['perms'][1][1])+')'
        elif e[2]['perms'][0][0] == e[1] and e[2]['perms'][1][0] == e[0]:
            gp = gp + '*Inverse(PermList('+str(e[2]['perms'][1][1])+'))'
            gp = gp + '*PermList('+str(e[2]['perms'][0][1])+')'
        else:
            # this case cannot happen.
            print('Error: something wrong happened in the function cycle',file=sys.stderr)
            sys.exit(1)
    symc.append(gp)

print("symc computation finished!!!")

fglist = args[4]
with open(fglist,mode='w') as f:
    # computing sym
    ispe = {}
    for sp in moleq[0].species:
        ispe[sp] = []
    for i, sp in enumerate(moleq[0].species):
        ispe[sp].append(i+1)
    f.write('sym:=Group([')
    for v in ispe.values():
       if len(v) >= 2:
           f.write('('+str(v[0])+','+str(v[1])+'),')
           f.write(str(tuple(v))+',')
    f.write(']);\n')

    f.write('symc:=Group([')
    for gp in symc:
        f.write(gp+',')
    # f.write('], GeneratorsOfGroup('+G.nodes[0]['group']+')));\n')
    f.write(']);\n')

    f.write('ur:=[')
    for n in G.nodes(data=True):
        f.write(n[1]['group']+',')
    f.write('];\n')

    f.write('urt:=[')
    for e in G.edges(data=True):
        f.write(e[2]['group']+',')
    f.write('];\n')

    f.write('ss:=[')
    for e in G.edges(data=True):
        f.write('[['+str(e[2]['perms'][0][0]+1)+','+'PermList('+str(e[2]['perms'][0][1])+')],')
        f.write('['+str(e[2]['perms'][1][0]+1)+','+'PermList('+str(e[2]['perms'][1][1])+')]],')
    f.write('];\n')

    f.write('org_eq:='+str([x+1 for x in org_eq])+';\n')
    f.write('org_ts:='+str([x+1 for x in org_ts])+';\n')

quit()
