# Synthetic RRM input -- NOT a GRRM output and not physically meaningful.
# Generator in the development fork (commit 570ad25):
# https://github.com/mochimochi-brian/reproduce_rrm/blob/570ad25/tools/make_synthetic_input.py
# Regenerate with that script:
#   python3 tools/make_synthetic_input.py --identical 6 --eqs 8 --ts 30 --seed 1 --large-eq-points 3 --large-eqs 1 --transposition-period 4 --light-period 3 --same-eq-period 5 -o <output path>
# sym = symc = S_6 on points 2..7 (no inversion extension), 8 EQs, 30 TSs
# light TSs (|sym|/|ur[1]| edges): 10, heavy TSs (|sym| edges): 20
sym:=SymmetricGroup([2..7]);;
symc:=sym;;
ur:=[Group([PermList([1,3,2,4,5,6,7]),PermList([1,3,4,2,5,6,7])]),Group([()]),Group([()]),Group([PermList([1,3,2,4,5,6,7])]),Group([()]),Group([()]),Group([()]),Group([PermList([1,3,2,4,5,6,7])])];;
ss:=[[[5,PermList([1,5,2,4,6,3,7])],[4,PermList([1,6,4,5,2,7,3])]],[[8,PermList([1,4,2,7,3,6,5])],[2,PermList([1,5,7,2,4,6,3])]],[[1,PermList([1,2,3,4,5,6,7])],[1,PermList([1,2,3,4,5,6,7])]],[[4,PermList([1,6,3,7,5,4,2])],[1,PermList([1,6,4,3,2,7,5])]],[[5,PermList([1,2,4,6,5,3,7])],[5,PermList([1,2,3,5,7,4,6])]],[[1,PermList([1,2,3,4,5,6,7])],[1,PermList([1,2,3,4,5,6,7])]],[[8,PermList([1,2,3,4,5,7,6])],[3,PermList([1,7,2,6,5,3,4])]],[[5,PermList([1,3,4,6,5,7,2])],[8,PermList([1,2,7,4,6,5,3])]],[[1,PermList([1,2,3,4,5,6,7])],[1,PermList([1,2,3,4,5,6,7])]],[[6,PermList([1,3,7,2,5,6,4])],[6,PermList([1,6,7,2,5,4,3])]],[[4,PermList([1,6,7,2,5,4,3])],[8,PermList([1,5,4,2,3,7,6])]],[[1,PermList([1,2,3,4,5,6,7])],[1,PermList([1,2,3,4,5,6,7])]],[[3,PermList([1,7,4,5,6,2,3])],[1,PermList([1,4,6,7,2,5,3])]],[[6,PermList([1,2,5,4,6,7,3])],[4,PermList([1,6,5,3,4,2,7])]],[[1,PermList([1,2,3,4,5,6,7])],[1,PermList([1,2,3,4,5,6,7])]],[[3,PermList([1,5,4,2,7,3,6])],[2,PermList([1,6,2,3,5,7,4])]],[[7,PermList([1,6,5,3,2,4,7])],[4,PermList([1,5,7,2,4,6,3])]],[[1,PermList([1,2,3,4,5,6,7])],[1,PermList([1,2,3,4,5,6,7])]],[[7,PermList([1,7,3,4,2,5,6])],[5,PermList([1,5,3,7,2,4,6])]],[[1,PermList([1,4,7,5,6,3,2])],[1,PermList([1,5,3,7,6,4,2])]],[[1,PermList([1,2,3,4,5,6,7])],[1,PermList([1,2,3,4,5,6,7])]],[[5,PermList([1,6,5,4,7,3,2])],[4,PermList([1,3,5,4,2,7,6])]],[[5,PermList([1,6,5,2,3,4,7])],[8,PermList([1,6,7,5,2,4,3])]],[[1,PermList([1,2,3,4,5,6,7])],[1,PermList([1,2,3,4,5,6,7])]],[[3,PermList([1,3,7,6,4,5,2])],[3,PermList([1,7,2,3,4,5,6])]],[[2,PermList([1,7,6,4,5,3,2])],[3,PermList([1,3,2,5,7,6,4])]],[[1,PermList([1,2,3,4,5,6,7])],[1,PermList([1,2,3,4,5,6,7])]],[[1,PermList([1,2,5,3,6,7,4])],[5,PermList([1,3,2,5,6,4,7])]],[[6,PermList([1,5,4,2,3,6,7])],[2,PermList([1,6,3,5,4,2,7])]],[[1,PermList([1,2,3,4,5,6,7])],[1,PermList([1,2,3,4,5,6,7])]]];;
urt:=List([1..Length(ss)], i -> Intersection(
        ConjugateGroup(ur[ss[i][1][1]], ss[i][1][2]),
        ConjugateGroup(ur[ss[i][1][1]], ss[i][1][2]^-1),
        ConjugateGroup(ur[ss[i][2][1]], ss[i][2][2]),
        ConjugateGroup(ur[ss[i][2][1]], ss[i][2][2]^-1)));;
org_eq:=[1..Length(ur)];;
org_ts:=[1..Length(ss)];;
