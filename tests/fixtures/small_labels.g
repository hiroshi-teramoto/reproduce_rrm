# Small labeled-graph fixture for byte comparison with v11.
# All transition-state groups are trivial; there are no Pechukas violations.
#
# Coverage:
#   EQ1  ur = <(1,2)>   -> 3 vertices, label "0"
#   EQ2  ur = 1         -> 6 vertices, label "1"
#   EQ3  ur = 1         -> 6 vertices, inversion isomer of EQ2, label "1*"
#   TS1  same EQ, s1 = s2            -> self-loops on EQ1
#   TS2  same EQ, s1 <> s2           -> EQ2 internal edges, multiplicity 2
#   TS3  EQ1 -- EQ2
#   TS4  identical connection to TS3 -> parallel edges with a distinct TS label,
#        and org_ts[4] = 3 makes it an inversion isomer, label "2*"
#   TS5  EQ3 -- EQ1, exercising the starred EQ endpoint
# 15 vertices, 5 x |S_3| = 30 edges.
sym:=Group([(1,2),(1,2,3)]);
symc:=sym;
ur:=[Group([(1,2)]),Group([()]),Group([()])];
urt:=[Group([()]),Group([()]),Group([()]),Group([()]),Group([()])];
ss:=[
  [[1,()],[1,()]],
  [[2,()],[2,(1,2)]],
  [[1,()],[2,()]],
  [[1,()],[2,()]],
  [[3,(1,2,3)],[1,()]]
];
org_eq:=[1,2,2];
org_ts:=[1,2,3,3,5];
