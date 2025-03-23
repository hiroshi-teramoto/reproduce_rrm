#! /bin/sh

# Graph node type
IS_POINT=true

# Specify the directory name 
MOL="AuCu4"
DIR="../Metal/Metal/${MOL}/"
NAME="${MOL}_AFIR"

VFILE="vertices_${NAME}.dat"
EFILE="edges_${NAME}.dat"
OFILE="rrm_${NAME}.dot"
GFILE="rrm_${NAME}.png"
CFILE="rrm_${NAME}.cli"

# Specify the following log files 
EQLIST="${DIR}/`ls ${DIR} | grep EQ_list.log`"
TSLIST="${DIR}/`ls ${DIR} | grep TS_list.log`"
TSHEAD="${TSLIST%_*}"

echo ${EQLIST}
echo ${TSLIST}
echo ${TSHEAD}

# output file to feed it into GAP
GLIST="data/${NAME}.g"

python3 rrm_reconstruction_v18.py $EQLIST $TSLIST $TSHEAD $GLIST 

/usr/local/gap-4.13.0/gap -b -q -m 12g generate_rrm_v11.g << EOI
Read("$GLIST");;
vlabel:=true;;
elabel:=true;;

Print("the index of symc in sym : ",IndexNC(sym,symc),"\n");
generate_rrm("$VFILE","$EFILE",symc,ur,urt,ss,org_eq,org_ts,vlabel,elabel);
QUIT;
EOI

echo "GAP computation done."

echo "graph $NAME {" > $OFILE
echo "        outputorder=\"edgesfirst\"" >> $OFILE
echo "        overlap=false" >> $OFILE
echo "        spline=true" >> $OFILE
echo "        frontname=\"Helvetica,Arial,sans-serif\"" >> $OFILE

if "$IS_POINT" ; then
    echo "        node [shape=point,fontname=\"Helvetica,Arial,sans-serif\",color=blue]" >> $OFILE
else
    echo "        node [fontname=\"Helvetica,Arial,sans-serif\"]" >> $OFILE
fi

echo "        edge [fontname=\"Helvetica,Arial,sans-serif\",color=lightgray]" >> $OFILE

cat $VFILE >> $OFILE
cat $EFILE >> $OFILE

echo "}" >> $OFILE

dot -Kfdp -Tpng $OFILE -o $GFILE
python3 check_number_of_edges_v3.py $OFILE 

exit 0

python3 enumerate_cliques_gap.py $OFILE > $CFILE
# python3 enumerate_cliques_gap.py $OFILE

# finding orbits of cliques
/usr/local/gap-4.13.0/gap -b -q -m 12g << EOI

Read("$GLIST");;
Read("$CFILE");
inv := Group((1,2));
cnpi := DirectProduct(sym,inv);
cnpic := DirectProduct(symc,inv);

Print("sinv:");
sinv;

OnEqRight:=function( eqright, g )
	local g1, g2, ind;
	g1 := Image(Projection(cnpic,1),g);
	g2 := Image(Projection(cnpic,2),g);
	if g2 = () then
		return [eqright[1], OnRight(eqright[2],g1)];
	else
		ind := Int(ReplacedString(eqright[1],"*",""))+1;
		if wchiral[ind] then
			if '*' in eqright[1] then
				return [ReplacedString(eqright[1],"*",""), OnRight(eqright[2],g1)];
			else
				return [Concatenation(eqright[1],"*"), OnRight(eqright[2],g1)];
			fi;
		else
			return [eqright[1], RightCoset(ActingDomain(eqright[2]),OnLeftInverse(Representative(OnRight(eqright[2],g1)),sinv[ind]))];
		fi;
	fi;
end;

OnEqRightFull:=function( eqright, g )
	local g1, g2, ind;
	g1 := Image(Projection(cnpi,1),g);
	g2 := Image(Projection(cnpi,2),g);
	if g2 = () then
		return [eqright[1], OnRight(eqright[2],g1)];
	else
		ind := Int(ReplacedString(eqright[1],"*",""))+1;
		if wchiral[ind] then
			if '*' in eqright[1] then
				return [ReplacedString(eqright[1],"*",""), OnRight(eqright[2],g1)];
			else
				return [Concatenation(eqright[1],"*"), OnRight(eqright[2],g1)];
			fi;
		else
			return [eqright[1], RightCoset(ActingDomain(eqright[2]),OnLeftInverse(Representative(OnRight(eqright[2],g1)),sinv[ind]))];
		fi;
	fi;
end;

eqcosets:=[];
for i in [1..Length(vx)] do
   uri:=ur[Int(ReplacedString(vx[i][1],"*",""))+1];
   eqcosets[i]:=[vx[i][1],RightCoset(uri,CanonicalRightCosetElement(uri,vx[i][2]^-1))];
od; 

hom:=ActionHomomorphism(cnpic,eqcosets,OnEqRight);
op:=Image(hom);

reps:=[];;
srep:=[];;
basis:=[];;

for i in [1..Length(cl[Length(cl)][1])] do
	Append(basis,[[]]);
od; 

for cli in cl do
	d := Length(cli[1]);
	if not IsTransitive(op,cli,OnSets) then
		orbs:=Orbits(op,cli,OnSets);
		for orb in orbs do
			Append(basis[d],[orb]);
			Append(reps,[orb[1]]);
			Append(srep,[Length(orb)]);
		od;
	else
		Append(basis[d],[cli]);
		Append(reps,[cli[1]]);
		Append(srep,[Length(cli)]);
	fi;
od;

Print("\n");
Print("Dimension of Simplicial Complexes in Each Dimension\n");
dims := List( basis, c -> Length(c) );
dims;

Print("\n");
Print("List of Representatives of Orbits\n");
for i in [1..Length(reps)] do
	Print("Orbit Representative: ",reps[i]); 
	Print("\n");
	Print("Orbit Length: ",srep[i]);
	Print("\n");
	stab:=Stabilizer(op,reps[i],OnSets);
	pg:=PreImages(hom,stab);
	Print("Generators of Point Group for corresponding saddle r:\n");
	for g in pg do
		g1 := Image(Projection(cnpic,1),g);
		g2 := Image(Projection(cnpic,2),g);
		if g2 = () then
			Print(g1,"*e\n");
		else
			Print(g1,"*i\n");
		fi;
	od;
	Print("\n");
od;
Print();

Print("Boundary Operators Represented by the Basis\n");
bops := [];;
for d in [1..Length(dims)-1] do
	bop := NullMat(dims[d],dims[d+1],GF(2));
	# compute boundary operator here !
	for b in basis[d+1] do
		bbp := [];
		bbm := [];
		for s in b do
			for i in [1..Length(s)] do
				bs := ShallowCopy(s);
				RemoveSet( bs, s[i] );
				if IsOddInt(i) then
					Append(bbp,[bs]);
				else
					Append(bbm,[bs]);
				fi;
			od;
		od;

		for bd in basis[d] do
			bop[Position(basis[d],bd)][Position(basis[d+1],b)] := (Number(bbp, x -> x = bd[1]) + Number(bbm, x -> x = bd[1])) mod 2;
		od;
	od;
	Append(bops,[bop]);
od;

Print("Boundary Operators:\n");
for d in [1..Length(bops)] do
	Print(d+1, "--->", d, "\n");
	Print(bops[d], "\n");
od;

Print("Dimensions of Homology:\n");
for i in [1..Length(dims)] do
	if i = 1 then
		Print("0-th Homology:", dims[1]-RankMatrix(bops[1]), "\n");
	elif i < Length(dims) then
		Print(i-1, "-th Homology:", dims[i]-RankMatrix(bops[i-1])-RankMatrix(bops[i]), "\n");
	else
		Print(i-1, "-th Homology:", dims[i]-RankMatrix(bops[i-1]), "\n");
	fi;
od;

EOI

echo "DONE"
