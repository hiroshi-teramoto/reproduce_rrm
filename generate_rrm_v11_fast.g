# Optional sequential speedup of generate_rrm_v11.g, with the same arguments.
# Vertex IDs retain RightTransversal order; a hash of canonical cosets replaces
# linear Position(vertices, ...) searches. Open streams avoid per-line file
# opens, and edges are written as generated instead of retaining the full list.
# The same-EQ Pechukas check reuses each EQ's transversal and action homomorphism.
# Pechukas diagnostics and continuation follow v11. No continue flag is needed.
# Invalid vertex lookups and failures to open output files abort with status 1.

# vfile : file to save vertex information
# efile : file to save edge information
# sym : sym_Omega
# ur, urt : list of U(r_R)s and U(r_T)s
# ss : list of [[i,s1],[j,s2]] (i, j: indices of equilibrium structures, s1 and s2 are permutation from the equilibrium structures)
# labels : vertex (edge) label is attached if true and is not otherwise

RrmNewCosetIndex:=function(sym)
	local h;
	h:=SparseIntKey(sym, One(sym));
	if h=fail then
		Print("Error, no permutation hash key for vertex index\n");
		FORCE_QUIT_GAP(1);
	fi;
	return SparseHashTable(h);
end;;

RrmVertexIndex:=function(offset, idx, ur, eq, g)
	local pos;
	pos:=LookupDictionary(idx[eq], CanonicalRightCosetElement(ur[eq], g));
	if pos=fail then
		# Error() returns under gap -T; FORCE_QUIT_GAP is the batch stop.
		Print("Error, vertex lookup failed for EQ",eq-1,"\n");
		FORCE_QUIT_GAP(1);
	fi;
	return offset[eq]+pos;
end;;

RrmBuildTransversals:=function(sym, ur)
	local rt, offset, nvert, idx, i, j, ht;
	rt:=[];
	offset:=[];
	idx:=[];
	nvert:=0;
	for i in [1..Length(ur)] do
		rt[i]:=RightTransversal(sym,ur[i]);
		ht:=RrmNewCosetIndex(sym);
		for j in [1..Length(rt[i])] do
			AddDictionary(ht, CanonicalRightCosetElement(ur[i], rt[i][j]), j);
		od;
		idx[i]:=ht;
		offset[i]:=nvert;
		nvert:=nvert+Length(rt[i]);
	od;
	return rec(rt:=rt, offset:=offset, idx:=idx);
end;;

RrmInversionSuffix:=function(original, index)
	if original<>index then return "*"; fi;
	return "";
end;;

RrmCheckPechukas:=function(sym, ur, urt, ss, org_eq, org_ts, rt)
	local homByEq, i, eq, actionHom, permGroup, pos1, pos2, stabPerm, stabEQs;
	homByEq:=[];
	# check if GRRM graph is consistent with the required symmetry.
	for i in [1..Length(ss)] do
		# In this case, the reactant and product are permutation isomers.
		if ss[i][1][1] = ss[i][2][1] then
			eq:=ss[i][1][1];
			if not IsBound(homByEq[eq]) then
				homByEq[eq]:=ActionHomomorphism(sym,rt[eq],OnRight,"surjective");
			fi;
			actionHom:=homByEq[eq];
			permGroup:=Image(actionHom);
			pos1:=PositionCanonical(rt[eq],ss[i][1][2]^-1);
			pos2:=PositionCanonical(rt[eq],ss[i][2][2]^-1);
			if pos1=fail or pos2=fail then
				Print("Error, Pechukas coset lookup failed for TS",i-1,"\n");
				FORCE_QUIT_GAP(1);
			fi;
			stabPerm:=Stabilizer(permGroup,Set([pos1,pos2]),OnSets);
			stabEQs:=PreImage(actionHom,stabPerm);
			if not IsSubgroup(stabEQs,urt[i]) then
				Assert(0,false,"Violation of Pechukus theorem: \n");
				Print("TS",i-1,RrmInversionSuffix(org_ts[i],i),":\n");

				Print("the reactant and product are permutation isomers\n");
				Print("EQ",org_eq[ss[i][1][1]]-1,RrmInversionSuffix(org_eq[ss[i][1][1]],ss[i][1][1]),"\n");

				Print("in what follows, minus 1 to convert to the atom labels\n");
				Print("U(r) (for reactant):\n");
				Print(ur[ss[i][1][1]],"\n");
				Print("U(r) (for product):\n");
				Print(ur[ss[i][2][1]],"\n");
				Print("stabEQs:\n");
				Print(stabEQs,"\n");
				Print("U(r) (for transition state):\n");
				Print(urt[i],"\n");
			fi;
		# In this case, the reactant and product are not.
		else
			if not IsSubgroup(ConjugateGroup(ur[ss[i][1][1]],ss[i][1][2]),urt[i]) or
				not IsSubgroup(ConjugateGroup(ur[ss[i][2][1]],ss[i][2][2]),urt[i]) then

				Assert(0,false,"Violation of Pechukus theorem: \n");
				Print("TS",i-1,RrmInversionSuffix(org_ts[i],i),":\n");

				Print("reactant:\n");
				Print("EQ",org_eq[ss[i][1][1]]-1,RrmInversionSuffix(org_eq[ss[i][1][1]],ss[i][1][1]),"\n");

				Print("product:\n");
				Print("EQ",org_eq[ss[i][2][1]]-1,RrmInversionSuffix(org_eq[ss[i][2][1]],ss[i][2][1]),"\n");

				Print("in what follows, minus 1 to convert to the atom labels\n");
				Print("U(r) (for reactant):\n");
				Print(ConjugateGroup(ur[ss[i][1][1]],ss[i][1][2]),"\n");
				Print("U(r) (for product):\n");
				Print(ConjugateGroup(ur[ss[i][2][1]],ss[i][2][2]),"\n");
				Print("U(r) (for transition state):\n");
				Print(urt[i],"\n");
			fi;
		fi;
	od;
end;;

RrmWriteVertices:=function(vfile, ur, org_eq, rt, vlabel)
	local vstream, vid, i, rturi, canon, suffix;
	vstream:=OutputTextFile(vfile,false);
	if vstream=fail then
		Print("Error, cannot write ",vfile,"\n");
		FORCE_QUIT_GAP(1);
	fi;
	SetPrintFormattingStatus(vstream,false);
	vid:=0;
	for i in [1..Length(ur)] do
		if i>1 then AppendTo(vstream,"}\n"); fi;
		suffix:=RrmInversionSuffix(org_eq[i],i);
		AppendTo(vstream,"subgraph cluster_",i-1," { label=\"",org_eq[i]-1,suffix,"\";\n");
		AppendTo(vstream,"fontsize=\"30pt\"\n");
		for rturi in rt[i] do
			vid:=vid+1;

			canon:=CanonicalRightCosetElement(ur[i],rturi);
			if vlabel then
				AppendTo(vstream,vid,"[label=\"",org_eq[i]-1,suffix," ",canon^-1,"\"]\n");
			else
				AppendTo(vstream,vid,"\n");
			fi;
		od;
	od;
	AppendTo(vstream,"}\n");
	CloseStream(vstream);
end;;

RrmWriteEdges:=function(efile, sym, urt, ss, org_ts, offset, idx, ur, elabel)
	local estream, i, rturti, id1, id2, canon, suffix;
	estream:=OutputTextFile(efile,false);
	if estream=fail then
		Print("Error, cannot write ",efile,"\n");
		FORCE_QUIT_GAP(1);
	fi;
	SetPrintFormattingStatus(estream,false);
	for i in [1..Length(ss)] do
		suffix:=RrmInversionSuffix(org_ts[i],i);
		for rturti in RightTransversal(sym,urt[i]) do
			id1:=RrmVertexIndex(offset,idx,ur,ss[i][1][1],ss[i][1][2]^-1*rturti);
			id2:=RrmVertexIndex(offset,idx,ur,ss[i][2][1],ss[i][2][2]^-1*rturti);
			canon:=CanonicalRightCosetElement(urt[i],rturti);
			if elabel then
				AppendTo(estream,id1,"--",id2,"[label=\"",org_ts[i]-1,suffix," ",canon^-1,"\"]\n");
			else
				AppendTo(estream,id1,"--",id2,"\n");
			fi;
		od;
	od;
	CloseStream(estream);
end;;

generate_rrm:=function(vfile, efile, sym, ur, urt, ss, org_eq, org_ts, labels...)
	local vlabel, elabel, built;

	if Length(labels)=0 then
		vlabel:=true;
		elabel:=true;
	elif Length(labels)=1 then
		vlabel:=labels[1];
		elabel:=true;
	else
		vlabel:=labels[1];
		elabel:=labels[2];
	fi;

	built:=RrmBuildTransversals(sym, ur);
	RrmCheckPechukas(sym, ur, urt, ss, org_eq, org_ts, built.rt);
	RrmWriteVertices(vfile, ur, org_eq, built.rt, vlabel);
	RrmWriteEdges(efile, sym, urt, ss, org_ts, built.offset, built.idx, ur, elabel);
	return;
end;;
