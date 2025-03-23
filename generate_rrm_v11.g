# vfile : file to save vertex information
# efile : file to save edge information
# sym : sym_Omega
# ur, urt : list of U(r_R)s and U(r_T)s
# ss : list of [[i,s1],[j,s2]] (i, j: indices of equilibrium structures, s1 and s2 are permutation from the equilibrium structures)  
# labels : vertex (edge) label is attached if true and is not otherwise
generate_rrm:=function(vfile, efile, sym, ur, urt, ss, org_eq, org_ts, labels...) 
	local rturi, rturti, pUri, vertices, v1, v2, edges, ind, p, vlabel, elabel, i, cosets, actionHom, permGroup, pos1, pos2, stabPerm, stabEQs;

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

    	vertices:=[];
        for i in [1..Length(ur)] do
            for rturi in RightTransversal(sym,ur[i]) do
                Append(vertices,[[i,CanonicalRightCosetElement(ur[i],rturi)]]); 
            od;
        od;

	# check if GRRM graph is consistent with the required symmetry.
	for i in [1..Length(ss)] do
		# In this case, the reactant and product are permutation isomers.
		if ss[i][1][1] = ss[i][2][1] then
			cosets := RightCosets(sym,ur[ss[i][1][1]]);
			actionHom := ActionHomomorphism(sym,cosets,OnRight);	
			permGroup := Image(actionHom);
			pos1 := Position(cosets,ur[ss[i][1][1]]*(ss[i][1][2]^-1));
			pos2 := Position(cosets,ur[ss[i][2][1]]*(ss[i][2][2]^-1));
			stabPerm := Stabilizer(permGroup,Set([pos1,pos2]),OnSets);
			stabEQs := PreImage(actionHom,stabPerm);
			if not IsSubgroup(stabEQs,urt[i]) then
				Assert(0,false,"Violation of Pechukus theorem: \n");
	        		if org_ts[i] = i then
					Print("TS",i-1,":\n");
				else
					Print("TS",i-1,"*:\n");
				fi;

				Print("the reactant and product are permutation isomers\n");
                        	if org_eq[ss[i][1][1]]=ss[i][1][1] then
					Print("EQ",org_eq[ss[i][1][1]]-1,"\n");
				else
					Print("EQ",org_eq[ss[i][1][1]]-1,"*\n");
				fi;

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
	        		if org_ts[i] = i then
					Print("TS",i-1,":\n");
				else
					Print("TS",i-1,"*:\n");
				fi;

				Print("reactant:\n");
                        	if org_eq[ss[i][1][1]]=ss[i][1][1] then
					Print("EQ",org_eq[ss[i][1][1]]-1,"\n");
				else
					Print("EQ",org_eq[ss[i][1][1]]-1,"*\n");
				fi;

				Print("product:\n");
                        	if org_eq[ss[i][2][1]]=ss[i][2][1] then
					Print("EQ",org_eq[ss[i][2][1]]-1,"\n");
				else
					Print("EQ",org_eq[ss[i][2][1]]-1,"*\n");
				fi;

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

	edges:=[];
        for i in [1..Length(ss)] do
		for rturti in RightTransversal(sym,urt[i]) do
			v1:=[ss[i][1][1],CanonicalRightCosetElement(ur[ss[i][1][1]],ss[i][1][2]^-1*rturti)];
			v2:=[ss[i][2][1],CanonicalRightCosetElement(ur[ss[i][2][1]],ss[i][2][2]^-1*rturti)];
			Append(edges,[[[i,CanonicalRightCosetElement(urt[i],rturti)],[Position(vertices,v1),Position(vertices,v2)]]]);
		od;
        od;

        PrintTo(vfile);
	ind:=0;
	for p in vertices do
		if ind <> p[1] then 
			if ind <> 0 then
				AppendTo(vfile,"}\n");
			fi;

                        if org_eq[p[1]]=p[1] then
				AppendTo(vfile,"subgraph cluster_",p[1]-1," { label=\"",org_eq[p[1]]-1,"\";\n");
			else
				AppendTo(vfile,"subgraph cluster_",p[1]-1," { label=\"",org_eq[p[1]]-1,"*\";\n");
			fi;

			AppendTo(vfile,"fontsize=\"30pt\"\n");
		fi;
		ind:=p[1];

		if vlabel then
                        if org_eq[p[1]]=p[1] then
    				AppendTo(vfile,Position(vertices,p),"[label=\"",org_eq[p[1]]-1," ",p[2]^-1,"\"]\n");
			else
    				AppendTo(vfile,Position(vertices,p),"[label=\"",org_eq[p[1]]-1,"* ",p[2]^-1,"\"]\n");
			fi;
		else
    			AppendTo(vfile,Position(vertices,p),"\n");
		fi;
	od;
	AppendTo(vfile,"}\n");

	PrintTo(efile);
	for p in edges do
		if elabel then
	        	if org_ts[p[1][1]]=p[1][1] then
				AppendTo(efile,p[2][1],"--",p[2][2],"[label=\"",org_ts[p[1][1]]-1," ",p[1][2]^-1,"\"]\n");
			else
				AppendTo(efile,p[2][1],"--",p[2][2],"[label=\"",org_ts[p[1][1]]-1,"* ",p[1][2]^-1,"\"]\n");
			fi;
		else	
			AppendTo(efile,p[2][1],"--",p[2][2],"\n");
		fi;
	od;
	return;
end;;
