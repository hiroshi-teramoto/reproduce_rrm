#! /bin/sh

# Graph node type
IS_POINT=true

# Specify the directory name 
# MOL="Au5Ag"
MOL="AuCu4"
DIR="Metal/${MOL}/"
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

echo "DONE"
