#! /bin/sh

# Graph node type
IS_POINT=true

# Specify the directory name 
MOL="pentane"
DIR="${MOL}_Restruct"
NAME="Re_${MOL}_AFIR"

# output files
VFILE="vertices_${NAME}.dat"
EFILE="edges_${NAME}.dat"
# resulting rrm in shape space in the dot format
OFILE="rrm_${NAME}.dot"
# png file in which the resulting rrm is visualized
GFILE="rrm_${NAME}.png"

# Specify the following log files 
EQLIST="${DIR}/`ls ${DIR} | grep EQ_list.log`"
TSLIST="${DIR}/`ls ${DIR} | grep TS_list.log`"
TSHEAD="${TSLIST%_*}"

echo ${EQLIST}
echo ${TSLIST}
echo ${TSHEAD}

# output file to feed it into GAP
GLIST="data/${NAME}.g"

python rrm_reconstruction_v12.py $EQLIST $TSLIST $TSHEAD $GLIST 

/usr/local/gap-4.11.1/bin/gap.sh -b -q -m 12g generate_rrm_v9.g << EOI
Read("$GLIST");;
vlabel:=false;;
elabel:=false;;

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

# use fdp if the resulting graph is too big
# dot -Kdot -Tpng $OFILE -o $GFILE
dot -Kfdp -Tpng $OFILE -o $GFILE

echo "DONE"
