#!/bin/bash

# start sampler with
# -b                        binary output 
# -l ODEmodel.so    shared library of model
# -c ODEmodel.cfg   configurations

# a simulation for model A may be run using different configurations B₁,… 
# prefix definitions
Model=CaMKII
ModelCFG=${Model}.cfg
ModelSO=${Model}.so
ModelXML=${Model}.xml
# number of sampling parameters:
SampleSize=$((2**19))
P=`grep '<Parameter ' ${ModelXML} | wc -l`
U=`egrep 'Description=".*Input' ${ModelXML} | wc -l`

echo "Model: ${ModelXML}; ${P} parameters total, ${U} of them are inputs."

ModelOPT="-l ./${ModelSO} -c ./${ModelCFG}"

cat<<EOF
$0
 output will be binary
 to load the sample in matlab or GNU Octave:

 NumOfParameters = $((P-U));
      SampleSize = ${SampleSize};
             fid = fopen('${SampleFile}','r');
          Sample = fread(fid,[NumOfParameters+1,SampleSize],'double');
                   fclose(fid);

 use «tailf $ModelCFG_….out» to see progress (press «Ctrl-C» to exit tailf)

 sampling now ...
EOF

for i in {1..8}; do
    SEED=$((1337+i))
    OPT="--seed $SEED --prior-start -b"
    OUT=${Model}_${i}.out
    ERR=${Model}_${i}.err
    SampleFile="sample/S${i}_seed${SEED}_${Model}_`date +%Y-%m-%dT%Hh%Mm`.double"
    echo -en "output: «tailf ${OUT}»\t"
    echo "«tailf ${ERR}» (stderr)"
    bin/ode_smmala ${OPT} -o ${SampleFile} -s ${SampleSize} -w $((SampleSize/4)) ${ModelOPT} 1> $OUT 2> $ERR &
done
