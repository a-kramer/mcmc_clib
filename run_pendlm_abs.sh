#!/bin/bash

# start sampler with
# -b                        binary output 
# -l ODEmodel11S26P4U.so    shared library of model
# -c ODEmodel11S26P4U.cfg   configurations
ModelCFG=pendlm_abs
ModelSO=pendulum
ModelXML=${ModelCFG}

SampleFile="sample/${ModelCFG}_`date +%Y-%m-%dT%Hh%Mm`.double"
SampleSize=$((2**16))
P=`grep '<Parameter ' ${ModelXML}.xml | wc -l`
U=`grep 'Input' ${ModelXML}.xml | wc -l`

cat<<EOF
$0
 redirecting standard output to $ModelCFG.out
 redirecting standard error  to $ModelCFG.err

 output will be binary (${SampleFile})
 to load the sample in matlab or GNU Octave:

 NumOfParameters = $((P-U));
      SampleSize = ${SampleSize};
             fid = fopen('${SampleFile}','r');
          Sample = fread(fid,[NumOfParameters+1,SampleSize],'double');
                   fclose(fid);

 use «tailf $ModelCFG.out» to see progress (press «Ctrl-C» to exit tailf)
 redirecting standard output to $ModelCFG.out
             standard error  to $ModelCFG.err
 sampling now ...
EOF

bin/ode_smmala -a 0.25 -b -s ${SampleSize} -w $((SampleSize/4)) -o $SampleFile -l ./$ModelSO.so -c ./$ModelCFG.cfg 1> $ModelCFG.out 2> $ModelCFG.err
