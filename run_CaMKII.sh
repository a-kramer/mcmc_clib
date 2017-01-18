#!/bin/bash

# start sampler with
# -b                        binary output 
# -l ODEmodel.so    shared library of model
# -c ODEmodel.cfg   configurations
ModelCFG=CaMKII
ModelSO=CaMKII
P=54
SampleFile="sample/${ModelCFG}_`date +%Y-%m-%dT%Hh%Mm`.double"
SampleSize=10000

cat<<EOF
$0
 redirecting standard output to $ModelCFG.out
 redirecting standard error  to $ModelCFG.err

 output will be binary (${SampleFile})
 to load the sample in matlab or GNU Octave:

 NumOfParameters = $P;
      SampleSize = ${SampleSize};
             fid = fopen('${SampleFile}','r');
          Sample = fread(fid,[NumOfParameters+1,SampleSize],'double');
                   fclose(fid);

 use «tailf $ModelCFG.out» to see progress (press «Ctrl-C» to exit tailf)

 sampling now ...
EOF

bin/ode_smmala  -b -o $SampleFile -s ${SampleSize} -w ${SampleSize} -l ./${ModelSO}.so -c ./${ModelCFG}.cfg 1> ${ModelCFG}.out 2> ${ModelCFG}.err
