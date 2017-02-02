#!/bin/bash

# start sampler with
# -b                        binary output 
# -l ODEmodel11S26P4U.so    shared library of model
# -c ODEmodel11S26P4U.cfg   configurations
ModelCFG=ODEmodel_norm_ft
ModelSO=ODEmodel11S26P4U

SampleFile="sample/${ModelCFG}_`date +%Y-%m-%dT%Hh%Mm`.double"
SampleSize=10000

cat<<EOF
$0
 redirecting standard output to $ModelCFG.out
 redirecting standard error  to $ModelCFG.err

 output will be binary (${SampleFile})
 to load the sample in matlab or GNU Octave:

 NumOfParameters = 26;
      SampleSize = ${SampleSize};
             fid = fopen('${SampleFile}','r');
          Sample = fread(fid,[NumOfParameters+1,SampleSize],'double');
                   fclose(fid);

 use «tailf $ModelCFG.out» to see progress (press «Ctrl-c» to exit tailf)

 sampling now ...
EOF

#echo "run -b -o $SampleFile -s ${SampleSize} -l ./$ModelSO.so -c ./$ModelCFG.cfg"
bin/ode_smmala  -b -o $SampleFile -s ${SampleSize} -w ${SampleSize} -l ./$ModelSO.so -c ./$ModelCFG.cfg > $ModelCFG.out 
