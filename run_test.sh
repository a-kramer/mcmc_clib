#!/bin/bash

# start sampler with
# -b                        binary output 
# -l ODEmodel11S26P4U.so    shared library of model
# -c ODEmodel11S26P4U.cfg   configurations
Model=ODEmodel11S26P4U
SampleFile="sample/${Model}_`date +%Y-%m-%dT%Hh%Mm`.double"
#SampleSize=$((1024**2))
SampleSize=$((2))
cat<<EOF
$0
 redirecting standard output to $Model.out
 redirecting standard error  to $Model.err

 output will be binary (${SampleFile})
 to load the sample in matlab or GNU Octave:

 NumOfParameters = 26;
      SampleSize = ${SampleSize};
             fid = fopen('${SampleFile}','r');
          Sample = fread(fid,[NumOfParameters+1,SampleSize],'double');
                   fclose(fid);

 use «tailf $Model.out» to see progress (press «Ctrl-C» to exit tailf)
 redirecting standard output to $Model.out
             standard error  to $Model.err
 sampling now ...
EOF

mpirun -np 2 bin/ode_smmala  -b -s ${SampleSize} -p -o $SampleFile -l ./$Model.so -c ./$Model.cfg 1> $Model.out 2> $Model.err
