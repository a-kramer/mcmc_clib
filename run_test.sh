#!/bin/bash

# start sampler with
# -b                        binary output 
# -l ODEmodel11S26P4U.so    shared library of model
# -c ODEmodel11S26P4U.cfg   configurations
Model=ODEmodel11S26P4U
SampleFile="`date +%Y-%m-%dT%Hh%Mm`.h5"
#SampleSize=$((1024**2))
DefaultSampleSize=$((2**10))
SampleSize=${2:-$DefaultSampleSize}
NP=${1:-2}
cat<<EOF
$0
 redirecting standard output to $Model.out
 redirecting standard error  to $Model.err

 output will be binary (${SampleFile})
 to load the sample in matlab or GNU Octave:

 NumOfParameters = 26;
      SampleSize = ${SampleSize};
          Sample = cell($NP,1);
  for i=1:$NP
             fid = fopen(sprintf("rank_%04i_of_$NP_${SampleFile}",i),'r');
       Sample{i} = fread(fid,[NumOfParameters+1,SampleSize],'double');
                   fclose(fid);
  endfor

 sampling now ...
EOF

mpirun -np $NP bin/ode_smmala -b -w ${SampleSize} -s ${SampleSize} -p -o $SampleFile -l ./$Model.so -c ./$Model.cfg 
