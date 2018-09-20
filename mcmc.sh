#!/bin/bash

ZENITY=`which zenity`
if [[ -z $ZENITY ]]; then
    echo "please install zenity"
    exit -1;
fi

function assert {
    if (($# < 2)); then
	echo "too few arguments to assert"
	exit -3;
    fi
    if (( $1 != 0 )); then
	echo "assertion failed, line $2."
	exit -3;
    fi
    return 0;
}

zenity --info --text="An MCMC run will require an ODE model suitable for integration by the CVODES solver (an <tt>.so</tt> file or a VFGEN file), a file with prior and data (specific hdf5 file or collection of <tt>.tsv</tt> files)"

Model=`zenity --file-selection --title="Select a Model File" --file-filter="*.so *.xml *.vf"`
assert $? $LINENO

## get the base name and strip the extension
ModelName=`basename $Model`
ModelName=${ModelName%.*}

## this block handles the model file of the mcmc run
if [[ $Model =~ .*[.](xml|vf) ]]; then
    ModelXML=$Model;
    NPvf=`grep '<Parameter' $Model | wc -l`
    ## with this we can determine the official model name
    ModelName=`grep '<VectorField' $ModelXML | sed -e 's/^.*Name="([^"]+)".*$/\1/'`
    VFGEN=`which vfgen`
    if [[ -z $VFGEN ]]; then
	zenity --info --text="Could not find <tt>vfgen</tt> to convert the <tt>vf</tt> model $Model to a shared object (<tt>.so</tt>)"
	exit -2;
    else
	vfgen_out=`vfgen cvodes:func=yes,sens=yes $Model`
	EC=$?;
	if ((EC==0)); then
	    ModelSO=${ModelName}.so	    
	    gcc -shared -fPIC -O2 -march=native -lm -o ${ModelSO} ${ModelName}_cvs.c
	else
	    zenity --error --text="vfgen error."
	    echo "$vfgen_out" | zenity --text-info --title="Error Text"
	    exit -2;
	fi
    fi
else
    ModelSO=$Model;
fi

## this block will handle the data part of the run

Data=`zenity --file-selection --title="Select a Data File" --multiple --file-filter="*.h5 *.tsv" --separator=" "`
assert $? $LINENO

if [[ "$Data" =~ \| ]]; then
    DataName=${ModelName:-"DataHDF5"}.h5
    sbtab $Data ${DataName}
    assert $? $LINENO
elif [[ -f "$Data" ]]; then
    DataName="$Data"
    ## one file
else
    zenity --error --text="Data file selection not understood"
fi

ResultSampleTag="`date +%Y-%m-%dT%Hh%Mm`.h5"

MCMCconf=( `zenity --forms --text="Parameters for this MCMC Run" --add-entry="Number of Processes" --add-entry="Sample Size" --add-entry="Burn in size" --add-entry="Temperature Schedule Factor [R⁺\0]" --separator=" "` )
assert $? $LINENO
declare -A Option
Option["NumProc"]=${MCMCconf[0]:-8};
Option["SampleSize"]=${MCMCconf[1]:-10240};
Option["WarmUp"]=${MCMCconf[2]:-1024};
Option["Gamma"]=${MCMCconf[3]:-"0.2"};

if [[ `which aprun` ]]; then
    MPICMD=aprun
else
    MPICMD=mpirun
fi

OutOpt="-b -o ${ResultSampleTag}"

ODEsmmalaCommand="$MPICMD -n ${Option['NumProc']} ode_smmala -s ${Option['SampleSize']} -w ${Option['WarmUp']} -d ${DataName} -l ${ModelSO} ${OutOpt}"
ODEsmmalaCommand=`echo "${ODEsmmalaCommand}" | zenity --text-info --editable --checkbox="That seems to be correct"`
echo ${ODEsmmalaCommand}
