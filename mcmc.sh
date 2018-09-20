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
    fi    
}

zenity --info --text="An MCMC run will require a model (an <tt>.so</tt> file or a VFGEN file), a file with prior and data (specific hdf5 file or collection of <tt>.tsv</tt> files)"

Model=`zenity --file-selection --title="Select a Model File" --file-filter="*.so *.xml *.vf"`
EC=$?

if (($EC!=0)); then
    
elif

## this block handles the model file of the mcmc run
if [[ $Model =~ .*[.](xml|vf) ]]; then
    NPvf=`grep '<Parameter' $Model | wc -l`
    ModelName=`grep '<VectorField' | sed -e 's/^.*Name="([^"]+)".*$/\1/'`
    VFGEN=`which vfgen`
    if [[ -z $VFGEN ]]; then
	zenity --info --text="Could not find <tt>vfgen</tt> to convert the <tt>vf</tt> model $Model to a shared object (<tt>.so</tt>)"
	exit -2;
    else
	vfgen_out=`vfgen cvodes:func=yes,sens=yes $Model`
	EC=$?;
	if ((EC==0)); then
	    gcc -shared -fPIC -O2 -march=native -lm -o ${ModelName}.so ${ModelName}_cvs.c
	else
	    zenity --error --text="vfgen error."
	    echo "$vfgen_out" | zenity --text-info --title="Error Text"
	    exit -2;
	fi
    fi
fi

## this block will handle the data part of the run

zenity
