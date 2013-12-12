#!/bin/bash

if (($#!=3)); then
    echo 'Usage: diag $rows $cols $value'
    echo 'returns a block of numbers with $value in diagonal entries (row=column)'
else
    for ((i=0;i<$1;i++)); do
	for ((j=0;j<$2;j++)); do
	    if ((i==j)); then echo -en "$3\t"; else echo -en "0\t"; fi
	done
	echo
    done
fi
