#!/bin/bash




if [ "$#" -lt 1 ]; then

	echo "Must specify the FASTA file"
	echo "Exiting"
	exit 1

fi

OPTIND=1         # Reset in case getopts has been used previously in the shell.



RATIO=""
PREFIX=""


# Initialize our own variables:

while getopts "r:p:" opt; do
    case "$opt" in

    r)  RATIO=$OPTARG
        ;;
    p)  PREFIX=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

FASTA=$@


if [ "$RATIO" = "" ]; then

	RATIO=7
fi

if [ "$PREFIX" = "" ]; then

	PREFIX=$FASTA

fi

#echo "RATIO=$RATIO, PREFIX=$PREFIX, FASTA=$@"

sed -i "s,\#define OCC_INTERVAL.*,\#define OCC_INTERVAL 0x80,g" ./bwa/bwt.h 

echo "Compiling bwa"

make -C ./bwa clean
make -C ./bwa

echo "Building the suffix position array with compression ratio = $RATIO"

./bwa/bwa index -b sa -r $RATIO -p $PREFIX $FASTA

rm ${PREFIX}.bwt


echo "Building BWT array"

sed -i "s,\#define OCC_INTERVAL.*,\#define OCC_INTERVAL 0x40,g" ./bwa/bwt.h 

make -C ./bwa clean
make -C ./bwa

./bwa/bwa index -b bwt -p $PREFIX ${FASTA} 
