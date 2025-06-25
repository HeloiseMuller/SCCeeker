#!/bin/bash

__usage="
Usage: $(basename $0) [OPTIONS]

Options:
  -d              	    absolute path to the database directory
  -o              	    absolute path to the output directory 
  -f           		    absolute path to a 2 tab deliminated colums file, where the first cololum contains the path of each fasta file and the second column contains the prefix to use for the output of the corresponding fasta. Do not use tilde
  -k 	                    'reference' or 'extended' depending on the database to use for the k-mer approach. No option skips the k-mer approach.
  -h, --help                Print USAGE; ignore all other parameters
  -v, --version             Print version number; ignore all other parameters
"

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  	echo "$__usage"
	exit 0
fi

if [[ "$1" == "-v" || "$1" == "--version" ]]; then
        echo "version v.0.1"
        exit 0
fi


#if [[ $# < 3 ||  $# > 4 ]]
#  then
#    echo "Wrong number of arguments: expect 3 or 4 arguments"
#    echo "$__usage"
#    exit 0
#fi

while getopts d:o:f:k: flag
do
    case "${flag}" in
        d) pathDB=${OPTARG};;
        o) pathOut=${OPTARG};;
        f) file=${OPTARG};;
        k) kmer=${OPTARG};;
    esac
done

echo "path of the database: $pathDB"
echo "path of the output: $pathOut"
echo "file where paths of the fasta are located: $file"

cd $(dirname $0) 

while read -r inp out
do
	echo ""
	bash SCCeeker.sh $pathDB $inp $pathOut $out $kmer

	#Add name prefix in output
#	awk -v var="$out" '{print $0, var}' OFS='\t' $pathOut/$out/${out}_blast.out > $pathOut/$out/${out}_blast.out2
#	mv $pathOut/$out/${out}_blast.out2 $pathOut/$out/${out}_blast.out
	
	if [[ $kmer == "reference" || $kmer == "extended" ]]
	then
		 awk -v var="$out" '{print $0, var}' OFS='\t' $pathOut/$out/${out}_kmer.txt > $pathOut/$out/${out}_kmer.txt2
		 mv $pathOut/$out/${out}_kmer.txt2 $pathOut/$out/${out}_kmer.txt
	 fi
	
done <$file

#Concatenare outputs all genomes
nFile=`wc $file | awk '{print $1}'`
cat ${pathOut}/*/*_blast.out > "${pathOut}/cat_all${nFile}_blast.out"

if [[ $kmer == "reference" || $kmer == "extended" ]]
then
	cat $pathOut/*/*_kmer.txt | grep -v Score > "$pathOut/cat_all${nFile}_kmer.tsv"
fi

exit 0