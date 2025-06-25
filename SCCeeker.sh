#!/bin/bash

#Arguments are
	#$1 absolute path DB
	#$2 absolute path input
	#$3 absolute path output
	#$4 name file output
	#$5 by default false, does not run kmer approch

cd $(dirname $0) 


__usage="
Usage: $(basename $0) [OPTIONS]

Options:
  -1st argument              absolute path to the database directory
  -2nd argument              absolute path to the fasta file
  -3rd argument              absolute path to the output directory 
  -4th argument              prefix to use in all output files
  -5th argument             'reference' or 'extended' depending on the dabase to use for the k-mer approach. No option skip the k-mer approach.
  -h, --help                 print help
  -v, --version              print version
"


if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  	echo "$__usage"
	exit 0
fi

if [[ "$1" == "-v" || "$1" == "--version" ]]; then
        echo "version v.0.1"
        exit 0
fi

if [[ $# < 4 ||  $# > 5 ]]
  then
    echo "Wrong number of arguments: expect 4 or 5 arguments"
    echo "$__usage"
    exit 0
fi
	
#Get names of DBs
DBs=`grep -v "#" $1/selectDB | grep fasta`

#Init i at 1 to start with the 1st DB
i=1

#If the genome is compressed, we need to uncompress t for blast
if [[ "$2" == *"gz"* ]]; then
	input=`echo $2 | sed s/".gz"//`
	gzip -dk $2 
       echo "Uncompressed fasta file"	
else   
	input=$2
	echo "No need for uncompression"
fi 

echo "Query for blast: $input"

#Go through each DB
for db in $DBs
do
	#Get db name (just remove pattern fasta for that)
	dbName=`echo $db | sed s/".fasta"//`
	echo "   Subject for blast $i: $dbName"	

	#Create a new directory for the output
	mkdir -p $3/$4
	#Run blast on db i
	blastn -query $input -subject $1/$db -out $3/$4/$4_tmp${i}.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore"

	#Add 2 columns (one for name output and one for name db)
	awk -v var="$4" '{print $0, var}' OFS='\t' $3/$4/$4_tmp${i}.out  >  $3/$4/$4_tmp${i}_named.out
	awk -v var="$dbName" '{print $0, var}' OFS='\t' $3/$4/$4_tmp${i}_named.out  >  $3/$4/$4_tmp${i}.out

	#Remove intermediate
	rm $3/$4/$4_tmp${i}_named.out
	
	#i = i+1
	i=$((i + 1))
done	

#Concatenate outputs of all db
cat $3/$4/$4_tmp*.out > $3/$4/$4_blast.out

#Delete files containing output of just a single db
rm $3/$4/$4_tmp*.out 


#####################################
## Running the kmer-based approach ##
#####################################

if [[ -z "$5" ]]
then
	echo "Do not proceed to the k-mer approach"
	echo "To processed to the k-mer approach, one has to set the argument 'kmer'"
elif  [[ ! -z "$5" && "$5" == "reference" ]]
then
	echo "Proceed to the k-mer approach with the reference database"
	python3 findtemplate.py -i $input  -t $1/template_db/MyKmerFinder_reference_template -o $3/$4/$4_kmer.txt
elif  [[ ! -z "$5" && "$5" == "extended" ]]
then
        echo "Proceed to the k-mer approach with the extended database"
        python3 findtemplate.py -i $input  -t $1/template_db/MyKmerFinder_extended_template -o $3/$4/$4_kmer.txt
else
	echo "Do not proceed to the k-mer approach"
        echo "The kmer argument allows one to run the k-mer approch and can only take one of the two following values: reference or extended"
fi


####################

#If we had to uncompress the genome, delete the uncompress one
if [[ "$2" == *"gz"* ]]; then
        #Remove the uncmpress genome
        rm $input
fi

exit 0
