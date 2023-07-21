
# Use agrep to find primer matches with allowed mismatches in given input fasta

mismatch=1

while getopts i:p:n:s: option
do 
    case "${option}"
        in
        i)input_file=${OPTARG};;
        p)primer=${OPTARG};;
        n)mismatch=${OPTARG};;
        s)strand=${OPTARG};;
    esac
done

if [ $strand == 'F' ]; then
    zcat $input_file | sed -n '2~4p' | awk '{print substr($0,1,50)}' > ./substrings.txt
elif [ $strand == 'R' ]; then
    zcat $input_file | sed -n '2~4p' | awk '{print substr($0,length($0)-49,50)}' > ./substrings.txt
fi

/hps/nobackup/rdf/metagenomics/service-team/users/chrisata/agrep/agrep -V0 -c -$mismatch $primer ./substrings.txt


# zcat $input_file | sed -n '2~2p' | awk '{print substr($0,1,50)substr($0,length($0)-49,50)}' > ./substrings.txt
