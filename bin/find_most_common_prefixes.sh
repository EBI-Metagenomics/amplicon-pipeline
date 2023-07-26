
# This following command is relatively knowledge-free (need approximate primer length and rank list size), 
# but generates all of the  unique strings in the first X characters along with their counts. 
# Checking the first couple (10) should help detect whether primers exist or not

# example run
# bash lib/find_most_common_prefixes.sh 
# -i test/PRJNA714811/SRR13978073_1.fastq.gz 
# -l 17 
# -c 10
# -o test/PRJNA714811/SRR13978073_1_mcp

beg=1
rev=0

while getopts i:c:l:b:r: option
do 
    case "${option}"
        in
        i)input_file=${OPTARG};;
        l)length=${OPTARG};;
        c)count=${OPTARG};;
        b)beg=${OPTARG};;
        r)rev=${OPTARG};;
    esac
done

IFS='.' read -r -a array <<< "$input_file"
file_type=${array[-1]}

if [ $rev == 0 ]
then
    zcat $input_file | sed -n '2~4p' | cut -c$beg-$length | sort | uniq -c | sort -rg
else
    zcat $input_file | sed -n '2~4p'| rev | cut -c$beg-$length | sort | uniq -c | sort -rg
fi


# | head -$count
# | tr -s ' ' | cut -d ' ' -f2 | paste -sd+ | bc
