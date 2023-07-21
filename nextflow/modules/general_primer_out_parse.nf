
process general_primer_out_parse {

    input:
    path general_primer_out

    output:
    stdout

    """
    primer_num="\$(grep 1 $general_primer_out | wc -l)"
    if [ \$primer_num -eq 2 ]; then
    	echo -n '2'
    elif [ \$primer_num -eq 1 ]; then
    	echo -n '1'
    elif [ \$primer_num -eq 0 ]; then
    	echo -n '0'
    fi
    """
}
    