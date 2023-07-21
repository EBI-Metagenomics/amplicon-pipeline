
process std_primer_out_parse {

    input:
    path std_primer_out

    output:
    stdout

    """
    primer_num="\$(wc -l $std_primer_out | cut -d' ' -f 1)"
    if [ \$primer_num -eq 2 ]; then
    	echo -n '2'
    elif [ \$primer_num -eq 1 ]; then
    	echo -n '1'
    elif [ \$primer_num -eq 0 ]; then
    	echo -n '0'
    fi
    """
}
