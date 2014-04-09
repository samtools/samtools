# allow cmd|md5sum to set $? to 1 if cmd fails but md5sum works
set -o pipefail

# Executes a single test and compares against the expected output
run_test() {
    p=$1; shift
    e=$1; shift
    test_iter=`expr $test_iter + 1`
    #echo "p=$p e=$e cmd=$@"
    result=`eval ${@+"$@"} 2>/dev/null > _out`
    #result=`eval ${@+"$@"} > _out`
    #result=`eval valgrind --error-exitcode=1 --leak-check=full ${@+"$@"}`
    if [ $? != 0 ]
    then
	echo "Error running $@"
    elif cmp -s _out expected/$e
    then
	if [ "$p" != "P" ]
	then
	    echo ""
	    echo "UNEXPECTED PASS: Task worked when we expected failure;"
	    echo "when running $@"
	    echo "See PASS-$e.${test_iter} vs expected/$e"
	    mv _out PASS-$e.${test_iter}
	    nupass=`expr $nupass + 1`
	    return 0
	else
	    nepass=`expr $nepass + 1`
	    return 1
	fi
    else
	if [ "$p" = "F" ]
	then
	    echo ""
	    echo "EXPECTED FAIL: Task failed, but expected to fail;"
	    echo "when running $@"
	    nefail=`expr $nefail + 1`
	    return 1
	else
	    echo ""
	    echo "UNEXPECTED FAIL: Output mismatch for $@"
	    echo "See FAIL-$e.${test_iter} vs expected/$e"
	    mv _out FAIL-$e.${test_iter}
	    nufail=`expr $nufail + 1`
	    return 0
	fi
    fi
}

# Process regresion file
regtest() {
    nupass=0; nepass=0
    nufail=0; nefail=0

    exec 9<"$1"
    while read line <&9
    do
	
	set -- $line
        case $1 in
        "#"*) # skip comments
            ;;
        "")   # skip blank lines too
            ;;
        *)
	    p=$1;shift
	    o=$1;shift
	    cmd=${@+"$@"}

	    if [ "$p" == "INIT" ]
	    then
		eval ${@+"$@"} 2>/dev/null > /dev/null
		continue
	    fi

	    test_iter=0
	    if [ "`expr \"$cmd\" : '.*\$fmt'`" != 0 -a "$samtools" != "./samtools-0.1.19" ]
	    then
		_cmd=`echo $cmd | sed 's/\$fmt/bam/'`
		run_test $p $o $_cmd
		#_cmd=`echo $cmd | sed 's/\$fmt/sam/'`
		#run_test $p $o $_cmd
		_cmd=`echo $cmd | sed 's/\$fmt/cram/'`
		run_test $p $o $_cmd
	    else
		_cmd=`echo $cmd | sed 's/\$fmt/bam/'`
		run_test $p $o $_cmd
	    fi
	    ;;
	esac
    done
    exec 9<&-

    echo ""
    echo "Expected   passes:   $nepass"
    echo "Unexpected passes:   $nupass"
    echo "Expected   failures: $nefail"
    echo "Unexpected failures: $nufail"

    if [ "$nupass" -gt 0 -o "$nufail" -gt 0 ]
    then
	echo "=> FAIL"
	return 1
    else
	echo "=> PASS"
	return 0
    fi
}

echo "Samtools mpileup tests:"

samtools="../../samtools"
bcftools="../../../bcftools/bcftools"
regtest mpileup.reg

# samtools="./samtools-0.1.19"
# bcftools="./bcftools-0.1.19"
# regtest mpileup.reg

exit $?
