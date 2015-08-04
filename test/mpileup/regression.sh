#!/bin/sh
#
#    Copyright (C) 2014 Genome Research Ltd.
#
#    Author: James Bonfield <jkb@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# allow cmd|md5sum to set $? to 1 if cmd fails but md5sum works
set +o | grep pipefail >/dev/null && set -o pipefail

# Use a local MD5 directory, which also has the benefit of overriding the
# inbuilt REF_PATH removing the external dependency on EBI.
#
# Note that these MD5 files are truncated versions of their full sequences;
# just enough to pass the tests.
REF_PATH=`pwd`/md5
export REF_PATH

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
        echo "$e: Error running $@"
        mv _out FAIL-$e.${test_iter}
        nufail=`expr $nufail + 1`
        return 0
    elif cmp -s _out expected/$e
    then
        if [ "$p" != "P" ]
        then
            echo ""
            echo "UNEXPECTED PASS: Task worked when we expected failure;" >&2
            echo "when running $@" >&2
            echo "See PASS-$e.${test_iter} vs expected/$e" >&2
            mv _out PASS-$e.${test_iter}
            nupass=`expr $nupass + 1`
            return 0
        else
            nepass=`expr $nepass + 1`
            rm _out
            return 1
        fi
    else
        if [ "$p" = "F" ]
        then
            echo ""
            echo "EXPECTED FAIL: Task failed, but expected to fail;"
            echo "when running $@"
            nefail=`expr $nefail + 1`
            rm _out
            return 1
        else
            echo ""
            echo "UNEXPECTED FAIL: Output mismatch for $@" >&2
            echo "See FAIL-$e.${test_iter} vs expected/$e" >&2
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
    while read -r line <&9
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

            if [ "$p" = "INIT" ]
            then
                #echo "p=$p cmd=$@"
                eval ${@+"$@"} > /dev/null || return 1
                continue
            fi

            test_iter=0
            if [ "`expr \"$cmd\" : '.*\$fmt'`" != 0 -a "$samtools" != "./samtools-0.1.19" ]
            then
                _cmd=`printf '%s' "$cmd" | sed 's/\$fmt/bam/'`
                run_test $p $o $_cmd
                #_cmd=`printf '%s' "$cmd" | sed 's/\$fmt/sam/'`
                #run_test $p $o $_cmd
                _cmd=`printf '%s' "$cmd" | sed 's/\$fmt/cram/'`
                run_test $p $o $_cmd
            else
                _cmd=`printf '%s' "$cmd" | sed 's/\$fmt/bam/'`
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

echo ""
echo "=== Testing $@ regressions ==="

samtools="../../samtools"
filter="../vcf-miniview -f"
regtest $@

# samtools="./samtools-0.1.19"
# filter="./bcftools-0.1.19 view - | sed etc"
# regtest mpileup.reg

exit $?
