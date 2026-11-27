# Run from the fuzz_data/runtime.$fuzzer directory

for i in `seq 1 10000`;do bash -c '(ASAN_OPTIONS=allocator_may_return_null=1:allow_addr2line=1:abort_on_error=1:detect_leaks=0 ../../fuzz_$fuzzer  -detect_leaks=0 corpus -rss_limit_mb=12560 -timeout=1 > /dev/null 2>&1 ) ' 2>/dev/null;done
