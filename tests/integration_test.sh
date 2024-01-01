#!/usr/bin/env bash
BIN=bin
PROG=../ms-index
TESTS=testfiles

mkdir -p $BIN

$PROG index -p $TESTS/integration_a.fa -k 3
$PROG index -p $TESTS/integration_b.fa -k 3

$PROG merge -p $TESTS/integration_a.fa -p $TESTS/integration_b.fa -r $BIN/merged.fa

$PROG query -k 3 -p $TESTS/integration_a.fa -q $TESTS/queries.txt > $BIN/a.txt 2> /dev/null
$PROG query -k 3 -p $TESTS/integration_a.fa -q $TESTS/queries.txt -f xor > $BIN/a_xor.txt 2> /dev/null
$PROG query -k 3 -p $TESTS/integration_b.fa -q $TESTS/queries.txt > $BIN/b.txt 2> /dev/null
$PROG query -k 3 -p $TESTS/integration_b.fa -q $TESTS/queries.txt -f xor > $BIN/b_xor.txt 2> /dev/null
$PROG query -k 3 -p $BIN/merged.fa -q $TESTS/queries.txt > $BIN/merged.txt 2> /dev/null
$PROG query -k 3 -p $BIN/merged.fa -q $TESTS/queries.txt -f xor > $BIN/merged_xor.txt 2> /dev/null




diff $TESTS/result_a_complements.txt $BIN/a.txt 
diff $TESTS/result_b_complements.txt $BIN/b.txt 
diff $TESTS/result_b_complements_xor.txt $BIN/b_xor.txt 
diff $TESTS/result_merged_complements.txt $BIN/merged.txt 



