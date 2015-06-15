#!/usr/bin/env bats

samtools=${BATS_TEST_DIRNAME}/../../samtools

setup() {
  # prepare bam with no issues
  cp ${BATS_TEST_DIRNAME}/../mpileup/ce#large_seq.bam ${BATS_TMPDIR}/ok.bam
}

teardown() {
  rm -f ${BATS_TMPDIR}/ok.bam
}

@test "samtools view on bam with no issues exits successfully" {
  run ${samtools} view ${BATS_TMPDIR}/ok.bam
  [ "$status" -eq 0 ]
}

@test "samtools index on bam with no issues exits successfully" {
  run ${samtools} index ${BATS_TMPDIR}/ok.bam
  [ "$status" -eq 0 ]
}

@test "samtools view on bam with no issues does not print message about truncation" {
  run ${samtools} view ${BATS_TMPDIR}/ok.bam
  echo "$output" | grep -vi 'truncat'
}

@test "samtools index on bam with no issues does not print message about truncation" {
  run ${samtools} index ${BATS_TMPDIR}/ok.bam
  echo "$output" | grep -vi 'truncat'
}

@test "samtools view on bam with no issues does not print message about corruption" {
  run ${samtools} view ${BATS_TMPDIR}/ok.bam
  echo "$output" | grep -vi 'corrupt'
}

@test "samtools index on bam with no issues does not print message about corruption" {
  run ${samtools} index ${BATS_TMPDIR}/ok.bam
  echo "$output" | grep -vi 'corrupt'
}

