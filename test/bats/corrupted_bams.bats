#!/usr/bin/env bats

samtools=${BATS_TEST_DIRNAME}/../../samtools

setup() {
  # prepare bam with bad internal bgzf block (corrupted)
  cp ${BATS_TEST_DIRNAME}/../mpileup/ce#large_seq.bam ${BATS_TMPDIR}/corrupted_bgzf_crc.bam
  (echo "corruption" | dd of=${BATS_TMPDIR}/corrupted_bgzf_crc.bam conv=notrunc seek=150000 bs=1 count=10) 2> /dev/null
}

teardown() {
  # remove corrupted bam
  rm -f ${BATS_TMPDIR}/corrupted_bgzf_crc.bam
}

@test "samtools view on corrupted bam exits with error" {
  run ${samtools} view ${BATS_TMPDIR}/corrupted_bgzf_crc.bam
  [ "$status" -ne 0 ]
}

@test "samtools index on corrupted bam exits with error" {
  run ${samtools} index ${BATS_TMPDIR}/corrupted_bgzf_crc.bam
  [ "$status" -ne 0 ]
}

@test "samtools view on corrupted bam prints message about corruption" {
  run ${samtools} view ${BATS_TMPDIR}/corrupted_bgzf_crc.bam
  echo "$output" | grep -i 'corrupt'
}

@test "samtools index on corrupted bam prints message about corruption" {
  run ${samtools} index ${BATS_TMPDIR}/corrupted_bgzf_crc.bam
  echo "$output" | grep -i 'corrupt'
}

@test "samtools view on corrupted bam does not print message about truncation" {
  run ${samtools} view ${BATS_TMPDIR}/corrupted_bgzf_crc.bam
  echo "$output" | grep -vi 'trunca'
}

@test "samtools index on corrupted bam does not print message about truncation" {
  run ${samtools} index ${BATS_TMPDIR}/corrupted_bgzf_crc.bam
  echo "$output" | grep -vi 'trunca'
}

