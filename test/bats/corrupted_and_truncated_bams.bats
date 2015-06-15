#!/usr/bin/env bats

samtools=${BATS_TEST_DIRNAME}/../../samtools

setup() {
  # prepare bam that is both corrupted and truncated
  cp ${BATS_TEST_DIRNAME}/../mpileup/ce#large_seq.bam ${BATS_TMPDIR}/corrupted_bgzf_crc_bad_eof.bam
  (echo "corruption" | dd of=${BATS_TMPDIR}/corrupted_bgzf_crc_bad_eof.bam conv=notrunc seek=150000 bs=1 count=10) 2> /dev/null
  (echo "ceci-nes-pas-un-block-d-eof" | dd of=${BATS_TMPDIR}/corrupted_bgzf_crc_bad_eof.bam conv=notrunc seek=$(echo $[$(stat -c "%s" ${BATS_TMPDIR}/corrupted_bgzf_crc_bad_eof.bam)-28]) bs=1 count=28) 2> /dev/null
}

teardown() {
  # remove corrupted bam
  rm -f ${BATS_TMPDIR}/corrupted_bgzf_crc_bad_eof.bam
}

@test "samtools view on corrupted and truncated bam exits with error" {
  run ${samtools} view ${BATS_TMPDIR}/corrupted_bgzf_crc_bad_eof.bam
  [ "$status" -ne 0 ]
}

@test "samtools index on corrupted and truncated bam exits with error" {
  run ${samtools} index ${BATS_TMPDIR}/corrupted_bgzf_crc_bad_eof.bam
  [ "$status" -ne 0 ]
}

@test "samtools view on corrupted and truncated bam prints message about truncation" {
  run ${samtools} view ${BATS_TMPDIR}/corrupted_bgzf_crc_bad_eof.bam
  echo "$output" | grep -i 'trunca'
}

@test "samtools view of header only on corrupted and truncated bam prints message about truncation" {
  run ${samtools} view -H ${BATS_TMPDIR}/corrupted_bgzf_crc_bad_eof.bam
  echo "$output" | grep -i 'trunca'
}

@test "samtools index on corrupted and truncated bam prints message about truncation" {
  run ${samtools} index ${BATS_TMPDIR}/corrupted_bgzf_crc_bad_eof.bam
  echo "$output" | grep -i 'trunca'
}

