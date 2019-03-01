

Dynamic read filters are dynamic C object loaded at runtime. The interface is  defined in `sam_dynreadfilter.h`  in the structure `SamDynReadFilter`.

Briefly a filter named `myfilter` should implement 3 functions:

  * `myfilter_initialize` that is used to initialize the filter
  * `myfilter_accept` a function returning `1` if the current reads is accepted
  * `myfilter_dispose` releases the resources associated to this filter

## Example

The source contains two examples

  * `misc/libfilterexample1.c` : https://github.com/lindenb/samtools/blob/pl_filter_hook/misc/filterexample1.c accepts reads if they have a clipped element with a length greater or equals to the env variable `${MIN_CIGAR_LEN}`.

  * `misc/libfilterexample2.c` : https://github.com/lindenb/samtools/blob/pl_filter_hook/misc/filterexample2.c accepts reads containing a EchoRI site



```
# update LD_LIBRARY_PATH

export LD_LIBRARY_PATH=/path/to/samtools/misc:${LD_LIBRARY_PATH}

# compile the filters

$ make misc/libfilterexample1.so misc/libfilterexample2.so
(...)
gcc -rdynamic -g -O2 -I. -I../htslib -I./lz4  -shared -fPIC  -o misc/libfilterexample1.so misc/filterexample1.c
gcc -rdynamic -g -O2 -I. -I../htslib -I./lz4  -shared -fPIC  -o misc/libfilterexample2.so misc/filterexample2.c

# invoke samtools without the filters

$ wget -q -O - "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19788/alignment/NA19788.chrom11.ILLUMINA.bwa.MXL.low_coverage.20120522.bam" |\
	head -n 2 | sed 's/GAATTC/.GAATTC./g'
SRR063241.26206328	163	11	60039	0	75M	=	60153	188	AGGCACTTAAATACACTGAAGCTGCCAAAACAATCTATCGTTTTGCCTACGTACTTATCAACTTCCTCATAGAAA	@@?8>BH;AGH6B?C?D?@:BCHCB17@@F964;;FEF=6<<D??D6H12.<:-<7347BJC<F>FA>;E=>6CB	X0:i:10	X1:i:0	MD:Z:72C2	RG:Z:SRR063241	AM:i:0	NM:i:1	SM:i:0	MQ:i:0	XT:A:R	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@O\[
SRR063352.22087164	163	11	60040	0	100M	=	60237	296	GGCACTTAAATACACTGAAGCTGCCAAAACAATCTATCGTTTTGCCTACGTACTTATCAACTTCCTCATAGCAAACTGGGAGAAAAAAGCAATGGAATGA	@DD=EGEFGGFFGHGIHGHIHIIHIIHHHGIIGHJGGI>GHGGIGJJGH@FFHJHDGIIHIIHIIJIHGFEHAGFIJHHJEHGGGIHBHAHGCGFFECDB	X0:i:10	X1:i:0	MD:Z:100	RG:Z:SRR063352	AM:i:0	NM:i:0	SM:i:0	MQ:i:0	XT:A:R	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@D

# invoke samtools with filter1

$ wget -q -O - "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19788/alignment/NA19788.chrom11.ILLUMINA.bwa.MXL.low_coverage.20120522.bam" |\
	MIN_CIGAR_LEN=30 ./samtools view -k filterexample1 - |\
	head -n 2 | sed 's/GAATTC/.GAATTC./g'
SRR063242.12606534	83	11	61196	0	40S35M	=	61010	-220	TGTGTAAGTTTTCCAAACAAAAAGGAACAGCATGAGCAAATGCAAGGAGGCCTAAAATAAAGAGATGTGTAAAGA	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	X0:i:10	X1:i:0	XC:i:35	MD:Z:35	RG:Z:SRR063242	AM:i:0	NM:i:0	SM:i:0	MQ:i:0	XT:A:R	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SRR063353.13443483	163	11	61534	0	35M65S	=	61647	212	CTATCCCTCCCCCCTCCCCCAACCCCACAACCCGCCGCGGGGGGTGATATTCCCCTTCCTGCGTCCTCGGGCTCTCATCGTTCTACCCCCCCCCCCGAGT	@DBCEGEHGEHHHHEGIIII8=<?CC0=########################################################################	X0:i:2	X1:i:88	XC:i:35	MD:Z:31A0G2	RG:Z:SRR063353	AM:i:0	NM:i:2	SM:i:0	MQ:i:0	XT:A:R	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# invoke samtools with filter2

$ wget -q -O - "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19788/alignment/NA19788.chrom11.ILLUMINA.bwa.MXL.low_coverage.20120522.bam" |\
	./samtools view -k filterexample2 - |\
	head -n 2 | sed 's/GAATTC/.GAATTC./g'
SRR063353.19306143	147	11	65662	0	100M	=	65523	-238	ACCTGCTATGTACCCACAAAAATTAAATTTAAAAACAATACATTGTTATCCACTATAGTCACCATATTGCACAATAGATCTGTT.GAATTC.ATTCCTCCTG	@@AD>BBEG8A:BB;CEFHFEB>DGGDDECEHFFHG>FDGH>@HDGHCFIIGJEHDJHGHFJIHHGEIIIGIGGGJHGHIIGHHGFGHGHFGFHGDFED@	X0:i:9	X1:i:1	MD:Z:100	RG:Z:SRR063353	AM:i:0	NM:i:0	SM:i:0	MQ:i:0	XT:A:R	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SRR063353.4519	163	11	65678	0	100M	=	65859	280	CAAAAATTAAATTTAAAAACAATACATTGTTATCCACTATAGTCACCATATTGCACAATAGATCTGTT.GAATTC.ATTCCTCCTGTACAATGCAATTTTGT	?CBDEACCADADEEDHHHEEAHGGGIGGFGG=GHIIHJGGHJHIGHHIHEHEJHIIIIHHKHHIKKIHIIIH?HHHFHIJJKKJHFDIIH?GH?DFECEC	X0:i:10	X1:i:0	MD:Z:100	RG:Z:SRR063353	AM:i:0	NM:i:0	SM:i:0	MQ:i:0	XT:A:R	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# invoke samtools with both filters

$ wget -q -O - "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19788/alignment/NA19788.chrom11.ILLUMINA.bwa.MXL.low_coverage.20120522.bam" |\
	MIN_CIGAR_LEN=30 ./samtools view -k filterexample2 -k filterexample1 - |\
	head -n 2 | sed 's/GAATTC/.GAATTC./g'
SRR063243.17728791	163	11	65732	0	41M34S	=	65815	157	ACAATAGATCTGTT.GAATTC.ATTCCTCCTGTACAATGCAATTTTGTACCCTTTGACCAACATCTACCCAATCCTC	@AE?EEEFF@EAAB<9<><5DA>56550@@GF?GDEC@DC###################################	X0:i:10	X1:i:0	XC:i:41	MD:Z:41	RG:Z:SRR063243	AM:i:0	NM:i:0	SM:i:0	MQ:i:0	XT:A:R	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
SRR063353.24813561	121	11	72794	0	65S35M	=	72794	0	CATAAAAGGAGGCCGCAGGGGTGAGGGGCTTGGTGGCT.GAATTC.CAACAAACACTTAGATGATTAACACCAATCCTTCCCAAACTCTTCCAAAAAAAATG	###########################################################################>:DD=CC<EEGIHJIGGGEGFD<<B	X0:i:5	X1:i:1	XC:i:35	MD:Z:35	RG:Z:SRR063353	AM:i:0	NM:i:0	SM:i:0	XT:A:R	BQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

```
