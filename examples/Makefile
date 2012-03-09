all:../libbam.a ../samtools ../bcftools/bcftools \
	ex1.glf ex1.pileup.gz ex1.bam.bai ex1f-rmduppe.bam ex1f-rmdupse.bam ex1.glfview.gz ex1.bcf calDepth
		@echo; echo \# You can now launch the viewer with: \'samtools tview ex1.bam ex1.fa\'; echo;

ex1.fa.fai:ex1.fa
		../samtools faidx ex1.fa
ex1.bam:ex1.sam.gz ex1.fa.fai
		../samtools import ex1.fa.fai ex1.sam.gz ex1.bam
ex1.bam.bai:ex1.bam
		../samtools index ex1.bam
ex1.pileup.gz:ex1.bam ex1.fa
		../samtools pileup -cf ex1.fa ex1.bam | gzip > ex1.pileup.gz
ex1.glf:ex1.bam ex1.fa
		../samtools pileup -gf ex1.fa ex1.bam > ex1.glf
ex1.glfview.gz:ex1.glf
		../samtools glfview ex1.glf | gzip > ex1.glfview.gz
ex1a.bam:ex1.bam
		../samtools view -h ex1.bam | awk 'BEGIN{FS=OFS="\t"}{if(/^@/)print;else{$$1=$$1"a";print}}' | ../samtools view -bS - > $@
ex1b.bam:ex1.bam
		../samtools view -h ex1.bam | awk 'BEGIN{FS=OFS="\t"}{if(/^@/)print;else{$$1=$$1"b";print}}' | ../samtools view -bS - > $@
ex1f.rg:
		(echo "@RG	ID:ex1	LB:ex1	SM:ex1"; echo "@RG	ID:ex1a	LB:ex1	SM:ex1"; echo "@RG	ID:ex1b	LB:ex1b	SM:ex1b") > $@
ex1f.bam:ex1.bam ex1a.bam ex1b.bam ex1f.rg
		../samtools merge -rh ex1f.rg $@ ex1.bam ex1a.bam ex1b.bam
ex1f-rmduppe.bam:ex1f.bam
		../samtools rmdup ex1f.bam $@
ex1f-rmdupse.bam:ex1f.bam
		../samtools rmdup -S ex1f.bam $@

ex1.bcf:ex1.bam ex1.fa.fai
		../samtools mpileup -gf ex1.fa ex1.bam > $@

../bcftools/bcftools:
		(cd ../bcftools; make bcftools)

../samtools:
		(cd ..; make samtools)

../libbam.a:
		(cd ..; make libbam.a)

calDepth:../libbam.a calDepth.c
		gcc -g -Wall -O2 -I.. calDepth.c -o $@ -L.. -lbam -lm -lz

clean:
		rm -fr *.bam *.bai *.glf* *.fai *.pileup* *~ calDepth *.dSYM ex1*.rg ex1.bcf

# ../samtools pileup ex1.bam|perl -ape '$_=$F[4];s/(\d+)(??{".{$1}"})|\^.//g;@_=(tr/A-Z//,tr/a-z//);$_=join("\t",@F[0,1],@_)."\n"'

# ../samtools pileup -cf ex1.fa ex1.bam|perl -ape '$_=$F[8];s/\^.//g;s/(\d+)(??{".{$1}"})|\^.//g;@_=(tr/A-Za-z//,tr/,.//);$_=join("\t",@F[0,1],@_)."\n"'
