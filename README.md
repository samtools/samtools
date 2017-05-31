samtools
========

This is **NOT** the official repository for samtools.<br>
The `mpileup` module has been updated to report a conversion matrix such as 
```
A	=	327690
A	A	0
A	C	1131
A	G	78
A	T	65
A	I	9
A	D	207
C	=	409710
C	A	454
C	C	0
C	G	31
C	T	116
C	I	505
C	D	630
G	=	377175
G	A	77
G	C	38
G	G	0
G	T	671
G	I	2
G	D	35
T	=	410496
T	A	79
T	C	109
T	G	123
T	T	0
T	I	16
T	D	39
```

Use it as 
```bash 
samtools mpileup -f /annotation/mm10.fa sample.bt2mm10.merged.bg.sorted.bam
```
