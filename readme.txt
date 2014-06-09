fastqz15.cpp is the source code for the latest version
of the FASTQ compressor. It compresses the common Sanger
variant. FASTQ is output by DNA sequencing machines.

fapack.cpp is a program to pack FASTA files into a format
suitable for input to fastqz as a reference genome for
better compression.
It packs 4 bases per byte and discards all but A,C,G,T.

fapacks.cpp works the same except that it does not ignore
lowercase a,c,g,t. Lowercase is used in hg19 to indicate
repeats. Generally it produces a larger reference but
gives better compression.

Other fastqz*.cpp are older versions. You don't need them.

Usage: fastqz {c[Q]|d|e[Q]|f} input output [reference]

Command c compresses input to output.fx?.zpaq (3 or 4 files)
Command d decompresses input.fx?.zpaq to output
Command e encodes input to output.fx?
Command f decodes input.fx? to output

Commands c and d are slow, require 1.5 GB memory, use 3 or
4 cores, but get very good compression. Commands e and f are much
faster, use little memory, and only one thread, but compression
ratio is not as good.

Commands cQ or eQ quantize the quality scores for lossy but
better compression. The default is c1 or e1, which is lossless.
Quality scores in the range 33..73 are rounded down to 35 plus
a multiple of Q.

You can supply a reference genome to improve compression.
If you use this, the same reference is needed to decompress.
It also increases the memory requirement to 1.2 GB for the
e command and 0.5 GB for the f command. c and d still need
1.5 GB.

You can prepare the reference genome from FASTA files like:

  fapacks hg19s *.fa

to produce the file hg19s. Then compress:

  fastqz c in.fastq arc hg19s

To decompress:

  fastqz d arc out.fastq hg19s

There are 4 compressed files:

  arc.fxh.zpaq - compressed headers
  arc.fxb.zpaq - compressed base calls
  arc.fxq.zpaq - compressed quality scores
  arc.fxa.zpaq - compressed alignments if a reference is used.

Commands e and f work the same way except the compressed
files do not have a .zpaq extension. If no reference is
used, then no .fxa or .fxa.zpaq file is produced or expected.

fastqz only works on the Sanger FASTQ variant. It assumes
that quality scores are Phred+33 (range ASCII 33 to 73).
Base calls must be A,C,G,T,N only. N must have a quality
score of 0, and all others 1 or higher. Maximum line length
is 4095. Lines must be terminated by linefeeds only (no
carriage returns). If a reference is used, it must be
smaller than 1 GB packed (4 billion bases).

To compile fastqz you will need the latest version of
libzpaq from https://sourceforge.net/projects/zpaq/
or http://mattmahoney.net/zpaq/
These programs will work in either Windows or Linux.
In Windows, you will also need Pthreads-Win32 from
http://sourceware.org/pthreads-win32/ to compile or run.
To compile (no Makefile, sorry):

  g++ -O3 -msse2 -s -lpthread fastqz.cpp libzpaq.cpp -o fastqz
  g++ -O3 -s fapack.cpp -o fapack

fastqz* and fapack* are written by Matt Mahoney, Dell Inc.
All are BSD-2 licensed. But note that libzpaq
is public domain and Pthreads-Win32 is LGPL.
