/* fapacks.cpp - pack FASTA 4 bases per byte
   includes lowercase a,c,g,t

  Copyright (C) 2012, Matt Mahoney, Dell Inc.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

This program produces packed DNA sequences from FASTA files.
The output may be used as a reference genome for the program fastqz.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

int main(int argc, char** argv) {
  if (argc<3)
    printf("To pack FASTA files: fapack output *.fa\n"), exit(1);
  FILE *out=fopen(argv[1], "wb");
  int b=1, c;
  for (int i=2; i<argc; ++i) {
    printf("%s\n", argv[i]);
    FILE *in=fopen(argv[i], "rb");
    if (!in) continue;
    bool dna=true;
    while ((c=getc(in))!=EOF) {
      if (c=='>') dna=false;
      else if (c==10) dna=true;
      if (islower(c)) c=toupper(c);
      if (dna) {
        if (c=='A') b=b*4;
        if (c=='C') b=b*4+1;
        if (c=='G') b=b*4+2;
        if (c=='T') b=b*4+3;
        if (b>=256) putc(b&255, out), b=1;
      }
    }
    if (in) fclose(in);
  }
  fclose(out);
  return 0;
}


