/*********************************************************************
  RIsearch2   --   RNA-RNA interaction search

  Copyright (c) 2016 by the contributors (see AUTHORS file)

  This file is part of RIsearch2.

  RIsearch2 is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  RIsearch2 is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with RIsearch, see file COPYING.
  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/
/* Simple API for FASTA file reading
 * for Bio5495/BME537 Computational Molecular Biology
 * SRE, Sun Sep  8 05:35:11 2002 [AA2721, transatlantic]
 * CVS $Id$
 */

/* Modified by Anders Frost Rudebeck to allow gzipped input */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "zlib.h"

#include "fasta.h"

#define FALLOC  1024

int file_exists(const char * filename)
{
	FILE *file;
	if ((file = fopen(filename, "r"))){
		fclose(file);
		return 1;
	}
	return 0;
}

fasta_t* fasta_open(const char *seqfile)
{
	fasta_t *ffp;

	if(!seqfile) {
		fprintf(stderr, "No file specified in fasta_open()\n");
		abort();
	}

	/*
	if (!file_exists(seqfile)) {
	fprintf(stderr, "File doesn't exist in fasta_open(). It's probably the query file.\n");
	abort();
	}
	*/

	ffp = malloc(sizeof(fasta_t));
	if (ffp == NULL) {
		fprintf(stderr, "Failed allocation in fasta open.\n");
		abort();
	}
	if (strcmp(seqfile, "-")) { /*returns 0/FALSE if they are same! */
		ffp->fp = gzopen(seqfile, "rb");              /* Assume seqfile exists & readable!   */
	} else {
		ffp->fp = gzdopen(STDIN_FILENO, "rb");
	}
	if (ffp->fp == NULL) {
		free(ffp);
		return NULL;
	}
	if ((gzgets(ffp->fp, ffp->buffer, FASTA_MAXLINE)) == NULL) {
		free(ffp);
		return NULL; 
	}
	return ffp;
}

char validrna(char c)
{
	switch(c) {
		case 't': return 'u';
		case 'T': return 'U';
		default : return 'N';
	}
}

int fasta_read(fasta_t *ffp, char **ret_seq, char **ret_name, saidx64_t *ret_L)
{
	char *s;
	char *name;
	char *seq;
	size_t   n;
	size_t   nalloc;
//	int   nonalpha = 0, modified = 0;
	const char *valid = "AaCcGgUuNn";

	/* Peek at the lookahead buffer; see if it appears to be a valid FASTA descline. */
	if (ffp->buffer[0] != '>')
		return 0;

	/* Parse out the name: the first non-whitespace token after the > */
	s  = strtok(ffp->buffer+1, " \t\r\n");
	name = malloc(sizeof(char) * (strlen(s)+1));
	if (name == NULL) {
		fprintf(stderr, "Failed allocation in fasta read.\n");
		abort();
	}
	strcpy(name, s);

	/* Everything else 'til the next descline is the sequence.
	 * Note the idiom for dynamic reallocation of seq as we
	 * read more characters, so we don't have to assume a maximum
	 * sequence length.
	 */
	seq = malloc(sizeof(char) * FALLOC);   /* allocate seq in blocks of residues */
	if (seq == NULL) {
		fprintf(stderr, "Failed allocation in fasta read.\n");
		abort();
	}
	nalloc = FALLOC;
	n = 0;
	while (gzgets(ffp->fp, ffp->buffer, FASTA_MAXLINE)) {
		if (ffp->buffer[0] == '>') break;    /* a-ha, we've reached the next descline */

		for (s = ffp->buffer; *s != '\0'; s++)  {
			if (! isalpha(*s)) {
				/* character is not is not alphabetic, skip it and flag for warning */
//				nonalpha = nonalpha + 1;
				continue;
			}
			/* validate alphabetic characters*/
			if (strchr(valid, *s) == NULL) {
				/* set Tt to Uu and all undefined to N */
				*s = validrna(*s);
//				modified = 1 ;
			}

			seq[n] = *s;                       /* store the character, bump length n */
			n++;
			if (nalloc == n) {                 /* are we out of room in seq? if so, expand */
			                                   /* (remember, need space for the final '\0')*/
				nalloc += FALLOC; 
				seq = realloc(seq, sizeof(char) * nalloc);
				if (seq == NULL) {
					fprintf(stderr, "%s\n", ffp->buffer);
					fprintf(stderr, "%ld - %ld\n", n, nalloc);
					fprintf(stderr, "Failed re-allocation in fasta read.\n");
					abort();
				}
			}
		}
	}
//	if (modified == 1) {
//		fprintf(stderr, "Nucleotide code T has been replaces with U; undefined codes are replaced by N\n");
//	}
//	if (nonalpha > 0) {
//		fprintf(stderr, "%d unexpected (non-alphabetic) character(s) in the fasta sequence skipped\n", nonalpha);
//	}
	seq = realloc(seq, sizeof(char) * (n+1));
	seq[n] = '\0';

	*ret_name = name;
	*ret_seq  = seq;
	*ret_L    = n;
	return 1;
}

void fasta_close(fasta_t *ffp)
{
	if(!ffp) {
		fprintf(stderr, "Attempted to close null-handle in fasta_close()\n");
		abort();
	}

	gzclose(ffp->fp);
	free(ffp);
}

