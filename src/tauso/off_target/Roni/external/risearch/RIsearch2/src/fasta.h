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
/* fasta.h
 * Declarations for simple FASTA i/o library
 * SRE, Sun Sep  8 05:37:38 2002 [AA2721, transatlantic]
 * CVS $Id$
 */

/* Modified by Anders Frost Rudebeck to allow gzipped input */

#ifndef __FASTA_H__
#define __FASTA_H__

#include <stdio.h>
#include <stdint.h>
#include <divsufsort64.h>
#include "zlib.h"

#define FASTA_MAXLINE 512	/* Requires FASTA file lines to be <512 characters */

typedef struct {
	gzFile fp;
	char  buffer[FASTA_MAXLINE];
} fasta_t;

extern fasta_t   *fasta_open(const char *seqfile);
extern int        fasta_read(fasta_t *fp, char **ret_seq, char **ret_name, saidx64_t *ret_L);
extern void       fasta_close(fasta_t *ffp);


#endif

