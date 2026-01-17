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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <zlib.h>

#include "main.h"

/**@brief parse seed argument and normalize to common representation
 * 
 * The user can define the seed in different formats, 
 * -s l (length), -s m:n (interval), -s m:n/l (interval+length)
 * (Also negative m/n are allowed to represent position from the end)
 * This function normalizes the representation to:
 * seed[0]=m=range start
 * seed[1]=n=range end
 * seed[2]=l=length
 *  
 * @param[out] seed normalized 
 * @param[in] length of the query sequence
 * 
 * @return Void.
 */
void seeds_setup(int *seed, saidx64_t *querylen)
{
	saidx64_t n = *querylen;

	seed[0] = seed_orig[0];
	seed[1] = seed_orig[1];
	seed[2] = seed_orig[2];

	/* length only */
	if(!seed_flag) {
		if(seed[0] <= 0) {
			fprintf(stderr, "Invalid seed length\n");
			abort();
		}
		seed[2] = MIN(seed[0], n);
		seed[0] = 1;
		seed[1] = n;
	}
	/* interval */
	else {
		/* there are many ways of fucking up intervals */
		if((seed[0] > 0 && (seed[1] < seed[0] || seed[1] > n))
				|| (seed[0] < 0 && (seed[1] < seed[0] || seed[0] < -n))
				|| ((seed[0] > 0 && seed[1] < 0) || (seed[0] < 0 && seed[1] > 0))) {
			fprintf(stderr, "Invalid seed interval\n");
			abort();
		}
		
		/* handle negative indexing (3'-end) */
		seed[0] += seed[0]<0 ? 1 : 0;
		seed[1] += seed[1]<0 ? 1 : 0;

		/* indexing should be in interval 1:n */
		seed[0] = (seed[0]+n)%n;
		seed[0] = seed[0] ? seed[0] : n;
		seed[1] = (seed[1]+n)%n;
		seed[1] = seed[1] ? seed[1] : n;

		if(seed[0] <= 0 || seed[1] <= 0 || seed[2] < 0) {
			fprintf(stderr, "Invalid seed interval\n");
			abort();
		}

		/* length within interval */
		if(seed[2] > 0) {
			if(seed[2] > seed[1]-seed[0]+1) {
				fprintf(stderr, "Invalid seed length (exceeds interval)\n");
				abort();
			}
		}
		else seed[2] = seed[1] - seed[0] + 1;
	}
}

