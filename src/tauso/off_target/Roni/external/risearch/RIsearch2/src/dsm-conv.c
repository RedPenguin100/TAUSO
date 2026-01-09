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
/* compiler complains about empty source files otherwise */
const int DSM_NOOP;

#ifndef RISEARCH2

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "dsm.h"

extern dsm_t dsm_t99;
extern dsm_t dsm_t04;
extern dsm_t dsm_extend;

void dsm_save(const char *name, dsm_t *dsm)
{
	int i,j,k,l;

	printf("const dsm_t %s = {\\\n", name);

	for(i = 0; i < DSM_N; i++) {
		printf("{\\\n");
		for(j = 0; j < DSM_N; j++) {
			printf("{");
			for(k = 0; k < DSM_N; k++) {
				printf("{");
				for(l = 0; l < DSM_N; l++) {
					printf("%d%s", (*dsm)[i][j][k][l], l==DSM_N-1?"":", ");
				}
				printf("}%s",k==DSM_N-1?" ":", ");
			}
			printf("}%s\\\n",j==DSM_N-1?"":",");
		}
		printf("}%s\\\n",i==DSM_N-1?"":",");
	}

	printf("};\n\n");
}

int main()
{
	int x,i,j,k,l,map[2][6] = {{1,3,2,4,5,0}, {4,2,3,1,5,0}};
	dsm_t dsm[2];


	for(x = 0; x < 2; x++)
		for(i = 0; i < DSM_N; i++)
			for(j = 0; j < DSM_N; j++)
				for(k = 0; k < DSM_N; k++)
					for(l = 0; l < DSM_N; l++)
						dsm[x][map[x][i]][map[x][j]][map[x][k]][map[x][l]] = dsm_extend[i][j][k][l];

	printf("#include \"dsm.h\"\n\n");

	dsm_save("dsm_extend_pos", &dsm[0]);
	dsm_save("dsm_extend_neg", &dsm[1]);

	return 0;
}

#endif

