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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "lists.h"

void sa_interval_add(sa_interval_list_t **l, saidx64_t start, saidx64_t end, saidx64_t qstart, saidx64_t qend, int slen)
{
	sa_interval_list_t *nl = (sa_interval_list_t*)malloc(sizeof(sa_interval_list_t));
	memset(nl, 0, sizeof(sa_interval_list_t));
	nl->start  = start;
	nl->end    = end;
	nl->qstart = qstart;
	nl->qend   = qend;
	nl->slen   = slen;
	nl->next   = *l;
	*l         = nl;
}

void risearch_list_add(risearch_list_t **l, saidx64_t seq_start, saidx64_t seq_end, saidx64_t query_start, saidx64_t query_end, double energy, const char *name, char **hits)
{
	risearch_list_t *nl = (risearch_list_t*)malloc(sizeof(risearch_list_t));
	nl->seq_start   = seq_start;
	nl->seq_end     = seq_end;
	nl->query_start = query_start;
	nl->query_end   = query_end;
	nl->energy      = energy;
	nl->name        = name;
	nl->hits        = hits;
	nl->next        = *l;
	*l              = nl;
}

/* merge risearch result lists */
risearch_list_t* risearch_list_merge(risearch_list_t *results[5])
{
	int i; risearch_list_t *p, *p2;

	for(i = 0, p = NULL; i < 5; i++) {
		while(results[i]) {
			p2 = results[i]; results[i] = results[i]->next;
			p2->next = p; p = p2;
		}
	}

	return p;
}

/* compare risearch results */
int risearch_list_compare(const list_t *v1, const list_t *v2)
{
	saidx64_t i; risearch_list_t *l1 = (risearch_list_t*)v1, *l2 = (risearch_list_t*)v2;

	/* always sort NULL last */
	if(!l1 || !l2) return (l1?-1:1);

	if(l1->name != l2->name) return strcmp(l1->name, l2->name); /* return (((size_t)l1->name < (size_t)l2->name)?-1:1); */
	if((i = l1->seq_start-l2->seq_start)) return (i<0?-1:1);
	if((i = l1->seq_end-l2->seq_end)) return (i<0?-1:1);

	/* remaining fields should be equal - no need to test */
	return 0;
}

/* make list unique - based on mergesort */
list_t* list_unique(list_t *list, int (*compare)(const list_t*, const list_t*), void (*elem_free)(void*))
{
	list_t *split[2] = {NULL, NULL}, *tmp, *res, *p = list;
	int i;

	/* boundary condition */
	if(!list || !list->next) return list;

	/* split list into two equal piles */
	for(i = 0; list; i=((i+1)&1)) {
		p = list; list = list->next;
		p->next = split[i]; split[i] = p;
	}

	/* sort and make sublists unique */
	split[0] = split[0]?list_unique(split[0], compare, elem_free):NULL;
	split[1] = split[1]?list_unique(split[1], compare, elem_free):NULL;

	/* pick head of result list */
	if(split[0] && split[1]) i = (compare(split[0], split[1]) <= 0)?0:1; else i = (split[0]?0:1); p = res = split[i]; split[i] = split[i]->next;

	/* merge unique results */
	while(split[0] || split[1]) {
		i = (compare(split[0], split[1]) <= 0)?0:1;

		/* check for uniqueness */
		if(!compare(p, split[i])) { tmp = split[i]; split[i] = split[i]->next; elem_free(tmp); }
		/* append unique element to result list */
		else { p->next = split[i]; split[i] = split[i]->next; p = p->next; }
	}

	/* null-terminate and return result list */
	p->next = NULL; return res;
}

