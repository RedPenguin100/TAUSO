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
#ifndef __LISTS_H__
#define __LISTS_H__

#include "sa.h"

#define LIST_FREE(node) while(node) { void *p = (void*)node; node = node->next; free(p); }

struct _list_t {
	struct _list_t *next;
};

struct _sa_interval_list_t {
	struct _sa_interval_list_t *next;
	saidx64_t start, end;
	saidx64_t qstart, qend;
	int slen;
};

struct _risearch_list_t {
	struct _risearch_list_t *next;
	saidx64_t seq_start, seq_end, query_start, query_end;
	double energy;
	const char *name;
	char **hits;
};

typedef struct _list_t list_t;
typedef struct _sa_interval_list_t sa_interval_list_t;
typedef struct _risearch_list_t risearch_list_t;

extern int risearch_list_compare(const list_t *v1, const list_t *v2);
extern risearch_list_t* risearch_list_merge(risearch_list_t *results[5]);
extern list_t* list_unique(list_t *list, int (*compare)(const list_t*, const list_t*), void (*elem_free)(void*));
extern void sa_interval_add(sa_interval_list_t **l, saidx64_t start, saidx64_t end, saidx64_t qstart, saidx64_t qend, int slen);
extern void risearch_list_add(risearch_list_t **l, saidx64_t seq_start, saidx64_t seq_end, saidx64_t query_start, saidx64_t query_end, double energy, const char *name, char **hits);

#endif

