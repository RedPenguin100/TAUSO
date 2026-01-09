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
#ifndef __SEARCH_H__
#define __SEARCH_H__ 1

#include "main.h"
#include "sa.h"
#include "lists.h"

extern saidx64_t sa_search_left(const saidx64_t *sa, saidx64_t start, saidx64_t end, saidx64_t offset, saint_t c);
extern void sa_search_interval(const saidx64_t *sa, saidx64_t start, saidx64_t end, saidx64_t offset, saidx64_t interval[6]);
extern void sa_evaluate_interval(sa_interval_list_t *results, query_t *query, const saidx64_t *sa, char strand);
extern void sa_parallel_match_neg(const saidx64_t *qsa, saidx64_t ql, saidx64_t qr, const saidx64_t *sa, saidx64_t sl, saidx64_t sr, saidx64_t offset, saidx64_t depth, int last_match_count, int mismatch_count, sa_interval_list_t **results);

#endif

