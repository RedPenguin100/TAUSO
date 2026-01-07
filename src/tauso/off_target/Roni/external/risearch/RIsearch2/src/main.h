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
#ifndef __MAIN_H__
#define __MAIN_H__ 1

#include "dsm.h"
#include "sa.h"
#include "lists.h"
#include "zlib.h"

typedef struct {
	int idx;

	saidx64_t k;
	saidx64_t j;

	int found;
	int score;
	double energy;

	saidx64_t best_left_i;
	saidx64_t best_left_j;

	saidx64_t best_right_i;
	saidx64_t best_right_j;

	saidx64_t new_left;
	saidx64_t new_right;

	saidx64_t q_new_left;
	saidx64_t q_new_right;

	char *id;
	char *result_str;
	char *query_str;
	char *aln_str;
	char *target_str;
	char *ia_str;

	char *temp_query_str;
	char *temp_aln_str;
	char *temp_target_str;
	char *temp_ia_str;
} aln_result_t;

typedef struct {
	char *name;
	char *seq;
	saidx64_t *qsa;
	saidx64_t length;
	sa_interval_list_t *intervals;
	int seed[3];
	gzFile out;
} query_t;

extern int comp[6];
extern int pair[6][6];
extern int pair_noGU[6][6];
extern saidx64_t *idx_lengths, num_idxs;
extern saidx64_t *sum_l;
extern char **name;
extern double min_energy;
extern const dsm_t *S;
extern int seed_flag, seed_orig[3], max_ext_len;
extern double min_seed_energy_per_length;
extern int seed_mismatch[3], seed_mismatch_flag, seed_threshold_flag;
//extern char* strdup(const char*);
extern int verbose;
extern int noGUseed;

extern void debug(const char *msg, ...);
extern void str_rev(char *s);

#endif

