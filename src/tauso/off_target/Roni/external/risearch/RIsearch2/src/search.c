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
#include <omp.h>
#include <limits.h>
#include <zlib.h>

#include "main.h"
#include "sa.h"
#include "lists.h"

#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
       _a < _b ? _a : _b; })

#define max(a,b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
       _a > _b ? _a : _b; })

//extern int do_extend_seed;
extern int show_alignment;
extern int extPen;

int print_debug = 0;
int print_mats_l = 0;
int print_mats_r = 0;

//saidx64_t *seed_lengths; //array storing how many maximal seeds of each size are found; counts are incremented in the extend_seed function

//reverse complements
const char scomp[256] = {
  ['a'] = 'u',
  ['u'] = 'a',
  ['c'] = 'g',
  ['g'] = 'c',
  ['n'] = 'n'
};

#define KEY_NOT_FOUND 0xffffffffffffff
#define NA INT_MIN/2

int
max3 (int a, int b, int c)
{
  int tmp = (a > b) ? a : b;
  return (c > tmp) ? c : tmp;
}

int
min3 (int a, int b, int c)
{
  int tmp = (a < b) ? a : b;
  return (c < tmp) ? c : tmp;
}

//!!! not used !!!
saidx64_t
binary_search (saidx64_t A[], saidx64_t key, saidx64_t imin, saidx64_t imax)
{
  // continue searching while [imin,imax] is not empty
  while (imax >= imin)
    {
      /* calculate the midpoint for roughly equal partition */
      int imid = imin + ((imax - imin) / 2);
      //fprintf(stderr, "imin: %lu imid: %lu imax: %lu A[imid]: %lu\n", imin, imid, imax, A[imid]);

      // determine which subarray to search
      if (A[imid] < key)
	// change min index to search upper subarray
	imin = imid + 1;
      else if (A[imid] > key)
	// change max index to search lower subarray
	imax = imid - 1;
      else
	// key found at index imid
	return imid;
    }
  // key not found
  return KEY_NOT_FOUND;
}

/**@brief Finds the leftmost position of where the character would have sorted
 *
 * Performs a binary search within the given range. 
 * Locates the leftmost position of a given character at the given offset.
 * Example: Already found the interval of all suffixes starting with agg, now search aggC at offset 3.
 * 
 * @param[in] sa The address to the suffix array (either query or target)
 * @param[in] start The first position (in SA) of the interval to examine/step into
 * @param[in] end The last position (exclusive) of that interval
 * @param[in] offset The number of steps already made into the suffix
 * @param[in] c The character that is searched
 * @return Leftmost position of character c
 */
saidx64_t
sa_search_left (const saidx64_t * sa, saidx64_t start, saidx64_t end,
		saidx64_t offset, saint_t c)
{
  saidx64_t half;

  half = (end - start) >> 1;
  while (start < end)
    {
      //fprintf(stderr, "start: %lu half: %lu offset: %lu XSA(sa[start+half]): %lu\n", start, half, offset, XSA(sa[start + half]));
      saidx64_t xsa_i = XSA (sa[start + half]) + offset;
      //fprintf(stderr, "xsa_i = XSA(sa[start + half]) + offset = %lu ; XSTRM(sa[xsa_i] = %c \n", xsa_i , XSTRM(sa[xsa_i]) );
      if (XSTRM (sa[xsa_i]) >= c)
	end = start + half;
      else
	start += half ? half : 1;
      half >>= 1;
    }
  //fprintf(stderr, "found pos=%lu; XSA(sa[start])+offset=%lu+%lu = %lu ; XSTRM(sa[XSA(sa[start]) + offset] = %c \n", start, XSA(sa[start]), offset, XSA(sa[start])+offset, XSTRM(sa[XSA(sa[start])+offset]));
  return start;
}

/**@brief Finds the intervals of where all characters sorted
 * @param[in] sa The address to the suffix array (either query or target)
 * @param[in] start The first position (in SA) of the interval to examine/step into 
 * @param[in] end The last position (exclusive) of that interval
 * @param[in] offset The number of steps already made into the suffix
 * @param[out] interval The positions of intervals in the range [a,c,g,n,u,$]
 * @return number of all suffixes inside the range (last - first pos)
 */
saidx64_t
sa_search_interval (const saidx64_t * sa, saidx64_t start, saidx64_t end,
		    saidx64_t offset, saidx64_t interval[6])
{
  int i;
  for (i = 0; i < 5; i++)
    {
      //fprintf(stderr, "find leftmost pos of %c:\n", "acgnu"[i]);
      interval[i] = sa_search_left (sa, start, end, offset, "acgnu"[i]);
    }
  interval[5] = end;
  return (end - interval[0]);
}


void
print_dsm_table ()
{
  int i, j, k, l;
  char *letters = "xagcun";

  //print the top row
  fprintf (stderr, "    ");
  for (i = 0; i < 6; i++)
    {
      for (j = 0; j < 6; j++)
	{
	  fprintf (stderr, "   %c%c ", letters[i], letters[j]);
	}
    }
  fprintf (stderr, "\n");

  for (i = 0; i < 6; i++)
    {
      for (j = 0; j < 6; j++)
	{
	  fprintf (stderr, "%c%c  ", letters[i], letters[j]);
	  for (k = 0; k < 6; k++)
	    {
	      for (l = 0; l < 6; l++)
		{
		  fprintf (stderr, "%5d ", (*S)[i][j][k][l]);
		}
	    }
	  fprintf (stderr, "\n");
	}
    }
}


/**
 * @brief Dynamic programming extension "left" of the seed (upstream on query)
 * 
 * Fills left DP matrices, finds max score with position
 *  
 * @return Max score found.
 */
int
DP_left (int *M,		//!<[out] DP matrix for match state
	 int *Bq,		//!<[out] DP matrix for query bulged state
	 int *Bt,		//!<[out] DP matrix for target bulged state
	 int lq,		//!<[in] Length of remaining query (height of DP matrices)
	 int lt,		//!<[in] Length of remaining target (width of DP matrices)
	 const saidx64_t * qsa,	//!<[in] Suffix array of the query
	 const saidx64_t * tsa,	//!<[in] Suffix array of the target
	 saidx64_t q_start,	//!<[in] Start within query sequence (first pos of seed)
	 saidx64_t t_start,	//!<[in] Start within target sequence (first pos of seed)
	 saidx64_t * best_i,	//!<[out] best query pos (as distance from q_start)
	 saidx64_t * best_j)	//!<[out] best target pos (as distance from t_start)
{

// macros to get the sequence letter at an offset from q_start or t_start
#define Q(ix) (XRIS(qsa[q_start-(ix)]))
#define T(ix) (comp[XRIS(tsa[t_start-(ix)])])

#define TI(q,t) ((q)*(lt) + (t))

  int i, j;
  int best_e;

  if (print_debug)
    fprintf (stderr, "dp_left, q_start: %lu t_start: %lu lq: %d lt: %d\n",
	     q_start, t_start, lq, lt);

  M[TI (0, 0)] = 0;
  best_e = (*S)[GAP][Q (0)][GAP][T (0)];
  *best_i = 0;
  *best_j = 0;
  if (lq <= 1 || lt <= 1)
    return best_e;

  Bq[TI (0, 0)] = NA;
  Bt[TI (0, 0)] = NA;
  M[TI (0, 1)] = NA;
  Bq[TI (0, 1)] = NA;
  M[TI (1, 0)] = NA;
  Bt[TI (1, 0)] = NA;
  Bq[TI (1, 1)] = NA;
  Bt[TI (1, 1)] = NA;

  Bt[TI (0, 1)] = (*S)[GAP][Q (0)][T (1)][T (0)];
  Bq[TI (1, 0)] = (*S)[Q (1)][Q (0)][GAP][T (0)];
  M[TI (1, 1)] = (*S)[Q (1)][Q (0)][T (1)][T (0)];
  if (M[TI (1, 1)] + (*S)[GAP][Q (1)][GAP][T (1)] > best_e)
    {
      best_e = M[TI (1, 1)] + (*S)[GAP][Q (1)][GAP][T (1)];
      *best_i = 1;
      *best_j = 1;
    }

  //1. row (i=0), only Bt (gap in query) possible
  for (j = 2; j < lt; ++j)
    {
      M[TI (0, j)] = NA;
      Bq[TI (0, j)] = NA;
      Bq[TI (1, j)] = NA;
      //fprintf(stderr, "t_start - j: %d\n", t_start - j);
      //fprintf(stderr, "TI(j,0): %d\n", TI(j,0));

      // gap here can only be extension from 0,1
      Bt[TI (0, j)] = Bt[TI (0, j - 1)] + (*S)[GAP][GAP][T (j)][T (j - 1)];
      // match can only be bulge closure from Bt, Bq and M are NA for i=0,j=2+
      M[TI (1, j)] = Bt[TI (0, j - 1)] + (*S)[Q (1)][GAP][T (j)][T (j - 1)];
      if (M[TI (1, j)] + (*S)[GAP][Q (1)][GAP][T (j)] > best_e)
	{
	  best_e = M[TI (1, j)] + (*S)[GAP][Q (1)][GAP][T (j)];
	  *best_i = 1;
	  *best_j = j;
	}
    }

  //1. col (j=0), only Bq (gap in target) possible
  for (i = 2; i < lq; ++i)
    {
      M[TI (i, 0)] = NA;
      Bt[TI (i, 0)] = NA;
      Bt[TI (i, 1)] = NA;
      Bq[TI (i, 0)] = Bq[TI (i - 1, 0)] + (*S)[Q (i)][Q (i - 1)][GAP][GAP];
      M[TI (i, 1)] = Bq[TI (i - 1, 0)] + (*S)[Q (i)][Q (i - 1)][T (1)][GAP];
      if (M[TI (i, 1)] + (*S)[GAP][Q (i)][GAP][T (1)] > best_e)
	{
	  best_e = M[TI (i, 1)] + (*S)[GAP][Q (i)][GAP][T (1)];
	  *best_i = i;
	  *best_j = 1;
	}
    }

  if (lq <= 2 || lt <= 2)
    return best_e;

  Bt[TI (1, 2)] = M[TI (1, 1)] + (*S)[GAP][Q (1)][T (2)][T (1)];
  Bq[TI (2, 1)] = M[TI (1, 1)] + (*S)[Q (2)][Q (1)][GAP][T (1)];
  M[TI (2, 2)] = M[TI (1, 1)] + (*S)[Q (2)][Q (1)][T (2)][T (1)];
  if (M[TI (2, 2)] + (*S)[GAP][Q (2)][GAP][T (2)] > best_e)
    {
      best_e = M[TI (2, 2)] + (*S)[GAP][Q (2)][GAP][T (2)];
      *best_i = 2;
      *best_j = 2;
    }
  Bq[TI (2, 2)] = M[TI (1, 2)] + (*S)[Q (2)][Q (1)][GAP][T (2)];
  Bt[TI (2, 2)] = M[TI (2, 1)] + (*S)[GAP][Q (2)][T (2)][T (1)];

  //initialize limited rows
  for (j = 3; j < lt; ++j)
    {
      //i=1: M already set; Bq is NA; Bt as in main recursion
      Bt[TI (1, j)] =
	max (M[TI (1, j - 1)] + (*S)[GAP][Q (1)][T (j)][T (j - 1)],
	     Bt[TI (1, j - 1)] + (*S)[GAP][GAP][T (j)][T (j - 1)]);
      // 3rd row (i=2)
      // Bq is not a possible source for M and Bq
      M[TI (2, j)] =
	max (M[TI (1, j - 1)] + (*S)[Q (2)][Q (1)][T (j)][T (j - 1)],
	     Bt[TI (1, j - 1)] + (*S)[Q (2)][GAP][T (j)][T (j - 1)]);
      Bq[TI (2, j)] = M[TI (1, j)] + (*S)[Q (2)][Q (1)][GAP][T (j)];
      //Bt as in main recursion
      Bt[TI (2, j)] =
	max (M[TI (2, j - 1)] + (*S)[GAP][Q (2)][T (j)][T (j - 1)],
	     Bt[TI (2, j - 1)] + (*S)[GAP][GAP][T (j)][T (j - 1)]);
      if (M[TI (2, j)] + (*S)[GAP][Q (2)][GAP][T (j)] > best_e)
	{
	  best_e = M[TI (2, j)] + (*S)[GAP][Q (2)][GAP][T (j)];
	  *best_i = 2;
	  *best_j = j;
	}
    }

  //initialize limited columns
  for (i = 3; i < lq; ++i)
    {
      //i=1: M already set; Bt is NA; Bq as in main recursion
      Bq[TI (i, 1)] =
	max (M[TI (i - 1, 1)] + (*S)[Q (i)][Q (i - 1)][GAP][T (1)],
	     Bq[TI (i - 1, 1)] + (*S)[Q (i)][Q (i - 1)][GAP][GAP]);
      // 3rd col (j=2)
      // Bt is nor a possible source for M and Bt
      M[TI (i, 2)] =
	max (M[TI (i - 1, 1)] + (*S)[Q (i)][Q (i - 1)][T (2)][T (1)],
	     Bq[TI (i - 1, 1)] + (*S)[Q (i)][Q (i - 1)][T (2)][GAP]);
      Bt[TI (i, 2)] = M[TI (i, 1)] + (*S)[GAP][Q (i)][T (2)][T (1)];
      // Bq as in main recursion
      Bq[TI (i, 2)] =
	max (M[TI (i - 1, 2)] + (*S)[Q (i)][Q (i - 1)][GAP][T (2)],
	     Bq[TI (i - 1, 2)] + (*S)[Q (i)][Q (i - 1)][GAP][GAP]);
      if (M[TI (i, 2)] + (*S)[GAP][Q (i)][GAP][T (2)] > best_e)
	{
	  best_e = M[TI (i, 2)] + (*S)[GAP][Q (i)][GAP][T (2)];
	  *best_i = i;
	  *best_j = 2;
	}
    }

  for (i = 3; i < lq; ++i)
    {
      for (j = 3; j < lt; ++j)
	{
	  M[TI (i, j)] =
	    max3 (M[TI (i - 1, j - 1)] +
		  (*S)[Q (i)][Q (i - 1)][T (j)][T (j - 1)],
		  Bq[TI (i - 1, j - 1)] + (*S)[Q (i)][Q (i - 1)][T (j)][GAP],
		  Bt[TI (i - 1, j - 1)] + (*S)[Q (i)][GAP][T (j)][T (j - 1)]);

	  Bq[TI (i, j)] =
	    max (M[TI (i - 1, j)] + (*S)[Q (i)][Q (i - 1)][GAP][T (j)],
		 Bq[TI (i - 1, j)] + (*S)[Q (i)][Q (i - 1)][GAP][GAP]);
	  Bt[TI (i, j)] =
	    max (M[TI (i, j - 1)] + (*S)[GAP][Q (i)][T (j)][T (j - 1)],
		 Bt[TI (i, j - 1)] + (*S)[GAP][GAP][T (j)][T (j - 1)]);

	  if (M[TI (i, j)] + (*S)[GAP][Q (i)][GAP][T (j)] > best_e)
	    {
	      best_e = M[TI (i, j)] + (*S)[GAP][Q (i)][GAP][T (j)];
	      *best_i = i;
	      *best_j = j;
	    }
	}
    }
#undef Q
#undef T

  if (print_mats_l)
    {
#define QS(ix) (XSTR(qsa[q_start-(ix)]))
#define TS(ix) (XSTR(tsa[t_start-(ix)]))
#define TSc(ix) (scomp[TS(ix)])
      fprintf (stderr, "\n#### DP left ###\n");
      fprintf (stderr, "=== M mat ===\n");
      fprintf (stderr, "%4c%c ", 'x', 'x');
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 1; j < lt; ++j)
	fprintf (stderr, "%4c%c ", TSc (j), TSc (j - 1));
      fprintf (stderr, "\n");
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 0; j < lt; ++j)
	{
	  if (M[TI (0, j)] == NA)
	    {
	      fprintf (stderr, " -NA- ");
	    }
	  else
	    {
	      fprintf (stderr, "%5d ", M[TI (0, j)]);
	    }
	}
      fprintf (stderr, "\n");
      for (i = 1; i < lq; ++i)
	{
	  fprintf (stderr, "%4c%c ", QS (i), QS (i - 1));
	  for (j = 0; j < lt; ++j)
	    {
	      if (M[TI (i, j)] == NA)
		{
		  fprintf (stderr, " -NA- ");
		}
	      else
		{
		  fprintf (stderr, "%5d ", M[TI (i, j)]);
		}
	    }
	  fprintf (stderr, "\n");
	}

      fprintf (stderr, "=== Bq mat ===\n");
      fprintf (stderr, "%4c%c ", 'x', 'x');
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 1; j < lt; ++j)
	fprintf (stderr, "%4c%c ", TSc (j), TSc (j - 1));
      fprintf (stderr, "\n");
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 0; j < lt; ++j)
	{
	  if (Bq[TI (0, j)] == NA)
	    {
	      fprintf (stderr, " -NA- ");
	    }
	  else
	    {
	      fprintf (stderr, "%5d ", Bq[TI (0, j)]);
	    }
	}
      fprintf (stderr, "\n");
      for (i = 1; i < lq; ++i)
	{
	  fprintf (stderr, "%4c%c ", QS (i), QS (i - 1));
	  for (j = 0; j < lt; ++j)
	    {
	      if (Bq[TI (i, j)] == NA)
		{
		  fprintf (stderr, " -NA- ");
		}
	      else
		{
		  fprintf (stderr, "%5d ", Bq[TI (i, j)]);
		}
	    }
	  fprintf (stderr, "\n");
	}

      fprintf (stderr, "=== Bt mat ===\n");
      fprintf (stderr, "%4c%c ", 'x', 'x');
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 1; j < lt; ++j)
	fprintf (stderr, "%4c%c ", TSc (j), TSc (j - 1));
      fprintf (stderr, "\n");
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 0; j < lt; ++j)
	{
	  if (Bt[TI (0, j)] == NA)
	    {
	      fprintf (stderr, " -NA- ");
	    }
	  else
	    {
	      fprintf (stderr, "%5d ", Bt[TI (0, j)]);
	    }
	}
      fprintf (stderr, "\n");
      for (i = 1; i < lq; ++i)
	{
	  fprintf (stderr, "%4c%c ", QS (i), QS (i - 1));
	  for (j = 0; j < lt; ++j)
	    {
	      if (Bt[TI (i, j)] == NA)
		{
		  fprintf (stderr, " -NA- ");
		}
	      else
		{
		  fprintf (stderr, "%5d ", Bt[TI (i, j)]);
		}
	    }
	  fprintf (stderr, "\n");
	}
#undef QS
#undef TS
#undef TSc
    }

  if (print_debug)
    {
      fprintf (stderr, "best_e: %d at ", best_e);
      fprintf (stderr, "best_i: %li, best_j: %li\n", *best_i, *best_j);
    }

#undef TI
  return best_e;
}



/**
 * @brief Dynamic programming extension "right" of the seed (downstream on query)
 * 
 * Fills right DP matrices, finds max score with position
 *  
 * @return Max score found.
 */
int
DP_right (int *M,		//!<[out] DP matrix for match state
	  int *Bq,		//!<[out] DP matrix for query bulged state
	  int *Bt,		//!<[out] DP matrix for target bulged state
	  int lq,		//!<[in] Length of remaining query (height of DP matrices)
	  int lt,		//!<[in] Length of remaining target (width of DP matrices)
	  const saidx64_t * qsa,	//!<[in] Suffix array of the query
	  const saidx64_t * tsa,	//!<[in] Suffix array of the target
	  saidx64_t q_start,	//!<[in] Start within query sequence (last pos of seed)
	  saidx64_t t_start,	//!<[in] Start within target sequence (last pos of seed)
	  saidx64_t * best_i,	//!<[out] best query pos (as distance from q_start)
	  saidx64_t * best_j)	//!<[out] best target pos (as distance from t_start)
{

// macros to get the sequence letter at an offset from q_start or t_start
#define Q(ix) (XRIS(qsa[q_start+(ix)]))
#define T(ix) (comp[XRIS(tsa[t_start+(ix)])])
#define TI(q,t) ((q)*(lt) + (t))

  int i, j;
  int best_e;

  if (print_debug)
    fprintf (stderr, "dp_right, q_start: %lu t_start: %lu lq: %d lt: %d\n",
	     q_start, t_start, lq, lt);

  M[TI (0, 0)] = 0;
  best_e = (*S)[Q (0)][GAP][T (0)][GAP];
  *best_i = 0;
  *best_j = 0;
  if (lq <= 1 || lt <= 1)
    return best_e;

  Bq[TI (0, 0)] = NA;
  Bt[TI (0, 0)] = NA;
  M[TI (0, 1)] = NA;
  Bq[TI (0, 1)] = NA;
  M[TI (1, 0)] = NA;
  Bt[TI (1, 0)] = NA;
  Bq[TI (1, 1)] = NA;
  Bt[TI (1, 1)] = NA;

  Bt[TI (0, 1)] = (*S)[Q (0)][GAP][T (0)][T (1)];
  Bq[TI (1, 0)] = (*S)[Q (0)][Q (1)][T (0)][GAP];
  M[TI (1, 1)] = (*S)[Q (0)][Q (1)][T (0)][T (1)];
  if (M[TI (1, 1)] + (*S)[Q (1)][GAP][T (1)][GAP] > best_e)
    {
      best_e = M[TI (1, 1)] + (*S)[Q (1)][GAP][T (1)][GAP];
      *best_i = 1;
      *best_j = 1;
    }

  //1. row (i=0), only Bt (gap in query) possible
  for (j = 2; j < lt; ++j)
    {
      M[TI (0, j)] = NA;
      Bq[TI (0, j)] = NA;
      Bq[TI (1, j)] = NA;
      //fprintf(stderr, "t_start + j: %d t_start: %d lt: %d\n", t_start + j, t_start, lt);
      //fprintf(stderr, "TI(0,j): %d\n", TI(0,j));
      Bt[TI (0, j)] = Bt[TI (0, j - 1)] + (*S)[GAP][GAP][T (j - 1)][T (j)];
      M[TI (1, j)] = Bt[TI (0, j - 1)] + (*S)[GAP][Q (1)][T (j - 1)][T (j)];
      if (M[TI (1, j)] + (*S)[Q (1)][GAP][T (j)][GAP] > best_e)
	{
	  best_e = M[TI (1, j)] + (*S)[Q (1)][GAP][T (j)][GAP];
	  *best_i = 1;
	  *best_j = j;
	}
    }
  //1. col (j=0), only Bq (gap in target) possible
  for (i = 2; i < lq; ++i)
    {
      M[TI (i, 0)] = NA;
      Bt[TI (i, 0)] = NA;
      Bt[TI (i, 1)] = NA;
      Bq[TI (i, 0)] = Bq[TI (i - 1, 0)] + (*S)[Q (i - 1)][Q (i)][GAP][GAP];
      M[TI (i, 1)] = Bq[TI (i - 1, 0)] + (*S)[Q (i - 1)][Q (i)][GAP][T (1)];
      if (M[TI (i, 1)] + (*S)[Q (i)][GAP][T (1)][GAP] > best_e)
	{
	  best_e = M[TI (i, 1)] + (*S)[Q (i)][GAP][T (1)][GAP];
	  *best_i = i;
	  *best_j = 1;
	}
    }

  if (lq <= 2 || lt <= 2)
    return best_e;

  Bt[TI (1, 2)] = M[TI (1, 1)] + (*S)[Q (1)][GAP][T (1)][T (2)];
  Bq[TI (2, 1)] = M[TI (1, 1)] + (*S)[Q (1)][Q (2)][T (1)][GAP];
  M[TI (2, 2)] = M[TI (1, 1)] + (*S)[Q (1)][Q (2)][T (1)][T (2)];
  if (M[TI (2, 2)] + (*S)[Q (2)][GAP][T (2)][GAP] > best_e)
    {
      best_e = M[TI (2, 2)] + (*S)[Q (2)][GAP][T (2)][GAP];
      *best_i = 2;
      *best_j = 2;
    }
  Bq[TI (2, 2)] = M[TI (1, 2)] + (*S)[Q (1)][Q (2)][T (2)][GAP];
  Bt[TI (2, 2)] = M[TI (2, 1)] + (*S)[Q (2)][GAP][T (1)][T (2)];

  //initialize limited rows
  for (j = 3; j < lt; ++j)
    {
      //i=1: M already set; Bq is NA; Bt as in main recursion
      Bt[TI (1, j)] =
	max (M[TI (1, j - 1)] + (*S)[Q (1)][GAP][T (j - 1)][T (j)],
	     Bt[TI (1, j - 1)] + (*S)[GAP][GAP][T (j - 1)][T (j)]);
      // 3rd row (i=2)
      // Bq is not a possible source for M and Bq
      M[TI (2, j)] =
	max (M[TI (1, j - 1)] + (*S)[Q (1)][Q (2)][T (j - 1)][T (j)],
	     Bt[TI (1, j - 1)] + (*S)[GAP][Q (2)][T (j - 1)][T (j)]);
      Bq[TI (2, j)] = M[TI (1, j)] + (*S)[Q (1)][Q (2)][T (j)][GAP];
      //Bt as in main recursion
      Bt[TI (2, j)] =
	max (M[TI (2, j - 1)] + (*S)[Q (2)][GAP][T (j - 1)][T (j)],
	     Bt[TI (2, j - 1)] + (*S)[GAP][GAP][T (j - 1)][T (j)]);
      if (M[TI (2, j)] + (*S)[Q (2)][GAP][T (j)][GAP] > best_e)
	{
	  best_e = M[TI (2, j)] + (*S)[Q (2)][GAP][T (j)][GAP];
	  *best_i = 2;
	  *best_j = j;
	}
    }

  //initialize limited columns
  for (i = 3; i < lq; ++i)
    {
      //i=1: M already set; Bt is NA; Bq as in main recursion
      Bq[TI (i, 1)] =
	max (M[TI (i - 1, 1)] + (*S)[Q (i - 1)][Q (i)][T (1)][GAP],
	     Bq[TI (i - 1, 1)] + (*S)[Q (i - 1)][Q (i)][GAP][GAP]);
      // 3rd col (j=2)
      // Bt is not a possible source for M and Bt:
      M[TI (i, 2)] =
	max (M[TI (i - 1, 1)] + (*S)[Q (i - 1)][Q (i)][T (1)][T (2)],
	     Bq[TI (i - 1, 1)] + (*S)[Q (i - 1)][Q (i)][GAP][T (2)]);
      Bt[TI (i, 2)] = M[TI (i, 1)] + (*S)[Q (i)][GAP][T (1)][T (2)];
      // Bq as in main recursion
      Bq[TI (i, 2)] =
	max (M[TI (i - 1, 2)] + (*S)[Q (i - 1)][Q (i)][T (2)][GAP],
	     Bq[TI (i - 1, 2)] + (*S)[Q (i - 1)][Q (i)][GAP][GAP]);
      if (M[TI (i, 2)] + (*S)[Q (i)][GAP][T (2)][GAP] > best_e)
	{
	  best_e = M[TI (i, 2)] + (*S)[Q (i)][GAP][T (2)][GAP];
	  *best_i = i;
	  *best_j = 2;
	}
    }

  for (i = 3; i < lq; ++i)
    {
      for (j = 3; j < lt; ++j)
	{
	  M[TI (i, j)] =
	    max3 (M[TI (i - 1, j - 1)] +
		  (*S)[Q (i - 1)][Q (i)][T (j - 1)][T (j)],
		  Bq[TI (i - 1, j - 1)] + (*S)[Q (i - 1)][Q (i)][GAP][T (j)],
		  Bt[TI (i - 1, j - 1)] + (*S)[GAP][Q (i)][T (j - 1)][T (j)]);

	  Bq[TI (i, j)] =
	    max (M[TI (i - 1, j)] + (*S)[Q (i - 1)][Q (i)][T (j)][GAP],
		 Bq[TI (i - 1, j)] + (*S)[Q (i - 1)][Q (i)][GAP][GAP]);
	  Bt[TI (i, j)] =
	    max (M[TI (i, j - 1)] + (*S)[Q (i)][GAP][T (j - 1)][T (j)],
		 Bt[TI (i, j - 1)] + (*S)[GAP][GAP][T (j - 1)][T (j)]);

	  if (M[TI (i, j)] + (*S)[Q (i)][GAP][T (j)][GAP] > best_e)
	    {
	      best_e = M[TI (i, j)] + (*S)[Q (i)][GAP][T (j)][GAP];
	      *best_i = i;
	      *best_j = j;
	    }
	}
    }
#undef Q
#undef T

  if (print_mats_r)
    {
#define QS(ix) (XSTR(qsa[q_start+(ix)]))
#define TS(ix) (XSTR(tsa[t_start+(ix)]))
#define TSc(ix) (scomp[TS(ix)])
      fprintf (stderr, "\n#### DP right ###\n");
      fprintf (stderr, "=== M mat ===\n");
      fprintf (stderr, "%4c%c ", 'x', 'x');
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 1; j < lt; ++j)
	fprintf (stderr, "%4c%c ", TSc (j - 1), TSc (j));
      fprintf (stderr, "\n");
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 0; j < lt; ++j)
	{
	  if (M[TI (0, j)] == NA)
	    {
	      fprintf (stderr, " -NA- ");
	    }
	  else
	    {
	      fprintf (stderr, "%5d ", M[TI (0, j)]);
	    }
	}
      fprintf (stderr, "\n");
      for (i = 1; i < lq; ++i)
	{
	  fprintf (stderr, "%4c%c ", QS (i - 1), QS (i));
	  for (j = 0; j < lt; ++j)
	    {
	      if (M[TI (i, j)] == NA)
		{
		  fprintf (stderr, " -NA- ");
		}
	      else
		{
		  fprintf (stderr, "%5d ", M[TI (i, j)]);
		}
	    }
	  fprintf (stderr, "\n");
	}

      fprintf (stderr, "=== Bq mat ===\n");
      fprintf (stderr, "%4c%c ", 'x', 'x');
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 1; j < lt; ++j)
	fprintf (stderr, "%4c%c ", TSc (j - 1), TSc (j));
      fprintf (stderr, "\n");
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 0; j < lt; ++j)
	{
	  if (Bq[TI (0, j)] == NA)
	    {
	      fprintf (stderr, " -NA- ");
	    }
	  else
	    {
	      fprintf (stderr, "%5d ", Bq[TI (0, j)]);
	    }
	}
      fprintf (stderr, "\n");
      for (i = 1; i < lq; ++i)
	{
	  fprintf (stderr, "%4c%c ", QS (i - 1), QS (i));
	  for (j = 0; j < lt; ++j)
	    {
	      if (Bq[TI (i, j)] == NA)
		{
		  fprintf (stderr, " -NA- ");
		}
	      else
		{
		  fprintf (stderr, "%5d ", Bq[TI (i, j)]);
		}
	    }
	  fprintf (stderr, "\n");
	}

      fprintf (stderr, "=== Bt mat ===\n");
      fprintf (stderr, "%4c%c ", 'x', 'x');
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 1; j < lt; ++j)
	fprintf (stderr, "%4c%c ", TSc (j - 1), TSc (j));
      fprintf (stderr, "\n");
      fprintf (stderr, "%4c%c ", '-', '-');
      for (j = 0; j < lt; ++j)
	{
	  if (Bt[TI (0, j)] == NA)
	    {
	      fprintf (stderr, " -NA- ");
	    }
	  else
	    {
	      fprintf (stderr, "%5d ", Bt[TI (0, j)]);
	    }
	}
      fprintf (stderr, "\n");
      for (i = 1; i < lq; ++i)
	{
	  fprintf (stderr, "%4c%c ", QS (i - 1), QS (i));
	  for (j = 0; j < lt; ++j)
	    {
	      if (Bt[TI (i, j)] == NA)
		{
		  fprintf (stderr, " -NA- ");
		}
	      else
		{
		  fprintf (stderr, "%5d ", Bt[TI (i, j)]);
		}
	    }
	  fprintf (stderr, "\n");
	}
#undef QS
#undef TS
#undef TSc
    }

  if (print_debug)
    {
      fprintf (stderr, "best_e: %d at ", best_e);
      fprintf (stderr, "best_i: %li, best_j: %li\n", *best_i, *best_j);
    }

#undef TI
  return best_e;
}

#define MA 0
#define BQ 1
#define BT 2

/**
 * @brief Backtrack to get the actual interacting string for output
 * 
 * Does backtrack on left and right DP matrices
 * Gets alignment-style interaction
 *  
 * @return Void.
 */
void
print_alignment (const saidx64_t * sa,	//!<[in] Suffix array of the target
		 const saidx64_t * qsa,	//!<[in] Suffix array of the query
		 aln_result_t * result,	//!<[in,out] result struct: gets positions (seed and best dp), sets query/aln/target_str and new_ pos
		 saidx64_t query_len,	//!<[in] Length of the query
		 int seed_len,	//!<[in] Length of the seed
		 int *M_left,	//!<[in] DP matrix for left extension, match state
		 int *Bq_left,	//!<[in] DP matrix for left extension, query bulged state
		 int *Bt_left,	//!<[in] DP matrix for left extension, target bulged state
		 int *M_right,	//!<[in] DP matrix for right extension, match state
		 int *Bq_right,	//!<[in] DP matrix for right extension, query bulged state
		 int *Bt_right)	//!<[in] DP matrix for right extension, target bulged state
{
  // MISMATCH UPDATE
  int q_seed_len = seed_len, t_seed_len = seed_len, qi, ti;


  int i, j, state = MA;

  char *query_str = result->query_str;
  char *temp_query_str = result->temp_query_str;
  char *aln_str = result->aln_str;
  char *temp_aln_str = result->temp_aln_str;
  char *target_str = result->target_str;
  char *temp_target_str = result->temp_target_str;
  char buf[4];
  saidx64_t idx;

  int rem_target_len_left = 0;
  int rem_target_len_right = 0;

  *temp_query_str = 0;
  *temp_aln_str = 0;
  *temp_target_str = 0;
  *query_str = 0;
  *aln_str = 0;
  *target_str = 0;

  saidx64_t q_start = XSA (qsa[result->k]);
  saidx64_t t_start = XSA (sa[result->j]);

  idx = XIDX (sa[result->j]);

  rem_target_len_left = min (	// general limiting parameter
			      max_ext_len,
			      // prevent run in previous seq
			      t_start - sum_l[idx] + 1);

  rem_target_len_right = min (	// general limiting parameter
			       max_ext_len,
			       // prevent overrun into next sequence / out of last
			       sum_l[idx + 1] - (t_start + t_seed_len - 1));

  int *M;
  int *Bq;
  int *Bt;

/*** backtrack left ***/

  M = M_left;
  Bq = Bq_left;
  Bt = Bt_left;

#define Q(ix) (XRIS(qsa[q_start-(ix)]))
#define T(ix) (XRIS(sa[t_start-(ix)]))
#define Tc(ix) (comp[T(ix)])
#define QS(ix) (XSTR(qsa[q_start-(ix)]))
#define TS(ix) (XSTR(sa[t_start-(ix)]))
#define TSc(ix) (scomp[TS(ix)])
#define TI(q,t) ((q)*(rem_target_len_left) + (t))

  i = result->best_left_i;
  j = result->best_left_j;

  result->new_left = t_start - j;
  result->q_new_left = q_start - i;

  if (print_debug)
    fprintf (stderr,
	     "printing, backtrack left from i: %d j: %d t_start: %lu\n", i, j,
	     t_start);

  state = MA;

  while (i > 0 || j > 0)
    {
      //fprintf(stderr, "M[TI(%d,%d)]: %d state: %d\n", i, j, M[TI(i,j)], state);
      if (state == MA)
	{
	  sprintf (buf, "%c", QS (i));
	  strcat (temp_query_str, buf);

	  if (Q (i) == T (j) && Q (i) != 5)	//same nt, not 'N'
	    strcat (temp_aln_str, "|");
	  else if (
	    /*(noGUseed==0) & this would be necessary if you don't want to show G-Us with ":" */
	    ((Q (i) == 2 && T (j) == 1) || (Q (i) == 4 && T (j) == 3)))
	    //ga or uc corresponds to wobble as T is reversed
	    strcat (temp_aln_str, ":");
	  else
	    strcat (temp_aln_str, " ");

	  sprintf (buf, "%c", TSc (j));
	  strcat (temp_target_str, buf);

	  if (M[TI (i, j)] ==
	      M[TI (i - 1, j - 1)] +
	      (*S)[Q (i)][Q (i - 1)][Tc (j)][Tc (j - 1)])
	    {
	      state = MA;
	    }
	  else if (M[TI (i, j)] ==
		   Bq[TI (i - 1, j - 1)] +
		   (*S)[Q (i)][Q (i - 1)][Tc (j)][GAP])
	    {
	      state = BQ;
	    }
	  else if (M[TI (i, j)] ==
		   Bt[TI (i - 1, j - 1)] +
		   (*S)[Q (i)][GAP][Tc (j)][Tc (j - 1)])
	    {
	      state = BT;
	    }
	  i--;
	  j--;
	}
      else if (state == BQ)
	{
	  sprintf (buf, "%c", QS (i));
	  strcat (temp_query_str, buf);
	  strcat (temp_aln_str, " ");
	  sprintf (buf, "-");
	  strcat (temp_target_str, buf);

	  if (Bq[TI (i, j)] ==
	      M[TI (i - 1, j)] + (*S)[Q (i)][Q (i - 1)][GAP][Tc (j)])
	    {
	      state = MA;
	      i--;
	    }
	  else if (Bq[TI (i, j)] ==
		   Bq[TI (i - 1, j)] + (*S)[Q (i)][Q (i - 1)][GAP][GAP])
	    {
	      state = BQ;
	      i--;
	    }
	}
      else if (state == BT)
	{
	  sprintf (buf, "-");
	  strcat (temp_query_str, buf);
	  strcat (temp_aln_str, " ");
	  sprintf (buf, "%c", TSc (j));
	  strcat (temp_target_str, buf);

	  if (Bt[TI (i, j)] ==
	      M[TI (i, j - 1)] + (*S)[GAP][Q (i)][Tc (j)][Tc (j - 1)])
	    {
	      state = MA;
	      j--;
	    }
	  else if (Bt[TI (i, j)] ==
		   Bt[TI (i, j - 1)] + (*S)[GAP][GAP][Tc (j)][Tc (j - 1)])
	    {
	      state = BT;
	      j--;
	    }
	}
    }

  strcat (query_str, temp_query_str);
  strcat (aln_str, temp_aln_str);
  strcat (target_str, temp_target_str);

  *temp_query_str = 0;
  *temp_aln_str = 0;
  *temp_target_str = 0;
#undef Q
#undef QS
#undef T
#undef Tc
#undef TS
#undef TSc
#undef TI

#define Q(ix) (XRIS(qsa[q_start+(ix)]))
#define QS(ix) (XSTR(qsa[q_start+(ix)]))
#define T(ix) (XRIS(sa[t_start+(ix)]))
#define Tc(ix) (comp[T(ix)])
#define TS(ix) (XSTR(sa[t_start+(ix)]))
#define TSc(ix) (scomp[TS(ix)])
#define TI(q,t) ((q)*(rem_target_len_right) + (t))
/*** seed interaction ***/

  if (print_debug)
    {
      strcat (query_str, "y");
      strcat (aln_str, "y");
      strcat (target_str, "y");
    }


  for (i = 0; i < seed_len; i++)
    {
      qi = i;
      ti = i;

      sprintf (buf, "%c", QS (qi));
      strcat (query_str, buf);
      if (Q (qi) == T (ti))
	strcat (aln_str, "|");
      else if (pair[Q (qi)][Tc (ti)])
	strcat (aln_str, ":");
      else
	strcat (aln_str, " ");
      sprintf (buf, "%c", TSc (ti));
      strcat (target_str, buf);
    }
  if (print_debug)
    {
      strcat (query_str, "x");
      strcat (aln_str, "x");
      strcat (target_str, "x");
    }

/*** backtrack right ***/
  M = M_right;
  Bq = Bq_right;
  Bt = Bt_right;

  t_start = t_start + t_seed_len - 1;
  q_start = q_start + q_seed_len - 1;

  i = result->best_right_i;
  j = result->best_right_j;

  result->new_right = t_start + j;
  result->q_new_right = q_start + i;

  if (print_debug)
    fprintf (stderr,
	     "printing, backtrack right from i: %d j: %d t_start: %lu \n", i,
	     j, t_start);

  state = MA;

  while (i > 0 || j > 0)
    {
      if (state == MA)
	{
	  sprintf (buf, "%c", QS (i));
	  strcat (temp_query_str, buf);

	  if (Q (i) == T (j) && Q (i) != 5)	// same nt, not 'N'
	    strcat (temp_aln_str, "|");
	  else if (
	    /*(noGUseed==0) & this would be necessary if you don't want to show G-Us with ":" */
	    ((Q (i) == 2 && T (j) == 1) || (Q (i) == 4 && T (j) == 3)))
	    //ga or uc corresponds to wobble as T is reversed
	    strcat (temp_aln_str, ":");
	  else
	    strcat (temp_aln_str, " ");

	  sprintf (buf, "%c", TSc (j));
	  strcat (temp_target_str, buf);

	  if (M[TI (i, j)] ==
	      M[TI (i - 1, j - 1)] +
	      (*S)[Q (i - 1)][Q (i)][Tc (j - 1)][Tc (j)])
	    {
	      state = MA;
	    }
	  else if (M[TI (i, j)] ==
		   Bq[TI (i - 1, j - 1)] +
		   (*S)[Q (i - 1)][Q (i)][GAP][Tc (j)])
	    {
	      state = BQ;

	    }
	  else if (M[TI (i, j)] ==
		   Bt[TI (i - 1, j - 1)] +
		   (*S)[GAP][Q (i)][Tc (j - 1)][Tc (j)])
	    {
	      state = BT;

	    }
	  i--;
	  j--;
	}
      else if (state == BQ)
	{
	  sprintf (buf, "%c", QS (i));
	  strcat (temp_query_str, buf);
	  strcat (temp_aln_str, " ");
	  sprintf (buf, "-");
	  strcat (temp_target_str, buf);

	  if (Bq[TI (i, j)] ==
	      M[TI (i - 1, j)] + (*S)[Q (i - 1)][Q (i)][Tc (j)][GAP])
	    {
	      state = MA;

	      i--;
	    }
	  else if (Bq[TI (i, j)] ==
		   Bq[TI (i - 1, j)] + (*S)[Q (i - 1)][Q (i)][GAP][GAP])
	    {
	      state = BQ;
	      i--;
	    }
	}
      else if (state == BT)
	{
	  sprintf (buf, "-");
	  strcat (temp_query_str, buf);
	  strcat (temp_aln_str, " ");
	  sprintf (buf, "%c", TSc (j));
	  strcat (temp_target_str, buf);

	  if (Bt[TI (i, j)] ==
	      M[TI (i, j - 1)] + (*S)[Q (i)][GAP][Tc (j - 1)][Tc (j)])
	    {
	      state = MA;
	      j--;
	    }
	  else if (Bt[TI (i, j)] ==
		   Bt[TI (i, j - 1)] + (*S)[GAP][GAP][Tc (j - 1)][Tc (j)])
	    {
	      state = BT;
	      j--;
	    }
	}
    }


  str_rev (temp_query_str);
  str_rev (temp_aln_str);
  str_rev (temp_target_str);

  strcat (query_str, temp_query_str);
  strcat (aln_str, temp_aln_str);
  strcat (target_str, temp_target_str);

#undef Q
#undef QS
#undef T
#undef Tc
#undef TS
#undef TSc
#undef TI

  /*
     result->query_str = query_str;
     result->aln_str = aln_str;
     result->target_str = target_str;
   */
}


/**
 * @brief Backtrack to get the condensed interacting string for output
 * 
 * Does backtrack on left and right DP matrices
 * Gets CIGAR-style interaction (for inline display)
 *  
 * @return Void.
 */
void
print_ali2 (const saidx64_t * sa,	//!<[in] Suffix array of the target
	    const saidx64_t * qsa,	//!<[in] Suffix array of the query
	    aln_result_t * result,	//!<[in,out] result struct: gets positions (seed and best dp), sets ia_str and new_ pos
	    saidx64_t query_len,	//!<[in] Length of the query
	    int seed_len,	//!<[in] Length of the seed
	    int *M_left,	//!<[in] DP matrix for left extension, match state
	    int *Bq_left,	//!<[in] DP matrix for left extension, query bulged state
	    int *Bt_left,	//!<[in] DP matrix for left extension, target bulged state
	    int *M_right,	//!<[in] DP matrix for right extension, match state
	    int *Bq_right,	//!<[in] DP matrix for right extension, query bulged state
	    int *Bt_right)	//!<[in] DP matrix for right extension, target bulged state
{
  // MISMATCH UPDATE
  int q_seed_len = seed_len, t_seed_len = seed_len, qi, ti;


  int i, j, state = MA;

  char *ia_str = result->ia_str;
  char *temp_ia_str = result->temp_ia_str;
  saidx64_t idx;

  int rem_target_len_left = 0;
  int rem_target_len_right = 0;

  *temp_ia_str = 0;
  *ia_str = 0;

  saidx64_t q_start = XSA (qsa[result->k]);
  saidx64_t t_start = XSA (sa[result->j]);

  idx = XIDX (sa[result->j]);

  rem_target_len_left = min (	// general limiting parameter
			      max_ext_len,
			      // prevent run in previous seq
			      t_start - sum_l[idx] + 1);

  rem_target_len_right = min (	// general limiting parameter
			       max_ext_len,
			       // prevent overrun into next sequence / out of last
			       sum_l[idx + 1] - (t_start + t_seed_len - 1));

  int *M;
  int *Bq;
  int *Bt;

/*** backtrack left ***/

  M = M_left;
  Bq = Bq_left;
  Bt = Bt_left;

#define Q(ix) (XRIS(qsa[q_start-(ix)]))
#define T(ix) (XRIS(sa[t_start-(ix)]))
#define Tc(ix) (comp[T(ix)])
//#define QS(ix) (XSTR(qsa[q_start-(ix)]))
//#define TS(ix) (XSTR(sa[t_start-(ix)]))
#define TI(q,t) ((q)*(rem_target_len_left) + (t))

  i = result->best_left_i;
  j = result->best_left_j;

  result->new_left = t_start - j;
  result->q_new_left = q_start - i;

  if (print_debug)
    fprintf (stderr,
	     "printing, backtrack left from i: %d j: %d t_start: %lu\n", i, j,
	     t_start);

  state = MA;

  while (i > 0 || j > 0)
    {
      //fprintf(stderr, "M[TI(%d,%d)]: %d state: %d\n", i, j, M[TI(i,j)], state);
      if (state == MA)
	{

	  if (Q (i) == T (j) && Q (i) != 5)	//same nt, not 'N'
	    strcat (ia_str, "P");
	  else if ((Q (i) == 2 && T (j) == 1) || (Q (i) == 4 && T (j) == 3))
	    //ga or uc corresponds to wobble as T is reversed
	    strcat (ia_str, "W");
	  else
	    strcat (ia_str, "U");

	  if (M[TI (i, j)] ==
	      M[TI (i - 1, j - 1)] +
	      (*S)[Q (i)][Q (i - 1)][Tc (j)][Tc (j - 1)])
	    {
	      state = MA;
	    }
	  else if (M[TI (i, j)] ==
		   Bq[TI (i - 1, j - 1)] +
		   (*S)[Q (i)][Q (i - 1)][Tc (j)][GAP])
	    {
	      state = BQ;
	    }
	  else if (M[TI (i, j)] ==
		   Bt[TI (i - 1, j - 1)] +
		   (*S)[Q (i)][GAP][Tc (j)][Tc (j - 1)])
	    {
	      state = BT;
	    }
	  i--;
	  j--;
	}
      else if (state == BQ)
	{
	  /* sprintf(buf, "%c", QS(i)) >> query_str ;  "-" >> target_str ; */
	  strcat (ia_str, "Q");

	  if (Bq[TI (i, j)] ==
	      M[TI (i - 1, j)] + (*S)[Q (i)][Q (i - 1)][GAP][Tc (j)])
	    {
	      state = MA;
	      i--;
	    }
	  else if (Bq[TI (i, j)] ==
		   Bq[TI (i - 1, j)] + (*S)[Q (i)][Q (i - 1)][GAP][GAP])
	    {
	      state = BQ;
	      i--;
	    }
	}
      else if (state == BT)
	{
	  /* sprintf(buf, "%c", TSc(j)) >> target_str ;  "-" >> query_str ; */
	  strcat (ia_str, "T");

	  if (Bt[TI (i, j)] ==
	      M[TI (i, j - 1)] + (*S)[GAP][Q (i)][Tc (j)][Tc (j - 1)])
	    {
	      state = MA;
	      j--;
	    }
	  else if (Bt[TI (i, j)] ==
		   Bt[TI (i, j - 1)] + (*S)[GAP][GAP][Tc (j)][Tc (j - 1)])
	    {
	      state = BT;
	      j--;
	    }
	}
    }


#undef Q
#undef T
#undef Tc
#undef TI

#define Q(ix) (XRIS(qsa[q_start+(ix)]))
//#define QS(ix) (XSTR(qsa[q_start+(ix)]))
#define T(ix) (XRIS(sa[t_start+(ix)]))
#define Tc(ix) (comp[T(ix)])
//#define TS(ix) (XSTR(sa[t_start+(ix)]))
//#define TSc(ix) (scomp[TS(ix)])
#define TI(q,t) ((q)*(rem_target_len_right) + (t))
/*** seed interaction ***/

  if (print_debug)
    {
      strcat (ia_str, "y");
    }

  for (i = 0; i < seed_len; i++)
    {
      qi = i;
      ti = i;

      if (Q (qi) == T (ti))
	strcat (ia_str, "P");
      else if (pair[Q (qi)][Tc (ti)])
	strcat (ia_str, "W");
      else
	strcat (ia_str, "U");
    }

  if (print_debug)
    {
      strcat (ia_str, "x");
    }

/*** backtrack right ***/
  M = M_right;
  Bq = Bq_right;
  Bt = Bt_right;

  t_start = t_start + t_seed_len - 1;
  q_start = q_start + q_seed_len - 1;

  i = result->best_right_i;
  j = result->best_right_j;

  result->new_right = t_start + j;
  result->q_new_right = q_start + i;

  if (print_debug)
    fprintf (stderr,
	     "printing, backtrack right from i: %d j: %d t_start: %lu\n", i,
	     j, t_start);

  state = MA;

  while (i > 0 || j > 0)
    {
      if (state == MA)
	{
	  if (Q (i) == T (j) && Q (i) != 5)	//same nt, not 'N'
	    strcat (temp_ia_str, "P");
	  else if ((Q (i) == 2 && T (j) == 1) || (Q (i) == 4 && T (j) == 3))
	    strcat (temp_ia_str, "W");
	  else
	    strcat (temp_ia_str, "U");

	  if (M[TI (i, j)] ==
	      M[TI (i - 1, j - 1)] +
	      (*S)[Q (i - 1)][Q (i)][Tc (j - 1)][Tc (j)])
	    {
	      state = MA;
	    }
	  else if (M[TI (i, j)] ==
		   Bq[TI (i - 1, j - 1)] +
		   (*S)[Q (i - 1)][Q (i)][GAP][Tc (j)])
	    {
	      state = BQ;

	    }
	  else if (M[TI (i, j)] ==
		   Bt[TI (i - 1, j - 1)] +
		   (*S)[GAP][Q (i)][Tc (j - 1)][Tc (j)])
	    {
	      state = BT;
	    }
	  i--;
	  j--;
	}
      else if (state == BQ)
	{
	  /* sprintf(buf, "%c", QS(i)) >> query_str ;  "-" >> target_str ; */
	  strcat (temp_ia_str, "Q");

	  if (Bq[TI (i, j)] ==
	      M[TI (i - 1, j)] + (*S)[Q (i - 1)][Q (i)][Tc (j)][GAP])
	    {
	      state = MA;
	      i--;
	    }
	  else if (Bq[TI (i, j)] ==
		   Bq[TI (i - 1, j)] + (*S)[Q (i - 1)][Q (i)][GAP][GAP])
	    {
	      state = BQ;
	      i--;
	    }
	}
      else if (state == BT)
	{
	  /* sprintf(buf, "%c", TSc(j)) >> target_str ;  "-" >> query_str ; */
	  strcat (temp_ia_str, "T");

	  if (Bt[TI (i, j)] ==
	      M[TI (i, j - 1)] + (*S)[Q (i)][GAP][Tc (j - 1)][Tc (j)])
	    {
	      state = MA;
	      j--;
	    }
	  else if (Bt[TI (i, j)] ==
		   Bt[TI (i, j - 1)] + (*S)[GAP][GAP][Tc (j - 1)][Tc (j)])
	    {
	      state = BT;
	      j--;
	    }
	}
    }

#undef Q
//#undef QS
#undef T
#undef Tc
//#undef TS
//#undef TSc
#undef TI

  str_rev (temp_ia_str);
  strcat (ia_str, temp_ia_str);
}

void
print_result (char *qname, const saidx64_t * sa, const saidx64_t * qsa,
	      aln_result_t * result)
{
  saidx64_t lpos, rpos;
  char strand;
  char *result_str = result->result_str;	//calloc(1024, sizeof(char));

  // iterate over the list of sequences to find out which one
  // contains the current result
  int imin = 0, imax = num_idxs, imid, idx = -1;
  while (imax >= imin)
    {
      imid = imin + ((imax - imin) / 2);

      //printf("imin: %d, imid: %d imax: %d sum_l[%d]: %lu sum_l[%d]: %lu\n", imin, imid, imax, imid, sum_l[imid], imid+1, sum_l[imid+1]);

      if (sum_l[imid] <= result->new_left &&
	  sum_l[imid + 1] > result->new_left)
	{
	  idx = imid;
	  break;
	}
      if (sum_l[imid + 1] <= result->new_left)
	imin = imid + 1;
      else if (sum_l[imid] > result->new_left)
	imax = imid - 1;
      else
	{
	  idx = imid;
	  break;
	}
    }

  if (idx == -1)
    {
      fprintf (stderr, "binary search failed\n");
      exit (1);
    }
  //fprintf(stderr, "idx: %d\n", idx);

  lpos = result->new_left - sum_l[idx];
  rpos = result->new_right - sum_l[idx];

  //printf("lpos: %lu rpos: %lu l[idx]: %lu idx: %d continuing\n", lpos, rpos, l[idx], idx);
  if (rpos > idx_lengths[idx])
    {
      return;
    }

  //fprintf(stderr, "lpos: %lu rpos: %lu \n", lpos, rpos);
  // odd-numbered sequences are the reverse-complement of the 
  // sequence before it
  if (idx % 2 == 1)
    {
      saidx64_t tpos = lpos;
      lpos = idx_lengths[idx] - rpos;
      rpos = idx_lengths[idx] - tpos;
      strand = '+';
    }
  else
    {
      strand = '-';
      lpos += 1;
      rpos += 1;
    }
  idx = (idx / 2) * 2;

  /*
     sprintf(result_str, "%s\t%s\t%c\t%lu\t%lu\t%lu\t%lu\t%.2f", 
     qname, name[idx],
     strand, lpos, rpos,
     result->q_new_left+1, result->q_new_right+1,
     result->energy);
   */
  sprintf (result_str, "%s\t%lu\t%lu\t%s\t%lu\t%lu\t%c\t%.2f",
	   qname, result->q_new_left + 1, result->q_new_right + 1,
	   name[idx], lpos, rpos, strand, result->energy);
}

void
really_print_alignment (gzFile qout, aln_result_t * result)
{
  // Actually print out the output strings stores in the result.
  // This is done so that all of the result strings are printed out
  // at once and the multi-threading code can print everythin in one
  // go

  if (show_alignment == 1)
    {
      gzprintf (qout, "%s\n", result->query_str);
      gzprintf (qout, "%s\n", result->aln_str);
      gzprintf (qout, "%s\n", result->target_str);
      gzprintf (qout, "%s\n", result->result_str);
    }
  else if (show_alignment == 2)
    {
      gzprintf (qout, "%s\t%s\n", result->result_str, result->ia_str);
    }
  else if (show_alignment == 3)
    {
      gzprintf (qout, "%s\t%s\t%s\t%s\t%s\n", result->result_str,
		result->ia_str, result->target_str, result->temp_target_str,
		result->temp_ia_str);
    }
  else
    {
      gzprintf (qout, "%s\n", result->result_str);
    }
}


/**
 * @brief Get final result after DP extension for maximal seed match
 *
 * Takes position of any seed match (within suffix array),
 * Tests whether seed match is maximal, by checking whether flanked to the left or right by a valid base pair,
 * Only performs DP extension if seed is maximal, because otherwise, it is contained into a maximal seed.
 * Finds boundaries for DP,
 * Calculates scores for left/right extension (calls)
 * Transform to energy and print results (call)
 * 
 * @return Void.
 */
void
extend_seed (saidx64_t j,	//!< Index in target suffix array (where seed match was found)
	     saidx64_t k,	//!< Index in query suffix array (where seed match was found)
	     const saidx64_t * sa,	//!< Suffix array of the target
	     const saidx64_t * qsa,	//!< Suffix array of the query
	     aln_result_t * result,	//!< result struct, fill in positions 
	     int seed_score,	//!< Score of the seed match
	     saidx64_t query_len,	//!< Length of the query, target length is in sum_l table
	     int seed_len,	//!< Length of the seed
	     saidx64_t best_spos,	//!< Start pos. of the most energy favourable subseed within any seed
	     int best_seed_len,	//!< Length of the most energy favourable subseed
	     int best_seed_score,	//!< Score of the most energy favourable subseed
	     char *qname,	//!< Name of the query sequence
	     int *seed,		//!< normalized (query-specific) seed info (m:n/l)
	     int *M_left,	//!< DP matrix for left extension, match state
	     int *Bq_left,	//!< DP matrix for left extension, query bulged state
	     int *Bt_left,	//!< DP matrix for left extension, target bulged state
	     int *M_right,	//!< DP matrix for right extension, match state
	     int *Bq_right,	//!< DP matrix for right extension, query bulged state
	     int *Bt_right)	//!< DP matrix for right extension, target bulged state
{
  int q_seed_len = seed_len, t_seed_len = seed_len;
  int total_score = seed_score;
  saidx64_t best_left_i = 0, best_left_j = 0, best_right_i = 0, best_right_j =
    0;

  result->found = 0;

  // get 0-based suffix numbers for query and target, and get seq id of the target region
  saidx64_t q_start = XSA (qsa[k]);
  saidx64_t t_start = XSA (sa[j]);
  saidx64_t idx = XIDX (sa[j]);

  //* get remaining lengths for left and right maximality check (min=0)
  int rem_query_len_left = q_start;
  int rem_target_len_left = t_start - sum_l[idx];
  int rem_query_len_right = query_len - (q_seed_len + q_start);
  int rem_target_len_right = sum_l[idx + 1] - (t_start + t_seed_len);

  if (print_debug)
    fprintf (stderr, "analyzing q_seed:%lu-%lu, t_start: %lu \n ",
	     q_start + 1, q_start + q_seed_len, t_start);

#define Q(ix) (XRIS(qsa[q_start+(ix)]))
#define T(ix) (comp[XRIS(sa[t_start+(ix)])])
  // test whether seed is flanked on left or right by a valid base pair
  // Maximal Test is different when seed is limited to some region
  // in that case, left maximality is important only if seed starts in mid positions. seed_flag determines if there is a positional constraint
  if (noGUseed)
    {
      if (seed_flag)
	{
	  if (q_start + 1 < seed[0] || q_start + q_seed_len > seed[1])
	    {
	      if (print_debug)
		{
		  fprintf (stderr,
			   "extension stopped due to start/end position for q_seed:%lu-%lu,  rem_query_len_left: %d, rem_query_len_right: %d, rem_target_len_left: %d, rem_target_len_right: %d \n",
			   q_start + 1, q_start + q_seed_len,
			   rem_query_len_left, rem_query_len_right,
			   rem_target_len_left, rem_target_len_right);
		}
	      // extension stops due to start or end position
	      return;
	    }

	  if ((rem_query_len_left - (seed[0] - 1) > 0
	       && rem_target_len_left > 0 && pair_noGU[Q (-1)][T (-1)])
	      || (rem_query_len_right - (query_len - seed[1]) > 0
		  && rem_target_len_right > 0
		  && pair_noGU[Q (q_seed_len)][T (t_seed_len)]))
	    {
	      if (print_debug)
		{
		  fprintf (stderr,
			   "extension stopped due to maximality for q_seed:%lu-%lu,  rem_query_len_left: %d, rem_query_len_right: %d, rem_target_len_left: %d, rem_target_len_right: %d \n",
			   q_start + 1, q_start + q_seed_len,
			   rem_query_len_left, rem_query_len_right,
			   rem_target_len_left, rem_target_len_right);
		}
	      // extension stops due to maximality within the seed
	      return;
	    }
	}
      else
	if ((rem_query_len_left > 0 && rem_target_len_left > 0
	     && pair_noGU[Q (-1)][T (-1)]) || (rem_query_len_right > 0
					       && rem_target_len_right > 0
					       &&
					       pair_noGU[Q (q_seed_len)][T
									 (t_seed_len)]))
	{
	  if (print_debug)
	    {
	      fprintf (stderr,
		       "stop extension whether maximality or str/end position! %lu %lu %d %d\n",
		       q_start, t_start, q_seed_len, t_seed_len);
	    }
	  return;
	}
    }
  else
    {
      if (seed_flag)
	{
	  if (q_start + 1 < seed[0] || q_start + q_seed_len > seed[1])
	    {
	      if (print_debug)
		{
		  fprintf (stderr,
			   "extension stopped due to start/end position for q_seed:%lu-%lu,  rem_query_len_left: %d, rem_query_len_right: %d, rem_target_len_left: %d, rem_target_len_right: %d \n",
			   q_start + 1, q_start + q_seed_len,
			   rem_query_len_left, rem_query_len_right,
			   rem_target_len_left, rem_target_len_right);
		}
	      // extension stops due to start or end position
	      return;
	    }

	  if ((rem_query_len_left - (seed[0] - 1) > 0
	       && rem_target_len_left > 0 && pair[Q (-1)][T (-1)])
	      || (rem_query_len_right - (query_len - seed[1]) > 0
		  && rem_target_len_right > 0
		  && pair[Q (q_seed_len)][T (t_seed_len)]))
	    {
	      if (print_debug)
		{
		  fprintf (stderr,
			   "extension stopped due to maximality for q_seed:%lu-%lu,  rem_query_len_left: %d, rem_query_len_right: %d, rem_target_len_left: %d, rem_target_len_right: %d \n",
			   q_start + 1, q_start + q_seed_len,
			   rem_query_len_left, rem_query_len_right,
			   rem_target_len_left, rem_target_len_right);
		}
	      // extension stops due to maximality within the seed
	      return;
	    }
	}
      else
	if ((rem_query_len_left > 0 && rem_target_len_left > 0
	     && pair[Q (-1)][T (-1)]) || (rem_query_len_right > 0
					  && rem_target_len_right > 0
					  &&
					  pair[Q (q_seed_len)][T
							       (t_seed_len)]))
	{
	  if (print_debug)
	    {
	      fprintf (stderr,
		       "stop extension whether maximality or str/end position! %lu %lu %d %d\n",
		       q_start, t_start, q_seed_len, t_seed_len);
	    }
	  return;
	}
    }
#undef Q
#undef T


  // WE DO NOT EXTEND THE MAXIMAL SEEDS IF there is a seed threshold
  // INSTEAD WE EXTEND THE MOST FAVOURABLE SUBSEED WITHIN MAXIMAL SEED
  if (seed_threshold_flag)
    {
      seed_len = best_seed_len;
      q_seed_len = best_seed_len;
      t_seed_len = best_seed_len;
      q_start += best_spos;
      t_start += best_spos;
      total_score = best_seed_score;
    }

  //Remaining sequence length definition has to change here since we are going to decide DP boundaries
  /* find boundaries for DP (remaining seq left/right of the seed)
   * all rem_{query,target}_len_{left,right} are + 1
   * i.e. nnnNNNNNnnnnn (with N for seed location) yields rem..left=4, rem..right=6
   * find value for query by actual boundaries (limited by max_ext_len)
   * for target its also limited by max_ext_len, possibly corrected by sequence boundaries
   */
  rem_query_len_left = min (max_ext_len, q_start + 1);
  rem_target_len_left = min (max_ext_len, t_start - sum_l[idx] + 1);
  rem_query_len_right =
    min (max_ext_len, query_len - q_seed_len - q_start + 1);
  rem_target_len_right =
    min (max_ext_len, sum_l[idx + 1] - (t_start + t_seed_len - 1));

  if (print_debug)
    {
      fprintf (stderr, "seed_score = seed_score : %d\n", total_score);
    }

  total_score +=
    DP_left (M_left, Bq_left, Bt_left, rem_query_len_left,
	     rem_target_len_left, qsa, sa, q_start, t_start, &best_left_i,
	     &best_left_j);

  if (print_debug)
    {
      fprintf (stderr, "extending the seed...\n");
      fprintf (stderr, "total_score += DP_left : %d\n", total_score);
    }

  total_score +=
    DP_right (M_right, Bq_right, Bt_right, rem_query_len_right,
	      rem_target_len_right, qsa, sa, q_start + q_seed_len - 1,
	      t_start + t_seed_len - 1, &best_right_i, &best_right_j);

  if (print_debug)
    {
      fprintf (stderr, "total_score += DP_right : %d\n", total_score);
    }

  int nt_count =
    best_left_j + best_left_i + best_right_i + best_right_j + q_seed_len +
    t_seed_len;

  if (print_debug)
    {
      fprintf (stderr,
	       "extend_seed lq_left: %d lt_left: %d lq_right: %d lt_right: %d\n",
	       rem_query_len_left, rem_target_len_left, rem_query_len_right,
	       rem_target_len_right);
      fprintf (stderr,
	       "nt_count is %d, using blj %li, bli %li, bri %li, brj %li, q_seedlen: %d, t_seedlen: %d \n",
	       nt_count, best_left_j, best_left_i, best_right_i, best_right_j,
	       q_seed_len, t_seed_len);
      fprintf (stderr, "**********************\n");
    }

  result->score = total_score + nt_count * extPen;
  /*Energy correction, BE AWARE of energy correction, if one energy matrix is accepted as input there will be some option here */
  result->energy = (result->score - 559.0f) / -100.0f;

  //ENERGY THRESHOLD for the interaction
  if (result->energy > min_energy)
    {
      if (print_debug)
	{
	  fprintf (stderr, "seed_score: %d \n", seed_score);
	  fprintf (stderr, "energy for this int: %f, thr: %f \n",
		   result->energy, min_energy);
	}
      return;
    }

  result->best_left_i = best_left_i;
  result->best_left_j = best_left_j;

  result->best_right_i = best_right_i;
  result->best_right_j = best_right_j;

  //exit(1);


  {
    if (show_alignment == 1)
      {
	print_alignment (sa, qsa, result, query_len, seed_len, M_left,
			 Bq_left, Bt_left, M_right, Bq_right, Bt_right);
      }
    else if (show_alignment == 2)
      {
	print_ali2 (sa, qsa, result, query_len, seed_len, M_left, Bq_left,
		    Bt_left, M_right, Bq_right, Bt_right);
      }
    else if (show_alignment == 3)
      {
	print_alignment (sa, qsa, result, query_len, seed_len, M_left,
			 Bq_left, Bt_left, M_right, Bq_right, Bt_right);
	print_ali2 (sa, qsa, result, query_len, seed_len, M_left, Bq_left,
		    Bt_left, M_right, Bq_right, Bt_right);

	/* Lets get the target & PAM sequence from target */
	char buf[4];
	/* Save the flanking sequences on already created temp strings */
	char *temp_fivend_str = result->temp_target_str;
	char *temp_thrend_str = result->temp_ia_str;
	*temp_fivend_str = 0;
	*temp_thrend_str = 0;

	saidx64_t lpos = result->new_left - sum_l[idx];
	saidx64_t rpos = result->new_right - sum_l[idx];
	if (idx % 2 == 1)
	  {
	    saidx64_t tpos = lpos;
	    lpos = idx_lengths[idx] - rpos;
	    rpos = idx_lengths[idx] - tpos;
	  }
	else
	  {
	    lpos += 1;
	    rpos += 1;
	  }
	saidx64_t trailing_left = lpos - 1;
	saidx64_t trailing_right = idx_lengths[idx] - rpos;

	// Get the flanking 5'end sequence of target
#define TSc(ix) (scomp[XSTR(sa[result->new_right+(ix+1)])])
	saidx64_t flank_length = 20;

	saidx64_t valid_length = min (flank_length, trailing_right);
	if (idx % 2 == 1)
	  valid_length = min (flank_length, trailing_left);

	saidx64_t i;
	for (i = 0; i < valid_length; i = i + 1)
	  {
	    sprintf (buf, "%c", TSc (i));
	    strcat (temp_fivend_str, buf);
	  }
#undef TSc

	// Get the flanking 3'end sequence of target
#define TSc(ix) (scomp[XSTR(sa[result->new_left-(ix+1)])])
	valid_length = min (flank_length, trailing_left);

	if (idx % 2 == 1)
	  valid_length = min (flank_length, trailing_right);

	for (i = 0; i < valid_length; i = i + 1)
	  {
	    sprintf (buf, "%c", TSc (i));
	    strcat (temp_thrend_str, buf);
	  }
#undef TSc

      }
    else
      {
	/* skip trackback, but still need to get the positions */
	q_start = XSA (qsa[result->k]);
	t_start = XSA (sa[result->j]);
	result->new_left = t_start - result->best_left_j;
	result->q_new_left = q_start - result->best_left_i;
	t_start = t_start + t_seed_len - 1;
	q_start = q_start + q_seed_len - 1;
	result->new_right = t_start + result->best_right_j;
	result->q_new_right = q_start + result->best_right_i;

      }
    result->found = 1;
    print_result (qname, sa, qsa, result);
  }

}


/**@brief Post-process intervals with risearch algorithm 
 * @return Void.
*/

void
sa_evaluate_interval (sa_interval_list_t * results, query_t * query,
		      const saidx64_t * sa, char strand)
{
  sa_interval_list_t *p;
  int seed_len, seed_score;
  saidx64_t ql, qr, sl, sr, j, k;
  saidx64_t q_start, t_start, xsa, idx;
  saidx64_t *qsa = query->qsa;
  int query_len = query->length;
  int dp_size, str_size;

  //For finding best subseed and seed energy filtering
  saidx64_t min_seed_length = query->seed[2], best_subseed_spos;;
  double best_perlength_score, temp_perlength_score;
  int temp_subseed_score, temp_subseed_len, best_subseed_len,
    best_subseed_score;

#ifdef DEBUG
  print_debug = 1;
#endif

#if VERBOSE>1
  print_mats_l = 1;
  print_mats_r = 1;
#endif

  if (max_ext_len > 0)
    {
      dp_size = (max_ext_len + 1) * (max_ext_len + 1);
      str_size = 4 * max_ext_len + query_len;	//should be enough to allow for x gaps on x nts extension
    }
  else
    {
      dp_size = (query_len) * (query_len);
      str_size = 2 * (query_len);
    }

  // the length of the dp matrix will be the length of the extension window
  // tables for DP to the left and right of the seed region
  int *M_left = (int *) calloc (dp_size, sizeof (int));
  int *Bq_left = (int *) calloc (dp_size, sizeof (int));
  int *Bt_left = (int *) calloc (dp_size, sizeof (int));

  int *M_right = (int *) calloc (dp_size, sizeof (int));
  int *Bq_right = (int *) calloc (dp_size, sizeof (int));
  int *Bt_right = (int *) calloc (dp_size, sizeof (int));

  if (!M_left || !Bq_left || !Bt_left || !M_right || !Bq_right || !Bt_right)
    {
      fprintf (stderr, "Failed allocating DP matrices\n");
      exit (EXIT_FAILURE);
    }

  //print_dsm_table();
  aln_result_t *result = (aln_result_t *) calloc (1, sizeof (aln_result_t));
  if (!result)
    {
      fprintf (stderr, "Failed allocating result struct\n");
      exit (EXIT_FAILURE);
    }
  if (show_alignment == 1)
    {
      result->query_str = calloc (str_size, sizeof (char));
      result->temp_query_str = calloc (str_size, sizeof (char));
      result->aln_str = calloc (str_size, sizeof (char));
      result->temp_aln_str = calloc (str_size, sizeof (char));
      result->target_str = calloc (str_size, sizeof (char));
      result->temp_target_str = calloc (str_size, sizeof (char));
      if (!result->query_str || !result->temp_query_str || !result->aln_str
	  || !result->temp_aln_str || !result->target_str
	  || !result->temp_target_str)
	{
	  fprintf (stderr, "Failed allocating result attributes\n");
	  exit (EXIT_FAILURE);
	}
    }
  else if (show_alignment == 2)
    {
      result->ia_str = calloc (str_size, sizeof (char));
      result->temp_ia_str = calloc (str_size, sizeof (char));
      if (!result->ia_str || !result->temp_ia_str)
	{
	  fprintf (stderr, "Failed allocating result attributes\n");
	  exit (EXIT_FAILURE);
	}
    }
  else if (show_alignment == 3)
    {
      result->query_str = calloc (str_size, sizeof (char));
      result->temp_query_str = calloc (str_size, sizeof (char));
      result->aln_str = calloc (str_size, sizeof (char));
      result->temp_aln_str = calloc (str_size, sizeof (char));
      result->target_str = calloc (str_size, sizeof (char));
      result->temp_target_str = calloc (str_size, sizeof (char));
      result->ia_str = calloc (str_size, sizeof (char));
      result->temp_ia_str = calloc (str_size, sizeof (char));
      if (!result->query_str || !result->temp_query_str || !result->aln_str
	  || !result->temp_aln_str || !result->target_str
	  || !result->temp_target_str || !result->ia_str
	  || !result->temp_ia_str)
	{
	  fprintf (stderr, "Failed allocating result attributes\n");
	  exit (EXIT_FAILURE);
	}
    }
  result->result_str = calloc (1024, sizeof (char));
  if (!result->result_str)
    {
      fprintf (stderr, "Failed allocating result string\n");
      exit (EXIT_FAILURE);
    }

  for (p = results; p; p = p->next)
    {
      // the locations that match in the index
      sl = p->start;
      sr = p->end;

      // the locations that match in the query
      ql = p->qstart;
      qr = p->qend;

      seed_len = p->slen;
      seed_score = 0;

      q_start = XSA (qsa[ql]);
      t_start = XSA (sa[sl]);

      idx = XIDX (sa[sl]);
      //#define Q(ix) (comp[XRIS(qsa[q_start+(ix)])])
#define Q(ix) (XRIS(qsa[q_start+(ix)]))
      //#define T(ix) (XRIS(sa[t_start+(ix)]))
#define T(ix) (comp[XRIS(sa[t_start+(ix)])])
      //fprintf(stderr, "orig q_start: %lu t_start: %lu\n", q_start, t_start);


      for (j = 0; j < seed_len - 1; j++)
	{
	  seed_score += (*S)[Q (j)][Q (j + 1)][T (j)][T (j + 1)];
	}

      //Minimum energy per length threshold for seeds
      best_perlength_score = (double) seed_score / seed_len;
      best_subseed_len = seed_len;
      best_subseed_score = seed_score;
      best_subseed_spos = 0;
      if (seed_threshold_flag)
	{
	  double threshold = min_seed_energy_per_length * (-100);
	  int hold_score;
	  temp_subseed_score = seed_score;
	  for (j = seed_len - 1; j >= (min_seed_length - 1); j--)
	    {
	      hold_score = temp_subseed_score;
	      for (k = 0; k <= (j - (min_seed_length - 1)); k++)
		{
		  temp_subseed_len = (j - k) + 1;
		  temp_perlength_score =
		    (double) temp_subseed_score / temp_subseed_len;
		  if (temp_perlength_score > best_perlength_score)
		    {
		      best_perlength_score = temp_perlength_score;
		      best_subseed_len = temp_subseed_len;
		      best_subseed_score = temp_subseed_score;
		      best_subseed_spos = k;
		    }
		  temp_subseed_score -=
		    (*S)[Q (k)][Q (k + 1)][T (k)][T (k + 1)];
		}
	      temp_subseed_score =
		hold_score - (*S)[Q (j - 1)][Q (j)][T (j - 1)][T (j)];
	    }
	  if (best_perlength_score < threshold)
	    {
	      continue;
	    }
	}

#undef Q
#undef T
      //printf("\n sl = %d   seed_len = %d \n ", sl, seed_len); 

      for (j = sl; j < sr; j++)
	{
	  xsa = XSA (sa[j]);	// = t_start
	  idx = XIDX (sa[j]);

	  if (xsa + seed_len > sum_l[idx + 1])
	    {
	      //fprintf(stderr, "seed match extends into next sequence, aborting\n");
	      //Since all valid seeds, whether maximal or not, are feeded into this stage
	      //Ignoring seeds that extends into next sequence is not a problem
	      continue;
	    }

	  for (k = ql; k < qr; k++)
	    {
	      result->j = j;
	      result->k = k;
	      //debug("\n seed len = %d, query len = %d, seed score = %d , t_start = %d, q_start = %d \n", seed_len, query_len, seed_score, XSA(sa[j]), XSA(qsa[k]));
	      if (print_debug)
		{
		  fprintf (stderr,
			   "\nseed len = %d, query len = %d, seed score = %d , t_start = %lu, q_start = %lu \n",
			   seed_len, query_len, seed_score, XSA (sa[j]),
			   XSA (qsa[k]));
		  fprintf (stderr, "extend seed: %d \n", seed_len);
		}
	      //Extend the current seed, (subseed info is given seperately)
	      extend_seed (j, k, sa, qsa, result, seed_score, query->length,
			   seed_len, best_subseed_spos, best_subseed_len,
			   best_subseed_score, query->name, query->seed,
			   M_left, Bq_left, Bt_left, M_right, Bq_right,
			   Bt_right);

//#pragma omp critical
//NOT critical as now all write to own file! reduces overhead that is spent waiting?
	      {
		if (result->found)
		  really_print_alignment (query->out, result);
	      }
	    }
	}
    }


  free (M_right);
  free (Bq_right);
  free (Bt_right);

  free (M_left);
  free (Bq_left);
  free (Bt_left);
  if (show_alignment == 1)
    {
      free (result->query_str);
      free (result->temp_query_str);
      free (result->aln_str);
      free (result->temp_aln_str);
      free (result->target_str);
      free (result->temp_target_str);
    }
  else if (show_alignment == 2)
    {
      free (result->ia_str);
      free (result->temp_ia_str);
    }
  else if (show_alignment == 3)
    {
      free (result->query_str);
      free (result->temp_query_str);
      free (result->aln_str);
      free (result->temp_aln_str);
      free (result->target_str);
      free (result->temp_target_str);
      free (result->ia_str);
      free (result->temp_ia_str);
    }
  free (result->result_str);
  free (result);
}

/** @brief parallel suffix array matching on sense-strand
 *
 * Recursive function to match query and target suffix in parallel.
 * Steps into suffixes until one of the trees has all its nodes visited.
 * Collects quadrupels (sl,sr,ql,qr) that define all seed match intervals, when the minimum seed length is reached.
 *  
 * @return Void. 
 */

void
sa_parallel_match_neg (const saidx64_t * qsa,	//!< Suffix array of the query
		       saidx64_t ql,	//!< The first position (in query SA) of the interval to examine/step into
		       saidx64_t qr,	//!< The last position (exclusive) of that interval
		       const saidx64_t * sa,	//!< Suffix array of the target
		       saidx64_t sl,	//!< The first position (in target SA) of the interval to examine/step into
		       saidx64_t sr,	//!< The last position (exclusive) of that interval
		       saidx64_t offset,	//!< The number of steps already made into the suffix
		       saidx64_t depth,	//!< The required minimum seed length
		       int last_match_count,	//!< Current number of matches after the last mismatch.
		       int mismatch_count,	//!< The number of mismatches in the current seed.
		       sa_interval_list_t ** results)	//!< List of intervals with seed matches
{
  saidx64_t qint[6], sint[6];

  if (depth <= offset)
    {
      // MISMATCH UPDATE
      if ((mismatch_count == 0)
	  || ((mismatch_count > 0) && (mismatch_count <= seed_mismatch[0])
	      && (last_match_count >= seed_mismatch[2])
	      && (last_match_count < depth)))
	{
	  sa_interval_add (results, sl, sr, ql, qr, offset);
	}
    }

  if (!sa_search_interval (qsa, ql, qr, offset, qint)
      || !sa_search_interval (sa, sl, sr, offset, sint))
    return;

  ++offset;

  /*
     a 1-0
     c 2-1
     g 3-2
     u 5-4
   */

  /* a */
  if (qint[1] - qint[0])
    {
      if (sint[1] - sint[0])
	sa_parallel_match_neg (qsa, qint[0], qint[1], sa, sint[0], sint[1],
			       offset, depth, last_match_count + 1,
			       mismatch_count, results);
    }
  /* c */
  if (qint[2] - qint[1])
    {
      if (sint[2] - sint[1])
	sa_parallel_match_neg (qsa, qint[1], qint[2], sa, sint[1], sint[2],
			       offset, depth, last_match_count + 1,
			       mismatch_count, results);
    }
  /*  g */
  if (qint[3] - qint[2])
    {
      if (sint[3] - sint[2])
	sa_parallel_match_neg (qsa, qint[2], qint[3], sa, sint[2], sint[3],
			       offset, depth, last_match_count + 1,
			       mismatch_count, results);
      /* wobble pair */
      /* g - comp('a') We need to match the complement of U: A */
      if (noGUseed == 0)
	if (sint[1] - sint[0])
	  sa_parallel_match_neg (qsa, qint[2], qint[3], sa, sint[0], sint[1],
				 offset, depth, last_match_count + 1,
				 mismatch_count, results);
    }

  /* u */
  if (qint[5] - qint[4])
    {
      if (sint[5] - sint[4])
	sa_parallel_match_neg (qsa, qint[4], qint[5], sa, sint[4], sint[5],
			       offset, depth, last_match_count + 1,
			       mismatch_count, results);
      /* wobble pair */
      /* source should be C because we're matching the complementary strand */
      if (noGUseed == 0)
	if (sint[2] - sint[1])
	  sa_parallel_match_neg (qsa, qint[4], qint[5], sa, sint[1], sint[2],
				 offset, depth, last_match_count + 1,
				 mismatch_count, results);
    }


	/********************** MISMATCH CASES ***********************/
  /* MISMATCH can only start in specified position and can have only specific length */
  if (seed_mismatch_flag)
    {

      if ((mismatch_count < seed_mismatch[0]) && (offset > seed_mismatch[1])
	  && (last_match_count < depth))
	{
	  /* MISMATCH for a */
	  if (qint[1] - qint[0])
	    {
	      if (sint[2] - sint[1])
		sa_parallel_match_neg (qsa, qint[0], qint[1], sa, sint[1],
				       sint[2], offset, depth, 0,
				       mismatch_count + 1, results);
	      if (sint[3] - sint[2])
		sa_parallel_match_neg (qsa, qint[0], qint[1], sa, sint[2],
				       sint[3], offset, depth, 0,
				       mismatch_count + 1, results);
	      if (sint[5] - sint[4])
		sa_parallel_match_neg (qsa, qint[0], qint[1], sa, sint[4],
				       sint[5], offset, depth, 0,
				       mismatch_count + 1, results);
	    }
	  /* MISMATCH for c */
	  if (qint[2] - qint[1])
	    {
	      if (sint[1] - sint[0])
		sa_parallel_match_neg (qsa, qint[1], qint[2], sa, sint[0],
				       sint[1], offset, depth, 0,
				       mismatch_count + 1, results);
	      if (sint[3] - sint[2])
		sa_parallel_match_neg (qsa, qint[1], qint[2], sa, sint[2],
				       sint[3], offset, depth, 0,
				       mismatch_count + 1, results);
	      if (sint[5] - sint[4])
		sa_parallel_match_neg (qsa, qint[1], qint[2], sa, sint[4],
				       sint[5], offset, depth, 0,
				       mismatch_count + 1, results);
	    }

	  /* MISMATCH for g */
	  if (qint[3] - qint[2])
	    {
	      if (sint[2] - sint[1])
		sa_parallel_match_neg (qsa, qint[2], qint[3], sa, sint[1],
				       sint[2], offset, depth, 0,
				       mismatch_count + 1, results);
	      if (sint[5] - sint[4])
		sa_parallel_match_neg (qsa, qint[2], qint[3], sa, sint[4],
				       sint[5], offset, depth, 0,
				       mismatch_count + 1, results);
	      if (noGUseed)
		if (sint[1] - sint[0])
		  sa_parallel_match_neg (qsa, qint[2], qint[3], sa, sint[0],
					 sint[1], offset, depth, 0,
					 mismatch_count + 1, results);

	    }

	  /* MISMATCH for u */
	  if (qint[5] - qint[4])
	    {
	      if (sint[1] - sint[0])
		sa_parallel_match_neg (qsa, qint[4], qint[5], sa, sint[0],
				       sint[1], offset, depth, 0,
				       mismatch_count + 1, results);
	      if (sint[3] - sint[2])
		sa_parallel_match_neg (qsa, qint[4], qint[5], sa, sint[2],
				       sint[3], offset, depth, 0,
				       mismatch_count + 1, results);
	      if (noGUseed)
		if (sint[2] - sint[1])
		  sa_parallel_match_neg (qsa, qint[4], qint[5], sa, sint[1],
					 sint[2], offset, depth, 0,
					 mismatch_count + 1, results);

	    }
	}
    }

}
