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
#include <ctype.h>

#include "main.h"
#include "fasta.h"
#include "sa.h"

/* maps nucleotides to 4bit representation */
saidx64_t sa_strmap(sauchar_t c)
{
	saidx64_t ret = 0;
	switch((char)c) {
		/* mappings 0x06 + 0x07 => 0x0E + 0x0F free for use */
		case  0 : ret = 0x00; break;
		case 'a': ret = 0x01; break;
		case 'g': ret = 0x02; break;
		case 'c': ret = 0x03; break;
		case 't': ret = 0x04; break;
		case 'u': ret = 0x04; break;
		case 'A': ret = 0x09; break;
		case 'G': ret = 0x0A; break;
		case 'C': ret = 0x0B; break;
		case 'T': ret = 0x0C; break;
		case 'U': ret = 0x0C; break;
		case 'X': ret = 0x0D; break;
		default : ret = 0x05; break;
	}

	return (ret<<34);
}

/* [Toru Kasai et al] Linear-Time Longest-Common-Prefix Computation in Suffix Arrays and Its Applications */
void sa_lcp(const sauchar_t *str, const saidx64_t *sa, uint8_t *lcp, saidx64_t *rank, saidx64_t n)
{
	saidx64_t h = 0, i, j, k; lcp[0] = 0;

	for(i = 0; i < n; i++)
		rank[sa[i]] = i;
	for(i = 0; i < n; i++) {
		k = rank[i];
		if(k > 0) {
			j = sa[k-1];
			for(;str[i+h] && str[j+h] && str[i+h] == str[j+h]; h++);
			/* lcp table is truncated to 8bits */
			lcp[k] = (h>255?255:(uint8_t)h);
			if(h > 0)
				h--;
		}
	}
}

void sa_length_to_index(saidx64_t **l, saidx64_t k)
{
	saidx64_t i,j,x;
	for(j = (*l)[0]+1, (*l)[0] = 1, i = 1; i <= k; i++) { 
		x = j+1;
		j += (*l)[i]+1;
		(*l)[i] = x;
	}
}

/* read suffix array from file */
void sa_read(const char *filename, saidx64_t **sa, saidx64_t *n, saidx64_t **l, saidx64_t **sum_l, saidx64_t *k, char ***name)
{
	FILE *in;
	saidx64_t i;
	saidx64_t tmp;
	uint16_t j;

	in = strcmp(filename,"-") ? fopen(filename, "rb") : stdin;

	if (!in) {
		fprintf(stderr, "file not found: %s\n", filename);
		exit(1);
	}

	if(1 != fread(n, sizeof(saidx64_t), 1, in)) {
		fprintf(stderr, "Erroneous index file given!\n");
		exit(1);
	}
	if(1 != fread(k, sizeof(saidx64_t), 1, in)) {
		fprintf(stderr, "Erroneous index file given!\n");
		exit(1);
	}
//	fprintf(stderr, "read : n=%"PRIdSAIDX64_T", k=%"PRIdSAIDX64_T"\n", *n, *k);

	*l = (saidx64_t*)malloc((*k)*sizeof(saidx64_t));
	*sum_l = (saidx64_t*)malloc((*k+1)*sizeof(saidx64_t));
	*name = (char**)malloc((*k)*sizeof(char*));
	if ((*l == NULL) || (*sum_l == NULL) || (*name == NULL)) {
		fprintf(stderr, "Failed allocating suffix array A\n");
		exit(EXIT_FAILURE);
	}

	tmp=fread(*l, sizeof(saidx64_t), *k, in);
	// check nitems is k , or feof / ferror !!! 
	if (tmp != *k) {
		fprintf(stderr, "Erroneous index file given\n");
//		fprintf(stderr, "read l: %"PRIdSAIDX64_T" not %"PRIdSAIDX64_T"\n", tmp,*k);
		exit(1);
	}

	(*sum_l)[0] = 0;

	for(i = 0; i < *k; i++) {
		if (1 != fread(&j, sizeof(uint16_t), 1, in)) {
			fprintf(stderr, "Erroneous index file given!\n");
//			fprintf(stderr, "no j at index %"PRIdSAIDX64_T"\n",i);
			exit(1);
		}
		(*name)[i] = (char*)malloc((j+1)*sizeof(char));
		if ((*name)[i] == NULL) {
			fprintf(stderr, "Failed allocating suffix array B\n");
			exit(EXIT_FAILURE);
		}
		if (j+1 != fread((*name)[i], sizeof(char), j+1, in)) {
			fprintf(stderr, "Erroneous index file given!\n");
//			fprintf(stderr, "no name of length %"PRIu16"  at index %"PRIdSAIDX64_T"\n",j+1, i);
			exit(EXIT_FAILURE);
		}
		(*sum_l)[i+1] = (*sum_l)[i] + (*l)[i];
//		fprintf(stderr, "read : l=%"PRIdSAIDX64_T", suml=%"PRIdSAIDX64_T", name=%s\n", (*l)[i], (*sum_l)[i+1], (*name)[i]);
	}

	*sa = (saidx64_t*)malloc((*n +1)*sizeof(saidx64_t));
	if (*sa == NULL) {
		fprintf(stderr, "Failed allocating suffix array C\n");
		exit(EXIT_FAILURE);
	}

	tmp = fread(*sa, sizeof(saidx64_t), *n, in);
	if (tmp != *n) {
		fprintf(stderr, "Erroneous index file given!\n");
//		fprintf(stderr, "read sa: %"PRIdSAIDX64_T" not %"PRIdSAIDX64_T"\n", tmp,*n);
		exit(1);
	}
	(*sa)[*n] = 0; /* preserve from reading unitilized values in sa_search_left */
	fclose(in);
}

/* write suffix array to file */
void sa_write(const char *filename, saidx64_t *sa, saidx64_t n, saidx64_t *l, saidx64_t k, char **name)
{
	saidx64_t i;
	uint16_t j;
	FILE *out;

	out = strcmp(filename,"-") ? fopen(filename, "wb") : stdout;

	if (out == NULL) {
		fprintf(stderr, "Failed writing to file %s, please make sure you have write permissions\n", filename);
		exit(EXIT_FAILURE);
	}

	fwrite(&n, sizeof(saidx64_t), 1, out);
	fwrite(&k, sizeof(saidx64_t), 1, out);
	fwrite(l, sizeof(saidx64_t), k, out);

	for(i = 0; i < k; i++) {
		j = strlen(name[i]);
		fwrite(&j, sizeof(uint16_t), 1, out);
		fwrite(name[i], sizeof(char), j+1, out);
	}

	fwrite(sa, sizeof(saidx64_t), n, out);

	fclose(out);
}

/* merge all tables into packed suffix array */
void sa_merge(const sauchar_t *str, saidx64_t *sa, saidx64_t n, saidx64_t *sum_l, saidx64_t k)
{
	saidx64_t i, imin, imax, imid, idx, tmp;

	for(i = 0; i < n; i++) {
		tmp = sa[i];
		imin=0;
		imax=k;
		idx=-1;
		while (imax >= imin) {
			imid = imin + ((imax - imin) / 2);
		//	if (sum_l[imid] <= tmp && sum_l[imid+1] > tmp) {
		//		idx = imid;
		//		break;
		//	}
			if (sum_l[imid+1] <= tmp)
				imin = imid + 1;
			else if (sum_l[imid] > tmp)
				imax = imid - 1;
			else {
				idx = imid;
				break;
			}
		}
		if (idx == -1) {
			fprintf(stderr, "binary search failed\n");
			exit(1);
		}
		//fprintf(stderr, "%ld\n", idx);
		sa[i] = MKSA(tmp) | MKSTR(str[i]) | MKIDX(idx);
	}
}

/* print partial sequence from suffix array */
void sa_debug_partial(const saidx64_t *sa, saidx64_t i, saidx64_t n)
{
	char c;
	for(; n<0 ? XSTR(sa[i]) : i<n; i++) { 
		c=XSTR(sa[i]); 
		fprintf(stderr, "%c", c ? c : '\n');
	} 
	fprintf(stderr, "\n");
}

/* print suffix array contents */
void sa_debug(const saidx64_t *sa, saidx64_t n)
{
	saidx64_t i;
	sa_debug_partial(sa, 1, n);
	for(i = 0; i < n; i++) {
		fprintf(stderr, "%.2li\t%.2u\t%.2u\t", XSA(sa[i]), XLCP(sa[i]), XIDX(sa[i]));
		sa_debug_partial(sa, XSA(sa[i]), -1);
	}
}

/* complement rna sequence */
char rna_compc(char c)
{
	switch(c) {
		case 'a': return 'u';
		case 'c': return 'g';
		case 'g': return 'c';
		case 't': return 'a';
		case 'u': return 'a';
		case 'A': return 'U';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		case 'U': return 'A';
		default : return c;
	}
}

/* reverse rna sequence */
void rna_reverse(char *s, saidx64_t n)
{
	saidx64_t i, n2;
	for(i = 0, n2 = n>>1; i < n2; i++) { 
		char c = s[i];
		s[i] = s[n-i-1];
		s[n-i-1] = c;
	}
}

/* complement rna sequence */
void rna_comp(char *s)
{
	for(; *s; s++) { 
		*s = rna_compc(*s);
	};
}

/* resort suffix array by index */
/* compare function for qsort */
int sa_sortidx(const void *v1, const void *v2)
{
	saidx64_t *s1, *s2;
	s1 = (saidx64_t*)v1;
	s2 = (saidx64_t*)v2;

	return (XIDX(*s1) - XIDX(*s2));
}

/* create partial suffix array in memory */
void sa_create_partial_reverse(const char *seq, saidx64_t start, saidx64_t end, saidx64_t cnt, saidx64_t **sa, saidx64_t *n)
{
	char *t;
	saidx64_t i, sai;

	*sa = (saidx64_t*)malloc((*n +1) * sizeof(saidx64_t));
	t = (char*)malloc((*n +1) * sizeof(char));
	if ((*sa == NULL) || (t == NULL)) {
		fprintf(stderr, "Failed allocating partial suffix array\n");
		exit(EXIT_FAILURE);
	}
	
	/* use lower case representation for lexicographic order and seed matching later */
	for (i = 0; i < *n; i++)
		t[i] = tolower(seq[i]);

	if ( divsufsort64((sauchar_t*)t, (*sa), *n) != 0) {
		fprintf(stderr, "Failed generating partial suffix array\n");
		exit(EXIT_FAILURE);
	}

	//fprintf(stderr, "i\tsa[i]\tt[i]\tseq[i]\n");
	//for (i = 0; i < *n; i++) fprintf(stderr, "%lu\t%lu\t%c\t%c\n", i, (*sa)[i], t[i], seq[i]);

	for(i = 0; i < *n; i++) {
		sai = (*sa)[i];	/* sai is position within original sequence, 0-based */
		if (sai < start || sai > (end-cnt+1)) {
			/* this suffix cannot fulfil seed criteria */
			(*sa)[i] = MKSA(*n);
			(*sa)[i] |= MKIDX(0);
		} else { /*valid suffix , set index for resorting (offset, as 0 is for invalid)*/
			(*sa)[i] = MKSA(sai);
			(*sa)[i] |= MKIDX(i+1);
		}
	}
	/*sort according to IDX*/
	qsort(*sa, *n, sizeof(saidx64_t), &sa_sortidx);
	/*nt as in original sequence*/
	for(i = 0; i < *n; i++) {
		(*sa)[i] |= MKSTR((sauchar_t)seq[i]);
	}

	(*sa)[*n] = 0; /* preserve from reading unitilized values in sa_search_left */
	free(t);

	//for (i = 0; i <= *n; i++) fprintf(stderr, "XSA(sa[i=%lu]) = %lu; XSTR(sa[i])=%c; XIDX(sa[i])=%u;\n", i, XSA((*sa)[i]), XSTR((*sa)[i]), XIDX((*sa)[i]) );
	//sa_debug(*sa, *n);
}

/* create partial suffix array in memory for reverse complement */
void sa_create_partial_revcomp(const char *seq, saidx64_t start, saidx64_t end, saidx64_t cnt, saidx64_t **sa, saidx64_t *n)
{
	char *s = strdup(seq);
	*n = strlen(s);
	rna_reverse(s, *n);
	//rna_comp(s);
	sa_create_partial_reverse(s, *n-end+2, *n-start+2, cnt, sa, n);
	free(s);
}

/* convert string to lower case */
void strtolower(char *s)
{
	for(; *s; s++) 
		*s = tolower(*s);
}

/* create suffix array from fasta file */
void sa_fasta_to_file(const char *fasta_file, const char *suffix_file)
{
	sauchar_t *T;
	saidx64_t *SA, N;
//	uint8_t   *LCP;

	fasta_t   *ffp;
	char      **seq, **name;
	saidx64_t *L;
	saidx64_t *sum_l;
	saidx64_t n_seqs = 1024;

	saidx64_t i, j, K;

	ffp = fasta_open(fasta_file);
	if (ffp == NULL) {
		fprintf(stderr, "Failed opening file %s\n", fasta_file);
		exit(EXIT_FAILURE);
	}
	debug("fasta\n");

	L    = (saidx64_t*)malloc(n_seqs*sizeof(saidx64_t));
	seq  = (char**)malloc(n_seqs*sizeof(char*));
	name = (char**)malloc(n_seqs*sizeof(char*));
	if ((L == NULL) || (seq == NULL) || (name == NULL)) {
		fprintf(stderr, "Failed allocating suffix array a\n");
		exit(EXIT_FAILURE);
	}

	for(i = 0, N = 0; fasta_read(ffp, &seq[i], &name[i], &L[i]); i++) {
		debug("alloc \"%s\" (%llu) i: %d\n", name[i], L[i], i);

		// concatenate a reversed complement of the sequence
		seq[i+1] = strdup(seq[i]);
		name[i+1] = strdup(name[i]);
		L[i+1] = L[i];

		rna_reverse(seq[i+1], L[i]);
		rna_comp(seq[i+1]);

		N += (saidx64_t)L[i];
		N += (saidx64_t)L[i];

		i += 1;

		if (i+2 >= n_seqs) {
			n_seqs += (1<<15);
			L    = (saidx64_t*)realloc(L, n_seqs*sizeof(saidx64_t));
			seq  = (char**)realloc(seq, n_seqs*sizeof(char*));
			name = (char**)realloc(name, n_seqs*sizeof(char*));
			//fprintf(stderr, "reallocating... n_seqs: %lu\n", n_seqs);
			if ((L == NULL) || (seq == NULL) || (name == NULL)) {
				fprintf(stderr, "Failed allocating suffix array b\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	K = i;

	//fprintf(stderr, "n_seqs: %d\n", n_seqs);
	fasta_close(ffp);

	debug("concatenate (%li)\n", N);
	T = (sauchar_t*)malloc((N+1)*sizeof(sauchar_t));
	if(!T) {
		fprintf(stderr, "failed to allocate memory for T\n");
		exit(1);
	}
	for(i = 0, j = 0; i < N; i += (L[j]), j++) {
		strtolower(seq[j]); 
		strcpy((char*)&T[i], (char*)seq[j]); 
		free(seq[j]);
	}
	free(seq);

	//fprintf(stderr, "N: %d\n", N);
	debug("suffix array\n");
	SA = (saidx64_t*)malloc(N*sizeof(saidx64_t));
	if(!SA) {
		fprintf(stderr, "failed to allocate memory for SA\n");
		exit(1);
	}
	if(divsufsort64(T, SA, N) != 0) {
		fprintf(stderr, "divsufsort64 failed\n");
		exit(1);
	}

/*
	for (i = 0; i < N; i++) {
		fprintf(stderr, "%ld %8s\n", i,  &T[SA[i]]);

	}

	for (i = 0; i < K; i++) {
		fprintf(stderr, "L[%ld]= %ld\n", i, L[i]);
	}
*/
	debug("merge\n");
	sum_l = (saidx64_t*)malloc((K+1)*sizeof(saidx64_t));
	if (sum_l == NULL) {
		fprintf(stderr, "Failed allocating sum_l\n");
		exit(EXIT_FAILURE);
	}
	sum_l[0] = 0;
	for(i = 0; i < K; i++) {
		sum_l[i+1] = sum_l[i] + L[i];
	}
	sa_merge(T, SA, N, sum_l, K);

	free(T);
	debug("write\n");
	sa_write(suffix_file, SA, N, L, K, name);
	free(L);

	debug("free\n");
	for(i = 0; i < K; i++) 
		free(name[i]);
	free(name);
	free(SA);
}

