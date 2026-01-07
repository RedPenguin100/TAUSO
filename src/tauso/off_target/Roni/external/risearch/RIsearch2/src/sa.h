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
#ifndef __SA_H__
#define __SA_H__ 1

#include <stdint.h>
#include <divsufsort64.h>

/* SA is stored in bits 34-0 */
#define XSA(_x) ((_x)&0x00000003ffffffff)
/* STR is stored in bits 38-34 */
#define XRIS(_x) (((_x)&0x0000001c00000000)>>34)
#define XSTR(_x) (sauchar_t)("\x00" "agcun  " "\x00" "AGCUN"[(((_x)&0x0000003c00000000)>>34)])
#define XSTRM(_x) (sauchar_t)("\x00" "agcun"[(((_x)&0x0000001c00000000)>>34)])
/* LCP is stored in bits 64-56 */
//#define XLCP(_x) (uint8_t)(((_x)&0xff00000000000000)>>56)
#define XLCP(_x) (uint8_t)(((_x)&0x0000000000000000)>>56)
/* IDX is stored in bits 56-38 */
// IDX is now stored in bits 64-38
//#define XIDX(_x) (uint32_t)(((_x)&0x00ffffc000000000)>>38)
#define XIDX(_x) (uint32_t)(((_x)&0xffffffc000000000)>>38)

#define MKSA(_x) ((_x)&0x00000003ffffffff)
#define MKLCP(_x) (((saidx64_t)((_x)&0x00))<<56)
#define MKSTR(_x) (sa_strmap(_x))
#define MKIDX(_x) ((((saidx64_t)(_x))<<38)&0xffffffc000000000);

extern void sa_read(const char *filename, saidx64_t **sa, saidx64_t *n, saidx64_t **l, saidx64_t **sum_l, saidx64_t *k, char ***name);
extern void sa_write(const char *filename, saidx64_t *sa, saidx64_t n, saidx64_t *l, saidx64_t k, char **name);
extern void sa_debug(const saidx64_t *sa, saidx64_t n);
extern void sa_fasta_to_file(const char *fasta_file, const char *suffix_file);
extern void sa_create_partial_reverse(const char *seq, saidx64_t start, saidx64_t end, saidx64_t cnt, saidx64_t **sa, saidx64_t *n);
extern void sa_create_partial_revcomp(const char *seq, saidx64_t start, saidx64_t end, saidx64_t cnt, saidx64_t **sa, saidx64_t *n);
extern void sa_length_to_index(saidx64_t **l, saidx64_t k);
#endif

