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
#ifndef __DSM_H__
#define __DSM_H__

#include <stdint.h>

/* dsm gap value */
#define GAP   0

/* dsm size */
#define DSM_N 6

/* dsm type */
typedef int dsm_t[DSM_N][DSM_N][DSM_N][DSM_N];

/* dsm RNA-RNA energy tables */
extern const dsm_t dsm_t99_pos;
extern const dsm_t dsm_t99_neg;
extern const dsm_t dsm_t04_pos;
extern const dsm_t dsm_t04_neg;
extern const dsm_t dsm_extend_pos;
extern const dsm_t dsm_extend_neg;


#endif
