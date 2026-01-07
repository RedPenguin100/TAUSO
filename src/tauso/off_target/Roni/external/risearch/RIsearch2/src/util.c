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
#include <stdarg.h>
#include <time.h>

#include "main.h"

int verbose = 0;

void debug(const char *msg, ...)
{
	char buffer[80];
	time_t t;
	va_list args;
	if(verbose) {
		t = time(NULL);
		strftime(buffer, 80, "%H:%M:%S", localtime(&t));
		fprintf(stderr, "RIsearch2 [%s]: ", buffer);
		va_start(args, msg);
		vfprintf(stderr, msg, args);
		va_end(args);
	}
}

