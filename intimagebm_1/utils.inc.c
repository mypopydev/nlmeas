/*
 * Copyright (c) 2013, Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under, at your option, the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version, or
 * the terms of the simplified BSD license.
 *
 * You should have received a copy of these licenses along this
 * program. If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

/**
 * @file utils.inc.c
 * @brief miscellaneous utlities
 *
 * This file should be #included in the source file using these
 * utilities, not linked with it.
 * 
 * The routines here are based on the work of Enric Meinhardt-Llopis.
 *
 * @author Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
 */

#ifndef _UTILS_INC_C
#define _UTILS_INC_C

#include <stdlib.h>
#include <stdio.h>

#if (defined(__STDC__) && defined(__STDC_VERSION__) \
     && (__STDC_VERSION__ >= 199901L))
#define I_CAN_HAS_C99
#endif

#ifdef I_CAN_HAS_C99
#define LOCAL inline static
#else
#define LOCAL static
#endif

/**
 * @brief Fatal error
 *
 * Print a message to standard-error output and exit.
 */
LOCAL void error(const char *msg)
{
   fprintf(stderr, "Error: %s\n", msg);
   exit(EXIT_FAILURE);
}

/**
 * @brief Safe memory allocation
 *
 * Print an error and exit if fails.
 */
LOCAL void *xmalloc(size_t size)
{
   void *p;
   if (0 == size)
      error("xmalloc: zero size");
   p = malloc(size);
   if (NULL == p)
      error("xmalloc: out of memory");
   return p;
}

/**
 * @brief Safe memory reallocation
 *
 * Print an error and exit if fails.
 */
LOCAL void *xrealloc(void *p, size_t size)
{
   if (0 == size)
      error("xrealloc: zero size");
   if (NULL == p)
      error("xrealloc: NULL pointer");
   p = realloc(p, size);
   if (NULL == p)
      error("xrealloc: out of memory");
   return p;
}

/**
 * @brief Safe memory array allocation
 *
 * Print an error and exit if fails.
 */
LOCAL void *xcalloc(size_t nobj, size_t size)
{
   void *p;
   if (0 == size)
      error("malloc: zero size");
   p = calloc(nobj, size);
   if (NULL == p)
      error("malloc: out of memory");
   return p;
}

#undef LOCAL

#ifdef I_CAN_HAS_C99
#undef I_CAN_HAS_C99
#endif

#endif				/* !_UTILS_INC_C */
