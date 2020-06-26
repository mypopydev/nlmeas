/**
 *
 * Copyright (c) 2020, Jose-Luis Lisani, joseluis.lisani@uib.es
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef IO_RGB_HEADER
#define IO_RGB_HEADER

void input2RGB(unsigned char *input, 
               unsigned char **RR, unsigned char **GG, unsigned char **BB, 
               int size);


void RGB2output(unsigned char *R, unsigned char *G, unsigned char *B, 
                unsigned char *output, int size);

void RGBtoI(unsigned char *R, unsigned char *G, unsigned char *B,
            unsigned char **II, int size);


#endif
