/* Provides memory management wrappers and safeguards.
 * Copyright (C) 2016-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include "memwrapper.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#endif


void *xmalloc(size_t size)
{
	if (size == 0) size = 1;
	void *ptr = malloc(size);
	if (!ptr) {
		fprintf(stderr, "[Libscientific] Memory exhausted allocating %zu bytes\n", size);
		abort();
	}
	return ptr;
}

void *xrealloc(void *ptr, size_t size)
{
	if (size == 0) { free(ptr); return NULL; }
	void *new_ptr = realloc(ptr, size);
	if (!new_ptr) {
		fprintf(stderr, "[Libscientific] Memory exhausted reallocating %zu bytes\n", size);
		abort();
	}
	return new_ptr;
}

void xfree(void *ptr)
{
	free(ptr);
}

void GetNProcessor(size_t *nprocs_online, size_t *nprocs_max)
{
	if (nprocs_online) *nprocs_online = 1;
	if (nprocs_max) *nprocs_max = 1;

#if defined(_WIN32)
	SYSTEM_INFO info;
	GetSystemInfo(&info);
	size_t n = (size_t)(info.dwNumberOfProcessors > 0 ? info.dwNumberOfProcessors : 1);
	if (nprocs_online) *nprocs_online = n;
	if (nprocs_max) *nprocs_max = n;
#else
#ifdef _SC_NPROCESSORS_ONLN
	long onln = sysconf(_SC_NPROCESSORS_ONLN);
	if (onln > 0 && nprocs_online) *nprocs_online = (size_t)onln;
#endif
#ifdef _SC_NPROCESSORS_CONF
	long conf = sysconf(_SC_NPROCESSORS_CONF);
	if (conf > 0 && nprocs_max) *nprocs_max = (size_t)conf;
#endif
#endif
}