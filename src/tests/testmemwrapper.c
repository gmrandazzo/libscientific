/* Unit tests for the memory wrapper module.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "memwrapper.h"

void test_xcalloc()
{
    puts("Testing xcalloc...");
    
    size_t nmemb = 10;
    size_t size = sizeof(int);
    int *ptr = xcalloc(nmemb, size);
    
    if (ptr == NULL) {
        puts("xcalloc returned NULL!");
        abort();
    }
    
    /* Verify zero initialization */
    for (size_t i = 0; i < nmemb; i++) {
        if (ptr[i] != 0) {
            printf("xcalloc memory not zeroed at index %zu!\n", i);
            abort();
        }
    }
    
    xfree(ptr);
    
    /* Test zero size cases (should return a valid pointer that can be freed) */
    ptr = xcalloc(0, size);
    if (ptr == NULL) {
        puts("xcalloc(0, size) returned NULL!");
        abort();
    }
    xfree(ptr);
    
    ptr = xcalloc(nmemb, 0);
    if (ptr == NULL) {
        puts("xcalloc(nmemb, 0) returned NULL!");
        abort();
    }
    xfree(ptr);
    
    puts("xcalloc tests passed.");
}

int main()
{
    test_xcalloc();
    return 0;
}
