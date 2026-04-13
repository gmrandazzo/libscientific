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
#include <stdint.h>
#include <unistd.h>
#include <sys/wait.h>
#include "memwrapper.h"

void testmemwrapper_xmalloc_success()
{
    puts("Testing xmalloc success path...");
    size_t size = 100;
    void *ptr = xmalloc(size);
    if (ptr == NULL) {
        puts("xmalloc returned NULL!");
        abort();
    }
    memset(ptr, 0xAA, size);
    xfree(ptr);
    puts("xmalloc success path passed.");
}

void testmemwrapper_xcalloc_success()
{
    puts("Testing xcalloc success path...");
    
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
    puts("xcalloc success path passed.");
}

void testmemwrapper_xrealloc_success()
{
    puts("Testing xrealloc success path...");
    size_t initial_size = 50;
    size_t new_size = 150;
    
    void *ptr = xmalloc(initial_size);
    memset(ptr, 0xBB, initial_size);
    
    void *new_ptr = xrealloc(ptr, new_size);
    if (new_ptr == NULL) {
        puts("xrealloc returned NULL!");
        abort();
    }
    
    /* Verify original data is preserved (first 50 bytes) */
    unsigned char *cptr = (unsigned char *)new_ptr;
    for (size_t i = 0; i < initial_size; i++) {
        if (cptr[i] != 0xBB) {
            printf("xrealloc data corruption at index %zu!\n", i);
            abort();
        }
    }
    
    xfree(new_ptr);
    puts("xrealloc success path passed.");
}

void testmemwrapper_xfree()
{
    puts("Testing xfree...");
    void *ptr = xmalloc(10);
    xfree(ptr);
    
    /* xfree(NULL) should be safe according to standard free() behavior */
    xfree(NULL);
    puts("xfree tests passed.");
}

void testmemwrapper_getnprocessor()
{
    puts("Testing GetNProcessor...");
    size_t n_online = 0;
    size_t n_max = 0;
    
    GetNProcessor(&n_online, &n_max);
    
    printf("Processors: Online=%zu, Max=%zu\n", n_online, n_max);
    
    if (n_online == 0) {
        puts("Warning: GetNProcessor returned 0 online processors.");
    }
    
    puts("GetNProcessor test passed.");
}

void testmemwrapper_malloc_fail()
{
    puts("Testing xmalloc failure...");
    void *ptr = xmalloc((size_t)-1);
    (void)ptr;
}

void testmemwrapper_calloc_fail()
{
    puts("Testing xcalloc failure...");
    void *ptr = xcalloc((size_t)-1, (size_t)-1);
    (void)ptr;
}

void testmemwrapper_realloc_fail()
{
    puts("Testing xrealloc failure...");
    void *ptr = xmalloc(10);
    void *new_ptr = xrealloc(ptr, (size_t)-1);
    (void)new_ptr;
    xfree(ptr);
}

int run_fail_test(void (*func)())
{
    pid_t pid = fork();
    if (pid == 0) {
        /* Child: run the failing function */
        func();
        exit(0); /* Should not reach here if func aborts */
    } else if (pid > 0) {
        /* Parent: check if child aborted */
        int status;
        waitpid(pid, &status, 0);
        if (WIFSIGNALED(status) && (WTERMSIG(status) == SIGABRT)) {
            return 0; /* Passed: child aborted as expected */
        }
        fprintf(stderr, "Child did not abort as expected. Status: %d\n", status);
        return 1; /* Failed: child did not abort */
    }
    perror("fork");
    return 1; /* Failed to fork */
}

int main(int argc, char **argv)
{
    if (argc > 1) {
        if (strcmp(argv[1], "--malloc-fail") == 0) {
            return run_fail_test(testmemwrapper_malloc_fail);
        } else if (strcmp(argv[1], "--calloc-fail") == 0) {
            return run_fail_test(testmemwrapper_calloc_fail);
        } else if (strcmp(argv[1], "--realloc-fail") == 0) {
            return run_fail_test(testmemwrapper_realloc_fail);
        }
    }

    testmemwrapper_xmalloc_success();
    testmemwrapper_xcalloc_success();
    testmemwrapper_xrealloc_success();
    testmemwrapper_xfree();
    testmemwrapper_getnprocessor();
    
    return 0;
}
