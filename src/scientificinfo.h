#ifndef SCIENTIFICINFO_H
#define SCIENTIFICINFO_H

#include <stdio.h>
#include <stdlib.h>

/* SIGNALS */
typedef int ssignal;

#define SIGSCIENTIFICSTOP 1
#define SIGSCIENTIFICRUN 0

void ScientificVersion();
const char *GetScientificVersion();

#endif
