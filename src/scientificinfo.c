#include "scientificinfo.h"

#define major_ 0
#define minor_ 7
#define patch_ 4

void ScientificVersion()
{
  printf("Scientific Library was writen by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>\nVersion: %d.%d.%d\n", major_, minor_, patch_);
}

const char *GetScientificVersion()
{
  static char c[10];
  sprintf(c, "%d.%d.%d", (int)major_, (int)minor_, (int)patch_);
  return c;
}
