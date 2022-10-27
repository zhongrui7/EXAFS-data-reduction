#include "stdlib.h"
#include "fcntl.h"
#include <io.h>
#include "dir.h"
#include "alloc.h"
#include "dos.h"
#include "string.h"
#include "iostream.h"
#include "stdio.h"
#include "conio.h"

#define pi 3.141592653589793

struct FIT
{float k;
 float phc;
 float mag;
 float phase;
 float factor;
 float lamda;
 float p;
 };
 struct FIT fitdat[4][420];

/*//End of FIT.h//*/