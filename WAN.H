// header file: Wan.h //
#include "alloc.h"
#include "stdio.h"
#include "conio.h"
#include "stdlib.h"
#include "string.h"
#include "io.h"
#include "math.h"
#include "fcntl.h"
#include "dir.h"
#include "alloc.h"
#include "dos.h"
#include "iostream.h"
#include "graphics.h"
#include "complex.h"

#define MAXNM 512
#define N 10
#define PI 3.141592653589793
#define C0 0.51210

 char filename[20], label[20], *string=NULL;
 int num, numb, deriv, ntag[4]={0,0,0,0}, flag, kweight;

 float b[N],c[N],d[2];
 float ev,es, ef, step;
 float X[MAXNM],Y[MAXNM];     //* working space *//
 float EX0[MAXNM],EY0[MAXNM], bkgd[MAXNM], EY[MAXNM];  //* raw absorption, background, bk-removed*//
 float EX2[MAXNM],EY2[MAXNM];  //* normalized absoprtion *//
 float KX[MAXNM],KY[MAXNM], KYW[MAXNM]; //* EXAFS and weighted *//
 float RX[MAXNM],RY[MAXNM], RYr[4*MAXNM], RYi[4*MAXNM]; //* Fourier Transform ampl. real imag *//
 float BX[MAXNM],BY[MAXNM], Amp[MAXNM], Pha[MAXNM]; //* Back Transform real,ampl,phase *//

// EXIT.cpp //
void exitwindow()
 { char c;
   c=getch();
   while(c!=27) { c=getch();}
   setfillstyle(1,0);
   bar(10,41,624,475);
  }
//* the end of EXIT.cpp *//

// * End of wan.h  *//
