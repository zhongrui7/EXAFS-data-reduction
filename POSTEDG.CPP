// postedg.cpp //

#include "spline.cpp"


float Power(float x,int n)
   {
    float p=1.0;
    int i;
    for(i=1; i<=n ;i++)  p = p*x;
    return(p);
    }


void PolSpl(float xx[512],float yy[512],int npts,float xl[9],float xh[9],int nr,float co[36],float nc[9])
  {
   int i,j,ibl,k,nk,ik,m,n,n1,ns,ne,ncol,i1,ni,ni1,nm1,ns1;
   float df[25],a[36][35], xk[9];
   float t,temp;
   int p,nbs[9];

   for(i=0; i< 36; i++){
      co[i] = 0.0;
      for(j=0; j < 35; j++)  a[i][j] = 0.0;
   }

   n = 0;
   nbs[0] = 0;

   for(i=0; i < nr; i++){
     n = n + nc[i];
     nbs[i+1] = n ;
     if(xl[i] >= xh[i]){
	t = xl[i];
	xl[i] = xh[i];
	xh[i] = t;
     }
   }


   n = n + 2 * (nr - 1);
   n1 = n + 1;
   xl[nr] = 0.0;
   xh[nr] = 0.0;

   for(ibl=0; ibl < nr ;ibl++)
     {xk[ibl] = 0.5 * (xh[ibl] + xl[ibl+1]);
      if(xl[ibl] > xl[ibl+1])
	 xk[ibl] = 0.5 * (xl[ibl] + xh[ibl+1]);
      ns = nbs[ibl];
      ne = nbs[ibl+1];
      for(i=0; i < npts ; i++)
	{ if((xx[i] >= xl[ibl]) && (xx[i] <= xh[ibl]))
          { df[ns] = 1.0;
	    ns1 = ns + 1;
	    for(j=ns1; j<ne ;j++) df[j] = df[j-1] * xx[i];
            
            for(j=ns; j<ne ;j++)
	      { for(k=j; k<ne ;k++)
		  a[k][j] = a[k][j] + df[j] * df[k] ;
               a[n1-1][j] = a[n1-1][j] + df[j] * yy[i];
               }
          }
       }
    }

    ncol = nbs[nr] - 1;
    nk = nr - 1;

    if(nk != 0 )
      {for(ik=0; ik<nk ;ik++)
	{ ncol ++;
	  ns = nbs[ik];
	  ne = nbs[ik+1] ;
	  a[ncol][ns] = -1.0;
	  ns ++;
	  for(i=ns; i<ne ; i++)
	     a[ncol][i] = a[ncol][i-1] * xk[ik];
          ncol ++;
	  a[ncol][ns] = -1.0;
	  ns++;
	  if(ns <= ne)
	   { for(i=ns; i<ne ;i++)
             {  p = i - ns + 1;
		a[ncol][i] = (ns - i -2) * Power((xk[ik]),p);
	     }
            }
	  ncol --;
	  ns = nbs[ik+1];
	  ne = nbs[ik+2] ;
	  a[ncol][ns] = 1.0;
	  ns++;
	  for(i=ns; i<ne ;i++)
	     a[ncol][i] = a[ncol][i-1] * xk[ik];
          ncol ++;
	  a[ncol][ns] = 1.0;
	  ns++;
	  if(ns <= ne)
	    { for(i=ns; i<ne ;i++)
		a[ncol][i] = (i -ns + 2) * Power((xk[ik]),(i - ns + 1));
             }
       }
    }

    for(i=0; i<n ;i++){
       i1 = i ;
       for(j=0; j<i1; j++) 
	   {a[j][i] = a[i][j];} 
    }

    nm1 = n -1;

  for(i=0; i<nm1; i++){
     i1 = i + 1;
     m = i;
     t = fabs(a[i][i]);
     for(j=i1; j<n ;j++){
        if(t < fabs(a[i][j])){
           m = j;
           t = fabs(a[i][j]);
        }
     }
     if(m != i){
        for(j=i; j<n1 ;j++){
	   t = a[j][i];
	   a[j][i] = a[j][m];
	   a[j][m] = t; 
        }
     }
     for(j=i1; j<n ;j++){
         t = a[i][j] / a[i][i];
         for(k=i1; k<n1 ;k++)
	     a[k][j] = a[k][j] - t * a[k][i];
     }
  }

       co[n-1] = a[n1-1][n-1] / a[n-1][n-1];

       for(i=0; i<=nm1 ;i++)
	{ ni = n - 1 - i;
	  t = a[n1-1][ni];
	  ni1 = ni + 1;
	  for(j=ni1; j<n ;j++){
	     t = t - co[j] * a[j][ni];
	     }
          co[ni] = t / a[ni][ni];
       }
 }


void normaliz(float U0[],float *step )
   {float a[5][5];
    float x0,x1,x2,x3,x4,x5,x6,u,u1,u2,u3,u4,t,postu0;
    int i,j;
    for(i=0;i<5;i++)
     {c[i]=0.0;
     for(j=0;j<5;j++) {a[i][j]=0.0;}
     }
    i=0;
    while((X[i]-ef)<100.0) {i=i+1;}
    j=0;
   while((X[i]-ef)<=300.0)
    {t=X[i];
     u=U0[i];
     x0=t*t;
     x1=x0*t;
     x2=x1*t;
     x3=x2*t;
     x4=x3*t;
     x5=x4*t;
     x6=x5*t;
     u1=u*t;
     u2=u*x0;
     u3=u*x1;
     u4=u*x2;
     i=i+1;
     j=j+1;
     a[0][0]=j;
     a[0][1]+=t;
     a[0][2]+=x0;
     a[0][3]+=x1;
     a[0][4]+=x2;
     a[1][0]+=t;
     a[1][1]+=x0;
     a[1][2]+=x1;
     a[1][3]+=x2;
     a[1][4]+=x3;
     a[2][0]+=x0;
     a[2][1]+=x1;
     a[2][2]+=x2;
     a[2][3]+=x3;
     a[2][4]+=x4;
     a[3][0]+=x1;
     a[3][1]+=x2;
     a[3][2]+=x3;
     a[3][3]+=x4;
     a[3][4]+=x5;
     a[4][0]+=x2;
     a[4][1]+=x3;
     a[4][2]+=x4;
     a[4][3]+=x5;
     a[4][4]+=x6;
     c[0]+=u;
     c[1]+=u1;
     c[2]+=u2;
     c[3]+=u3;
     c[4]+=u4;
    }
    gaussj2(a,5);
    postu0=c[0]+c[1]*ef+c[2]*ef*ef+c[3]*ef*ef*ef+c[4]*ef*ef*ef*ef;
    *step=postu0;  //-preu0;
 }


void mu0(int wmax)
  {  float *U0, *X0, xl[N], xh[N], nc[N], co[4*N];
     float  delta, xval, yval, cstart;
     int i,j,k,menu,nr=3,d=3,l;
 	char n;
     X0=new float[wmax];
     U0=new float[wmax];

     for(i=0; i<numb; i++){X[i]=EX2[i]; Y[i]=EY2[i];}
     for(i=numb; i<MAXNM; i++){X[i]=0.0; Y[i]=0.0;}

   sprintf(label,"Energy [eV]");
   plot(X, Y, 0, 0, numb);

outtextxy(80, 450,"Select method:");
outtextxy(250, 450,"0)Equal_interval; 1)Manual");
     n=getch();
     while(n!=13)
      {  setcolor(RED);
	while(n<48||n>49)
	{
	 sprintf(string,"%c",n);
	 outtextxy(400,450,string);
	 n=getch();
	 }

	 sprintf(string,"%c",n);
	 outtextxy(400,450,string);
       if(n>47 && n<50) menu=n-48;
       n=getch();
      }
	setcolor(15);	
     setfillstyle(1,0);
     bar(10,430,620,470);
switch(menu)
  {
   case 0:
    {outtextxy(150, 450,"Spline intervals [2-6]:");
     n=getch();
     while(n!=13)
      {  
	while(n<50||n>55)
	{
	 sprintf(string,"%c",n);
	 outtextxy(300,450,string);
	 n=getch();
	 }
	 sprintf(string,"%c",n);
	 outtextxy(300,450,string);
       if(n>49 && n<56) nr=n-48;
       n=getch();
      }
  	xl[0] = X[0];
   	xh[nr-1] = X[numb-2];
   	delta = (xh[nr-1] - xl[0])/float(nr);
     setwritemode(XOR_PUT);
     for(k=0; k < nr-1; k++)
 	{
	xl[k+1] = xl[k] + delta;
	xh[k] = xl[k+1];
       j=80+int(unit_x*(xh[k]-min_x));
       line(j,400,j,80);
       	}

     outtextxy(150, 465,"Spline degree [2-4]:");
     n=getch();
     while(n!=13)
      {  
	while(n<50||n>52)
	{
	 sprintf(string,"%c",n);
	 outtextxy(300,465,string);
	 n=getch();
	 }

	 sprintf(string,"%c",n);
	 outtextxy(300,465,string);
       if(n>49 && n<53) d=n-48;
       n=getch();
      }
	if(nr==6&&d==4) d=3;
   	for(j=0;j<nr; j++) nc[j]= d+1; 
     }	
	break;

   case 1:
     {outtextxy(200,460,"'Space' spline: A<-S|D->F");
      while(n!=32) {n=getch();}
       setcolor(15);
       k=numb/3; j=1;
      setwritemode(XOR_PUT);
      while(n!=27&&j<7)
      { 
	sprintf(string,"%-4.2f",EX2[k]);
	moveto(300,70);
	outtext("E=");
	outtextxy(350,70,string);
	delta=16;
     while(n!=13)
	 {
       if(n=='A'||n=='a')	 { k=k-delta;  }
       if(n=='F'||n=='f')	 { k=k+delta;  }
       if(n=='S'||n=='s')	 { k=k-1;  }
       if(n=='D'||n=='d')	 { k=k+1;  }
       if(k<9) k=9;
       if(k>numb-10) k=numb-10;
       i=80+int(unit_x*(EX2[k]-min_x));
       line(i,400,i,80);
	sprintf(string,"%-8.4f",EX2[k]);
	moveto(300,70);
	outtext("E=");
	outtextxy(450,70,string);
	n=getch();
	line(i,400,i,80);
	setfillstyle(1,1);
	bar(280,60,420,80);
       }
       xl[j]=EX2[k];
       j=j+1; nr=j;
       n=getch();
	line(i,400,i,80);
      }

	xl[0] = X[0]; xh[nr-1]=X[numb-2];
        for(i=0; i < nr;i++){
	  xh[i] = xl[i+1];
	  nc[i]=4;
          }
   }
	break;
 }

	PolSpl(X,Y,numb,xl,xh,nr,co,nc);

  /** change xl[1] and xh[nr] to extrapolate the fit.   */
   xl[0] = xh[nr-1] = X[0];
   for(i=0; i < numb; i++){
       xl[0]  = (xl[0] >=  X[i]? X[i]:xl[0]);
       xh[nr-1] = (xh[nr-1] >= X[i]? xh[nr-1]:X[i]);
   }
  

   /** calculate the first points.    */
   xval = X[0];
   yval = 0.0;
   for(k=0; k < nc[1]; k++) yval += co[k] * Power((xval),k);

   X0[0] = xval;
   U0[0] = yval;

   /* now the rest of the points.    */
   for(i=0; i < numb;i++){
      for(j=0; j < nr; j++){
	 if((X[i] > xl[j]) && (X[i] <= xh[j])){
	    cstart = 0;
	    for(l=0; l <= j-1; l++)
	        cstart += nc[l];
	    xval = X[i];
	    yval = 0.0;
	    for(k=0; k < nc[j]; k++)
		yval += co[cstart + k] * Power((xval),k);
            X0[i] = xval;
	    U0[i] = yval;
           }
	 }
	}

 getch();
 plotD(X,Y, X0, U0, 0, 0, numb);
    outtextxy(540,460,"'Esc'exit");
    getch();
 normaliz(U0,&step);
 for(i=0;i<numb;i++)
      {
	Y[i]=(Y[i]-U0[i])/step;
       }

     delete[] X0;
     delete[] U0;
  }



void E2k(int pnum,int wmax)
    { char c1,c,string[4];
      float yp1,ypn, deta;
      float *X1,*Y1, *Y2;
      X1=new float[pnum];
      Y1=new float[pnum];
      int i,j;
      Y2=new float[pnum];

     for(i=0;i<pnum;i++)  Y2[i]=0.0;
     for(i=0;i<pnum;i++)
      {
        X[i]=sqrt(X[i]-ef)*C0;
	X1[i]=X[i];
	Y1[i]=Y[i];
       }

       yp1=0.0;
       ypn=(Y[pnum-2]-Y[pnum-1])/(X[pnum-2]-X[pnum-1]);
  spline(X1,Y1,pnum,yp1,ypn,Y2);

      deta=0.05;
  for(j=0;j<wmax;j++)
   { X[j]=deta*float(j); 
   if(X[j]>=X1[0]&&X[j]<=X1[pnum-1])
    {Y[j]=splint(X1,Y1,Y2,pnum,X[j]);}
   else
    {Y[j]=0.0;}
      }

//* store EXAFS data into KX[],KY[] *//
    for(j=0;j<MAXNM;j++){KX[j]=X[j]; KY[j]=Y[j];}  

 	exitwindow();   
	sprintf(label,"k [A^-1]");
  plot(X, Y, 0.0, 20.0, wmax);
   ntag[2]=wmax;
	delete[] Y2;
	delete[] X1;
	delete[] Y1;
  }
/*///End of POSTEDG.cpp//*/