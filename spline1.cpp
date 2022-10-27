
//* Spline1.cpp *//

void spline(float X1[],float Y1[],int n,float yp1,float ypn,float Y2[])
 {float *u;
  float p,qn,sig,un;
  int i,k;

   u=new float[n];
  if(yp1>0.99e30)
   Y2[0]=u[0]=0.0;
  else
   {Y2[0]=-0.5;
    u[0]=(3.0/(X1[1]-X1[0]))*((Y1[1]-Y1[0])/(X1[1]-X1[0])-yp1);
   }
  for(i=1;i<n-1;i++)
   { 
    sig=(X1[i]-X1[i-1])/(X1[i+1]-X1[i-1]);
    p=sig*Y2[i-1]+2.0;
    Y2[i]=(sig-1.0)/p;
    u[i]=(Y1[i+1]-Y1[i])/(X1[i+1]-X1[i])-(Y1[i]-Y1[i-1])/(X1[i]-X1[i-1]);
    u[i]=(6.0*u[i]/(X1[i+1]-X1[i-1])-sig*u[i-1])/p;
    }
   
  if(ypn>0.99e30)  qn=un=0.0;
  else  {
      qn=0.5;
      un=(3.0/(X1[n-1]-X1[n-2]))*(ypn-(Y1[n-1]-Y1[n-2])/(X1[n-1]-X1[n-2]));
      }
  Y2[n-1]=(un-qn*u[n-2])/(qn*Y2[n-2]+1.0);
  for(k=n-2;k>=0;k--)  Y2[k]=Y2[k]*Y2[k+1]+u[k];
  delete[] u;
  }


 float splint(float Xa[],float Ya[],float Y2[],int n,float x)
  {
    float y,h,b,a;
    int klo,khi,k;
    klo=0;
    khi=n-1;
   while(khi-klo>1)
    {k=(khi+klo)>>1;
     if(Xa[k]>x) khi=k;
     else
      klo=k;
     }
   h=Xa[khi]-Xa[klo];
   if(h==0.0)
     return(Ya[klo]);
   a=(Xa[khi]-x)/h;
   b=(x-Xa[klo])/h;
   y=a*Ya[klo]+b*Ya[khi]+((a*a*a-a)*Y2[klo]+(b*b*b-b)*Y2[khi])*(h*h)/6.0;

   return(y);
   }

// * End of spline1.cpp  *//
