// SPLINE.cpp //

#define NMAX 50

void gaussj2(float a[][5],int n)
   {
     int i,j,k,l,ll,icol,irow,indxc[NMAX],indxr[NMAX],ipiv[NMAX];
     float big,dum,pivinv;
     for(j=0;j<n;j++) ipiv[j]=0;
     for(i=0;i<n;i++)
       {big=0;
        for(j=0;j<n;j++)
         {if(ipiv[j]!=1)
          {for(k=0;k<n;k++)
            {if(ipiv[k]==0)
              {if(fabs(*(*(a+j)+k))>=big)
                {big=fabs(*(*(a+j)+k));
                 irow=j;
                 icol=k;}
              }
              if(ipiv[k]>1)
             { cout<<"singlar matrix in gauss"<<endl;
              exit(1);}
              }
              }
              }
        ipiv[icol]=ipiv[icol]+1;
        if(irow!=icol)
          {for(l=0;l<n;l++)
            {dum=*(*(a+irow)+l);
             *(*(a+irow)+l)=*(*(a+icol)+l);
             *(*(a+icol)+l)=dum;
             }

             dum=c[irow];
             c[irow]=c[icol];
             c[icol]=dum;
            }
         indxr[i]=irow;
         indxc[i]=icol;
         if (*(*(a+icol)+icol)==0)
            {cout<<"singular matrix in gaussj"<<endl;
             exit(1);
            }
            pivinv=1.0/(*(*(a+icol)+icol));
            *(*(a+icol)+icol)=1;
            for(l=0;l<n;l++)
            {*(*(a+icol)+l)=*(*(a+icol)+l)*pivinv;}

             {c[icol]=c[icol]*pivinv;}
            for(ll=0;ll<n;ll++)
             {if(ll!=icol)
               {dum=*(*(a+ll)+icol);
                *(*(a+ll)+icol)=0;
                for(l=0;l<n;l++)
                *(*(a+ll)+l)=*(*(a+ll)+l)-*(*(a+icol)+l)*dum;
                 c[ll]=c[ll]-c[icol]*dum;
              }}
              }
        for(l=n-1;l>=0;l--)
          {if(indxr[l]!=indxc[l])
            {for(k=0;k<n;k++)
             {dum=*(*(a+k)+indxr[l]);
             *(*(a+k)+indxr[l])= *(*(a+k)+indxc[l]);
              *(*(a+k)+indxc[l])=dum;
             }}
             }
       }

//*Given a tabulated function y of x (unordered), and Given the values of the first derivatives at the end points This routine returns an array y2, that contains the second derivatives of the function at the tabulated points.*//
void spline(float X1[],float Y1[],int n,float yp1,float ypn,float Y2[])
 { float *u, p,qn,sig,un;
   int i,k;
   u=new float[n];
   if(yp1>0.99e30) Y2[0]=u[0]=0.0;
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
   
    if(ypn>0.99e30) qn=un=0.0;
    else  
      {qn=0.5;
       un=(3.0/(X1[n-1]-X1[n-2]))*(ypn-(Y1[n-1]-Y1[n-2])/(X1[n-1]-X1[n-2]));
       } 
    Y2[n-1]=(un-qn*u[n-2])/(qn*Y2[n-2]+1.0);
    for(k=n-2;k>=0;k--)  Y2[k]=Y2[k]*Y2[k+1]+u[k];
      delete[] u;
  }


//*Given the arrays xa(ordered, ya of length n, which tabulate a function and given the array y2a which is the output of spline and an unordered array xq, this routine returns a cubic-spline interpolated array yq.*//
  float splint(float Xa[],float Ya[],float Y2[],int n,float x)
    {float y,h,b,a;
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
    printf("\nbad Xa input to routine splint");
   a=(Xa[khi]-x)/h;
   b=(x-Xa[klo])/h;
   y=a*Ya[klo]+b*Ya[khi]+((a*a*a-a)*Y2[klo]+(b*b*b-b)*Y2[khi])*(h*h)/6.0;

   return(y);
   }

/*////End of SPLINE.cpp/////*/