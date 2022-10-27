// * amebsa1.cpp  *//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>

#define FABS(a) ((a)>0.0?(a):-(a))
#define SWAP(a,b) {float temp=(a);(a)=(b);(b)=temp;}
#define NR_END 1
#define FREE_ARG char*
#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

 static long idum=100000000;
 float tt;

 /*allocate a float vector with subscript range v[nl...nh]*/
float *vector(long nl,long nh)
  { float *v;
   v=(float *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(float)));
   return v-nl+NR_END;
   }


 /*allocate a float matrix with subscript range m[nrl...nrh][ncl...nch]*/
float **matrix(long nrl,long nrh,long ncl,long nch)
 {long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;
  m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  m+=NR_END;
  m-=nrl;
  m[nrl]=(float *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  m[nrl]+=NR_END;
  m[nrl]-=ncl;
  for(i=nrl+1;i<=nrh;i++)
   m[i]=m[i-1]+ncol;
  return m;
 }

 void free_vector(float *v,long nl,long nh)
 { free((FREE_ARG)(v+nl-NR_END));
  }

 void free_matrix(float **m,long nrl,long nrh,long ncl,long nch)
  {
   free((FREE_ARG)(m[nrl]+ncl-NR_END));
   free((FREE_ARG)(m+nrl-NR_END));
   }
/*//End of NRUTIL.H//*/



float ran1(long *idum)
   {
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;
   if((*idum<=0)||(!iy))
     {
     if(-(*idum)<1) *idum=1;
      else *idum=-(*idum);
      for(j=NTAB+7;j>=0;j--)
       {k=(*idum)/IQ;
        *idum=IA*(*idum-k*IQ)-IR*k;
        if(*idum<0) *idum+=IM;
        if(j<NTAB) iv[j]=*idum;
       }
       iy=iv[0];
     }
     k=(*idum)/IQ;
     *idum=IA*(*idum-k*IQ)-IR*k;
     if(*idum<0)  *idum+=IM;
     j=iy/NDIV;
     iy=iv[j];
     iv[j]=*idum;
     if((temp=(1.0/IM)*iy)>RNMX) return RNMX;
     else
      return temp;
    }

float Camotsa(float **p,float fvertex[],float psum[],int dim,float best_P[],float *best_f,
        float (*func)(float[],float [],float [],float ,struct FIT [][420],int ,int ,int[],int ,float ,float),int ihi,float *yhi,float fac ,float S0 ,struct FIT fitdat[][420],int fitma,int ma,int flag[],int coeff,float ekl,float ekh,float x[],float y[],float **init)
         {int j,actor=1;
          float fac1,fac2,fflu,ftry,*ptry;
          ptry=vector(1,dim);
          fac1=(1.0-fac)/dim;
          fac2=fac1-fac;
         for(j=1;j<=dim;j++)
           {ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
            if(init[j][1]<=ptry[j]&&init[j][2]>=ptry[j])
              actor=1;
            else
              {actor=0;
               break;
               }
              }
         if(actor==1)
            {ftry=(*func)( x,y,ptry,S0,fitdat,fitma,ma,flag,coeff,ekl,ekh);
             if(ftry<*best_f){
                for(j=1;j<=dim;j++) best_P[j]=ptry[j];
                *best_f=ftry;
                     }

         fflu=ftry-tt*log(ran1(&idum));

         if(fflu<*yhi){
            fvertex[ihi]=ftry;
            *yhi=fflu;
            for(j=1;j<=dim;j++)
               { psum[j]+=ptry[j]-p[ihi][j];
                 p[ihi][j]=ptry[j];
                }
                  }
            free_vector(ptry,1,dim);
            return fflu;
              }
          else
            {free_vector(ptry,1,dim);
             return *yhi;
             }
            }



void Camebsa(float **p,float fvertex[],int dim,float best_P[],float *best_f, float ftol,
              float (*func)(float[],float [],float [],float  ,struct FIT [][420],int ,int ,int[],int ,float ,float),int *iter,float temptr,float S0 ,struct FIT fitdat[][420],int fitma,int ma,int flag[],int coeff,float ekl,float ekh,float x[],float y[],float **init)
      {
        int i,ihi,ilo,j,m,n,actor=1,mpts=dim+1;
        float rtol,sum,temp,fhi,flo,fnhi,fsave,ft,ftry,*psum;

        psum=vector(1,dim);
        tt=-temptr;
        for(n=1;n<=dim;n++){
        for(sum=0.0,m=1;m<=mpts;m++) sum+=p[m][n];
        psum[n]=sum;}
        for(;;){
           ilo=1;
           ihi=2;
           fnhi=flo=fvertex[1]+tt*log(ran1(&idum));
           fhi=fvertex[2]+tt*log(ran1(&idum));
          if(flo>fhi){
             ilo=2;
             ihi=1;
             fnhi=fhi;
             fhi=flo;
             flo=fnhi;
              }
          for(i=3;i<=mpts;i++){
              ft=fvertex[i]+tt*log(ran1(&idum));
            if(ft<=flo){
               ilo=i;
               flo=ft;
               }
            if(ft>fhi){
               fnhi=fhi;
               ihi=i;
               fhi=ft;
               }
            else if (ft>fnhi){
               fnhi=ft;
               }
             }
        rtol=2.0*FABS(fhi-flo)/(FABS(fhi)+FABS(flo));
        if(rtol<ftol||*iter<0){
          SWAP(fvertex[1],fvertex[ilo])
          for(n=1;n<=dim;n++) SWAP(p[1][n],p[ilo][n])
          break;
          }
    *iter-=2;
    ftry=Camotsa(p,fvertex,psum,dim,best_P,best_f,func,ihi,&fhi,-1.0, S0 ,fitdat, fitma, ma, flag, coeff, ekl, ekh,x,y,init);
   if(ftry<=flo)
    ftry=Camotsa(p,fvertex,psum,dim,best_P,best_f,func,ihi,&fhi,2.0 ,S0 ,fitdat, fitma, ma, flag, coeff, ekl, ekh,x,y,init);
   else if(ftry>=fnhi){
       fsave=fhi;
       ftry=Camotsa(p,fvertex,psum,dim,best_P,best_f,func,ihi,&fhi,0.5, S0 ,fitdat, fitma, ma, flag, coeff, ekl, ekh,x,y,init);
       if(ftry>=fsave){
          for(i=1;i<=mpts;i++){
              if(i!=ilo){
                 for(j=1;j<=dim;j++){
                    psum[j]=0.5*(p[i][j]+p[ilo][j]);
                    if(init[j][1]>=psum[j]||init[j][2]<=psum[j])
                      {actor=0;
                       break;
                       }
                   p[i][j]=psum[j];
                         }
                if(actor==1)
                   fvertex[i]=(*func)( x,y,psum,S0,fitdat,fitma,ma,flag,coeff,ekl,ekh);
                    }
                      }
            *iter-=dim;
         for(n=1;n<=dim;n++){
            for(sum=0.0,m=1;m<=mpts;m++) sum+=p[m][n];
            psum[n]=sum;}
                  }
                }
             else
                ++(*iter);
          }
         free_vector(psum,1,dim);
     }

// * The END of amebsa1.cpp  *//

// End of amebsa1.cpp //