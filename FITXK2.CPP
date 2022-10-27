// * Fitxk2.cpp * //
// * this main program of FITXK2, which inludes files xfit.h, fitxk2.cpp, amebsa1.cpp, spline1.cpp, IO.cpp, and plot.cpp *//

#include "math.h"
#include "graphics.h"

#include "xfit.h"
#include "spline1.cpp"
#include "amebsa1.cpp"
#include "Plot.cpp"
#include "IO.cpp"

float XkGender(float x[],float Y[],float parm[],float S0 ,struct FIT temp[][420],int fitma,int ma,int flag[],int coeff,float ekl,float ekh)
 /*Q(k)=2KR+2@(k)+A(k)+0.263E0*R/k"*/
  {int i,j,m;
   float z=0.0, y=0.0;
   float Ak,Fk,Qk,Mk,Nk,Sk;
   j=1;
   while(x[j]<ekl) {j=j+1;}
   while(x[j]<=ekh)
      { z=0.0;
	Nk=1.0;
       for(m=0;m<coeff;m++) Nk=Nk*x[j];
       for(i=0;i<fitma;i++)
	 { switch(flag[i]){
	   case -1:
        	{//x[j]=sqrt(x[j]*x[j]-0.263*parm[i*ma+4]);
		  Qk= 2*x[j]*parm[i*ma+1]+temp[i][j-1].phc+temp[i][j-1].phase+0.263*parm[i*ma+4]*parm[i*ma+1]/x[j];
		  Fk=(temp[i][j-1].mag*temp[i][j-1].factor*S0/x[j])*exp(-2*parm[i*ma+3]*parm[i*ma+3]*x[j]*x[j]);
		  Mk=Fk*exp(-2*parm[i*ma+1]/temp[i][j-1].lamda)/(parm[i*ma+1]*parm[i*ma+1]);
		  if(x[j]!=0.0)
		     z+=Nk*parm[i*ma+2]*Mk*sin(Qk);
		  else
		     z+=0.0;
		 }
	       break;
           case 2:
	    {// x[j]=sqrt(x[j]*x[j]-0.263*parm[i*ma+4]);
             Sk=4*x[j]*x[j]*parm[i*ma+5]*parm[i*ma+5];
	    Qk= 2*x[j]*parm[i*ma+1]+temp[i][j-1].phc+temp[i][j-1].phase+0.263*parm[i*ma+4]*parm[i*ma+1]/x[j]+atan(2*x[j]*parm[i*ma+5]*(3.0-Sk)/(1.0-3*Sk));
	    Fk=(temp[i][j-1].mag*temp[i][j-1].factor*S0/x[j])*exp(-2*parm[i*ma+3]*parm[i*ma+3]*x[j]*x[j]);
	    Mk=Fk*exp(-2*parm[i*ma+1]/temp[i][j-1].lamda)/(parm[i*ma+1]*parm[i*ma+1]);
	    Ak=Mk/((1+Sk)*sqrt(1+Sk));
	    if(x[j]!=0.0)
	      z+=Nk*parm[i*ma+2]*Ak*sin(Qk);
            else
            z+=0.0;
	    }
	    break;
           case 1:
	    {
	     Sk=4*x[j]*x[j]*parm[i*ma+5]*parm[i*ma+5];
	    Qk= 2*x[j]*parm[i*ma+1]+temp[i][j-1].phc+temp[i][j-1].phase+0.263*parm[i*ma+4]*parm[i*ma+1]/x[j]+atan(2*x[j]*parm[i*ma+5]*(2.0-Sk)/(1.0-2*Sk));
	    Fk=(temp[i][j-1].mag*temp[i][j-1].factor*S0/x[j])*exp(-2*parm[i*ma+3]*parm[i*ma+3]*x[j]*x[j]);
	    Mk=Fk*exp(-2*parm[i*ma+1]/temp[i][j-1].lamda)/(parm[i*ma+1]*parm[i*ma+1]);
	    Ak=Mk/(1.0+Sk);
            if(x[j]!=0.0)
	      z+=Nk*parm[i*ma+2]*Ak*sin(Qk);
	    else
	      z+=0.0;
	    }
	    break;

           case 0:
            {
	     Sk=4*x[j]*x[j]*parm[i*ma+5]*parm[i*ma+5];
	    Qk= 2*x[j]*parm[i*ma+1]+temp[i][j-1].phc+temp[i][j-1].phase+0.263*parm[i*ma+4]*parm[i*ma+1]/x[j]+atan(2*x[j]*parm[i*ma+5]);
            Fk=(temp[i][j-1].mag*temp[i][j-1].factor*S0/x[j])*exp(-2*parm[i*ma+3]*parm[i*ma+3]*x[j]*x[j]);
	    Mk=Fk*exp(-2*parm[i*ma+1]/temp[i][j-1].lamda)/(parm[i*ma+1]*parm[i*ma+1]);
            Ak=Mk/(sqrt(1.0+Sk));
	    if(x[j]!=0.0)
	      z+=Nk*parm[i*ma+2]*Ak*sin(Qk);
	    else
	      z+=0.0;
	    }
	    break;
	   }
     }
      y+=(Y[j]-z)*(Y[j]-z);
     j=j+1;
    }
    return(y);
     }


 float **p, **init, *temp1, *fvertex, *tem, *best;
 float best_f, ftol,temptr,min,max,T,atm,temp;
 float S0,ekl,ekh;
 int dim ,jk,iter,MAX, flag[10];
 int i,j,num,ma,fitma,flag1,handle,coeff,nl,nh;
 FILE *file;
 char c,ch,*ptr,*path;
 char string[40],buff[80],ikk[40],stri[20], feff[10][40];

void main()
{
/*
 float **p, **init, *temp1, *fvertex, *tem, *best;
 float best_f, ftol,temptr,min,max,T,atm,temp;
 float S0,ekl,ekh;
 int dim ,jk,iter,MAX, flag[10];
 int i,j,num,ma,fitma,flag1,handle,coeff,nl,nh;

 FILE *file;
 char c,ch,*ptr,*path;
 char string[40],buff[80],ikk[40],stri[20], feff[10][40];
*/
 float *x,*y,*z;
  x=vector(1,420);
  y=vector(1,420);
  z=vector(1,420);

  ReadInit();
    for(j=1;j<=dim;j++)
      temp1[j]=(init[j][2]-init[j][1]);
    for(i=1;i<=dim+1;i++)
      for(j=1;j<=dim;j++)
      p[i][j]=init[j][1]+ran1(&idum)*temp1[j];
     ReadSource(feff,fitdat,fitma);
     num=ReadOrigin(ikk,x,y);
     printf("Files Checking Passed !\n");getch();
    select_es(y,num,&ekl,&ekh);
     i=1;
   while(x[i]<ekl)    {i=i+1;}
    nl=i;
   while(x[i]<=ekh)   {i=i+1;}
    nh=i;
   for(i=1;i<=dim+1;i++)
    { for(j=1;j<=dim;j++)
       tem[j]=p[i][j];
    fvertex[i]=XkGender(x,y,tem,S0,fitdat,fitma,ma,flag,coeff,ekl,ekh);
    }
    min=fvertex[1];
    max=fvertex[1];
    jk=1;
    for(i=2;i<=dim+1;i++)
     {if(fvertex[i]<min)
       { min=fvertex[i];
	 jk=i;
       }
      if(fvertex[i]>max)  {max=fvertex[i];}
     }
    best_f=1.0e10;
    for(i=1;i<=dim;i++)
     *(best+i)=p[jk][i];
    ftol=EPS;
    T=max;
    atm=0.01;
    MAX=500;
    temptr=T/(1.0-atm);
    flag1=1;
    printf("\n Do you want to fit[Y/N]:");
    ch=getch();
    if(ch=='Y'||ch=='y')
    {
     do{temptr=temptr*(1.0-atm);
	iter=MAX;
	for(i=1;i<=dim+1;i++)
	 tem[i]=fvertex[i];
     Camebsa(p,fvertex,dim,best,&best_f,ftol, XkGender,&iter,temptr, S0 ,fitdat,fitma,ma,flag,coeff,ekl,ekh,x,y,init);
	printf("%f\n",best_f);


   for(i=1;i<=dim+1;i++)
      {//printf("%f\n",FABS(tem[i]-fvertex[i]));
	printf("%d, %f, ",i, fvertex[i]);
	if(FABS(tem[i]-fvertex[i])==0)
	     flag1=1;
	else
	 { flag1=0;break;}
	  }
	printf("%d\n",flag1);
     if(flag1)
	{ for(j=1;j<=dim;j++)
	   p[1][j]=best[j];
       for(i=2;i<=dim+1;i++)
	{for(j=1;j<=dim;j++)
	do
	 {p[i][j]=best[j]+(ran1(&idum)-0.5)*temp1[j];
	 }while(init[j][1]>=p[i][j]||init[j][2]<=p[i][j]);
	}
       for(i=1;i<=dim+1;i++)
      { for(j=1;j<=dim;j++)
	 tem[j]=p[i][j];
      fvertex[i]=XkGender(x,y,tem,S0,fitdat,fitma,ma,flag,coeff,ekl,ekh);
      }
    min=fvertex[1];
    max=fvertex[1];
    jk=1;
    for(i=2;i<=dim+1;i++)
     {if(fvertex[i]<min)
       { min=fvertex[i];
	 jk=i;
       }
      if(fvertex[i]>max)
       {max=fvertex[i];}
     }

     T=max;
     temptr=T/(1.0-atm);
     }
      }while(iter<=0&&!kbhit());
    }
    else
   {ekl=2.0;ekh=18.0;}
    printf("Fitting is Done !\n");
	temp=0.0;
    while(x[i]<ekl)      i=i+1;
    while(x[i]<=ekh)
     {temp+=y[i]*y[i];
      i=i+1;
     }
   temp=sqrt(best_f/temp);
   printf("\n---------#The Related Error :%f  #---------\n",100.0*temp);getch();

    getch();
   save_dat(x,y,z,best,S0,fitdat,fitma,ma,flag,coeff,ekl,ekh);
  do{
    clrscr();
   printf("Imput the sourse file(*.par):");
   gets(stri);
      c = '.';
      ptr = strchr(stri, c);
      if (!ptr)
      strcat(stri,".par");
      path=searchpath(stri);
    }while((file=fopen(stri,"w+"))==NULL);
    fprintf(file,"fit %d series of paraments\n",fitma);
   for(i=0;i<fitma;i++)
    {fprintf(file,"the feff file name:%s\n",feff[2*i]);
     fprintf(file,"the feff file name:%s\n",feff[2*i+1]);
      fprintf(file,"the fit domain(k):  %f    %f \n",ekl,ekh);
     fprintf(file," the paraments:\n");
     if(ma==5)
       {fprintf(file,"Rj:  %f\n",best[i*ma+1]+(1.0+flag[i])*best[i*ma+5]);
	fprintf(file,"Rm:  %f\n",best[i*ma+1]+flag[i]*best[i*ma+5]);
       }
     else
       fprintf(file,"Rj:  %f\n",best[i*ma+1]);
     fprintf(file,"Nj:  %f\n",best[i*ma+2]);
     fprintf(file,"Sigma_t:  %f\n",best[i*ma+3]);
     fprintf(file,"E0:  %f\n",best[i*ma+4]);
     if(ma==5)
      fprintf(file,"Sigma_s:  %f\n",best[i*ma+5]);
    }
 fclose(file);
 plot(x,y,z,nl,nh,1);
 getch();
 closegraph();
 }


// * End of Fitxk2.cpp * //

