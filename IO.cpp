//* IO.cpp *//


// void ReadInit(float S0, int fitma, int ma, char ikk[40], int coeff, float **init)
void ReadInit()
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
  do{
   clrscr();
   printf("********************************************************************** \n");
   printf("*                  USTC EXAFS FITTING PROGRAM 0.1a                    * \n");
   printf("*           ( Written by Xiaohong Wan, Zr Li, 2009/01/19 )           * \n");
   printf("* ------------------------------------------------------------------ * \n");
   printf("*  This program needs an Input file (ASCII.ini) like following       * \n");
   printf("*  S0      0.7-1.0 ; Amplitude Reduction Factor                      * \n");
   printf("*  fitma   1 - 5   ; Shell Number to be fitted                       * \n");
   printf("*  ma      4 - 5   ; Parameter Number                                * \n");
   printf("*  ikk     XXX.ikk ; Data in k-Space to be fitted [N / K Chi]        * \n");
   printf("*  coeff   0 - 3   ; k-weight in Fourier transform                   * \n");
   printf("*  feff    AAA.dat ; FEFF file for Absorption Atom                   * \n");
   printf("*  feff    BBB.dat ; FEFF file for Scattering Atom                   * \n");
   printf("*  flag    -1 - 3  ; EXAFS Function [-1:Gaussian, 0-3:Weighted Exp]  * \n");
   printf("*  R       x.x y.y ; Distance Range (Angstrom)                       * \n");
   printf("*  N       x.x y.y ; Coordination Number Range                       * \n");
   printf("*  sigma   x.x y.y ; Thermal Disorder (.000-0.500 Angstrom)          * \n");
   printf("*  E0      x.x y.y ; Energy Shift Range (-20.~20. eV)                * \n");
   printf("*  s       x.x y.y ; Static Disorder (.000-0.500 Angstrom)           * \n");
   printf("********************************************************************** \n");
   printf("Imput the sourse file(*.ini):");
   gets(string);
      c = '.';
      ptr = strchr(string, c);
      if (!ptr)
      strcat(string,".ini");
      path=searchpath(string);
    }while(((handle=open(path,O_RDONLY))==-1)||(file=fdopen(handle,"r+"))==NULL);

   fscanf(file,"%s",buff);
   fscanf(file,"%s",buff);
   S0=atof(buff);
   fscanf(file,"%s",buff);
   fscanf(file,"%s",buff);
   fitma=atoi(buff);
   fscanf(file,"%s",buff);
   fscanf(file,"%s",buff);
   ma=atoi(buff);
   dim=ma*fitma;

   p=matrix(1,ma*fitma+1,1,ma*fitma);
   init=matrix(1,ma*fitma,1,2);
   temp1=vector(1,ma*fitma);
   fvertex=vector(1,ma*fitma+1);
   best=vector(1,ma*fitma);
   tem=vector(1,ma*fitma+1);
   fscanf(file,"%s",buff);
   fscanf(file,"%s",ikk);
   fscanf(file,"%s",buff);
   fscanf(file,"%s",buff);
   coeff=atoi(buff);
   for(j=0,i=0;i<fitma;i++)
    {fscanf(file,"%s",buff);
     fscanf(file,"%s",feff[j]);
      j=j+1;
     fscanf(file,"%s",buff);
     fscanf(file,"%s",feff[j]);
      j=j+1;
     fscanf(file,"%s",buff);
     fscanf(file,"%s",buff);
     flag[i]=atof(buff);
     fscanf(file,"%s",buff);
     fscanf(file,"%s",buff);
     init[i*ma+1][1]=atof(buff);
     fscanf(file,"%s",buff);
     init[i*ma+1][2]=atof(buff);
     fscanf(file,"%s",buff);
     fscanf(file,"%s",buff);
     init[i*ma+2][1]=atof(buff);
     fscanf(file,"%s",buff);
     init[i*ma+2][2]=atof(buff);
     fscanf(file,"%s",buff);
     fscanf(file,"%s",buff);
     init[i*ma+3][1]=atof(buff);
     fscanf(file,"%s",buff);
     init[i*ma+3][2]=atof(buff);
     fscanf(file,"%s",buff);
     fscanf(file,"%s",buff);
     init[i*ma+4][1]=atof(buff);
     fscanf(file,"%s",buff);
     init[i*ma+4][2]=atof(buff);
     if(ma==5)
      {fscanf(file,"%s",buff);
       fscanf(file,"%s",buff);
       init[i*ma+5][1]=atof(buff);
       fscanf(file,"%s",buff);
       init[i*ma+5][2]=atof(buff);
      }
    }
}


void nihe(struct FIT y[],int pnum)
    { float *Y1,*Y2,*Y3,*Y4,*Y5,*X1,*X2,*X3,*X4,*X5,*X6;
      float yp1,ypn, xmax,deta;
      int i,j;
       Y1=new float[pnum];
       Y2=new float[pnum];
       Y3=new float[pnum];
       Y4=new float[pnum];
       Y5=new float[pnum];
       X1=new float[pnum];
       X2=new float[pnum];
       X3=new float[pnum];
       X4=new float[pnum];
       X5=new float[pnum];
       X6=new float[pnum];

      for(i=0;i<pnum;i++)
       {X1[i]=y[i].k;
	X2[i]=y[i].phc;
	X3[i]=y[i].mag;
	X4[i]=y[i].phase;
	X5[i]=y[i].factor;
	X6[i]=y[i].lamda;

	Y1[i]=0.0;
	Y2[i]=0.0;
	Y3[i]=0.0;
	Y4[i]=0.0;
	Y5[i]=0.0;

       }
      xmax=y[pnum-1].k;
      deta=0.05;
      yp1=1.0e30;
      ypn=1.0e30;
      spline(X1,X2,pnum,yp1,ypn,Y1);
      spline(X1,X3,pnum,yp1,ypn,Y2);
      spline(X1,X4,pnum,yp1,ypn,Y3);
      spline(X1,X5,pnum,yp1,ypn,Y4);
      spline(X1,X6,pnum,yp1,ypn,Y5);

       j=-1;
       do
	 { j=j+1;
	   y[j].k=deta*double(j);
	   y[j].phc=splint(X1,X2,Y1,pnum,y[j].k);
	   y[j].mag=splint(X1,X3,Y2,pnum,y[j].k);
	   y[j].phase=splint(X1,X4,Y3,pnum,y[j].k);
	   y[j].factor=splint(X1,X5,Y4,pnum,y[j].k);
	   y[j].lamda=splint(X1,X6,Y5,pnum,y[j].k);

	}while(y[j].k!=xmax);
	delete[] Y1;
	delete[] Y2;
	delete[] Y3;
	delete[] Y4;
	delete[] Y5;
	delete[] X1;
	delete[] X2;
	delete[] X3;
	delete[] X4;
	delete[] X5;
	delete[] X6;

	}


int comp(char *buf)
  {char string1[15]="real[p]",string2[15]="real[p]@#";
   int l1,l2;
   l1=strcmp(string1,buf);
   l2=strcmp(string2,buf);
   l1=l1*l2;
   return l1;
  }


void  ReadSource(char feff[][40],struct FIT fitdat[][420],int fitma)
 {FILE *file;
  int i,j,m,handle,num=0;
  char c,*ptr,*path, buff[80];
  for(m=0,i=0;i<fitma;i++)
   {
      c = '.';
      ptr = strchr(feff[m], c);
      if (!ptr)
      strcat(feff[m],".inp");
      printf("%s\n",feff[m]);
      path=searchpath(feff[m]);
      if((handle=open(path,O_RDONLY))==-1||(file=fdopen(handle,"r+"))==NULL)
       { printf("Open file:%s error",feff[m]);
	exit(1);
       }
    do{
     fscanf(file,"%s",buff);
      }while(comp(buff)!=0);
     j=-1;
   do
    { j=j+1;
    fscanf(file,"%s",buff);
    fitdat[i][j].k=atof(buff);
    fscanf(file,"%s",buff);
    fitdat[i][j].phc=atof(buff);
    fscanf(file,"%s",buff);
    fscanf(file,"%s",buff);
    fscanf(file,"%s",buff);
    fscanf(file,"%s",buff);
    fscanf(file,"%s",buff);
    } while(fitdat[i][j].k!=20.0);
   num=j+1;
   fclose(file);
      m=m+1;
      c = '.';
      ptr = strchr(feff[m], c);
      if (!ptr)
      strcat(feff[m],".inp");
      printf("%s\n",feff[m]);
      path=searchpath(feff[m]);
      if((handle=open(path,O_RDONLY))==-1||(file=fdopen(handle,"r+"))==NULL)
       { printf("open file:%s error",feff[m]);
	exit(1);
       }
    do{
     fscanf(file,"%s",buff);
      }while(comp(buff)!=0);
     j=-1;
   do
    { j=j+1;
    fscanf(file,"%s",buff);
    fscanf(file,"%s",buff);
    fscanf(file,"%s",buff);
    fitdat[i][j].mag=atof(buff);
    fscanf(file,"%s",buff);
    fitdat[i][j].phase=atof(buff);
    fscanf(file,"%s",buff);
    fitdat[i][j].factor=atof(buff);
    fscanf(file,"%s",buff);
    fitdat[i][j].lamda=atof(buff);
    fscanf(file,"%s",buff);
    } while(fitdat[i][j].k!=20.0);
   fclose(file);
   m=m+1;
  nihe(fitdat[i],num);
  }
}



int ReadOrigin(char ikk[],float x[],float xkk[])
  {FILE *file;
  int j,handle,num;
  char c,*ptr,*path, buf1[40],buf2[40],buf3[40];
   c = '.';
   ptr = strchr(ikk, c);
   if (!ptr)
      strcat(ikk,".ikk");
   printf("%s\n",ikk);
   path=searchpath(ikk);
   if((handle=open(path,O_RDONLY))==-1||(file=fdopen(handle,"r+"))==NULL)
     { printf("Open file:%s error\n",ikk);
	exit(1);
      }

   fscanf(file,"%s", buf1);
   num=atoi(buf1);
   j=0;
   do
   {j=j+1;
   fscanf(file,"%s",buf2);
   fscanf(file,"%s",buf3);
   x[j]=atof(buf2);
   xkk[j]=atof(buf3);
   }while(x[j]<20.0);
  fclose(file);
  num=j;
  return(num);
 }


void save_dat(float x[],float y[],float z[],float parm[],float S0,struct FIT temp[][420],int fitma,int ma,int flag[],int coeff,float ekl,float ekh)
      {int i,j,m;
	float Ak,Fk,Qk,Mk,Nk,Sk;
	char c,*ptr, string[20];
	FILE *fout;
           
	for(j=1;j<=420;j++)
	   z[j]=0.0;
	   j=1;
       while(x[j]<ekl) {j=j+1;}
       while(x[j]<=ekh)
	 { z[j]=0.0;
	     Nk=1.0;
	    for(m=0;m<coeff;m++) Nk=Nk*x[j];
	    for(i=0;i<fitma;i++)
	  {
	 switch(flag[i]){
           case -1:
		 {//x[j]=sqrt(x[j]*x[j]-0.263*parm[i*ma+4]);
		  Qk= 2*x[j]*parm[i*ma+1]+temp[i][j-1].phc+temp[i][j-1].phase+0.263*parm[i*ma+4]*parm[i*ma+1]/x[j];
		  Fk=(temp[i][j-1].mag*temp[i][j-1].factor*S0/x[j])*exp(-2*parm[i*ma+3]*parm[i*ma+3]*x[j]*x[j]);
		  Mk=Fk*exp(-2*parm[i*ma+1]/temp[i][j-1].lamda)/(parm[i*ma+1]*parm[i*ma+1]);
		  z[j]+=Nk*parm[i*ma+2]*Mk*sin(Qk);
		 }
	       break;
           case 2:
 	    {// x[j]=sqrt(x[j]*x[j]-0.263*parm[i*ma+4]);
	     Sk=4*x[j]*x[j]*parm[i*ma+5]*parm[i*ma+5];
	    Qk= 2*x[j]*parm[i*ma+1]+temp[i][j-1].phc+temp[i][j-1].phase+0.263*parm[i*ma+4]*parm[i*ma+1]/x[j]+atan(2*x[j]*parm[i*ma+5]*(3.0-Sk)/(1.0-3*Sk));
	    Fk=(temp[i][j-1].mag*temp[i][j-1].factor*S0/x[j])*exp(-2*parm[i*ma+3]*parm[i*ma+3]*x[j]*x[j]);
	    Mk=Fk*exp(-2*parm[i*ma+1]/temp[i][j-1].lamda)/(parm[i*ma+1]*parm[i*ma+1]);
	    Ak=Mk/((1+Sk)*sqrt(1+Sk));
	    z[j]+=Nk*parm[i*ma+2]*Ak*sin(Qk);
	    }
	    break;
           case 1:
	    {
	     Sk=4*x[j]*x[j]*parm[i*ma+5]*parm[i*ma+5];
	    Qk= 2*x[j]*parm[i*ma+1]+temp[i][j-1].phc+temp[i][j-1].phase+0.263*parm[i*ma+4]*parm[i*ma+1]/x[j]+atan(2*x[j]*parm[i*ma+5]*(2.0-Sk)/(1.0-2*Sk));
	    Fk=(temp[i][j-1].mag*temp[i][j-1].factor*S0/x[j])*exp(-2*parm[i*ma+3]*parm[i*ma+3]*x[j]*x[j]);
	    Mk=Fk*exp(-2*parm[i*ma+1]/temp[i][j-1].lamda)/(parm[i*ma+1]*parm[i*ma+1]);
	    Ak=Mk/(1.0+Sk);
	    z[j]+=Nk*parm[i*ma+2]*Ak*sin(Qk);
	    }
	    break;
           case 0:
	    {
	     Sk=4*x[j]*x[j]*parm[i*ma+5]*parm[i*ma+5];
	    Qk= 2*x[j]*parm[i*ma+1]+temp[i][j-1].phc+temp[i][j-1].phase+0.263*parm[i*ma+4]*parm[i*ma+1]/x[j]+atan(2*x[j]*parm[i*ma+5]);
	    Fk=(temp[i][j-1].mag*temp[i][j-1].factor*S0/x[j])*exp(-2*parm[i*ma+3]*parm[i*ma+3]*x[j]*x[j]);
	    Mk=Fk*exp(-2*parm[i*ma+1]/temp[i][j-1].lamda)/(parm[i*ma+1]*parm[i*ma+1]);
	    Ak=Mk/(sqrt(1.0+Sk));
	   z[j]+=Nk*parm[i*ma+2]*Ak*sin(Qk);
	   }
	    break;
	   }
	 z[j]=z[j];
	}
	j=j+1;
      }
	 do{
   clrscr();
   printf("Imput the out file(*.dat):");
   gets(string);
      c = '.';
      ptr = strchr(string, c);
      if (!ptr)
      strcat(string,".dat");
    }while((fout=fopen(string,"w+"))==NULL);

	i=1;
    while(x[i]<ekl)   {i=i+1;}
    while(x[i]<=ekh)
   {
   fprintf(fout,"% -13.8f     %-13.8f    %-13.8f\n",x[i],y[i],z[i]);
     i=i+1;
   }
   fclose(fout);
  }

// End of IO.cpp //
