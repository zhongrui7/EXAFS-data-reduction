// FFT.cpp //
//ll=1, forward FT; ll=-1, inverse transform //
void zrfft(float *re, float *im, int nn, int ll)
{
 int j, ij, m, l, istep;
 float sc, actep, bctep, bcarg, aw, bw;
 float atemp, btemp;
 j = 0;
 sc = sqrt(1.0 / nn);

for (ij=0; ij<=nn-1; ij++)
 {
  if (ij <= j)
   {
   actep = re[j] * sc;
   bctep = im[j] * sc;
   re[j] = re[ij] * sc;
   im[j] = im[ij] * sc;
   re[ij] = actep;
   im[ij] = bctep;
   }
 m = nn/2;
 do
  {  if ((j+1)<=m)   break;
   j = j - m;
   m = m/2;
   }while (m>=1);
   j = j + m;
  }

 l = 1;
 do
  {  istep = 2 * l;
  for (m=1; m<=l; m++)
   {
    bcarg = -(PI * ll * (m-1)) / l;
    aw = cos(bcarg);
    bw = sin(bcarg);
    for (ij=m-1; ij<=nn-1; ij=ij+istep)
    { atemp = aw * re[ij+l] - bw * im[ij+l];
     btemp = aw * im[ij+l] + bw * re[ij+l];
     re[ij+l] = re[ij] - atemp;
     im[ij+l] = im[ij] - btemp;
     re[ij] = re[ij] + atemp;
     im[ij] = im[ij] + btemp;
     }
   }
 l = istep;
 }while( l < nn );
}


// * End of fft.cpp *//


 void transform(int wmax)
  {  complex *temp;
     int i,j,wmax1, delta;
     float rd,yd,yp1,ypn, tmp;
     float *X4, *Y4, *Y5, *Y6, *F;
     char c,c1, string1[32]; 
     temp=new complex[4*wmax];
     X4=new float[4*wmax];
     Y4=new float[4*wmax];
     Y5=new float[4*wmax];
     Y6=new float[4*wmax];
     F=new float[wmax];
 
//* restore EXAFS data from KX[], KY[] *//
    for(i=0; i<wmax; i++){X4[i]=KX[i]; Y4[i]=KY[i]; Y6[i]=0.0;}  
    for(i=wmax;i<4*MAXNM;i++)	{ Y4[i]=0.0; Y6[i]=0.0;}
	
	sprintf(label,"k [A^-1]");
      plot(X4, Y4, 0, 20, wmax);
     setfillstyle(1,0);
     bar(10,430,620,470);
      setcolor(15);
     outtextxy(150, 450,"Enter K-weght [n=0/1/2/3]:");
     outtextxy(300, 470,"'Esc'exit");
     c1=getch();
     while(c1!=13)
      {  c=c1;
	while(c<48||c>57)
	{setcolor(14);
	bar(480,445,520,453);
	setcolor(4);
	 sprintf(string,"%c",c);
	 outtextxy(350,450,string);
	 c=getch();
	 }
	 sprintf(string,"%c",c);
	 outtextxy(300,450,string);
       if(c=='1'||c=='2'||c=='3'||c=='0')
       kweight=c-48;
       c1=getch();
      }
  
	for(i=0;i<wmax;i++)
	 {for(j=0;j<kweight;j++) {Y4[i]=Y4[i]*X4[i]; }
          KYW[i]=Y4[i];
          }
       exitwindow();
       plot(X4, Y4, 0, 20, wmax);
       sprintf(string,"%d",kweight);
       outtextxy(300,450,string);
     setfillstyle(1,0);
     bar(10,430,620,470);
       setcolor(15);     
       outtextxy(150, 450,"FT range [k=2.0,15.0]:");
     ev=2.0; es=15.0;
     c=getch();
       string=NULL;
   while(c!=13){
     j=0;
      while(c!=',')
	{
	 sprintf(string,"%c",c);
	 outtextxy(300+j*8,450,string);
         strcat(string1,string);
	 c=getch();j++;
	 }
        ev=atof(string1);
       
      string=NULL; string1[0]='\0';
      j=0;
      while(c!=13)
      {
        sprintf(string,"%c",c);
	outtextxy(350+j*8,450,string);
        strcat(string1,string);
        c=getch();j++;
      }
        es=atof(string1);

      if(ev>es)	{tmp=es; es=ev; ev=tmp;}
   ev=2.0; es=15;
    }

     for(i=0;i<wmax;i++)
      {if(KX[i]<ev-1) F[i]=0.0; 
       if(KX[i]>=ev-1&&KX[i]<ev)F[i]=cos((ev*20-i)/40*PI)*cos((ev*20-i)/40*PI);
       if(KX[i]>=ev&&KX[i]<=es)   F[i]=1.0;
       if(KX[i]>es&&KX[i]<=es+1)F[i]=cos((i-es*20)/40*PI)*cos((i-es*20)/40*PI);
       if(KX[i]>es+1) F[i]=0.0;
       
     Y5[i]=F[i]*Y4[i]; 
     for(j=1;j<4;j++) {Y5[i+j*wmax]=0.0;  }
     }
     getch();
     plotD(KX, Y4, KX, Y5, 0, 20, wmax);

	
   zrfft(Y5, Y6, 4*wmax, 1);

//*backup Fourier Transform result*/
    for(i=0;i<4*wmax;i++)
     {RYr[i]=Y5[i]; RYi[i]=Y6[i]; 
      temp[i]=complex(Y5[i], Y6[i]);
      Y4[i]=abs(temp[i]);
      X4[i]=arg(temp[i]);
     }

    for(i=0;i<wmax;i++) 
     {RX[i]=(float(i)/1024.0)*31.4; RY[i]=Y4[i]; }
    sprintf(label,"R [A]");
     getch();
     plot(RX, RY, 0, 10, wmax);
     outtextxy(400, 100,"K-weight:");
	 sprintf(string,"%d",kweight);
	 outtextxy(450,100,string);
     getch();
     plotD(RX, RYr, RX, RYi, 0, 10, wmax);
 int x, y;
 for(i=0;i<wmax/2;i++)
      { 
	x=80+int((RX[i]-min_x)*unit_x+0.5);
	y=400-int((RY[i]-min_y)*unit_y+0.5);
	putpixel(x,y,15);
      }
ntag[3]=wmax;
    delete[] X4;     
    delete[] Y4;  
    delete[] Y5;     
    delete[] Y6;  
    delete[] temp;
  }


void inverse(int wmax)
     { char c,c1;
       int i,j,k,p, wax, wv,ws,qa;
       float max_dot,rd,tempwan,tmp;
       float *Re, *Im,*Y7,*Y8;
      complex *temp;
      temp=new complex[MAXNM];
       Re=new float[4*wmax];
       Im=new float[4*wmax];
       Y7=new float[4*wmax];
       Y8=new float[4*wmax];
      for(i=0;i<4*wmax;i++) {Re[i]=RYr[i]; Im[i]=RYi[i];}

      sprintf(label, "R [A]");
      setcolor(COL);
      plot(RX, RY, 0, 10, MAXNM);

   tmp=RY[0];
   for(i=1;i<MAXNM/4;i++)
      {	if(RY[i]>tmp) {tmp=RY[i]; p=i;}
       }

     setwritemode(XOR_PUT);
     setcolor(15);
     i=80+int(RX[p]*unit_x+0.5); 
     line(i,400,i,80);
     sprintf(string,"%-4.2f",RX[p]);
     moveto(280,70);
     outtext("r=");
     outtextxy(320,70,string);
     setcolor(15);
     setfillstyle(1,0);
     bar(10,430,620,475);
     setcolor(15);
     outtextxy(340,450,"To Back FT this peak [Y/N]:");
     outtextxy(200, 470,"'Esc'exit");
     c=getch();
      if(c=='Y'||c=='y'||c=='N'||c=='n')
      {sprintf(string,"%c",c);
       outtextxy(460,450,string);
       }
  if(c=='Y'||c=='y')
    {
     setfillstyle(1,0);
     bar(10,430,620,475);
     setcolor(15);
     outtextxy(340,450,"Put 'Space' to select :[<-|->]");
       c=getch();
     while(c!=32) {c=getch();}  //press SPACE
     setcolor(7);
     k=p;      
     //Begin of peak selection zone
     c=getch();
     while(c!=13)  //press ENTER
     {
       if(c==75)	 { k=k-1;  }
       if(c==77)	 { k=k+1;  }
       if(k<1)k=1;
       if(k>MAXNM/2)k=MAXNM/2;
       i=80+int(RX[k]*unit_x+0.5); 
       line(i,400,i,80);
       wv=k;
       sprintf(string,"%-8.4f",RX[wv]);
     setfillstyle(1,0);
     bar(10,430,620,470);
       moveto(280,70);
       outtext("r=");
       outtextxy(320,70,string);
	c=getch();
       }
       sprintf(string,"%-4.2f",RX[wv]);
       moveto(120,70);
       outtext("Left=");
       outtextxy(160,70,string);
       
       k=p;
       c=getch();
      while(c!=13)
     { if(c==75)	 { k=k-1;  }
       if(c==77)	 { k=k+1;  }
       if(k<1)k=1;
       if(k>MAXNM/2)k=MAXNM/2;
       j=80+int(RX[k]*unit_x+0.5); 
       line(j,400,j,80);
       ws=k;
       sprintf(string,"%-8.4f",RX[ws]);
	setfillstyle(1,0);
        bar(250,60,350,80);
	moveto(280,70);
       outtext("r=");
       outtextxy(320,70,string);
	  c=getch();
        }
       sprintf(string,"%-4.2f",RX[ws]);
       moveto(450,80);
       outtext("Righ=");
       outtextxy(500,80,string);

   setfillstyle(EMPTY_FILL,7);
   rectangle(i,80,j,400);

      if(wv>ws)	{tmp=ws; ws=wv; wv=tmp;}
      //End of peak selection zone

       qa=5; 	
      for(i=0;i<4*wmax;i++)
	{ 
	tempwan=0.0;
	if(i>=wv&&i<=ws)
	  {if(i<=wv+qa)  tempwan=0.5*(1.0+cos((i-wv-qa)*PI/qa));
	   if(i>=ws-qa)  tempwan=0.5*(1.0+cos((ws-qa-i)*PI/qa));
	   else	       tempwan=1;
	    }

	if(i>=(4*wmax-ws)&&i<=(4*wmax-wv))
	   {if(i<=(4*wmax-ws)+qa)
	       tempwan=0.5*(1.0+cos((i-(4*wmax-ws)-qa)*PI/qa)); 
	    if(i>=(4*wmax-wv)-qa)
	       tempwan=0.5*(1.0+cos(((4*wmax-wv)-qa-i)*PI/qa)); 
	    else tempwan=1;   
           }
	 Re[i]=Re[i]*tempwan; Y7[i]=Re[i]; 
         Im[i]=Im[i]*tempwan; Y8[i]=Im[i];
         }
    
       for(k=0; k<wmax; k++)
           {BX[k]= 25.6*(float(k)/float(wmax)); }

        zrfft(Y7, Y8, 4*wmax, -1);
	exitwindow(); 
	sprintf(label,"k [A^-1]");
          plotD(BX, Y7, BX, Y8, 1, 20, wmax);
        outtextxy(500, 100,"K-weight:");
	 sprintf(string,"%d",kweight);
	 outtextxy(550,100,string);
      
       for(i=10;i<wmax;i++)
	 { 
          for(j=1; j<=kweight; j++){ Y7[i]=Y7[i]/BX[i]; Y8[i]=Y8[i]/BX[i];}
          BY[i]=Y7[i];  
          }

	getch();
   plot(BX, BY, 1, 21, wmax);
   outtextxy(400, 100,"Chi(k)");
	getch();
   plotD(BX, BY, KX, KY, 1, 21, wmax);
	setcolor(15);
	outtextxy(500, 100,"RAW"); 

	getch();
      for(i=0;i<2*wmax;i++) {Y7[i]=Re[i]; Y8[i]=Im[i];}
      for(i=2*wmax;i<4*wmax;i++) {Y7[i]=0.0; Y8[i]=0.0;}
        zrfft(Y7, Y8, 4*wmax, -1);

	tmp=0.0; j=0.0;
 for(i=1; i<MAXNM; i++)
   { 
	for(k=1; k<=kweight; k++){ Y7[i]=Y7[i]/BX[i]; Y8[i]=Y8[i]/BX[i];}
	temp[i]=complex(Y7[i], Y8[i]);
     Amp[i]=2*abs(temp[i]);
     Pha[i]=arg(temp[i]);
     if(tmp>0&&Pha[i]<0) j++;  
     tmp=Pha[i];
     Pha[i]+=j*2*PI;

//    Pha[i]=180/PI*(Pha[i]-2*BX[i]*RX[p]); 
//    Pha[i]=Pha[i]-2*BX[i]*RX[p];
    }

     plotD(BX, BY, BX, Amp, 1, 21, wmax); 
	setcolor(14);
	outtextxy(500, 100,"Amplitude");
	setcolor(15);
	outtextxy(500, 120,"Chi(k)");
	getch();
     plot(BX, Pha, 1, 21, wmax); 
	setcolor(15);
	outtextxy(500, 100,"Phase"); 
    }


   else
     { setfillstyle(1,8);
       bar(70,60,170,80);
       setcolor(7);
       line(i,400,i,80);
      }

     setwritemode(COPY_PUT);
     setcolor(7);
     delete[] temp;
     delete[] Re;
     delete[] Im;
     delete[] Y7;
     delete[] Y8;
 
    }
/*///End of FFT.cpp////*/



