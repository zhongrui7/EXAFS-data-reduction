// Preedg.cpp //

// #include "exit.cpp"

 void Linear( )
 {
  float t,x1,u,u1,x2=0,x3=0,u3=0,u2=0,n=0;
  int i=0,j=0;

   while((X[i]-ef)<-300.0) {i=i+1;}
   while((X[i]-ef)<=-50.0)
    {t =X[i];
     u =Y[i];
     x1=t*t;
     u1=u*t;
     j =j+1;
     i =i+1;
     n =float(j);
     x2+=t;
     x3+=x1;
     u2+=u;
     u3+=u1;
    }
  d[0]=(u3*x2-x3*u2)/(x2*x2-n*x3);
  d[1]=(u2*x2-n*u3)/(x2*x2-n*x3);
  for(i=0;i<num;i++)
    {bkgd[i]=d[0]+d[1]*X[i]; }
   }


 void Victorian( )
 {float t,x1,x2,y1=0,y2=0,y3=0,y4=0,y5=0,DT=0;
  float t2, t3;
  int i=0;

   while((X[i]-ef)<-300.0) {i=i+1;}
   while(X[i]<ev)
    {
     t=10000.0/X[i];
     x1=t*t*t;
     x2=x1*t;
     y1+=x1*x1;
     y2+=x1*x2;
     y3+=x2*x2;
     y4+=x1*Y[i];
     y5+=x2*Y[i];
     i=i+1;
    }
  DT=y1*y3-y2*y2;
  if(DT==0.0){exit(1);}
  b[0]=(y4*y3-y5*y2)/DT;
  b[1]=(y5*y1-y4*y2)/DT;
 
  for(i=0;i<num;i++)
    {t=X[i];
     t=10000.0/t;
     t2=t*t*t;
     t3=t2*t;
     bkgd[i]=b[0]*t2+b[1]*t3;
     }
 }


 void getE0()
  {float d1, d2;
   int i, j;
  for(i=0; i<num; i++){EX0[i]=X[i]; EY0[i]=Y[i];}
  sprintf(label,"Energy [eV]");
    plot(X, Y, 0, 0, num);
    getch();
   d1=0.0;
   for(j=6;j<num/2;j++)
   {d2=(Y[j+3]-Y[j-3])/(X[j+3]-X[j-3]);
    if(d2>d1){deriv=j;  d1=d2; }
   }
    ef=X[deriv];
  }


void select_es()
  { float temp;
    char c;
    int i, j;

  sprintf(label,"Energy [eV]");
//* restore pre-edge background remved data *//
    for(i=0; i<num; i++){X[i]=EX0[i]; Y[i]=EY0[i];}  
    plot(X, Y, 0, 0, num);
    setfillstyle(1,0);
    setcolor(15);
    bar(10,430,620,473);
    outtextxy(120,450,"Press 'Space' to  <-|-> ");
    sprintf(string,"Select Pre-END & Post-START");
    outtextxy(350,450,string);
    outtextxy(550, 460,"'Enter'to process");
    c=getch();
    while(c!=32) {c=getch();}
    setcolor(15);
    setwritemode(XOR_PUT);
    sprintf(string,"%-7.2f",X[deriv]);
	moveto(250,70);
    outtext("E0 =");
    outtextxy(340,70,string);
    j=deriv;
    c=getch();
    while(c!=13)
     {
       if(c==75)	 { j=j-1;  }
       if(c==77)	 { j=j+1;  }
       if(j<6) j=6;
       if(j>num/2) j=num/2;
       i=80+int(unit_x*(X[j]-min_x)+0.5);
	line(i,400,i,80);
       ev=X[j];
       sprintf(string,"%-7.2f",ev-ef);
     setfillstyle(1,0);
     bar(100,70,180,90);
       moveto(110,80);
       outtext("Ev=");
       outtextxy(150,80,string);
	 c=getch();
     line(i,400,i,80);
      } 
     line(i,400,i,80);

        j=deriv;
	c=getch();
      while(c!=13)
     {
       if(c==75)	 { j=j-1;  }
       if(c==77)	 { j=j+1;  }
       if(j<9) j=9;
       if(j>num) j=num;
       i=80+int(unit_x*(X[j]-min_x)+0.5);
	line(i,400,i,80);
       es=X[j];
	sprintf(string,"%-7.2f",es-ef);
     setfillstyle(1,0);
     bar(400,70,500,90);
     moveto(420,80);
     outtext("ES=");
     outtextxy(460,80,string);
	c=getch();
     line(i,400,i,80);
      } 
  
    if(ev>es)	{temp=es;    es=ev;   ev=temp; }
  setfillstyle(1,0);
  setcolor(15);
  bar(10,430,620,473);
  outtextxy(500, 470,"press 'Esc' to exit");
 }


void sub()
{   char c, c1;
    int i, j, js, type=1;
    
    plot(X, Y, 0, 0, num);
     setfillstyle(1,0);
     setcolor(15);
     bar(10,430,620,473);
     outtextxy(150, 450,"Bkgrd: 0-Linear; 1-Victoreen");

     c1=getch();
     while(c1!=13)
      {  c=c1;
	while(c<48||c>51)
	{setcolor(14);
	bar(480,445,520,453);
	setcolor(14);
	 sprintf(string,"%c",c);
	 outtextxy(350,450,string);
	 c=getch();
	 }
     	 setcolor(15);
	 sprintf(string,"%c",c);
	 outtextxy(300,450,string);
    if(c=='1'||c=='0') type=c-48;
       c1=getch();
      }
  setfillstyle(1,0);
  bar(10,430,620,473);
  outtextxy(500, 470,"'Esc'exit");
 switch(type)
  {case 0:   
     Linear();
     break;
   case 1:  
     Victorian();
     break;
   }

 //* plot the raw data along with the pre-edge background *//
 plotD(EX0, EY0, EX0, bkgd, 0, 0, num); 
 for(i=0; i<num; i++)
      {
       EY[i]=Y[i]-bkgd[i];
      }
 getch();
 plot(EX0,EY, 0,0, num);

  i=0;
  while(es>X[i])  { i=i+1; }
  js= i;
  numb=num-js;
  for(i=0;i<numb;i++)
    { EX2[i]=X[i+js];
      EY2[i]=Y[i+js]-bkgd[i+js];
     }
  getch();
  setfillstyle(1,0);
  bar(10,430,620,473);
  plot(EX2, EY2, 0, 0, numb);
  outtextxy(500, 470,"'Esc'exit");
  ntag[1]=numb;
}
/*///End of PREEDG.cpp///*/