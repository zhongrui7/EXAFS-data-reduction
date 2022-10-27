
//* Plot.cpp *//

 void iniwindow ()
  {int gdriver = DETECT, gmode;
   initgraph(&gdriver, &gmode, "");
    setbkcolor(0);
    setfillstyle(2,6);
    setcolor(14);
    bar(10,50,620,430);
    settextstyle(10,0,3);
    outtextxy(150,40,"EXAFS Fitting 0.1c--*--(NSRL/USTC, 2018/08/19)");
     setcolor(15);
     settextstyle(0,0,0);
     moveto(60,400);
     line(60,400,620,400);
     line(60,400,60,50);
   }


  void plot2( float Y5[],int n,float *unit_x)
    { float max_a,min_a,unit_y, digit1,digit2;
       char *string=NULL;

     int i,x,y,x0,y0, Y;
     iniwindow();
     max_a=min_a=Y5[1];
     setwritemode(COPY_PUT);
     setcolor(14);
     for(i=2;i<=n;i++)
     {if(Y5[i]>max_a)  max_a=Y5[i];
      if(Y5[i]<min_a)  min_a=Y5[i];
     }

      unit_y=320.0/(max_a-min_a);
     *unit_x=500.0/n;
	i=1;
	x0=60+int(i*(*unit_x)+0.5);
	y0=400-int((Y5[i]-min_a)*(unit_y)+0.5);
		Y=400-int((-min_a)*(unit_y)+0.5);
		line(60,Y,600,Y);
		outtextxy(25,Y,"0");
      for(i=2;i<=n;i++)
      {
	x=60+int(i*(*unit_x)+0.5);
	y=400-int((Y5[i]-min_a)*(unit_y)+0.5);
	 line(x0,y0,x,y);
	 x0=x;
	 y0=y;
	}

     for(i=0;i<6;i++)
     {  digit1=(max_a-min_a)*(i-2.5)/5;
	digit2=i*(0.05*double(n)/5);
	sprintf(string,"%5.1f",digit1);
       x=12;
       y=Y-64*(i-2.5);
     outtextxy(60,y,"-");
     outtextxy(x,y,string);
     sprintf(string,"%-5.1f",digit2);
     x=56+100*i;
     y=415;
     outtextxy(x-1,400,"|");
     outtextxy(x,y,string);
     }
   }


void plot(float x[],float y[],float z[],int nl,int nh ,int l)
    {
     float digit1,digit2;
     char *string=NULL;
     float max_a,min_a,max_e,min_e,unit_x,unit_y;
    int i,x1,y1,x0,y0,y2;
    iniwindow();
     setwritemode(XOR_PUT);
     max_a=min_a=y[nl];
     max_e=min_e=x[nl];

     for(i=nl+1;i<=nh;i++)
     {if(y[i]>max_a)  max_a=y[i];
      if(y[i]<min_a)  min_a=y[i];
      if(x[i]>max_e)  max_e=x[i];
      if(x[i]<min_e)  min_e=x[i];
      }

     unit_y=320.0/(max_a-min_a);
     unit_x=500.0/(max_e-min_e);

	i=nl;
	x0=60+int((x[i]-min_e)*unit_x+0.5);
	y0=400-int((y[i]-min_a)*unit_y+0.5);
       setcolor(15);
      for(i=nl;i<=nh;i++)
      {
	x1=60+int((x[i]-min_e)*unit_x+0.5);
	y1=400-int((y[i]-min_a)*unit_y+0.5);
	y2=400-int((z[i]-min_a)*unit_y+0.5);
	 line(x0,y0,x1,y1);
	 x0=x1;
	 y0=y1;
	 putpixel(x1,y2,1);
	circle(x1,y2,1);
	}
     if(l==1)
      {y1=400-int((0-min_a)*unit_y+0.5);
       line(60,y1,600,y1);
       outtextxy(25,y1,"0");
      }
     digit1=min_a;
     digit2=min_e;
     setcolor(14);
     settextstyle(0,0,0);
     for(i=0;i<6;i++)
     {  digit1=(max_a-min_a )*(i-2.5)/5;
	digit2=(max_e-min_e)*i/5+min_e;
       sprintf(string,"%5.1f",digit1);
       x0=12;
       y0=y1-64*(i-2.5);
     outtextxy(60,y0,"-");
     outtextxy(x0,y0,string);
     sprintf(string,"%-5.1f",digit2);
     x0=60+100*i;
     y0=415;
     outtextxy(x0-1,400,"|");
     outtextxy(x0,y0,string);
     }
      }


  void select_es(float y[],int wmax,float *ekl,float *ekh)
    {
    float unit,ek[2], ke, temp;
    char c, string[40];
    int i,j=0, delta=16;
     plot2(y,wmax,&unit);
     setcolor(15);
     outtextxy(100,450,"1. 'Space' to select Fitting Regime:[USE <-A,S]|[D,F->]");
     outtextxy(100,470,"2. 'Esc' to fit");

      c=getch();
      while(c!=32)
      {c=getch();}
       setcolor(15);
      setwritemode(XOR_PUT);

  do{
      i=60+250;
      line(i,400,i,80);
       ke=float(i-60)*0.05/unit;
       sprintf(string,"%-5.2f",ke);
       moveto(80,70);
       outtext("k=");
       outtextxy(110,70,string);

      c=getch();
      while(c!=13)
      {
       if(c=='A'||c=='a')	 { i=i-delta;  }
       if(c=='F'||c=='f')	 { i=i+delta;  }
       if(c=='S'||c=='s')	 { i=i-1;  }
       if(c=='D'||c=='d')	 { i=i+1;  }
       if(i<60)i=60;
       if(i>600)i=600;
       line(i,400,i,80);
       ek[j]=float(i-58)*0.05/unit;
       sprintf(string,"%-5.2f",ek[j]);
       moveto(80,70);
       outtext("K=");
       outtextxy(110,70,string);
	c=getch();
       line(i,400,i,80);
       setfillstyle(1,2);
       bar(70,60,170,80);
       }
       line(i,400,i,80);
    j++;
    } while (j<2);


       if(ek[0]>ek[1])
	{temp=ek[0];
	 ek[0]=ek[1];
	 ek[1]=temp;
	}
       *ekl=ek[0];
       *ekh=ek[1];
      getch();
      closegraph();
    }


// * End of plot.cpp   *//