// PLOT.cpp //
#define COL 14
 float unit_x,unit_y;  
 float max_y,max_x,min_y,min_x;

 void iniwindow ()
 { setbkcolor(0);
   setcolor(14);
   setwritemode(COPY_PUT);
   rectangle(5,41,625,475);
   setfillstyle(2,8);
   bar(7,43,624,430);
   moveto( 60,400);
   line(80,400,620,400); // X-axis 520 pixels
   line(80,400,80,60);   // Y-axis 340 pixels 
  }


void frame()
{ int i, x,y0, y; 
  float tag1,tag2;
  setcolor(15);
  settextstyle(0,0,0);
  for(i=0;i<6;i++)
    { 
     tag2=(max_x-min_x)*i/5+min_x;
     sprintf(string,"%-6.2f",tag2);
     x=80+100*i;
     y=415;    
     outtextxy(x,y,string);   
     outtextxy(x-1,400,"|");
     }

  if(max_y>0.000001&&min_y<-0.000001)
   {y0=400-int((0-min_y)*unit_y+0.5);
    setcolor(YELLOW); line(90,y0,600,y0);  
    for(i=-6;i<6;i++)
     { tag1=(max_y-min_y)*i/5;
     if(tag1<max_y&&tag1>min_y) 
      {sprintf(string,"%-6.3f",tag1);
       x=50;
       y=y0-tag1*unit_y;
      setcolor(YELLOW); outtextxy(80,y,"-");
      setcolor(WHITE); outtextxy(x,y,string);}
      }
     }

  else
   {for(i=0;i<6;i++)
     {tag1=(max_y-min_y )*i/5+min_y;
      sprintf(string,"%-6.3f",tag1);
      x=50;
      y=400-64*i;
     outtextxy(80,y,"-");
     outtextxy(x,y,string);
      }
    }
 outtextxy(550,65,filename);
 outtextxy(550,390,label);
 }


void plot(float X[], float Y[],float minx, float maxx, int n)
  {
   int i, j, k, x,y,x0,y0;
   iniwindow();
   if(minx==maxx)
   {
    max_x=min_x=X[0];
    for(i=1;i<n;i++)
     {
      if(X[i]>max_x)  max_x=X[i];
      if(X[i]<min_x)  min_x=X[i];
      }
    j=1; k=n;
    }

  else
    { j=1;
     for(i=1;i<n;i++)
     {if(X[i]<=minx)  j=i; 
      if(X[i]<=maxx)  k=i;
      }
       min_x=minx; max_x=maxx;
     }

  max_y=min_y=Y[j];
  for(i=j;i<k;i++)
   {max_y=((Y[i]>max_y)?Y[i]:max_y);
    min_y=((Y[i]<min_y)?Y[i]:min_y);
    }
  unit_y=320.0/(max_y-min_y);
  unit_x=500.0/(max_x-min_x);
  i=j-1;
  x0=80+int((X[i]-min_x)*unit_x+0.5);
  y0=400-int((Y[i]-min_y)*unit_y+0.5);
  for(i=j;i<k;i++)
   {
    x=80+int((X[i]-min_x)*unit_x+0.5);
    y=400-int((Y[i]-min_y)*unit_y+0.5);
    line(x0,y0,x,y);
	x0=x;
	y0=y;
	putpixel(x,y,4);
       }

 frame();
  }




void plotD(float X1[], float Y1[], float X2[], float Y2[],float minx,float maxx, int n)
   {
    float tag1,tag2;
    int i, j, k,x,y,x0,y0;
     iniwindow();

  if (minx==maxx)
    {
    max_x=min_x=X[0];
    for(i=1;i<n;i++)
     { 
      if(X1[i]>max_x)  max_x=X1[i]; if(X2[i]>max_x)  max_x=X2[i];
      if(X1[i]<min_x)  min_x=X1[i]; if(X2[i]<min_x)  min_x=X2[i];
      }
    j=1; k=n;
    }

  else
    {j=1;
    for(i=1;i<n;i++)
     {if(X[i]<=minx)  j=i;
      if(X[i]<=maxx)  k=i;
      }
     min_x=minx; max_x=maxx;
     }

   max_y=min_y=Y1[j];
   for(i=j;i<k;i++)
     {if(Y1[i]>max_y)  max_y=Y1[i]; if(Y2[i]>max_y)  max_y=Y2[i];
      if(Y1[i]<min_y)  min_y=Y1[i]; if(Y2[i]<min_y)  min_y=Y2[i];
      }
   unit_y=320.0/(max_y-min_y);
   unit_x=500.0/(max_x-min_x);

	setcolor(14);
	i=j-1;
	x0=80+int((X1[i]-min_x)*unit_x+0.5);
	y0=400-int((Y1[i]-min_y)*unit_y+0.5);
   for(i=j;i<k;i++)
     {
	x=80+int((X1[i]-min_x)*unit_x+0.5);
	y=400-int((Y1[i]-min_y)*unit_y+0.5);
	line(x0,y0,x,y);
	x0=x;
	y0=y;
	putpixel(x,y,4);
       }

	setcolor(15);
	i=j-1;
	x0=80+int((X2[i]-min_x)*unit_x+0.5);
	y0=400-int((Y2[i]-min_y)*unit_y+0.5);
    for(i=j;i<k;i++)
      {
	x=80+int((X2[i]-min_x)*unit_x+0.5);
	y=400-int((Y2[i]-min_y)*unit_y+0.5);
	 line(x0,y0,x,y);
	 x0=x;
	 y0=y;
	putpixel(x,y,4);
       }

 frame();
 }

/*///End of PLOT.cpp//*/