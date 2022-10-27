/// EXF.cpp //

void savedata()
  { FILE *fout;
    char name[30];
    char *string1,*ptr;
    int i,wmax;
    char c,c1, ch[2];

     ch[1]='\0';
     setfillstyle(1,0);
     bar(6,45,620,474);
     setcolor(15);

  settextjustify(LEFT_TEXT, CENTER_TEXT);
  outtextxy(100,100,"Choose data set to Save [0-4]:");
  outtextxy(130,140,"[0]: Normalized Absorption data");
  outtextxy(130,160,"[1]: EXAFS function data Chi(k)");
  outtextxy(130,180,"[2]: k-Weighted EXAFS and FT Ampl." );
  outtextxy(130,200,"[3]: FT Ampl. + REAL + IMAG. parts");
  outtextxy(130,220,"[4]: Filtered EXAFS data in k-space");
  outtextxy(130,240,"[5]: EXIT");
  do{c=getch();} while(c<48||c>53);
    sprintf(name,"%c",c);
    outtextxy(500,260,name);  

  if(c>47&&c<53)
    {
     outtextxy(100,300,"Input file Name:");
     setcolor(4);
     i=0;
     ch[0]=getch();
     while(ch[0]!=13||i==0)
       {if(ch[0]!=8&&ch[0]>45&&ch[0]<123)
	  { outtextxy(250+i*8,300,ch);
	    name[i]=ch[0];
	    i=i+1;
	  }
	else
	  { if(i!=0)
	    {i=i-1;
	    bar(496+i*8,236,504+i*8,244);
	     }
	   } 
	  ch[0]=getch();
	 }
    name[i]='\0';
    c1 = '.';
    ptr = strchr(name, c1);
    if (!ptr) strcat(name,string1);
		fout=fopen(name,"w+");
    }
	
  switch(c){
    case '0':  
        for(i=0;i<num;i++)
	fprintf(fout,"%-13.8f   %-13.8f\n",EX0[i],EY[i]/step);
	break;
    case '1': 
        for(i=0;i<wmax;i++)
	fprintf(fout,"%-13.8f   %-13.8f\n",KX[i],KY[i]);
	break;
    case '2':  fprintf(fout,"FT range from %-13.8f to  %-13.8f\n",ev,es); 
        for(i=0;i<wmax;i++)
	fprintf(fout,"%-13.8f   %-13.8f   %-13.8f   %-13.8f\n",KX[i],KYW[i], RX[i],RY[i]);
	break;
    case '3':  fprintf(fout,"R   FT   Real   Imagin\n");
        for(i=0;i<wmax;i++)
        fprintf(fout,"%-13.8f   %-13.8f   %-13.8f   %-13.8f\n",RX[i],RY[i], RYr[i], RYi[i]);
        break;
    case '4':  fprintf(fout,"k   Real   Amplitude  Phase\n");
        for(i=0;i<wmax;i++)
        fprintf(fout,"%-13.8f   %-13.8f   %-13.8f   %-13.8f\n",BX[i],BY[i], Amp[i], Pha[i]);
    case '5':     break;
	}
   fclose(fout);
   outtextxy(380, 460,"press 'Esc' to exit");
   exitwindow();
  settextjustify(CENTER_TEXT, CENTER_TEXT);
  }



void rd_fl(int proces)
{ FILE *file;
  int i,handle;
  char name[20],*string1,*path,*buf1,*buf2,*buf3,*ptr;
  char  c,ch[2];
  char  extention[5][8]={"ra.dat","bk.dat","pk.dat","xk.dat","sp.dat"};
  setcolor(15);
  ch[1]='\0';
       buf1=new char(20);
       buf2=new char(40);
       buf3=new char(40);
	switch(proces){
	 case 1:  string1=extention[0];break;
	 case 2:  string1=extention[1];break;
	 case 3:  string1=extention[2];break;
	 case 4:  string1=extention[3];break;
	 case 5:  string1=extention[4];break;
	 default: break;
	 }
     do{
       setfillstyle(1,0);
       bar(10,60,620,473);
       outtextxy(100,240,"Input the File Name:");
       outtextxy(194,240,"(*");
       outtextxy(230,240,string1);
       outtextxy(260,240,")");
       i=0;
       ch[0]=getch();
       while(ch[0]!=13||i==0)
	{if(ch[0]!=8)
	  { outtextxy(300+i*8,240,ch);
	    name[i]=ch[0];
	    i=i+1;
	  }
	 else
	   { if(i!=0)
	     {i=i-1;
	    bar(496+i*8,236,504+i*8,244);
	     }
	   }
	  ch[0]=getch();

	 }
	name[i]='\0';
      c = '.';
      ptr = strchr(name, c);
      if (!ptr)
      strcat(name,string1);
    stpcpy(filename, name);
      path=searchpath(name);
    }while(((handle=open(path,O_RDONLY))==-1)||(file=fdopen(handle,"r+"))==NULL);
    for(i=0;i<MAXNM;i++)
    {Y[i]=0.0;
     X[i]=0.0;
    }
   fscanf(file,"%s", name);
   num=atoi(name);
   ntag[0]=num;
   for(i=0;i<num;i++)
   {fscanf(file,"%s",buf1);
   fscanf(file,"%s",buf2);
   fscanf(file,"%s",buf3);
      X[i]=atof(buf2);
      Y[i]=atof(buf3);
   }
   if(proces>=2)
    {fscanf(file,"%s",buf1);
     fscanf(file,"%s",buf2);
     ef=atof(buf2);
    }
     fclose(file);
  }

/*///End of EXF.cpp////*/