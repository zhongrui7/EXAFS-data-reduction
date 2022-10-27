// Main program for EXAFS analysis: XAFS0.cpp //

#include "wan.h"
#include "plot.cpp"
#include "preedg.cpp"
#include "postedg.cpp"
#include "FFT.cpp"
#include "IMEX.cpp"

void routine(int proces)
  {char ch;
   setfillstyle(SOLID_FILL,1);
   switch(proces)
    {case 1:    rd_fl(proces);
	 	getE0();
		break;
    case 2:     if(ntag[0]>1){select_es();
		exitwindow();
		sub();}
		break;
    case 3:     if(ntag[1]>1){mu0(numb);     //* post-edge mu0 fitting *//
 		E2k(numb,MAXNM);}
		break;
    case 4:     if(ntag[2]>1){transform(MAXNM);}  //* Fourier transform *//
                break;
    case 5:     if(ntag[3]>1){inverse(MAXNM);}  //* inverse(2*MAXNM); Inversed FT *// 
	        break;
    case 7:     if(ntag[2]>1){savedata();}
		break;
    }
  }


void main(void)
{
//  int gdriver = VGA;
//   GraphDriver = DETECT;
//   initgraph(&GraphDriver, &GraphMode, "" );
//   int gmode=VGAHI;
//  registerbgidriver(EGAVGA_driver);
   int errorcode, gdriver = DETECT, gmode;
   initgraph(&gdriver, &gmode, " ");
   if ((errorcode = graphresult() ) != grOk)
     {
      printf ("Error : %s\n", grapherrormsg(errorcode));
      exit (1);
      };
   setbkcolor(0);
   clearviewport();
   rectangle(0,0,639,479);
   void *image;
   setfillstyle(SOLID_FILL,15);
   bar(20,10,90,40);
   int imsize=imagesize(20,10,90,40);
   image=malloc(imsize);
   getimage(20,10,90,40,image);
   cleardevice();
   settextstyle(DEFAULT_FONT,HORIZ_DIR,1);
   settextjustify(CENTER_TEXT,CENTER_TEXT);
   setfillstyle(SOLID_FILL,6);
   bar(20,10,580,40);
   setfillstyle(EMPTY_FILL,7);
   rectangle(20,10,90,40);
   rectangle(90,10,160,40);
   rectangle(160,10,230,40);
   rectangle(230,10,300,40);
   rectangle(300,10,370,40);
   rectangle(370,10,440,40);
   rectangle(440,10,510,40);
   rectangle(510,10,580,40);

   setcolor(10);

   outtextxy(50,25,"OPEN ");
   outtextxy(125,25,"RmBg");
   outtextxy(200,25,"Spline");
   outtextxy(260,25,"FFT ");
   outtextxy(330,25,"IFT");
   outtextxy(400,25,"Util.");
   outtextxy(470,25,"SAVE");
   outtextxy(540,25,"EXIT");
	outtextxy(300,200,"Welcome to use XAFSpak v0.11e, NSRL/USTC");
	outtextxy(320,260,"Last update 2018-08-26");
   putimage(20,10,image,1);
   int flag=1, post=20;
   char ch;
   ch=getch();
   while(1){
	    if((ch==75)&&(flag!=1))
	     {
		exitwindow();
		putimage(post,10,image,1);
	        post=post-70;
	        putimage(post,10,image,1);
		flag=flag-1;
	      }
	     if((ch==77)&&(flag!=8))
	     {
	        putimage(post,10,image,1);
	        post=post+70;
	        putimage(post,10,image,1);
	        flag=flag+1;
	      }
	      if(ch==13){
		 switch(flag)
		     {
			       case 8:closegraph();clrscr();exit(0);
			       case 1:
			       case 2:
			       case 3:
			       case 4:
			       case 5:
			       case 6:
			       case 7:  routine(flag);break;
		       }
			}
	    ch=getch();
	    }

      }
/*//End of XAFS0.cpp//*/