/*
 * lineread.c - output adequate line information for Transit.
 *              Part of Transit program.
 *
 * Copyright (C) 2003 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

/* TD: replace gabby_read in dbread_pands and lineread.h by verblevel */

#include <transit.h>
#include <lineread.h>
/*#include "../transit/transitstd.c"*/
#include <math.h>

/* Version history:
   0.5:  First light of procedure. Oct 7 2003. Patricio Rojo
   0.8:  Added 'dbid' and its effect on '.isoid' field. Added a starting
         value of 0 for 'dindex' marks in infofile. Added transitprint().
 	 100703. PMR
   0.9:  Changed 'dbid' function, now we store a correlative number in
         '.isoid' that starts at 0 in the first isotope of the first
         database and increase from there. 100903. PMR
   0.11: added option '-n' for no file output. 102603. PMR
   0.13: Verbose improved. 110903. PMR
   1.2:  Isotope information in the Infofile is now given in a
         sequential order and not per database as before. 012504. PMR
   1.3:  Version of TWII is now stored in the .inf file. 021204. PMR
   1.4:  Bug found and corrected that wrote a lot of masses instead of
         1 per isotope. 041203. PMR.
   1.5:  Change to directory structure and linking (several files
         instead of one). 040904. PMR.
*/
static int lineread_ver=1;
static int lineread_rev=5;



short gabby_dbread=0;
/* Add sources and name for every new reader driver */
#define dbread_nfcn 1
PREC_NREC (*linefcn[dbread_nfcn])(char *,struct linedb **,float,float
				    ,char *, PREC_ZREC ***,PREC_ZREC **
				    ,PREC_ZREC **, int *, int *, char ***)={
				      dbread_pands
				    };
char *dname[dbread_nfcn]={
  "Partridge & Schwenke. Water"
};


/*
  synhelp: Help on syntax
*/
static void synhelp(float deltw,char *datafile,char *infofile,
		    float deltwmark, int verblevel)
{
  fprintf(stderr,
	  "Syntax:\n\tlineread [..options..] <wavl_i> <wavl_f>\n\n"
	  " Where available options and their defaults are:\n"
	  "  -d<delt_wavl>   : range of wavelength in nanometers to be read\n"
	  "                    and write at once from each database (%.2f).\n"
	  "  -m<mark_wavl>   : Generate an index in the info file every this\n"
	  "                    much nanometers (%.2f)\n\n"
	  "  -n              : Dummy run!, do not write anything to files\n"
	  "  -o              : output data file (\"%s\").\n"
	  "  -O              : output info file (\"%s\").\n"
	  "        For any of the previous two options, you can specify"
	  " a dash '-' to\n"
	  "     indicate output to be sent to the standard output."
	  "  Note however, that\n"
	  "     the output is in binary format  and that at least the"
	  " output of one of\n"
	  "     files should be directed to file to be compatible with"
	  " Transit.\n\n"
	  "  -v and -q       : increase or decrease verbose level(%i).\n"
	  "  -h              : show this help.\n\n"
	  ,deltw,deltwmark,datafile,infofile,verblevel);
  exit(EXIT_FAILURE);
}

int main(int argc,char *argv[])
{
  int i,j,rc;
  int left; //Number of databases with lines available for reorder.
  float deltw;
  double iniw,finw,parw;
  struct linedb **lineread,**crnt;
#ifndef NODEBUG
  struct line_transition *ltest;
#endif
  PREC_ZREC ***Z,**T,**mass;
  PREC_LNDATA tmin;
  int *nT,*nIso,*dbid,totaliso,adb;
  PREC_NREC *nlines,pmin;
  FILE *fpdout,*fpiout;
  char *datafile,*infofile;
  PREC_NREC dindex;
  void *tempp;
  int dummy;
  int verblevel;
  float prevmark,deltwmark;    //Mark every that much nanometers in info
			       //file
  char ***isonames;


  if(dbread_nfcn<1)
    transiterror(TERR_CRITICAL,
		 "No drivers for reading database selected or found!!");


  crnt=    (struct linedb **)calloc(dbread_nfcn,sizeof(struct linedb *));
  lineread=(struct linedb **)calloc(dbread_nfcn,sizeof(struct linedb *));
  Z=       (PREC_ZREC ***)   calloc(dbread_nfcn,sizeof(PREC_ZREC **)   );
  T=       (PREC_ZREC **)    calloc(dbread_nfcn,sizeof(PREC_ZREC *)    );
  mass=    (PREC_ZREC **)    calloc(dbread_nfcn,sizeof(PREC_ZREC *)    );
  nlines=  (PREC_NREC *)     calloc(dbread_nfcn,sizeof(PREC_NREC)      );
  isonames=(char ***)        calloc(dbread_nfcn,sizeof(char **)        );
  dbid=    (int *)           calloc(dbread_nfcn,sizeof(int)            );
  nIso=    (int *)           calloc(dbread_nfcn,sizeof(int)            );
  nT=      (int *)           calloc(dbread_nfcn,sizeof(int)            );

  deltw=40;
  verblevel=1;
  datafile=(char *)calloc(19,sizeof(char));
  infofile=(char *)calloc(19,sizeof(char));
  strcpy(datafile,"./res/lineread.dat");
  strcpy(infofile,"./res/lineread.inf");
  deltwmark=1;
  fpdout=fpiout=NULL;
  dummy=0;

  if(argc<3)
    synhelp(deltw,datafile,infofile,deltwmark,verblevel);

  while(1){
    rc=getopt(argc,argv,"vhqd:o:O:m:n");
    if (rc==-1)
      break;

    switch(rc){
    case 'n':
      dummy=1;
      break;
    case 'h':
      synhelp(deltw,datafile,infofile,deltwmark,verblevel);
      break;
    case 'm':
      deltwmark=atof(optarg);
      break;
    case 'd':
      deltw=atof(optarg);
      break;
    case 'q':
      verblevel--;
      break;
    case 'v':
      verblevel++;
      break;
    case 'o':
      datafile=(char *)realloc(datafile,strlen(optarg)*sizeof(char));
      strcpy(datafile,optarg);
      break;
    case 'O':
      infofile=(char *)realloc(infofile,strlen(optarg)*sizeof(char));
      strcpy(infofile,optarg);
      break;
    default:
      synhelp(deltw,datafile,infofile,deltwmark,verblevel);
      break;
    }
  }
  iniw=atof(argv[optind++]);
  finw=atof(argv[optind]);
  gabby_dbread=verblevel;

  transitprint(1,verblevel,
	    "           LINEREAD v%i.%i. Part of Transit package\n"
	    "--------------------------------------------------------------\n",
	    lineread_ver,lineread_rev);
  transitprint(1,verblevel,
	       "Reading %i line database(s)\n\n",dbread_nfcn);

  rc=0;
  if(strcmp("-",datafile)==0){
    rc=1;
    fpdout=stdout;
  }
  if(strcmp("-",infofile)==0){
    rc+=2;
    fpiout=stdout;
    if(fpdout==stdout)
      transiterror(TERR_WARNING,
		   "Both information and data will go to the standard\n"
		   "output. Be warned, that this won't be compatible\n"
		   "with main Transit program\n");
  }

  transitprint(2,verblevel,
	       "Extra blah, blah enabled. Be prepared!\n\n"
	       "Total wavelength range is %.2f to %.2f Nanometers.\n"
	       "The databases are going to be read and written in the"
	       " standard TWII (Transit\n"
	       "wavelength information  interchange)  format in smaller"
	       " wavelength ranges of\n"
	       "%.1f nanometers each.\n\n"
	       "Indexing for speed-up of reading will be done every %.1f"
	       " nanometers\n\n"
	       ,iniw,finw,deltw,deltwmark);

  if(dummy)
    transitprint(1,verblevel,
		 "Dummy run: No output. However it simulates\n");

  transitprint(1,verblevel,
	       "TWIIf requires two output files that were chosen to be:\n"
	       " '%s' as the info file, and\n"
	       " '%s' as the data file\n\n"
	       ,rc&2?"standard output":infofile
	       ,rc&1?"standard output":datafile);


  if(!dummy){
    if(fpdout==NULL&&(fpdout=fopen(datafile,"w"))==NULL){
      transiterror(TERR_SERIOUS,
		   "Data file '%s' cannot be opened for writing.\n"
		   ,datafile);
    }

    if(fpiout==NULL&&(fpiout=fopen(infofile,"w"))==NULL){
      transiterror(TERR_SERIOUS,
		   "Information file '%s' cannot be opened for writing.\n"
		   ,infofile);
    }
  }

  if(finw<iniw)
    transiterror(TERR_SERIOUS,
		 "Final wavelength (%.2f) has to be greater than\n"
		 "initial wavelength (%.2f)\n",finw,iniw);

  if(!dummy){
    fwrite(&lineread_ver, sizeof(int),1,fpiout);
    fwrite(&lineread_rev, sizeof(int),1,fpiout);
    rc=strlen(datafile);
    fwrite(&rc,sizeof(int),1,fpiout);
    fwrite(datafile,sizeof(char),rc,fpiout);
    fwrite(&deltwmark,sizeof(float),1,fpiout);
    fwrite(&iniw,sizeof(double),1,fpiout);
    fwrite(&finw,sizeof(double),1,fpiout);
    rc=dbread_nfcn;
    fwrite(&rc,sizeof(int),1,fpiout);
  }

  dindex=0;
  prevmark=iniw;
  while(iniw<finw){
    parw=iniw+deltw;
    if(parw>finw)
      parw=finw;

    transitprint(1,verblevel,
		 "*******    Wavelength range %8.2f - %8.2f    *******\n"
		 ,iniw,parw);

    left=0;
    /* Reading of data in the range [iniw,parw] */
    for (i=0;i<dbread_nfcn;i++){
      transitprint(1,verblevel,"    Database %i (%s): \n",i+1,dname[i]);
      tempp=dindex?NULL:Z+left;	/* Only read Z if first time */
      if((nlines[left]=(linefcn[i])(NULL,lineread+left,iniw,parw,NULL,
				    tempp,T+left,mass+left,nT+left,
				    nIso+left,isonames+left))>0){
	crnt[left]=lineread[left];
	transitDEBUG(12,verblevel,"wlprev:%f\nwlnew:%f\n"
		     ,lineread[left]->wl,
		     (ltest=(void *)&(crnt[left]->wl))->wl);
	if(left)
	  dbid[left]=dbid[left-1]+nIso[left-1];
	else
	  dbid[left]=0;
	left++;
      }
      else
	transiterror(TERR_WARNING,
		     "Database %i didn't have any line in the wavelength\n"
		     "range %f - %f, or there was an error.\n"
		     ,iniw,parw);

    }
    totaliso=dbid[left-1]+nIso[left-1];

    /* Sorting and output of infofile is done here only while looking at
       the first range */
    if(!dindex&&!dummy){
      for(i=0;i<dbread_nfcn;i++){
	rc=strlen(dname[i]);
	fwrite(&rc,sizeof(int),1,fpiout);
	fwrite(dname[i],sizeof(char),rc,fpiout);
	fwrite(nT+i,sizeof(int),1,fpiout);
	fwrite(nIso+i,sizeof(int),1,fpiout);
	fwrite(T[i],sizeof(PREC_ZREC),nT[i],fpiout);
      }
      fwrite(&totaliso,sizeof(int),1,fpiout);
      for(adb=0;adb<dbread_nfcn;adb++){
	for(j=0;j<nIso[adb];j++){
	  fwrite(&adb,sizeof(int),1,fpiout);
	  fwrite(mass[adb]+j,sizeof(PREC_ZREC),1,fpiout);
	  transitDEBUG(22,verblevel,
		       "Just wrote mass %g at %li\n"
		       ,mass[adb][j],ftell(fpiout));
	  rc=strlen(isonames[adb][j]);
	  fwrite(&rc,sizeof(int),1,fpiout);
	  transitDEBUG(22,verblevel,
		       "Just wrote name's length %i at %li\n"
		       ,rc,ftell(fpiout));
	  fwrite(isonames[adb][j],sizeof(char),rc,fpiout);
	  fwrite(Z[adb][j],sizeof(PREC_ZREC),nT[adb],fpiout);
	}
      }
      fwrite(&dindex,sizeof(PREC_NREC),1,fpiout);
    }

    transitprint(1,verblevel,"sorting... ");

    /* Merge and output of sorted data line list */
    while(left){
      tmin=crnt[0]->wl;
      pmin=0;
      for(i=1;i<left;i++){
	if(crnt[i]->wl<tmin){
	  tmin=crnt[i]->wl;
	  pmin=i;
	}
      }

      /* Instead we are leaving isoid as a correlative list and we are
	 storing its database information in the array 'dbinfo' */
      crnt[pmin]->isoid=dbid[pmin]+crnt[pmin]->isoid;

      if(!dummy)
	fwrite(&(crnt[pmin]->wl),sizeof(struct line_transition),1,fpdout);
      dindex++;
      if(crnt[pmin]->wl - prevmark > deltwmark){
	if(!dummy)
	  fwrite(&dindex,sizeof(PREC_NREC),1,fpiout);
	prevmark=crnt[pmin]->wl;
	transitDEBUG(15,verblevel,"Mark set: %li\n",dindex);
      }

      if(++crnt[pmin]-lineread[pmin]>=nlines[pmin]){
	for(i=pmin;i<left-1;i++){
	  dbid[i]=dbid[i+1];
	  crnt[i]=crnt[i+1];
	}
	left--;
      }
    }
    //    Pprintf(2,"DD:read %f - %f\n",iniw,parw);
    iniw+=deltw;
    for (i=0;i<dbread_nfcn;i++){
      free(lineread[i]);
      free(Z[i]);
      free(T[i]);
      free(mass[i]);
    }
    transitprint(1,verblevel,"done\n");
    //      Pprintf(2,"DD:Freed\n",left);
  }
  rc=0;
  if(!dummy){
    fwrite(&rc,sizeof(int),1,fpiout);

    fclose(fpiout);
    fclose(fpdout);
  }
  transitprint(1,verblevel,"\n");

  return EXIT_SUCCESS;
}


