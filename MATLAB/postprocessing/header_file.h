/* THIS IS THE HEADER FILE WHERE SPECIFIC HEADER FILES
 * PROTOTYPES, STRUCTURES, AND FUNCTIONS IN ARE DEFINED */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>

/*  PROTOTYPE DEFINITIONS */

/*#define max(m, n) ((m) > (n)  ? (m) : (n))*/
/*#define min(m, n) ((m) < (n)  ? (m) : (n))*/

#define MAXLINE 300 /* Maximum number of characters in a line */
#define MAXWORD 50 /* Maximum number of characters in a word */
#define INF 1.e20 /* Maximum function value */
#define mpar 21 /* Maximum number of parameters */
#define mpt 1980 /* Maximum number of points in population*/
                  /* 60 groups x 37 points (2*mpar+1) */
#define mps 50 /* Maximum number of points in a subcomplex */
#define MAX 50 
#define MH 721 
#define MX 1976 

#define MAXDAY 1 /*15000*/  /* Maximum number of days in simulation */
/*#define MAXTSTEP 3650 /*4081 13214 */ /* 4081 */ /* Maximum numberoftimesteps */
#define MAXFLUX  20 /* Maximum number of fluxes in data set */
#define MAXSTATS 15 /*  Maximum number of statistics */
#define MAXRESULTS 1  /* Maximum number of best objective results */
#define MAXFLOWS 5  /* Maximum number of outflows from model */
#define MAXGAGE 1 /* Maximum number of gages */
#define MAXPAR 25  /* Maximum number of parameters */
#define MINVOLUME 0.001
#define AREA 1944.0  /* km^2 */
#define STARTDAY 0
/*#define TIMESPERDAY 4*/

/* STRUCTURE DEFINITIONS */

struct SMA { /* I/O passed in function FLAND1 (THE SMA
	      * accounting operation module */
  double uztwm,uzfwm,uzk,pctim,adimp,riva,zperc,rexp;   /* Parameters */
  double lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,rserv,side;  /* Parameters */
  double pxmlt,pemlt,mpemlt[12];                        /* Parameters */
  double etmin,etshp,etrng;                    /* et curve parameters */
  double uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc;               /* States */
  double dt,pxv,ep,epdist,tlci;                   /* Other Parameters */
  double UH[19],UH_nord,tlci_flows[1];  
};

struct FSUM1 { /* I/O passed in function FLAND1 (THE SMA
		* accounting operation module */
  double srot,simpvt,srodt,srost,sintft,sgwfp,sgwfs;
  double srecht,sett,se1,se3,se4,se5;
};

struct CRITERION {
  char name [10];
  int calflag;
  double current;
};



struct DATE {
     int day;
     int month;
     int year;
};

struct MODELTYPE {
     char name[10];
     int Qperday;
     int Precipperday;
     int PETperday;
};


/* SCE STRUCTURES */



/* FUNCTION DECLARATION */

/*FILE *fopen(); */
/*FILE *in,*out1,*out2,*out3,*out4;*/

/*int sac_sma(double Pars[], double States[], double Obs[], double tlci_out[], double SimStates[], double basefcc[]);*/
                   
/*int fland1(struct SMA *sma,struct FSUM1 *fsum1);*/

/* int print(struct NUM *pnum,struct PARM *pparm,struct OPT_PARM *popt,int m,int n,int ngs1);*/

/* void read_model_pars(struct PARAMETERS *parameters,struct OBS *obs,char filename_parms [MAXLINE] ){ */

/*double func(double *x,struct NUM *pnum, struct PARM *pparm,struct PARAMETERS *parameters,struct OBS *obs,struct MODELTYPE *modeltype,struct OBJECTIVES *objectives,struct STATS *stats,struct OUTPUT *output,struct OPT_PARM *popt);*/

/* END OF FILE */ 
