/* Conceptual rainfall-runoff model with adaptive-step explicit 
 * Runge-Kutta integrator 
 ---written by JA Vrugt, Feb. 2022
*/

#include "mex.h"                                                                                                                                                 
#include <math.h>
#include <stdlib.h>                                                                                                                                              
#include <stdio.h>       
#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)

static void runge_kutta(int nvar, int nt, double *tout, double *y0, const mxArray *data, const mxArray *options, double *y);
static void rk2(int nvar, int s, double h, double *u, double *LTE, double Precip, double EvapP, double Sumax, double beta, double alfa, double Ks, double Kf, double m);
static int  fRhs(int s, double *u, double *udot, double Precip, double EvapP, double Sumax, double beta, double alfa, double Ks, double Kf, double m);
/*static double expFlux(double Sr, double a); */
/*static double exponen(double x); */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
  double *tout, *y0, *y;
  const mxArray *options, *data;
  mwSize nvar, nt;
  
  /* Create C pointers to input arguments */
  tout    = mxGetPr(prhs[0]);
  y0      = mxGetPr(prhs[1]);
  data    = prhs[2];
  options = prhs[3];
    
  /* Allocate memory for output arguments */
  nvar = mxGetM(prhs[1]);
  nt   = mxGetN(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(nvar,nt,mxREAL);

  /* Create C pointers to output arguments */
  y = mxGetPr(plhs[0]);

  /* Call the C subroutine */
  runge_kutta(nvar,nt,tout,y0,data,options,y);
  
  return;
}

static void runge_kutta(int nvar, int nt, double *tout, double *y0, const mxArray *data, const mxArray *options, double *y)
{
  double hin, hmax_, hmin_, reltol, *abstol, order;
  double h, t, t1, t2, *LTE, *ytmp, wrms, *w;
  int i, s, ns, accept;
  double *P, *Ep, Sumax, beta, alfa, Ks, Kf, m, Precip, EvapP;
  
  /* Initialize */
  ns = nt - 1;
  for (i=1; i<=nvar; i++) y[i-1] = y0[i-1];
  
  /* Memory allocation for temporary vectors */
  LTE    = malloc(nvar*sizeof(double));
  ytmp   = malloc(nvar*sizeof(double));
  w      = malloc(nvar*sizeof(double));
    
  /* Extract integration options */
  hin    = *mxGetPr(mxGetField(options,0,"InitialStep"));
  hmax_   = *mxGetPr(mxGetField(options,0,"MaxStep"));
  hmin_   = *mxGetPr(mxGetField(options,0,"MinStep"));
  reltol = *mxGetPr(mxGetField(options,0,"RelTol"));
  abstol =  mxGetPr(mxGetField(options,0,"AbsTol"));
  order  = *mxGetPr(mxGetField(options,0,"Order"));
  
  P     =  mxGetPr(mxGetField(data,0,"P"));
  Ep    =  mxGetPr(mxGetField(data,0,"Ep"));
  Sumax = *mxGetPr(mxGetField(data,0,"Sumax"));
  beta  = *mxGetPr(mxGetField(data,0,"beta"));
  alfa  = *mxGetPr(mxGetField(data,0,"alfa"));
  Ks    = *mxGetPr(mxGetField(data,0,"Ks"));
  Kf    = *mxGetPr(mxGetField(data,0,"Kf"));
  m     = *mxGetPr(mxGetField(data,0,"m"));  
/*  mexPrintf("%s %5i %5f %5f\n", "nt",nt,order,abstol); */
  
  /* Integrate from tout[0] to tout[end] */
  for (s=1; s<=ns; s++) {
      /* Set start and end times */
      t1 = tout[s-1];
      t2 = tout[s];
      /* Set initial step */
      h = hin;
      h = max(hmin_,min(h,hmax_));
      h = min(h,t2-t1);
      /* Set initial y */
      for (i=1; i<=nvar; i++) y[i-1+nvar*s] = y[i-1+nvar*(s-1)];
      /* Integrate from t1 to t2 */
      t = t1;
      while (t < t2) {
        /* Advance solution by step h */
        for (i=1; i<=nvar; i++) ytmp[i-1] = y[i-1+nvar*s];
        Precip = P[s-1]; EvapP = Ep[s-1];   
        rk2(nvar,s,h,ytmp,LTE,Precip,EvapP,Sumax,beta,alfa,Ks,Kf,m);
        /* Decide whether to accept current step */
    	accept = 0;
        wrms = 0;
        for (i=1; i<=nvar; i++) {
            w[i-1] = 1.0/(reltol*fabs(ytmp[i-1]) + abstol[i-1]);
            wrms = wrms + pow(w[i-1]*LTE[i-1],2);
        }
        wrms = pow(wrms/nvar,0.5);
        if (wrms <= 1) accept = 1;
        /* Next line added by JAV to avoid getting stuck */
        if (h <= hmin_) { 
            wrms = 0.5;
            accept = 1;
        }
    	/* If accepted, update y and t */
        if (accept > 0) {
    		for (i=1; i<=nvar; i++) y[i-1+nvar*s] = ytmp[i-1];
            t = t + h;
    	}
        /* Compute new step */
        h = h*max(0.2,min(5.0,0.9*pow(wrms,-1.0/order)));
        h = max(hmin_,min(h,hmax_));
        h = min(h,t2-t); 
      /*  mexPrintf("%s %6.3f\n", "TimeStep",h);  
        mexPrintf("%s %d\n", "Accept",accept);  
        mexPrintf("%s %5.2f\n", "MaxTStep",hmax_);  */
      }
  }
  
  /* Free memory for temporary vectors */
  free(LTE);
  free(ytmp);
  free(w);

  return;
}

static void rk2(int nvar, int s, double h, double *u, double *LTE, double Precip, double EvapP, double Sumax, double beta, double alfa, double Ks, double Kf, double m)
{
	int i, flag;
    double *udotE, *uE, *udot;
    
	/* Memory allocation for temporary vectors */
    udotE = malloc(nvar*sizeof(double));
    uE    = malloc(nvar*sizeof(double));
    udot  = malloc(nvar*sizeof(double));
  
	/* Euler solution */
	flag = fRhs(s,u,udotE,Precip,EvapP,Sumax,beta,alfa,Ks,Kf,m);
    for (i=1; i<=nvar; i++) uE[i-1] = u[i-1] + h*udotE[i-1];
	/* Heun solution */
	flag = fRhs(s,uE,udot,Precip,EvapP,Sumax,beta,alfa,Ks,Kf,m);
    for (i=1; i<=nvar; i++) u[i-1] = u[i-1] + 0.5*h*(udotE[i-1] + udot[i-1]);    
	/* Compute estimate of LTE */
    for (i=1; i<=nvar; i++) LTE[i-1] = fabs(uE[i-1]-u[i-1]);

	/* Free memory for temporary vectors */
    free(udotE);
    free(uE);
    free(udot);
   
    return;
}

/* Conceptual rainfall-runoff model */
static int fRhs(int s, double *u, double *udot, double Precip, double EvapP, double Sumax, double beta, double alfa, double Ks, double Kf, double m)
{
    double Su, Sf1, Sf2, Sf3, Ss, Su_r;
    double Ea, Perc, qs, qs_out, qf, qf_out1, qf_out2, qf_out3;
    double EvapI, P_e, Ep_e;

    Su = u[0];
    Ss = u[1];
    Sf1 = u[2];
    Sf2 = u[3];
    Sf3 = u[4];

    Su_r = min(Su/Sumax,1);           
   /* Precip = P[s-1]; EvapP = Ep[s-1];  */
    Perc = Precip * (1 - pow(1-Su_r,beta));
    Ea = EvapP * Su_r * (1 + m)/(Su_r + m);  
    udot[0] = Precip - Perc - Ea; 

    qs = (1-alfa) * Perc; 
    qs_out = Ks * Ss;
    udot[1] = qs - qs_out;

    qf = alfa * Perc;
    qf_out1 = Kf * Sf1;         
    qf_out2 = Kf * Sf2;         
    qf_out3 = Kf * Sf3;         
    udot[2] = qf - qf_out1;
    udot[3] = qf_out1 - qf_out2;
    udot[4] = qf_out2 - qf_out3;
    udot[5] = qf_out3 + qs_out;

return(0);
}
