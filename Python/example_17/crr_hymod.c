#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)

void runge_kutta(int nvar, int nt, double *tout, double *y0, double *data, double *options, double *y);
void rk2(int nvar, int s, double h, double *u, double *LTE, double Precip, double EvapP, double Sumax, double beta, double alfa, double Ks, double Kf, double m);
int fRhs(int s, double *u, double *udot, double Precip, double EvapP, double Sumax, double beta, double alfa, double Ks, double Kf, double m);

void runge_kutta(int nvar, int nt, double *tout, double *y0, double *data, double *options, double *y) {
    double hin, hmax_, hmin_, reltol, abstol, order;
    double h, t, t1, t2, *LTE, *ytmp, wrms, *w;
    int i, s, ns;
    double *P, *Ep, Sumax, beta, alfa, Ks, Kf, m, Precip, EvapP;
  
    ns = nt - 1;
    for (i = 1; i <= nvar; i++) y[i - 1] = y0[i - 1];
  
    LTE = malloc(nvar * sizeof(double));
    ytmp = malloc(nvar * sizeof(double));
    w = malloc(nvar * sizeof(double));
    
    hin = options[0];
    hmax_ = options[1];
    hmin_ = options[2];
    reltol = options[3];
    abstol = options[4];
    order = options[5];
  
    P = data; 
    Ep = data + nvar;  // Assuming data is structured to store these values
    Sumax = data[6];
    beta = data[7];
    alfa = data[8];
    Ks = data[9];
    Kf = data[10];
    m = data[11];
  
    for (s = 1; s <= ns; s++) {
        t1 = tout[s - 1];
        t2 = tout[s];
        h = hin;
        h = max(hmin_, min(h, hmax_));
        h = min(h, t2 - t1);
      
        for (i = 1; i <= nvar; i++) y[i - 1 + nvar * s] = y[i - 1 + nvar * (s - 1)];
      
        t = t1;
        while (t < t2) {
            for (i = 1; i <= nvar; i++) ytmp[i - 1] = y[i - 1 + nvar * s];
            Precip = P[s - 1];
            EvapP = Ep[s - 1];
            rk2(nvar, s, h, ytmp, LTE, Precip, EvapP, Sumax, beta, alfa, Ks, Kf, m);
            wrms = 0;
            for (i = 1; i <= nvar; i++) {
                w[i - 1] = 1.0 / (reltol * fabs(ytmp[i - 1]) + abstol);
                wrms = wrms + pow(w[i - 1] * LTE[i - 1], 2);
            }
            wrms = pow(wrms / nvar, 0.5);
            if (h <= hmin_) wrms = 0.5;
            if ((wrms <= 1) || (h <= hmin_)) {
                for (i = 1; i <= nvar; i++) y[i - 1 + nvar * s] = ytmp[i - 1];
                t = t + h;
            }
            h = h * max(0.2, min(5.0, 0.9 * pow(wrms, -1.0 / order)));
            h = max(hmin_, min(h, hmax_));
            h = min(h, t2 - t);
        }
    }
  
    free(LTE);
    free(ytmp);
    free(w);
}

void rk2(int nvar, int s, double h, double *u, double *LTE, double Precip, double EvapP, double Sumax, double beta, double alfa, double Ks, double Kf, double m) {
    int i, flag;
    double *udotE, *uE, *udot;
    udotE = malloc(nvar * sizeof(double));
    uE = malloc(nvar * sizeof(double));
    udot = malloc(nvar * sizeof(double));
  
    flag = fRhs(s, u, udotE, Precip, EvapP, Sumax, beta, alfa, Ks, Kf, m);
    for (i = 1; i <= nvar; i++) uE[i - 1] = u[i - 1] + h * udotE[i - 1];
    flag = fRhs(s, uE, udot, Precip, EvapP, Sumax, beta, alfa, Ks, Kf, m);
    for (i = 1; i <= nvar; i++) u[i - 1] = u[i - 1] + 0.5 * h * (udotE[i - 1] + udot[i - 1]);
    for (i = 1; i <= nvar; i++) LTE[i - 1] = fabs(uE[i - 1] - u[i - 1]);

    free(udotE);
    free(uE);
    free(udot);
}

int fRhs(int s, double *u, double *udot, double Precip, double EvapP, double Sumax, double beta, double alfa, double Ks, double Kf, double m) {
    double Su, Sf1, Sf2, Sf3, Ss, Su_r;
    double Ea, Perc, qs, qs_out, qf, qf_out1, qf_out2, qf_out3;
    double EvapI, P_e, Ep_e;

    Su = u[0];
    Ss = u[1];
    Sf1 = u[2];
    Sf2 = u[3];
    Sf3 = u[4];

    Su_r = min(Su / Sumax, 1);
    Perc = Precip * (1 - pow(1 - Su_r, beta));
    Ea = EvapP * Su_r * (1 + m) / (Su_r + m);
    udot[0] = Precip - Perc - Ea;

    qs = (1 - alfa) * Perc;
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

    return 0;
}
