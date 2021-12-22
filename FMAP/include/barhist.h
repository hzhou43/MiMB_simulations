/**
 *    @file    barhist.h
 *    @author  Sanbo Qin
 *    @brief   Histgram for BAR.
 */
#ifndef _BarHist_H_
#define _BarHist_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct BarHistST *BarHist;
BarHist BarHistInit(double drt, int n, double low, double bw);
void BarHistFree(BarHist *b);
void BarHistZero(BarHist b);
void BarHistUpdate(BarHist b, double value);
void BarHistUpdateWght(BarHist b, double value, double w);
double BarHistSum(BarHist b, double beta, double parC, double (*barfwbw)(double, double, double, double), double *nsum);
void BarHistSave(BarHist b, char* fname, bool txt);
BarHist BarHistLoad(char* fname, bool txt);
void BarHistOpt(BarHist forw, BarHist back, double beta, double tol, double boe[6]);
void BarHistOptMore(BarHist forw, BarHist back, double beta, double tol, int n, double boe[]);
void BarHistOptLoad(char* fnforw, char* fnback, bool txt, double beta, double tol, double boe[6]);
void BoeRep(FILE *fp, int n, double res[],double beta);

#ifdef __cplusplus
}
#endif

#endif  /* _BarHist_H_ */
