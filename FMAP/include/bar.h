#ifndef _BAR_H_
#define _BAR_H_

#include<stdlib.h>
#include<stdio.h>
#include<stddef.h>
#include<ctype.h>
#include<float.h>
#include<string.h>
#include<math.h>
#include<assert.h>

#ifdef __cplusplus
extern "C" {
#endif

double BarIns(double mu, double beta, double parC);
double BarDel(double mu, double beta, double parC);
double OsIns(double mu, double beta);
double OsDel(double mu, double beta);
double ExpIns(double mu, double beta);
double ExpDel(double mu, double beta);
double BarForBack(double mu, double beta, double parC, double drt);
double OsForBack(double mu, double beta, double parC, double drt);
double ExpForBack(double mu, double beta, double parC, double drt);

#ifdef __cplusplus
}
#endif

#endif  /* _BAR_H_ */
