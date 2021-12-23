#ifndef ljs_hpp__
#define ljs_hpp__ 1

#ifdef __cplusplus
extern "C" {
#endif

double LjLrc(double rc,double rho);
double LjesLrc(double eps, double sig, double rc,double rho);
double LJ(double r2);
double LJes(double eps, double sig, double r2);
double LJForce(double r);
double LJesForce(double eps, double sig, double r);
double LjLfs(double r2, double rc2, double rc, double vRc, double ljFRc);
double LjesLfs(double eps, double sig, double r2, double rc2, double rc, double vRc, double ljFRc);
double LjLfsSlow(double r2, double rc2);
double LjesLfsSlow(double eps, double sig, double r2, double rc2);
double LjRc(double r2, double rc2);
double LjesRc(double eps, double sig,double r2, double rc2);
double LjesWRc(double eps, double sig,double r2, double rc2);
double LjWRc(double r2, double rc2);
double LjesLfsWRc(double eps, double sig,double r2, double rc2);
double LjLfsWRc(double r2, double rc2);
double LjShfSlow(double r2, double rc2);
double LjShfWRc(double r2, double rc2);
double LjesShfSlow(double eps, double sig, double r2, double rc2);

#ifdef __cplusplus
}
#endif

#endif
