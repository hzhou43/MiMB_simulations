/**
 *	@file    ljs.c
 *	@author  Sanbo Qin
 *	@brief   Lennard-Jones potential
 */ 
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/// Lennard-Jones potential Long-range correction per particle in reduced unit.
double LjLrc(double rc,double rho){
	const double pi8_3=8.0*3.14159265358979/3.0;
	double ern=pi8_3*rho*(1.0/pow(rc,9)/3.0-1.0/pow(rc,3));
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f rho:%f\n",__FILE__,__LINE__,__func__,rc,ern,rho);
#endif
	return ern;
}

/**
 *	@brief Lennard-Jones potential Long-range correction per particle.
 *
 *      \f[ U_{LRC} = 2\pi\rho\int_{rc}^{\inf} r^2 V_{LJ}(r) dr \f]
 *      Mol. Phy. 100, 2025 (2002) Eq. 3 
 *      \f[ \Phi_{LRC}=\frac{8}{3}\pi\rho\epsilon\sigma^{3}[\frac{1}{3}(\frac{\sigma}{r_c})^9-(\frac{\sigma}{r_c})^3] \f]
 *      \f[ \rho=N/V \f]
 */
double LjesLrc(double eps, double sig, double rc,double rho){
        const double pi8_3=8.0*3.14159265358979/3.0;
	double sr=sig/rc;
        double ern=pi8_3*rho*eps*pow(sig,3)*(pow(sr,9)/3.0-pow(sr,3));
#ifdef DEBUG
        printf("%s:%d %s\tr2: %f ern:%f rho:%f\n",__FILE__,__LINE__,__func__,rc,ern,rho);
#endif
        return ern;
}

/// Lennard-Jones potential in reduced unit.
double LJ(double r2){
	double r6=r2*r2*r2;
	double r6v=1.0/r6;
	double ern=4.0*r6v*(r6v-1.0);
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
	return ern;
}

/**
 *      @brief Lennard-Jones potential.
 *
 *	\f[ V_{LJ}(r)=4\epsilon[(\frac{\sigma}{r})^{12}-(\frac{\sigma}{r})^6] \f]
 */
double LJes(double eps, double sig, double r2){
	double sr2=sig*sig/r2;
	double sr6=sr2*sr2*sr2;
	double ern=4.0*eps*sr6*(sr6-1.0);
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
	return ern;
}

/// Lennard-Jones force in reduced unit.
double LJForce(double r){
	double r2=r*r;
	double r6=r2*r2*r2;
	double r6v=1.0/r6;
	//double ern=4.0*r6v*(-12.0*r6v+6.0)/r;
	double ern=-48.0*r6v*(r6v-0.5)/r;
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif 
	return ern;
}

/**
 *	@brief Lennard-Jones force.
 *	
 *	\f[ \frac{\delta V_{LJ}(r)}{\delta r}=4\epsilon[-12(\frac{\sigma}{r})^{12}+6(\frac{\sigma}{r})^6]/r \f]
 */
double LJesForce(double eps, double sig, double r){
	double sr=sig/r;
	double sr2=sr*sr;
        double sr6=sr2*sr2*sr2;
	//double ern=4.0*eps*sr6*(-12.0*sr6+6.0)/r;
	double ern=-48.0*eps*sr6*(sr6-0.5)/r;
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r*r,ern);
#endif
	return ern;
}

/// Lennard-Jones potential with Linear-force shift in reduced unit.
double LjLfs(double r2, double rc2, double rc, double vRc, double ljFRc){
	double ern=0.0;
	if (r2<=rc2){
		double r=sqrt(r2);
		ern=LJ(r2)-vRc-ljFRc*(r-rc);
	}
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
	return ern;
}

/// Lennard-Jones potential with Linear-force shift.
double LjesLfs(double eps, double sig, double r2, double rc2, double rc, double vRc, double ljFRc){
	double ern=0.0;
	if (r2<=rc2){
		double r=sqrt(r2);
		ern=LJes(eps,sig,r2)-vRc-ljFRc*(r-rc);
	}
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
	return ern;
}

/// Lennard-Jones potential with Linear-force shift in reduced unit.
double LjLfsSlow(double r2, double rc2){
	double ern=0.0;
	if (r2<=rc2){
		double rc=sqrt(rc2);
		double vRc=LJ(rc2);
		double ljFRc=LJForce(rc);
		double r=sqrt(r2);
		ern=LJ(r2)-vRc-ljFRc*(r-rc);
	}
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
	return ern;
}

/// Lennard-Jones potential shifted in reduced unit.
double LjShfSlow(double r2, double rc2){
        double ern=0.0;
        if (r2<=rc2){
                double vRc=LJ(rc2);
                ern=LJ(r2)-vRc;
        }
#ifdef DEBUG
        printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
        return ern;
}

/**
 *	@brief	Lennard-Jones potential with Linear-force shift.
 *
 * 	For \f$ r \leq r_c \f$
 * 	\f[ V(r)=V_{LJ}(r)-V_{LJ}(r_c)-\frac{\delta V_{LJ}}{\delta r}|_{r_c}(r-rc) \f]
 * 	Slow ~2 folds, since recalculating force and potential at cutoff.
 */
double LjesLfsSlow(double eps, double sig, double r2, double rc2){
	double ern=0.0;
	if (r2<=rc2){
		double rc=sqrt(rc2);
		double vRc=LJes(eps,sig,rc2);
		double ljFRc=LJesForce(eps,sig,rc);
		double r=sqrt(r2);
		ern=LJes(eps,sig,r2)-vRc-ljFRc*(r-rc);
	}
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
	return ern;
}

/// Lennard-Jones potential shifted
double LjesShfSlow(double eps, double sig, double r2, double rc2){
        double ern=0.0;
        if (r2<=rc2){
                double rc=sqrt(rc2);
                double vRc=LJes(eps,sig,rc2);
                ern=LJes(eps,sig,r2)-vRc;
        }
#ifdef DEBUG
        printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
        return ern;
}

///Lennard-Jones potential with cutoff in reduced unit.
double LjRc(double r2, double rc2){
	double ern=0.0;
	if (r2<=rc2){
		ern=LJ(r2);
	}
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
	return ern;
}

///Lennard-Jones potential with cutoff.
double LjesRc(double eps, double sig,double r2, double rc2){
	double ern=0.0;
	if (r2<=rc2){
		ern=LJes(eps,sig,r2);
	}
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
	return ern;
}

/**
 *      @brief  Lennard-Jones instantaneous pair virial.
 *
 *      For \f$ r \leq r_c \f$
 *      \f[ W_{ij} = r_{ij}f_{V_{LJ}} = -r_{ij}\frac{\delta V_{LJ}}{\delta r}|_{r_{ij}}
 *                 = 4\epsilon[12(\frac{\sigma}{r})^{12}-6(\frac{\sigma}{r})^6] \f]
 *      
 */
double LjesWRc(double eps, double sig,double r2, double rc2){
	double ern=0.0;
	if (r2<=rc2){
		double sr2=sig*sig/r2;
       		double sr6=sr2*sr2*sr2;
		//double ern=-4.0*eps*sr6*(-12.0*sr6+6.0)/r*r;
		ern=48.0*eps*sr6*(sr6-0.5);
	}
#ifdef DEBUG
	printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
	return ern;
}


///Lennard-Jones instantaneous pair virial in reduced unit.
double LjWRc(double r2, double rc2){
	double ern=0.0;
	if (r2<=rc2){
        	double r6=r2*r2*r2;
        	double r6v=1.0/r6;
        	//double ern=-4.0*r6v*(-12.0*r6v+6.0)/r*r;
        	ern=48.0*r6v*(r6v-0.5);
	}
#ifdef DEBUG
        printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
        return ern;	
}

/**
 *      @brief  Lennard-Jones potential with Linear-force shift instantaneous pair virial.
 *
 *      For \f$ r \leq r_c \f$
 *      \f[ W_{ij} = {LjesWRc} + {LjesFRc}*r \f]	
 */
double LjesLfsWRc(double eps, double sig,double r2, double rc2){
	double ern=0.0;
	if (r2<=rc2){
		ern=LjesWRc(eps,sig,r2,rc2);
		double rc=sqrt(rc2);
		double r=sqrt(r2);
		double ljFRc=LJesForce(eps,sig,rc);
		ern += ljFRc*r;
	}
#ifdef DEBUG
        printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
        return ern;
}

///Lennard-Jones potential with Linear-force shift instantaneous pair virial in reduced unit
double LjLfsWRc(double r2, double rc2){
        double ern=0.0;
        if (r2<=rc2){
                ern=LjWRc(r2,rc2);
                double rc=sqrt(rc2);
                double r=sqrt(r2);
                double ljFRc=LJForce(rc);
                ern += ljFRc*r;
	}
#ifdef DEBUG
        printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
        return ern;
}

///Lennard-Jones potentiial shifted instantaneous pair virial in reduced unit
double LjShfWRc(double r2, double rc2){
        double ern=0.0;
        if (r2<=rc2){
                ern=LjWRc(r2,rc2);
        }
#ifdef DEBUG
        printf("%s:%d %s\tr2: %f ern:%f\n",__FILE__,__LINE__,__func__,r2,ern);
#endif
        return ern;
}

#ifdef __cplusplus
}
#endif
