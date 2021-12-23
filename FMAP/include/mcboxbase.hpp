#ifndef mcboxbase_hpp__
#define mcboxbase_hpp__ 1

#include "prng.hpp"
///Base class for MC box.
class MCBoxBase {
  public:
    ///Constructor
    MCBoxBase(int n, double l, double t, bool r, PrngBase *p):num(n), length(l), temperature(t), reduced(r),prng(p){}
    virtual ~MCBoxBase(){};
    
    ///Return volume.
    inline double vol(){
      return length*length*length;
    }

    ///Return beta from reduce temperature or in K. 
    inline double beta(){
      if (reduced){
	return 1.0/temperature;
      }else{
	//const double N_A = 6.022045000e+23;
        //SI units, 2010 CODATA value, J/K = m2·kg/(s2·K) in SI base units
        //const double k_B = 1.380648813e-23;
        //cal/K         1 steam table calorie = 4.1868 J
        //const double k_B = 3.297623030e-24;
        //return N_A*k_B*Temp/1000.0;
	return 0.001985843*temperature;
      }
    }

    ///Return number density.
    inline double numdensity(){
      return num/this->vol();
    }
  
    ///move n steps
    virtual void steps(int n);
    
    ///Move one step.
    virtual void step();

    ///Report status or set status parameter.
    virtual void stat(bool verb=true);

    virtual double acceptratio(bool init=false);

  protected:
    ///Number of particle.
    int num;
    ///Box Length.
    double length;
    ///reduced unit or not?
    bool reduced;
    ///Temperature.
    double temperature;
    ///Poniter to random number generator.
    PrngBase *prng;
    ///beta 1/kBT.
    double mbeta;
    /// accepted and total steps.
    long nacc,ntot;
};
#endif //mcboxbase_hpp__
