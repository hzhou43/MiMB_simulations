#ifndef mcnvtbase_hpp__
#define mcnvtbase_hpp__ 1

#include <assert.h>
#include <iostream>
#include "mcboxbase.hpp"
#include "atoms3d.hpp"
#include "ernbase.hpp"

///Base class for MC nvt.
class MCNVTBase:public MCBoxBase, public Atoms3D {
  public:
    MCNVTBase(int n, double l, double t, bool r, PrngBase *p, int c, ErnBase *e):MCBoxBase(n,l,t,r,p),Atoms3D(c),ff(e){
      mbeta=beta();
    }
    ~MCNVTBase(){}

    ///set dr, the max displacement.
    void setdr(double x){
      mdr=x;
    }

    ///Move one step.
    void step();
 
    void init_cube(){
      Atoms3D::init_cube(length);
    }
    
    void init_fcc(){
      Atoms3D::init_fcc(length);
    }

    ///Report status.
    virtual void stat(bool verb=true);
    
    ///Return ernergy between i and j.
    virtual double ern(int i, int j);

    ///Return ernergy between i and others.
    virtual double ern(int i);

    ///Return total ernergy.
    virtual double ern();

    ///Return total ernergy, skip recalculation if false.
    virtual double ern(bool redo);

    ///Return insertion ernergy.
    virtual double ern(const std::valarray<double> &probe);

    ///Return widom insertion ernergy.
    virtual double widom();

    ///Return widom insertion ernergy, with Prng.
    virtual double widom(PrngBase *p);

    ///Return widom insertion ernergy, with Prng, and ErnBase.
    virtual double widom(PrngBase *p, ErnBase *f);

    ///Return sum of multiple widom insertion ernergies ernergy in boltzman space.
    virtual double widomsbltz(int n);
    
    ///Return sum of multiple widom insertion ernergies ernergy in boltzman space, with Prng.
    virtual double widomsbltz(PrngBase *p, int n);

    ///Return random deletion ernergy.
    virtual double pickdel();

    ///Return random deletion ernergy,with Prng.
    virtual double pickdel(PrngBase *p);

    ///Return sum of multiple random deletion ernergies in boltzman space.
    virtual double pickdelsbltz(int n);

    ///Return sum of multiple random deletion ernergies in boltzman space,with Prng.
    virtual double pickdelsbltz(PrngBase *p, int n);
    
    ///Return virial presure between i and j.
    virtual double vir(int i, int j);

    ///Return virial presure between i and others.
    virtual double vir(int i);

    ///Return virial presure ernergy.
    virtual double vir();

    ///Return virial presure of insertion
    virtual double vir(const std::valarray<double> &probe);
  
  protected:
    //Pointer to ernergy object.
    ErnBase *ff;
    double mern,mvir;
    double mbeta;
    double mdr;
};
#endif //mcnvtbase_hpp__
