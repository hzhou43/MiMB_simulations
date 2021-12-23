#ifndef atomsbase_hpp__
#define atomsbase_hpp__ 1

#include <cstdio>
#include <valarray>

///Base class for Atoms.
class AtomsBase {
  public:
    AtomsBase(int c, int d):cap(c),dim(d){
      cords.resize(cap*dim);
      dump.resize(dim);
    }
    virtual ~AtomsBase(){}
    
    ///Return distance between i and j
    virtual double distance(int i, int j);

    ///Return distance^2, should be faster than distance
    virtual double distance2(int i, int j);

    ///Return distance to i 
    virtual double distance(int i, const std::valarray<double> &probe);

    ///Return distance^2 to i
    virtual double distance2(int i, const std::valarray<double> &probe);
    
    ///Return PBC distance between i and j
    virtual double distance(int i, int j, double boxlen);

    ///Return PBC distance^2, should be faster than distance
    virtual double distance2(int i, int j, double boxlen);

    ///Return PBC distance to i 
    virtual double distance(int i, const std::valarray<double> &probe, double boxlen);

    ///Return PBC distance^2 to i
    virtual double distance2(int i, const std::valarray<double> &probe, double boxlen);

    ///print out cordinations.
    virtual void showcords(bool txt);

    ///print out cordinations to file object
    virtual void showcords(FILE *fp, bool txt);

    ///save i to dump. 
    virtual void save(int i);

    ///retrieve i from dump.
    virtual void retrieve(int i);

    ///set i from x.
    virtual void set(int i, const std::valarray<double> &x);

    ///dispalce i by x.
    virtual void displace(int i, const std::valarray<double> &x);
    
    ///dispalce i by x.
    virtual void displace(int i, double x[]);

    ///dispalce i by x and apply PBC.
    virtual void displace(int i, const std::valarray<double> &x, double boxlen);

    ///dispalce i by x and apply PBC.
    virtual void displace(int i, double x[],double boxlen);

  protected:
    ///dimesion of Cordinators.
    int dim;
    ///Cap of container.
    int cap;
    ///Cordinators.
    std::valarray<double> cords;
    ///Dump site. 
    std::valarray<double> dump;

  private:
    double pbcdst(double d, double boxlen, double hl);
};
#endif //atomsbase_hpp__
