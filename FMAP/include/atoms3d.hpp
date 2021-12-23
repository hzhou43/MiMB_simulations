#ifndef atoms3d_hpp__
#define atoms3d_hpp__ 1

#include <vector>
#include "atomsbase.hpp"

struct Point{
  double x;
  double y;
  double z;
};

struct PointInt{
  int x;
  int y;
  int z;
};

///3D class for Atoms.
class Atoms3D:public AtomsBase{
  public:
    Atoms3D(int c):AtomsBase(c,3){}
    virtual ~Atoms3D(){}
 
    void init_fcc(double boxlen);
    
    void init_cube(double boxlen);   
    
    void init(int n, double xyz[][3]);

    void init(std::vector<Point> xyz);
};
#endif //atoms3d_hpp__
