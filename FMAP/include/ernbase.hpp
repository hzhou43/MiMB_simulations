#ifndef ernbase_hpp__
#define ernbase_hpp__ 1

///Base class for Ern.
class ErnBase {
  public:
    ///Constructor
    ErnBase(double c):low(0.0),cutoff(c){}
    ErnBase(double l, double c):low(l),cutoff(c){}
    virtual ~ErnBase(){}
    
    ///Return ern with distance d.
    double ern(double d){
      return ern2(d*d);
    }
    
    ///Return ern with distance d^2.
    virtual double ern2(double d2)=0;

    ///Return virial presure.
    double vir(double d){
      return vir2(d*d);
    }

    ///Return virial presure with distance d^2.
    virtual double vir2(double d2)=0;

  protected:
    ///Low limit, flatten when below.
    double low;
    ///Cutoff, zero when beyond.
    double cutoff;
    ///Value at low limit.
};
#endif //ernbase_hpp__
