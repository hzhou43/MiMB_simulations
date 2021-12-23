#ifndef prng_hpp__
#define prng_hpp__ 1


///pure abstract class psudeorandom number generator.
class PrngBase {
  public:
    ///Constructor.
    PrngBase(int s): seed(s){}
    virtual ~PrngBase(){};
   
    ///Return uniform random number between [0,1).
    virtual double ran()=0;
 
    ///Return uniform random number integer between [0,up)
    virtual int ranI(int up);   
 
    ///check ran
    double check(int n);
  
  protected:
    ///seed for random number genrator.
    int seed;
    ///pointer for random number genrator state.
    void *prng=NULL;
    /// Initilize random number genrator.
    virtual void init()=0;
};

/// Prng use srand from system, not safe with multiple objects.
class PrngSys:public PrngBase{
  public:
    PrngSys(int s):PrngBase(s){init();}
    ~PrngSys(){};

    ///Return uniform random number between [0,1).
    double ran();
    
  protected:
    void init(); 
};

/// Prng from gsl, safe with multiple objects, should link with -lgsl.
class PrngGSL:public PrngBase{
  public:
    PrngGSL(int s):PrngBase(s){init();}
    ~PrngGSL();

    ///Return uniform random number between [0,1).
    double ran();
    
    ///Return uniform random number integer between [0,up)
    //int ranI(int up);
  
  protected:
    void init();  
};

/// Prng from Numerical Recipe, safe with multiple objects.
class PrngNR:public PrngBase{
  public:
    PrngNR(int s):PrngBase(s){init();}
    ~PrngNR();

    ///Return uniform random number between [0,1).
    double ran();
     
  protected:
    void init();
}; 

#endif //prng_hpp__
