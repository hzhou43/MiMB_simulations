#ifndef timebase_hpp__
#define timebase_hpp__ 1

#include <ctime>
#include <cstdarg>
#include <string>
#include <cstdio>

class TimeBase {
  public:
    ///Constructor
    TimeBase(int n, std::string *jobs);
    ~TimeBase();
    
    void init(); 
    ///Start counting on i
    void start(int i);
    ///Start counting on ...
    void starts(int n, ...);
    
    ///Stop counting on i
    void stop(int i);
    ///Stop counting on ...
    void stops(int n, ...);

    ///Report.
    void report(FILE *fp);
  protected:
    int njob;
    std::string *names;
  private:
    double *mstart;
    double *msum;
#ifdef _OPENMP
    double *cstart;
    double *csum;
#endif
};
#endif //timebase_hpp__
