#ifndef FUNCTIONOBJECT__H
#define FUNCTIONOBJECT__H
class FunctionObject {
  public:
    FunctionObject(){};
    virtual ~FunctionObject(){};
    virtual double operator() (double *x, double *p) {return 0;};
    virtual Bool_t isValid() { return kFALSE; };
    virtual void Randomize() { printf("Warning! Randomize() has not been implemented in your derived class. Either implement it, or don't call it!\n");};
    virtual void Restore() { printf("Warning! Restore() has not been implemented in your derived class. Either implement it, or don't call it!\n");};
};
#endif
