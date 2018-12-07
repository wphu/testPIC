#ifndef INTERPOLATOR1D1ORDER_TEST_H
#define INTERPOLATOR1D1ORDER_TEST_H

class Parameters;

//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D1Order_test
{

public:
    Interpolator1D1Order_test(Parameters& params);
    ~Interpolator1D1Order_test(){};

    void operator() (double *Ex, double *x, int ipart, double &Ex_local);

private:
    double dx_inv;
    double xjn, xjmxi;
    int ip;

};//END class

#endif
