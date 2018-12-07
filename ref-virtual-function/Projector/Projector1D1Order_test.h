#ifndef PROJECTOR1D1ORDER_TEST_H
#define PROJECTOR1D1ORDER_TEST_H

class Parameters;

class Projector1D1Order_test {
public:
    Projector1D1Order_test(Parameters& params);
    ~Projector1D1Order_test();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    //! Not used for now
    void operator() (double *rho, double *x, int ipart);


private:
    double dx_inv;
    double dx_ov_dt;
    double xjn, xjmxi;
    int ip;
};

#endif
