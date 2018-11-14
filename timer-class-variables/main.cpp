#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>

#include "Timer.h"
#include "Parameters.h"

using namespace std;

int n_particle = 10000;
int n_timestep = 10000;
double dt = 1.0e-5;

void empty_function()
{


}

int main(int argc, char* argv[])
{
    Parameters params;

    double time_clock, time1, time2, time_temp;
    double *x;
    double *vx;

    Timer timer;
    timer.init("time1");
    time1 = 0.0;
    time2 = 0.0;

    x = new double[n_particle];
    vx = new double[n_particle];
    for(int i = 0; i < n_particle; i++)
    {
        x[i] = 12.23124e-2;
        vx[i] = 2424232.34242;
    }  

    timer.restart();
    for(int i_time = 0; i_time < n_timestep; i_time++)
    {
        //time_temp = clock();
        //timer.restart();
        //timer.update();

        //empty_function();
        //empty_function();
        //time1 = time1 + clock() - time_temp;

        //time_temp = clock();
        for(int i = 0; i < n_particle; i++)
        {
            //x[i] = x[i] + vx[i] * dt;
            x[i] = x[i] + vx[i] * params.dt;
        }
        //time2 = time2 + clock() - time_temp;
    }
    timer.update();
    timer.print_clock();

}
