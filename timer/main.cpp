#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>

#include "Timer.h"

using namespace std;

int n_particle = 1;
int n_timestep = 10000;
double dt = 1.0e-5;

void empty_function()
{


}

int main(int argc, char* argv[])
{
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

    //time_clock = clock();
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
            x[i] = x[i] + vx[i] * dt;
        }
        //time2 = time2 + clock() - time_temp;
    }
    //time_clock = clock() - time_clock;
    timer.update();
    cout<<setw(15)<<"total time: " << "\t" <<time_clock<<endl;
    //cout<<setw(15)<<"time1: " << "\t" <<time1<<endl;
    //cout<<setw(15)<<"time2: " << "\t" <<time2<<endl;
    timer.print_clock();

}
