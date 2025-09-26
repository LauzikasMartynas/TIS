#include "../globals.h"

#define m_p_cgs 1.6726e-24
#define k_B_cgs 1.3807e-16

// Adjust all these to your liking
void setup_isothermal_torus ( )
{
    Problem.Periodic[0] = Problem.Periodic[1] = Problem.Periodic[2] = false;
    Problem.Boxsize[0] = 1.0;
    Problem.Boxsize[1] = 1.0;
    Problem.Boxsize[2] = 1.0;
    
    sprintf ( Problem.Name, "isothermal_torus" );
    
    Problem.Rho_Max = 1e6 * m_p_cgs / UDensity;

    Density_Func_Ptr = &isothermal_torus_model;
}

float isothermal_torus_model ( const int ipart , const double bias )
{
    double x = P[ipart].Pos[0] - 0.5 * Problem.Boxsize[0];
    double y = P[ipart].Pos[1] - 0.5 * Problem.Boxsize[1];
    double z = P[ipart].Pos[2] - 0.5 * Problem.Boxsize[2];

    double r = sqrt(x*x+y*y+z*z) + 0.01;

    double v = 220 * 1e5 / UVel;
    double c_s = 0.58 * v;
    double v_0 = 0.42 * v;
    double v_c2 = p2(v_0)/p2(c_s);
    double A = 2 + v_c2;
    double B = 1;
    double g = 6.67259e-8 / ULength / p2(UVel) * UMass; // G in internal units

    double theta = atan2(sqrt(x*x+y*y), z);
    double rho = 0;

        if (x != 0 && y != 0)
        {
          rho = 1. * p2(c_s )*A *A *B*pow(tan(theta/2.),A )/(2.*pi*g*r*r*(p2(sin(theta))*pow(1.+B*pow(tan(theta/2.),A ),2.)));

        if (Problem.Rho_Max < rho)
          printf("UDensity:%.5e rho_max:%.5e rho:%.5e\n", UDensity, Problem.Rho_Max, rho * UDensity / m_p_cgs);

        }
         //Clouds

    return rho;
}

void isothermal_torus_Velocity ( const int ipart, float out[3] )
{
    double x = P[ipart].Pos[0] - 0.5 * Problem.Boxsize[0];
    double y = P[ipart].Pos[1] - 0.5 * Problem.Boxsize[1];
    double z = P[ipart].Pos[2] - 0.5 * Problem.Boxsize[2];

    double r = sqrt(x*x + y*y + z*z);

    double vel_gas_abs = 1/3. * 220 * 1e5 / 3.086e21;
    out[0] = - y / r * vel_gas_abs;
    out[1] = x / r * vel_gas_abs;
    out[2] = 0.0;
}

