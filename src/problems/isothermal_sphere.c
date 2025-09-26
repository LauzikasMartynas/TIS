#include "../globals.h"

#define m_p_cgs 1.6726e-24
#define k_B_cgs 1.3807e-16

// Adjust all these to your liking
void setup_isothermal_sphere ( )
{
    
    Problem.Periodic[0] = Problem.Periodic[1] = Problem.Periodic[2] = true;

    Problem.Boxsize[0] = 20.0;
    Problem.Boxsize[1] = 20.0;
    Problem.Boxsize[2] = 10.0;
    
    sprintf ( Problem.Name, "isothermal_Sphere" );
    
    Problem.Rho_Max = 1e2 * m_p_cgs / UDensity; //404 * m_p_cgs / UDensity;

    Density_Func_Ptr = &isothermal_model;
    
}

float isothermal_model ( const int ipart , const double bias )
{
    P[ipart].CMZ_type = 2.0;

    double x = P[ipart].Pos[0] - 0.5 * Problem.Boxsize[0];
    double y = P[ipart].Pos[1] - 0.5 * Problem.Boxsize[1];
    double z = P[ipart].Pos[2] - 0.5 * Problem.Boxsize[2];

    double r = sqrt(x*x+y*y+z*z);

    double rho0 = 1e2 * m_p_cgs / UDensity;
    double r_0 = 1;
    double R_over_r_0 = r / r_0;
    double rho = rho0 / p2(R_over_r_0);
    if (r < 0.1)
        rho = rho0;
    if (fabs(z) < 1)
        return rho0;
    return rho0;
}

void isothermal_Velocity ( const int ipart, float out[3] )
{
    
    //unused
    double x = P[ipart].Pos[0] - 0.5 * Problem.Boxsize[0];
    double y = P[ipart].Pos[1] - 0.5 * Problem.Boxsize[1];
    double z = P[ipart].Pos[2] - 0.5 * Problem.Boxsize[2];

    double r = sqrt(x*x + y*y + z*z);

    double rho_0 = 1.0e2 * m_p_cgs / UDensity;
    double r_0 = 1;
    double g = 6.67259e-8 * UDensity * p2(ULength/UVel); // G in internal units
    double vel_gas_abs = sqrt(4.0 * pi * g * rho_0 * p2(r_0));
    out[0] = - y / r * vel_gas_abs;
    out[1] = x / r * vel_gas_abs;
    out[2] = 0.0;
}

