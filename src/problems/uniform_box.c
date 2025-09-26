#include "../globals.h"

#define m_p_cgs 1.6726e-24
#define k_B_cgs 1.3807e-16

// Adjust all these to your liking
void setup_uniform_box ( )
{
    
    Problem.Periodic[0] = Problem.Periodic[1] = Problem.Periodic[2] = true;

    Problem.Boxsize[0] = 1.0;
    Problem.Boxsize[1] = 1.0;
    Problem.Boxsize[2] = 1.0;

    sprintf ( Problem.Name, "Uniform_box" );
    
    Problem.Rho_Max = 1e-4; //404 * m_p_cgs / UDensity;

    Density_Func_Ptr = &uniform_Density;
    U_Func_Ptr = &uniform_U;
    Velocity_Func_Ptr = &uniform_Velocity;
}

float uniform_Density ( const int ipart , const double bias )
{
    return 1e-4;
}

float uniform_U ( const int ipart )
{
    float u_to_temp = 0.63 * m_p_cgs / k_B_cgs * (5.0 / 3.0 - 1.0) * UEnergy;
    return 100/u_to_temp;
}

void uniform_Velocity ( const int ipart, float out[3] )
{
    out[0]=0.0;
    out[1]=0.0;
    out[2]=0.0;
}

