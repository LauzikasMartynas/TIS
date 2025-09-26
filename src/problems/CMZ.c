#include "../globals.h"

#define m_p_cgs 1.6726e-24
#define k_B_cgs 1.3807e-16
#define G_cgs   6.67259e-8

// Adjust all these to your liking
void setup_CMZ ( )
{
    
    Problem.Periodic[0] = Problem.Periodic[1] = Problem.Periodic[2] = false;

    Problem.Boxsize[0] = 4.0;
    Problem.Boxsize[1] = 4.0;
    Problem.Boxsize[2] = 4.0;

    sprintf ( Problem.Name, "IC_CMZ" );
    
    Problem.Rho_Max = 2e3 * m_p_cgs / UDensity;

    Density_Func_Ptr = &CMZ_Density;
    U_Func_Ptr = &CMZ_U;
    Velocity_Func_Ptr = &CMZ_Velocity;
}

double Beta_model ( float r, const int ipart  )
{
    P[ipart].CMZ_type = 2.0;

    double rho0 = 404.0 * m_p_cgs / UDensity;
    double beta = 1.0;
    double rc = 0.1;

    return rho0 * pow ( 1 + p2 ( r / rc ), -1.5 * beta );
}

//https://ned.ipac.caltech.edu/level5/Sept16/Sofue/Sofue4.html
double Iso_model ( double r , const int ipart )
{
    P[ipart].CMZ_type = 0.0;

    double rho0 = 1e1 * m_p_cgs / UDensity;
    double R_over_h = r/0.3;
    double rho = rho0 / p2(R_over_h);
  
    return rho;
}

double Uniform_model ( double r , const int ipart )
{
    P[ipart].CMZ_type = 0.0;
    double rho = 1e0 * m_p_cgs / UDensity;
    return rho;
}

double NFW_model ( double r, const int ipart )
{
    P[ipart].CMZ_type = 2.0;

    double rho0 = 1e4 * m_p_cgs / UDensity;
    double R_over_h = r/0.1;
    double rho = rho0 * (1 + p2(R_over_h))/ R_over_h;
  
    return rho;
}

double Torus_model (float x, float y, float z, const int ipart )
{
    double rho0 = Problem.Rho_Max;
    double R_in = 0.15;
    double R_out = 0.25;
    double alpha = 5. * pi / 180;
    double r_in = R_in * tan(alpha);
    double r_out = R_out * tan(alpha);
    double flow_R = 0.25;
    double flow_r = 0.05;
    
    double rho = 0;
    
    double R = sqrt(x*x+y*y+z*z);
    double R_flat = sqrt(x*x+y*y);

    if (R < R_in-r_in)
        return 0;

    if  (((sqrt(z * z + pow(R-R_in, 2)) < r_in) && (R_flat < R_in)) ||
        ((sqrt(z * z + pow(R-R_out, 2)) < r_out) && (R_flat > R_out)) ||
        ((fabs(z) < R_flat * tan(alpha)) && (R_flat >= R_in) && (R_flat <= R_out)))
        {
            P[ipart].CMZ_type = 1.0;
            rho += rho0 * exp(-fabs(z)/0.05)* exp(-R_flat/0.1);
            return rho;
        }

    double r_flow = sqrt(pow(fabs(y) - flow_R, 2) + pow(z, 2));

    if  ((r_flow <= flow_r) && (x * y > 0) && (fabs(x)< 1))
        {
            P[ipart].CMZ_type = 2.0;
            if (fabs(x)<0.1)
                rho += rho0/8* exp(-r_flow/0.02) * exp(-(1-fabs(x)/0.1));
            else
                rho += rho0/8* exp(-r_flow/0.02);
            return rho;
        }

    P[ipart].CMZ_type = 0.0;

    double rho_iso = 1e0 * m_p_cgs / UDensity;
    double R_over_h = R/0.5;
    rho += rho_iso / p2(R_over_h);
  
    return rho;

}

double Inflow_model ( float x, float y, float z, const int ipart )
{
    //P[ipart].CMZ_type = 3.0;

    double mu_cgs = 4.77764784e+33; //3.7e7 M_sun*G

    double b = 0.2;
    double e = 1.02;
    double r_p = 0.2;
    double a = -r_p / (e-1);
    double i = pi /180 * 5;
    double Om = 3/2 * pi;
    double w = pi/3;
    /*
    double EA_array = np.linspace(-np.pi/8, 0, 20);

    double nu = 2*np.arctan(np.sqrt((e+1)/(e-1)) * np.tanh(EA/2));

    double r = a*(1 - e*np.cosh(EA));

    double h = np.sqrt(-mu_cgs*a * (e**2-1));

    double X = r*(cos(Om)*cos(w+nu) - sin(Om)*sin(w+nu)*cos(i));
    double Y = r*(sin(Om)*cos(w+nu) + cos(Om)*sin(w+nu)*cos(i));
    double Z = r*(sin(i)*sin(w+nu));

    return [X,Y,Z]
    */
    return 0;
}

float CMZ_Density ( const int ipart , const double bias )
{
    double x = P[ipart].Pos[0] - Problem.Boxsize[0] / 2.0;
    double y = P[ipart].Pos[1] - Problem.Boxsize[1] / 2.0;
    double z = P[ipart].Pos[2] - Problem.Boxsize[2] / 2.0;

    double r = sqrt(x*x+y*y+z*z);

    double rho = 0;

    //rho += Uniform_model(r, ipart);
   
    //rho += Iso_model(r, ipart);

    rho += Torus_model(x, y, z, ipart);
    
    return (float)rho;
}

float CMZ_U ( const int ipart )
{
    double T = 1e4;
    double u_to_temp = 0.63 * m_p_cgs / k_B_cgs * (5.0/3.0-1.0) * UEnergy;
    return T/u_to_temp;
}

void CMZ_Velocity ( const int ipart, float out[3] )
{

    
    double x = P[ipart].Pos[0] / Problem.Boxsize[0] - 0.5;
    double y = P[ipart].Pos[1] / Problem.Boxsize[0] - 0.5;
    double z = P[ipart].Pos[2] / Problem.Boxsize[0] - 0.5;
    double r = sqrt(x*x+y*y+z*z);

    double v = sqrt(G_cgs * 5e7 * 1.988435e33/(r*ULength));

    out[2] = 0;
    out[1] = 0; //x/r * v/UVel;
    out[0] = 0; //-y/r * v/UVel;

}

