#ifndef LBM_FEQ_H
#define LBM_FEQ_H

inline double Domain::Feq(size_t k, double Rho, Vec3_t & V)
{
    double VdotC = dot(V,C[k]);
    double VdotV = dot(V,V);
    return W[k]*Rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
}

#endif