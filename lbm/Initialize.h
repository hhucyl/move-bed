#ifndef LBM_INITIALIZE_H
#define LBM_INITIALIZE_H

inline void Domain::Initialize(iVec3_t idx, double TheRho, Vec3_t & TheVel)
{
    size_t ix = idx(0);
    size_t iy = idx(1);
    size_t iz = idx(2);

    BForce[ix][iy][iz] = OrthoSys::O;

    for (size_t k=0;k<Nneigh;k++)
    {
        F[ix][iy][iz][k] = Feq(k,TheRho,TheVel);
    }

    if (!IsSolid[ix][iy][iz])
    {
        Vel[ix][iy][iz] = TheVel;
        Rho[ix][iy][iz] = TheRho;
    }
    else
    {
        Vel[ix][iy][iz] = OrthoSys::O;
        Rho[ix][iy][iz] = 0.0;
    }
}

inline void Domain::Initial(double rho, Vec3_t &v0,  Vec3_t &g0)
{
    
    for(size_t ix=0; ix<Ndim(0); ix++)
    for(size_t iy=0; iy<Ndim(1); iy++)
    for(size_t iz=0; iz<Ndim(2); iz++)
    {
        Rho[ix][iy][iz] = rho;
        Vel[ix][iy][iz] = 0.0, 0.0, 0.0;
        BForce[ix][iy][iz] = g0;
        for(size_t k=0; k<Nneigh; ++k)
        {
            F[ix][iy][iz][k] = Feq(k,rho,v0);            
            Ftemp[ix][iy][iz][k] = Feq(k,rho,v0);            
        }
    // std::cout<<F[ix][iy][iz][18]<<std::endl;
        
    }
    Rho0 = rho;//very important
}

#endif