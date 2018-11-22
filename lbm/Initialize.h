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

#endif