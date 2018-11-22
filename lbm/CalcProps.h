#ifndef LBM_CalcProps_H
#define LBM_CalcProps_H

inline void Domain::CalcProps()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        //BForce[ix][iy][iz] = OrthoSys::O;
        Vel   [ix][iy][iz] = OrthoSys::O;
        Rho   [ix][iy][iz] = 0.0;
        if (!IsSolid[ix][iy][iz])
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Rho[ix][iy][iz] +=  F[ix][iy][iz][k];
                Vel[ix][iy][iz] +=  F[ix][iy][iz][k]*C[k];
            }
            Vel[ix][iy][iz] *= Cs/Rho[ix][iy][iz];
            // Vel[ix][iy][iz] += 0.5*dt*BForce[ix][iy][iz];
        }
    }
}

inline void Domain::CalcPropsForCell(Vec3_t &idx)
{
    size_t ix = (size_t) idx(0);
    size_t iy = (size_t) idx(1);
    size_t iz = (size_t) idx(2);
    Vel   [ix][iy][iz] = OrthoSys::O;
    Rho   [ix][iy][iz] = 0.0;
    if (!IsSolid[ix][iy][iz])
    {
        for (size_t k=0;k<Nneigh;k++)
        {
            Rho[ix][iy][iz] +=  F[ix][iy][iz][k];
            Vel[ix][iy][iz] +=  F[ix][iy][iz][k]*C[k];
        }
        Vel[ix][iy][iz] *= Cs/Rho[ix][iy][iz];
        // Vel[ix][iy][iz] += 0.5*dt*BForce[ix][iy][iz];
    }
}

#endif