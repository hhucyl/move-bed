#ifndef LBM_STREAM_H
#define LBM_STREAM_H

inline void Domain::Stream()
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
        for (size_t k=0;k<Nneigh;k++)
        {
            size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)Ndim(0))%Ndim(0);
            size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)Ndim(1))%Ndim(1);
            size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)Ndim(2))%Ndim(2);
            F[nix][niy][niz][k] = Ftemp[ix][iy][iz][k];
        }
    }

}

#endif