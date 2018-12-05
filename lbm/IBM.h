#ifndef LBM_IBM_H
#define LBM_IBM_H

#include <mechsys/dem/special_functions.h>


inline double Domain::KernelIBM(double r, double x)
{
    double xx = (r-x)/dx;
    double ll = std::sqrt(xx*xx);
    if(ll>=2.0)
    {
        return 0;
    }else if(ll>=1.0)
    {
        return 0.125*(5.0 - 2.0*ll - std::sqrt(-7.0 + 12.0*ll - 4.0*ll*ll));
    }else{
        return 0.125*(3.0 - 2.0*ll + std::sqrt(1.0 + 4.0*ll - 4.0*ll*ll));
        
    }
}

inline double Domain::KernelIBM1(double r, double x)
{
    double xx = (r-x)/dx;
    double ll = std::sqrt(xx*xx);
    if(ll>1.0)
    {
        return 0;
    }else{
        return 1-std::fabs(xx);
        
    }
}


inline void Domain::ApplyIBM2D(Vec3_t &pos, double R)
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    // size_t nz = Ndim(2);
    int N = points.size();
    Vec3_t VelIBM[N];
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int im=0; im<N; im++)
    {
        Vec3_t r = points[im];
        int ixs = std::max(std::floor(r(0) - 3*dx),0.0);
        int ixe = std::min(std::ceil(r(0) + 3*dx),(double) nx);
        int iys = std::max(std::floor(r(1) - 3*dx),0.0);
        int iye = std::min(std::ceil(r(1) + 3*dx),(double) ny);
        VelIBM[im] = 0.0,0.0,0.0;
        
        for(int ix= ixs; ix<ixe; ix++)
        for(int iy= iys; iy<iye; iy++)
        {
            VelIBM[im] += Vel[ix][iy][0]*KernelIBM(r(0),ix)*KernelIBM(r(1),iy); 
        }
        
        
    }
    int ixs = std::max(std::floor(pos(0) - (R+3)*dx),0.0);
    int ixe = std::min(std::ceil(pos(0) + (R+3)*dx),(double) nx);
    int iys = std::max(std::floor(pos(1) - (R+3)*dx),0.0);
    int iye = std::min(std::ceil(pos(1) + (R+3)*dx),(double) ny);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int ix= ixs; ix<ixe; ix++)
    for(int iy= iys; iy<iye; iy++)
    {
        Flbm[ix][iy][0] = 0.0,0.0,0.0;
        for(int im=0; im<N; im++)
        {
            Vec3_t r = points[im];        
            Vec3_t FIBM = 2.0*Rho[ix][iy][0]*(0.0-VelIBM[im])/dt;
            Flbm[ix][iy][0] += FIBM*KernelIBM(r(0),ix)*KernelIBM(r(1),iy)/(dx*dx)*dS[im]; 
        }
        
    }
    
}

inline void Domain::GenPts(Vec3_t &pos, double R, int N)
{
    if(Ndim(2)==1)
    {
        double alpha = 2*M_PI/((double) N);
        for(int im=0; im<N; im++)
        {
            Vec3_t r(R*std::cos(im*alpha)+pos(0),R*std::sin(im*alpha)+pos(1),0.0);
            points.push_back(r);
            dS.push_back(alpha*R);
        }
    }else{
        // Rb = std::pow(0.5*(R*R*R + (R-dx)*(R-dx)*(R-dx)),1.0/3.0);
        // int ns = std::round(M_PI*Rb+1);
        
    }

}

inline void Domain::ApplyIBM3D(Vec3_t &pos, double R)
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    int N = points.size();
    Vec3_t VelIBM[N];
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int im=0; im<N; im++)
    {
        Vec3_t r = points[im];
        int ixs = std::max(std::floor(r(0) - 3*dx),0.0);
        int ixe = std::min(std::ceil(r(0) + 3*dx),(double) nx);
        int iys = std::max(std::floor(r(1) - 3*dx),0.0);
        int iye = std::min(std::ceil(r(1) + 3*dx),(double) ny);
        int izs = std::max(std::floor(r(2) - 3*dx),0.0);
        int ize = std::min(std::ceil(r(2) + 3*dx),(double) nz);
        VelIBM[im] = 0.0,0.0,0.0;
        
        for(int ix= ixs; ix<ixe; ix++)
        for(int iy= iys; iy<iye; iy++)
        for(int iz= izs; iz<ize; iz++)
        {
            VelIBM[im] += Vel[ix][iy][iz]*KernelIBM1(r(0),ix)*KernelIBM1(r(1),iy)*KernelIBM1(r(2),iz); 
        }
        
        
    }
    int ixs = std::max(std::floor(pos(0) - (R+3)*dx),0.0);
    int ixe = std::min(std::ceil(pos(0) + (R+3)*dx),(double) nx);
    int iys = std::max(std::floor(pos(1) - (R+3)*dx),0.0);
    int iye = std::min(std::ceil(pos(1) + (R+3)*dx),(double) ny);
    int izs = std::max(std::floor(pos(2) - (R+3)*dx),0.0);
    int ize = std::min(std::ceil(pos(2) + (R+3)*dx),(double) nz);
    /*int ixs = 0;
    int ixe = nx;
    int iys = 0;
    int iye = ny;
    int izs = 0;
    int ize = nz;*/
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int ix= ixs; ix<ixe; ix++)
    for(int iy= iys; iy<iye; iy++)
    for(int iz= izs; iz<ize; iz++)
    // for(size_t ix=0; ix<nx; ++ix)
    // for(size_t iy=0; iy<ny; ++iy)
    // for(size_t iz=0; iz<nz; ++iz)  
    {
        // Flbm[ix][iy][iz] = 0.0,0.0,0.0;
        for(int im=0; im<N; im++)
        {
            Vec3_t r = points[im];
            if(std::fabs(KernelIBM1(r(0),ix)*KernelIBM1(r(1),iy)*KernelIBM1(r(2),iz))<1e-9) continue;       
            Vec3_t FIBM = 2.0*Rho[ix][iy][iz]*(0.0-VelIBM[im])/dt;
            Flbm[ix][iy][iz] += FIBM*KernelIBM1(r(0),ix)*KernelIBM1(r(1),iy)*KernelIBM1(r(2),iz)/(dx*dx*dx)*dS[im]; 
        }
        Vel[ix][iy][iz] += dt/(2.0*Rho[ix][iy][iz])*Flbm[ix][iy][iz];
        
    }
    
}

inline void Domain::ApplyIBM3D()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    int N = points.size();
    Vec3_t VelIBM[N];
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int im=0; im<N; im++)
    {
        Vec3_t r = points[im];
        int ixs = std::max(std::floor(r(0) - 3*dx),0.0);
        int ixe = std::min(std::ceil(r(0) + 3*dx),(double) nx);
        int iys = std::max(std::floor(r(1) - 3*dx),0.0);
        int iye = std::min(std::ceil(r(1) + 3*dx),(double) ny);
        int izs = std::max(std::floor(r(2) - 3*dx),0.0);
        int ize = std::min(std::ceil(r(2) + 3*dx),(double) nz);
        VelIBM[im] = 0.0,0.0,0.0;
        
        for(int ix= ixs; ix<ixe; ix++)
        for(int iy= iys; iy<iye; iy++)
        for(int iz= izs; iz<ize; iz++)
        {
            VelIBM[im] += Vel[ix][iy][iz]*KernelIBM1(r(0),ix)*KernelIBM1(r(1),iy)*KernelIBM1(r(2),iz); 
        }
        
        
    }
    // int ixs = std::max(std::floor(pos(0) - (R+3)*dx),0.0);
    // int ixe = std::min(std::ceil(pos(0) + (R+3)*dx),(double) nx);
    // int iys = std::max(std::floor(pos(1) - (R+3)*dx),0.0);
    // int iye = std::min(std::ceil(pos(1) + (R+3)*dx),(double) ny);
    // int izs = std::max(std::floor(pos(2) - (R+3)*dx),0.0);
    // int ize = std::min(std::ceil(pos(2) + (R+3)*dx),(double) nz);
    int ixs = 0;
    int ixe = nx;
    int iys = 0;
    int iye = ny;
    int izs = 0;
    int ize = nz;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int ix= ixs; ix<ixe; ix++)
    for(int iy= iys; iy<iye; iy++)
    for(int iz= izs; iz<ize; iz++)
    // for(size_t ix=0; ix<nx; ++ix)
    // for(size_t iy=0; iy<ny; ++iy)
    // for(size_t iz=0; iz<nz; ++iz)  
    {
        // Flbm[ix][iy][iz] = 0.0,0.0,0.0;
        for(int im=0; im<N; im++)
        {
            Vec3_t r = points[im];
            if(std::fabs(KernelIBM1(r(0),ix)*KernelIBM1(r(1),iy)*KernelIBM1(r(2),iz))<1e-9) continue;       
            Vec3_t FIBM = 2.0*Rho[ix][iy][iz]*(0.0-VelIBM[im])/dt;
            Flbm[ix][iy][iz] += FIBM*KernelIBM1(r(0),ix)*KernelIBM1(r(1),iy)*KernelIBM1(r(2),iz)/(dx*dx*dx)*dS[im]; 
        }
        Vel[ix][iy][iz] += dt/(2.0*Rho[ix][iy][iz])*Flbm[ix][iy][iz];
        
    }
    
}

inline void Domain::ApplyIBM3DIM(Vec3_t &pos, double R, size_t IT)
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    int N = points.size();
    Vec3_t VelIBM[N];
    int ixss = std::max(std::floor(pos(0) - (R+3)*dx),0.0);
    int ixee = std::min(std::ceil(pos(0) + (R+3)*dx),(double) nx);
    int iyss = std::max(std::floor(pos(1) - (R+3)*dx),0.0);
    int iyee = std::min(std::ceil(pos(1) + (R+3)*dx),(double) ny);
    int izss = std::max(std::floor(pos(2) - (R+3)*dx),0.0);
    int izee = std::min(std::ceil(pos(2) + (R+3)*dx),(double) nz);
    for(size_t it=0; it<IT; it++)
    {
        #ifdef USE_OMP
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        #endif
        for(int im=0; im<N; im++)
        {
            Vec3_t r = points[im];
            int ixs = std::max(std::floor(r(0) - 3*dx),0.0);
            int ixe = std::min(std::ceil(r(0) + 3*dx),(double) nx);
            int iys = std::max(std::floor(r(1) - 3*dx),0.0);
            int iye = std::min(std::ceil(r(1) + 3*dx),(double) ny);
            int izs = std::max(std::floor(r(2) - 3*dx),0.0);
            int ize = std::min(std::ceil(r(2) + 3*dx),(double) nz);
            VelIBM[im] = 0.0,0.0,0.0;
            
            for(int ix= ixs; ix<ixe; ix++)
            for(int iy= iys; iy<iye; iy++)
            for(int iz= izs; iz<ize; iz++)
            {
                VelIBM[im] += Vel[ix][iy][iz]*KernelIBM1(r(0),ix)*KernelIBM1(r(1),iy)*KernelIBM1(r(2),iz); 
            }
            
            
        }
        
        
        
        #ifdef USE_OMP
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        #endif
        for(int ix= ixss; ix<ixee; ix++)
        for(int iy= iyss; iy<iyee; iy++)
        for(int iz= izss; iz<izee; iz++)
        // for(size_t ix=0; ix<nx; ++ix)
        // for(size_t iy=0; iy<ny; ++iy)
        // for(size_t iz=0; iz<nz; ++iz)  
        {
            // Flbm[ix][iy][iz] = 0.0,0.0,0.0;
            Vec3_t Flbmt(0.0,0.0,0.0);
            for(int im=0; im<N; im++)
            {
                Vec3_t r = points[im];
                if(std::fabs(KernelIBM1(r(0),ix)*KernelIBM1(r(1),iy)*KernelIBM1(r(2),iz))<1e-9) continue;       
                Vec3_t FIBM = 2.0*Rho[ix][iy][iz]*(0.0-VelIBM[im])/dt;
                Flbmt += FIBM*KernelIBM1(r(0),ix)*KernelIBM1(r(1),iy)*KernelIBM1(r(2),iz)/(dx*dx*dx)*dS[im]; 
                
            }
            // if(ix == 14&&iy ==14&&iz ==1) std::cout<<Time<<" "<<Flbmt<<" "<<Flbm[14][14][1]<<std::endl;
            Vel[ix][iy][iz] += dt/(2.0*Rho[ix][iy][iz])*Flbmt;
            Flbm[ix][iy][iz] += Flbmt;
        }
        // std::cout<<"IT "<<IT<<" "<<Flbm[14][14][1]<<std::endl;
    }
}
    

























#endif