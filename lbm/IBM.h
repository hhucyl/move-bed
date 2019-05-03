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
    

inline void Domain::adddiskIBM_sub(DEM::Disk *Pa)
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    // size_t nz = Ndim(2);
    int N = std::ceil(2.0*M_PI*Pa->Rh/dx);
    double alpha = 2*M_PI/((double) N);
    double dS = 2.0*M_PI*Pa->Rh*dx/N;
    Vec3_t points[N];
    Vec3_t VelIBM[N];
    Vec3_t FIBM[N];
    for(int im=0; im<N; im++)
    {
        Vec3_t r(Pa->Rh*std::cos(im*alpha)+Pa->X(0),Pa->Rh*std::sin(im*alpha)+Pa->X(1),0.0);
        points[im] = r;
        int ixs = std::max(std::floor(r(0) - 3*dx),0.0);
        int ixe = std::min(std::ceil(r(0) + 3*dx),(double) nx);
        int iys = std::max(std::floor(r(1) - 3*dx),0.0);
        int iye = std::min(std::ceil(r(1) + 3*dx),(double) ny);
        VelIBM[im] = 0.0,0.0,0.0;
        for(int ix= ixs; ix<ixe; ix++)
        for(int iy= iys; iy<iye; iy++)
        {
            VelIBM[im] += Vel[ix][iy][0]*KernelIBM1(r(0),ix)*KernelIBM1(r(1),iy); 
        }
        Vec3_t tmp;
        Rotation(Pa->W,Pa->Q,tmp);     
        Vec3_t VelPt   = Pa->V + cross(tmp,r);
        Vec3_t FIBMt = -(VelPt-VelIBM[im])/dt*dS;
        FIBM[im] = -FIBMt;
        Vec3_t T,Tt;
        Vec3_t B = r-Pa->X;
        Tt =           cross(B,FIBMt);
        Quaternion_t q;
        Conjugate    (Pa->Q,q);
        Rotation     (Tt,q,T);
            //std::cout << "1" << std::endl;
    #ifdef USE_OMP
        omp_set_lock      (&Pa->lck);
    #endif
        Pa->Fh          += FIBMt;
        Pa->Th          += T;
    #ifdef USE_OMP
        omp_unset_lock    (&Pa->lck);
    #endif
        
    // }
    // int ixs = std::max(std::floor(Pa->X(0) - (Pa->Rh+3)*dx),0.0);
    // int ixe = std::min(std::ceil(Pa->X(0) + (Pa->Rh+3)*dx),(double) nx);
    // int iys = std::max(std::floor(Pa->X(1) - (Pa->Rh+3)*dx),0.0);
    // int iye = std::min(std::ceil(Pa->X(1) + (Pa->Rh+3)*dx),(double) ny);
    // for(int im=0; im<N; im++)
    // {
        for(int ix= ixs; ix<ixe; ix++)
        for(int iy= iys; iy<iye; iy++)
        {
        
        
        // Flbm[ix][iy][0] = 0.0,0.0,0.0;
        
            Vec3_t r = points[im];
            // Vec3_t FIBM = 2.0*Rho[ix][iy][0]*(VelPt-VelIBM[im])/dt;
            Flbm[ix][iy][0] += 2.0*Rho[ix][iy][0]*FIBM[im]*KernelIBM1(r(0),ix)*KernelIBM1(r(1),iy)/(dx*dx)*dS; 
        }
        
       
    }
    
}

inline void Domain::FindNeighNodes(Vec3_t &r, Vec3_t &gr, int depth, std::vector<std::pair<int,int>> &NodeList, std::vector<bool> &NodeType, int ip)
{
    int nx = (int) Ndim(0);
    int ny = (int) Ndim(1);
    
    int ixs = std::max(std::floor(r(0)-depth),0.0);
    int ixe = std::min(std::ceil(r(0)+depth),(double) nx);
    int iys = std::max(std::floor(r(1)-depth),0.0);
    int iye = std::min(std::ceil(r(1)+depth),(double) ny);
    int ixsg = std::max(std::floor(gr(0)-depth),0.0);
    int ixeg = std::min(std::ceil(gr(0)+depth),(double) nx);
    int iysg = std::max(std::floor(gr(1)-depth),0.0);
    int iyeg = std::min(std::ceil(gr(1)+depth),(double) ny);
    // Vec3_t IS(ixs,iys,0);
    // Vec3_t IE(ixe,iye,0);
    // bool flag = r(modexy)>=Box(0) && r(modexy)<=Box(1);
    std::pair<int,int> node;
    
   
    for(int ix=ixs; ix<ixe; ++ix)
    for(int iy=iys; iy<iye; ++iy)
    {
        node.first = ix;
        node.second = iy;
        NodeList.push_back(node);
        NodeType.push_back(0==0);
        // CheckIBM[ix][iy][0].insert(ip);
    }

    for(int ix=ixsg; ix<ixeg; ++ix)
    for(int iy=iysg; iy<iyeg; ++iy)
    {
        node.first = ix;
        node.second = iy;
        NodeList.push_back(node);
        NodeType.push_back(0==1);
        // CheckIBM[ix][iy][0].insert(ip);
    }
    
    
}


inline void Domain::adddiskIBM_sub_periodic(DEM::Disk *Pa, DEM::Disk *GPa, int ip)
{
    // int nx = (int) Ndim(0);
    // int ny = (int) Ndim(1);
    // size_t nz = Ndim(2);
    int N = std::ceil(2.0*M_PI*Pa->Rh/dx);
    // N = 240;
    double alpha = 2*M_PI/((double) N);
    double dS = 2.0*M_PI*Pa->Rh*dx/N;
    // dS = 0.785;
    // Vec3_t points[N];
    Vec3_t VelIBM[N];
    Vec3_t FIBM[N];

    for(int im=0; im<N; im++)
    {
        Vec3_t r(Pa->Rh*std::cos(im*alpha)+Pa->X(0),Pa->Rh*std::sin(im*alpha)+Pa->X(1),0.0);
        Vec3_t gr(GPa->Rh*std::cos(im*alpha)+GPa->X(0),GPa->Rh*std::sin(im*alpha)+GPa->X(1),0.0);
        // Pa->points[im] = r;
        // GPa->points[im] = gr;
        
        VelIBM[im] = 0.0,0.0,0.0;
        std::vector<std::pair<int,int>> NodeList;
        std::vector<bool> NodeType;
        FindNeighNodes(r, gr, 3, NodeList, NodeType, ip);
        // Pa->NodeType.assign(NodeType.begin(),NodeType.end());
        Vec3_t r_temp(0.0,0.0,0.0);
        for(int i=0; i<(int) NodeList.size(); ++i)
        {
            int ix = NodeList[i].first;
            int iy = NodeList[i].second;
            if(NodeType[i])
            {
                r_temp = r;
            }else{
                r_temp = gr;
            }
            VelIBM[im] += Vel[ix][iy][0]*KernelIBM1(r_temp(0),ix)*KernelIBM1(r_temp(1),iy); 
        }
        Vec3_t B = r-Pa->X;
        Vec3_t tmp;
        Rotation(Pa->W,Pa->Q,tmp);
        Vec3_t VelPt   = Pa->V + cross(tmp,B);
        // Vec3_t VelPt   = Pa->V;
        Vec3_t FIBMt = -(VelPt-VelIBM[im])/dt*dS;
        FIBM[im] = -FIBMt;
        Vec3_t T,Tt;
        Tt =           cross(B,FIBMt);
        Quaternion_t q;
        Conjugate    (Pa->Q,q);
        Rotation     (Tt,q,T);
            //std::cout << "1" << std::endl;
    // #ifdef USE_OMP
    //     omp_set_lock      (&Pa->lck);
    // #endif
        Pa->Fh          += FIBMt;
        Pa->Th          += Tt;
    // #ifdef USE_OMP
    //     omp_unset_lock    (&Pa->lck);
    // #endif
      
    
        for(int i=0; i<(int) NodeList.size(); ++i)
        {
            int ix = NodeList[i].first;
            int iy = NodeList[i].second;
            if(NodeType[i])
            {
                r_temp = r;
            }else{
                r_temp = gr;
            }
            Flbm[ix][iy][0] += 1.0*Rho[ix][iy][0]*FIBM[im]*KernelIBM1(r_temp(0),ix)*KernelIBM1(r_temp(1),iy)/(dx*dx)*dS;//may cause trouble when kernel overlap 
        }
    }
    
}





inline void Domain::AddDisksIBM()
{
    // if(Time<0.5) std::cout<<"--- "<<"PSM"<<" ---"<<std::endl;    
    // std::cout<<1<<std::endl;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t ip=0; ip<Particles.size(); ++ip)
    {
        Particles[ip].Fh = 0.0,0.0,0.0;
        Particles[ip].Th = 0.0,0.0,0.0;
        GhostParticles[ip].Fh = 0.0,0.0,0.0;
        GhostParticles[ip].Th = 0.0,0.0,0.0;
        // adddiskIBM_sub(&Particles[ip]);
        adddiskIBM_sub_periodic(&Particles[ip],&GhostParticles[ip],(int) ip);
        // std::cout<<ip<<std::endl;
        
    }
    // std::cout<<2<<std::endl;

    // #ifdef USE_OMP
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    // #endif
    // for(size_t ip=0; ip<GhostParticles.size(); ++ip)
    // {
    //     DEM::Disk *Pa = &GhostParticles[ip];
    //     Pa->Fh = 0.0,0.0,0.0;
    //     Pa->Th = 0.0,0.0,0.0;
    //     if(!Pa->Ghost) continue;
    //     // std::cout<<ip<<std::endl;
    //     adddiskIBM_sub(&GhostParticles[ip]);
        
    // }
    // // std::cout<<3<<std::endl;

    
    // #ifdef USE_OMP
    // #pragma omp parallel for schedule(static) num_threads(Nproc)
    // #endif
    // for(size_t ip=0; ip<Particles.size(); ++ip)
    // {
    //     DEM::Disk *Pa = &Particles[ip];
    //     DEM::Disk *Pa_ghost = &GhostParticles[ip];
    //     if(!Pa_ghost->Ghost) continue;
    //     Pa->Fh += Pa_ghost->Fh;
    //     Pa->Th += Pa_ghost->Th;
    // }
    
}


















#endif