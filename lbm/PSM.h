#ifndef LBM_PSM_H
#define LBM_PSM_H





inline void Domain::AddSphereG(Vec3_t &pos, double R)
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy)
    for(size_t iz=0; iz<nz; ++iz)  
    {
        Vec3_t CC(ix,iy,iz);
        double len = DEM::SphereCube(pos,CC,R,dx);
        if (std::fabs(len)<1.0e-12) continue;
        double gamma  = len/(12.0*dx);
        // std::cout<<ix<<" "<<iy<<" "<<iz<<" "<<gamma<<std::endl;
                
        Gamma[ix][iy][iz] = std::min(gamma,1.0);
        
        //cell->Gamma   = std::max(gamma,cell->Gamma);
        //cell->Gamma   = std::min(gamma+cell->Gamma,1.0);
        //if (fabs(cell->Gamma-1.0)<1.0e-12)
        //if (fabs(cell->Gamma-1.0)<1.0e-12&&(fabs(Lat[0].G)>1.0e-12||Gmix>1.0e-12)) 
        
        // Vec3_t B      = CC - pos;
        // Vec3_t tmp;
        // Rotation(Pa->w,Pa->Q,tmp);
        // Vec3_t VelP   = Pa->v + cross(tmp,B);
        Vec3_t VelPt(0.0,0.0,0.0);
        double rho = Rho[ix][iy][iz];
        double Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
        //double Bn  = gamma;
        //double Bn  = floor(gamma);
        
        Vec3_t Flbmt(0.0,0.0,0.0);
        // double *FF = F[ix][iy][iz];
        
        for (size_t k=0;k<Nneigh;k++)
        {
            double Fvpp     = Feq(Op[k],rho,VelPt);
            double Fvp      = Feq(k    ,rho,VelPt);
            double Omega    = F[ix][iy][iz][Op[k]] - Fvpp - (F[ix][iy][iz][k] - Fvp);
            //cell->Omeis[k] += Omega;
            //cell->Omeis[k] += gamma*Omega;
            // cell->Omeis[k] = Omega;
            Omeis[ix][iy][iz][k] = Omega;
            Flbmt += -Bn*Omega*C[k];
        }
        
        Flbm[ix][iy][iz] = Flbmt;
    }
    
    
}




inline void Domain::AddDiskG(Vec3_t &pos, double R)
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    // size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy) 
    {
        Vec3_t CC(ix,iy,0);
        double len = DEM::DiskSquare(pos,CC,R,dx);
        if (std::fabs(len)<1.0e-12) continue;
        double gamma  = len/(4.0*dx);
        // std::cout<<ix<<" "<<iy<<" "<<iz<<" "<<gamma<<std::endl;
                
        Gamma[ix][iy][0] = std::min(gamma,1.0);
        
        //cell->Gamma   = std::max(gamma,cell->Gamma);
        //cell->Gamma   = std::min(gamma+cell->Gamma,1.0);
        //if (fabs(cell->Gamma-1.0)<1.0e-12)
        //if (fabs(cell->Gamma-1.0)<1.0e-12&&(fabs(Lat[0].G)>1.0e-12||Gmix>1.0e-12)) 
        
        // Vec3_t B      = CC - pos;
        // Vec3_t tmp;
        // Rotation(Pa->w,Pa->Q,tmp);
        // Vec3_t VelP   = Pa->v + cross(tmp,B);
        Vec3_t VelPt(0.0,0.0,0.0);
        double rho = Rho[ix][iy][0];
        double Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
        //double Bn  = gamma;
        //double Bn  = floor(gamma);
        
        Vec3_t Flbmt(0.0,0.0,0.0);
        // double *FF = F[ix][iy][iz];
        
        for (size_t k=0;k<Nneigh;k++)
        {
            double Fvpp     = Feq(Op[k],rho,VelPt);
            double Fvp      = Feq(k    ,rho,VelPt);
            double Omega    = F[ix][iy][0][Op[k]] - Fvpp - (F[ix][iy][0][k] - Fvp);
            //cell->Omeis[k] += Omega;
            //cell->Omeis[k] += gamma*Omega;
            // cell->Omeis[k] = Omega;
            Omeis[ix][iy][0][k] = Omega;
            Flbmt += -Bn*Omega*C[k];
        }
        
        Flbm[ix][iy][0] = Flbmt;
    }
    
    
}



inline void Domain::adddiskG_sub(DEM::Disk *Pa)
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    // size_t nz = Ndim(2);
    int ixs = std::max(std::floor(Pa->X(0) - Pa->Rh - 3*dx),0.0);
    int ixe = std::min(std::ceil(Pa->X(0) + Pa->Rh + 3*dx),(double) nx);
    int iys = std::max(std::floor(Pa->X(1) - Pa->Rh - 3*dx),0.0);
    int iye = std::min(std::ceil(Pa->X(1) + Pa->Rh + 3*dx),(double) ny);
    // std::cout<<ixs<<" "<<ixe<<" "<<iys<<" "<<iye<<std::endl; 
    for(int ix=ixs; ix<ixe; ++ix)
    for(int iy=iys; iy<iye; ++iy) 
    {
        double x = (double) ix;
        double y = (double) iy;
        Vec3_t CC(x,y,0);
        double len = DEM::DiskSquare(Pa->X,CC,Pa->Rh,dx);
        // std::cout<<ix<<" "<<iy<<" "<<len<<std::endl;
        if (std::fabs(len)<1.0e-12) continue;
        double gamma  = len/(4.0*dx);
        // std::cout<<ix<<" "<<iy<<" "<<iz<<" "<<gamma<<std::endl;
                
        Gamma[ix][iy][0] = std::min(gamma+Gamma[ix][iy][0],1.0);
        // CheckRh[ix][iy][0] = ip;
        //cell->Gamma   = std::max(gamma,cell->Gamma);
        //cell->Gamma   = std::min(gamma+cell->Gamma,1.0);
        //if (fabs(cell->Gamma-1.0)<1.0e-12)
        //if (fabs(cell->Gamma-1.0)<1.0e-12&&(fabs(Lat[0].G)>1.0e-12||Gmix>1.0e-12)) 
        
        Vec3_t B      = CC - Pa->X;
        Vec3_t tmp;
        Rotation(Pa->W,Pa->Q,tmp);
        Vec3_t VelPt   = Pa->V + cross(tmp,B);
        VelP[ix][iy][0] = VelPt;
        double rho = Rho[ix][iy][0];
        double Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
        //double Bn  = gamma;
        //double Bn  = floor(gamma);
        
        Vec3_t Flbmt(0.0,0.0,0.0);
        // double *FF = F[ix][iy][iz];
        
        for (size_t k=0;k<Nneigh;k++)
        {
            double Fvpp     = Feq(Op[k],rho,VelPt);
            double Fvp      = Feq(k    ,rho,VelPt);
            double Omega    = F[ix][iy][0][Op[k]] - Fvpp - (F[ix][iy][0][k] - Fvp);
            //cell->Omeis[k] += Omega;
            //cell->Omeis[k] += gamma*Omega;
            // cell->Omeis[k] = Omega;
            Omeis[ix][iy][0][k] = Omega;
            Flbmt += -Bn*Omega*C[k];
        }
        
        Flbm[ix][iy][0] = Flbmt;
        Vec3_t T,Tt;
        Tt =           cross(B,Flbmt);
        Quaternion_t q;
        Conjugate    (Pa->Q,q);
        Rotation     (Tt,q,T);
            //std::cout << "1" << std::endl;
    #ifdef USE_OMP
        omp_set_lock      (&Pa->lck);
    #endif
        Pa->Fh          += Flbmt;
        Pa->Th          += T;
    #ifdef USE_OMP
        omp_unset_lock    (&Pa->lck);
    #endif
    }
}

inline void Domain::AddDisksG()
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
        adddiskG_sub(&Particles[ip]);
        // std::cout<<ip<<std::endl;
        
    }
    // std::cout<<2<<std::endl;

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t ip=0; ip<GhostParticles.size(); ++ip)
    {
        DEM::Disk *Pa = &GhostParticles[ip];
        Pa->Fh = 0.0,0.0,0.0;
        Pa->Th = 0.0,0.0,0.0;
        if(!Pa->Ghost) continue;
        // std::cout<<ip<<std::endl;
        adddiskG_sub(&GhostParticles[ip]);
        
    }
    // std::cout<<3<<std::endl;

    
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t ip=0; ip<Particles.size(); ++ip)
    {
        DEM::Disk *Pa = &Particles[ip];
        DEM::Disk *Pa_ghost = &GhostParticles[ip];
        if(!Pa_ghost->Ghost) continue;
        Pa->Fh += Pa_ghost->Fh;
        Pa->Th += Pa_ghost->Th;
    }
    
}


inline void Domain::BoundaryGamma()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy)
    for(size_t iz=0; iz<nz; ++iz)  
    {
        
        double gamma  = Gamma[ix][iy][iz];
        if(gamma<1e-12) continue;
        
        Vec3_t VelPt = VelP[ix][iy][iz];
        double rho = Rho[ix][iy][iz];
        double Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
        //double Bn  = gamma;
        //double Bn  = floor(gamma);
        
        Vec3_t Flbmt(0.0,0.0,0.0);
        
        for (size_t k=0;k<Nneigh;k++)
        {
            double Fvpp     = Feq(Op[k],rho,VelPt);
            double Fvp      = Feq(k    ,rho,VelPt);
            double Omega    = F[ix][iy][iz][Op[k]] - Fvpp - (F[ix][iy][iz][k] - Fvp);
            //cell->Omeis[k] += Omega;
            //cell->Omeis[k] += gamma*Omega;
            // cell->Omeis[k] = Omega;
            Omeis[ix][iy][iz][k] = Omega;
            Flbmt += -Bn*Omega*C[k];
        }
        
        Flbm[ix][iy][iz] = Flbmt;
    }
    
    
}

























#endif