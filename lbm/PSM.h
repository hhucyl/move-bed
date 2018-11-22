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
    size_t nz = Ndim(2);
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

inline void Domain::UpdateParticlesContacts()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    ListofContacts.clear();
    std::set<std::pair<int,int>> myset;
    std::set<std::pair<int,int>>* myset_private = new std::set<std::pair<int,int>>[Nproc];
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for(int ip=0; ip<Particles.size(); ++ip)
    {
        DEM::Disk *Pa = &Particles[ip];
        int ixs = std::max(std::floor(Pa->X(0) - Pa->R - 3*dx),0.0);
        int ixe = std::min(std::ceil(Pa->X(0) + Pa->R + 3*dx),(double) nx);
        int iys = std::max(std::floor(Pa->X(1) - Pa->R - 3*dx),0.0);
        int iye = std::min(std::ceil(Pa->X(1) + Pa->R + 3*dx),(double) ny);
        for(int ix=ixs; ix<ixe; ++ix)
        for(int iy=iys; iy<iye; ++iy) 
        {
            double x = (double) ix;
            double y = (double) iy;
            Vec3_t CC(x,y,0);
            double len = DEM::DiskSquare(Pa->X,CC,Pa->R,dx);
            if (std::fabs(len)<1.0e-12) continue;
            if(Check[ix][iy][0]<0)
            {
                Check[ix][iy][0] = ip;
            }else{
                // std::cout<<"Collide!!!!!"<<std::endl;
                int ip1 = std::min(Check[ix][iy][0],ip);
                int ip2 = std::max(Check[ix][iy][0],ip);
                std::pair<int,int> temp(ip1,ip2);
                myset_private[omp_get_thread_num()].insert(temp);
                
            }
        }
    }
    std::set<std::pair<int,int>>::iterator it;
    for(size_t i=0; i<Nproc; ++i)
    {
        //std::set_union(myset_private[i].begin(),myset_private[i].end(),myset.begin(),myset.end(),std::insert_iterator<std::set<std::pair<int,int>>>(myset,myset.begin()));
        for(it=myset_private[i].begin();it!=myset_private[i].end();++it)
        {
           myset.insert(*it);
        }    
    }
    
    for(it=myset.begin();it!=myset.end();++it)
    {
        std::pair<int,int> temp((*it).first,(*it).second);
        ListofContacts.push_back(temp);    
    }
}

inline void Domain::AddDisksG()
{
    if(Time<0.5) std::cout<<"--- "<<"PSM"<<" ---"<<std::endl;    
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(size_t ip=0; ip<Particles.size(); ++ip)
    {
        DEM::Disk *Pa = &Particles[ip];
        // std::cout<<ip<<std::endl;
        int ixs = std::max(std::floor(Pa->X(0) - Pa->Rh - 3*dx),0.0);
        int ixe = std::min(std::ceil(Pa->X(0) + Pa->Rh + 3*dx),(double) nx);
        int iys = std::max(std::floor(Pa->X(1) - Pa->Rh - 3*dx),0.0);
        int iye = std::min(std::ceil(Pa->X(1) + Pa->Rh + 3*dx),(double) ny);
        for(size_t ix=ixs; ix<ixe; ++ix)
        for(size_t iy=iys; iy<iye; ++iy) 
        {
            double x = (double) ix;
            double y = (double) iy;
            Vec3_t CC(x,y,0);
            double len = DEM::DiskSquare(Pa->X,CC,Pa->Rh,dx);
            if (std::fabs(len)<1.0e-12) continue;
            double gamma  = len/(4.0*dx);
            // std::cout<<ix<<" "<<iy<<" "<<iz<<" "<<gamma<<std::endl;
                    
            Gamma[ix][iy][0] = std::min(gamma+Gamma[ix][iy][0],1.0);
            
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
            Pa->F          += Flbmt;
            Pa->T          += T;
        #ifdef USE_OMP
            omp_unset_lock    (&Pa->lck);
        #endif
        }
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