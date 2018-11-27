#ifndef LBM_CoupleDEM_H
#define LBM_CoupleDEM_H

inline void Domain::join_contactlist_sub(std::set<std::pair<int,int>> *myset_private, std::vector<std::pair<int,int>> &ListofContacts)
{
    std::set<std::pair<int,int>> myset;
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

inline void Domain::UpdateParticlesContacts()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    ListofContactsPP.clear();
    std::set<std::pair<int,int>> myset_privatepp[Nproc];
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
                myset_privatepp[omp_get_thread_num()].insert(temp);
                
            }
        }
    }
    join_contactlist_sub(myset_privatepp,ListofContactsPP);

    ListofContactsPG.clear();
    std::set<std::pair<int,int>> myset_privatepg[Nproc];
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for(int ip=0; ip<GhostParticles.size(); ++ip)
    {
        DEM::Disk *Pa = &GhostParticles[ip];
        if(!Pa->Ghost) continue;
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
            if(Check[ix][iy][0]>0)
            {
                // std::cout<<"Collide!!!!!"<<std::endl;
                int ip1 = std::min(Check[ix][iy][0],ip);
                int ip2 = std::max(Check[ix][iy][0],ip);
                if(std::fabs(ip2-ip1)<1e-6) continue;
                std::pair<int,int> temp(ip1,ip2);
                myset_privatepg[omp_get_thread_num()].insert(temp);
                
            }
        }
    }
    join_contactlist_sub(myset_privatepg,ListofContactsPG);

    
}

inline void Domain::update_pair_sub(DEM::DiskPair &pair, DEM::Disk* P1, DEM::Disk* P2)
{
    if(pair.delta>0)
    {
        omp_set_lock  (&P1->lck);
            P1->Fc += pair.F1;
            P1->Tc += pair.T1;
        omp_unset_lock(&P1->lck);
        omp_set_lock  (&P2->lck);
            P2->Fc += pair.F2;
            P2->Tc += pair.T2;
        omp_unset_lock(&P2->lck);
    }
    
}

inline void Domain::UpdateParticlePairForce()
{
    //ordinary particle
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for(size_t i=0; i<ListofContactsPP.size();++i)
    {
        int ip1 = ListofContactsPP[i].first;
        int ip2 = ListofContactsPP[i].second;
        if(!Particles[ip1].IsFree() && !Particles[ip2].IsFree()) continue;
        
        DEM::DiskPair pair(&Particles[ip1],&Particles[ip2]);
        pair.CalcForce(dtdem);
        
        update_pair_sub(pair,&Particles[ip1],&Particles[ip2]);    
        
        
    }
    //ghost particle
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for(size_t i=0; i<ListofContactsPG.size();++i)
    {
        int ip1 = ListofContactsPG[i].first;
        int ip2 = ListofContactsPG[i].second;
        if(!Particles[ip1].IsFree() && !Particles[ip2].IsFree()) continue;
        bool flag1 = GhostParticles[ip1].Ghost;
        bool flag2 = GhostParticles[ip2].Ghost;
        if(flag1&&!flag2)
        {
            DEM::DiskPair pair(&GhostParticles[ip1],&Particles[ip2]);
            pair.CalcForce(dtdem);
            
            update_pair_sub(pair,&Particles[ip1],&Particles[ip2]);    
            
        }
        if(!flag1&&flag2)
        {
            DEM::DiskPair pair(&Particles[ip1],&GhostParticles[ip2]);
            pair.CalcForce(dtdem);
            
            update_pair_sub(pair,&Particles[ip1],&Particles[ip2]); 
        }

    }
}

inline void Domain::MoveParticles()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0; i<Particles.size(); i++)
    {
        Particles[i].Translate(dtdem);
        Particles[i].Rotate(dtdem);
    }
}


inline void Domain::LeaveAndForcedForce()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for(size_t i=0;i<Particles.size();i++)
    {
        Particles[i].Fh = 0.0,0.0,0.0;
        Particles[i].Fc = 0.0,0.0,0.0;
        Particles[i].Th = 0.0,0.0,0.0;
        Particles[i].Tc = 0.0,0.0,0.0;
        
        if(modexy<0) continue;
        Particles[i].Leave(modexy,Box);
    }
}

inline void Domain::GhostPeriodic()
{   
    
    if(modexy<0) return;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for(size_t i=0; i<GhostParticles.size();++i)
    {
        DEM::Disk *Pa = &GhostParticles[i];
        Pa->Periodic(modexy,Box);
        Pa->Ff = 0.0,0.0,0.0;
        Pa->F = 0.0,0.0,0.0;
        Pa->Tf = 0.0,0.0,0.0;
        Pa->T = 0.0,0.0,0.0;

    }
    
    
}


#endif