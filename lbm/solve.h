#ifndef LBM_SOLVE_H
#define LBM_SOLVE_H

inline void Domain::Solve(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    std::cout<<"dt of LBM "<<dt<<" dt of DEM "<<dtdem<<std::endl;
    double tout = 0;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    GhostParticles.assign(Particles.begin(), Particles.end()); 
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, idx_out);
            
            WriteXDMF(fn.CStr());
            idx_out++;
            
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            
            

            tout += dtout;
        }
        

        


        for(int i=0; i<std::floor(dt/dtdem); ++i)
        {
            // std::cout<<i<<std::endl;
            bool flag = i==0 || i==(std::floor(dt/dtdem)-1);
        
            if(flag){
                SetZero();
            }
            //set added force and check leave particles
            LeaveAndForcedForce();

            GhostParticles.clear();
            GhostParticles.assign(Particles.begin(), Particles.end()); 
        
            GhostPeriodic();
        
         
            //set fluid force
            if(flag){
                AddDisksG();
            }

            //update particles contact
            if(flag){
                UpdateParticlesContacts();
            }
            
            UpdateParticlePairForce();
        
            //move
            MoveParticles();
        }
        

        //collide and streaming
        (this->*ptr2collide)();
        Stream();
        BounceBack(false);
        CalcProps();

        Time += 1;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

inline void Domain::SolveRW(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    std::cout<<"dt of LBM "<<dt<<" dt of DEM "<<dtdem<<std::endl;
    double tout = 0;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    GhostParticles.assign(Particles.begin(), Particles.end()); 
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, idx_out);
            
            WriteXDMF(fn.CStr());
            idx_out++;
            
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            
            

            tout += dtout;
        }
        

        


        for(int i=0; i<std::floor(dt/dtdem); ++i)
        {
            // std::cout<<i<<std::endl;
            bool flag = i==0 || i==(std::floor(dt/dtdem)-1);
        
            if(flag){
                SetZero();
            }
            //set added force and check leave particles
            LeaveAndForcedForce();

            GhostParticles.clear();
            GhostParticles.assign(Particles.begin(), Particles.end()); 
        
            GhostPeriodic();
        
         
            //set fluid force
            if(flag){
                AddDisksG();
            }

            //update particles contact
            if(flag){
                UpdateParticlesContacts();
            }
            
            UpdateParticlePairForce();
        
            //move
            MoveParticles();
        }
        

        //collide and streaming
        (this->*ptr2collide)();
        Stream();
        BounceBack(false);
        CalcProps();

        //trace particle
        rwsolve_sub(dt);

        Time += 1;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

inline void Domain::rwsolve_sub(double dt)
{
    // std::cout<<1<<std::endl;
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int i=0;i<RWParticles.size();++i)
    {
        RW::Particle *RWP = &RWParticles[i]; 
        int x1 = std::floor(RWP->X(0));
        int y1 = std::floor(RWP->X(1));
        int x2 = x1+1;
        int y2 = y1+1;
        if(x1<0 || y1<0 || x2>Ndim(0)-1 || y2>Ndim(1)-1)
        {
            RWP->Leave(modexy,Box);
            x1 = std::floor(RWP->X(0));
            y1 = std::floor(RWP->X(1));
            x2 = x1+1;
            y2 = y1+1;
        }
        std::vector<Vec3_t> VV{Vel[x1][y1][0],Vel[x2][y1][0],Vel[x1][y2][0],Vel[x2][y2][0]};
        std::vector<int> idx{x1,x2,y1,y2};
        RWP->Move(VV,idx,dt);
        RWP->Leave(modexy,Box);
        int ix = std::round(RWP->X(0));
        int iy = std::round(RWP->X(1));
        if(Gamma[ix][iy][0]>1e-9 && Check[ix][iy][0]>0) 
        {
            int ip = Check[ix][iy][0];
            if(Norm(Particles[ip].X-RWP->X)<=Particles[ip].Rh && Norm(Particles[ip].Xb-RWP->Xb)>Particles[ip].Rh)
            {
                RWP->Reflect(Particles[ip].X,Particles[ip].Rh);
                RWP->Leave(modexy,Box);

            }
        }
    }

}

inline void Domain::CheckInside()
{
    for(size_t i=0; i<RWParticles.size(); ++i)
    {
        RW::Particle *RWP = &RWParticles[i]; 
        for(size_t ip=0; ip<Particles.size(); ++ip)
        {
            DEM::Disk *Pa = &Particles[ip];
            if(Norm(Pa->X-RWP->X)<Pa->Rh)
            {
                RWParticles.erase(RWParticles.begin()+i);
                break;
            }
        }

    }
}

#endif