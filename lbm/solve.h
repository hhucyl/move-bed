#ifndef LBM_SOLVE_H
#define LBM_SOLVE_H

inline void Domain::Solve(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport)
{
    StartSolve();
    std::cout<<"Box "<<Box<<std::endl;
    std::cout<<"modexy "<<modexy<<std::endl;
    double tout = 0;
    while(Time<Tf)
    {
        if (Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, idx_out);
            
            WriteXDMF(fn.CStr());
            idx_out++;
            // std::cout<<"--- Time = "<<Time<<" "<<Tf<<" ---"<<std::endl;
            if (ptReport!=NULL) (*ptReport) ((*this), UserData); 
            

            tout += dtout;
        }
        
        
        for(int i=0; i<std::floor(dt/dtdem); ++i)
        {
            // std::cout<<i<<std::endl;
            SetZero();
            //set added force and check leave particles
            LeaveAndForcedForce();

            GhostParticles.clear();
            GhostParticles.assign(Particles.begin(), Particles.end()); 
        
            GhostPeriodic();
        
         
            //set fluid force
            AddDisksG();

            //update particles contact
            UpdateParticlesContacts();
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
    EndSolve();
}


#endif