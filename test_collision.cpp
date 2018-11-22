#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
    double rhos;
};

void Report(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    size_t nz = dom.Ndim(2);
    double dx = dom.dx;
    if(dom.Time <1e-6)
    {
        String fs;
        fs.Printf("%s.out","2");
        // fs.Printf("%s_%d_%d_%g.out","Permeability",dat.bbtype,nx,dom.Tau);
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"X1"<<Util::_8s<<"X2"<<Util::_8s<<"DX"<<Util::_8s<<"Collision\n";
    }else{
        // double Cd = 8*dom.Particles[0].R*dat.g*(dat.rhos-1.0)/(3.0*dom.Particles[0].V(1)*dom.Particles[0].V(1));
        double F = 0;
        for(size_t ix=0; ix<nx; ix++)
        for(size_t iy=0; iy<ny; iy++)
        {
            F+=dom.Flbm[ix][iy][0](1);
        }
        
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<dom.Particles[0].X(1)<<Util::_8s<<dom.Particles[1].X(1)<<Util::_8s<<dom.Particles[1].X(1)-dom.Particles[0].X(1)<<Util::_8s<<dom.ListofContacts.size()<<std::endl;
        
    }
}

void Initial(LBM::Domain &dom, double rho, Vec3_t &v0,  Vec3_t &g0)
{
    
    for(size_t ix=0; ix<dom.Ndim(0); ix++)
    for(size_t iy=0; iy<dom.Ndim(1); iy++)
    for(size_t iz=0; iz<dom.Ndim(2); iz++)
    {
        dom.Rho[ix][iy][iz] = rho;
        dom.Vel[ix][iy][iz] = 0.0, 0.0, 0.0;
        dom.BForce[ix][iy][iz] = g0;
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            dom.F[ix][iy][iz][k] = dom.Feq(k,rho,v0);            
            dom.Ftemp[ix][iy][iz][k] = dom.Feq(k,rho,v0);            
        }
    // std::cout<<dom.F[ix][iy][iz][18]<<std::endl;
        
    }
    dom.Rho0 = rho;//very important
}


  
int main (int argc, char **argv) try
{
    
    
    size_t Nproc = 8;
    size_t h = 400;
    double nu = 0.01;
    int N=8;
    if(argc>=2) nu = atof(argv[1]);     
    if(argc>=3) N = atoi(argv[2]);
    if(argc>=4) Nproc = atoi(argv[3]); 

    size_t nx = h;
    size_t ny = h;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double R = 10;
    std::cout<<"R = "<<R<<std::endl;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = 2e-5;
    my_dat.R = R;
    Vec3_t g0(0.0,0.0,0.0);
    dom.Nproc = Nproc;       

    //dom.Isq = true;
    // dom.IsF = false;
    // dom.IsFt = false;
   
    //initial
    double rho = 1.0;
    double rhos = 2.0;
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    Initial(dom,rho,v0,g0);
    
    Vec3_t pos(nx*0.1+10.0*R,0.4*ny,0.0);
    Vec3_t pos1(nx*0.1+10.0*R+1,0.6*ny,0.0);
    Vec3_t dxp(2.5*R,0.0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    for(size_t ip=0; ip<N/2; ++ip)
    {
        // std::cout<<pos<<std::endl;
        dom.Particles.push_back(DEM::Disk(-ip, pos, v, w, rhos, R, dt));
        dom.Particles.push_back(DEM::Disk(-ip, pos1, v, w, rhos, R, dt));
        pos = pos+dxp;
        pos1 = pos1+dxp;
    }
    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(int ip=0; ip<N; ++ip)
    {
        dom.Particles[ip].Ff = 0.0, std::pow(-1,ip)*M_PI*R*R*rhos*my_dat.g, 0.0;
        dom.Particles[ip].Kn = 10;
        dom.Particles[ip].Gn = 0.0;
        dom.Particles[ip].Kt = 0.0;
        dom.Particles[ip].Mu = 0.0;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = R;
        // dom.Particles[ip].FixVeloc();

    }
    for(size_t ix=0; ix<nx; ix++)
    {
        dom.IsSolid[ix][0][0] = true;
        dom.IsSolid[ix][ny-1][0] = true;
    }
    for(size_t iy=0; iy<ny; iy++)
    {
        dom.IsSolid[0][iy][0] = true;
        dom.IsSolid[nx-1][iy][0] = true;
    }


    double Tf = 1e4;
    
    double dtout = 1e2;
    char const * TheFileKey = "test_collision";
    //solving
    dom.StartSolve();
    double tout = 0;
    dom.dtdem = dt;
    while(dom.Time<Tf)
    {
        if (dom.Time>=tout)
        {
            
            String fn;
            fn.Printf("%s_%04d", TheFileKey, dom.idx_out);
            
            dom.WriteXDMF(fn.CStr());
            dom.idx_out++;
            // std::cout<<"--- Time = "<<dom.Time<<" "<<Tf<<" ---"<<std::endl;
            Report(dom,&my_dat); 
            tout += dtout;
        }
        dom.SetZero();
        //set added force
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for(size_t i=0;i<dom.Particles.size();i++)
        {
            dom.Particles[i].F = dom.Particles[i].Ff;
            dom.Particles[i].T = dom.Particles[i].Tf;
        }
        //set fluid force 
        dom.AddDisksG();
        dom.UpdateParticlesContacts();
        #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
        for(size_t i=0; i<dom.ListofContacts.size();++i)
        {
            int ip1 = dom.ListofContacts[i].first;
            int ip2 = dom.ListofContacts[i].second;
            DEM::DiskPair pair(&dom.Particles[ip1],&dom.Particles[ip2]);
            pair.CalcForce(dom.dtdem);
            omp_set_lock  (&pair.P1->lck);
            pair.P1->F += pair.F1;
            pair.P1->T += pair.T1;
            omp_unset_lock(&pair.P1->lck);
            omp_set_lock  (&pair.P2->lck);
            pair.P2->F += pair.F2;
            pair.P2->T += pair.T2;
            omp_unset_lock(&pair.P2->lck);
        }
        //move
        #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
        for (size_t i=0; i<N; i++)
        {
		    dom.Particles[i].Translate(dt);
		    dom.Particles[i].Rotate(dt);
        }
        //collide and streaming
        (dom.*dom.ptr2collide)();
        dom.Stream();
        dom.BounceBack(false);
        dom.CalcProps();
        dom.Time += 1;
    }
    dom.EndSolve();
    return 0;
}MECHSYS_CATCH