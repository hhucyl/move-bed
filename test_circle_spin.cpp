#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double w;
    double vb;
    double nu;
    double R;
    double rhos;
};

void inout(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(0);
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t iy=0; iy<ny; ++iy)
    {
        double *f = dom.F[0][iy][0];
        double *f1 = dom.F[1][iy][0];
        double yy = (double)iy+0.5;
        double vx = -dat.vb*4.0/(((double) ny)*((double) ny))*yy*(yy-(double)ny);
        Vec3_t vel(vx,0.0,0.0);
        double rho = dom.Rho[1][iy][0];
        Vec3_t vel1 = dom.Vel[1][iy][0];
        for(size_t k=0;k<dom.Nneigh; ++k)
        {
            f[k] = dom.Feq(k,rho,vel) + f1[k] - dom.Feq(k,rho,vel1); 
        }
        Vec3_t idx(0,iy,0);
        dom.CalcPropsForCell(idx);
    }

    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t iy=0; iy<ny; ++iy)
    {
        double *f = dom.F[nx-1][iy][0];
        double *f1 = dom.F[nx-2][iy][0];
        for(size_t k=0;k<dom.Nneigh; ++k)
        {
            f[k] =  f1[k] ; 
        }
        Vec3_t idx(nx-1,iy,0);
        dom.CalcPropsForCell(idx);
    }
    

}

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
        fs.Printf("%s.out","spinning");
        // fs.Printf("%s_%d_%d_%g.out","Permeability",dat.bbtype,nx,dom.Tau);
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"w"<<Util::_8s<<"R"<<Util::_8s<<"Spa"<<Util::_8s<<"Cl\n";
    }else{
        // double Cd = 8*dom.Particles[0].R*dat.g*(dat.rhos-1.0)/(3.0*dom.Particles[0].V(1)*dom.Particles[0].V(1));
        double F = 0;
        for(size_t ix=0; ix<nx; ix++)
        for(size_t iy=0; iy<ny; iy++)
        {
            F+=dom.Flbm[ix][iy][0](1);
        }
        double spa = dat.w*dom.Particles[0].R/dat.vb;
        double Cl = F/(dat.vb*dat.vb*dom.Particles[0].R);
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<dom.Particles[0].W(2)<<Util::_8s<<dom.Particles[0].R<<Util::_8s<<spa<<Util::_8s<<Cl<<std::endl;
        
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



//something to make compare simple
  
int main (int argc, char **argv) try
{
    
    
    size_t Nproc = 1;
    size_t h = 100;
    double nu = 0.1;
    int N=1;
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
    my_dat.w = 2e-4;
    my_dat.vb = 1e-3;
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
    
    Vec3_t pos(nx*0.5,ny*0.5,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,my_dat.w);
    for(size_t ip=0; ip<N; ++ip)
    {

        dom.Particles.push_back(LBM::Disk(-ip, pos, v, w, rhos, R, dt));
    }
    dom.Particles[0].FixVeloc();
    std::cout<<dom.Particles[0].W(2)<<std::endl;
    for(size_t ix=0;ix<nx; ++ix)
    {
        dom.IsSolid[ix][0][0] = true;
        dom.IsSolid[ix][ny-1][0] = true;
    }

    double Tf = 1e3;
    
    double dtout = 1;
    char const * TheFileKey = "test_spin";
    //solving
    dom.StartSolve();
    double tout = 0;
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
        inout(dom,&my_dat);
        std::cout<<1<<std::endl;
        //collide and streaming
        dom.AddDisksG();
                
        (dom.*dom.ptr2collide)();

        //set added force
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for(size_t i=0;i<dom.Particles.size();i++)
        {
            dom.Particles[i].F = dom.Particles[i].Ff;
            dom.Particles[i].T = dom.Particles[i].Tf;
        }
        //set fluid force 
        
        //move
        #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
        for (size_t i=0; i<N; i++)
        {
		    dom.Particles[i].Translate(dt);
		    dom.Particles[i].Rotate(dt);
        }

        dom.Stream();
        dom.BounceBack(false);
        dom.CalcProps();
        dom.Time += 1;
    }
    dom.EndSolve();
    return 0;
}MECHSYS_CATCH