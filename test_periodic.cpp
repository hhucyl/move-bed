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
        fs.Printf("%s.out","periodic");
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
        
        // dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<F<<Util::_8s<<dom.Particles[0].F(1)-dom.Particles[0].Ff(1) <<Util::_8s<<dom.Particles[0].X(1)<<Util::_8s<<dom.Particles[0].F<<Util::_8s<<dom.ListofContacts.size()<<dom.GhostParticles[0].Ghost<<std::endl;
        
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
    size_t h = 200;
    double nu = 0.01;
    if(argc>=2) Nproc = atoi(argv[1]); 

    size_t nx = h+5;
    size_t ny = h;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double R = 20;
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
    dom.IsF = false;
    // dom.IsFt = false;
   
    //initial
    double rho = 1.0;
    double rhos = 2.0;
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    
    Vec3_t pos(nx-1-(2*R+0.1),0.5*ny+0.1,0.0);
    // Vec3_t pos(0.5*nx+0.1,ny-1-(R+0.1),0.0);
    // Vec3_t pos(nx*0.5,ny-1-3.0*R,0.0);
    // Vec3_t pos1(nx*0.5,2.5*R,0.0);
    Vec3_t dxp(2*R+1,0.0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t v1(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    //DEM
    for(size_t ip=0; ip<1; ++ip)
    {
        // std::cout<<pos<<std::endl;
    // dom.Particles.push_back(DEM::Disk(-1, pos, v, w, rhos, R, dt));
    // dom.Particles.push_back(DEM::Disk(-2, pos1, v1, w, rhos, R, dt));
    // dom.Particles[1].FixVeloc();
        dom.Particles.push_back(DEM::Disk(-ip, pos, v, w, rhos, R, dt));
        pos = pos+dxp;
        // pos1 = pos1+dxp;
    }

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(int ip=0; ip<dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff =  M_PI*R*R*rhos*my_dat.g, 0.0 , 0.0;
        // dom.Particles[ip].Ff =  0.0, -M_PI*R*R*rhos*my_dat.g , 0.0;
        dom.Particles[ip].Kn = 100;
        dom.Particles[ip].Gn = 0.0;
        dom.Particles[ip].Kt = 0.0;
        dom.Particles[ip].Mu = 0.0;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = R;
        // dom.Particles[ip].FixVeloc();

    }
    
    // for(size_t ix=0; ix<nx; ix++)
    // {
    //     dom.IsSolid[ix][0][0] = true;
    //     dom.IsSolid[ix][ny-1][0] = true;
    // }
    // for(size_t iy=0; iy<ny; iy++)
    // {
    //     dom.IsSolid[0][iy][0] = true;
    //     dom.IsSolid[nx-1][iy][0] = true;
    // }

    Initial(dom,rho,v0,g0);
    // dom.InitialFromH5("test_periodic1_0020.h5",g0);
    double Tf = 1e4;
    
    double dtout = 1e2;
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    dom.dtdem = dt;
    dom.IsF = true;
    //solving
    dom.SolveIBM( Tf, dtout, "test_periodic1", NULL, Report);    
    
    return 0;

}MECHSYS_CATCH