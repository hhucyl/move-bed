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
    // size_t nz = dom.Ndim(2);
    // double dx = dom.dx;
    if(dom.Time <1e-6)
    {
        String fs;
        fs.Printf("%s.out","periodic2");
        // fs.Printf("%s_%d_%d_%g.out","Permeability",dat.bbtype,nx,dom.Tau);
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"X1"<<Util::_8s<<"X2"<<Util::_8s<<"DX"<<Util::_8s<<"Collision\n";
    }else{
        // double Cd = 8*dom.Particles[0].R*dat.g*(dat.rhos-1.0)/(3.0*dom.Particles[0].V(1)*dom.Particles[0].V(1));
        double F = 0;
        for(size_t ix=0; ix<nx; ix++)
        for(size_t iy=0; iy<ny; iy++)
        {
            F+=dom.Flbm[ix][iy][0](0);
        }
        
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<dom.Particles[0].F(0)-dom.Particles[0].Ff(0) <<Util::_8s<<dom.Particles[0].X(0)<<Util::_8s<<dom.Particles[0].F<<Util::_8s<<dom.ListofContacts.size()<<dom.GhostParticles[0].Ghost<<std::endl;
        
    }
}


  
int main (int argc, char **argv) try
{
    
    
    size_t Nproc = 8;
    size_t h = 400;
    double nu = 0.01;
    if(argc>=2) Nproc = atoi(argv[1]); 

    size_t nx = h;
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
    my_dat.g = 2e-7;
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
    dom.Initial(rho,v0,g0);
    
    Vec3_t pos(R+5,ny/2,0.0);
    // Vec3_t pos(nx*0.5,ny-1-3.0*R,0.0);
    // Vec3_t pos1(nx*0.5,2.5*R,0.0);
    Vec3_t dxp(2*R,0.0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t v1(0.1,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    //DEM
    dom.dtdem = 0.1*dt;
    dom.Particles.push_back(DEM::Disk(-1, pos, v, w, rhos, R, dom.dtdem));
    // dom.Particles[0].FixVeloc();
    pos(0) = nx/2;
    dom.Particles.push_back(DEM::Disk(-2, pos, v1, w, rhos, R, dom.dtdem));
    pos(0) = nx-1-R-5;
    dom.Particles.push_back(DEM::Disk(-3, pos, v, w, rhos, R, dom.dtdem));

    pos = R+5, ny/2 - 3*R, 0.0;
    dom.Particles.push_back(DEM::Disk(-4, pos, v, w, rhos, R, dom.dtdem));
    // dom.Particles[0].FixVeloc();
    pos(0) = nx/2;
    dom.Particles.push_back(DEM::Disk(-5, pos, v1, w, rhos, R, dom.dtdem));
    pos(0) = nx-1-R-5;
    dom.Particles.push_back(DEM::Disk(-6, pos, v, w, rhos, R, dom.dtdem));

    pos = R+5, ny/2 + 3*R, 0.0;
    dom.Particles.push_back(DEM::Disk(-7, pos, v, w, rhos, R, dom.dtdem));
    // dom.Particles[0].FixVeloc();
    pos(0) = nx/2;
    dom.Particles.push_back(DEM::Disk(-8, pos, v1, w, rhos, R, dom.dtdem));
    pos(0) = nx-1-R-5;
    dom.Particles.push_back(DEM::Disk(-9, pos, v, w, rhos, R, dom.dtdem));
    

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        // dom.Particles[ip].Ff = M_PI*R*R*rhos*my_dat.g, 0.0 , 0.0;
        dom.Particles[ip].Kn = 10;
        dom.Particles[ip].Gn = 0.0;
        dom.Particles[ip].Kt = 0.0;
        dom.Particles[ip].Mu = 0.0;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = 0.8*R;
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


    double Tf = 1e4;
    
    double dtout = 10;
    dom.Box = 0.0, (double) ny-1, 0.0;
    dom.modexy = 0;
    //solving
    std::cout<<std::floor(dom.dt/dom.dtdem)<<std::endl;

    dom.Solve( Tf, dtout, "test_periodic1", NULL, Report);
    return 0;

}MECHSYS_CATCH