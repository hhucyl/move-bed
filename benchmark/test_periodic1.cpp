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
        fs.Printf("%s.out","periodic1_6");
        // fs.Printf("%s_%d_%d_%g.out","Permeability",dat.bbtype,nx,dom.Tau);
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"Y"<<Util::_8s<<"Fx"<<Util::_8s<<"Fy\n";
    }else{
        
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<dom.Particles[0].X(1)<<Util::_8s<<dom.Particles[0].Fh(0)<<Util::_8s<<dom.Particles[0].Fh(1)<<std::endl;
        
    }
}


  
int main (int argc, char **argv) try
{
    
    
    size_t Nproc = 8;
    size_t h = 100;
    double nu = 0.002;
    if(argc>=2) Nproc = atoi(argv[1]); 

    size_t nx = h;
    size_t ny = h;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double R = 5;
    std::cout<<"R = "<<R<<std::endl;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = 2e-7;
    my_dat.R = R;
    Vec3_t g0(2e-7,0.0,0.0);
    dom.Nproc = Nproc;       

    
   
    //initial
    double rho = 1.0;
    double rhos = 1.0;
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    
    Vec3_t pos(nx/2,ny-1-40,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    //DEM
    dom.dtdem = 1.0*dt;
    dom.Particles.push_back(DEM::Disk(-1, pos, v, w, rhos, R, dom.dtdem));
    

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff = 0.0 , 0.0 , 0.0;
        dom.Particles[ip].Kn = 100;
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
    // for(size_t iy=0; iy<ny; iy++)
    // {
    //     dom.IsSolid[0][iy][0] = true;
    //     dom.IsSolid[nx-1][iy][0] = true;
    // }

    dom.Initial(rho,v0,g0);
    dom.IsF = true;
    // dom.InitialFromH5("test_periodic1_0001.h5",g0);
    double Tf = 1e6;
    
    double dtout = 1e4;
    dom.Box = 0.0, (double) nx-1, 0.0;
    dom.modexy = 0;
    
    //solving
    dom.SolveIBM( Tf, dtout, "test_periodic2", NULL, Report);
    return 0;

}MECHSYS_CATCH