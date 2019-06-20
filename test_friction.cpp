#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double x;
    double R;
    double rhos;
};

void Report(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    if(dom.Time <1e-6)
    {
        String fs;
        fs.Printf("%s.out","2");
        // fs.Printf("%s_%d_%d_%g.out","Permeability",dat.bbtype,nx,dom.Tau);
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"X\n";
    }else{
        double ddd = dom.Particles[1].X(0)-dat.x;
        if(ddd<0) ddd = dom.Particles[1].X(0) + dom.Ndim(0)-dat.x;
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<ddd<<std::endl;
        
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
    if(argc>=2) nu = atof(argv[1]);

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
    my_dat.x = 2.0*R;
    my_dat.g = -0.5e-7;
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
    
    
    Vec3_t pos(my_dat.x,0.4*ny,0.0);
    // Vec3_t pos1(nx*0.1+10.0*R,0.4*ny+4*R+0.5*R,0.0);
    Vec3_t pos1(my_dat.x,0.4*ny + 10,0.0);
    Vec3_t pos2(my_dat.x,0.4*ny + 30,0.0);
    Vec3_t pos3(my_dat.x,0.4*ny + 40,0.0);
    Vec3_t dxp(2.5*R,0.0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t v1(1e-3,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    dom.dtdem = 0.01*dt;

        // std::cout<<pos<<std::endl;
    dom.Particles.push_back(DEM::Disk(-1, pos, v, w, rhos, R, dom.dtdem));
    dom.Particles.push_back(DEM::Disk(-2, pos1, v1, w, rhos, R, dom.dtdem));
    
    dom.Particles.push_back(DEM::Disk(-3, pos2, v, w, rhos, R, dom.dtdem));
    dom.Particles.push_back(DEM::Disk(-4, pos3, v1, w, rhos, R, dom.dtdem));

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        //dom.Particles[ip].Ff = 0.0, std::pow(-1,ip)*M_PI*R*R*rhos*my_dat.g, 0.0;
        dom.Particles[ip].Ff = 0.0, M_PI*R*R*rhos*my_dat.g, 0.0;
        dom.Particles[ip].Kn = 5;
        dom.Particles[ip].Gn = 0.0;
        dom.Particles[ip].Kt = 2.5;
        dom.Particles[ip].Mu = 0.4;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = R;
        dom.Particles[ip].nu = nu;
        
        // dom.Particles[ip].FixVeloc();

    }
    dom.Particles[0].FixVeloc();
    dom.Particles[2].FixVeloc();
    for(size_t ix=0; ix<nx; ix++)
    {
        dom.IsSolid[ix][0][0] = true;
        dom.IsSolid[ix][ny-1][0] = true;
    }
    
    Initial(dom,rho,v0,g0);
    // dom.InitialFromH5("test_2_1_0063.h5",g0);

    double Tf = 5e4;
    
    double dtout = 1e2;
    dom.IsF = true;
    dom.Box = 0.0, nx-1, 0.0;
    dom.modexy = 0;
    //solving
    dom.SolveIBM( Tf, dtout, "test_2", NULL, Report);
    
    
    
    return 0;
}MECHSYS_CATCH
