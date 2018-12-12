#include "./lbm/Domain.h"
#include <time.h>

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
        fs.Printf("%s.out","mvbed");
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"Collision"<<Util::_8s<<"GhostNum\n";
    }else{
        
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<dom.ListofContacts.size()<<Util::_8s<<dom.GhostParticles[0].Ghost<<std::endl;
        
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

double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}
  
int main (int argc, char **argv) try
{
    
    std::srand((unsigned)time(NULL));    
    size_t Nproc = 12;
    size_t h = 400;
    double nu = 0.01;
    if(argc>=2) Nproc = atoi(argv[1]); 

    size_t nx = h+10;
    size_t ny = 2*h;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double R = 20;
    double Ga = 2.4;
    double rho = 1.0;
    double rhos = 2.0;
    double gy = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
    std::cout<<"gy = "<<gy<<std::endl;

    std::cout<<"R = "<<R<<std::endl;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = -2e-9;
    my_dat.R = R;
    Vec3_t g0(my_dat.g,0.0,0.0);
    dom.Nproc = Nproc;       

    //dom.Isq = true;
    // dom.IsF = false;
    // dom.IsFt = false;
   
    //initial
    
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    Initial(dom,rho,v0,g0);
    
    Vec3_t pos(R+1,R,0.0);
    Vec3_t dxp(0.0,0.0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    //DEM
    dom.dtdem = 0.1*dt;
    //fixed
    for(size_t ip=0; ip<10; ++ip)
    {
        // std::cout<<pos<<std::endl;
        dom.Particles.push_back(DEM::Disk(-ip, pos, v, w, rhos, R, dom.dtdem));
        dom.Particles[ip].FixVeloc();
        dxp = 2.0*R+1,0.0,0.0;
        pos = pos+dxp;
    }
    //move
    int num = 0;
    for(int ipy=0; ipy<9; ++ipy)
    {
        pos = R+1,(2*ipy+2)*R+R+1,0.0;
        for(int ipx=0; ipx<10; ++ipx)
        {
            Vec3_t dxr(random(-0.1,0.1),random(-0.1,0.1));
            dom.Particles.push_back(DEM::Disk(-num, pos+dxr, v, w, rhos, R, dom.dtdem));
            // std::cout<<pos(0)<<" "<<pos(1)<<std::endl;
            dxp = 2.0*R+1,0.0,0.0;
            pos = pos+dxp;
            num++;
        }   
    }

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(int ip=0; ip<dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff = 0.0, -M_PI*R*R*rhos*gy, 0.0;
        dom.Particles[ip].Kn = 5;
        dom.Particles[ip].Gn = 0.8;
        dom.Particles[ip].Kt = 2.5;
        // dom.Particles[ip].Mu = 0.0;
        // dom.Particles[ip].Eta = 0.0;
        // dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = 0.8*R;
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


    double Tf = 2e2;
    
    double dtout = 1e2;;
    dom.Box = 0, nx-1, 0;
    dom.modexy = 0;
    //solving
    dom.Solve( Tf, dtout, "test_mvbed", NULL, NULL);
    
    return 0;
}MECHSYS_CATCH