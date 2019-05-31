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




double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}
  
int main (int argc, char **argv) try
{
    
    std::srand((unsigned)time(NULL));    
    size_t Nproc = 12;
    size_t h = 3200;
    double nu = 0.01;
    int ratio = 2;
    bool initfromfile = false;
    char *h5name = NULL;
    if(argc>=2) Nproc = atoi(argv[1]);
    if(argc>=3) ratio = atoi(argv[2]);
    if(argc>=4)
    {
        initfromfile = true;
        h5name = argv[3];
    }

    size_t nx = h;
    size_t ny = h/2;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double R = 10;
    double Ga = 20.0;
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
    my_dat.g = -2e-9*ratio;
    my_dat.R = R;
    Vec3_t g0(my_dat.g,0.0,0.0);
    std::cout<<"gx = "<<my_dat.g<<std::endl;
    dom.Nproc = Nproc;       

    //initial
    
    my_dat.rhos = rhos;
    
    
    Vec3_t pos(R,R,0.0);
    Vec3_t dxp(0.0,0.0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    //DEM
    dom.dtdem = 0.01*dt;
    int pnum = 0;
    //fixed
    for(size_t ip=0; ip<160; ++ip)
    {
        // std::cout<<pos<<std::endl;
        dom.Particles.push_back(DEM::Disk(pnum, pos, v, w, rhos, R, dom.dtdem));
        dom.Particles[ip].FixVeloc();
        dxp = 2.0*R,0.0,0.0;
        pos = pos+dxp;
        pnum++;
    }
    //move
    for(int ipy=0; ipy<40; ++ipy)
    {
        pos = R,(2*ipy+2)*R+R,0.0;
        for(int ipx=0; ipx<160; ++ipx)
        {
            Vec3_t dxr(random(-0.1,0.1),random(-0.1,0.1));
            dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, rhos, R, dom.dtdem));
            // std::cout<<pos(0)<<" "<<pos(1)<<std::endl;
            dxp = 2.0*R,0.0,0.0;
            pos = pos+dxp;
            pnum++;
        }   
    }

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    
    for(int ip=0; ip<(int) dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff = 0.0, -M_PI*R*R*(rhos/rho-1)*gy, 0.0;
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

    Vec3_t v0(0.0,0.0,0.0);
    dom.IsF = true;
  
    dom.InitialFromH5("test_mvbed_2_0523.h5",g0);
    

        //dom.Initial(rho,v0,g0);
    for(int i=0;i<nx-1;++i)
    {
        Vec3_t xt(0.5+i,189.5,0);
        int ix = std::round(xt(0));
        int iy = std::round(xt(1));
        // std::cout<<ix<<" "<<iy<<" "<<dom.Gamma[ix][iy][0]<<std::endl;
        if(dom.Gamma[ix][iy][0]<1e-9)
        {
            dom.RWParticles.push_back(RW::Particle(xt));
        }
    }
    // std::cout<<"Checking Inside"<<std::endl;
    // dom.CheckInside();

    double Tf = 1e6;
    
    double dtout = 1e3;
    dom.Box = 0.0,(double) nx-1, 0.0;
    dom.modexy = 0;
    //solving
    dom.SolveRW( Tf, dtout, "test_mvbedrw_1", NULL, NULL);
    
    return 0;
}MECHSYS_CATCH
