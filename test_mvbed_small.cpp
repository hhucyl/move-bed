#include "./lbm/Domain.h"
#include <time.h>

struct myUserData
{
    std::ofstream oss_ss;
    double vb;
    double nu;
    double R;
    double rhos;
};

void Report(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    // size_t nx = dom.Ndim(0);
    // size_t ny = dom.Ndim(1);
    // size_t nz = dom.Ndim(2);
    // double dx = dom.dx;
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

void Setup(LBM::Domain &dom, void *UD)
{
    size_t nx = dom.Ndim(0);
    size_t ny = dom.Ndim(1);
    myUserData &dat = (*static_cast<myUserData *> (UD));
    Vec3_t vvb(dat.vb,0.0,0.0);
    // std::cout<<vvb<<std::endl;
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for(size_t ix=0;ix<nx;++ix)
    {
        double *f = dom.F[ix][ny-1][0];
        double *f1 = dom.F[ix][ny-2][0];
        double rho1 = dom.Rho[ix][ny-2][0];
        Vec3_t vel1 = dom.Vel[ix][ny-2][0];
        
        for(size_t k=0; k<dom.Nneigh; ++k)
        {
            f[k] = dom.Feq(k,rho1,vvb) + f1[k] - dom.Feq(k,rho1,vel1);
        }
        Vec3_t idx(ix,ny-1,0);
        dom.CalcPropsForCell(idx);
    }

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
    int Nx = 40;
    int Ny = 40;
    int Nduny = 10;
    int Ndunx = 30;
    size_t Rn = 10;
    double gap = 0.3;
    double nu = 0.01;
    double ratio = 0.5;
    double vb = 0.01;
    bool initfromfile = false;
    char *h5name = NULL;
    if(argc>=2) Nproc = atoi(argv[1]);
    if(argc>=3) ratio = atof(argv[2]);
    if(argc>=4)
    {
        initfromfile = true;
        h5name = argv[3];
    }
    int gapn = std::ceil(gap*Nx);
    std::cout<<"extra gap n "<<gapn<<std::endl;
    size_t nx = 2*Rn*Nx+gapn;
    size_t ny = 2*Rn*(Ny+1)+4*Nduny*Rn;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double R = (double) Rn;
    double Ga = 10.0;
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
    my_dat.vb = vb;
    my_dat.R = R;
    Vec3_t g0(0.0,0.0,0.0);
    std::cout<<"vb = "<<my_dat.vb<<std::endl;
    dom.Nproc = Nproc;       

    //initial
    
    my_dat.rhos = rhos;
    
    
    Vec3_t pos(R+gap,R,0.0);
    Vec3_t dxp(0.0,0.0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    //DEM
    dom.dtdem = 1*dt;
    int pnum = 0;
    //fixed
    for(int ip=0; ip<Nx; ++ip)
    {
        // std::cout<<pos<<std::endl;
        dom.Particles.push_back(DEM::Disk(pnum, pos, v, w, rhos, R, dom.dtdem));
        dom.Particles[ip].FixVeloc();
        dxp = 2.0*R+gap,0.0,0.0;
        pos = pos+dxp;
        pnum++;
    }
    //move
    for(int ipy=0; ipy<Ny; ++ipy)
    {
        pos = R+gap,(2*ipy+3)*R+gap,0.0;
        for(int ipx=0; ipx<Nx; ++ipx)
        {
            Vec3_t dxr(random(-0.1,0.1),random(-0.1,0.1),0.0);
            // Vec3_t dxr(0.0,0.0,0.0);
            dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, rhos, R, dom.dtdem));
            // std::cout<<pos(0)<<" "<<pos(1)<<std::endl;
            

            dxp = 2.0*R+gap,0.0,0.0;
            pos = pos+dxp;
            pnum++;
        }   
    }
    Vec3_t posd1(R+gap,(2*Ny+3)*R+gap,0.0);
    Vec3_t posd(R+2*R*(Ndunx-1)+gap*Ndunx,(2*Ny+3)*R+2*R*(Nduny-1)+gap*Nduny,0.0);
    Vec3_t posd2(R+2*R*(Nx-1)+gap*Nx,(2*Ny+3)*R+gap,0.0);
    double k1 = (posd(1)-posd1(1))/(posd(0)-posd1(0));
    double k2 = (posd(1)-posd2(1))/(posd(0)-posd2(0));
    //dune
    for(int idy=0; idy<Nduny; ++idy)
    {
        double yt = posd1(1)+2*idy*R+idy*gap;
        Vec3_t posts((yt-posd1(1))/k1+posd1(0),yt,0.0);
        Vec3_t poste((yt-posd2(1))/k2+posd2(0),yt,0.0);
        int Nt = std::ceil((poste(0)-posts(0))/(2*R));
        for (int idx=0; idx<Nt; ++idx)
        {
            Vec3_t post(posts(0)+2*R*idx+gap*idx,yt,0.0);
            Vec3_t dxr(random(-0.1,0.1),random(-0.1,0.1),0.0);
            // Vec3_t dxr(0.0,0.0,0.0);
            dom.Particles.push_back(DEM::Disk(-pnum, post+dxr, v, w, rhos, R, dom.dtdem));
            
            pnum++;
        }
        if(Nt==0)
        {
            Vec3_t dxr(random(-0.1,0.1),random(-0.1,0.1),0.0);

            // Vec3_t dxr(0.0,0.0,0.0);
            dom.Particles.push_back(DEM::Disk(-pnum, posd, v, w, rhos, R, dom.dtdem));

        }
        
    }

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    
    for(int ip=0; ip<(int) dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff = 0.0, -M_PI*R*R*rhos*gy, 0.0;
        dom.Particles[ip].Kn = 5;
        dom.Particles[ip].Gn = 0.8;
        dom.Particles[ip].Kt = 2.5;
        // dom.Particles[ip].Mu = 0.0;
        // dom.Particles[ip].Eta = 0.0;
        // dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = 0.8*R;
        dom.Particles[ip].FixVeloc();
        // dom.Particles[ip].FixVeloc();

    }
    
    for(size_t ix=0; ix<nx; ix++)
    {
        dom.IsSolid[ix][0][0] = true;
        // dom.IsSolid[ix][ny-1][0] = true;
    }
    // for(size_t iy=0; iy<ny; iy++)
    // {
    //     dom.IsSolid[0][iy][0] = true;
    //     dom.IsSolid[nx-1][iy][0] = true;
    // }

    Vec3_t v0(0.0,0.0,0.0);
    dom.IsF = true;
    if(initfromfile)
    {
        dom.InitialFromH5(h5name,g0);

    }else{
        dom.Initial(rho,v0,g0);

    }


    double Tf = 1e6;
    
    double dtout = 1e3;
    dom.Box = 0.0,(double) nx-1, 0.0;
    dom.modexy = 0;
    //solving
    dom.SolveIBM( Tf, dtout, "test_mvbed_c", Setup, NULL);
    
    return 0;
}MECHSYS_CATCH