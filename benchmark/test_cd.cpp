#include "./lbm/Domain.h"


struct myUserData
{
    std::ofstream oss_ss;
    double w;
    double nu;
    double R;
    double rhos;
    double U;
};

void Report(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    
    if(dom.Time <1e-6)
    {
        String fs;
        fs.Printf("%s.out","drag");
        // fs.Printf("%s_%d_%d_%g.out","Permeability",dat.bbtype,nx,dom.Tau);
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"Cd"<<Util::_8s<<"Cl\n";
    }else{
        
        double Fd = dom.Particles[0].Fh(0);
        double Fl = dom.Particles[0].Fh(1);
        double Cd = Fd/(dat.U*dat.U*dat.R);
        double Cl = Fl/(dat.U*dat.U*dat.R);
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<Cd<<Util::_8s<<Cl<<std::endl;
        
    }
}

void Setup(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    size_t ny = dom.Ndim(1);
    size_t nx = dom.Ndim(0);
    
    for(size_t iy=0; iy<ny; ++iy)
    {
        size_t index = 0;
        Vec3_t idx(index,iy,0);
        double *f = dom.F[index][iy][0];
        double *f1= dom.F[1][iy][0];
        double rho1 = dom.Rho[1][iy][0];
        Vec3_t vel1 = dom.Vel[1][iy][0];
        double L = (double) (ny-1);
        double yy = (double) iy-0.5;
        // Vec3_t vel(-4.0*dat.U/(L*L)*yy*(yy-L),0.0,0.0);
        Vec3_t vel(dat.U,0.0,0.0);

        for(size_t k=0;k<dom.Nneigh; ++k)
        {
            f[k] = dom.Feq(k,rho1,vel) + f1[k] - dom.Feq(k,rho1,vel1);
        }
        dom.CalcPropsForCell(idx);

        index = nx-1;
        idx = index,iy,0;
        f = dom.F[index][iy][0];
        f1 = dom.F[index-1][iy][0];
        for(size_t k=0;k<dom.Nneigh; ++k)
        {
            f[k] = f1[k];
        }
        dom.CalcPropsForCell(idx);

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
    
    
    size_t Nproc = 8;
    if(argc>=2) Nproc = atoi(argv[1]); 

    size_t nx = 800;
    size_t ny = 800;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double R = 10.0;
    double Re = 20.0;
    double U = 0.03;
    double nu = 2.0*R*(2.0/3.0*U)/Re;
    std::cout<<"nu =  "<<nu<<std::endl;
    std::cout<<"R = "<<R<<std::endl;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.U = U;
    my_dat.R = R;
    Vec3_t g0(0.0,0.0,0.0);
    dom.Nproc = Nproc;       

    
   
    //initial
    double rho = 1.0;
    double rhos = 2.0;
    my_dat.rhos = rhos;
    Vec3_t v0(0.0,0.0,0.0);
    
    // Vec3_t pos(300.0,120.0,0.0);
    Vec3_t pos(0.5*nx,0.5*ny,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    //DEM
    dom.dtdem = 1.0*dt;
    dom.Particles.push_back(DEM::Disk(-1, pos, v, w, rhos, R, dom.dtdem));
    

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    for(size_t ip=0; ip<dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff = 0.0, 0.0 , 0.0;
        // dom.Particles[ip].Ff = 0.0, 0.0 , 0.0;
        dom.Particles[ip].Kn = 100;
        dom.Particles[ip].Gn = 0.0;
        dom.Particles[ip].Kt = 0.0;
        dom.Particles[ip].Mu = 0.0;
        dom.Particles[ip].Eta = 0.0;
        dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = 1.0*R;
        dom.Particles[ip].FixVeloc();

    }
    // for(size_t ix=0; ix<nx; ix++)
    // {
    //     dom.IsSolid[ix][0][0] = true;
    //     dom.IsSolid[ix][ny-1][0] = true;
    // }
    

    dom.Initial(rho,v0,g0);
    double Tf = 1e6;
    
    double dtout = 1e3;
    dom.Box = 0, nx-1, 0;
    dom.modexy = 0;
    //solving
    dom.SolveIBM( Tf, dtout, "test_cd", Setup, Report);
    
    return 0;
}MECHSYS_CATCH