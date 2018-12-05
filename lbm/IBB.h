#ifndef LBM_IBB_H
#define LBM_IBB_H

#include <mechsys/dem/special_functions.h>


inline void Domain::AddDiskQ(Vec3_t &pos, double R)
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    for(size_t iz=0; iz<nz; ++iz)
    for(size_t ix=0; ix<nx; ++ix)
    for(size_t iy=0; iy<ny; ++iy)
    {
        Vec3_t Cell(ix,iy,iz);
        double L = std::sqrt(dot(Cell - pos,Cell - pos));
        double Lx = Cell(0) - pos(0);
        double Ly = Cell(1) - pos(1);
        for(size_t k=0; k<Nneigh; ++k)
        {
            size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)nx)%nx;
            size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)ny)%ny;
            size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)nz)%nz;
            Vec3_t Cell1(nix,niy,niz);
            double L1 = std::sqrt(dot(Cell1 - pos,Cell1 - pos));
            double l = std::sqrt(dot(Cell1 - Cell,Cell1 - Cell));
            double lx = Cell1(0) - Cell(0);
            double ly = Cell1(1) - Cell(1);
            double delta = 4.0*((lx*Lx+ly*Ly)*(lx*Lx+ly*Ly) - (l)*(l)*(L*L-R*R));
            if(delta<0)
            {
                q[ix][iy][iz][k] = -1;
                continue;
            }
            double q1 = (-(2.0*lx*Lx+2.0*ly*Ly) + std::sqrt(delta))/(2.0*l*l);
            double q2 = (-(2.0*lx*Lx+2.0*ly*Ly) - std::sqrt(delta))/(2.0*l*l);
            bool flag1 = q1>=0 && q1-1< 1e-5;
            // double flag1 = q1>=0;
            bool flag2 = q2>=0 && q2-1< 1e-5;
            bool flag3 = L>R;
            bool flag4 = L1<R;
            // double flag2 = q2>=0;
            if(flag3)
            {
                if(flag4)
                {
                    // std::cout<<"!!! "<<q1<<" "<<q2<<std::endl;
                    if(flag1)
                    {
                        q[ix][iy][iz][k] = q1;
                        IsSolid[nix][niy][niz] = true;
                        if(flag2)
                        {
                            if(std::abs(q1-q2)>1e-6)
                            std::cout<<"!!! "<<q1<<" "<<q2<<std::endl;
                        }    
                    }else{
                        if(flag2)
                        {
                            q[ix][iy][iz][k] = q2;
                            IsSolid[nix][niy][niz] = true;
                                
                        }else{
                            q[ix][iy][iz][k] = -1;
                        }
                    }
                }
            }else{
            //    IsSolid[ix][iy][iz] = true; 
            }
        }
        // if(L<1.1*R) IsSolid[ix][iy][iz] = true;
    }
}

inline void Domain::AddSphereQ(Vec3_t &pos, double R)
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    double px = pos(0);
    double py = pos(1);
    double pz = pos(2); 
    int xstart = std::max(px-R-2,0.0);
    int ystart = std::max(py-R-2,0.0);
    int zstart = std::max(pz-R-2,0.0);
    int xend = std::min(px+R+2,(double) nx);
    int yend = std::min(py+R+2,(double) ny);
    int zend = std::min(pz+R+2,(double) nz);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for(int ix= xstart; ix<xend; ++ix)
    for(int iy= ystart; iy<yend; ++iy)
    for(int iz= zstart; iz<zend; ++iz)  
    {
        double x = (double) ix;
        double y = (double) iy;
        double z = (double) iz;
        Vec3_t Cell(x,y,z);
        double L = dot(Cell-pos,Cell-pos);
        if(L<R*R) 
        {
            IsSolid[ix][iy][iz] = true;
            continue;
        }
        for(size_t k=0; k<Nneigh; ++k)
        {
            size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)nx)%nx;
            size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)ny)%ny;
            size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)nz)%nz;
            x = (double) nix;
            y = (double) niy;
            z = (double) niz;
            Vec3_t Cell1(x,y,z);
            double L1 = dot(Cell1-pos,Cell1-pos);
            if(L1>R*R) continue;
            // Gamma[ix][iy][iz] = 1.0;                
            
            double K = dot(Cell-pos,Cell1-pos);
            double A = (L + L1- 2.0*K);
            double B = (2.0*K - 2.0*L);
            double C = L - R*R;
            double delta = B*B - 4.0*A*C;
            if(delta<0)
            {
                q[ix][iy][iz][k] = -1;
                continue;
            }
            double q1 = (-B + std::sqrt(delta))/(2.0*A);
            double q2 = (-B - std::sqrt(delta))/(2.0*A);
            bool flag1 = q1>=0 && q1-1 <1e-6;
            bool flag2 = q2>=0 && q2-1 <1e-6;
            if(flag1)
            {
                q[ix][iy][iz][k] = q1;
                IsSolid[nix][niy][niz] = true;
                if(std::abs(q1-q2)>1e-6)
                {
                     std::cout<<"!!!!!"<<ix<<" "<<iy<<" "<<iz<<" "<<k<<" "<<q1<<" "<<q2<<std::endl;
                }
            }else{
                if(flag2)
                {
                    q[ix][iy][iz][k] = q2;
                    IsSolid[nix][niy][niz] = true;
                
                    
                }else{
                    q[ix][iy][iz][k] = -1.0;
                    
                }

            }
            
                
        }
        
    }
}

inline void Domain::BounceBack(bool calcF)
{
    if(Time<0.5) std::cout<<"--- "<<"SBB"<<" ---"<<std::endl;
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx; ix++)
	for(size_t iy=0; iy<ny; iy++)
	for(size_t iz=0; iz<nz; iz++)
    {
    // for(size_t i=0; i<Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        

        for(size_t k=0; k<Nneigh; k++)
        {
            size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)nx)%nx;
            size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)ny)%ny;
            size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)nz)%nz;
            
            if(IsSolid[nix][niy][niz])
            {   
                F[ix][iy][iz][Op[k]] = Ftemp[ix][iy][iz][k];
                
            }
            
        }
        if(!calcF) continue;
        
        Vec3_t Flbmt(0,0,0);
        for(size_t k=0; k<Nneigh; ++k)
        {
            Flbmt += (Ftemp[ix][iy][iz][k]+F[ix][iy][iz][Op[k]])*C[k] - 2.0*W[k]*1.0*3.0*dot(C[Op[k]],VelP[ix][iy][iz])/Cs/Cs*VelP[ix][iy][iz];
        }
        Flbm[ix][iy][iz] = Flbmt;

	    
    }
}

inline void Domain::BounceBackLIBB(bool calcF)
{
                    
    if(Time<0.5) std::cout<<"--- "<<"LIBB"<<" ---"<<std::endl;    
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx; ix++)
	for(size_t iy=0; iy<ny; iy++)
	for(size_t iz=0; iz<nz; iz++)
    {
    // for(size_t i=0; i<Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        // double Temp[Nneigh];
        // double Temp1[Nneigh];
        for(size_t k=0; k<Nneigh; k++)
        {
            double qt = q[ix][iy][iz][k];
            if(qt<0) continue;
            double fw[Nneigh];
            size_t oix = (size_t)((int)ix + (int)C[Op[k]](0) + (int)nx)%nx;
            size_t oiy = (size_t)((int)iy + (int)C[Op[k]](1) + (int)ny)%ny;
            size_t oiz = (size_t)((int)iz + (int)C[Op[k]](2) + (int)nz)%nz;
            double *f = F[ix][iy][iz];
            double *fo = F[oix][oiy][oiz];
            double *ft = Ftemp[ix][iy][iz];
            double *fto = Ftemp[oix][oiy][oiz];
            
            fw[k] = qt*ft[k] + (1-qt) * fto[k];
            if(IsSolid[oix][oiy][oiz])
            {
                f[Op[k]] = ft[k];
            }else{
                //Peng et al 2016 Yu's method
                double fu = 2.0*W[k]*1.0*3.0*dot(C[Op[k]],VelP[ix][iy][iz]);
                f[Op[k]] = 1.0/(1+qt)*(fw[k]+fu) + qt/(1+qt)*fo[Op[k]];
            }
        
            
            
            
        }
        if(!calcF) continue;
        Vec3_t Flbmt(0,0,0);
        for(size_t k=0; k<Nneigh; ++k)
        {
            double qt = q[ix][iy][iz][k];
            if(qt<0) continue;
            Flbmt += (Ftemp[ix][iy][iz][k]+F[ix][iy][iz][Op[k]])*C[k] - 2.0*W[k]*1.0*3.0*dot(C[Op[k]],VelP[ix][iy][iz])/Cs/Cs*VelP[ix][iy][iz];
        }
        Flbm[ix][iy][iz] = Flbmt;

        

    }   
}

inline void Domain::BounceBackQIBB(bool calcF)
{
    if(Time<0.5) std::cout<<"--- "<<"QIBB"<<" ---"<<std::endl;    
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx; ix++)
	for(size_t iy=0; iy<ny; iy++)
	for(size_t iz=0; iz<nz; iz++)
    {
    // for(size_t i=0; i<Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        

        for(size_t k=0; k<Nneigh; k++)
        {
            double qt = q[ix][iy][iz][k];
            if(qt<0) continue;
            double fw[Nneigh];
            size_t oix = (size_t)((int)ix + (int)C[Op[k]](0) + (int)nx)%nx;
            size_t oiy = (size_t)((int)iy + (int)C[Op[k]](1) + (int)ny)%ny;
            size_t oiz = (size_t)((int)iz + (int)C[Op[k]](2) + (int)nz)%nz;
        
            size_t ooix = (size_t)((int)oix + (int)C[Op[k]](0) + (int)nx)%nx;
            size_t ooiy = (size_t)((int)oiy + (int)C[Op[k]](1) + (int)ny)%ny;
            size_t ooiz = (size_t)((int)oiz + (int)C[Op[k]](2) + (int)nz)%nz;
            double *f = F[ix][iy][iz];
            double *foo = F[ooix][ooiy][ooiz];
            double *fo = F[oix][oiy][oiz];
            double *ft = Ftemp[ix][iy][iz];
            double *fto = Ftemp[oix][oiy][oiz];
            double *ftoo = Ftemp[ooix][ooiy][ooiz];
            if(IsSolid[ooix][ooiy][ooiz])
            {
                fw[k] = qt*ft[k] + (1-qt) * fto[k];
                if(IsSolid[oix][oiy][oiz])
                {
                    f[Op[k]] = ft[k];
                }else{
                    //Peng et al 2016 Yu's method
                    double fu = 2.0*W[k]*1.0*3.0*dot(C[Op[k]],VelP[ix][iy][iz]);
                    f[Op[k]] = 1.0/(1+qt)*(fw[k]+fu) + qt/(1+qt)*fo[Op[k]];
                }
            }else{
                fw[k] = qt*(1.0+qt)/2.0*ft[k] + (1.0-qt)*(1.0+qt)*fto[k] - qt*(1-qt)/2.0*ftoo[k];
                //Peng et al 2016 Yu's method
                double fu = 2.0*W[k]*1.0*3.0*dot(C[Op[k]],VelP[ix][iy][iz]);
                f[Op[k]] = 2.0/((1+qt)*(2+qt))*(fw[k]+fu) + 2.0*qt/(1+qt)*fo[Op[k]] - qt/(2.0+qt)*foo[Op[k]];
            }
            
                
            
        }
        if(!calcF) continue;
        Vec3_t Flbmt(0,0,0);
        for(size_t k=0; k<Nneigh; ++k)
        {
            double qt = q[ix][iy][iz][k];
            if(qt<0) continue;
            Flbmt += (Ftemp[ix][iy][iz][k]+F[ix][iy][iz][Op[k]])*C[k] - 2.0*W[k]*1.0*3.0*dot(C[Op[k]],VelP[ix][iy][iz])/Cs/Cs*VelP[ix][iy][iz];
        }
        Flbm[ix][iy][iz] = Flbmt;

        

    }   
}

inline void Domain::BounceBackMR(bool calcF)
{
    if(Time<0.5) std::cout<<"--- "<<"MR"<<" ---"<<std::endl;    
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx; ix++)
	for(size_t iy=0; iy<ny; iy++)
	for(size_t iz=0; iz<nz; iz++)
    {
    // for(size_t i=0; i<Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        for(size_t k=0; k<Nneigh; k++)
        {
            double qt = q[ix][iy][iz][k];
            if(qt<0) continue;
            double fw[Nneigh];
            size_t oix = (size_t)((int)ix + (int)C[Op[k]](0) + (int)nx)%nx;
            size_t oiy = (size_t)((int)iy + (int)C[Op[k]](1) + (int)ny)%ny;
            size_t oiz = (size_t)((int)iz + (int)C[Op[k]](2) + (int)nz)%nz;
        
            size_t ooix = (size_t)((int)oix + (int)C[Op[k]](0) + (int)nx)%nx;
            size_t ooiy = (size_t)((int)oiy + (int)C[Op[k]](1) + (int)ny)%ny;
            size_t ooiz = (size_t)((int)oiz + (int)C[Op[k]](2) + (int)nz)%nz;
            double *f = F[ix][iy][iz];
            // double *foo = F[ooix][ooiy][ooiz];
            double *fo = F[oix][oiy][oiz];
            double *ft = Ftemp[ix][iy][iz];
            double *fto = Ftemp[oix][oiy][oiz];
            double *ftoo = Ftemp[ooix][ooiy][ooiz];
            if(IsSolid[ooix][ooiy][ooiz])
            {
                fw[k] = qt*ft[k] + (1-qt) * fto[k];
                if(IsSolid[oix][oiy][oiz])
                {
                    // f[Op[k]] = ft[k];
                    // Vec3_t v0(0.0,0.0,0.0); 
                    // f[k] = Feq(k,1.0,v0);
                }else{
                    //Peng et al 2016 Yu's method
                    double fu = 2.0*W[k]*1.0*dot(C[Op[k]],VelP[ix][iy][iz]);
                    f[Op[k]] = 1.0/(1+qt)*(fw[k]+fu) + qt/(1+qt)*fo[Op[k]];
                }
            }else{
                double k1 = (1-2.0*qt-2.0*qt*qt)/((1.0+qt)*(1.0+qt));
                double k2 = qt*qt/((1.0+qt)*(1.0+qt));
                // double fu = 4.0/((1.0+qt)*(1.0+qt))*W[k]*1.0*3.0*dot(C[Op[k]],VelP[ix][iy][iz]);
                // double FF =  -(t[ix][iy][iz][k]/Nu)*1.0/(4.0*(1+qt)*(1+qt));
                // std::cout<<k<<" "<<t[ix][iy][iz][k]<<" "<<FF<<" "<<f[Op[k]]<<std::endl;
                // double FF =  -t[ix][iy][iz][k]*(4.0/(3.0*Nu*(1+qt)*(1+qt)))*((Tau-0.5)/(1.0/S(k)-0.5));
                f[Op[k]] = ft[k] + k1*fto[k] + k2*ftoo[k] - k1*ft[Op[k]] - k2*fto[Op[k]] ;
                if(f[Op[k]]<0) throw new Fatal("Negative distribution function value!!!!!");   
            }
            
                
            
        }
       
        if(!calcF) continue;
        Vec3_t Flbmt(0,0,0);
        for(size_t k=0; k<Nneigh; ++k)
        {
            double qt = q[ix][iy][iz][k];
            if(qt<0) continue;
            Flbmt += (Ftemp[ix][iy][iz][k]+F[ix][iy][iz][Op[k]])*C[k] - 2.0*W[k]*1.0*3.0*dot(C[Op[k]],VelP[ix][iy][iz])/Cs/Cs*VelP[ix][iy][iz];
        }
        Flbm[ix][iy][iz] = Flbmt;

        

    }   
}

inline void Domain::BounceBackCLI(bool calcF)
{
    if(Time<0.5) std::cout<<"--- "<<"CLI"<<" ---"<<std::endl;    
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx; ix++)
	for(size_t iy=0; iy<ny; iy++)
	for(size_t iz=0; iz<nz; iz++)
    {
    // for(size_t i=0; i<Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        
        
        for(size_t k=0; k<Nneigh; k++)
        {
            double qt = q[ix][iy][iz][k];
            if(qt<0) continue;
            // double fw[Nneigh];
            size_t oix = (size_t)((int)ix + (int)C[Op[k]](0) + (int)nx)%nx;
            size_t oiy = (size_t)((int)iy + (int)C[Op[k]](1) + (int)ny)%ny;
            size_t oiz = (size_t)((int)iz + (int)C[Op[k]](2) + (int)nz)%nz;
            double *f = F[ix][iy][iz];
            // double *fo = F[oix][oiy][oiz];
            double *ft = Ftemp[ix][iy][iz];
            double *fto = Ftemp[oix][oiy][oiz];
            
            
            if(IsSolid[oix][oiy][oiz])
            {
                f[Op[k]] = ft[k];
            }else{
                
                // double fu = 2.0*W[k]*1.0*dot(C[Op[k]],VelP[ix][iy][iz]);
                f[Op[k]] = (1.0-2.0*qt)/(1+2.0*qt)*fto[k] + ft[k] - (1.0-2.0*qt)/(1.0+2.0*qt)*ft[Op[k]];
                
            }
            
        }
        if(!calcF) continue;
        Vec3_t Flbmt(0,0,0);
        for(size_t k=0; k<Nneigh; ++k)
        {
            double qt = q[ix][iy][iz][k];
            if(qt<0) continue;
            Flbmt += (Ftemp[ix][iy][iz][k]+F[ix][iy][iz][Op[k]])*C[k] - 2.0*W[k]*1.0*3.0*dot(C[Op[k]],VelP[ix][iy][iz])/Cs/Cs*VelP[ix][iy][iz];
        }
        Flbm[ix][iy][iz] = Flbmt;

        

    }   
}


inline void Domain::CalcForceQ()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
	for(size_t ix=0; ix<nx;ix++)
	for(size_t iy=0; iy<ny;iy++)
	for(size_t iz=0; iz<nz;iz++)
    {
    // for(size_t i=0; i<Ncells; ++i)
	// {
    //     iVec3_t iv(0,0,0);
    //     idx2Pt(i,iv);
    //     size_t ix = iv(0);
    //     size_t iy = iv(1);
    //     size_t iz = iv(2);
        Vec3_t Flbmt(0,0,0);
        for(size_t k=0; k<Nneigh; ++k)
        {
            double qt = q[ix][iy][iz][k];
            if(qt<0) continue;
            // size_t oix = (size_t)((int)ix + (int)C[Op[k]](0) + (int)nx)%nx;
            // size_t oiy = (size_t)((int)iy + (int)C[Op[k]](1) + (int)ny)%ny;
            // size_t oiz = (size_t)((int)iz + (int)C[Op[k]](2) + (int)nz)%nz;
            Flbmt += (Ftemp[ix][iy][iz][k]+F[ix][iy][iz][Op[k]])*C[k];
        }
        Flbm[ix][iy][iz] = Flbmt;
    }

}



    

























#endif