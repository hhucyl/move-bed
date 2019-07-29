#ifndef LBM_INITIALIZE_H
#define LBM_INITIALIZE_H

inline void Domain::Initialize(iVec3_t idx, double TheRho, Vec3_t & TheVel)
{
    size_t ix = idx(0);
    size_t iy = idx(1);
    size_t iz = idx(2);

    BForce[ix][iy][iz] = OrthoSys::O;

    for (size_t k=0;k<Nneigh;k++)
    {
        F[ix][iy][iz][k] = Feq(k,TheRho,TheVel);
    }

    if (!IsSolid[ix][iy][iz])
    {
        Vel[ix][iy][iz] = TheVel;
        Rho[ix][iy][iz] = TheRho;
    }
    else
    {
        Vel[ix][iy][iz] = OrthoSys::O;
        Rho[ix][iy][iz] = 0.0;
    }
}

inline void Domain::Initial(double rho, Vec3_t &v0,  Vec3_t &g0)
{
    
    for(size_t ix=0; ix<Ndim(0); ix++)
    for(size_t iy=0; iy<Ndim(1); iy++)
    for(size_t iz=0; iz<Ndim(2); iz++)
    {
        Rho[ix][iy][iz] = rho;
        Vel[ix][iy][iz] = 0.0, 0.0, 0.0;
        BForce[ix][iy][iz] = g0;
        for(size_t k=0; k<Nneigh; ++k)
        {
            F[ix][iy][iz][k] = Feq(k,rho,v0);            
            Ftemp[ix][iy][iz][k] = Feq(k,rho,v0);            
        }
    // std::cout<<F[ix][iy][iz][18]<<std::endl;
        
    }
    Rho0 = rho;//very important
}

inline void Domain::InitialFromH5(char const * TheFileKey, Vec3_t &g0)
{
    String fn(TheFileKey);
    std::cout<<"Initializing From "<<fn.CStr()<<std::endl;
    hid_t file_id = H5Fopen(fn.CStr(),H5F_ACC_RDONLY,H5P_DEFAULT);
    double *Ff = new double[Nneigh*Ndim(0)*Ndim(1)*Ndim(2)];
    double *Ga = new double[Ndim(0)*Ndim(1)*Ndim(2)];
    double *BF = new double[3*Ndim(0)*Ndim(1)*Ndim(2)];
    H5LTread_dataset_double(file_id,"/F_0",Ff);
    H5LTread_dataset_double(file_id,"/Gamma",Ga);
    H5LTread_dataset_double(file_id,"/BForce_0",BF);
    size_t nn=0;
	//fluid
    for (size_t iz=0; iz<Ndim(2); ++iz)
    for (size_t iy=0; iy<Ndim(1); ++iy)
    for (size_t ix=0; ix<Ndim(0); ++ix)
    {
		double * f  = F[ix][iy][iz];
        Vec3_t bf(BF[3*nn],BF[3*nn+1],BF[3*nn+2]);
        Vel   [ix][iy][iz] = OrthoSys::O;
        Rho   [ix][iy][iz] = 0.0;
        BForce[ix][iy][iz] = bf;
        Gamma [ix][iy][iz] = Ga[nn];
		for(size_t k=0; k<Nneigh; k++)
		{
			
			f[k] = Ff[Nneigh*nn + k];
            Rho[ix][iy][iz] += f[k];
            Vel[ix][iy][iz] += f[k]*C[k];
		}
        Vel[ix][iy][iz] *= Cs/Rho[ix][iy][iz];
        if(IsSolid[ix][iy][iz])
        {
            Vel   [ix][iy][iz] = OrthoSys::O;
            Rho   [ix][iy][iz] = 0.0;
        }
		nn++;
	}
    std::cout<<"Fluid Complete"<<std::endl;
    //particle
    int NP = Particles.size(); 
    if(NP>0)
    {
        double *Ppos = new double[3*NP*2];
        H5LTread_dataset_double(file_id,"/Pposition",Ppos);
        double *Pposb = new double[3*NP*2];
        H5LTread_dataset_double(file_id,"/Ppositionb",Pposb);
        double *Pvec = new double[3*NP*2];
        H5LTread_dataset_double(file_id,"/PVeloc",Pvec);
        double *PW = new double[3*NP*2];
        H5LTread_dataset_double(file_id,"/PW",PW);
        double *PWb = new double[3*NP*2];
        H5LTread_dataset_double(file_id,"/PWb",PWb);
        int N[1];
        // std::cout<<1<<std::endl;
        H5LTread_dataset_int(file_id,"/PListRNum",N);
        int NLR = N[0];
        int *PlistR = new int[NLR*2];
        // std::cout<<2<<std::endl;

        H5LTread_dataset_int(file_id,"/PListR",PlistR);
        double *SFR = new double[NLR*3];
        // std::cout<<3<<std::endl;

        H5LTread_dataset_double(file_id,"/SFR",SFR);
        double *FDR = new double[NLR*3];
        // std::cout<<4<<std::endl;

        H5LTread_dataset_double(file_id,"/FDR",FDR);


        for(int i=0; i<NP; ++i)
        {
            DEM::Disk *Pa = &Particles[i];
            Pa->X = Ppos[3*i], Ppos[3*i+1], Ppos[3*i+2]; 
            Pa->Xb = Pposb[3*i], Pposb[3*i+1], Pposb[3*i+2]; 
            Pa->V = Pvec[3*i], Pvec[3*i+1], Pvec[3*i+2];
            Pa->W = PW[3*i], PW[3*i+1], PW[3*i+2]; 
            Pa->Wb = PWb[3*i], PWb[3*i+1], PWb[3*i+2];
        }

        std::map<std::pair<int,int>,Vec3_t> fmap;
        std::map<std::pair<int,int>,Vec3_t> rmap;
        for(int i=0; i<NLR; ++i)
        {
            std::pair<int,int> p;
            p.first = PlistR[2*i];
            p.second = PlistR[2*i+1];
            Vec3_t sfr(SFR[3*i],SFR[3*i+1],SFR[3*i+2]);
            Vec3_t fdr(FDR[3*i],FDR[3*i+1],FDR[3*i+2]);
            fmap[p] = sfr;
            rmap[p] = fdr;
        }
        Friction = fmap;
        Rolling = rmap;

        delete[] Ppos;
        delete[] Pposb;
        delete[] Pvec;
        delete[] PW;
        delete[] PWb;
		delete[] PlistR;
		delete[] SFR;
		delete[] FDR;
    }


	delete[] Ff;
	H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
}

#endif
