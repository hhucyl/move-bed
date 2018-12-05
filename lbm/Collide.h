#ifndef LBM_Collide_H
#define LBM_Collide_H

inline void Domain::CollideSRT()
{
    if(Time<0.5) std::cout<<"--- "<<"SRT"<<" ---"<<std::endl;    
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[ix][iy][iz])
        {
            double NonEq[Nneigh];
            double Q = 0.0;
            double tau = Tau;
            double rho = Rho[ix][iy][iz];
            Vec3_t vel = Vel[ix][iy][iz];
            double VdotV = dot(vel,vel);
            for (size_t k=0;k<Nneigh;k++)
            {
                double VdotC = dot(vel,C[k]);
                double Feq   = W[k]*rho*(1.0 + 3.0*VdotC + 4.5*VdotC*VdotC - 1.5*VdotV);
                NonEq[k] = F[ix][iy][iz][k] - Feq;
                Q +=  NonEq[k]*NonEq[k]*EEk[k];
            }
            Q = sqrt(2.0*Q);
            tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*Sc/rho));

            // bool valid = true;
            // double alpha = 1.0;
            // while (valid)
            // {
                // valid = false;
            double Bn = (Gamma[ix][iy][iz]*(Tau-0.5))/((1.0-Gamma[ix][iy][iz])+(Tau-0.5));
            // Bn = Gamma[ix][iy][iz];    
                for (size_t k=0;k<Nneigh;k++)
                {
                    double ForceTerm = dt*3.0*W[k]*dot(BForce[ix][iy][iz],C[k]);
                    // Vec3_t BFt(0.0, 0.0, 0.0);
                    // BFt = 3.0*(C[k] - vel)/(Cs*Cs) + 9.0*dot(C[k],vel)/(Cs*Cs*Cs*Cs)*C[k]; 
                    // double ForceTerm = dt*(1 - 1.0/(2.0*Tau))*W[k]*dot(BFt,BForce[ix][iy][iz]); 
                    Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k] - (1-Bn)*NonEq[k]/tau + Bn*Omeis[ix][iy][iz][k]+ForceTerm;
                    // if (Ftemp[ix][iy][iz][k]<-1.0e-12)
                    // {
                    //     //std::cout << Ftemp[ix][iy][iz][k] << std::endl;
                    //     double temp =  tau*F[ix][iy][iz][k]/NonEq[k];
                    //     if (temp<alpha) alpha = temp;
                    //     valid = true;
                    // }
                    // if (std::isnan(Ftemp[ix][iy][iz][k]))
                    // {
                    //     std::cout << "CollideSC: Nan found, resetting" << std::endl;
                    //     std::cout << " " << alpha << " " << iVec3_t(ix,iy,iz) << " " << k << " " << std::endl;
                    //     throw new Fatal("Domain::CollideSC: Distribution funcitons gave nan value, check parameters");
                    // }
                }
            // }
        }
        else
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][Op[k]];
            }
        }
    }
}

inline void Domain::CollideMRT()
{
    if(Time<0.5) std::cout<<"--- "<<"MRT"<<" ---"<<std::endl;    
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[ix][iy][iz])
        {
            //double NonEq[Nneigh];
            //double Q = 0.0;
            //double tau = Tau;
            double rho = Rho[ix][iy][iz];
            Vec3_t vel = Vel[ix][iy][iz];
            // vel += 0.5*dt*BForce[ix][iy][iz];
            //double VdotV = dot(vel,vel);
            //for (size_t k=0;k<Nneigh;k++)
            //{
                //double VdotC = dot(vel,C[k]);
                //double Feq   = W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                //NonEq[k] = F[ix][iy][iz][k] - Feq;
                //Q +=  NonEq[k]*NonEq[k]*EEk[k];
            //}
            //Q = sqrt(2.0*Q);
            //tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*Sc/rho));
            double *f = F[ix][iy][iz];
            double *m = Ftemp[ix][iy][iz];
            double fneq[Nneigh];
            int n=Nneigh,mm=1;
            double a = 1.0,b = 0.0;
            dgemv_("N",&n,&n,&a,M.data,&n,f,&mm,&b,m,&mm);
            // for(size_t i=0; i<Nneigh; i++)
            // {
            //     m[i] = 0.0;
            //     for(size_t j=0; j<Nneigh; j++)
            //     {
            //         m[i] += M(i,j)*f[j];                   
            //     }
            // }
            
            
            (this->*ptr2meq)(m,rho,vel);
            // MeqD3Q19(m,rho,vel);
            
            dgemv_("N",&n,&n,&a,Minv.data,&n,m,&mm,&b,fneq,&mm);
            // for(size_t i=0; i<Nneigh; i++)
            // {
            //     fneq[i] = 0.0;
            //     for(size_t j=0; j<Nneigh; j++)
            //     {
            //         fneq[i] += Minv(i,j)*m[j];                   
            //     }
            // }
            double Bn = (Gamma[ix][iy][iz]*(Tau-0.5))/((1.0-Gamma[ix][iy][iz])+(Tau-0.5));
            //Bn = Gamma[ix][iy][iz];
            for (size_t k=0; k<Nneigh; k++)
            {
                double ForceTerm = dt*3.0*W[k]*dot(BForce[ix][iy][iz],C[k]);
                // Vec3_t BFt(0.0, 0.0, 0.0);
                // BFt = 3.0*(C[k] - vel)/(Cs*Cs) + 9.0*dot(C[k],vel)/(Cs*Cs*Cs*Cs)*C[k]; 
                // double ForceTerm = dt*(1 - 1.0/(2.0*Tau))*W[k]*dot(BFt,BForce[ix][iy][iz]); 
                // Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k]-fneq[k]+ForceTerm;
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k] - (1-Bn)*fneq[k] + Bn*Omeis[ix][iy][iz][k] + ForceTerm;
                if(Ftemp[ix][iy][iz][k]<0) Ftemp[ix][iy][iz][k] = 0.0;

                // if(Ftemp[ix][iy][iz][k]<0)
                // {
                    // std::cout<<ix<<" "<<iy<<" "<<iz<<std::endl;
                    // throw new Fatal("Negative distribution function value!!!!!");   
                // }
                // if(ix == 1 &&iy==10&&iz ==0)
                // {
                //     std::cout<<"F "<<k<<"="<<F[ix][iy][iz][k]<<std::endl;
                //     std::cout<<"Ft "<<k<<"="<<Ftemp[ix][iy][iz][k]<<std::endl;
                // }
            }

            //bool valid = true;
            //double alpha = 1.0, alphat = 1.0;
            //while (valid)
            //{
                //alpha = alphat;
                //valid = false;
              //  for (size_t k=0;k<Nneigh;k++)
            //    {
          //          Vec3_t BFt(0.0, 0.0, 0.0);
        //            BFt = 3.0*(C[k] - vel)/(Cs*Cs) + 9.0*dot(C[k],vel)/(Cs*Cs*Cs*Cs)*C[k]; 
      //              double ForceTerm = dt*(1 - 1.0/(2.0*Tau))*W[k]*dot(BFt,BForce[ix][iy][iz]); 
                    //Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k] - alpha*(fneq[k]-ForceTerm);
                    //if (Ftemp[ix][iy][iz][k]<-1.0e-12)
                    //{
                        //std::cout << Ftemp[ix][iy][iz][k] << std::endl;
                      //  valid = true;
                    //    double temp =  F[ix][iy][iz][k]/(fneq[k]-ForceTerm);
                  //      if (temp<alphat) alphat = temp;
                //    }
                    //if (std::isnan(Ftemp[ix][iy][iz][k]))
                    //{
                        //std::cout << "CollideSC: Nan found, resetting" << std::endl;
                        //std::cout << " " << alpha << " " << iVec3_t(ix,iy,iz) << " " << k << " " << std::endl;
                        //throw new Fatal("Domain::CollideSC: Distribution funcitons gave nan value, check parameters");
                    //}
              //  }
            //}
        }
        else
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][Op[k]];
            }
        }
    }

}

inline void Domain::CollideMRTLES()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[ix][iy][iz])
        {
            // double NonEq[Nneigh];
            
            
            double tau = Tau;
            double rho = Rho[ix][iy][iz];
            Vec3_t vel = Vel[ix][iy][iz];
            // vel += 0.5*dt*BForce[ix][iy][iz];
            // double VdotV = dot(vel,vel);
            // for (size_t k=0;k<Nneigh;k++)
            // {
            //     double VdotC = dot(vel,C[k]);
            //     double Feq   = W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
            //     NonEq[k] = F[ix][iy][iz][k] - Feq;
                //Q +=  NonEq[k]*NonEq[k]*EEk[k];
            //}
            //Q = sqrt(2.0*Q);
            //tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*Sc/rho));
            double *f = F[ix][iy][iz];
            double *m = Ftemp[ix][iy][iz];
            double fneq[Nneigh];
            int n=Nneigh,mm=1;
            double a = 1.0,b = 0.0;
            dgemv_("N",&n,&n,&a,M.data,&n,f,&mm,&b,m,&mm);
            
            double Q = 0.0;
            double drho = rho - Rho0;
            double P[3][3] = {0};
            P[0][0] = ((m[1]+2.0*drho) + m[9])/3.0;
            P[1][1] = P[0][0] + 0.5*(m[11]-m[9]);
            P[2][2] = P[1][1] - m[11];
            P[1][0] = P[0][1] = m[13];
            P[2][1] = P[1][2] = m[14];
            P[0][2] = P[2][0] = m[15];
            for(size_t im=0; im<3; im++)
            for(size_t in=0; in<3; in++)
            {
                double dd = im==in ? 1.0: 0.0;
                double QQ = dd*drho/3.0 + m[2*im+3]*m[2*in+3]-P[im][in];
                Q += QQ*QQ;
            }
            Q = std::sqrt(Q);
            tau = 0.5*(std::sqrt(tau*tau + 18.0*Sc*Sc*std::sqrt(Q))+tau);
            S(9) = S(11) = S(13) = S(14) = S(15) = 1/tau;

            (this->*ptr2meq)(m,rho,vel);
            // MeqD3Q19(m,rho,vel);
            
            dgemv_("N",&n,&n,&a,Minv.data,&n,m,&mm,&b,fneq,&mm);
            
            double Bn = (Gamma[ix][iy][iz]*(Tau-0.5))/((1.0-Gamma[ix][iy][iz])+(Tau-0.5));
            for (size_t k=0; k<Nneigh; k++)
            {
                double ForceTerm = dt*3.0*W[k]*dot(BForce[ix][iy][iz],C[k]);
                // Vec3_t BFt(0.0, 0.0, 0.0);
                // BFt = 3.0*(C[k] - vel)/(Cs*Cs) + 9.0*dot(C[k],vel)/(Cs*Cs*Cs*Cs)*C[k]; 
                // double ForceTerm = dt*(1 - 1.0/(2.0*Tau))*W[k]*dot(BFt,BForce[ix][iy][iz]); 
                
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k] - (1-Bn)*fneq[k] + Bn*Omeis[ix][iy][iz][k] + ForceTerm;
               
            }
        }else{
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][Op[k]];
            }
        }
    }
}



inline void Domain::CollideMRTMR()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[ix][iy][iz])
        {
            
            double rho = Rho[ix][iy][iz];
            Vec3_t vel = Vel[ix][iy][iz];
            double *f = F[ix][iy][iz];
            double *tt = t[ix][iy][iz];
            double *m = Ftemp[ix][iy][iz];
            double fneq[Nneigh];
            int n=Nneigh,mm=1;
            double a = 1.0,b = 0.0;
            dgemv_("N",&n,&n,&a,M.data,&n,f,&mm,&b,m,&mm);
            
            (this->*ptr2meq)(m,rho,vel);
            
            dgemv_("N",&n,&n,&a,Minv.data,&n,m,&mm,&b,fneq,&mm);

            tt[4] = fneq[4];
            tt[6] = fneq[6];
            tt[8] = fneq[8];
            tt[16] = fneq[16];
            tt[17] = fneq[17];
            tt[18] = fneq[18];
            
            
            double Bn = (Gamma[ix][iy][iz]*(Tau-0.5))/((1.0-Gamma[ix][iy][iz])+(Tau-0.5));
            for (size_t k=0; k<Nneigh; k++)
            {
                double ForceTerm = dt*3.0*W[k]*dot(BForce[ix][iy][iz],C[k]);
                
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k] - (1-Bn)*fneq[k] + Bn*Omeis[ix][iy][iz][k] + ForceTerm;
                if(Ftemp[ix][iy][iz][k]<0) Ftemp[ix][iy][iz][k] = 0; 
                
            }

        }
        else
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][Op[k]];
            }
        }
    }

}

inline void Domain::MeqD2Q9(double *m, double rho, Vec3_t &vel)
{
    m[0] = S(0)*(m[0] - rho);
    m[1] = S(1)*(m[1] - rho*(-2.0 + 3.0*dot(vel,vel)));
    m[2] = S(2)*(m[2] - rho*( 1.0 - 3.0*dot(vel,vel)));
    m[3] = S(3)*(m[3] - rho*vel(0));
    m[4] = S(4)*(m[4] + rho*vel(0));
    m[5] = S(5)*(m[5] - rho*vel(1));
    m[6] = S(6)*(m[6] + rho*vel(1));
    m[7] = S(7)*(m[7] - rho*(vel(0)*vel(0)-vel(1)*vel(1)));
    m[8] = S(8)*(m[8] - rho*(vel(0)*vel(1)));
}

inline void Domain::MeqD3Q15(double *m, double rho, Vec3_t &vel)
{
    m[ 0] = S(0)*(m[0] - rho); 
    m[ 1] = S( 1)*(m[ 1] + rho - rho*dot(vel,vel));
    m[ 2] = S( 2)*(m[ 2] + rho);
    m[ 3] = 0.0;
    m[ 4] = S( 4)*(m[ 4] + 7.0/3.0*rho*vel(0)); 
    m[ 5] = 0.0;
    m[ 6] = S( 6)*(m[ 6] + 7.0/3.0*rho*vel(1)); 
    m[ 7] = 0.0;
    m[ 8] = S( 8)*(m[ 8] + 7.0/3.0*rho*vel(2)); 
    m[ 9] = S( 9)*(m[ 9] - rho*(3.0*vel(0)*vel(0)-dot(vel,vel)));
    m[10] = S(10)*(m[10] - rho*(vel(1)*vel(1)-vel(2)*vel(2)));
    m[11] = S(11)*(m[11] - rho*(vel(0)*vel(1)));
    m[12] = S(12)*(m[12] - rho*(vel(1)*vel(2)));
    m[13] = S(13)*(m[13] - rho*(vel(0)*vel(2)));
    m[14] = S(14)* m[14];
}

inline void Domain::MeqD3Q19(double *m, double rho, Vec3_t &vel)
{
    m[ 0] = S( 0)*( m[ 0] - rho );
    m[ 1] = S( 1)*( m[ 1] - (-11.0*rho + 19.0 * dot(vel,vel)) );
    m[ 2] = S( 2)*( m[ 2] - (we*rho + wej*rho*dot(vel,vel)) );
    m[ 3] = 0.0;
    m[ 4] = S( 4)*( m[ 4] - (-2.0/3.0*rho*vel(0)) );
    m[ 5] = 0.0;
    m[ 6] = S( 6)*( m[ 6] - (-2.0/3.0*rho*vel(1)) );
    m[ 7] = 0.0;
    m[ 8] = S( 8)*( m[ 8] - (-2.0/3.0*rho*vel(2)) );
    m[ 9] = S( 9)*( m[ 9] - (2.0*rho*vel(0)*vel(0)-rho*vel(1)*vel(1)-rho*vel(2)*vel(2)) );
    m[10] = S(10)*( m[10] - (wxx*(2.0*rho*vel(0)*vel(0)-rho*vel(1)*vel(1)-rho*vel(2)*vel(2))) );
    m[11] = S(11)*( m[11] - rho*(vel(1)*vel(1)-vel(2)*vel(2)) );
    m[12] = S(12)*( m[12] - (-0.5*rho*(vel(1)*vel(1)-vel(2)*vel(2))) );
    m[13] = S(13)*( m[13] - rho*vel(0)*vel(1) );
    m[14] = S(14)*( m[14] - rho*vel(1)*vel(2) );
    m[15] = S(15)*( m[15] - rho*vel(0)*vel(2) );
    m[16] = S(16)*  m[16];
    m[17] = S(17)*  m[17];
    m[18] = S(18)*  m[18];

    // m[ 0] = S( 0)*( m[ 0] - rho );
    // m[ 1] = S( 1)*( m[ 1] - (-11.0*rho) );
    // m[ 2] = S( 2)*( m[ 2] - (3.0*rho) );
    // m[ 3] = 0.0;
    // m[ 4] = S( 4)*( m[ 4] - (-2.0/3.0*rho*vel(0)) );
    // m[ 5] = 0.0;
    // m[ 6] = S( 6)*( m[ 6] - (-2.0/3.0*rho*vel(1)) );
    // m[ 7] = 0.0;
    // m[ 8] = S( 8)*( m[ 8] - (-2.0/3.0*rho*vel(2)) );
    // m[ 9] = S( 9)*( m[ 9] );
    // m[10] = S(10)*( m[10] );
    // m[11] = S(11)*( m[11]  );
    // m[12] = S(12)*( m[12] );
    // m[13] = S(13)*( m[13]  );
    // m[14] = S(14)*( m[14]  );
    // m[15] = S(15)*( m[15]  );
    // m[16] = S(16)*  m[16];
    // m[17] = S(17)*  m[17];
    // m[18] = S(18)*  m[18];
}

inline void Domain::CollideSRTIBM()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[ix][iy][iz])
        {
            
            double tau = Tau;
            double rho = Rho[ix][iy][iz];
            Vec3_t bf = Flbm[ix][iy][iz];
            Vec3_t vel = Vel[ix][iy][iz] + 0.5*bf/rho*dt;
            double NonEq[Nneigh];
            for (size_t k=0;k<Nneigh;k++)
            {
                NonEq[k] = F[ix][iy][iz][k] - Feq(k,rho,vel);
            }
                
            for (size_t k=0;k<Nneigh;k++)
            {
                double ForceTerm1 = dt*3.0*W[k]*dot(BForce[ix][iy][iz],C[k]);
                Vec3_t BFt(0.0, 0.0, 0.0);
                BFt = 3.0*(C[k] - vel)/(Cs*Cs) + 9.0*dot(C[k],vel)/(Cs*Cs*Cs*Cs)*C[k]; 
                double ForceTerm = dt*(1 - 1.0/(2.0*Tau))*W[k]*dot(BFt,bf); 
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k] - NonEq[k]/tau + ForceTerm + ForceTerm1;
              
            }
        }
        else
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][Op[k]];
            }
        }
    }
}

inline void Domain::CollideMRTIBM()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[ix][iy][iz])
        {
            
            double tau = Tau;
            double rho = Rho[ix][iy][iz];
            Vec3_t bf = Flbm[ix][iy][iz];
            // Vec3_t vel = Vel[ix][iy][iz] + 0.5*bf/rho*dt;
            Vec3_t vel = Vel[ix][iy][iz];
            double *f = F[ix][iy][iz];
            double *m = Ftemp[ix][iy][iz];
            double fneq[Nneigh];
            int n=Nneigh,mm=1;
            double a = 1.0,b = 0.0;
            dgemv_("N",&n,&n,&a,M.data,&n,f,&mm,&b,m,&mm);
            
            (this->*ptr2meq)(m,rho,vel);
            
            dgemv_("N",&n,&n,&a,Minv.data,&n,m,&mm,&b,fneq,&mm);

                
            for (size_t k=0;k<Nneigh;k++)
            {
                double ForceTerm1 = dt*3.0*W[k]*dot(BForce[ix][iy][iz],C[k]);
                Vec3_t BFt(0.0, 0.0, 0.0);
                BFt = 3.0*(C[k] - vel)/(Cs*Cs) + 9.0*dot(C[k],vel)/(Cs*Cs*Cs*Cs)*C[k]; 
                double ForceTerm = dt*(1 - 1.0/(2.0*tau))*W[k]*dot(BFt,bf); 
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][k] - fneq[k]  + ForceTerm + ForceTerm1;
              
            }
        }
        else
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[ix][iy][iz][k] = F[ix][iy][iz][Op[k]];
            }
        }
    }
}

#endif