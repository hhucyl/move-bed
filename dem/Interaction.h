/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/


#ifndef DEM_INTERACTION_H
#define DEM_INTERACTION_H

// Mechsys
#include "Disk.h"
#include <mechsys/dem/basic_functions.h>
namespace DEM
{
class DiskPair
{
public:
    //Constructor
    DiskPair () {};
    DiskPair(Disk * D1, Disk * D2);

    //Methods
    void CalcForce      (double dt);
    void CalcD      ();
    bool UpdateContacts (double Alpha);

#ifdef USE_THREAD
    pthread_mutex_t lck;   ///< Lock to protect variables from race conditions.
#endif

    //Data
    Disk * P1;       ///< Pointer to first particle
    Disk * P2;       ///< Pointer to second particle
    double     Kn;       ///< Normal Spring constant 
    double     Kt;       ///< Tangential Spring constant 
    double     Gn;       ///< dissipation constant
    double     Gt;       ///< dissipation constant
    double     Mu;       ///< Friction coefficient
    double     Beta;     ///< Rolling stiffness coeffcient
    double     Eta;      ///< Plastic moment coefficient
    double     delta;
    Vec3_t    SFr;       ///< Vector of static friction
    Vec3_t    Fdr;       ///< Vector of rolling resistance
    Vec3_t     F1;       ///< net force over particle 1
    Vec3_t     F2;       ///< net force over particle 2
    Vec3_t     T1;       ///< net torque over particle 1
    Vec3_t     T2;       ///< net torque over particle 2
};

DiskPair::DiskPair(Disk * Dp1, Disk * Dp2)
{
    P1 = Dp1;
    P2 = Dp2;
    Kn = 2.0*ReducedValue(P1->Kn,P2->Kn);
    Kt = 2.0*ReducedValue(P1->Kt,P2->Kt);
    Gn = 2.0*ReducedValue(P1->Gn,P2->Gn);
    Gt = 2.0*ReducedValue(P1->Gt,P2->Gt);
    Mu = 2.0*ReducedValue(P1->Mu,P2->Mu);
    Eta= 2.0*ReducedValue(P1->Eta,P2->Eta);
    Beta= 2.0*ReducedValue(P1->Beta,P2->Beta);
    SFr= OrthoSys::O;
    Fdr= OrthoSys::O;
    double me = ReducedValue(P1->M,P2->M);
    // std::cout<<Gn<<std::endl;
    if (Gn < 0.0)
    {
        if (fabs(Gn)>1.0) throw new Fatal("CInteractonSphere the restitution coefficient is greater than 1");
        Gn = 2.0*sqrt((pow(log(-Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-Gn),2.0)));
        Gt = 0.0;
    }
    Gn *= me;
    Gt *= me;
#ifdef USE_THREAD
    pthread_mutex_init(&lck,NULL);
#endif
}

void DiskPair::CalcD()
{   
    double dist  = norm(P2->X - P1->X);
    delta = P1->R + P2->R - dist;
    if(P1->X(0)<(-10*P1->R)||P2->X(0)<(-10*P2->R)) delta = -1000;//couple with periodic
}

void DiskPair::CalcForce(double dt)
{
    F1      = OrthoSys::O;
    F2      = OrthoSys::O;
    T1      = OrthoSys::O;
    T2      = OrthoSys::O;
    
    double dist  = norm(P2->X - P1->X);
    delta = P1->R + P2->R - dist;
    if(P1->X(0)<(-10*P1->R)||P2->X(0)<(-10*P2->R)) delta = -1000;//couple with periodic
    if (delta>0)
    {

        //Force
        Vec3_t n    = (P2->X - P1->X)/dist;
        Vec3_t x    = P1->X+n*((P1->R*P1->R-P2->R*P2->R+dist*dist)/(2*dist));
        // Vec3_t Fn   = Kn*std::pow(delta,1.5)*n;
        Vec3_t Fn   = Kn*delta*n;
        Vec3_t t1,t2,x1,x2;
        Rotation(P1->W,P1->Q,t1);
        Rotation(P2->W,P2->Q,t2);
        x1 = x - P1->X;
        x2 = x - P2->X;
        Vec3_t Vrel = -((P2->V-P1->V)+cross(t2,x2)-cross(t1,x1));
        Vec3_t vt = Vrel - dot(n,Vrel)*n;
        SFr      += dt*vt;
        SFr      -= dot(SFr,n)*n;
        Vec3_t tan= SFr;
        if(norm(tan)>0.0) tan/=norm(tan);
        if(norm(SFr)>Mu*norm(Fn)/Kt)
        {
            SFr = Mu*norm(Fn)/Kt*tan;
        }
        Vec3_t F    = Fn + Gn*dot(n,Vrel)*n + Kt*SFr + Gt*vt;
        // Vec3_t F    = Fn + Gn*dot(n,Vrel)*n;
        // std::cout<<F(1)<<std::endl;
        if (dot(F,n)<0) F-=dot(F,n)*n;

        
        //Rolling resistance
        Vec3_t Vr = P1->R*P2->R*cross(Vec3_t(t1 - t2),n)/(P1->R+P2->R);
        Fdr += Vr*dt;
        Fdr -= dot(Fdr,n)*n;
        tan = Fdr;
        if (norm(tan)>0.0) tan/=norm(tan);
        double Kr = Beta*Kt;
        if (norm(Fdr)>Eta*Mu*norm(Fn)/Kr)
        {
            Fdr = Eta*Mu*norm(Fn)/Kr*tan;
        }
        Vec3_t Ft = -Kr*Fdr;

        //Assigning torque values
        Vec3_t T1 = -cross(x1,F) + P1->R*cross(n,Ft);
        Vec3_t T2 =  cross(x2,F) - P2->R*cross(n,Ft);
        Quaternion_t q;
        Conjugate (P1->Q,q);
        Rotation  (T1,q,t1);
        Conjugate (P2->Q,q);
        Rotation  (T2,q,t2);

        F1      = -F;
        F2      =  F;
        T1      = t1;
        T2      = t2;

        //std::cout << F << " " << t1 << P1->IsFree() << P2->IsFree() << std::endl;

    }else{
        // std::cout<<1<<std::endl;
        // double e1 = P1->e1;
        // double eal = P1->eal;//pay attention to 0.125 this may vary from case to case see Brandle de Motta et al 2014
        // auto lambda = [](double d) {return (0.5/d - 9.0/20.0*std::log(d) - 3.0/56.0*d*std::log(d) + 1.346);};
        // double RR = 2*P1->R*P2->R/(P1->R+P2->R);
        // double ee = -delta/RR;
        // double mu = 1.0*P1->nu;
        // if(ee<e1 && ee>0)
        // {
        //     Vec3_t n    = (P1->X - P2->X)/dist;
        //     double Ft = std::fabs(-6*M_PI*mu*RR*dot((P1->V-P2->V),n)*(lambda(e1)-lambda(eal)));
        //     F1      =  Ft*n;
        //     F2      =  -Ft*n; 
        // }else{
        //     if(ee<=eal && ee>0)
        //     {
        //         Vec3_t n    = (P1->X - P2->X)/dist;
                
        //         double Ft = std::fabs(6*M_PI*mu*RR*dot((P1->V-P2->V),n)*(lambda(ee)-lambda(eal)));
        //         F1      =  Ft*n;
        //         F2      =  -Ft*n; 
        //     }

        // }
    }
}

bool DiskPair::UpdateContacts(double Alpha)
{
    if (norm(P1->X-P2->X) <= P1->R + P2->R + 2*Alpha) return true;
    else                                              return false;
}
}
#endif
