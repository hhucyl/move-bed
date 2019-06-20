#ifndef DEM_DISK_H
#define DEM_DISK_H

// STD
#include <algorithm>

// Mechsys
#include <mechsys/linalg/quaternion.h>

namespace DEM
{

class Disk
{
public:
    Disk(int Tag, Vec3_t const & X, Vec3_t const & V, Vec3_t const & W, double rho, double R, double dt);

    //Methods
    void Translate   (double dt);
    void Rotate      (double dt);
    void FixVeloc    () {vf=true,true,true; wf=true,true,true;};
    bool IsFree () {return !vf(0)&&!vf(1)&&!vf(2)&&!wf(0)&&!wf(1)&&!wf(2);}; ///< Ask if the particle has any constrain in its movement
    bool Ghost;
    void Periodic    (int modexy, Vec3_t &Box);
    void Leave    (int modexy, Vec3_t &Box);
#ifdef USE_THREAD
    pthread_mutex_t lck;   ///< Lock to protect variables from race conditions.
#elif USE_OMP
    omp_lock_t      lck;             ///< to protect variables in multithreading
#endif
    // Data
    int  Tag;              ///< Id of the particle
    Vec3_t X0;             ///< initial position of the particle
    Vec3_t X;              ///< position of the particle
    Vec3_t Xb;             ///< position of the particle before
    Vec3_t Xbt;
    Vec3_t V;              ///< velocity of the CM
    Vec3_t W;              ///< angular velocity
    Vec3_t Wb;             ///< angular velocity
    Vec3_t F;              ///< Force
    Vec3_t Fh;             ///< Force From hydraulics
    Vec3_t Flb;            ///< Force from lubrication
    Vec3_t Fc;             ///< Force From collision
    Vec3_t Ff;             ///< fixed Force
    Vec3_t T;              ///< Torque
    Vec3_t Th;             ///< Torque from hydraulics
    Vec3_t Tc;             ///< Torque from collision
    Vec3_t Tf;             ///< fixed Torque
    Quaternion_t Q;        ///< The quaternion representing the rotation
    double R;              ///< Disk radius
    double Rh;              ///< Disk radius
    double M;              ///< mass of the disk
    double I;              ///< inertia moment of the particle
    double Gn;             ///< dissipation constant for collision
    double Gt;             ///< dissipation constant for collision
    double Kn;             ///< Stiffness constant
    double Kt;             ///< Tangential stiffness constant
    double Mu;             ///< Friction coefficient
    double Beta;           ///< Rolling stiffness coeffcient
    double Eta;            ///< Plastic moment coefficient
    double nu;             ///< Water kinetic viscosity the dem in
    double eal;            ///< For lubrification 
    double e1;             ///< For lubrification
    bVec3_t vf;            ///< prescribed velocities
    bVec3_t wf;            ///< prescribed angular velocities

    // std::vector<bool> NodeType;
    // std::vector<Vec3_t> points;


};

inline Disk::Disk(int TheTag, Vec3_t const & TheX, Vec3_t const & TheV, Vec3_t const & TheW, double Therho, double TheR, double dt)
{
    Tag = TheTag;
    X   = TheX;
    V   = TheV;
    W   = TheW;
    R   = TheR;
    M   = M_PI*R*R*Therho;
    I   = 0.5*M*R*R;
    Xb  = X - dt*V;
    Xbt = Xb;
    Wb  = W;
    F   = 0.0,0.0,0.0;
    Fh   = 0.0,0.0,0.0;
    Fc   = 0.0,0.0,0.0;
    Flb   = 0.0,0.0,0.0;
    T   = 0.0,0.0,0.0;
    Ff  = 0.0,0.0,0.0;
    Tf  = 0.0,0.0,0.0;
    vf  = false,false,false;
    wf  = false,false,false;
    Gn  = 8.0;
    Gt  = 0.0;
    Kn  = 1.0e3;
    Kt  = 5.0e2;
    Mu  = 0.4;
    nu  = 0.0;
    eal = 0.0;
    e1 = 0.0;
    Eta = 1.0;  
    Beta = 0.12; 
    Q    = 1.0,0.0,0.0,0.0;
    Ghost = false;
#ifdef USE_THREAD
    pthread_mutex_init(&lck,NULL);
#elif USE_OMP
    omp_init_lock(&lck);
#endif
}

inline void Disk::Periodic(int modexy, Vec3_t &Box)
{
    //modexy is to distingusih perodic in xy direction, Box is where the periodic boundary was
    //remember that in periodic 0 is not equal to nx-1. It is actual at the location of nx since the definition of periodic boundary in LBM. This is accounting for + -1 in Periodic as well as Leave. This is perticularly important in IBM Periodic Boundary.
        
    // double dist1 = std::fabs(X(modexy) - Box(0));
    // double dist2 = std::fabs(X(modexy) - Box(1));
    double dist1 = X(modexy) - Box(0);
    double dist2 = X(modexy) - Box(1);
    if(std::fabs(dist1)<std::fabs(dist2))
    {
        if(dist1<2*R)
        {
            Ghost = true;
            double distb = Xb(modexy)-Box(0);
            X(modexy) = Box(1) + dist1 + 1.0;
            Xb(modexy) = Box(1) + distb + 1.0;
        }
        
    }else{
        if(dist2<2*R)
        {
            Ghost = true;
            double distb = Xb(modexy)-Box(1);
            X(modexy) = Box(0) + dist2 - 1.0;
            Xb(modexy) = Box(0) + distb - 1.0;

        }
        
}
    
      
}

inline void Disk::Leave(int modexy, Vec3_t &Box)
{
    if(IsFree())
    {
        if(X(modexy)<(Box(0)-1.5*R))
        {   
            // double dist = std::fabs(X(modexy)-Box(0));
            // X(modexy) = Box(1)-dist+1.0;        
            // double distb = Xb(modexy)-Box(0);
            // Xb(modexy) = Box(1)+distb+1.0;
            double dist = X(modexy)-Box(0);
            X(modexy) = Box(1)+dist+1.0;        
            double distb = Xb(modexy)-Box(0);
            Xb(modexy) = Box(1)+distb+1.0;

        }
        if(X(modexy)>Box(1)+1.5*R)
        {
            // double dist = std::fabs(X(modexy)-Box(1));
            // X(modexy) = Box(0)+dist-1.0;

            // double distb = Xb(modexy)-Box(1);
            // Xb(modexy) = Box(0)+distb-1.0;
            double dist = X(modexy) - Box(1);
            double distb = Xb(modexy)-Box(1);
            X(modexy) = Box(0) + dist - 1.0;
            Xb(modexy) = Box(0) + distb - 1.0;
        }
    }
    

}


inline void Disk::Translate(double dt)
{
    //std::cout << F(0) << " " << M << " " << V(0) << std::endl;
    // F = Fh+Fc+Ff+Flb;
    F = Fh+Fc+Ff;
    // F = Fc+Ff;
    Vec3_t Ft = F;
    if (vf(0)) Ft(0) = 0.0;
    if (vf(1)) Ft(1) = 0.0;
    if (vf(2)) Ft(2) = 0.0;
    //if (isnan(norm(F))) 
    //{
        //std::cout << Tag << std::endl;
    //}

    Vec3_t Xa = 2*X - Xb + Ft*(dt*dt/M);
    Vec3_t tp = Xa - X;
    // if(Tag==-396) std::cout<<tp(0)<<std::endl;
    V         = 0.5*(Xa - Xb)/dt;
    Xbt       = Xb;
    Xb        = X;
    X         = Xa;
}

inline void Disk::Rotate (double dt)
{
    double q0,q1,q2,q3,wx,wy,wz;
    q0 = 0.5*Q(0);
    q1 = 0.5*Q(1);
    q2 = 0.5*Q(2);
    q3 = 0.5*Q(3);
    T = Th+Tc+Tf;
    // T = Tc+Tf;
    Vec3_t Tt = T;
    if (wf(0)) Tt(0) = 0.0;
    if (wf(1)) Tt(1) = 0.0;
    if (wf(2)) Tt(2) = 0.0;

    Vec3_t Td = Tt/I;
    W = Wb+0.5*dt*Td;
    wx = W(0);
    wy = W(1);
    wz = W(2);
    Quaternion_t dq(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz),qm;

    Wb  = Wb+Td*dt;
    qm  = Q+dq*(0.5*dt);
    q0  = 0.5*qm(0);
    q1  = 0.5*qm(1);
    q2  = 0.5*qm(2);
    q3  = 0.5*qm(3);
    wx  = Wb(0);
    wy  = Wb(1);
    wz  = Wb(2);
    dq  = Quaternion_t(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz);
    Quaternion_t Qd = (qm+dq*0.5*dt),temp;
    Q  = Qd/norm(Qd);
}

}
#endif
