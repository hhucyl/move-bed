

#ifndef LBM_DOMAIN_H
#define LBM_DOMAIN_H

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#ifdef USE_OCL
#include <mechsys/oclaux/cl.hpp>
#endif

// Std lib
#ifdef USE_OMP
#include <omp.h>
#endif

//STD
#include<iostream>
#include<set>
#include<iterator>
#include<utility>
#include<algorithm>

// Mechsys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/util/numstreams.h>

#include <mechsys/dem/special_functions.h>
#include "../dem/Disk.h"
#include "../dem/Interaction.h"

#include "../rw/particle.h"
enum LBMethod
{
    D2Q5,     ///< 2D 5 velocities
    D2Q9,     ///< 2D 9 velocities
    D3Q15,    ///< 3D 15 velocities
    D3Q19,    ///< 3D 19 velocities
    //D3Q27     ///< 3D 27 velocities
};

enum CollideMethod
{
    SRT,
    MRT
};

enum BounceBackMethod
{
    SBB,
    LIBB,
    QIBB,
    MR,
    CLI
};


namespace LBM
{



inline size_t Pt2idx(iVec3_t iv, iVec3_t & Dim) // Calculates the index of the cell at coordinates iv for a cubic lattice of dimensions Dim
{
    return iv(0) + iv(1)*Dim(0) + iv(2)*Dim(0)*Dim(1);
}

inline void   idx2Pt(size_t n, iVec3_t & iv, iVec3_t & Dim) // Calculates the coordinates from the index
{
    iv(0) = n%Dim(0);
    iv(1) = (n/Dim(0))%(Dim(1));
    iv(2) = n/(Dim(0)*Dim(1));
}



class Domain
{
public:
	static const double   WEIGHTSD2Q5   [ 5]; ///< Weights for the equilibrium distribution functions (D2Q5)
	static const double   WEIGHTSD2Q9   [ 9]; ///< Weights for the equilibrium distribution functions (D2Q9)
	static const double   WEIGHTSD3Q15  [15]; ///< Weights for the equilibrium distribution functions (D3Q15)
	static const double   WEIGHTSD3Q19  [19]; ///< Weights for the equilibrium distribution functions (D3Q19)
	//static const double   WEIGHTSD3Q27  [27]; ///< Weights for the equilibrium distribution functions (D3Q27)
	static const Vec3_t   LVELOCD2Q5    [ 5]; ///< Local velocities (D2Q5) 
	static const Vec3_t   LVELOCD2Q9    [ 9]; ///< Local velocities (D2Q9) 
	static const Vec3_t   LVELOCD3Q15   [15]; ///< Local velocities (D3Q15)
	static const Vec3_t   LVELOCD3Q19   [19]; ///< Local velocities (D3Q19)
	//static const Vec3_t   LVELOCD3Q27   [27]; ///< Local velocities (D3Q27)
	static const size_t   OPPOSITED2Q5  [ 5]; ///< Opposite directions (D2Q5) 
	static const size_t   OPPOSITED2Q9  [ 9]; ///< Opposite directions (D2Q9) 
	static const size_t   OPPOSITED3Q15 [15]; ///< Opposite directions (D3Q15)
	static const size_t   OPPOSITED3Q19 [19]; ///< Opposite directions (D3Q19)
	//static const size_t   OPPOSITED3Q27 [27]; ///< Opposite directions (D3Q27)
    static const double   MD2Q5       [5][5]; ///< MRT transformation matrix (D2Q5)
    static const double   MD2Q9       [9][9]; ///< MRT transformation matrix (D2Q9)
    static const double   MD3Q15    [15][15]; ///< MRT transformation matrix (D3Q15)
    static const double   MID3Q15    [15][15]; ///< MRT transformation matrix (D3Q15)
    static const double   MD3Q19    [19][19]; ///< MRT transformation matrix (D3Q19)
    //static const size_t   MD3Q27    [27][27]; ///< MRT transformation matrix (D3Q27)
    //static const double   SD2Q5          [5]; ///< MRT relaxation time vector (D2Q5)
    //static const double   SD2Q9          [9]; ///< MRT relaxation time vector (D2Q9)
    //static const double   SD3Q15        [15]; ///< MRT relaxation time vector (D3Q15)
    //static const double   SD3Q19        [19]; ///< MRT relaxation time vector (D3Q19)
    //static       double   SD3Q19        [27]; ///< MRT relaxation time vector (D3Q27)
    
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);
    typedef void (LBM::Domain::*ptFM_t)  (double *m, double rho, Vec3_t &vel); 
    typedef void (LBM::Domain::*ptFC_t)  ();
    typedef void (LBM::Domain::*ptFB_t)  (bool calcf);
    
   
    
    //Constructors
    

    //Special constructor with only one component, the parameters are the same as above
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    CollideMethod         MethodC, //< Type of collision
    double                nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step
    // ~Domain();
    //Methods
    
    void   CollideMRT();///< The collide step of LBM with MRT                                  
    void   CollideMRTMR();///< The collide step of LBM with MRT                                  
    void   CollideMRTLES();///< The collide step of LBM with MRT                        
    
    void   CollideSRTIBM();
    void   CollideMRTIBM();

    void   MeqD2Q9(double *m, double rho, Vec3_t &vel);                      
    void   MeqD3Q15(double *m, double rho, Vec3_t &vel);             
    void   MeqD3Q19(double *m, double rho, Vec3_t &vel);             
    
                           ///< The collide step of LBM with MRT
    void   CollideSRT();                                                           ///< The collide step of LBM for single component simulations
    void   Stream();
    void   BounceBack(bool calcF);
    void   BounceBackLIBB(bool calcF);
    void   BounceBackQIBB(bool calcF);
    void   BounceBackCLI(bool calcF);
    void   BounceBackMR(bool calcF);
    void   BoundaryGamma();
    double KernelIBM(double r, double x);
    double KernelIBM1(double r, double x);
    void   ApplyIBM2D(Vec3_t &pos, double R);
    void   ApplyIBM3D(Vec3_t &pos, double R);
    void   ApplyIBM3D();
    void   ApplyIBM3DIM(Vec3_t &pos, double R, size_t IT);
    void   GenPts(Vec3_t &pos, double R, int N);
    void   CalcProps();
    void   CalcPropsForCell(Vec3_t &idx);
    void   Initialize(iVec3_t idx, double Rho, Vec3_t & Vel);           ///< Initialize each cell with a given density and velocity
    double Feq(size_t k, double Rho, Vec3_t & Vel);                               ///< The equilibrium function
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);            ///< Solve the Domain dynamics
    
    //Writing Methods
    void WriteXDMF         (char const * FileKey);                                ///< Write the domain data in xdmf file
    void AddDiskQ(Vec3_t &pos, double R);
    void AddSphereQ(Vec3_t &pos, double R);
    void AddDiskG(Vec3_t &pos, double R);
    void AddDisksG();
    void adddiskG_sub(DEM::Disk *Pa);
    void AddSphereG(Vec3_t &pos, double R);
    void CalcForceQ();
    void SetZero();
    void StartSolve();
    void EndSolve();
    void SetBounceBack();
    void update_pair_sub(DEM::DiskPair &pair, DEM::Disk* P1, DEM::Disk* P2);
    void join_contactlist_sub(std::set<std::pair<int,int>> *myset_private, std::vector<std::pair<int,int>> &ListofContacts);
    void UpdateParticlesContacts();
    void UpdateParticlePairForce();
    void LeaveAndForcedForce();
    void GhostPeriodic();
    void MoveParticles();
    void Initial(double rho, Vec3_t &v0,  Vec3_t &g0);
    void InitialFromH5(char const * TheFileKey, Vec3_t &g0);


    void Solve(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport);
    void SolveRW(double Tf, double dtout, char const * TheFileKey, ptDFun_t ptSetup, ptDFun_t ptReport);
    void rwsolve_sub(double dt);
    void CheckInside();
    #ifdef USE_OMP
    omp_lock_t      lck;                      ///< to protect variables in multithreading
    #endif
    ptFM_t ptr2meq;
    ptFC_t ptr2collide;
    ptFB_t ptr2bb;
    //Data
    //time
    std::clock_t StartTime;
    std::clock_t EndTime;
    double ****q;
    double ****t;
    double ****Omeis;
    double **** F;                           ///< The array containing the individual functions with the order of the lattice, the x,y,z coordinates and the order of the function.
    double **** Ftemp;                       ///< A similar array to hold provitional data
    bool   ***  IsSolid;                     ///< An array of bools with an identifier to see if the cell is a solid cell
    Vec3_t ***  Vel;                         ///< The fluid velocities
    Vec3_t ***  BForce;                      ///< Body Force for each cell
    double ***  Rho;                         ///< The fluid densities
    double      Tau;                         ///< The characteristic time of the lattice
                          ///< The mixing constant for multicomponent simulations
    size_t const * Op;                        ///< An array containing the indexes of the opposite direction for bounce back conditions
    double const *  W;                        ///< An array with the direction weights
    double *     EEk;                         ///< Diadic product of the velocity vectors
    Vec3_t const *  C;                        ///< The array of lattice velocities
    Mat_t        M;                           ///< Transformation matrix to momentum space for MRT calculations
    Mat_t        Minv;                        ///< Inverse Transformation matrix to momentum space for MRT calculations
    Vec_t        S;                               ///< Vector of relaxation times for MRT
    size_t       Nneigh;                        ///< Number of Neighbors, depends on the scheme
    double       dt;                          ///< Time Step
    double       dx;                          ///< Grid size
    double       Cs;                          ///< Lattice Velocity
    bool         IsFirstTime;                 ///< Bool variable checking if it is the first time function Setup is called
    iVec3_t      Ndim;                        ///< Lattice Dimensions
    size_t       Ncells;                      ///< Number of cells
    size_t       Nproc;                       ///< Number of processors for openmp
    size_t       idx_out;                     ///< The discrete time step for output
    String       FileKey;                     ///< File Key for output files
    void *       UserData;                    ///< User Data
    size_t       Step;                        ///< Lenght of averaging cube to save data
    double       Time;                        ///< Simulation time variable
    size_t       Nl;                          ///< Number of lattices (fluids)
    double       Sc;                          ///< Smagorinsky constant
    double       Nu;
    double       Rho0;
    LBMethod     Method;
    CollideMethod MethodC;
    BounceBackMethod MethodB;
    bool IsMR;
    bool IsF;
    bool IsFt;
    bool Isq;
    Vec3_t ***VelP;
    Vec3_t ***Flbm;
    double ***Gamma;
    int ***Check;
    // int ***CheckRh;
    double we;
    double wej;
    double wxx;
    std::vector<Vec3_t> points; ///IBM
    std::vector<double> dS;     ///IBM
    // std::vector<LBM::Disk> Particles; ///DEM
    std::vector<DEM::Disk> Particles; ///DEM
    std::vector<DEM::Disk> GhostParticles; ///DEM
    std::vector<std::pair<int,int>> ListofContacts;
    // std::vector<std::pair<int,int>> ListofContactsPP;
    // std::vector<std::pair<int,int>> ListofContactsPG;
    std::vector<RW::Particle> RWParticles;

    double dtdem;
    Vec3_t Box;
    int modexy;

    
};





inline Domain::Domain(LBMethod TheMethod, CollideMethod TheMethodC,  double Thenu, iVec3_t TheNdim, double Thedx, double Thedt)
{
    Array<double> nu(1);
    nu[0] = Thenu;
    
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    if (nu.Size()==0) throw new Fatal("LBM::Domain: Declare at leat one fluid please");
    if (TheNdim(2) >1&&(TheMethod==D2Q9 ||TheMethod==D2Q5 ))  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (TheNdim(2)==1&&(TheMethod==D3Q15||TheMethod==D3Q19))  throw new Fatal("LBM::Domain: Ndim(2) is 1. Either change the method to D2Q9 or increase the z-dimension");
    IsF = false;
    IsFt = false;
    Isq = false;
    IsMR = false;
    ptr2meq = NULL;
    ptr2collide = NULL;
    ptr2bb = NULL;
    Time        = 0.0;
    dt          = Thedt;
    dx          = Thedx;
    Cs          = dx/dt;
    Step        = 1;
    Sc          = 0.17;
    Nl          = 1;
    Ndim        = TheNdim;
    Rho0        = 1.0;
    Ncells      = Ndim(0)*Ndim(1)*Ndim(2);
    IsFirstTime = true;
    Method = TheMethod;
    MethodC = TheMethodC;
    Nu = nu[0];
    Box = 0.0, 0.0, 0.0;
    modexy = -1;
    Tau         = 3.0*nu[0]*dt/(dx*dx)+0.5;
    if (TheMethod==D2Q5)
    {
        Nneigh = 5;
        W      = WEIGHTSD2Q5;
        C      = LVELOCD2Q5;
        Op     = OPPOSITED2Q5;
    }
    if (TheMethod==D2Q9)
    {
        Nneigh = 9;
        W      = WEIGHTSD2Q9;
        C      = LVELOCD2Q9;
        Op     = OPPOSITED2Q9;
        ptr2meq = &LBM::Domain::MeqD2Q9;
        ptr2collide = &LBM::Domain::CollideSRT; 
        if(TheMethodC == MRT)
        {
            M.Resize(Nneigh,Nneigh);
            Minv.Resize(Nneigh,Nneigh);
            for (size_t n=0;n<Nneigh;n++)
            for (size_t m=0;m<Nneigh;m++)
            {
                M(n,m) = MD2Q9[n][m];
                // Minv(n,m) = MID3Q15[n][m];
            }
            Inv(M,Minv);
            
            double s = 1.0/Tau;
            S.Resize(Nneigh);
            S = 1.0,1.4,1.4,1.0,1.2,1.0,1.2,s,s;
            ptr2collide = &LBM::Domain::CollideMRT;
            
        }
       
        
    }
    if (TheMethod==D3Q15)
    {
        Nneigh = 15;
        W      = WEIGHTSD3Q15;
        C      = LVELOCD3Q15;
        Op     = OPPOSITED3Q15;
        ptr2meq = &LBM::Domain::MeqD3Q15;
        ptr2collide = &LBM::Domain::CollideSRT; 
        if(TheMethodC == MRT)
        {
            M.Resize(Nneigh,Nneigh);
            Minv.Resize(Nneigh,Nneigh);
            for (size_t n=0;n<Nneigh;n++)
            for (size_t m=0;m<Nneigh;m++)
            {
                M(n,m) = MD3Q15[n][m];
                // Minv(n,m) = MID3Q15[n][m];
            }
            Inv(M,Minv);
            double tau = 3.0*Thenu*dt/(dx*dx)+0.5;
            double s   = 8.0*(2.0-1.0/tau)/(8.0-1.0/tau);
            S.Resize(Nneigh);
            S = 0.0,1.0/tau,1.0/tau,0.0,s,0.0,s,0.0,s,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,s;
            // S = 0.0,1.6,1.2,0.0,1.6,0.0,1.6,0.0,1.6,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.2;
            ptr2collide = &LBM::Domain::CollideMRT;
            
        }
        
        
    }
    if (TheMethod==D3Q19)
    {
        Nneigh = 19;
        W      = WEIGHTSD3Q19;
        C      = LVELOCD3Q19;
        Op     = OPPOSITED3Q19;
        ptr2meq = &LBM::Domain::MeqD3Q19;
        ptr2collide = &LBM::Domain::CollideSRT; 
        if(TheMethodC == MRT)
        {
            M.Resize(Nneigh,Nneigh);
            Minv.Resize(Nneigh,Nneigh);
            for (size_t n=0;n<Nneigh;n++)
            for (size_t m=0;m<Nneigh;m++)
            {
                M(n,m) = MD3Q19[n][m];
                // Minv(n,m) = MID3Q15[n][m];
            }
            Inv(M,Minv);
            double tau = 3.0*Thenu*dt/(dx*dx)+0.5;
            double s   = 8.0*(2.0-1.0/tau)/(8.0-1.0/tau);
            S.Resize(Nneigh);
            // S = 0.0,1.0/tau,1.0/tau,0.0,s,0.0,s,0.0,s,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,s;
            S = 0.0,s,s,0.0,s,0.0,s,0.0,s,1.0/tau,s,1.0/tau,s,1.0/tau,1.0/tau,1.0/tau,s,s,s;
            ptr2collide = &LBM::Domain::CollideMRT;
            we = 3.0;
            wej = -11.0/2.0;
            wxx = -0.5;
        }
       
        
    }




    t           = new double *** [Ndim(0)];    
    q           = new double *** [Ndim(0)];
	Omeis       = new double *** [Ndim(0)];
    VelP        = new Vec3_t **  [Ndim(0)];
    Flbm        = new Vec3_t **  [Ndim(0)];
    F           = new double *** [Ndim(0)];
    Ftemp       = new double *** [Ndim(0)];
    Vel         = new Vec3_t **  [Ndim(0)];
    BForce      = new Vec3_t **  [Ndim(0)];
    Rho         = new double **  [Ndim(0)];
    IsSolid     = new bool   **  [Ndim(0)];
    Gamma       = new double **  [Ndim(0)];
    Check       = new int **  [Ndim(0)];
    // CheckRh       = new int **  [Ndim(0)];
    for (size_t nx=0;nx<Ndim(0);nx++)
    {
        F       [nx]    = new double ** [Ndim(1)];
        Ftemp   [nx]    = new double ** [Ndim(1)];
        Vel     [nx]    = new Vec3_t *  [Ndim(1)];
        BForce  [nx]    = new Vec3_t *  [Ndim(1)];
        Rho     [nx]    = new double *  [Ndim(1)];
        IsSolid [nx]    = new bool   *  [Ndim(1)];
        q       [nx]    = new double ** [Ndim(1)];
        t       [nx]    = new double ** [Ndim(1)];
	    Omeis   [nx]    = new double ** [Ndim(1)];
        VelP    [nx]    = new Vec3_t *  [Ndim(1)];
        Flbm    [nx]    = new Vec3_t *  [Ndim(1)];
        Gamma   [nx]    = new double *  [Ndim(1)];
        Check   [nx]    = new int *  [Ndim(1)];
        // CheckRh   [nx]    = new int *  [Ndim(1)];
        for (size_t ny=0;ny<Ndim(1);ny++)
        {
            F       [nx][ny]    = new double * [Ndim(2)];
            Ftemp   [nx][ny]    = new double * [Ndim(2)];
            Vel     [nx][ny]    = new Vec3_t   [Ndim(2)];
            BForce  [nx][ny]    = new Vec3_t   [Ndim(2)];
            Rho     [nx][ny]    = new double   [Ndim(2)];
            IsSolid [nx][ny]    = new bool     [Ndim(2)];
            q       [nx][ny]    = new double * [Ndim(2)];
            t       [nx][ny]    = new double * [Ndim(2)];
	        Omeis   [nx][ny]    = new double * [Ndim(2)];
            VelP    [nx][ny]    = new Vec3_t   [Ndim(2)];
            Flbm    [nx][ny]    = new Vec3_t   [Ndim(2)];
            Gamma   [nx][ny]    = new double   [Ndim(2)];        
            Check   [nx][ny]    = new int   [Ndim(2)];        
            // CheckRh   [nx][ny]    = new int   [Ndim(2)];        
            for (size_t nz=0;nz<Ndim(2);nz++)
            {
                F    [nx][ny][nz]    = new double [Nneigh];
                Ftemp[nx][ny][nz]    = new double [Nneigh];
                q    [nx][ny][nz]    = new double [Nneigh];
                t    [nx][ny][nz]    = new double [Nneigh];
	            Omeis[nx][ny][nz]    = new double [Nneigh];
                IsSolid[nx][ny][nz]  = false;
                VelP[nx][ny][nz] = 0.0, 0.0, 0.0;
                Flbm[nx][ny][nz] = 0.0, 0.0, 0.0;
                BForce[nx][ny][nz] = 0.0, 0.0, 0.0;
                Gamma[nx][ny][nz] = 0.0;
                Check[nx][ny][nz] = -1;
                // CheckRh[nx][ny][nz] = -1;
                for (size_t nn=0;nn<Nneigh;nn++)
                {
                    F    [nx][ny][nz][nn] = 0.0;
                    Ftemp[nx][ny][nz][nn] = 0.0;
                    q    [nx][ny][nz][nn] = -1.0;
                    t    [nx][ny][nz][nn] = 0.0;
                    Omeis[nx][ny][nz][nn] = 0.0;
                }
            }
        }
    }
    


    EEk = new double [Nneigh];
    for (size_t k=0;k<Nneigh;k++)
    {
        EEk[k]    = 0.0;
        for (size_t n=0;n<3;n++)
        for (size_t m=0;m<3;m++)
        {
            EEk[k] += fabs(C[k][n]*C[k][m]);
        }
    }

    

    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Nl*Ncells,TERM_RST);
#ifdef USE_OMP
    omp_init_lock(&lck);
#endif
}

#include "Feq.h"
#include "Output.h"
#include "Collide.h"
#include "Stream.h"
#include "CalcProps.h"
#include "Initialize.h"
#include "ApplySolidBoundary.h"
#include "CoupleDEM.h"
#include "solve.h"


inline void Domain::SetBounceBack()
{
    if(MethodB == SBB)
    {
        ptr2bb = &Domain::BounceBack;
    }else if(MethodB == LIBB){
        ptr2bb = &Domain::BounceBackLIBB;        
    }else if(MethodB == QIBB){
        ptr2bb = &Domain::BounceBackQIBB;
    }else if(MethodB == CLI){
        ptr2bb = &Domain::BounceBackCLI;
    }else if(MethodB == MR){
        ptr2bb = &Domain::BounceBackMR;         
    }else{
        throw new Fatal("Collide Type is NOT RIGHT!!!!!");    
    }
}

inline void Domain::SetZero()
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
        Gamma[ix][iy][iz] = 0.0;
        Check[ix][iy][iz] = -1;
        // CheckRh[ix][iy][iz] = -1;
        VelP[ix][iy][iz] = 0.0, 0.0, 0.0;
        Flbm[ix][iy][iz] = 0.0, 0.0, 0.0;
        for(size_t k=0; k<Nneigh; ++k)
        {
            Omeis[ix][iy][iz][k] = 0.0;
        }
    }
}

inline void Domain::StartSolve()
{
    if(ptr2collide == NULL)
    {
        throw new Fatal("Collsion is not initialized, check the pointer");
    }else{
        if(MethodC == MRT)
        {
            if(ptr2meq == NULL)
                throw new Fatal("Matrix for MRT is not initialized, check the pointer");
            
        }
    }
    // printf("\033[01;33m--- Solving ---\033[0m\n");
    printf("--- Tau = %g ---\n",Tau);
    // StartTime = std::clock();
    idx_out = 0;
}
inline void Domain::EndSolve()
{
    // using namespace std;
    // EndTime = std::clock();
    // double ttime = (double)(EndTime - StartTime)/CLOCKS_PER_SEC;
    // printf("\033[01;34m---Elapsed Time = %f s---\033[0m\n",ttime);
    
}

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{
    


    idx_out     = 0;
    FileKey.Printf("%s",TheFileKey);
    Nproc = TheNproc;
    
    if(ptr2collide == NULL)
    {
        throw new Fatal("Collsion is not initialized, check the pointer");
    }else{
        if(MethodC == MRT)
        {
            if(ptr2meq == NULL)
                throw new Fatal("Matrix for MRT is not initialized, check the pointer");
            
        }
    }

    
    // if()
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving -------------------------------------------------------------------%s\n",TERM_CLR1    , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                    , TERM_RST);
    
    printf("%s  Tau of Lattice                  =  %g%s\n"       ,TERM_CLR2, Tau                             , TERM_RST);
    
  
    

    
    
    //std::cout << "2" << std::endl;

    double tout = Time;
    while (Time < Tf)
    {
        if (Time >= tout)
        {
            //std::cout << "3" << std::endl;
            
            //std::cout << "4" << std::endl;
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                if ( RenderVideo) 
                {
                    #ifdef USE_HDF5
                    WriteXDMF(fn.CStr());
                    #else
                    //WriteVTK (fn.CStr());
                    #endif
                }
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
                std::cout<<"-------------Time = "<<Time<<" ----------------"<<std::endl;
            }
            tout += dtOut;
            idx_out++;
        }

        //SetZero();
        Vec3_t pos(49.5,49.5,49.5);
        AddSphereG(pos,5.0);
        
        // CollideSRT(ptr2meq);
        CollideSRT();
        // (this->*ptr2collide)(ptr2meq);
        
        Stream();
        // (this->*ptr2bounceback)();
        // BounceBackLIBB();
        // CalcForceQ();
        // if (ptr2calcforce!=NULL)(this->*ptr2calcforce)();
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        
        CalcProps();
        
        

        Time += dt;
        //std::cout << Time << std::endl;
    }
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);
    
}

const double Domain::WEIGHTSD2Q5   [ 5] = { 2./6., 1./6., 1./6., 1./6., 1./6 };
const double Domain::WEIGHTSD2Q9   [ 9] = { 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36. };
const double Domain::WEIGHTSD3Q15  [15] = { 2./9., 1./9., 1./9., 1./9., 1./9.,  1./9.,  1./9., 1./72., 1./72. , 1./72., 1./72., 1./72., 1./72., 1./72., 1./72.};
const double Domain::WEIGHTSD3Q19  [19] = { 1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};
const size_t Domain::OPPOSITED2Q5  [ 5] = { 0, 3, 4, 1, 2 };                                                       ///< Opposite directions (D2Q5) 
const size_t Domain::OPPOSITED2Q9  [ 9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };                                           ///< Opposite directions (D2Q9) 
const size_t Domain::OPPOSITED3Q15 [15] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};                     ///< Opposite directions (D3Q15)
const size_t Domain::OPPOSITED3Q19 [19] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};     ///< Opposite directions (D3Q19)
const Vec3_t Domain::LVELOCD2Q5  [ 5] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0} };
const Vec3_t Domain::LVELOCD2Q9  [ 9] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}, {1,1,0}, {-1,1,0}, {-1,-1,0}, {1,-1,0} };
const Vec3_t Domain::LVELOCD3Q15 [15] =
{
	{ 0, 0, 0}, { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, 
	{ 0, 0, 1}, { 0, 0,-1}, { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, 
	{-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} 
};

const Vec3_t Domain::LVELOCD3Q19 [19] =
{
	{ 0, 0, 0}, 
    { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1}, 
    { 1, 1, 0}, {-1,-1, 0}, { 1, 0, 1}, {-1, 0,-1}, { 0, 1, 1}, { 0,-1,-1},
    { 1,-1, 0}, {-1, 1, 0}, { 1, 0,-1}, {-1, 0, 1}, { 0, 1,-1}, { 0,-1, 1}
};

const double   Domain::MD2Q9 [ 9][ 9] = { { 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0},
                                          {-4.0, -1.0, -1.0, -1.0, -1.0,  2.0,  2.0,  2.0,  2.0},{ 4.0, -2.0, -2.0, -2.0, -2.0,  1.0,  1.0,  1.0,  1.0},
                                          { 0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0},
                                          { 0.0, -2.0,  0.0,  2.0,  0.0,  1.0, -1.0, -1.0,  1.0},
                                          { 0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0},
                                          { 0.0,  0.0, -2.0,  0.0,  2.0,  1.0,  1.0, -1.0, -1.0},
                                          { 0.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0},
                                          { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0}};

const double Domain::MD3Q15 [15][15]  =  /*{ { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                            {-2.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                            {16.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                            { 0.0, 1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
                                            { 0.0,-4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
                                            { 0.0, 0.0, 0.0, 1.0,-1.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0},
                                            { 0.0, 0.0, 0.0,-4.0, 4.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0,-4.0, 4.0, 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0},
                                            { 0.0, 2.0, 2.0,-1.0,-1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                            { 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0},
                                            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0, 1.0,-1.0} };*/
                                        { { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                          {-2.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                          {16.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                          { 0.0, 1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
                                          { 0.0,-4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
                                          { 0.0, 0.0, 0.0, 1.0,-1.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0},
                                          { 0.0, 0.0, 0.0,-4.0, 4.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0,-4.0, 4.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0},
                                          { 0.0, 2.0, 2.0,-1.0,-1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                          { 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0, 1.0,-1.0} };

const double Domain::MD3Q19 [19][19]= { {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                        {-30.0, -11.0, -11.0, -11.0, -11.0, -11.0, -11.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0},
                                        {12.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                        {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0},
                                        {0.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 1.0, -1.0},
                                        {0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 1.0, -1.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, -1.0, 1.0, -1.0, 1.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, -1.0, 1.0, -1.0, 1.0},
                                        {0.0, 2.0, 2.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0},
                                        {0.0, -4.0, -4.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0},
                                        {0.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, -2.0, -2.0, 2.0, 2.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, -1.0, 1.0, 1.0, -1.0}};

}
#endif
