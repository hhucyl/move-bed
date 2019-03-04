#ifndef RW_PARTICLE_H
#define RW_PARTICLE_H

#include <mechsys/linalg/matvec.h>


namespace RW
{
class Particle
{
    public:
    Particle(Vec3_t &X);
    Vec3_t X;
    Vec3_t Xb;
    void Move(std::vector<Vec3_t> &VV, std::vector<int> &idx, double dt);
    void Leave(int modexy, Vec3_t &Box);
    void Reflect(Vec3_t &C, double R);
};
inline Particle::Particle(Vec3_t &X0)
{
    X = X0;
    Xb = X;
}

inline void Particle::Move(std::vector<Vec3_t> &VV, std::vector<int> &idx, double dt)
{
    //VV 0=>11 1=>21 2=>12 3=>22
    //*12   *22
    //*11   *21
    //idx 0=>x1 1=>x2 2=>y1 3=>y2
    double x1 = (double) idx[0];
    double x2 = (double) idx[1];
    double y1 = (double) idx[2];
    double y2 = (double) idx[3];
    double x = X(0);
    double y = X(1);
    Vec3_t V(0,0,0);
    V = 1.0/((x2-x1)*(y2-y1))*(VV[0]*(x2-x)*(y2-y)+VV[1]*(x-x1)*(y2-y)+VV[2]*(x2-x)*(y-y1)+VV[3]*(x-x1)*(y-y1));
    // std::cout<<VV[0]<<" "<<VV[1]<<" "<<VV[2]<<" "<<VV[3]<<std::endl;
    // std::cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<std::endl;
    // std::cout<<V(0)<<"    "<<V(1)<<std::endl;
    Xb = X;
    X = X+V*dt;    
}

inline void Particle::Leave(int modexy, Vec3_t &Box)
{
    if(1000*X(modexy)<1000*Box(0))
    {   
        double dist = std::fabs(X(modexy)-Box(0));
        X(modexy) = Box(1)-dist;
        double distb = Xb(modexy)-Box(0);
        Xb(modexy) = Box(1)+distb;;   

    }
    if(X(modexy)>Box(1))
    {
        double dist = std::fabs(X(modexy)-Box(1));
        X(modexy) = Box(0)+dist;
        double distb = Xb(modexy)-Box(1);
        Xb(modexy) = Box(0)+distb;

    }
    


}

inline void Particle::Reflect(Vec3_t &C, double R)
{
    double L1 = Norm(X-Xb);
    double L2 = Norm(Xb-C);
    // std::cout<<"L1 = "<<L1<<" L2 ="<<L2<<std::endl;
    double AA = L1*L1;
    double BB = -2.0*((Xb(1)-C(1))*(Xb(1)-X(1))+(Xb(0)-C(0))*(Xb(0)-X(0)));
    double CC = L2*L2 - R*R;
    double Delta = BB*BB-4.0*AA*CC;
    // std::cout<<"Delta = "<<Delta<<std::endl;
    
    if(Delta>0)
    {
        double q;
        double q1 = (-BB+std::sqrt(Delta))/(2.0*AA);
        double q2 = (-BB-std::sqrt(Delta))/(2.0*AA);
        bool flag1 = q1>=0 && q1-1<1e-9;
        bool flag2 = q2>=0 && q2-1<1e-9;
        if(flag1)
        {
            q = q1;
        }else{
            if(flag2)
            {
                q = q2;
            }else{
                std::cout<<q2<<" "<<X<<" "<<Xb<<" "<<"ERROR IN RWPARTICLE REFLECT!!!"<<std::endl;
            }
        }
        // std::cout<<"q = "<<q<<std::endl;
        Vec3_t Xi(Xb(0)-q*(Xb(0)-X(0)),Xb(1)-q*(Xb(1)-X(1)),0.0);
        // std::cout<<"D = "<<Xi<<std::endl;
        double y = -((X(0)-2.*Xi(0))*(C(0)-Xi(0))*(C(1)-Xi(1)) + (X(1)-2.*Xi(1))*(C(1)-Xi(1))*(C(1)-Xi(1)) + X(0)*(C(0)-Xi(0))*(C(1)-Xi(1)) - X(1)*(C(0)-Xi(0))*(C(0)-Xi(0)))/(R*R);
        double x;
        if(std::fabs(C(1)-Xi(1))>0)
        {
            x = ((C(0)-Xi(0))*(y-X(1)))/(C(1)-Xi(1)) + X(0);
        }else{
            x = -((C(1)-Xi(1))*y + (X(0)-2.*Xi(0))*(C(0)-Xi(0)) + (X(1)-2.*Xi(1))*(C(1)-Xi(1)))/(C(0)-Xi(0));
        }
        // std::cout<<x<<" "<<y<<std::endl;
        X = x,y,0.0;
    }
    
}


}


#endif