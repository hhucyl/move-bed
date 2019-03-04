#include "./rw/particle.h"

int main (int argc, char **argv) try
{
    Vec3_t xt(0.5,0,0);
    RW::Particle RWP(xt);
    RWP.Xb = 1.5 ,0.0, 0.0;
    Vec3_t c(0,0,0);
    RWP.Reflect(c,1);
    std::cout<<RWP.X<<std::endl;
}MECHSYS_CATCH