#include "./lbm/Domain.h"

int main (int argc, char **argv) try
{    
    Vec3_t w(0,0,2);
    Vec3_t t1;
    Quaternion_t Q(1.0,0,0,0);
    Rotation(w,Q,t1);
    return 0;

}MECHSYS_CATCH