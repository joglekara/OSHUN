/*! \brief time-integrators for various methods
 * \author PICKSC
 * \file   steppers.cpp
 *
 * Includes spatial advection, electric field advection, and electric field update routines
 *
 */
//--------------------------------------------------------------
//  Standard libraries
#include <mpi.h>
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>
#include <algorithm>
#include <cstdlib>
#include <mpi.h>

#include <math.h>
#include <map>

//  My libraries
#include "lib-array.h"
#include "lib-algorithms.h"

//  Declerations
#include "gpu.h"
#include "state.h"
#include "input.h"
#include "formulary.h"
#include "collisions.h"
#include "parallel.h"
#include "vlasov.h"
#include "setup.h"
#include "functors.h"
#include "stepper.h"


//**************************************************************
//**************************************************************
//--------------------------------------------------------------
ARK32::ARK32(State1D& Yin): Yhv1(Yin), Yhv2(Yin), Yhv3(Yin), Yhv4(Yin),
                        Yhc1(Yin),Yhc2(Yin), Yhc3(Yin), Yhc4(Yin), 
                        Yt(Yin),
        ae21(1767732205903./2027836641118.),
        ae31(5535828885825./10492691773637.), ae32(788022342437./10882634858940.),
        ae41(6485989280629./16251701735622.), ae42(-4246266847089./9704473918619), ae43(10755448449292./10357097424841),
        ai21(1767732205903./4055673282236), ai22(1767732205903./4055673282236.),
        ai31(2746238789719./10658868560708.), ai32(-640167445237./6845629431997), ai33(1767732205903./4055673282236.),
        ai41(1471266399579./7840856788654.), ai42(-4482444167858./7529755066697), ai43(11266239266428./11593286722821.), ai44(1767732205903./4055673282236.),

        b1(1471266399579./7840856788654.), b2(-4482444167858./7529755066697), b3(11266239266428./11593286722821.), b4(1767732205903./4055673282236.),
        b1_LO(2756255671327./12835298489170.), b2_LO(-10771552573575./22201958757719), b3_LO(9247589265047./10645013368117), b4_LO(2193209047091./5459859503100)
    {
        
        Yhv1 = static_cast<complex<double> >(0.0);  Yhc1 = static_cast<complex<double> >(0.0);
        Yhv2 = static_cast<complex<double> >(0.0);  Yhc2 = static_cast<complex<double> >(0.0);
        Yhv3 = static_cast<complex<double> >(0.0);  Yhc3 = static_cast<complex<double> >(0.0);
        Yhv4 = static_cast<complex<double> >(0.0);  Yhc4 = static_cast<complex<double> >(0.0);
    }
//--------------------------------------------------------------
ARK32:: ~ARK32(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
}
void ARK32::take_step(State1D& Y2, State1D& Y3, double time, double h, VlasovFunctor1D_explicitE& vF, collisions_1D& coll, Parallel_Environment_1D& PE) 
{
//      Take a step using ARK32

//      Yh1, Stage 1
        // z1 = Y2;

        vF(Y3,Yhv1); PE.Neighbor_Communications(Yhv1);
        Yhv1 *= ae21*h; Yhc1 *= ai21*h;
        
        Yt = Y3;
        Yt += Yhv1; Yt += Yhc1;     
        
        coll(Yt,Yhc2,time,(ai22*h));
        Yhc2 *= (ai22*h);    Yt += Yhc2;
        
        // z2 = Yt;

        vF(Yt,Yhv2); PE.Neighbor_Communications(Yhv2);
        Yhv1 *= ae31/ae21;  Yhc1 *= ai31/ai21;  
        Yhv2 *= ae32*h;     Yhc2 *= ai32/ai22;  

        Yt = Y2;
        Yt += Yhv1; Yt += Yhc1;
        Yt += Yhv2; Yt += Yhc2;
        
        coll(Yt,Yhc3,time,(ai33*h));
        Yhc3 *= (ai33*h);   Yt += Yhc3;
        
        // z3 = Yt;

        vF(Yt,Yhv3); PE.Neighbor_Communications(Yhv3);
        Yhv1 *= ae41/ae31;  Yhc1 *= ai41/ai31;  
        Yhv2 *= ae42/ae32;  Yhc2 *= ai42/ai32;  
        Yhv3 *= ae43*h;     Yhc3 *= ai43/ai33;

        Yt = Y2;
        Yt += Yhv1; Yt += Yhc1;
        Yt += Yhv2; Yt += Yhc2;
        Yt += Yhv3; Yt += Yhc3;

        coll(Yt,Yhc4,time,(ai44*h));
        Yhc4 *= (ai44*h);   Yt += Yhc4;
        
        // z4 = Yt;
        vF(Yt,Yhv4); PE.Neighbor_Communications(Yhv4);

        //  Assemble 3rd order solution
        Y2 = Y3;
        Yhv1 *= b1/ae41;    Y3 += Yhv1;
        Yhv2 *= b2/ae42;    Y3 += Yhv2;
        Yhv3 *= b3/ae43;    Y3 += Yhv3;
        Yhv4 *= b4*h;       Y3 += Yhv4;

        Yhc1 *= b1/ai41;    Y3 += Yhc1;
        Yhc2 *= b2/ai42;    Y3 += Yhc2;
        Yhc3 *= b3/ai43;    Y3 += Yhc3;
        Yhc4 *= b4/ai44;    Y3 += Yhc4;

        //  Assemble 2nd order solution
        Yhv1 *= b1_LO/b1;   Y2 += Yhv1;
        Yhv2 *= b2_LO/b2;   Y2 += Yhv2;
        Yhv3 *= b3_LO/b3;   Y2 += Yhv3;
        Yhv4 *= b4_LO/b4;   Y2 += Yhv4;

        Yhc1 *= b1_LO/b1;   Y2 += Yhc1;
        Yhc2 *= b2_LO/b2;   Y2 += Yhc2;
        Yhc3 *= b3_LO/b3;   Y2 += Yhc3;
        Yhc4 *= b4_LO/b4;   Y2 += Yhc4;

        Yhc4 *= 1./(b4_LO*h);
        Yhc1 = Yhc4;

//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
//--------------------------------------------------------------
ARK43::ARK43(State1D& Yin): Yhv1(Yin), Yhv2(Yin), Yhv3(Yin), Yhv4(Yin), Yhv5(Yin), Yhv6(Yin),
                        Yhc1(Yin), Yhc2(Yin), Yhc3(Yin), Yhc4(Yin), Yhc5(Yin), Yhc6(Yin), 
                        Yt(Yin),
        
        ae21(0.5),
        ae31(13861.0/62500.0), ae32(6889.0/62500.0),
        ae41(-116923316275.0/2393684061468.0), ae42(-2731218467317.0/15368042101831.0), ae43(9408046702089.0/11113171139209.0),
        ae51(-451086348788.0/2902428689909.0), ae52(-2682348792572.0/7519795681897.0), ae53(12662868775082.0/11960479115383.0), ae54(3355817975965.0/11060851509271.0),
        ae61(647845179188.0/3216320057751.0), ae62(73281519250.0/8382639484533.0), ae63(552539513391.0/3454668386233.0), ae64(3354512671639.0/8306763924573.0), ae65(4040.0/17871.0),
        
        ai21(1.0/4.0), ai22(1.0/4.0),
        ai31(8611.0/62500.0), ai32(-1743.0/31250.0), ai33(1.0/4.0),
        ai41(5012029.0/34652500.0), ai42(-654441.0/2922500.0), ai43(174375.0/388108.0), ai44(1.0/4.0),
        ai51(15267082809.0/155376265600.0), ai52(-71443401.0/120774400.0), ai53(730878875.0/902184768.0), ai54(2285395.0/8070912.0), ai55(1.0/4.0),
        ai61(82889.0/524892.0), ai63(15625.0/83664.0), ai64(69875.0/102672.0), ai65(-2260.0/8211.0), ai66(1.0/4.0),

        b1(82889.0/524892.0), b3(15625.0/83664.0), b4(69875.0/102672.0), b5(-2260.0/8211.0), b6(1.0/4.0),
        b1_LO(4586570599.0/29645900160.0), b3_LO(178811875.0/945068544.0), b4_LO(814220225.0/1159782912.0), b5_LO(-3700637.0/11593932.0), b6_LO(61727.0/225920.0)
    {
        
        Yhv1 = static_cast<complex<double> >(0.0);  Yhc1 = static_cast<complex<double> >(0.0);
        Yhv2 = static_cast<complex<double> >(0.0);  Yhc2 = static_cast<complex<double> >(0.0);
        Yhv3 = static_cast<complex<double> >(0.0);  Yhc3 = static_cast<complex<double> >(0.0);
        Yhv4 = static_cast<complex<double> >(0.0);  Yhc4 = static_cast<complex<double> >(0.0);
        Yhv5 = static_cast<complex<double> >(0.0);  Yhc5 = static_cast<complex<double> >(0.0);
        Yhv6 = static_cast<complex<double> >(0.0);  Yhc6 = static_cast<complex<double> >(0.0);
    }
//--------------------------------------------------------------
ARK43:: ~ARK43(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
}
void ARK43::take_step(State1D& Y3, State1D& Y4, double time, double h, VlasovFunctor1D_explicitE& vF, collisions_1D& coll, Parallel_Environment_1D& PE) 
{
//      Take a step using ARK43

//      Yh1, Stage 1
    // z1 = Y2;

    vF(Y4,Yhv1); PE.Neighbor_Communications(Yhv1);
    Yhv1 *= ae21*h; Yhc1 *= ai21*h;
    
    Yt = Y4;
    Yt += Yhv1; Yt += Yhc1;     

    coll(Yt,Yhc2,time,(ai22*h));   
    Yhc2 *= (ai22*h);    Yt += Yhc2;
    
    // z2 = Yt;

    vF(Yt,Yhv2); PE.Neighbor_Communications(Yhv2);
    Yhv1 *= ae31/ae21;  Yhc1 *= ai31/ai21;  
    Yhv2 *= ae32*h;     Yhc2 *= ai32/ai22;  

    Yt = Y4;
    Yt += Yhv1; Yt += Yhc1;
    Yt += Yhv2; Yt += Yhc2;
    coll(Yt,Yhc3,time,(ai33*h));

    Yhc3 *= (ai33*h);   Yt += Yhc3;
    
    // z3 = Yt;

    vF(Yt,Yhv3); PE.Neighbor_Communications(Yhv3);
    Yhv1 *= ae41/ae31;  Yhc1 *= ai41/ai31;  
    Yhv2 *= ae42/ae32;  Yhc2 *= ai42/ai32;  
    Yhv3 *= ae43*h;     Yhc3 *= ai43/ai33;

    Yt = Y4;
    Yt += Yhv1; Yt += Yhc1;
    Yt += Yhv2; Yt += Yhc2;
    Yt += Yhv3; Yt += Yhc3;

    coll(Yt,Yhc4,time,(ai44*h));

    Yhc4 *= (ai44*h);   Yt += Yhc4;
    
    // z4 = Yt;
    vF(Yt,Yhv4); PE.Neighbor_Communications(Yhv4);
    Yhv1 *= ae51/ae41;  Yhc1 *= ai51/ai41;  
    Yhv2 *= ae52/ae42;  Yhc2 *= ai52/ai42;  
    Yhv3 *= ae53/ae43;  Yhc3 *= ai53/ai43;  
    Yhv4 *= ae54*h;     Yhc4 *= ai54/ai44;

    Yt = Y4;
    Yt += Yhv1; Yt += Yhc1;
    Yt += Yhv2; Yt += Yhc2;
    Yt += Yhv3; Yt += Yhc3;
    Yt += Yhv4; Yt += Yhc4;

    coll(Yt,Yhc5,time,(ai55*h));

    Yhc5 *= (ai55*h);   Yt += Yhc5;
    
    // z5 = Yt;
    vF(Yt,Yhv5); PE.Neighbor_Communications(Yhv5);
    Yhv1 *= ae61/ae51;  Yhc1 *= ai61/ai51;  
    Yhv2 *= ae62/ae52;  Yhc2 *= ai62/ai52;  
    Yhv3 *= ae63/ae53;  Yhc3 *= ai63/ai53;  
    Yhv4 *= ae64/ae54;  Yhc4 *= ai64/ai54;  
    Yhv5 *= ae65*h;     Yhc5 *= ai65/ai55;

    Yt = Y4;
    Yt += Yhv1; Yt += Yhc1;
    Yt += Yhv2; Yt += Yhc2;
    Yt += Yhv3; Yt += Yhc3;
    Yt += Yhv4; Yt += Yhc4;
    Yt += Yhv5; Yt += Yhc5;

    coll(Yt,Yhc6,time,(ai66*h));

    Yhc6 *= (ai66*h);   Yt += Yhc6;

    // z6 = Yt;
    vF(Yt,Yhv6);  PE.Neighbor_Communications(Yhv6);

    //  Assemble 4th order solution
    Y3 = Y4;
    Yhv1 *= b1/ae61;    Y4 += Yhv1;
    Yhv3 *= b3/ae63;    Y4 += Yhv3;
    Yhv4 *= b4/ae64;    Y4 += Yhv4;
    Yhv5 *= b5/ae65;    Y4 += Yhv5;
    Yhv6 *= b6*h;       Y4 += Yhv6;

    Yhc1 *= b1/ai61;    Y4 += Yhc1;
    Yhc3 *= b3/ai63;    Y4 += Yhc3;
    Yhc4 *= b4/ai64;    Y4 += Yhc4;
    Yhc5 *= b5/ai65;    Y4 += Yhc5;
    Yhc6 *= b6/ai66;    Y4 += Yhc6;

    //  Assemble 3rd order solution
    Yhv1 *= b1_LO/b1;   Y3 += Yhv1;
    Yhv3 *= b3_LO/b3;   Y3 += Yhv3;
    Yhv4 *= b4_LO/b4;   Y3 += Yhv4;
    Yhv5 *= b5_LO/b5;   Y3 += Yhv5;
    Yhv6 *= b6_LO/b6;   Y3 += Yhv6;

    Yhc1 *= b1_LO/b1;   Y3 += Yhc1;
    Yhc3 *= b3_LO/b3;   Y3 += Yhc3;
    Yhc4 *= b4_LO/b4;   Y3 += Yhc4;
    Yhc5 *= b5_LO/b5;   Y3 += Yhc5;
    Yhc6 *= b6_LO/b6;   Y3 += Yhc6;

    Yhc6 *= 1./(b6_LO*h);
    Yhc1 = Yhc6;

//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
//--------------------------------------------------------------
ARK54::ARK54(State1D& Yin): Yhv1(Yin), Yhv2(Yin), Yhv3(Yin), Yhv4(Yin), Yhv5(Yin), Yhv6(Yin), Yhv7(Yin), Yhv8(Yin),
                        Yhc1(Yin), Yhc2(Yin), Yhc3(Yin), Yhc4(Yin), Yhc5(Yin), Yhc6(Yin), Yhc7(Yin), Yhc8(Yin),
                        Yt(Yin),

    ae21(41.0/100.0),
    ae31(367902744464.0/2072280473677.0), ae32(677623207551.0/8224143866563.0),
    ae41(1268023523408.0/10340822734521.0), ae43(1029933939417.0/13636558850479.0),
    ae51(14463281900351.0/6315353703477.0), ae53(66114435211212.0/5879490589093.0), ae54(-54053170152839.0/4284798021562.0),
    ae61(14090043504691.0/34967701212078.0), ae63(15191511035443.0/11219624916014.0), ae64(-18461159152457.0/12425892160975.0), ae65(-281667163811.0/9011619295870.0),
    ae71(19230459214898.0/13134317526959.0), ae73(21275331358303.0/2942455364971.0), ae74(-38145345988419.0/4862620318723.0), ae75(-1.0/8.0), ae76(-1.0/8.0),
    ae81(-19977161125411.0/11928030595625.0), ae83(-40795976796054.0/6384907823539.0), ae84(177454434618887.0/12078138498510.0), ae85(782672205425.0/8267701900261.0), ae86(-69563011059811.0/9646580694205.0), ae87(7356628210526.0/4942186776405.0),
        
    ai21(41.0/200.0), ai22(41.0/200.0),
    ai31(41.0/400.0), ai32(-567603406766.0/11931857230679.0), ai33(41.0/200.0),
    ai41(683785636431.0/9252920307686.0), ai43(-110385047103.0/1367015193373.0), ai44(41.0/200.0),
    ai51(3016520224154.0/10081342136671.0), ai53(30586259806659.0/12414158314087.0), ai54(-22760509404356.0/11113319521817.0), ai55(41.0/200.0),
    ai61(218866479029.0/1489978393911.0), ai63(638256894668.0/5436446318841.0), ai64(-1179710474555.0/5321154724896.0), ai65(-60928119172.0/8023461067671.0), ai66(41.0/200.0),
    ai71(1020004230633.0/5715676835656.0), ai73(25762820946817.0/25263940353407.0), ai74(-2161375909145.0/9755907335909.0), ai75(-211217309593.0/5846859502534.0), ai76(-4269925059573.0/7827059040749.0), ai77(41.0/200.0),
    ai81(-872700587467.0/9133579230613.0), ai84(22348218063261.0/9555858737531.0), ai85(-1143369518992.0/8141816002931.0), ai86(-39379526789629.0/19018526304540.0), ai87(32727382324388.0/42900044865799.0), ai88(41.0/200.0),

    b1(-872700587467.0/9133579230613.0), b4(22348218063261.0/9555858737531.0), b5(-1143369518992.0/8141816002931.0), b6(-39379526789629.0/19018526304540.0), b7(32727382324388.0/42900044865799.0), b8(41.0/200.0),

    b1_LO(-975461918565.0/9796059967033.0), b4_LO(78070527104295.0/32432590147079.0), b5_LO(-548382580838.0/3424219808633.0), b6_LO(-33438840321285.0/15594753105479.0), b7_LO(3629800801594.0/4656183773603.0), b8_LO(4035322873751.0/18575991585200.0)
    {
        
        Yhv1 = static_cast<complex<double> >(0.0);  Yhc1 = static_cast<complex<double> >(0.0);
        Yhv2 = static_cast<complex<double> >(0.0);  Yhc2 = static_cast<complex<double> >(0.0);
        Yhv3 = static_cast<complex<double> >(0.0);  Yhc3 = static_cast<complex<double> >(0.0);
        Yhv4 = static_cast<complex<double> >(0.0);  Yhc4 = static_cast<complex<double> >(0.0);
        Yhv5 = static_cast<complex<double> >(0.0);  Yhc5 = static_cast<complex<double> >(0.0);
        Yhv6 = static_cast<complex<double> >(0.0);  Yhc6 = static_cast<complex<double> >(0.0);
    }
//--------------------------------------------------------------
ARK54:: ~ARK54(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
}
void ARK54::take_step(State1D& Y4, State1D& Y5, double time, double h, VlasovFunctor1D_explicitE& vF, collisions_1D& coll, Parallel_Environment_1D& PE) 
{
//      Take a step using ARK43

//      Yh1, Stage 1
    // z1 = Y2;

    vF(Y5,Yhv1); PE.Neighbor_Communications(Yhv1);
    Yhv1 *= ae21*h; Yhc1 *= ai21*h;
    
    Yt = Y5;
    Yt += Yhv1; Yt += Yhc1;     

    coll(Yt,Yhc2,time,(ai22*h));   
    Yhc2 *= (ai22*h);    Yt += Yhc2;
    
    // z2 = Yt;

    vF(Yt,Yhv2); PE.Neighbor_Communications(Yhv2);
    Yhv1 *= ae31/ae21;  Yhc1 *= ai31/ai21;  
    Yhv2 *= ae32*h;     Yhc2 *= ai32/ai22;  

    Yt = Y5;
    Yt += Yhv1; Yt += Yhc1;
    Yt += Yhv2; Yt += Yhc2;
    coll(Yt,Yhc3,time,(ai33*h));

    Yhc3 *= (ai33*h);   Yt += Yhc3;
    
    // z3 = Yt;

    vF(Yt,Yhv3); PE.Neighbor_Communications(Yhv3);
    Yhv1 *= ae41/ae31;  Yhc1 *= ai41/ai31;  
    Yhv3 *= ae43*h;     Yhc3 *= ai43/ai33;

    Yt = Y5;
    Yt += Yhv1; Yt += Yhc1;
    Yt += Yhv3; Yt += Yhc3;

    coll(Yt,Yhc4,time,(ai44*h));

    Yhc4 *= (ai44*h);   Yt += Yhc4;
    
    // z4 = Yt;
    vF(Yt,Yhv4); PE.Neighbor_Communications(Yhv4);
    Yhv1 *= ae51/ae41;  Yhc1 *= ai51/ai41;  
    Yhv3 *= ae53/ae43;  Yhc3 *= ai53/ai43;  
    Yhv4 *= ae54*h;     Yhc4 *= ai54/ai44;

    Yt = Y5;
    Yt += Yhv1; Yt += Yhc1;
    Yt += Yhv3; Yt += Yhc3;
    Yt += Yhv4; Yt += Yhc4;

    coll(Yt,Yhc5,time,(ai55*h));

    Yhc5 *= (ai55*h);   Yt += Yhc5;
    
    // z5 = Yt;
    vF(Yt,Yhv5); PE.Neighbor_Communications(Yhv5);
    Yhv1 *= ae61/ae51;  Yhc1 *= ai61/ai51;  
    Yhv3 *= ae63/ae53;  Yhc3 *= ai63/ai53;  
    Yhv4 *= ae64/ae54;  Yhc4 *= ai64/ai54;  
    Yhv5 *= ae65*h;     Yhc5 *= ai65/ai55;

    Yt = Y5;
    Yt += Yhv1; Yt += Yhc1;
    Yt += Yhv3; Yt += Yhc3;
    Yt += Yhv4; Yt += Yhc4;
    Yt += Yhv5; Yt += Yhc5;

    coll(Yt,Yhc6,time,(ai66*h));

    Yhc6 *= (ai66*h);   Yt += Yhc6;

    // z6 = Yt;
    vF(Yt,Yhv6);  PE.Neighbor_Communications(Yhv6);
    Yhv1 *= ae71/ae61;  Yhc1 *= ai71/ai61;  
    Yhv3 *= ae73/ae63;  Yhc3 *= ai73/ai63;  
    Yhv4 *= ae74/ae64;  Yhc4 *= ai74/ai64;  
    Yhv5 *= ae75/ae65;  Yhc5 *= ai75/ai65;
    Yhv6 *= ae76*h;     Yhc6 *= ai76/ai66;

    Yt = Y5;
    Yt += Yhv1; Yt += Yhc1;
    Yt += Yhv3; Yt += Yhc3;
    Yt += Yhv4; Yt += Yhc4;
    Yt += Yhv5; Yt += Yhc5;
    Yt += Yhv6; Yt += Yhc6;

    coll(Yt,Yhc7,time,(ai77*h));

    Yhc7 *= (ai77*h);   Yt += Yhc7;

    // z6 = Yt;
    vF(Yt,Yhv7);  PE.Neighbor_Communications(Yhv7);
    Yhv1 *= ae81/ae71;  Yhc1 *= ai81/ai71;  
    Yhv3 *= ae83/ae73;  
    Yhv4 *= ae84/ae74;  Yhc4 *= ai84/ai74;  
    Yhv5 *= ae85/ae75;  Yhc5 *= ai85/ai75;
    Yhv6 *= ae86/ae76;  Yhc6 *= ai86/ai76;
    Yhv7 *= ae87*h;     Yhc7 *= ai87/ai77;

    Yt = Y5;
    Yt += Yhv1; Yt += Yhc1;
    Yt += Yhv3;
    Yt += Yhv4; Yt += Yhc4;
    Yt += Yhv5; Yt += Yhc5;
    Yt += Yhv6; Yt += Yhc6;
    Yt += Yhv7; Yt += Yhc7;

    coll(Yt,Yhc8,time,(ai88*h));

    Yhc8 *= (ai88*h);   Yt += Yhc8;

    // z6 = Yt;
    vF(Yt,Yhv8);  PE.Neighbor_Communications(Yhv8);

    //  Assemble 4th order solution
    Y4 = Y5;
    Yhv1 *= b1/ae81;    Y5 += Yhv1;
    Yhv4 *= b4/ae84;    Y5 += Yhv4;
    Yhv5 *= b5/ae85;    Y5 += Yhv5;
    Yhv6 *= b6/ae86;    Y5 += Yhv6;
    Yhv7 *= b7/ae87;    Y5 += Yhv7;
    Yhv8 *= b8*h;       Y5 += Yhv8;

    Yhc1 *= b1/ai81;    Y5 += Yhc1;
    Yhc4 *= b4/ai84;    Y5 += Yhc4;
    Yhc5 *= b5/ai85;    Y5 += Yhc5;
    Yhc6 *= b6/ai86;    Y5 += Yhc6;
    Yhc7 *= b7/ai87;    Y5 += Yhc7;
    Yhc8 *= b8/ai88;    Y5 += Yhc8;

    //  Assemble 3rd order solution
    Yhv1 *= b1_LO/b1;   Y4 += Yhv1;
    Yhv4 *= b4_LO/b4;   Y4 += Yhv4;
    Yhv5 *= b5_LO/b5;   Y4 += Yhv5;
    Yhv6 *= b6_LO/b6;   Y4 += Yhv6;
    Yhv7 *= b7_LO/b7;   Y4 += Yhv7;
    Yhv8 *= b8_LO/b8;   Y4 += Yhv8;

    Yhc1 *= b1_LO/b1;    Y4 += Yhc1;
    Yhc4 *= b4_LO/b4;    Y4 += Yhc4;
    Yhc5 *= b5_LO/b5;    Y4 += Yhc5;
    Yhc6 *= b6_LO/b6;    Y4 += Yhc6;
    Yhc7 *= b7_LO/b7;    Y4 += Yhc7;
    Yhc8 *= b8_LO/b8;    Y4 += Yhc8;

    Yhc8 *= 1./(b8_LO*h);
    Yhc1 = Yhc8;

//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
//--------------------------------------------------------------
RKCK45::RKCK45(State1D& Yin): Yh1(Yin), Yh3(Yin), Yh4(Yin), Yh5(Yin), Yh6(Yin), Yt(Yin),
        
        a21(0.2), 
        a31(3./40.), a32(9./40.),
        a41(.3), a42(-.9), a43(1.2),
        a51(-11./54.), a52(2.5) ,a53(-70./27.), a54(35./27.),
        a61(1631./55296.), a62(175./512.), a63(575./13824.), a64(44275./110592.), a65(253./4096.),
        b1_5(37./378.), b3_5(250./621.), b4_5(125./594.), b6_5(512./1771.),
        b1_4(2825./27648.), b3_4(18575./48384.), b4_4(13525./55296.), b5_4(277./14336.), b6_4(0.25)
    {}
//--------------------------------------------------------------
RKCK45:: ~RKCK45(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
}
void RKCK45::take_step(State1D& Y5, State1D& Y4, double time, double h, VlasovFunctor1D_explicitE& vF, collisions_1D& coll, Parallel_Environment_1D& PE) 
{
//  Take a step using RKCK

//      Yh1, Stage 1
    // z1 = Y2;
    vF(Y4,Yh1,time,1.); Yh1 *= h;
    Yh1 *= a21;
    PE.Neighbor_Communications(Yh1);
    Yt = Y4;    Yt  += Yh1;                              // Y1 = Y1 + (h/5)*Yh

    //      Step 2
    // vF(Yt,Yh2);                                   // f(Y1)
    vF(Yt,Y5,time,1.); Y5 *= h;                                   // f(Y1)
    Yh1 *= a31/a21; Y5 *= a32; 
    PE.Neighbor_Communications(Y5); 
    Yt = Y4;    Yt += Yh1;  Yt += Y5;

    //      Step 3
    vF(Yt,Yh3,time,1.); Yh3 *= h;
    Yh1 *= a41/a31; Y5 *= a42/a32;   Yh3 *= a43;
    PE.Neighbor_Communications(Yh3);
    Yt = Y4;    Yt += Yh1;  Yt += Y5; Yt += Yh3;
    
    //      Step 4
    vF(Yt,Yh4,time,1.); Yh4 *= h;
    Yh1 *= a51/a41; Y5 *= a52/a42;   Yh3 *= a53/a43;    Yh4 *= a54;
    PE.Neighbor_Communications(Yh4);
    Yt = Y4;    Yt += Yh1;  Yt += Y5; Yt += Yh3; Yt += Yh4;
    
    //      Step 5
    vF(Yt,Yh5,time,1.); Yh5 *= h;
    Yh1 *= a61/a51; Y5 *= a62/a52;   Yh3 *= a63/a53;    Yh4 *= a64/a54; Yh5 *= a65;
    PE.Neighbor_Communications(Yh5);
    Yt = Y4;    Yt += Yh1;  Yt += Y5; Yt += Yh3; Yt += Yh4; Yt += Yh5;
    
    //      Step 6
    vF(Yt,Yh6,time,1.); Yh6 *= h;
    PE.Neighbor_Communications(Yh6);


    //      Assemble 5th order solution
    Y5 = Y4;
    Yh1 *= b1_5/a61;    Y5 += Yh1;
    Yh3 *= b3_5/a63;    Y5 += Yh3;
    Yh4 *= b4_5/a64;    Y5 += Yh4;
    Yh6 *= b6_5;        Y5 += Yh6;

    //      Assemble 4th order solution
    Yh1 *= b1_4/b1_5;   Y4 += Yh1;
    Yh3 *= b3_4/b3_5;   Y4 += Yh3;
    Yh4 *= b4_4/b4_5;   Y4 += Yh4;
    Yh5 *= b5_4/a65;    Y4 += Yh5;
    Yh6 *= b6_4/b6_5;   Y4 += Yh6;
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
//--------------------------------------------------------------
RK4C::RK4C(State1D& Yin): Y0(Yin), Y1(Yin), Y2(Yin), Yh(Yin), //Y0_2D(), Y1_2D(), Y2_2D(), Yh_2D(),
                            onethird(1./3.), twothird(2./3.)
    {}

// RK4C::RK4C(State2D& Yin): Y0_2D(Yin), Y1_2D(Yin), Y2_2D(Yin), Yh_2D(Yin), Y0(), Y1(), Y2(), Yh(),
                            // onethird(1./3.), twothird(2./3.)
    // {}    
//--------------------------------------------------------------
RK4C:: ~RK4C(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
}
void RK4C::take_step(State1D& Ystar, State1D& Y, double time, double h, VlasovFunctor1D_explicitE& vF, collisions_1D& coll, Parallel_Environment_1D& PE) 
{
//  Take a step using RKCK

    // Initialization
        Y0 = Y; Y1 = Y;

//      Step 1
        vF(Y1,Yh,time,h);                    // slope in the beginning
        PE.Neighbor_Communications(Yh);
        Yh *= (0.5*h);   Y1 += Yh;      // Y1 = Y1 + (h/2)*Yh  Yhc = (*CF)(Y1,time,0.5*h)
        Yh *= (onethird); Y  += Yh;      // Y  = Y  + (h/6)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        vF(Y1,Yh,time,h);     Y1  = Y0;      // slope in the middle
        PE.Neighbor_Communications(Yh);
        Yh *= (0.5*h);   Y1 += Yh;      // Y1 = Y0 + (h/2)*Yh
        Yh *= (twothird); Y  += Yh;      // Y  = Y  + (h/3)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 3
        vF(Y1,Yh,time,h);                    // slope in the middle again
        PE.Neighbor_Communications(Yh);
        Yh *= h;          Y0 += Yh;     // Y0 = Y0 + h*Yh
        Yh *= (onethird);  Y  += Yh;     // Y  = Y  + (h/3)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 4
        vF(Y0,Yh,time,h);                    // slope at the end
        PE.Neighbor_Communications(Yh);
        Yh *= (h/6.0);    Y += Yh;      // Y  = Y  + (h/6)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
//--------------------------------------------------------------
// void RK4C::take_step(State2D& Ystar, State2D& Y, double time, double h, VlasovFunctor2D_explicitE& vF, collisions_2D& coll, Parallel_Environment_2D& PE) 
// {
// //  Take a step using RKCK

//     // Initialization
//         Y0_2D = Y; Y1_2D = Y;

// //      Step 1
//         vF(Y1_2D,Yh_2D,time,h);                    // slope in the beginning
//         PE.Neighbor_Communications(Yh_2D);
//         Yh_2D *= (0.5*h);   Y1_2D += Yh_2D;      // Y1 = Y1 + (h/2)*Yh  Yhc = (*CF)(Y1,time,0.5*h)
//         Yh_2D *= (onethird); Y  += Yh_2D;      // Y  = Y  + (h/6)*Yh
// //      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// //      Step 2
//         vF(Y1_2D,Yh_2D,time,h);     Y1_2D  = Y0_2D;      // slope in the middle
//         PE.Neighbor_Communications(Yh_2D);
//         Yh_2D *= (0.5*h);   Y1_2D += Yh_2D;      // Y1 = Y0 + (h/2)*Yh
//         Yh_2D *= (twothird); Y  += Yh_2D;      // Y  = Y  + (h/3)*Yh
// //      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// //      Step 3
//         vF(Y1_2D,Yh_2D,time,h);                    // slope in the middle again
//         PE.Neighbor_Communications(Yh_2D);
//         Yh_2D *= h;          Y0_2D += Yh_2D;     // Y0 = Y0 + h*Yh
//         Yh_2D *= (onethird);  Y  += Yh_2D;     // Y  = Y  + (h/3)*Yh
// //      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// //      Step 4
//         vF(Y0_2D,Yh_2D,time,h);                    // slope at the end
//         PE.Neighbor_Communications(Yh_2D);
//         Yh_2D *= (h/6.0);    Y += Yh_2D;      // Y  = Y  + (h/6)*Yh
// //      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// //      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// }
//--------------------------------------------------------------
/*RKT54::RKT54(State1D& Yin): Yh1(Yin), Yh2(Yin), Yh3(Yin), Yh4(Yin), Yh5(Yin), Yh6(Yin), Yh7(Yin), Yt(Yin),
        
        a21(0.161),
        a31(-0.008480655492356989),   a32(0.335480655492357),
        a41(2.8971530571054935),    a42(-6.359448489975075),    a43(4.3622954328695815),
        a51(5.325864828439257),   a52(-11.748883564062828),  a53(7.4955393428898365),  a54(-0.09249506636175525),
        a61(5.86145544294642),    a62(-12.92096931784711),    a63(8.159367898576159), a64(-0.071584973281401), a65(-0.028269050394068383),
        a71(0.09646076681806523), a72(0.01),   a73(0.4798896504144996),   a74(1.379008574103742),    a75(-3.290069515436081),    a76(2.324710524099774),

        b1_5(0.09468075576583945),    b2_5(0.009183565540343254),   b3_5(0.4877705284247616),    b4_5(1.234297566930479),  b5_5(-2.7077123499835256),  b6_5(1.866628418170587),    b7_5(0.015151515151515152),
        btilde1(-0.00178001105222577714),   btilde2(-0.0008164344596567469),   btilde3(0.007880878010261995), btilde4(-0.1447110071732629),    btilde5(0.5823571654525552),    btilde6(-0.45808210592918697),  btilde7(0.015151515151515152)
    {}
//--------------------------------------------------------------
RKT54:: ~RKT54(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
}
void RKT54::take_step(State1D& Y4, State1D& Y5, double time, double h, VlasovFunctor1D_explicitE& vF, collisions_1D& coll, Parallel_Environment_1D& PE) 
{
//      Step 1
    vF(Y5,Yh1); Yh1 *= h;
    Yh1 *= a21;
    PE.Neighbor_Communications(Yh1);
    Yt = Y5;    Yt  += Yh1;                              // Y1 = Y1 + (h/5)*Yh

    //      Step 2
    vF(Yt,Yh2); Yh2 *= h;                                   // f(Y1)
    Yh1 *= a31/a21; Yh2 *= a32;
    PE.Neighbor_Communications(Yh2);  
    Yt = Y5;    Yt += Yh1;  Yt += Yh2;

    //      Step 3
    vF(Yt,Yh3); Yh3 *= h;
    Yh1 *= a41/a31; Yh2 *= a42/a32;   Yh3 *= a43;
    PE.Neighbor_Communications(Yh3);
    Yt = Y5;    Yt += Yh1; Yt += Yh2; Yt += Yh3;
    
    //      Step 4
    vF(Yt,Yh4); Yh4 *= h;
    Yh1 *= a51/a41; Yh2 *= a52/a42;   Yh3 *= a53/a43;    Yh4 *= a54;
    PE.Neighbor_Communications(Yh4);
    Yt = Y5;    Yt += Yh1; Yt += Yh2; Yt += Yh3; Yt += Yh4;
    
    //      Step 5
    vF(Yt,Yh5); Yh5 *= h;
    Yh1 *= a61/a51; Yh2 *= a62/a52;   Yh3 *= a63/a53;    Yh4 *= a64/a54; Yh5 *= a65;
    PE.Neighbor_Communications(Yh5);
    Yt = Y5;    Yt += Yh1; Yt += Yh2; Yt += Yh3; Yt += Yh4; Yt += Yh5;
        
    //      Step 6
    vF(Yt,Yh6); Yh6 *= h;
    Yh1 *= a71/a61; Yh2 *= a72/a62;   Yh3 *= a73/a63;    Yh4 *= a74/a64; Yh5 *= a75/a65; Yh6 *= a76;
    PE.Neighbor_Communications(Yh6);
    Yt = Y5;    Yt += Yh1;  Yt += Yh2;   Yt += Yh3;  Yt += Yh4;  Yt += Yh5;  Yt += Yh6;

    //      Step 7
    vF(Yt,Yh7); Yh7 *= h;
    PE.Neighbor_Communications(Yh7);

    //      Assemble 5th order solution
    Y4 = Y5;
    Yh1 *= b1_5/a71;    Y5 += Yh1;
    Yh2 *= b2_5/a72;    Y5 += Yh2;
    Yh3 *= b3_5/a73;    Y5 += Yh3;
    Yh4 *= b4_5/a74;    Y5 += Yh4;
    Yh5 *= b5_5/a75;    Y5 += Yh5;
    Yh6 *= b6_5/a76;    Y5 += Yh6;
    Yh7 *= b7_5;        Y5 += Yh7;

    //      Assemble 4th order solution
    Yh1 *= btilde1/b1_5;  Y4 += Yh1;
    Yh2 *= btilde2/b2_5;  Y4 += Yh2;
    Yh3 *= btilde3/b3_5;  Y4 += Yh3;
    Yh4 *= btilde4/b4_5;  Y4 += Yh4;
    Yh5 *= btilde5/b5_5;  Y4 += Yh5;
    Yh6 *= btilde6/b6_5;  Y4 += Yh6;
    Yh7 *= btilde7/b7_5;  Y4 += Yh7;

//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}*/
//--------------------------------------------------------------
RKDP85::RKDP85(State1D& Yin): Yh1(Yin), Yh2(Yin), Yh3(Yin), Yh4(Yin), Yh5(Yin), Yh6(Yin),
                                Yh7(Yin), Yh8(Yin), Yh9(Yin), Yh10(Yin), Yt(Yin),

  a0201(0.05260015195876773),
  a0301(0.0197250569845379), a0302(0.0591751709536137),
  a0401(0.02958758547680685), a0403(0.08876275643042054),
  a0501(0.2413651341592667), a0503(-0.8845494793282861), a0504(0.924834003261792),
  a0601(0.037037037037037035), a0604(0.17082860872947386), a0605(0.12546768756682242),
  a0701(0.037109375), a0704(0.17025221101954405), a0705(0.06021653898045596), a0706(-0.017578125),
  a0801(0.03709200011850479), a0804(0.17038392571223998), a0805(0.10726203044637328), a0806(-0.015319437748624402), a0807(0.008273789163814023),
  a0901(0.6241109587160757), a0904(-3.3608926294469414), a0905(-0.868219346841726), a0906(27.59209969944671), a0907(20.154067550477894), a0908(-43.48988418106996),
  a1001(0.47766253643826434), a1004(-2.4881146199716677), a1005(-0.590290826836843), a1006(21.230051448181193), a1007(15.279233632882423), a1008(-33.28821096898486), a1009(-0.020331201708508627),
  a1101(-0.9371424300859873), a1104(5.186372428844064), a1105(1.0914373489967295), a1106(-8.149787010746927), a1107(-18.52006565999696), a1108(22.739487099350505), a1109(2.4936055526796523), a1110(-3.0467644718982196),
  a1201(2.273310147516538), a1204(-10.53449546673725), a1205(-2.0008720582248625), a1206(-17.9589318631188), a1207(27.94888452941996), a1208(-2.8589982771350235), a1209(-8.87285693353063), a1210(12.360567175794303), a1211(0.6433927460157636),

  b1(0.054293734116568765), b6(4.450312892752409), b7(1.8915178993145003), b8(-5.801203960010585), b9(0.3111643669578199), b10(-0.1521609496625161), b11(0.20136540080403034), b12(0.04471061572777259),
  bhh1(0.2440944881889764), bhh2(0.7338466882816118), bhh3(0.022058823529411766),
  er1(0.01312004499419488), er6(-1.2251564463762044), er7(-0.4957589496572502), er8(1.6643771824549864), er9(-0.35032884874997366), er10(0.3341791187130175), er11(0.08192320648511571), er12(-0.022355307863886294)
{}
//--------------------------------------------------------------
RKDP85:: ~RKDP85(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
}
void RKDP85::take_step(State1D& Y5, State1D& Y8, double time, double h, VlasovFunctor1D_explicitE& vF, collisions_1D& coll, Parallel_Environment_1D& PE) 
{
//      Step 1
    vF(Y8,Yh1,time,1.); Yh1 *= h;
    Yh1 *= a0201;
    PE.Neighbor_Communications(Yh1);
    Yt = Y8;    Yt  += Yh1;                              // Y1 = Y1 + (h/5)*Yh

    //      Step 2
    vF(Yt,Yh2,time,1.); Yh2 *= h;                                   // f(Y1)
    Yh1 *= a0301/a0201; Yh2 *= a0302;
    PE.Neighbor_Communications(Yh2);  
    Yt = Y8;    Yt += Yh1;  Yt += Yh2;

    //      Step 3
    vF(Yt,Yh3,time,1.); Yh3 *= h;
    Yh1 *= a0401/a0301; Yh3 *= a0403;
    PE.Neighbor_Communications(Yh3);
    Yt = Y8;    Yt += Yh1; Yt += Yh3;
    
    //      Step 4
    vF(Yt,Yh4,time,1.); Yh4 *= h;
    Yh1 *= a0501/a0401;    Yh3 *= a0503/a0403;    Yh4 *= a0504;
    PE.Neighbor_Communications(Yh4);
    Yt = Y8;    Yt += Yh1;  Yt += Yh3; Yt += Yh4;
    
    //      Step 5
    vF(Yt,Yh5,time,1.); Yh5 *= h;
    Yh1 *= a0601/a0501; Yh4 *= a0604/a0504; Yh5 *= a0605;
    PE.Neighbor_Communications(Yh5);
    Yt = Y8;    Yt += Yh1; Yt += Yh4; Yt += Yh5;
        
    //      Step 6
    vF(Yt,Yh6,time,1.); Yh6 *= h;
    Yh1 *= a0701/a0601; Yh4 *= a0704/a0604; Yh5 *= a0705/a0605; Yh6 *= a0706;
    PE.Neighbor_Communications(Yh6);
    Yt = Y8;    Yt += Yh1;  Yt += Yh4;  Yt += Yh5;  Yt += Yh6;

    //      Step 7
    vF(Yt,Yh7,time,1.); Yh7 *= h;
    Yh1 *= a0801/a0701; Yh4 *= a0804/a0704; Yh5 *= a0805/a0705; Yh6 *= a0806/a0706; Yh7 *= a0807;
    PE.Neighbor_Communications(Yh7);
    Yt = Y8;    Yt += Yh1;  Yt += Yh4;  Yt += Yh5;  Yt += Yh6;  Yt += Yh7;

    //      Step 8
    vF(Yt,Yh8,time,1.); Yh8 *= h;
    Yh1 *= a0901/a0801; Yh4 *= a0904/a0804; Yh5 *= a0905/a0805; Yh6 *= a0906/a0806; Yh7 *= a0907/a0807; Yh8 *= a0908;
    PE.Neighbor_Communications(Yh8);
    Yt = Y8;    Yt += Yh1;  Yt += Yh4;  Yt += Yh5;  Yt += Yh6;  Yt += Yh7;  Yt += Yh8;

    //      Step 9
    vF(Yt,Yh9,time,1.); Yh9 *= h;
    Yh1 *= a1001/a0901; Yh4 *= a1004/a0904; Yh5 *= a1005/a0905; Yh6 *= a1006/a0906; Yh7 *= a1007/a0907; Yh8 *= a1008/a0908; Yh9 *= a1009;
    PE.Neighbor_Communications(Yh9);
    Yt = Y8;    Yt += Yh1;  Yt += Yh4;  Yt += Yh5;  Yt += Yh6;  Yt += Yh7;  Yt += Yh8;  Yt += Yh9;  

    //      Step 10
    vF(Yt,Yh10,time,1.); Yh10 *= h;
    Yh1 *= a1101/a1001; Yh4 *= a1104/a1004; Yh5 *= a1105/a1005; Yh6 *= a1106/a1006; Yh7 *= a1107/a1007; Yh8 *= a1108/a1008; Yh9 *= a1109/a1009; Yh10 *= a1110;
    PE.Neighbor_Communications(Yh10);
    Yt = Y8;    Yt += Yh1;  Yt += Yh4;  Yt += Yh5;  Yt += Yh6;  Yt += Yh7;  Yt += Yh8;  Yt += Yh9;  Yt += Yh10;  

    //      Step 12
    vF(Yt,Yh2,time,1.); Yh2 *= h;
    Yh1 *= a1201/a1101; Yh4 *= a1204/a1104; Yh5 *= a1205/a1105; Yh6 *= a1206/a1106; Yh7 *= a1207/a1107; Yh8 *= a1208/a1108; Yh9 *= a1209/a1109; Yh10 *= a1210/a1110; Yh2 *= a1211;
    PE.Neighbor_Communications(Yh2);
    Yt = Y8;    Yt += Yh1;  Yt += Yh4;  Yt += Yh5;  Yt += Yh6;  Yt += Yh7;  Yt += Yh8;  Yt += Yh9;  Yt += Yh10; Yt += Yh2;

    //      Step 13
    vF(Yt,Yh3,time,1.); Yh3 *= h;
    PE.Neighbor_Communications(Yh3);
        
    //      Assemble 5th order solution
    // Y5 = Y8;
    
    Yh1 *= b1/a1201;    Y8 += Yh1;
    Yh6 *= b6/a1206;    Y8 += Yh6;
    Yh7 *= b7/a1207;    Y8 += Yh7;
    Yh8 *= b8/a1208;    Y8 += Yh8;
    Yh9 *= b9/a1209;    Y8 += Yh9;
    Yh10*= b10/a1210;   Y8 += Yh10;
    Yh2 *= b11/a1211;   Y8 += Yh2;
    Yh3 *= b12;         Y8 += Yh3;

    Y5 = Y8;

    //      Assemble 4th order solution
    // Yh1 *= er1/b1;      Y5 += Yh1;
    // Yh6 *= er6/b6;      Y5 += Yh6;
    // Yh7 *= er7/b7;      Y5 += Yh7;
    // Yh8 *= er8/b8;      Y5 += Yh8;
    // Yh9 *= er9/b9;      Y5 += Yh9;
    // Yh10*= er10/b10;    Y5 += Yh10;
    // Yh2 *= er11/b11;    Y5 += Yh2;
    // Yh3 *= er12/b12;    Y5 += Yh3;
    // 
    Yh1 *= (-1.*bhh1)/b1;   Y5 += Yh1;
    Yh9 *= (-1.*bhh2)/b9;   Y5 += Yh9;
    Yh3 *= (-1.*bhh3)/b12;  Y5 += Yh3;


//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
//--------------------------------------------------------------
