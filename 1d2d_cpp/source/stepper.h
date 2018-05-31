/*!\brief  Functors for various time-integration methodds - Declarations
 * \author PICKSC
 * \file   stepper.h
 *
 * Includes functors for fully explicit, implicit B, implicit E
 *
 */
#ifndef OSHUN_STEPPERS_H
#define OSHUN_STEPPERS_H

//**************************************************************

class ARK32 {
public:
//      Constructor
    ARK32(State1D& Yin);
    ~ARK32();

    void take_step(State1D& Y3, State1D& Y2, double time, double h, 
    	VlasovFunctor1D_explicitE& vF, collisions_1D& cF, Parallel_Environment_1D& PE);
private:

    State1D  Yhv1, Yhv2, Yhv3, Yhv4, Yhc1, Yhc2, Yhc3, Yhc4, Yt;

    double ae21;
    double ae31,ae32;
    double ae41,ae42,ae43;
    double ai21, ai22;
    double ai31,ai32,ai33;
    double ai41,ai42,ai43,ai44;

    double b1,b2,b3,b4;
    double b1_LO,b2_LO,b3_LO,b4_LO;
};
//--------------------------------------------------------------
class ARK43 {
public:
//      Constructor
    ARK43(State1D& Yin);
    ~ARK43();

    void take_step(State1D& Y3, State1D& Y4, double time, double h, 
        VlasovFunctor1D_explicitE& vF, collisions_1D& cF, Parallel_Environment_1D& PE);
private:

    State1D  Yhv1, Yhv2, Yhv3, Yhv4, Yhv5, Yhv6;
    State1D  Yhc1, Yhc2, Yhc3, Yhc4, Yhc5, Yhc6, Yt;

    double ae21;
    double ae31,ae32;
    double ae41,ae42,ae43;
    double ae51,ae52,ae53,ae54;
    double ae61,ae62,ae63,ae64,ae65;
    
    double ai21, ai22;
    double ai31,ai32,ai33;
    double ai41,ai42,ai43,ai44;
    double ai51,ai52,ai53,ai54,ai55;
    double ai61,ai62,ai63,ai64,ai65,ai66;

    double b1,b2,b3,b4,b5,b6;
    double b1_LO,b3_LO,b4_LO,b5_LO,b6_LO;
};
//--------------------------------------------------------------
//--------------------------------------------------------------
class ARK54 {
public:
//      Constructor
    ARK54(State1D& Yin);
    ~ARK54();

    void take_step(State1D& Y4, State1D& Y5, double time, double h, 
        VlasovFunctor1D_explicitE& vF, collisions_1D& cF, Parallel_Environment_1D& PE);
private:

    State1D  Yhv1, Yhv2, Yhv3, Yhv4, Yhv5, Yhv6, Yhv7, Yhv8;
    State1D  Yhc1, Yhc2, Yhc3, Yhc4, Yhc5, Yhc6, Yhc7, Yhc8, Yt;

    double ae21;
    double ae31,ae32;
    double ae41,ae42,ae43;
    double ae51,ae52,ae53,ae54;
    double ae61,ae62,ae63,ae64,ae65;
    double ae71,ae72,ae73,ae74,ae75, ae76;
    double ae81,ae78,ae83,ae84,ae85, ae86, ae87;

    
    double ai21, ai22;
    double ai31,ai32,ai33;
    double ai41,ai42,ai43,ai44;
    double ai51,ai52,ai53,ai54,ai55;
    double ai61,ai62,ai63,ai64,ai65,ai66;
    double ai71,ai72,ai73,ai74,ai75,ai76,ai77;
    double ai81,ai82,ai83,ai84,ai85,ai86,ai87,ai88;

    double b1,b4,b5,b6,b7,b8;
    double b1_LO,b4_LO,b5_LO,b6_LO,b7_LO,b8_LO;
};
//--------------------------------------------------------------
class RKCK45 {
public:
//      Constructor
    RKCK45(State1D& Yin);
    ~RKCK45();

    void take_step(State1D& Y5, State1D& Y4, double time, double h, 
        VlasovFunctor1D_explicitE& vF, collisions_1D& cF, Parallel_Environment_1D& PE);
private:

    State1D  Yh1, Yh3, Yh4, Yh5, Yh6, Yt;

    double a21;
    double a31,a32;
    double a41,a42,a43;
    double a51,a52,a53,a54;
    double a61,a62,a63,a64,a65;

    double b1_5,b3_5,b4_5,b6_5;
    double b1_4,b3_4,b4_4,b5_4,b6_4;

};
//--------------------------------------------------------------
class RK4C {
public:
//      Constructor
    RK4C(State1D& Yin);
    // RK4C(State2D& Yin);
    ~RK4C();

    void take_step(State1D& Y5, State1D& Y4, double time, double h, 
        VlasovFunctor1D_explicitE& vF, collisions_1D& cF, Parallel_Environment_1D& PE);

    // void take_step(State2D& Y5, State2D& Y4, double time, double h, 
    //     VlasovFunctor2D_explicitE& vF, collisions_2D& cF, Parallel_Environment_2D& PE);
private:

    State1D  Y0, Y1, Y2, Yh;

    // State2D  Y0_2D, Y1_2D, Y2_2D, Yh_2D;

    complex<double> onethird, twothird;
};
//--------------------------------------------------------------
class RKDP85 {
public:
//      Constructor
    RKDP85(State1D& Yin);
    ~RKDP85();

    void take_step(State1D& Y5, State1D& Y4, double time, double h, 
        VlasovFunctor1D_explicitE& vF, collisions_1D& cF, Parallel_Environment_1D& PE);
private:

    State1D Yh1, Yh2, Yh3, Yh4, Yh5, Yh6;
    State1D Yh7, Yh8, Yh9, Yh10, Yt;

    double a0201;
    double a0301, a0302;
    double a0401, a0403;
    double a0501, a0503, a0504;
    double a0601, a0604, a0605;
    double a0701, a0704, a0705, a0706;
    double a0801, a0804, a0805, a0806, a0807;
    double a0901, a0904, a0905, a0906, a0907, a0908;
    double a1001, a1004, a1005, a1006, a1007, a1008, a1009;
    double a1101, a1104, a1105, a1106, a1107, a1108, a1109, a1110;
    double a1201, a1204, a1205, a1206, a1207, a1208, a1209, a1210, a1211;
    
    double b1, b6, b7, b8, b9, b10, b11, b12;
    double bhh1, bhh2, bhh3;
    double er1, er6, er7, er8, er9, er10, er11, er12;

};
//--------------------------------------------------------------
/*class RKT54 {
public:
//      Constructor
    RKT54(State1D& Yin);
    ~RKT54();

    void take_step(State1D& Y4, State1D& Y5, double time, double h, 
        VlasovFunctor1D_explicitE& vF, collisions_1D& cF, Parallel_Environment_1D& PE);
private:

    State1D  Yh1, Yh2, Yh3, Yh4, Yh5, Yh6, Yh7, Yt;

    double a21;
    double a31,a32;
    double a41,a42,a43;
    double a51,a52,a53,a54;
    double a61,a62,a63,a64,a65;
    double a71,a72,a73,a74,a75,a76;
    

    double b1_5,b2_5,b3_5,b4_5,b5_5,b6_5, b7_5;
    double btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7;

};*/
//--------------------------------------------------------------
//**************************************************************
#endif //OSHUN_STEPPERS_H

