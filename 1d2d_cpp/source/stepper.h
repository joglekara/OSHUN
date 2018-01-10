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

