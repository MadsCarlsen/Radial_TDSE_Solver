//
// Created by au642261 on 8/3/25.
//

#ifndef RADIALPROPAGATOR_H
#define RADIALPROPAGATOR_H

#include "Types.h"
#include "Fields.h"

enum class Gauge {
    velocity,
    length
};

struct settingsStruct {
    // Grid settings
    int Nr;
    double dr;
    int lMax;

    // Imag time propagation

    // Real time propagation settings
    Gauge gauge;
    double dt;
    double tFin;  // In times of total field duration

    // Field settings
    int NCycles;
    double lambda_nm;
    double CEP;
    double I0_Wcm2;
};

class RadialPropagator {
public:
    const double sqrt2Inv = 1./std::sqrt(2.);
    settingsStruct settings;
    Fields laserFields;
    int Nr;
    double dr, dr_inv2;
    int lMax;
    dVec rVec;

    Gauge gauge;

    // CN matrices
    dMat D2, M2;  // Second derivative and Numerov boost matrix
    dMat D1, M1;  // First derivative and Numerov boost matrix
    std::vector<cdMat> AtomCNForward;  // All the atomic forward CN matrices (for each l)
    std::vector<cdMat> AtomCNBackward;  // All the atomic backwards CN matrices
    double deltaCoulomb;  // Fixes for l=m=0 to match Coulomb solution

    dMat potMat;  // Matrix with potential vector for each angular momentum
    dVec clmVec;  // Vector with Clebsch-Gordan coeffs

    // Various arrays used in calculations
    cdVec phiTemp_l;  // Temporary wavefunction arrays
    cdVec phiTemp_l1;
    dVec FrVec;  // Terms in the R^lm matrix
    dVec slmVec;  // Terms in the R^lm matrix
    dVec slmFac1Vec;  // Terms in the R^lm matrix
    dVec slmFac2Vec;  // Terms in the R^lm matrix
    dMat Glm_p;  // Matrix for performing Mixing step
    dMat Glm_m;


    RadialPropagator(settingsStruct settings);
    double potential(double r, int l);  // Radial SAE potential
    void prepareDerivMats();
    void prepareAtomicCNMats(double dt);

    // Imaginary time propagation and eigenstate determination
    void groundStateImagProp(dVec& phi0, double dt, int NSteps);

    // Real time propagation
    void realTimeProp(cdMat& psi, double dt, double tFin);  // The time loop itself, calls propagators
    void velocityPropagator(cdMat& psi, double Ai, double dt);
    void lengthPropagator(cdMat& psi, double Ei, double dt);





};


#endif //RADIALPROPAGATOR_H
