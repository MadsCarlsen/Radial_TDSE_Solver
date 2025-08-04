//
// Created by au642261 on 8/3/25.
//

#include "RadialPropagator.h"
#include <numeric>
#include <iostream>
#include "LinalgRoutines.h"

RadialPropagator::RadialPropagator(settingsStruct settings) : settings(settings) {
    laserFields = Fields(settings.lambda_nm, settings.I0_Wcm2, settings.CEP, settings.NCycles);

    // Prepare the radial grid
    Nr = settings.Nr;
    dr = settings.dr;
    dr_inv2 = 1./(dr*dr);
    rVec = std::vector<double>(Nr);
    for (int i=0; i<Nr; i++) {
        rVec[i] = (i+1.)*dr;
    }
    lMax = settings.lMax;

    // Prepare various matrices representing the derivatives on the grid
    prepareDerivMats();

    // Build matrix with potential vectors
    potMat = dMat(lMax+1, dVec(Nr));
    for (int li=0; li<lMax+1; li++) {
        for (int i=0; i<Nr; i++) {
            potMat[li][i] = potential(rVec[i], li);
        }
    }

    // Build vector with Clebsch-Gordan coeffs
    clmVec = dVec(lMax+1);
    for (int i=0; i<lMax+1; i++) {
        auto li = static_cast<double>(i);
        clmVec[i] = std::sqrt((li+1.)*(li+1.) / ((2.*li+1.)*(2.*li+3.)));
    }

    // Various arrays used in propagation - just here to prevent new initialization for each step...
    phiTemp_l = cdVec(Nr);  // Temporary wavefunction arrays
    phiTemp_l1 = cdVec(Nr);
    FrVec = dVec(Nr);  // Terms in the R^lm matrix
    slmVec = dVec(Nr);  // Terms in the R^lm matrix
    slmFac1Vec = dVec(Nr);  // Terms in the R^lm matrix
    slmFac2Vec = dVec(Nr);  // Terms in the R^lm matrix
    Glm_p = dMat(3, dVec(Nr));  // Matrix for performing Mixing step
    Glm_m = dMat(3, dVec(Nr));
}

double RadialPropagator::potential(double r, int l) {
    return -1./r + l*(l+1.)/(2.*r*r);
}

void RadialPropagator::prepareDerivMats() {
    // DOUBLE DERIV MAT - NB [0] is the lower diag, [1] is diag, [2] is upper diag!
    D2 = dMat(3, std::vector<double>(Nr));
    D2[0] = std::vector<double>(Nr, dr_inv2);
    D2[1] = std::vector<double>(Nr, -2. * dr_inv2);
    D2[2] = std::vector<double>(Nr, dr_inv2);

    // Numerov boost mat
    M2 = dMat(3, std::vector<double>(Nr));
    M2[0] = std::vector<double>(Nr, -1./6.);
    M2[1] = std::vector<double>(Nr, -10./6.);
    M2[2] = std::vector<double>(Nr, -1./6.);

    // Fixes for m=l=0
    deltaCoulomb = -2./(dr*dr) * (1. - dr/(12.-10.*dr));
    D2[1][0] = deltaCoulomb;
    M2[1][0] = -2.*(1.+dr*dr/12. * deltaCoulomb);


    // SINGLE DERIV MAT
    D1 = dMat(3, std::vector<double>(Nr));
    D1[0] = std::vector<double>(Nr, -1./(2.*dr));
    D1[1] = std::vector<double>(Nr, 0.);
    D1[2] = std::vector<double>(Nr, 1./(2.*dr));
    double gamma = std::sqrt(3.) - 2.;
    D1[1][0] = gamma; D1[1][Nr-1] = -gamma;  // Corrections to make Hermitian

    // Single deriv numerov boost mat
    M1 = dMat(3, std::vector<double>(Nr));
    M1[0] = std::vector<double>(Nr, 1./6.);
    M1[1] = std::vector<double>(Nr, 4.);
    M1[2] = std::vector<double>(Nr, 1./6.);
    M1[0][0] = 4.+gamma; M1[0][Nr-1] = 4.+gamma;  // Corrections to make Hermitian
}

void RadialPropagator::prepareAtomicCNMats(double dt) {
    AtomCNForward = std::vector<cdMat> (lMax+1, cdMat(3, cdVec(Nr)));
    AtomCNBackward = std::vector<cdMat>(lMax+1, cdMat(3, cdVec(Nr)));

    for (int li=0; li<lMax+1; li++) {
        cdMat& CNForward = AtomCNForward[li];
        cdMat& CNBackward = AtomCNBackward[li];
        dMat MVprod = dTridiagDiagMult(M2, potMat[li]);
        for (int i=0; i<3; i++) {
            for (int j=0; j<Nr; j++) {
                CNForward[i][j] = M2[i][j] - Cdouble(0.,1.) * dt/2.*(D2[i][j] + MVprod[i][j]);
                CNBackward[i][j] = M2[i][j] + Cdouble(0.,1.) * dt/2.*(D2[i][j] + MVprod[i][j]);
            }
        }
        if (li == 0) {
            // Remove the Coulomb fix as only for l=m=0
            D2[1][0] = -2. * dr_inv2;
            M2[1][0] = -10./6.;
        }
    }
}

void RadialPropagator::velocityPropagator(cdMat& psi, double Ai, double dt) {
    // Prebuild some of the R^lm matrix
    for (int i=0; i<Nr; i++) {
        FrVec[i] = Ai/rVec[i] * dt/4.;
    }

    for (int li=0; li<lMax; li++) {  // NB last li term has no mat mult - included only in previous
        cdVec& phi_l = psi[li];
        cdVec& phi_l1 = psi[li+1];

        // PURE ANGULAR TERM  TODO CHECK THAT LOOPS HERE ARE VECTORIZED!
        // Build the R^lm matrix for each r
        for (int i=0; i<Nr; i++) {
            slmVec[i] = FrVec[i] * (li+1.)*clmVec[li];
        }
        for (int i=0; i<Nr; i++) {
            slmFac1Vec[i] = (1.-slmVec[i]*slmVec[i]) / (1.+slmVec[i]*slmVec[i]);
            slmFac2Vec[i] = 2.*slmVec[i] / (1. + slmVec[i]*slmVec[i]);
        }

        // Apply the 2x2 R^lm matrix for each r
        for (int i=0; i<Nr; i++) {  // TODO Probably need temporary vecs here if I want vecoriation!
            Cdouble temp_l = slmFac1Vec[i]*phi_l[i] - slmFac2Vec[i]*phi_l1[i];
            Cdouble temp_l1 = slmFac2Vec[i]*phi_l[i] + slmFac1Vec[i]*phi_l1[i];
            phi_l[i] = temp_l;
            phi_l1[i] = temp_l1;
        }

        // MIXING TERM
        // Transform using B (COULD MAYBE JOIN THIS LOOP WITH THE ABOVE TO ENABLE SMOOOTH VECTORIZATION?)
        for (int i=0; i<Nr; i++) {
            phiTemp_l[i] = sqrt2Inv * (phi_l[i] + phi_l1[i]);
            phiTemp_l1[i] = sqrt2Inv * (-phi_l[i] + phi_l1[i]);
        }

        // Apply G^lm operator on two new phi
        double tempFac = dt/4.*Ai*clmVec[li];
        for (int i=0; i<3; i++) {
            for (int j=0; j<Nr; j++) {
                Glm_p[i][j] = M1[i][j] + tempFac * D1[i][j];
                Glm_m[i][j] = M1[i][j] - tempFac * D1[i][j];
            }
        }
        tridiagMatVecMult(phi_l, Glm_m, phiTemp_l);
        tridiagMatVecMult(phi_l1, Glm_p, phiTemp_l1);
        tridiagSolver(phiTemp_l, Glm_p, phi_l);
        tridiagSolver(phiTemp_l1, Glm_m, phi_l1);

        // Transform back using B^T - also brings data back to phi_l and phi_l1
        for (int i=0; i<Nr; i++) {
            phi_l[i] = sqrt2Inv * (phiTemp_l[i] - phiTemp_l1[i]);
            phi_l1[i] = sqrt2Inv * (phiTemp_l[i] + phiTemp_l1[i]);
        }
    }

    // ATOMIC HAMILTONIAN
    for (int li=0; li<lMax+1; li++) {
        cdVec& phi_l = psi[li];
        cdMat& CNForward = AtomCNForward[li];
        cdMat& CNBackward = AtomCNBackward[li];
        tridiagMatVecMult(phiTemp_l, CNForward, phi_l);
        tridiagSolver(phi_l, CNBackward, phiTemp_l);
    }

    // FIELD OPERATORS AGAIN - now reversed in order and l
    for (int li=lMax-1; li >= 0; li--) {  // NB last li term has no mat mult - included only in previous
        cdVec& phi_l = psi[li];
        cdVec& phi_l1 = psi[li+1];

        // MIXING TERM
        // Transform using B (COULD MAYBE JOIN THIS LOOP WITH THE ABOVE TO ENABLE SMOOOTH VECTORIZATION?)
        for (int i=0; i<Nr; i++) {
            phiTemp_l[i] = sqrt2Inv * (phi_l[i] + phi_l1[i]);
            phiTemp_l1[i] = sqrt2Inv * (-phi_l[i] + phi_l1[i]);
        }

        // Apply G^lm operator on two new phi
        double tempFac = dt/4.*Ai*clmVec[li];
        for (int i=0; i<3; i++) {
            for (int j=0; j<Nr; j++) {
                Glm_p[i][j] = M1[i][j] + tempFac * D1[i][j];
                Glm_m[i][j] = M1[i][j] - tempFac * D1[i][j];
            }
        }
        tridiagMatVecMult(phi_l, Glm_m, phiTemp_l);
        tridiagMatVecMult(phi_l1, Glm_p, phiTemp_l1);
        tridiagSolver(phiTemp_l, Glm_p, phi_l);
        tridiagSolver(phiTemp_l1, Glm_m, phi_l1);

        // Transform back using B^T - also brings data back to phi_l and phi_l1
        for (int i=0; i<Nr; i++) {
            phi_l[i] = sqrt2Inv * (phiTemp_l[i] - phiTemp_l1[i]);
            phi_l1[i] = sqrt2Inv * (phiTemp_l[i] + phiTemp_l1[i]);
        }


        // PURE ANGULAR TERM  TODO CHECK THAT LOOPS HERE ARE VECTORIZED!
        // Build the R^lm matrix for each r
        for (int i=0; i<Nr; i++) {
            slmVec[i] = FrVec[i] * (li+1.)*clmVec[li];
        }
        for (int i=0; i<Nr; i++) {
            slmFac1Vec[i] = (1.-slmVec[i]*slmVec[i]) / (1.+slmVec[i]*slmVec[i]);
            slmFac2Vec[i] = 2.*slmVec[i] / (1. + slmVec[i]*slmVec[i]);
        }

        // Apply the 2x2 R^lm matrix for each r
        for (int i=0; i<Nr; i++) {  // TODO Probably need temporary vecs here if I want vecoriation!
            Cdouble temp_l = slmFac1Vec[i]*phi_l[i] - slmFac2Vec[i]*phi_l1[i];
            Cdouble temp_l1 = slmFac2Vec[i]*phi_l[i] + slmFac1Vec[i]*phi_l1[i];
            phi_l[i] = temp_l;
            phi_l1[i] = temp_l1;
        }
    }
}

void RadialPropagator::lengthPropagator(cdMat &psi, double Ei, double dt) {
    // Prebuild some of the R^lm matrix
    for (int i=0; i<Nr; i++) {
        FrVec[i] = Ei*rVec[i] * dt/4.;
    }

    for (int li=0; li<lMax; li++) {  // NB last li term has no mat mult - included only in previous
        cdVec& phi_l = psi[li];
        cdVec& phi_l1 = psi[li+1];

        // PURE ANGULAR TERM  TODO CHECK THAT LOOPS HERE ARE VECTORIZED!
        // Build the R^lm matrix for each r
        for (int i=0; i<Nr; i++) {
            slmVec[i] = FrVec[i] * clmVec[li];
        }
        for (int i=0; i<Nr; i++) {
            slmFac1Vec[i] = (1.-slmVec[i]*slmVec[i]) / (1.+slmVec[i]*slmVec[i]);
            slmFac2Vec[i] = 2.*slmVec[i] / (1. + slmVec[i]*slmVec[i]);
        }

        // Apply the 2x2 R^lm matrix for each r
        for (int i=0; i<Nr; i++) {  // TODO Probably need temporary vecs here if I want vecoriation!
            Cdouble temp_l = slmFac1Vec[i]*phi_l[i] - Cdouble(0.,1.)*slmFac2Vec[i]*phi_l1[i];
            Cdouble temp_l1 = -Cdouble(0.,1.)*slmFac2Vec[i]*phi_l[i] + slmFac1Vec[i]*phi_l1[i];
            phi_l[i] = temp_l;
            phi_l1[i] = temp_l1;
        }
    }

    // ATOMIC HAMILTONIAN
    for (int li=0; li<lMax+1; li++) {
        cdVec& phi_l = psi[li];
        cdMat& CNForward = AtomCNForward[li];
        cdMat& CNBackward = AtomCNBackward[li];
        tridiagMatVecMult(phiTemp_l, CNForward, phi_l);
        tridiagSolver(phi_l, CNBackward, phiTemp_l);
    }

    // FIELD OPERATORS AGAIN - now reversed in order and l
    for (int li=lMax-1; li >= 0; li--) {  // NB last li term has no mat mult - included only in previous
        cdVec& phi_l = psi[li];
        cdVec& phi_l1 = psi[li+1];

        // PURE ANGULAR TERM  TODO CHECK THAT LOOPS HERE ARE VECTORIZED!
        // Build the R^lm matrix for each r
        for (int i=0; i<Nr; i++) {
            slmVec[i] = FrVec[i] * clmVec[li];
        }
        for (int i=0; i<Nr; i++) {
            slmFac1Vec[i] = (1.-slmVec[i]*slmVec[i]) / (1.+slmVec[i]*slmVec[i]);
            slmFac2Vec[i] = 2.*slmVec[i] / (1. + slmVec[i]*slmVec[i]);
        }

        // Apply the 2x2 R^lm matrix for each r
        for (int i=0; i<Nr; i++) {  // TODO Probably need temporary vecs here if I want vecoriation!
            Cdouble temp_l = slmFac1Vec[i]*phi_l[i] - Cdouble(0.,1.)*slmFac2Vec[i]*phi_l1[i];
            Cdouble temp_l1 = -Cdouble(0.,1.)*slmFac2Vec[i]*phi_l[i] + slmFac1Vec[i]*phi_l1[i];
            phi_l[i] = temp_l;
            phi_l1[i] = temp_l1;
        }
    }
}

void RadialPropagator::realTimeProp(cdMat& psi, double dt, double tFin) {
    // Prepare the atomic CN mats with the given time step
    prepareAtomicCNMats(dt);

    int NSteps = static_cast<int>(tFin/dt + 1);
    int printEvery = NSteps / 20;

    std::cout << std::endl << "STARTING REAL TIME PROPAGATION..." << std::endl;
    for (int ti=0; ti<NSteps; ti++) {
        if (ti % printEvery == 0) {std::cout << ti << " / " << NSteps << std::endl;}

        double tEval = (ti+0.5)*dt;  // Time at which to evaluate the time-dependent part of the Hamiltonian

        if (gauge == Gauge::velocity) {
            // Velocity gauge
            double Ai = laserFields.AFieldSin2(tEval);
            velocityPropagator(psi, Ai, dt);
        } else {
            // Length gauge
            double Ei = laserFields.EFieldSin2(tEval);
            lengthPropagator(psi, Ei, dt);
        }

        // IMPLEMENT VARIOUS CALCULATIONS ON PSI(t) HERE
    }

    // Last check the norm of each l component of the WP
    double totNorm = 0.;
    dVec normVec(Nr);
    for (int li=0; li<lMax+1; li++) {
        cdVec& phi_l = psi[li];
        for (int i=0; i<Nr; i++) {
            normVec[i] = std::pow(std::abs(phi_l[i]),2);
        }
        double norm = dr*std::accumulate(normVec.begin(), normVec.end(), 0.);
        std::cout << li << ": " << norm << std::endl;
        totNorm += norm;
    }
    std::cout << "Total norm: " << totNorm << std::endl;
}

void RadialPropagator::groundStateImagProp(dVec &phi0, double dt, int NSteps) {
    // Make sure that the D2 and M2 have Coulomb corrections
    D2[1][0] = deltaCoulomb;
    M2[1][0] = -2.*(1.+dr*dr/12. * deltaCoulomb);

    // Perform the M * V matrix product
    dMat MVprod = dTridiagDiagMult(M2, potMat[0]);

    // Build the two CN matrices
    dMat CNmatForward(3, dVec(phi0.size()));
    dMat CNmatBackward(3, dVec(phi0.size()));
    dMat ECalMat(3, dVec(phi0.size()));  // Matrix for calculating energy
    for (int i=0; i<3; i++) {
        for (int j=0; j<phi0.size(); j++) {
            CNmatForward[i][j] = M2[i][j] - dt/2.*(D2[i][j] + MVprod[i][j]);
            CNmatBackward[i][j] = M2[i][j] + dt/2.*(D2[i][j] + MVprod[i][j]);
            ECalMat[i][j] = D2[i][j] + MVprod[i][j];
        }
    }

    // Now perform imag time steps
    std::cout << "IMAGINARY TIME PROPAGATION FOR THE GROUND STATE:" << std::endl;
    dVec phiCal(phi0.size()), phiEcal(phi0.size());
    for (int i=0; i<NSteps; i++) {
        // Mult the forward CN mat on the wavefunction
        tridiagMatVecMult(phiCal, CNmatForward, phi0);

        // Do backwards solve
        tridiagSolver(phi0, CNmatBackward, phiCal);

        // Renormalize wavefunction to 1
        for (int j=0; j<phi0.size(); j++) {phiCal[j] = phi0[j] * phi0[j];}
        double norm = std::sqrt(dr*std::accumulate(phiCal.begin(), phiCal.end(), 0.));
        for (int j=0; j<phi0.size(); j++) {phi0[j] = phi0[j]/norm;}

        // Do energy calculation every x iteration

        if (i%(NSteps/10) == 0) {
            tridiagSolver(phiCal, M2, phi0);  // Get <\phi| * M^{-1}
            tridiagMatVecMult(phiEcal, ECalMat, phi0);  // Act with Hamiltonian on state
            for (int j=0; j<phi0.size(); j++) {phiCal[j] = phiCal[j] * phiEcal[j];}  // Inner product to get energy
            double E = dr*std::accumulate(phiCal.begin(), phiCal.end(), 0.);

            std::cout << "Energy at step " << i+1 << ": " << E << std::endl;
        }

    }
    tridiagSolver(phiCal, M2, phi0);  // Get <\phi| * M^{-1}
    tridiagMatVecMult(phiEcal, ECalMat, phi0);  // Act with Hamiltonian on state
    for (int j=0; j<phi0.size(); j++) {phiCal[j] = phiCal[j] * phiEcal[j];}  // Inner product to get energy
    double E = std::accumulate(phiCal.begin(), phiCal.end(), 0.);;

    std::cout << "Energy: " << E << std::endl;
}


