//
// Created by au642261 on 8/3/25.
//

#ifndef FIELDS_H
#define FIELDS_H

#include <cmath>

class Fields {
public:
    double omega, CEP;
    int NCycles;
    double tFinSin2, Up, sqrtUp, A0, Eamp;

    Fields(){}

    Fields(double lambda_nm, double I0_Wcm2, double CEP, int NCycles) : CEP(CEP), NCycles(NCycles) {
        // Convert to atomic units
        Eamp = std::sqrt(I0_Wcm2 / 3.50945e16);  // Maximum electric field amplitude in a.u.
        omega = 2.*M_PI * 137.036 / (lambda_nm * 1.0e-9 / 5.29177e-11);

        // Calculate various variables
        Up = Eamp*Eamp / (4. * omega*omega);  // Ponderomotive energy in a.u.
        sqrtUp = std::sqrt(Up);
        A0 = 2.*std::sqrt(Up);  // Field amplitude in a.u.
        tFinSin2 = 2.*M_PI*NCycles/omega;
    };

    double AFieldSin2(double t) {
        if (t >= 0. && t <= tFinSin2) {
            return A0 * std::pow(std::sin(omega*t/(2.*NCycles)),2) * std::cos(omega*t + CEP);
        } else {
            return 0.;
        }
    }

    // TODO IMPLEMENT EFieldSin2
    double EFieldSin2(double t) {}
};

#endif //FIELDS_H
