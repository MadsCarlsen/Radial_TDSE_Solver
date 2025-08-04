#include "RadialPropagator.h"
#include "Types.h"

int main() {
    settingsStruct settings{};

    settings.gauge = Gauge::velocity;

    // Grid settings
    settings.Nr = 6000;
    settings.dr = 0.1;
    settings.lMax = 25;

    // Field settings
    settings.NCycles = 2;
    settings.I0_Wcm2 = 1.3e14;
    settings.lambda_nm = 800.;
    settings.CEP = M_PI/2.;

    RadialPropagator radProp(settings);

    // Get the ground state
    dVec phi0(settings.Nr, 1.);
    radProp.groundStateImagProp(phi0, 2., 200);

    // Transfer and propagate
    cdMat psi(settings.lMax+1, cdVec(settings.Nr, Cdouble(0.,0.)));
    for (int i=0; i<settings.Nr; i++) {psi[0][i] = phi0[i];}
    radProp.realTimeProp(psi, 0.05, radProp.laserFields.tFinSin2);

}