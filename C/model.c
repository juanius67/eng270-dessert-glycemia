#include "model.h"
#include <math.h>
#include <stddef.h>

#ifdef _WIN32
#define API __declspec(dllexport)
#else
#define API
#endif

// y = [G, X, I]
void derivatives(double t, const double y[], double dydt[], const Params *p) {
    (void)t; // not used except in D(t)
    const double G = y[0], X = y[1], I = y[2];
    const double Dt = p->A * exp(-p->k * t);     // D(t) = A e^{-kt}

    dydt[0] = -(p->p1 + X) * G + p->p1 * p->Gb + Dt;             // dG/dt
    dydt[1] = -p->p2 * X + p->p3 * (I - p->Ib);                  // dX/dt
    const double secretion = (G > p->p5) ? p->p4 * (G - p->p5) : 0.0;
    dydt[2] = -p->p6 * (I - p->Ib) + secretion;                  // dI/dt
}

// Simple fixed-step RK4. Outputs G and I only (X is internal).
API void simulate(double *G_out, double *I_out,
                  int nsteps, double dt, const Params *p)
{
    double t = 0.0;
    double y[3] = { p->Gb, 0.0, p->Ib }; // start at basal, X=0

    // scratch
    double k1[3], k2[3], k3[3], k4[3], ytmp[3], dydt[3];

    for (int i = 0; i < nsteps; ++i) {
        // store outputs (include t=0 sample at index 0)
        G_out[i] = y[0];
        I_out[i] = y[2];

        // RK4
        derivatives(t, y, dydt, p);
        for (int j=0;j<3;++j) k1[j] = dt * dydt[j];

        for (int j=0;j<3;++j) ytmp[j] = y[j] + 0.5*k1[j];
        derivatives(t + 0.5*dt, ytmp, dydt, p);
        for (int j=0;j<3;++j) k2[j] = dt * dydt[j];

        for (int j=0;j<3;++j) ytmp[j] = y[j] + 0.5*k2[j];
        derivatives(t + 0.5*dt, ytmp, dydt, p);
        for (int j=0;j<3;++j) k3[j] = dt * dydt[j];

        for (int j=0;j<3;++j) ytmp[j] = y[j] + k3[j];
        derivatives(t + dt, ytmp, dydt, p);
        for (int j=0;j<3;++j) k4[j] = dt * dydt[j];

        for (int j=0;j<3;++j)
            y[j] += (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]) / 6.0;

        t += dt;
    }
}
