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

    dydt[0] = -(p->S_G_min1 + X) * G + p->S_G_min1 * p->Gb_mg_dL + Dt;             // dG/dt
    dydt[1] = -p->p2_min1 * X + p->p3_min1_per_uU_per_mL * (I - p->Ib_uU_mL);      // dX/dt
    const double secretion = (G > p->G_thr_mg_dL)
        ? p->phi_G_uU_mL_min1_per_mg_dL * (G - p->G_thr_mg_dL)
        : 0.0;
    dydt[2] = -p->n_min1 * (I - p->Ib_uU_mL) + secretion;                          // dI/dt
}

static void derivatives_ex(
    double t, const double y[], double dydt[], const Params *p,
    double Afast, double kfast,
    double Aslow, double kslow,
    double Aprot, double kprot)
{
    const double G = y[0];
    const double X = y[1];
    const double I = y[2];

    const double fast = Afast * exp(-kfast * t);
    const double slow = Aslow * exp(-kslow * t);
    const double D = fast + slow;

    const double secretion = (G > p->G_thr_mg_dL)
        ? p->phi_G_uU_mL_min1_per_mg_dL * (G - p->G_thr_mg_dL)
        : 0.0;
    const double Iprot = Aprot * exp(-kprot * t);

    dydt[0] = -(p->S_G_min1 + X) * G + p->S_G_min1 * p->Gb_mg_dL + D;
    dydt[1] = -p->p2_min1 * X + p->p3_min1_per_uU_per_mL * (I - p->Ib_uU_mL);
    dydt[2] = -p->n_min1 * (I - p->Ib_uU_mL) + secretion + Iprot;
}

// Simple fixed-step RK4. Outputs G and I only (X is internal).
API void simulate(double *G_out, double *I_out,
                  int nsteps, double dt, const Params *p)
{
    double t = 0.0;
    double y[3] = { p->Gb_mg_dL, 0.0, p->Ib_uU_mL }; // start at basal, X=0

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

API void simulate_ex(
    double* G_out, double* I_out, int nsteps, double dt,
    const Params* prm,
    double Afast, double kfast,
    double Aslow, double kslow,
    double Aprot, double kprot)
{
    double t = 0.0;
    double y[3] = { prm->Gb_mg_dL, 0.0, prm->Ib_uU_mL };

    double k1[3], k2[3], k3[3], k4[3], ytmp[3], dydt[3];

    for (int i = 0; i < nsteps; ++i) {
        G_out[i] = y[0];
        I_out[i] = y[2];

        derivatives_ex(t, y, dydt, prm, Afast, kfast, Aslow, kslow, Aprot, kprot);
        for (int j = 0; j < 3; ++j) {
            k1[j] = dt * dydt[j];
            ytmp[j] = y[j] + 0.5 * k1[j];
        }

        derivatives_ex(t + 0.5 * dt, ytmp, dydt, prm, Afast, kfast, Aslow, kslow, Aprot, kprot);
        for (int j = 0; j < 3; ++j) {
            k2[j] = dt * dydt[j];
            ytmp[j] = y[j] + 0.5 * k2[j];
        }

        derivatives_ex(t + 0.5 * dt, ytmp, dydt, prm, Afast, kfast, Aslow, kslow, Aprot, kprot);
        for (int j = 0; j < 3; ++j) {
            k3[j] = dt * dydt[j];
            ytmp[j] = y[j] + k3[j];
        }

        derivatives_ex(t + dt, ytmp, dydt, prm, Afast, kfast, Aslow, kslow, Aprot, kprot);
        for (int j = 0; j < 3; ++j) {
            k4[j] = dt * dydt[j];
        }

        for (int j = 0; j < 3; ++j) {
            y[j] += (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]) / 6.0;
        }

        t += dt;
    }
}
