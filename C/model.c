#include "model.h"

#include <math.h>

static double S_G_ = 0.025;
static double p2_  = 0.025;
static double p3_  = 1.3e-3;
static double n_   = 0.14;
static double Gb_  = 90.0;
static double Ib_  = 7.0;

EXPORT void set_params(BergmanParams p) {
    S_G_ = p.S_G;
    p2_  = p.p2;
    p3_  = p.p3;
    n_   = p.n;
    Gb_  = p.Gb;
    Ib_  = p.Ib;
}

static void compute_derivatives(double G, double X, double I,
                                double D, double u,
                                double* dG, double* dX, double* dI) {
    *dG = -(S_G_ + X) * (G - Gb_) + D;
    *dX = -p2_ * X + p3_ * (I - Ib_);
    *dI = -n_  * (I - Ib_) + u;
}

EXPORT void step(double* G, double* X, double* I,
                 double Gb, double Ib,
                 double dt, double D, double u) {
    (void)Gb;
    (void)Ib;

    double g = *G;
    double x = *X;
    double ins = *I;

    double k1_g, k1_x, k1_i;
    double k2_g, k2_x, k2_i;
    double k3_g, k3_x, k3_i;
    double k4_g, k4_x, k4_i;

    double dG, dX, dI;

    compute_derivatives(g, x, ins, D, u, &dG, &dX, &dI);
    k1_g = dt * dG;
    k1_x = dt * dX;
    k1_i = dt * dI;

    compute_derivatives(g + 0.5 * k1_g, x + 0.5 * k1_x, ins + 0.5 * k1_i, D, u, &dG, &dX, &dI);
    k2_g = dt * dG;
    k2_x = dt * dX;
    k2_i = dt * dI;

    compute_derivatives(g + 0.5 * k2_g, x + 0.5 * k2_x, ins + 0.5 * k2_i, D, u, &dG, &dX, &dI);
    k3_g = dt * dG;
    k3_x = dt * dX;
    k3_i = dt * dI;

    compute_derivatives(g + k3_g, x + k3_x, ins + k3_i, D, u, &dG, &dX, &dI);
    k4_g = dt * dG;
    k4_x = dt * dX;
    k4_i = dt * dI;

    g   += (k1_g + 2.0 * k2_g + 2.0 * k3_g + k4_g) / 6.0;
    x   += (k1_x + 2.0 * k2_x + 2.0 * k3_x + k4_x) / 6.0;
    ins += (k1_i + 2.0 * k2_i + 2.0 * k3_i + k4_i) / 6.0;

    *G = g;
    *X = x;
    *I = ins;
}

static double appearance(double A, double k, double t) {
    if (A == 0.0) {
        return 0.0;
    }
    if (k <= 0.0) {
        return A;
    }
    return A * exp(-k * t);
}

EXPORT void simulate_dual(double* G_out, double* I_out,
                          int nsteps, double dt,
                          double Afast, double kfast,
                          double Aslow, double kslow,
                          double Aprot, double kprot) {
    double G = Gb_;
    double X = 0.0;
    double I = Ib_;
    double t = 0.0;

    for (int i = 0; i < nsteps; ++i) {
        G_out[i] = G;
        I_out[i] = I;
        if (i == nsteps - 1) {
            break;
        }

        double D1 = appearance(Afast, kfast, t) + appearance(Aslow, kslow, t);
        double u1 = appearance(Aprot, kprot, t);

        double g = G;
        double x = X;
        double ins = I;
        double dG, dX, dI;

        compute_derivatives(g, x, ins, D1, u1, &dG, &dX, &dI);
        double k1_g = dt * dG;
        double k1_x = dt * dX;
        double k1_i = dt * dI;

        double t_half = t + 0.5 * dt;
        double D2 = appearance(Afast, kfast, t_half) + appearance(Aslow, kslow, t_half);
        double u2 = appearance(Aprot, kprot, t_half);
        compute_derivatives(g + 0.5 * k1_g, x + 0.5 * k1_x, ins + 0.5 * k1_i,
                            D2, u2, &dG, &dX, &dI);
        double k2_g = dt * dG;
        double k2_x = dt * dX;
        double k2_i = dt * dI;

        compute_derivatives(g + 0.5 * k2_g, x + 0.5 * k2_x, ins + 0.5 * k2_i,
                            D2, u2, &dG, &dX, &dI);
        double k3_g = dt * dG;
        double k3_x = dt * dX;
        double k3_i = dt * dI;

        double t_full = t + dt;
        double D4 = appearance(Afast, kfast, t_full) + appearance(Aslow, kslow, t_full);
        double u4 = appearance(Aprot, kprot, t_full);
        compute_derivatives(g + k3_g, x + k3_x, ins + k3_i,
                            D4, u4, &dG, &dX, &dI);
        double k4_g = dt * dG;
        double k4_x = dt * dX;
        double k4_i = dt * dI;

        G += (k1_g + 2.0 * k2_g + 2.0 * k3_g + k4_g) / 6.0;
        X += (k1_x + 2.0 * k2_x + 2.0 * k3_x + k4_x) / 6.0;
        I += (k1_i + 2.0 * k2_i + 2.0 * k3_i + k4_i) / 6.0;

        t += dt;
    }
}
