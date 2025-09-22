#include "model.h"
#include <math.h>

// state vector y = [G, X, I]
void derivatives(double t, const double y[], double dydt[], const Params *p) {
    double G = y[0], X = y[1], I = y[2];

    // Dessert input D(t) = A * exp(-k*t)
    double Dt = p->A * exp(-p->k * t);

    // Glucose dynamics
    dydt[0] = -(p->p1 + X) * G + p->p1 * p->Gb + Dt;
    // Insulin action dynamics
    dydt[1] = -p->p2 * X + p->p3 * (I - p->Ib);
    // Insulin dynamics with secretion
    double secretion = (G > p->p5) ? p->p4 * (G - p->p5) : 0.0;
    dydt[2] = -p->p6 * (I - p->Ib) + secretion;
}
