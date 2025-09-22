#ifndef MODEL_H
#define MODEL_H

typedef struct {
    double p1, p2, p3, p4, p5, p6;
    double Gb, Ib;
    double A, k; // dessert input params
} Params;

// Derivative function: y = [G, X, I]
void derivatives(double t, const double y[], double dydt[], const Params *params);

#endif
