#ifndef MODEL_H
#define MODEL_H
typedef struct {
    double p1, p2, p3, p4, p5, p6;
    double Gb, Ib;
    double A, k; // dessert input params
} Params;
void derivatives(double t, const double y[], double dydt[], const Params *params);
#ifdef _WIN32
#define API __declspec(dllexport)
#else
#define API
#endif
API void simulate(double *G_out, double *I_out,
                  int nsteps, double dt, const Params *p);

// Extended dual-exponential appearance with protein-driven insulin pulse.
// D(t) = Afast*exp(-kfast t) + Aslow*exp(-kslow t)
// extra insulin drive from protein: Iprot(t) = Aprot*exp(-kprot t)
API void simulate_ex(
    double* G_out, double* I_out, int nsteps, double dt,
    const Params* prm,
    double Afast, double kfast,
    double Aslow, double kslow,
    double Aprot, double kprot);
#endif
