#ifndef MODEL_H
#define MODEL_H
typedef struct {
    double S_G_min1;
    double p2_min1;
    double p3_min1_per_uU_per_mL;
    double phi_G_uU_mL_min1_per_mg_dL;
    double G_thr_mg_dL;
    double n_min1;
    double Gb_mg_dL;
    double Ib_uU_mL;
    double A;
    double k;
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
