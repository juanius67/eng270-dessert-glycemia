#ifndef MODEL_H
#define MODEL_H

#ifdef _WIN32
  #define EXPORT __declspec(dllexport)
#else
  #define EXPORT
#endif

typedef struct {
    double S_G;  /* min^-1 */
    double p2;   /* min^-1 */
    double p3;   /* (min^-1)/(uU mL^-1) */
    double n;    /* min^-1 */
    double Gb;   /* mg/dL */
    double Ib;   /* uU/mL */
} BergmanParams;

EXPORT void step(double* G, double* X, double* I,
                 double Gb, double Ib,
                 double dt, double D, double u);

EXPORT void set_params(BergmanParams p);

EXPORT void simulate_dual(double* G_out, double* I_out,
                          int nsteps, double dt,
                          double Afast, double kfast,
                          double Aslow, double kslow,
                          double Aprot, double kprot);

#endif /* MODEL_H */
