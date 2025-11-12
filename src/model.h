#ifndef MODEL_H
#define MODEL_H

#ifdef _WIN32
  #ifdef BUILDING_MODEL
    #define EXPORT __declspec(dllexport)
  #else
    #define EXPORT __declspec(dllimport)
  #endif
#else
  #define EXPORT __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double S_G;
    double p2;
    double p3;
    double n;
    double Gb;
    double Ib;
} BergmanParams;

EXPORT void set_params(BergmanParams p);
EXPORT void step(double* G, double* X, double* I,
                 double Gb, double Ib,
                 double dt, double D, double u);

#ifdef __cplusplus
}
#endif

#endif /* MODEL_H */
