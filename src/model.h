#ifndef MODEL_H
#define MODEL_H

/*
 * This header declares the small C interface that Python loads via ctypes.
 * This document was written with the help of generative AI tools
 * (ChatGPT / ChatGPT Codex and “Jules” from Google AI) to:
 *   - standardise cross-platform export macros across windows and POSIX,
 *   - generate clear, consistent documentation comments, and
 *   - reduce time spent on compatability issues so that human effort could focus on
 *     the mathematical model design, parameter choices, and analysis.
 *
 */

#ifdef _WIN32
  #ifdef BUILDING_MODEL
    /* When building the DLL on Windows */
    #define EXPORT __declspec(dllexport)
  #else
    /* When consuming the DLL on Windows */
    #define EXPORT __declspec(dllimport)
  #endif
#else
  /*
   * On POSIX systems (Linux, macOS), make the symbols visible from shared
   * libraries by default so ctypes can locate them.
   */
  #define EXPORT __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Parameters for the Bergman minimal model.
 *
 * This struct groups all physiological constants that characterise a
 * nominal subject.  The units are chosen to match the Python layer:
 *
 *   - Glucose in mg/dL
 *   - Insulin in µU/mL
 *   - Time in minutes
 *
 * These values are typically loaded from configs/params.yaml on the
 * Python side and passed once to set_params() at the start of a run.
 *
 * @note Struct definition and documentation text generated with AI
 *       (ChatGPT / Jules) and then reviewed by the project author.
 */
typedef struct {
    double S_G; /**< Glucose effectiveness (min^-1). Controls how fast G decays toward Gb in the absence of insulin. */
    double p2;  /**< Decay rate of remote insulin effect X (min^-1). */
    double p3;  /**< Gain from insulin above baseline (I - Ib) to remote effect X (min^-1 per µU/mL). */
    double n;   /**< Insulin clearance rate (min^-1). */
    double Gb;  /**< Basal (fasting) glucose concentration (mg/dL). */
    double Ib;  /**< Basal (fasting) insulin concentration (µU/mL). */
} BergmanParams;

/**
 * @brief Set global parameters for the Bergman minimal model.
 *
 * The implementation stores the contents of @p p in internal static
 * variables.  Subsequent calls to step() use these stored values and
 * therefore do not need to pass parameters on every time step.
 *
 * Typical usage (from C):
 * @code
 *   BergmanParams p = { ... };  // filled from YAML on the Python side
 *   set_params(p);
 * @endcode
 *
 * @param p  Structure containing all model parameters.
 *
 * @note Function signature and comments generated with AI and kept
 *       small on purpose so that ctypes bindings remain simple.
 */
EXPORT void set_params(BergmanParams p);

/**
 * @brief Advance the Bergman minimal model by one RK4 time step.
 *
 * This function integrates the coupled ODE system for glucose G,
 * remote insulin effect X, and insulin I over a single step of size
 * @p dt using a classical Runge–Kutta 4 scheme.
 *
 * All state variables are updated in place so that the caller can run
 * a full simulation with a simple loop:
 *
 * @code
 *   for (int k = 0; k < n_steps; ++k) {
 *       double t = k * dt;
 *       double D = appearance(t);   // gut appearance term (mg/dL/min)
 *       double u = insulin_input(t); // insulin stimulus (µU/mL/min)
 *       step(&G, &X, &I, Gb, Ib, dt, D, u);
 *   }
 * @endcode
 *
 * @param G   Pointer to current glucose concentration (mg/dL). Updated in place.
 * @param X   Pointer to current remote insulin effect (min^-1). Updated in place.
 * @param I   Pointer to current insulin concentration (µU/mL). Updated in place.
 * @param Gb  Basal glucose concentration (mg/dL) used in the ODEs.
 * @param Ib  Basal insulin concentration (µU/mL) used in the ODEs.
 * @param dt  Time step (minutes).
 * @param D   Glucose appearance rate from the gut (mg/dL/min).
 * @param u   Insulin secretion stimulus (µU/mL/min).
 *
 * @note Documentation block drafted with AI; numerical behaviour and
 *       use in the project were specified and checked by the author.
 */
EXPORT void step(double* G, double* X, double* I,
                 double Gb, double Ib,
                 double dt, double D, double u);

#ifdef __cplusplus
}
#endif

#endif /* MODEL_H */
