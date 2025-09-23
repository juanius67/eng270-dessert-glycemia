import os, sys
import ctypes as C
import numpy as np

HERE = os.path.dirname(__file__)
LIB_DIR = os.path.join(HERE, "..", "C")

def _try_load_lib():
    # try Linux/Mac .so, then Windows .dll
    for name in ("libmodel.so", "model.dll"):
        path = os.path.join(LIB_DIR, name)
        if os.path.exists(path):
            lib = C.CDLL(path)
            return lib
    return None

# numpy-based RK4 fallback (so you can run without compiling C yet)
def _python_rk4(params, A, k, dt, nsteps):
    G = np.empty(nsteps); I = np.empty(nsteps)
    Gb, Ib = params["Gb"], params["Ib"]
    p1,p2,p3,p4,p5,p6 = [params[k_] for k_ in ("p1","p2","p3","p4","p5","p6")]
    y = np.array([Gb, 0.0, Ib], dtype=float)
    t = 0.0
    def deriv(t, y):
        G_, X_, I_ = y
        Dt = A * np.exp(-k*t)
        dG = -(p1 + X_) * G_ + p1*Gb + Dt
        dX = -p2 * X_ + p3 * (I_ - Ib)
        sec = p4 * (G_ - p5) if G_ > p5 else 0.0
        dI = -p6 * (I_ - Ib) + sec
        return np.array([dG,dX,dI])
    for i in range(nsteps):
        G[i], I[i] = y[0], y[2]
        k1 = dt*deriv(t, y)
        k2 = dt*deriv(t+0.5*dt, y+0.5*k1)
        k3 = dt*deriv(t+0.5*dt, y+0.5*k2)
        k4 = dt*deriv(t+dt, y+k3)
        y = y + (k1 + 2*k2 + 2*k3 + k4)/6.0
        t += dt
    return G, I

def run_simulation(params, A, k):
    dt = float(params["dt"]); t_end = float(params["t_end"])
    nsteps = int(t_end/dt) + 1
    lib = _try_load_lib()
    if lib is None:
        # fallback
        G, I = _python_rk4(params, A, k, dt, nsteps)
        t = np.linspace(0.0, t_end, nsteps)
        return {"t": t, "G": G, "I": I, "source": "python"}
    # set signature
    class ParamsC(C.Structure):
        _fields_ = [("p1", C.c_double),("p2", C.c_double),("p3", C.c_double),
                    ("p4", C.c_double),("p5", C.c_double),("p6", C.c_double),
                    ("Gb", C.c_double),("Ib", C.c_double),
                    ("A", C.c_double),("k", C.c_double)]
    lib.simulate.argtypes = [C.POINTER(C.c_double), C.POINTER(C.c_double),
                             C.c_int, C.c_double, C.POINTER(ParamsC)]
    lib.simulate.restype = None

    G = (C.c_double * nsteps)()
    I = (C.c_double * nsteps)()
    pc = ParamsC(params["p1"],params["p2"],params["p3"],
                 params["p4"],params["p5"],params["p6"],
                 params["Gb"],params["Ib"], A, k)
    lib.simulate(G, I, nsteps, dt, C.byref(pc))
    G = np.frombuffer(G, dtype=np.float64, count=nsteps)
    I = np.frombuffer(I, dtype=np.float64, count=nsteps)
    t = np.linspace(0.0, t_end, nsteps)
    return {"t": t, "G": G, "I": I, "source": "C"}
