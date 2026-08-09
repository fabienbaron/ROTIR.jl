import numpy as np, libphoebe, os, sys
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
ld = lambda f: np.ascontiguousarray(np.loadtxt(os.path.join(D, f)), dtype=np.float64)
V  = [ld("V1.txt"), ld("V2.txt")]
Tr = [np.ascontiguousarray(np.loadtxt(os.path.join(D,f)), dtype=np.int32) for f in ("Tr1.txt","Tr2.txt")]
N  = [ld("N1.txt"), ld("N2.txt")]
A  = [ld("A1.txt"), ld("A2.txt")]
R  = [ld("R1.txt"), ld("R2.txt")]
F0 = [ld("F0_1.txt"), ld("F0_2.txt")]
print("triangles:", [len(t) for t in Tr], "vertices:", [len(v) for v in V]); sys.stdout.flush()

arr = lambda *v: np.ascontiguousarray(v, dtype=np.float64)
cases = (("uniform", [(b"uniform", arr()), (b"uniform", arr())]),
         ("linear",  [(b"linear",  arr(0.6)), (b"linear", arr(0.6))]))
for tag, LD in cases:
    for meth in ("Horvat", "Wilson"):
        F = libphoebe.mesh_radiosity_problem_nbody_convex(
                V, Tr, N, A, R, F0, LD, meth.encode(), b"triangles",
                epsF=1e-12, max_iter=100)
        for k in (0, 1):
            jl = ld("jl_%s_%s_%d.txt" % (tag, meth.lower(), k+1))
            ph = np.asarray(F[k], dtype=np.float64)
            rj, rp = jl - F0[k], ph - F0[k]        # the reflected part
            rel = np.abs(rj - rp).max() / max(np.abs(rp).max(), 1e-30)
            print("  %-8s %-7s body%d: max|Δ(reflected)|/max = %9.3e | "
                  "max Fout/F0  julia %.9f  phoebe %.9f"
                  % (tag, meth, k+1, rel, (jl/F0[k]).max(), (ph/F0[k]).max()))
            sys.stdout.flush()
