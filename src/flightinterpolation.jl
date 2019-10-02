import PyCall; const py = PyCall
import PyPlot; const plt = PyPlot
import Dierckx; const spl = Dierckx
# Import interpolations from scipy for PCHIP interpolations
# Conda.add("pyplot")
const itp = py.PyNULL()
copy!(itp, py.pyimport_conda("scipy.interpolate", "scipy"))

pchip = itp.PchipInterpolator([1,2,3,4,5,6,7], [0,0,0,3,6,5.7,6])
spline = spl.Spline1D([1,2,3,4,5,6,7], [0,0,0,3,6,5.7,6], bc="extrapolate")
plt.clf()
plt.scatter([1,2,3,4,5,6,7], [0,0,0,3,6,5.7,6], marker="x", color="k", label="data")
plt.plot(0:0.01:8,[pchip(i) for i in 0:0.01:8], label="pchip")
plt.plot(0:0.01:8,spline(0:0.01:8), label="spline")
plt.legend(); plt.grid(ls=":")
plt.gcf()

pwd()
plt.savefig("TrackMatcher.jl/data/flightInterpolation.pdf")
