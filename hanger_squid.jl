using JosephsonCircuits
using Plots

# Hanger resonator (feedline P1→P2, resonator Lr||Cr coupled via Cc)
# capacitively coupled via Crj to a flux-tunable SQUID (two JJs in a loop).
# DC flux bias applied through mutual inductance K on port 3.
# Sweeps S21 vs signal frequency and DC bias current.
#
# Node map:
#   1 — feedline left  (P1, R1)
#   2 — feedline right (P2, R2)
#   3 — resonator      (Cc from 1, Lr||Cr to gnd, Crj to 4)
#   4 — SQUID top      (Crj from 3, Lj1||Cj1 to gnd, Ll to 5)
#   5 — SQUID bottom   (Lj2||Cj2 to gnd)
#   6 — flux bias coil (Ldc to gnd, K with L2, P3, R3=1kΩ)

@variables R Cc Lr Cr Lf Crj Lj Cj Ll Ldc K Lg

circuit = [
    # Feedline
    ("P1","1","0",1),
    ("R1","1","0",R),
    ("Lg1","1","0",Lg),   # large shunt L: anchors DC node flux at node 1 (resistors alone don't)
    ("P2","2","0",2),
    ("R2","2","0",R),
    ("Lg2","2","0",Lg),   # same for node 2
    ("Lf","1","2",Lf),

    # Hanger resonator
    ("Cc","1","3",Cc),
    ("Lr","3","0",Lr),
    ("Cr","3","0",Cr),

    # Coupling cap: resonator → SQUID
    ("Crj","3","4",Crj),

    # SQUID
    ("Lj1","4","0",Lj),
    ("Cj1","4","0",Cj),
    ("L2","4","5",Ll),
    ("Lj2","5","0",Lj),
    ("Cj2","5","0",Cj),

    # Flux bias
    ("L3","6","0",Ldc),
    ("K1","L2","L3",K),
    ("P3","6","0",3),
    ("R3","6","0",1000.0),
]

# Convert bias current array to flux in units of Φ₀.
# M = K*sqrt(Ll*Ldc) is the mutual inductance between the bias coil and SQUID arm.
function current_to_flux(I; K, Ll, Ldc)
    Phi0 = 2.067833848e-15   # magnetic flux quantum h/(2e), Wb
    M    = 0.5 * K * sqrt(Ll * Ldc)
    return M .* I ./ Phi0
end

circuitdefs = Dict(
    R    => 50.0,
    Lf   => 10e-12,       # near-ideal feedline short
    Cc   => 5.2e-15,      # feedline→resonator coupling
    Lr   => 0.738e-9,     # resonator inductance (~6.84 GHz bare)
    Cr   => 728e-15,      # resonator capacitance
    Crj  => 7.5e-15,      # resonator→SQUID coupling
    Lj   => 1e-9,   # junction inductance (each)
    Cj   => 10.0e-15,     # junction capacitance (each)
    Ll   => 34e-12,       # SQUID arm inductance
    K    => 0.999,        # mutual coupling (K=1.0 makes inv-L matrix singular)
    Ldc  => 0.74e-12,     # flux bias coil inductance
    Lg   => 100e-9,       # large shunt L (~4 kΩ at 6.8 GHz, negligible loading)
)

# S21 heatmap: frequency × DC bias current
# Fine sweep over the tuning range found with Lj=1nH: 6.72–6.82 GHz (~93 MHz span)
ws           = 2*pi*(6.75:0.0005:6.85)*1e9   # 1 MHz steps, 181 freq pts
wp           = (2*pi*13.7*1e9,)              # placeholder pump freq; Ip=0 throughout
currentvals  = (-100:0.5:100)*1e-5             # ±400 μA in 5 μA steps, 161 bias pts

Npumpharmonics       = (1,)
Nmodulationharmonics = (1,)

outvals = zeros(Complex{Float64}, length(ws), length(currentvals))

println("Sweeping $(length(currentvals)) bias points × $(length(ws)) freq points …")
@time for (k, Idc) in enumerate(currentvals)
    sources = [
        (mode=(0,), port=3, current=Idc),
        (mode=(1,), port=3, current=0.0),
    ]
    sol = hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        circuit, circuitdefs; dc=true, threewavemixing=true, fourwavemixing=true)
    outvals[:,k] = sol.linearized.S((0,),2,(0,),1,:)
end

fluxvals = current_to_flux(currentvals;
    K   = circuitdefs[K],
    Ll  = circuitdefs[Ll],
    Ldc = circuitdefs[Ldc],
)

p = plot(
    fluxvals,
    ws/(2*pi*1e9),
    10*log10.(abs2.(outvals)),
    seriestype = :heatmap,
    xlabel     = "DC flux bias (Φ/Φ₀)",
    ylabel     = "Signal frequency (GHz)",
    title      = "S21 (dB): hanger resonator + flux-tunable SQUID",
    clims      = (-15, 0),
    size       = (800, 500),
)
