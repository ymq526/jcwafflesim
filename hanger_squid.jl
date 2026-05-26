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

@variables R Cc Lr Cr Lf Crj Lj Cj Ll Ldc K Lg Rr Rj

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
    ("Rr","3","0",Rr),     # resonator internal loss, Q_int ≈ 80k at f_r

    # Coupling cap: resonator → SQUID
    ("Crj","3","4",Crj),

    # SQUID
    ("L2","4","5",Ll),
    ("L3","4","6",Ll),
    ("Lj1","5","0",Lj),
    ("Cj1","5","0",Cj),
    ("Rj1","5","0",Rj),    # junction 1 internal loss, Q_int ≈ 80k at plasma freq
    ("Lj2","6","0",Lj),
    ("Cj2","6","0",Cj),
    ("Rj2","6","0",Rj),    # junction 2 internal loss

    # Flux bias
    ("L4","7","0",Ldc),
    ("K1","L2","L4",K),
    ("K2","L3","L4",-K),
    ("P3","7","0",3),
    ("R3","7","0",1000.0),
]

# Convert bias current array to flux in units of Φ₀.
# M = K*sqrt(Ll*Ldc) is the mutual inductance between the bias coil and SQUID arm.
function current_to_flux(I; K, Ll, Ldc)
    Phi0 = 2.067833848e-15   # magnetic flux quantum h/(2e), Wb
    M    = K * sqrt( Ll * Ldc)
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
    Cj   => 100e-15,      # junction capacitance (each)
    Ll   => 34e-12,       # SQUID arm inductance
    K    => 0.999,        # mutual coupling (K=1.0 makes inv-L matrix singular)
    Ldc  => 0.74e-12,     # flux bias coil inductance
    Lg   => 100e-9,       # large shunt L (~4 kΩ at 6.8 GHz, negligible loading)
    Rr   => 2.55e6,       # resonator internal Q ≈ 80k at f_r ≈ 6.87 GHz
    Rj   => 8.0e6,        # junction internal Q ≈ 80k at plasma freq ≈ 15.9 GHz
)

# S21 heatmap: frequency × DC bias current
# Fine sweep over the tuning range found with Lj=1nH: 6.72–6.82 GHz (~93 MHz span)
ws           = 2*pi*(6.76:0.0002:6.85)*1e9   # 1 MHz steps, 181 freq pts
wp           = (2*pi*13.7*1e9,)              # placeholder pump freq; Ip=0 throughout
currentvals  = (-40:0.1:40)*1e-5             # ±400 μA in 5 μA steps, 161 bias pts

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
savefig(p, raw"C:\Users\Carissa K\Documents\jcwafflesim\hanger_squid_s21.png")

# ── Two-tone spectroscopy ──────────────────────────────────────────────────
# Park signal at the notch; sweep pump through port 3.
# For each (flux, pump_freq), compute |S21_pump - S21_ref| in the complex
# S21 plane, where S21_ref is taken from the pump-off sweep above.

notch_idx   = [argmin(abs.(outvals[:,k])) for k in 1:size(outvals,2)]
notch_freqs = ws[notch_idx]

# Reference S21 at the notch per flux point — reuse the Ip=0 sweep
S21_ref_all = [outvals[notch_idx[k], k] for k in 1:length(currentvals)]

# Subsample flux axis — every 4th of 801 pts → 201 flux points
step_2t = 4
idx_2t  = collect(1:step_2t:length(currentvals))
nf_2t   = notch_freqs[idx_2t]
flux_2t = fluxvals[idx_2t]

wp_2t_GHz = 0:0.002:12.0             # 2 MHz steps → 6001 pump pts
wp_2t     = 2*pi * collect(wp_2t_GHz) * 1e9
Ip_2t     = 800e-9

out_2t = fill(NaN, length(wp_2t), length(idx_2t))  # |ΔS21| linear

println("Two-tone (|ΔS21|): $(length(idx_2t)) flux pts × $(length(wp_2t)) pump pts …")
@time for (ki, k) in enumerate(idx_2t)
    Idc     = currentvals[k]
    ws_park = [nf_2t[ki]]
    S21_ref = S21_ref_all[k]
    for (j, wp_j) in enumerate(wp_2t)
        src = [
            (mode=(0,), port=3, current=Idc),
            (mode=(1,), port=3, current=Ip_2t),
        ]
        try
            sol = hbsolve(ws_park, (wp_j,), src, (1,), (1,),
                circuit, circuitdefs; dc=true, threewavemixing=true, fourwavemixing=true)
            out_2t[j,ki] = abs(sol.linearized.S((0,),2,(0,),1,1) - S21_ref)
        catch
            # singular at this pump freq — leave NaN
        end
    end
end

p2 = plot(
    flux_2t,
    collect(wp_2t_GHz),
    out_2t,
    seriestype  = :heatmap,
    xlabel      = "DC flux bias (Φ/Φ₀)",
    ylabel      = "Pump frequency (GHz)",
    title       = "Two-tone |ΔS21|: complex distance from pump-off reference",
    colorbar_title = "|ΔS21|",
    clims       = (0, 1),
    size        = (800, 500),
)
savefig(p2, raw"C:\Users\Carissa K\Documents\jcwafflesim\hanger_squid_2tone.png")
