using JosephsonCircuits
using Plots

# Hanger- (notch-) type resonator measured in transmission (S21).
#
# A 50 Ohm feedline runs between port 1 (node 1) and port 2 (node 2). The package
# allows only one port per branch, so the two ports sit on distinct nodes joined
# by Lf, a tiny series inductance that models the electrically-short feedline as
# a near-ideal short (omega*Lf << 50 Ohm over the band). A lambda/2 CPW resonator
# (Z0 = 50 Ohm) is represented by its effective lumped parallel tank Lr || Cr and
# couples to the feedline through Cc. Near resonance the series combination
# Cc + (Lr||Cr) shorts the feedline to ground, producing a notch in S21 at
# f0 = 1/(2*pi*sqrt(Lr*(Cr+Cc))).

@variables R Cc Lr Cr Lf
circuit = [
    ("P1","1","0",1),     # port 1 on the feedline
    ("R1","1","0",R),     # 50 Ohm reference/source impedance
    ("P2","2","0",2),     # port 2 on the feedline
    ("R2","2","0",R),     # 50 Ohm reference/load impedance
    ("Lf","1","2",Lf),    # feedline modeled as ideal short (tiny series inductance)
    ("Cc","1","3",Cc),    # coupling capacitor feedline <-> resonator
    ("Lr","3","0",Lr),    # resonator effective inductance
    ("Cr","3","0",Cr),    # resonator effective capacitance
]

circuitdefs = Dict(
    R  => 50.0,
    Lf => 10e-12,        # negligible feedline inductance (ideal short)
    Cc => 5.2e-15,       # coupling capacitance (sets external Q) - adjust as needed
    Lr => 0.738e-9,       # resonator effective inductance
    Cr => 728e-15,        # resonator effective capacitance
)

# signal frequency sweep around the resonance (~6.82 GHz)
ws = 2*pi*(6.5:0.0002:7.2)*1e9

sol = hblinsolve(ws, circuit, circuitdefs)

# complex S21 (output port 2, input port 1, signal mode (0,))
S21 = sol.S((0,),2,(0,),1,:)

p1 = plot(
    ws/(2*pi*1e9),
    20*log10.(abs.(S21)),
    xlabel="Signal frequency (GHz)",
    ylabel="|S21| (dB)",
    title="Hanger resonator transmission",
    legend=false,
)

p2 = plot(
    ws/(2*pi*1e9),
    rad2deg.(angle.(S21)),
    xlabel="Signal frequency (GHz)",
    ylabel="arg(S21) (deg)",
    legend=false,
)

plot(p1, p2, layout=(2,1), size=(700,600))
