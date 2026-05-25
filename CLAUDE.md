# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Julia scripts for simulating superconducting microwave circuits (JPAs, resonators) with [JosephsonCircuits.jl](https://josephsoncircuits.org). Runs on the global Julia v1.12 environment — there is no `Project.toml`.

## Running scripts and viewing plots

- Run with `julia path/to/script.jl`. The first call in a fresh process recompiles packages (~5 s); this is normal.
- Plots do NOT display when Julia runs from the terminal. To see a figure, add `savefig(p, "name.png")` and then open the PNG (read the image file). The scripts here build a plot object but don't save it by default.
- Frequencies are angular throughout: `ws = 2*pi*f`. Convert back for plotting with `ws/(2*pi*1e9)` (GHz).

## JosephsonCircuits.jl modeling rules

These are easy to get wrong and cause errors or silently bad results:

- **Allowed netlist components only:** `Lj, NL, L, C, K, I, R, P`. There is no transmission-line element — model feedlines and distributed (λ/2, λ/4) resonators as lumped equivalents (e.g. a parallel `L`/`C` tank, a near-ideal short as a tiny series `L`).
- **One port `P` per branch (node pair).** Two ports cannot share the same node pair ("Only one port allowed per branch"). Put them on distinct nodes and join with a through element.
- **Each port `P` needs a co-located resistor `R`** on the same node pair; that resistor's value is the S-parameter reference impedance (and the source/load impedance), default 50 Ω.
- **Every node needs a DC path to ground.** A purely capacitively-coupled node floats; add a large inductor to ground if needed.
- **Solver choice:** `hblinsolve(ws, circuit, circuitdefs)` for purely linear circuits (no junctions); `hbsolve(...)` for pumped/nonlinear circuits with Josephson junctions. `hblinsolve` returns the linearized object directly (`sol.S(...)`); `hbsolve` returns it nested under `sol.linearized.S(...)`.
- **S-parameter indexing:** `sol.S(outputmode, outputport, inputmode, inputport, freqindex)`. The signal mode is `(0,)`. Example — S21: `sol.S((0,),2,(0,),1,:)`.
