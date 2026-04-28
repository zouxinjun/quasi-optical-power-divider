# THz Power Divider Reflector Optimization

MATLAB implementation of a Physical Optics (PO)-based electromagnetic field solver and Differential Evolution (DE) optimizer for an 8-way terahertz (380 GHz) power divider reflector surface.

## Overview

The system models a reflector-based power divider where a horn antenna illuminates a perturbed parabolic reflector (H2 plane), and the resulting field is evaluated at an observation plane (H3 plane). The optimizer adjusts the reflector surface perturbation to maximize the correlation between the computed H3 field and a target 8-beam field pattern.

**Computation flow:**

```
Horn aperture (H1)
      │  free-space propagation (Helmholtz-Kirchhoff)
Reflector surface (H2)  ← parabolic baseline + perturbation
      │  free-space propagation
Observation plane (H3)  ← compare with target pattern
```

## Repository Structure

```
├── Generate_Target_Field/
│   └── gen_target_field.m          # Generate the 8-beam target H-field pattern
│
├── Generate_Window_Function/
│   └── gen_weight_window.m         # Generate the weighting window for fitness evaluation
│
├── Generate_Reflector_Perturbation/
│   └── generate_divider_surface.m  # Generate the initial reflector perturbation file
│
├── Generate_STL/
│   ├── generate_reflector_stl.m    # Export the reflector surface as an STL file
│   └── stlwrite.m                  # STL file writer utility
│
├── Single_Run/
│   └── compute_H_from_surface.m    # Standalone H-field solver (single evaluation + plots)
│
└── Optimize_Reflector/
    ├── main_optimize.m             # Top-level optimization script (run this)
    ├── optimize_reflector_surface.m# Differential Evolution optimizer
    ├── calc_H_field.m              # H-field solver (called by optimizer)
    └── eval_field_correlation.m    # Fitness function (normalized cross-correlation)
```

## Requirements

- MATLAB R2018b or later
- No additional toolboxes required

## Quick Start

### Step 1 — Prepare input files

Run the following scripts **in order**, each in its own folder:

```matlab
% 1. Generate the target field pattern
run('Generate_Target_Field/gen_target_field.m')

% 2. Generate the weighting window
run('Generate_Window_Function/gen_weight_window.m')

% 3. Generate the initial reflector surface (all-zero = standard paraboloid)
run('Generate_Reflector_Perturbation/generate_divider_surface.m')
```

This produces three data files: `Divider_C.txt`, `Divider_mag.txt`, `Divider_phase.txt`, `Weight.txt`, and `divider_surface.txt`.

### Step 2 — (Optional) Single evaluation

To verify the setup and inspect the H3-plane field before optimizing:

```matlab
% Copy divider_surface.txt and the Divider_*.txt / Weight.txt files into Single_Run/
run('Single_Run/compute_H_from_surface.m')
```

### Step 3 — Run optimization

Copy all data files into `Optimize_Reflector/`, then:

```matlab
run('Optimize_Reflector/main_optimize.m')
```

The script displays a real-time convergence curve and writes the optimized surface back to `divider_surface.txt` when finished.

### Step 4 — (Optional) Export STL

To export the optimized reflector surface for EM simulation (e.g., HFSS):

```matlab
% Copy the optimized divider_surface.txt into Generate_STL/
run('Generate_STL/generate_reflector_stl.m')
```

## Key Parameters

| Parameter | Location | Description |
|---|---|---|
| `scale` | all scripts | Geometric scaling factor (default `3.8`) |
| `freq_GHz` | `calc_H_field.m` | Operating frequency (default `380` GHz) |
| `beam_weights` | `gen_target_field.m` | Per-beam amplitude weights (default all ones) |
| `n_generations` | `main_optimize.m` | DE generations (default `300`) |
| `pop_size` | `main_optimize.m` | Population size (default `10`) |
| `Fe` | `main_optimize.m` | DE crossover weight (default `0.1`) |
| `Fm` | `main_optimize.m` | Mutation strength (default `0.1`) |

## Algorithm

The optimizer uses **Differential Evolution (DE/best/1)** with roulette-wheel parent selection:

1. Evaluate fitness of all individuals via PO field computation
2. Elite individual is carried over unchanged
3. Crossover: `child = best + Fe × (r1 − r2)`
4. Mutation: randomly perturb `mut_num` nodes by `Fm × surf_depth`
5. Repeat for `n_generations` generations

Fitness is defined as `100 × Re(k)`, where `k` is the normalized cross-correlation between the computed H3 field and the weighted target pattern.
