# Magnetic Field on Seawater Pipe Flow — Ansys Fluent UDF

## Problem description

Seawater (electrically conducting) flows axially through a straight pipe of
length **L = 1 m**.  An external permanent magnet or electromagnet is placed
around the pipe at **x = 0.5 m** and imposes a transverse magnetic field.

The interaction of the magnetic field with the moving conducting fluid
produces a Lorentz body force that opposes the flow (Hartmann braking) and
can redirect it in the cross-sectional plane.

---

## Physical model (low magnetic Reynolds number)

| Quantity | Symbol | Value used |
|---|---|---|
| Seawater electrical conductivity | σ | 4.0 S/m |
| Peak magnetic flux density | B₀ | 0.5 T |
| Field direction (transverse) | y | (configurable) |
| Magnet centre | x_c | 0.5 m |
| Gaussian half-width | σ_mag | 0.05 m |

### Field profile
```
B(x) = B0 · exp[ −(x − x_c)² / (2 σ_mag²) ]
```

### Induced current density (low-Rm approximation)
```
J = σ · (V × B)
```
The electric potential gradient ∇φ is neglected, which is valid when the
Hartmann number Ha is moderate and the walls are insulating.

### Lorentz force density (source term injected into Navier–Stokes)
```
F = J × B = σ · (V × B) × B   [N/m³]
```
For **B ∥ ŷ** this reduces to:
```
Fx = −σ B² vx
Fy =  0
Fz = −σ B² vz
```

---

## Files

| File | Purpose |
|---|---|
| `magnetic_field_seawater.c` | Fluent UDF source |

---

## Configurable parameters (`#define` block in the UDF)

| Macro | Default | Description |
|---|---|---|
| `B0` | 0.5 T | Peak flux density |
| `B_DIR` | 1 | Field direction: 0=x, 1=y, 2=z |
| `X_CENTER` | 0.5 m | Axial centre of the magnet |
| `SIGMA_MAG` | 0.05 m | Gaussian width of the field zone |
| `SIGMA_EL` | 4.0 S/m | Seawater electrical conductivity |

---

## Step-by-step Fluent setup

### 1. Compile the UDF
**Define → User-Defined → Functions → Compiled…**
- Add `magnetic_field_seawater.c` to the source list
- Click **Build**, then **Load**

### 2. Allocate UDM storage (optional — for post-processing)
**Define → User-Defined → Memory…**  →  set to **4** slots

Slot | Stored quantity
-----|----------------
0 | Local \|B\| (T)
1 | Fx (N/m³)
2 | Fy (N/m³)
3 | Fz (N/m³)

Hook `mhd_store_udm` under **Define → User-Defined → Function Hooks →
Execute at End of Time Step / Iteration**.

### 3. Assign source terms to the fluid zone
**Cell Zone Conditions → (fluid zone) → Source Terms**

Equation | UDF to select
---------|---------------
X Momentum | `mhd_x_mom`
Y Momentum | `mhd_y_mom`
Z Momentum | `mhd_z_mom`

### 4. Seawater material properties
Go to **Materials → Fluid** and set (or create a new material "seawater"):

| Property | Value |
|---|---|
| Density | 1025 kg/m³ |
| Dynamic viscosity | 1.08 × 10⁻³ Pa·s |
| Thermal conductivity | 0.6 W/(m·K) |
| Electrical conductivity | *handled inside UDF via `SIGMA_EL`* |

> Fluent's MHD add-on module is **not** required here; the Lorentz force
> is injected directly as a volumetric source term.

### 5. Boundary conditions (example)
- **Inlet** (x = 0): velocity-inlet, u = 0.1–1.0 m/s, v = w = 0
- **Outlet** (x = 1 m): pressure-outlet, gauge pressure = 0
- **Pipe wall**: no-slip

### 6. Solution methods (recommended)
- Solver: pressure-based, steady
- Scheme: SIMPLE or SIMPLEC
- Spatial discretisation: Second-Order Upwind for momentum
- Convergence: residuals < 1 × 10⁻⁵

---

## Expected physics

| Hartmann number | `Ha = B0·R·√(σ/μ)` |
|---|---|
| Typical value (R=0.01 m) | Ha ≈ 0.5·0.01·√(4/1.08×10⁻³) ≈ **9.6** |

At Ha ≈ 10 the flow profile flattens noticeably in the magnetic field region
(Hartmann layers form on the walls perpendicular to B) and the axial pressure
drop increases compared to the purely hydrodynamic case.

---

## Extending the model

- **Non-uniform B field in the cross-section**: add radial or angular
  dependence inside `gaussian_B()` using `x[1]` and `x[2]`.
- **Electric potential equation**: add a scalar transport equation for φ
  and include the `∇φ` term in the Lorentz force for higher accuracy at
  large Ha.
- **Joule heating**: add an energy source term `Q = J·J/σ` via an
  additional `DEFINE_SOURCE` hooked to the energy equation.
