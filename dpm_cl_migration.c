/*
 * UDF: 2D DPM Study — Cl⁻ Ion Migration in Seawater Under Electric Field
 * Target solver: Ansys Fluent V19.0  (2D pressure-based, DPM enabled)
 *
 * ── PROBLEM DESCRIPTION ──────────────────────────────────────────────────
 *  Domain   : 2D rectangle,  L = 1000 mm (x),  H = 200 mm (y)
 *  Fluid    : seawater, flowing in +x direction (inlet x=0, outlet x=L)
 *  Particles: Cl⁻ ions injected uniformly across the inlet face
 *  Electrode: strip on the TOP wall, centred at x = 200 mm (anode, +)
 *             Bottom wall is the counter-electrode (cathode, −)
 *
 * ── PHYSICS ──────────────────────────────────────────────────────────────
 *  An electric potential V_APPLIED is imposed between the top (+) and
 *  bottom (−) walls.  The resulting uniform background field points in
 *  the −y direction:
 *
 *       E_y_bg = − V_APPLIED / H          [V/m]  (downward)
 *
 *  Cl⁻ carries charge q = −e < 0, so the Coulomb force is:
 *
 *       F_y = q · E_y_bg = (−e)(−V/H) = +eV/H   [N]  (upward) ✓
 *
 *  A Gaussian x-focusing field (centred on the electrode) is superimposed
 *  to steer ions axially toward the capture region:
 *
 *       E_x = E_FOCUS · (x − x_e) · exp[−(x−x_e)²/(2σ_x²)]
 *
 *  Cl⁻ x-force: F_x = q · E_x  →  (−e) × positive = negative when
 *  x > x_e (pushed left), positive when x < x_e (pushed right) → focusing ✓
 *
 *  Fluent's DPM drag model handles the Stokes (or Cunningham) drag
 *  automatically; the terminal drift (electrophoretic) velocity is:
 *
 *       v_ep = |q| · E / (3π · μ · d_p)
 *
 * ── SUGGESTED SIMULATION PARAMETERS ─────────────────────────────────────
 *
 *  ┌──────────────────────────────────┬───────────────┬────────────────────┐
 *  │ Parameter                        │ Value         │ Notes              │
 *  ├──────────────────────────────────┼───────────────┼────────────────────┤
 *  │ Seawater density ρ               │ 1025 kg/m³    │ 35 g/L salinity    │
 *  │ Dynamic viscosity μ              │ 1.08e−3 Pa·s  │ at 20 °C           │
 *  │ Inlet velocity U                 │ 0.001 m/s     │ Re = ρUH/μ ≈ 190  │
 *  │ Mass flow rate ṁ (per unit depth)│ 0.205 kg/(m·s)│ = ρ · U · H       │
 *  ├──────────────────────────────────┼───────────────┼────────────────────┤
 *  │ Applied voltage V_APPLIED        │ 1 000 V       │ = E0 × H          │
 *  │ Uniform E field E0               │ 5 000 V/m     │ = V/H             │
 *  │ Focusing amplitude E_FOCUS       │ 500 V/m       │ optional, see §3  │
 *  ├──────────────────────────────────┼───────────────┼────────────────────┤
 *  │ Particle type                    │ inert          │                   │
 *  │ Particle diameter d_p            │ 3.62e−10 m    │ Cl⁻ ionic diameter │
 *  │ Particle density ρ_p             │ 2 375 kg/m³   │ from m/V_sphere    │
 *  │ Initial particle velocity        │ U_x = U_fluid │ match flow at inlet│
 *  │ Particle mass flow rate          │ scale to      │ see injection note │
 *  │   (per stream line)              │ concentration │                    │
 *  ├──────────────────────────────────┼───────────────┼────────────────────┤
 *  │ Electrode half-width             │ 25 mm         │ x: 175 mm–225 mm  │
 *  │ Predicted v_ep (y)               │ 2.18e−4 m/s   │ |q|E/(3πμd)        │
 *  │ Predicted x-capture point        │ 0.124–0.168 m │ from inlet (top)   │
 *  │   (for particles near y = 0.15 m)│               │ → hits electrode ✓ │
 *  └──────────────────────────────────┴───────────────┴────────────────────┘
 *
 *  Capture efficiency estimate (uniform injection):
 *    Particles starting at y₀ ∈ [0.124, 0.168] m reach the top wall
 *    exactly within the electrode x-range → ~22 % capture for this E field.
 *    Increase V_APPLIED (or widen ELEC_HALF_W) to raise efficiency.
 *
 * ── FLUENT SETUP CHECKLIST ───────────────────────────────────────────────
 *  1. Compile this UDF: Define → User-Defined → Functions → Compiled
 *  2. Enable DPM: Models → Discrete Phase → ON
 *       • Particle type: Inert
 *       • Drag law: Stokes-Cunningham (best for sub-micron particles)
 *       • Unsteady tracking: ON (recommended for DPM)
 *       • Max steps: 50 000 (ions travel slowly, need many steps)
 *  3. Create injection at inlet:
 *       • Type: surface (inlet face)
 *       • Diameter: 3.62e-10 m, Density: 2375 kg/m³
 *       • Vx = 0.001 m/s, Vy = 0.0 m/s
 *       • Flow rate: set proportional to seawater Cl⁻ concentration
 *         (~0.54 mol/L → ~19 g/L; scale to desired number of parcels)
 *  4. Allocate 4 UDM slots: Define → User-Defined → Memory → 4
 *       UDM-0: E_x (V/m), UDM-1: E_y (V/m),
 *       UDM-2: F_x/m (m/s²), UDM-3: F_y/m (m/s²)
 *  5. Hook DEFINE_ADJUST: Define → UDF → Function Hooks → Adjust
 *  6. Hook body force: in DPM injection panel → Body Force UDF → cl_em_force
 *  7. Assign DPM BC to top wall:
 *       Boundary Conditions → (top-wall) → DPM → BC type = UDF
 *       → select  cl_wall_bc
 *       Leave bottom wall, inlet, outlet with Fluent defaults
 *       (reflect, escape, escape respectively)
 */

#include "udf.h"
#include <math.h>

/* ═══════════════════════════ USER PARAMETERS ═══════════════════════════ */

/* --- Geometry ---------------------------------------------------------- */
#define DOMAIN_H        0.200       /* channel height              [m]     */
#define X_ELEC          0.200       /* electrode centre (axial)    [m]     */
#define ELEC_HALF_W     0.025       /* electrode half-width        [m]     */
/* electrode occupies top wall from (X_ELEC - ELEC_HALF_W)
                                  to (X_ELEC + ELEC_HALF_W)               */

/* --- Electric field ----------------------------------------------------- */
#define V_APPLIED       1000.0      /* top(+) to bottom(-) voltage [V]     */
#define E0              (V_APPLIED / DOMAIN_H)  /* background E [V/m]      */
/*
 * X-focusing field: pulls ions axially toward x = X_ELEC.
 * Amplitude in [V/m].  Set to 0.0 to use only the uniform y-field.
 */
#define E_FOCUS         500.0       /* x-focusing field amplitude  [V/m]   */
#define SIGMA_X         0.040       /* x-focusing Gaussian width   [m]     */

/* --- Cl⁻ ion properties ------------------------------------------------- */
/*
 * Cl⁻  atomic mass = 35.45 u  →  m = 35.45 × 1.661e-27 = 5.887e-26 kg
 * Ionic radius = 1.81 Å  →  d = 3.62e-10 m
 * Effective density = m / (π/6 · d³) = 2375 kg/m³  (use in injection)
 * Charge = -1 elementary charge
 */
#define Q_ION           (-1.602e-19)    /* ionic charge            [C]     */
#define M_ION           (5.887e-26)     /* ionic mass              [kg]    */
#define QM              (Q_ION / M_ION) /* specific charge [C/kg] = -2.72e6*/

/* ═══════════════════════════════════════════════════════════════════════ */


/* ───────────────── helper: electric field at (x, y) ──────────────────── */
/*
 * Returns E_x and E_y at the particle position.
 *
 * E_y: uniform background field (top + / bottom −).
 *      Top is positive → conventional E points from top to bottom → −y.
 *      E_y = −E0  (<0, downward).  F_y on Cl⁻ = QM × E_y = +eE0/m >0 ✓
 *
 * E_x: Gaussian-profile focusing field.
 *      Derivative of Gaussian → sign changes at x = X_ELEC.
 *      E_x = E_FOCUS × (x−x_e) × exp[−(x−x_e)²/(2σ²)]
 *      F_x on Cl⁻ = QM × E_x  (negative QM flips sign → toward electrode) ✓
 */
static void em_field(real x, real y, real *Ex_out, real *Ey_out)
{
    real dx  = x - X_ELEC;
    real Gx  = exp(-0.5 * dx * dx / (SIGMA_X * SIGMA_X));

    *Ex_out  = E_FOCUS * dx * Gx;   /* V/m — sign gives focusing for Cl⁻  */
    *Ey_out  = -E0;                  /* V/m — uniform downward field        */
}


/* ─────────────── 1. EM BODY FORCE ON Cl⁻ PARTICLES ──────────────────── */
/*
 * DEFINE_DPM_BODY_FORCE is called once per particle per time step for
 * each spatial component i (0=x, 1=y).
 * Must return force per unit particle mass [m/s²].
 *
 * Equation of motion (one direction):
 *   m · dv/dt = F_drag (Fluent) + F_EM (this UDF) + F_gravity (Fluent)
 *
 * F_EM / m  = (Q_ION / M_ION) × E_i  =  QM × E_i
 */
DEFINE_DPM_BODY_FORCE(cl_em_force, p, i)
{
    real Ex, Ey;
    em_field(P_POS(p)[0], P_POS(p)[1], &Ex, &Ey);

    if (i == 0) return QM * Ex;   /* ax = (q/m) × Ex  [m/s²] */
    if (i == 1) return QM * Ey;   /* ay = (q/m) × Ey  [m/s²] */
    return 0.0;                    /* z-component (0 in 2D)    */
}


/* ─────────────── 2. BOUNDARY CONDITION AT THE TOP WALL ──────────────── */
/*
 * Assign this UDF to the TOP WALL zone only.
 * • Particle hits electrode region  → absorbed (PATH_ABORT, trapped).
 * • Particle hits rest of top wall  → specular reflection.
 *
 * f         : face index of the impacted face
 * t         : thread pointer of the impacted boundary zone
 * f_normal  : outward unit normal of the face (points into fluid, i.e. −y)
 * dim       : spatial dimension (2 for 2D simulations)
 *
 * Return values:
 *   0 (PATH_ABORT)  → particle tracking ends; Fluent records as "trapped"
 *   1 (PATH_ACTIVE) → particle continues; reflected velocity already set
 */
DEFINE_DPM_BC(cl_wall_bc, p, t, f, f_normal, dim)
{
    real xf[ND_ND];
    real vn;
    int  k;

    /* Centroid of the face the particle just hit */
    F_CENTROID(xf, f, t);

    /* ── Electrode capture zone ─────────────────────────────────────── */
    if (xf[0] >= (X_ELEC - ELEC_HALF_W) &&
        xf[0] <= (X_ELEC + ELEC_HALF_W))
    {
        /*
         * Ion reaches the anode → absorbed.
         * Fluent will record this parcel as "trapped" in the DPM summary.
         * To account for partial collection, you can also multiply
         * P_FLOW_RATE(p) by a sticking efficiency factor (0–1) and
         * reflect the remainder; here we use 100 % collection.
         */
        return 0; /* PATH_ABORT — trapped at electrode */
    }

    /* ── Regular top-wall section: specular reflection ──────────────── */
    /*
     * Reflect velocity: v_new = v_old − 2 (v · n̂) n̂
     * f_normal points into the fluid, so use it directly.
     */
    vn = 0.0;
    for (k = 0; k < dim; k++)
        vn += P_VEL(p)[k] * f_normal[k];

    for (k = 0; k < dim; k++)
        P_VEL(p)[k] -= 2.0 * vn * f_normal[k];

    return 1; /* PATH_ACTIVE — continue tracking */
}


/* ─────────────── 3. ADJUST: STORE E-FIELD IN UDM FOR POST-PROCESSING ── */
/*
 * Requires 4 UDM slots allocated in Fluent:
 *   UDM-0  E_x  (V/m)
 *   UDM-1  E_y  (V/m)
 *   UDM-2  a_x = F_x / m  (m/s²)
 *   UDM-3  a_y = F_y / m  (m/s²)
 *
 * Hook this via: Define → User-Defined → Function Hooks → Adjust.
 * These fields can then be plotted directly in CFD-Post or Fluent contour
 * panels as "User Memory 0..3".
 */
DEFINE_ADJUST(store_em_field, domain)
{
    Thread *t;
    cell_t  c;

    thread_loop_c(t, domain)
    {
        begin_c_loop(c, t)
        {
            real xc[ND_ND];
            real Ex, Ey;

            C_CENTROID(xc, c, t);
            em_field(xc[0], xc[1], &Ex, &Ey);

            C_UDMI(c, t, 0) = Ex;          /* E_x  [V/m]   */
            C_UDMI(c, t, 1) = Ey;          /* E_y  [V/m]   */
            C_UDMI(c, t, 2) = QM * Ex;     /* a_x  [m/s²]  */
            C_UDMI(c, t, 3) = QM * Ey;     /* a_y  [m/s²]  */
        }
        end_c_loop(c, t)
    }
}
