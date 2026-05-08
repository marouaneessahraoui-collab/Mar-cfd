/*
 * UDF: External Magnetic Field on Seawater Pipe Flow (MHD Lorentz Force)
 *
 * Physics:
 *   Seawater (electrically conducting) flows through a 1-meter pipe.
 *   An external magnetic field B is applied transversely (y-direction by
 *   default) and localised around x = 0.5 m via a Gaussian envelope.
 *
 *   The induced Lorentz body force density is:
 *       J  = sigma * (V x B)          [A/m^2]  (low-Rm approximation)
 *       F  = J  x B                   [N/m^3]
 *
 *   Source terms are injected into the three momentum equations.
 *
 * Coordinate conventions:
 *   x  – axial (flow) direction,  0 … 1 m
 *   y  – transverse (direction of B by default)
 *   z  – remaining transverse direction
 *
 * Parameters (edit the #define block below):
 *   B0        – peak magnetic flux density [T]
 *   B_DIR     – field direction: 0=x, 1=y, 2=z
 *   X_CENTER  – axial centre of the magnet [m]
 *   SIGMA_MAG – Gaussian half-width of the field region [m]
 *   SIGMA_EL  – electrical conductivity of seawater [S/m]
 */

#include "udf.h"
#include <math.h>

/* ──────────────────────────── USER PARAMETERS ──────────────────────────── */
#define B0         0.5      /* peak magnetic flux density [T]                 */
#define B_DIR      1        /* field direction: 0=x, 1=y, 2=z                */
#define X_CENTER   0.5      /* axial position of magnet centre [m]            */
#define SIGMA_MAG  0.05     /* Gaussian width of field region [m]             */
#define SIGMA_EL   4.0      /* seawater electrical conductivity [S/m]        */
/* ─────────────────────────────────────────────────────────────────────────── */

/*
 * gaussian_B:
 *   Returns local magnetic flux density magnitude based on axial position x.
 *   B(x) = B0 * exp( -(x - x_c)^2 / (2*sigma^2) )
 */
static real gaussian_B(real x)
{
    real dx = x - X_CENTER;
    return B0 * exp(-0.5 * (dx * dx) / (SIGMA_MAG * SIGMA_MAG));
}

/*
 * lorentz_force:
 *   Computes the Lorentz force density vector F = sigma*(VxB)xB
 *   and stores the three components in f[0..2].
 *
 *   With B = B_mag * e_dir, the closed-form result for each direction is:
 *
 *   B in y (B_DIR=1):
 *     Fx = -sigma * B^2 * vx
 *     Fy =  0
 *     Fz = -sigma * B^2 * vz
 *
 *   B in z (B_DIR=2):
 *     Fx = -sigma * B^2 * vx
 *     Fy = -sigma * B^2 * vy
 *     Fz =  0
 *
 *   B in x (B_DIR=0):
 *     Fx =  0
 *     Fy = -sigma * B^2 * vy
 *     Fz = -sigma * B^2 * vz
 */
static void lorentz_force(real vx, real vy, real vz,
                          real B_mag, real f[3])
{
    real sB2 = SIGMA_EL * B_mag * B_mag;

#if   B_DIR == 1   /* B along y */
    f[0] = -sB2 * vx;
    f[1] =  0.0;
    f[2] = -sB2 * vz;

#elif B_DIR == 2   /* B along z */
    f[0] = -sB2 * vx;
    f[1] = -sB2 * vy;
    f[2] =  0.0;

#else              /* B along x (B_DIR == 0) */
    f[0] =  0.0;
    f[1] = -sB2 * vy;
    f[2] = -sB2 * vz;
#endif
}

/* ──────────────────── MOMENTUM SOURCE TERMS ─────────────────────────────── */

DEFINE_SOURCE(mhd_x_mom, cell, thread, dS, eqn)
{
    real x[ND_ND];
    C_CENTROID(x, cell, thread);

    real B   = gaussian_B(x[0]);
    real vx  = C_U(cell, thread);
    real vy  = C_V(cell, thread);
    real vz  = C_W(cell, thread);
    real f[3];

    lorentz_force(vx, vy, vz, B, f);

    /* Linearise: dF/du = -sigma*B^2 (only when B_DIR != 0) */
#if B_DIR == 0
    dS[eqn] = 0.0;
#else
    dS[eqn] = -SIGMA_EL * B * B;
#endif

    return f[0];
}

DEFINE_SOURCE(mhd_y_mom, cell, thread, dS, eqn)
{
    real x[ND_ND];
    C_CENTROID(x, cell, thread);

    real B   = gaussian_B(x[0]);
    real vx  = C_U(cell, thread);
    real vy  = C_V(cell, thread);
    real vz  = C_W(cell, thread);
    real f[3];

    lorentz_force(vx, vy, vz, B, f);

#if B_DIR == 1
    dS[eqn] = 0.0;
#else
    dS[eqn] = -SIGMA_EL * B * B;
#endif

    return f[1];
}

DEFINE_SOURCE(mhd_z_mom, cell, thread, dS, eqn)
{
    real x[ND_ND];
    C_CENTROID(x, cell, thread);

    real B   = gaussian_B(x[0]);
    real vx  = C_U(cell, thread);
    real vy  = C_V(cell, thread);
    real vz  = C_W(cell, thread);
    real f[3];

    lorentz_force(vx, vy, vz, B, f);

#if B_DIR == 2
    dS[eqn] = 0.0;
#else
    dS[eqn] = -SIGMA_EL * B * B;
#endif

    return f[2];
}

/* ──────────────────── OPTIONAL: UDM FIELD STORAGE ───────────────────────
 *
 * If you allocate 4 User-Defined Memory slots (UDM 0..3) in Fluent
 * (Define → User-Defined → Memory), this macro will write:
 *   UDM-0  : local |B| [T]
 *   UDM-1  : Fx    [N/m^3]
 *   UDM-2  : Fy    [N/m^3]
 *   UDM-3  : Fz    [N/m^3]
 * at every flow iteration so the fields are visible in post-processing.
 * ────────────────────────────────────────────────────────────────────────── */

DEFINE_ADJUST(mhd_store_udm, domain)
{
    Thread *t;
    cell_t  c;

    thread_loop_c(t, domain)
    {
        begin_c_loop(c, t)
        {
            real x[ND_ND];
            C_CENTROID(x, c, t);

            real B  = gaussian_B(x[0]);
            real vx = C_U(c, t);
            real vy = C_V(c, t);
            real vz = C_W(c, t);
            real f[3];

            lorentz_force(vx, vy, vz, B, f);

            C_UDMI(c, t, 0) = B;
            C_UDMI(c, t, 1) = f[0];
            C_UDMI(c, t, 2) = f[1];
            C_UDMI(c, t, 3) = f[2];
        }
        end_c_loop(c, t)
    }
}
