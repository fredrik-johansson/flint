/* Benchmark: integer relations among real numbers.

   Replicates the lattice built by _qqbar_acb_lindep (src/qqbar/acb_lindep.c):
   an n x (n+1) matrix, identity in the left block, and one column holding
   round(2^scale * x_i).  Highly unbalanced: tiny left block, huge last column.

   Planted relation: log(2) + log(3) - log(6) = 0, placed at three positions
   among n otherwise-random reals (the same setup as FLINT's p-lindep).     */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "arb.h"
#include "arf.h"
#include "mag.h"

/* Build the lindep lattice for the real vector u (length n) at precision prec. */
static void
build_lindep_lattice(fmpz_mat_t A, arb_srcptr u, slong n, slong prec)
{
    arf_t tmpr, halfr;
    mag_t max_size, tmpmag;
    fmpz_t scale_exp;
    slong i;

    arf_init(tmpr);
    arf_init(halfr);
    mag_init(max_size);
    mag_init(tmpmag);
    fmpz_init(scale_exp);

    arf_set_d(halfr, 0.5);

    for (i = 0; i < n; i++)
    {
        arf_get_mag(tmpmag, arb_midref(u + i));
        mag_max(max_size, max_size, tmpmag);
    }

    if (mag_is_zero(max_size))
        fmpz_zero(scale_exp);
    else
    {
        fmpz_neg(scale_exp, MAG_EXPREF(max_size));
        fmpz_add_ui(scale_exp, scale_exp, prec);
    }
    fmpz_sub_ui(scale_exp, scale_exp, FLINT_MAX(10, prec * 0.05));

    fmpz_mat_zero(A);
    for (i = 0; i < n; i++)
        fmpz_one(fmpz_mat_entry(A, i, i));

    for (i = 0; i < n; i++)
    {
        arf_mul_2exp_fmpz(tmpr, arb_midref(u + i), scale_exp);
        arf_add(tmpr, tmpr, halfr, prec, ARF_RND_NEAR);
        arf_floor(tmpr, tmpr);
        arf_get_fmpz(fmpz_mat_entry(A, i, n), tmpr, ARF_RND_NEAR);
    }

    arf_clear(tmpr);
    arf_clear(halfr);
    mag_clear(max_size);
    mag_clear(tmpmag);
    fmpz_clear(scale_exp);
}

/* Does row 0 encode the planted relation (+-1 at i1,i2,i3, zero elsewhere)? */
static int
check_relation(const fmpz_mat_t A, slong n, slong i1, slong i2, slong i3)
{
    slong i;
    int ok;
    const fmpz *r = fmpz_mat_entry(A, 0, 0);

    ok = (fmpz_is_one(r + i1) && fmpz_is_one(r + i2) && fmpz_equal_si(r + i3, -1))
      || (fmpz_equal_si(r + i1, -1) && fmpz_equal_si(r + i2, -1) && fmpz_is_one(r + i3));

    for (i = 0; i < n; i++)
        if (i != i1 && i != i2 && i != i3)
            ok = ok && fmpz_is_zero(r + i);

    return ok;
}

static double
wall(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

int
main(int argc, char **argv)
{
    slong n, prec, i, i1, i2, i3;
    slong nmax = 120, K = 16;
    arb_ptr u;
    flint_rand_t state;
    fmpz_lll_t fl;

    if (argc > 1)
        nmax = atol(argv[1]);
    if (argc > 2)
        K = atol(argv[2]);

    flint_rand_init(state);
    /* qqbar uses delta = 0.75, eta = 0.51 for lindep */
    fmpz_lll_context_init(fl, 0.75, 0.51, 1, 0);

    u = _arb_vec_init(nmax);

    flint_printf("Integer relations: log(2) + log(3) - log(6) = 0 among n reals\n");
    flint_printf("lattice is n x (n+1), last column ~prec bits\n\n");
    flint_printf("%-5s %-7s %-13s %-13s %-8s %-7s %-7s\n",
                 "n", "prec", "fmpz_lll(s)", "blaster(s)", "ratio", "lll_ok", "bl_ok");

    for (n = 20; n <= nmax; n += 20)
    {
        fmpz_mat_t A0, A1, A2;
        double t0, t_lll, t_bl;
        int ok1, ok2, r1, r2;

        for (i = 0; i < n; i++)
            arb_urandom(u + i, state, 4096);

        i1 = n / 4;
        i2 = n / 2;
        i3 = (3 * n) / 4;

        /* Precision must grow with n or spurious short vectors beat the
           planted relation.  Use prec = K*n (K from argv[2], default 16). */
        prec = FLINT_MAX(64, K * n);

        arb_log_ui(u + i1, 2, prec);
        arb_log_ui(u + i2, 3, prec);
        arb_log_ui(u + i3, 6, prec);
        arb_mul_2exp_si(u + i1, u + i1, -1);
        arb_mul_2exp_si(u + i2, u + i2, -1);
        arb_mul_2exp_si(u + i3, u + i3, -1);

        fmpz_mat_init(A0, n, n + 1);
        build_lindep_lattice(A0, u, n, prec);

        /* Time both algorithms on identical copies of the same lattice. */
        fmpz_mat_init(A1, n, n + 1);
        fmpz_mat_set(A1, A0);
        t0 = wall();
        fmpz_lll(A1, NULL, fl);
        t_lll = wall() - t0;
        ok1 = check_relation(A1, n, i1, i2, i3);

        fmpz_mat_init(A2, n, n + 1);
        fmpz_mat_set(A2, A0);
        t0 = wall();
        r2 = fmpz_lll_blaster(A2, NULL, fl);
        t_bl = wall() - t0;
        ok2 = check_relation(A2, n, i1, i2, i3);
        (void) r2; (void) r1;

        flint_printf("%-5wd %-7wd %-13.4f %-13.4f %-8.2f %-7s %-7s\n",
                     n, prec, t_lll, t_bl,
                     (t_bl > 0) ? t_lll / t_bl : 0.0,
                     ok1 ? "yes" : "NO", ok2 ? "yes" : "NO");
        fflush(stdout);

        fmpz_mat_clear(A0);
        fmpz_mat_clear(A1);
        fmpz_mat_clear(A2);
    }

    _arb_vec_clear(u, nmax);
    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
