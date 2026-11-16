// simulate_hr.cpp
//
// Purpose
//   Compiled simulation engine of the R package hrsim. It integrates a model
//   of heteroresistance with compensatory mutations under serial transfer:
//   twenty bacterial states (R0..R9 carry 0 to 9 copies of a beta-lactamase
//   gene; R0C..R9C are the same with a compensatory mutation), a limiting
//   resource, an antibiotic and beta-lactamase. Three functions are exported
//   to R: simulate_hr() runs one replicate; debug_step1() and
//   debug_step1_fired() expose the first step so that the tests can compare
//   it with the R reference.
//
// Paper
//   Berryhill, B. A. and Andersson, D. I. (2026). Unstable Heteroresistance
//   as a Stepping Stone to Stable Resistance via Compensatory Mutations.
//   bioRxiv. doi:10.64898/2026.10.08.756632
//   Supplementary Equations (S. Eq), Supplementary Tables and the
//   Supplementary Text cited below are those of this paper.
//
// Numerical scheme
//   Time is in hours. A run takes n_steps = round(horizon_h / dt) steps of
//   fixed length dt. The scheme has three components.
//   Explicit Euler: the rates of S. Eq 1 to 5 are evaluated from the state at
//   the start of the step and applied with step dt.
//   Monte Carlo transition layer: every transition between states
//   (recombination between copy numbers, or compensation) fires with a
//   probability that is also evaluated from the start-of-step state. A
//   transition that fires moves one cell per mL from its source state to its
//   target state. These are the M terms of S. Eq 6 to 25; the procedure is
//   introduced in the Supplementary Text, "Stochastic Transitions Between
//   States".
//   Serial transfer: a clock advances by dt each step. When it reaches
//   transfer_interval_h the state is recorded, the populations and the
//   antibiotic are diluted, a fresh antibiotic dose is added, the resource is
//   restored and beta-lactamase is set to zero.
//   Within a step the order is: evaluate all rates and probabilities; draw
//   which transitions fire; apply the Euler increments and the fired
//   transitions together and floor A and the populations at zero; advance
//   the transfer clock and, if a transfer is due, record a row and carry out
//   the transfer; record a grid row if one is due.
//
// Equations
//   i runs over the 20 states; sums are over all 20.
//     psi     = r / (r + K)                                           S. Eq 2
//     Pi_i    = (v_i - (v_i - vmin_i) * (A / MIC_i)^kappa_i
//                      / ((A / MIC_i)^kappa_i - vmin_i / v_i)) * psi  S. Eq 1
//     dr/dt   = -e * psi * sum(v_i * R_i)                             S. Eq 3
//     dA/dt   = -d * A - KB * B * A                                   S. Eq 4
//     dB/dt   =  e * psi * sum(p_i * R_i)                             S. Eq 5
//     dR_i/dt =  R_i * Pi_i + transitions in - transitions out        S. Eq 6 to 25
//   One step, with every right-hand side taken at the start of the step:
//     R_i <- max(0, R_i + dt * R_i * Pi_i + net_i)
//     r   <- r + dt * dr/dt
//     A   <- max(0, A + dt * dA/dt)
//     B   <- B + dt * dB/dt
//   net_i is the number of transitions that fired into state i in the step
//   minus the number that fired out of it. Transition k, from state f to
//   state t with rate constant rate_k, fires with probability
//     P_k = rate_k * dt * R_f * psi
//   (never when P_k <= 0, always when P_k >= 1).
//   Serial transfer:
//     R_i <- dil * R_i;   A <- dil * A + AMAX;   r <- C;   B <- 0
//
// Symbols
//   Units of the states and parameters that appear in Supplementary Tables
//   1 and 2 are those listed there.
//
//   Symbol      Meaning                           Unit               Code variable     S. Eq
//   R_i         density of state i                cells per mL       pop[i]            6 to 25
//   r           limiting resource                 ug per mL          r                 3
//   A           antibiotic                        ug per mL          A                 4
//   B           beta-lactamase                    units per mL       B                 5
//   psi         resource limitation term          dimensionless      Snapshot::psi     2
//   Pi_i        net growth or death rate          per cell per hour  Q                 1
//   v_i         maximum growth rate               per cell per hour  Vv[i]  (v)        1, 3
//   vmin_i      minimum growth or death rate      per cell per hour  Vm[i]  (vm)       1
//   kappa_i     Hill coefficient                  dimensionless      Kh[i]  (kH)       1
//   MIC_i       minimum inhibitory concentration  ug per mL          Mm[i]  (M)        1
//   p_i         beta-lactamase production rate    units per ug       Pp[i]  (Pb)       5
//                                                 per hour
//   K           Monod constant                    ug per mL          kmon              2
//   e           conversion efficiency             ug per cell        eR                3, 5
//   d           antibiotic decay rate             per hour           d0                4
//   KB          beta-lactamase association        mL per unit        KB                4
//               constant                          per hour
//   mu, gamma,  recombination and compensation    per cell per hour  rate[k]           6 to 25
//   omega       rates; rate[k] is one of them                        (rate_vec)        (M terms)
//               times the transition's
//               multiplicity
//   dt          Euler step                        hours              dt                -
//   AMAX        antibiotic dose at each transfer  ug per mL          AMAX              -
//   dil         dilution factor at each transfer  dimensionless      dil               -
//   C           resource after each transfer      ug per mL          Cres              -
//
//   The name in brackets is the argument of the exported functions that
//   supplies the value.
//
//   K is written k in Supplementary Table 2. In the equations above, K is
//   added to r, in ug per mL, and d * A and KB * B * A are terms of dA/dt,
//   so that d is per hour and, with B in units per mL, KB is in mL per unit
//   per hour.
//
// What this file does not do
//   It does not run replicates in parallel: one call runs one replicate on
//   one thread. It makes no draw from R's random number generator; the
//   generator defined below is its only source of random numbers. It does
//   not classify outcomes, and it does not check the lengths of its input
//   vectors or the range of the transition indices (see the contract of
//   simulate_hr()).
//
// Floating point
//   The package is compiled with -ffp-contract=off (src/Makevars), which
//   stops the compiler contracting a multiplication and an addition into one
//   fused multiply-add, which rounds once where the two operations performed
//   separately round twice. The agreement between the two engines reported
//   in the Methods of the paper was obtained with this flag.
//
// Provenance. This file was translated from the R reference implementation
// (`R/simulation/03_engine_R.R`) to C++ with the assistance of a generative
// AI model (Claude Opus 4.8, Anthropic), as described in the Methods. The
// translation was verified against the R reference before any results were
// generated from it: deterministically, with the random layer disabled, and
// statistically, on outcome frequencies. The checks are in `tests/`. The
// authors reviewed the code and take responsibility for it.

#include <Rcpp.h>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>

static const int N_POP  = 20;   // cell states, in the fixed order R0..R9, R0C..R9C
static const int N_COLS = 25;   // output columns: time, 20 states, r, A, B, cycle_end

// ---- Power function for the Hill term ---------------------------------------
//
// Computes x^e for the term (A / MIC)^kappa of S. Eq 1. x^1 is x exactly, so
// for an exponent of 1 the base is returned as it is: no rounding is
// introduced and no library call is made. Every state has kappa = 1 in the
// paper's parameter set (Supplementary Table 2, Hill Coefficient), so this is
// the branch taken in each of the 20 evaluations per step. Exponents 2 and 0
// are also evaluated directly (x * x and 1); any other exponent goes to
// std::pow.
static inline double pow_fast(double x, double e) {
    if (e == 1.0) return x;
    if (e == 2.0) return x * x;
    if (e == 0.0) return 1.0;
    return std::pow(x, e);
}

// ---- Random number generator ------------------------------------------------
//
// Algorithm: xoshiro256++ (Blackman and Vigna, https://prng.di.unimi.it/).
// Seeding: the four 64-bit state words are four successive outputs of
// SplitMix64 started from the run's seed.
// Uniform on [0, 1): the top 53 bits of a 64-bit output, multiplied by 2^-53.
// Draw order: within a step, one draw for each eligible transition (one whose
// probability lies strictly between 0 and 1), in transition-table order.
// The generator is independent of R's: no draw is taken from R's generator,
// and set.seed() has no effect on the stream. It is a simulation generator
// and is not cryptographically secure.
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

struct SplitMix64 {
    uint64_t x;
    uint64_t next() {
        uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
};
struct Xoshiro256pp {
    uint64_t s[4];
    void seed(uint64_t seed64) { SplitMix64 sm{seed64}; for (int i = 0; i < 4; ++i) s[i] = sm.next(); }
    uint64_t next() {
        const uint64_t r = rotl(s[0] + s[3], 23) + s[0];
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t;    s[3] = rotl(s[3], 45);
        return r;
    }
    double unif() { return (next() >> 11) * 0x1.0p-53; }
};

// ---- Rates and probabilities of one step ------------------------------------
//
// Snapshot holds everything a step needs from the state at its start.
// simulate_hr() and both debug functions fill it through the single function
// snapshot_kernel(), so the tests exercise the same arithmetic as a full run.
struct Snapshot {
    double psi = 0.0, dr = 0.0, dB = 0.0, dA = 0.0;   // psi (S. Eq 2) and the time derivatives of r, B and A (S. Eq 3, 5, 4)
    std::vector<double> growth;   // Euler increment dt * R_i * Pi_i of each state, cells per mL; length N_POP
    std::vector<double> probs;    // probability that each transition fires in this step; length NTR
};

// snapshot_kernel: evaluates S. Eq 1 to 5 and the transition probabilities
// for the state it is given. It reads the state and writes only to s.
//
// Arguments
//   pop          densities of the 20 states, cells per mL; length N_POP
//   r, A, B      resource (ug per mL), antibiotic (ug per mL) and
//                beta-lactamase (units per mL)
//   Vv, Vm, Kh,  per-state v, vmin, kappa, MIC and p of the symbol table;
//   Mm, Pp       length N_POP each
//   f            source state of each transition, 0-based; length NTR
//   rate         rate constant of each transition, per cell per hour;
//                length NTR, which is taken from this vector
//   dt           Euler step, hours
//   kmon, eR,    K, e, d and KB of the symbol table
//   d0, KB
//   s            output; s.growth must already have length N_POP and s.probs
//                length NTR
//
// Preconditions are not checked here. The three exported functions size the
// vectors and convert the transition indices to 0-based before they call.
static inline void snapshot_kernel(
    const std::vector<double>& pop, double r, double A, double B,
    const std::vector<double>& Vv, const std::vector<double>& Vm,
    const std::vector<double>& Kh, const std::vector<double>& Mm,
    const std::vector<double>& Pp,
    const std::vector<int>& f, const std::vector<double>& rate,
    double dt, double kmon, double eR, double d0, double KB,
    Snapshot& s)
{
    s.psi = r / (r + kmon);                       // S. Eq 2. r is not floored by the update, and no guard is applied here.
    for (int i = 0; i < N_POP; ++i) {
        const double AM = pow_fast(A / Mm[i], Kh[i]);
        const double HH = ((Vv[i] - Vm[i]) * AM) / (AM - Vm[i] / Vv[i]);
        const double Q  = (Vv[i] - HH) * s.psi;   // Pi_i of S. Eq 1
        s.growth[i] = dt * pop[i] * Q;
    }
    double sv = 0.0, sp = 0.0;
    for (int i = 0; i < N_POP; ++i) { sv += Vv[i] * pop[i]; sp += Pp[i] * pop[i]; }
    s.dr = -eR * s.psi * sv;                      // S. Eq 3
    s.dB =  eR * s.psi * sp;                      // S. Eq 5
    s.dA = -A * d0 - KB * B * A;                  // S. Eq 4
    const int NTR = (int) rate.size();
    for (int k = 0; k < NTR; ++k)
        s.probs[k] = rate[k] * dt * pop[f[k]] * s.psi;   // factors in the order the R reference multiplies them
}

// ---- simulate_hr ------------------------------------------------------------
//
// Runs one replicate from time 0 to horizon_h and returns the recorded rows.
//
// Arguments
//   init          initial densities of the 20 states, cells per mL; length 20,
//                 in the order R0..R9, R0C..R9C
//   r0, A0, B0    initial resource (ug per mL), antibiotic (ug per mL) and
//                 beta-lactamase (units per mL)
//   v, vm, kH,    per-state v (per cell per hour), vmin (per cell per hour),
//   M, Pb         kappa (dimensionless), MIC (ug per mL) and p (units per ug
//                 per hour); length 20 each, in the same state order
//   from, to      source and target state of each transition; length NTR.
//                 They arrive 1-based (1..20), as R supplies them, and are
//                 converted to 0-based here.
//   rate_vec      rate constant of each transition, per cell per hour (mu,
//                 gamma or omega times the transition's multiplicity);
//                 length NTR. NTR is taken from this vector.
//   dt            Euler step, hours
//   horizon_h     simulated time, hours
//   transfer_interval_h
//                 time between serial transfers, hours
//   kmon, eR,     K, e, d and KB of the symbol table
//   d0, KB
//   dil           factor by which the populations and the antibiotic are
//                 multiplied at each transfer, dimensionless
//   AMAX          antibiotic added at each transfer, ug per mL
//   Cres          resource concentration set at each transfer, ug per mL
//   n_out         approximate number of evenly spaced grid rows to record in
//                 addition to the transfer rows (the exact stride is given
//                 under Return value); a value <= 0 records only time 0 and
//                 the transfer rows
//   seed          seed of this run's random stream: a non-negative whole
//                 number. It is passed as a double because R has no 64-bit
//                 integer type, and is converted to uint64_t once. The stream
//                 depends on this value alone.
//   stochastic    false disables the Monte Carlo layer: no transition fires
//                 and no random number is drawn
//
// Preconditions
//   Not checked here: the six per-state vectors have length 20; from and to
//   have the length of rate_vec and values in 1..20; seed is a non-negative
//   whole number. The R wrapper run_sim() does not check them either: it
//   passes on what it finds in its parameter list, transition table, initial
//   state and seed. The caller must ensure them.
//   Checked here: n_out must not exceed n_steps, and n_rec must not exceed
//   1e8. n_rec is the number of rows allocated for the output, an upper
//   bound on the number of rows returned: the grid rows, plus
//   ceil(horizon_h / transfer_interval_h), plus 5. Either violation stops
//   with an R error.
//
// Return value
//   A numeric matrix with 25 named columns, in this order: time (hours),
//   R0..R9, R0C..R9C (cells per mL), r, A, B, cycle_end. Rows are in time
//   order. The first row is the initial state at time 0. Each transfer adds
//   one row holding the state just before dilution, with cycle_end = 1. If
//   n_out > 0, a grid row with cycle_end = 0 is added every out_every steps,
//   where out_every = max(1, floor(n_steps / n_out)); on a step that also has
//   a transfer the grid row follows the transfer row and holds the state
//   after the transfer. Two attributes are attached: "final_state" (length
//   23: the 20 densities, then r, A and B, after the last step) and
//   "n_steps".
// [[Rcpp::export]]
Rcpp::NumericMatrix simulate_hr(
    Rcpp::NumericVector init,
    double r0, double A0, double B0,
    Rcpp::NumericVector v,
    Rcpp::NumericVector vm,
    Rcpp::NumericVector kH,
    Rcpp::NumericVector M,
    Rcpp::NumericVector Pb,
    Rcpp::IntegerVector from,
    Rcpp::IntegerVector to,
    Rcpp::NumericVector rate_vec,
    double dt, double horizon_h, double transfer_interval_h,
    double kmon, double eR, double d0, double KB,
    double dil, double AMAX, double Cres,
    int n_out,
    double seed,
    bool stochastic = true)
{
    // The inputs are copied into plain C++ vectors once, so that the step loop
    // does not read them through R objects. Transition indices become 0-based
    // here.
    std::vector<double> pop(N_POP), Vv(N_POP), Vm(N_POP), Kh(N_POP), Mm(N_POP), Pp(N_POP);
    for (int i = 0; i < N_POP; ++i) {
        pop[i] = init[i]; Vv[i] = v[i]; Vm[i] = vm[i]; Kh[i] = kH[i]; Mm[i] = M[i]; Pp[i] = Pb[i];
    }
    double r = r0, A = A0, B = B0;
    const int NTR = rate_vec.size();
    std::vector<int> f(NTR), t(NTR);
    std::vector<double> rate(NTR);
    for (int k = 0; k < NTR; ++k) { f[k] = from[k] - 1; t[k] = to[k] - 1; rate[k] = rate_vec[k]; }

    const long long n_steps = llround(horizon_h / dt);
    // Grid rows are triggered by the step counter, every out_every steps, and
    // not by the transfer clock. The stride max(1, floor(n_steps / n_out)) is
    // the one the R reference uses, so both engines record the same steps.
    // n_rec is the number of rows allocated, an upper bound on the number
    // returned: the grid rows, ceil(horizon_h / transfer_interval_h) for the
    // transfer rows, and a margin of 5 that covers the row at time 0.
    const long long out_every = (n_out > 0)
        ? std::max(1LL, (long long) std::floor((double) n_steps / n_out)) : 0;
    const long long n_rec = (out_every > 0 ? n_steps / out_every : 0)
                          + (long long) std::ceil(horizon_h / transfer_interval_h) + 5;

    // An n_out above n_steps would give a stride of 1, that is one row per
    // step. It is rejected before the output matrix is allocated, as is an
    // n_rec above 1e8; n_rec is the number of rows that would be allocated.
    if (n_out > 0 && (long long) n_out > n_steps)
        Rcpp::stop("n_out (%d) exceeds n_steps (%lld): the recording grid would not thin. "
                   "Use 0 < n_out <= n_steps (R uses 400), or n_out <= 0 for transfer rows only.",
                   n_out, (long long) n_steps);
    if (n_rec > 100000000LL)
        Rcpp::stop("computed n_rec (%lld) is implausibly large; check n_out.", (long long) n_rec);

    Rcpp::NumericMatrix out(n_rec, N_COLS);
    long long rec = 0;
    auto record = [&](double now, int cycle_end) {
        if (rec >= n_rec) return;            // a row beyond the allocated bound is dropped, as in the R reference
        out(rec, 0) = now;
        for (int i = 0; i < N_POP; ++i) out(rec, i + 1) = pop[i];
        out(rec, 21) = r; out(rec, 22) = A; out(rec, 23) = B; out(rec, 24) = cycle_end;
        ++rec;
    };

    Xoshiro256pp rng; rng.seed((uint64_t) seed);
    std::vector<int> net(N_POP);
    Snapshot snap; snap.growth.assign(N_POP, 0.0); snap.probs.assign(NTR, 0.0);

    record(0.0, 0);
    // tq is the transfer clock, in hours: the sum of dt over the steps taken
    // so far, less transfer_interval_h for each transfer carried out. See the
    // transfer block below.
    double tq = 0.0;

    for (long long step = 1; step <= n_steps; ++step) {
        // Every rate and probability is computed from the state at the start
        // of the step, before pop, r, A or B changes. The updates further
        // down therefore do not feed into one another within a step.
        snapshot_kernel(pop, r, A, B, Vv, Vm, Kh, Mm, Pp, f, rate, dt, kmon, eR, d0, KB, snap);
        for (int i = 0; i < N_POP; ++i) net[i] = 0;

        // Monte Carlo layer. Transitions are visited in table order. A
        // probability <= 0 (for example an empty source state or a zero rate)
        // can never fire and a probability >= 1 always fires, so neither
        // consumes a random number; every other transition consumes exactly
        // one. Each firing moves one cell per mL from source to target.
        if (stochastic) {
            for (int k = 0; k < NTR; ++k) {
                const double pk = snap.probs[k];
                if (pk <= 0.0) continue;
                const bool fire = (pk >= 1.0) ? true : (rng.unif() < pk);
                if (fire) { net[t[k]] += 1; net[f[k]] -= 1; }
            }
        }

        // All increments are applied together. The antibiotic and the
        // populations are floored at zero: a fired transition removes a whole
        // cell per mL even from a state that holds less than one, and the
        // explicit antibiotic step goes below zero if d + KB * B exceeds
        // 1 / dt. The resource and beta-lactamase are not floored; they are
        // left exactly as computed, as in the R reference. The populations
        // are floored only after growth and transitions have both been added.
        for (int i = 0; i < N_POP; ++i) pop[i] += snap.growth[i] + net[i];
        r += dt * snap.dr;
        A  = std::max(0.0, A + dt * snap.dA);
        B += dt * snap.dB;
        for (int i = 0; i < N_POP; ++i) if (pop[i] < 0.0) pop[i] = 0.0;

        const double now = step * dt;                // time at the end of this step, as a product and not a running sum

        // Serial transfer. The clock tq is advanced by adding dt at every
        // step, as the R reference does, and is not derived from the step
        // counter, so that both engines compare the same rounded sum with the
        // interval. With dt = 1e-5 h and a 24 h interval the sum after
        // 2,400,000 additions is just below 24, so the first transfer falls
        // at step 2,400,001, at a recorded time of 24.00001 h, and a 24 h run
        // contains no transfer. The interval is subtracted from the clock,
        // which is not set to zero, so the remainder carries into the next
        // cycle.
        tq += dt;
        if (tq >= transfer_interval_h) {
            record(now, 1);                          // state before dilution: the end-of-cycle row
            for (int i = 0; i < N_POP; ++i) pop[i] *= dil;
            A = A * dil + AMAX;  r = Cres;  B = 0.0;
            tq -= transfer_interval_h;
        }
        // A grid row on a transfer step is recorded after the transfer.
        if (out_every > 0 && step % out_every == 0) record(now, 0);
    }

    Rcpp::NumericMatrix res = out(Rcpp::Range(0, rec - 1), Rcpp::_);   // out has n_rec rows, an upper bound; only the first rec were filled
    Rcpp::CharacterVector cn(N_COLS);
    cn[0] = "time";
    for (int i = 0; i < 10; ++i) {
        cn[1  + i] = "R" + std::to_string(i);          // 0-based columns 1..10: R0..R9
        cn[11 + i] = "R" + std::to_string(i) + "C";    // 0-based columns 11..20: R0C..R9C
    }
    cn[21] = "r"; cn[22] = "A"; cn[23] = "B"; cn[24] = "cycle_end";
    res.attr("dimnames") = Rcpp::List::create(R_NilValue, cn);

    // The state after the last step, and the number of steps taken, are
    // attached as attributes for inspection.
    Rcpp::NumericVector final_state(23);
    for (int i = 0; i < N_POP; ++i) final_state[i] = pop[i];
    final_state[20] = r; final_state[21] = A; final_state[22] = B;
    res.attr("final_state") = final_state;
    res.attr("n_steps") = (double) n_steps;
    return res;
}

// ---- debug_step1 ------------------------------------------------------------
//
// Evaluates snapshot_kernel() once for the supplied state and returns what it
// computed. The state is not advanced and no random number is drawn.
//
// Arguments
//   As for simulate_hr(), with the same units, lengths and 1-based from and
//   to; horizon_h, transfer_interval_h, dil, AMAX, Cres, n_out, seed and
//   stochastic are absent. to does not enter the calculation and is only
//   returned.
//
// Preconditions
//   As for simulate_hr(); none is checked here.
//
// Return value
//   A named list: psi (scalar), growth (length 20, cells per mL), dr, dB and
//   dA (scalars: the time derivatives of r, B and A), probs (length NTR), and
//   from_1based and to_1based (the from and to vectors as supplied, still
//   1-based, so that R can compare them with its own).
// [[Rcpp::export]]
Rcpp::List debug_step1(
    Rcpp::NumericVector init, double r0, double A0, double B0,
    Rcpp::NumericVector v, Rcpp::NumericVector vm, Rcpp::NumericVector kH,
    Rcpp::NumericVector M, Rcpp::NumericVector Pb,
    Rcpp::IntegerVector from, Rcpp::IntegerVector to, Rcpp::NumericVector rate_vec,
    double dt, double kmon, double eR, double d0, double KB)
{
    std::vector<double> pop(N_POP), Vv(N_POP), Vm(N_POP), Kh(N_POP), Mm(N_POP), Pp(N_POP);
    for (int i = 0; i < N_POP; ++i) {
        pop[i] = init[i]; Vv[i] = v[i]; Vm[i] = vm[i]; Kh[i] = kH[i]; Mm[i] = M[i]; Pp[i] = Pb[i];
    }
    const int NTR = rate_vec.size();
    std::vector<int> f(NTR);
    std::vector<double> rate(NTR);
    for (int k = 0; k < NTR; ++k) { f[k] = from[k] - 1; rate[k] = rate_vec[k]; }

    Snapshot snap; snap.growth.assign(N_POP, 0.0); snap.probs.assign(NTR, 0.0);
    snapshot_kernel(pop, r0, A0, B0, Vv, Vm, Kh, Mm, Pp, f, rate, dt, kmon, eR, d0, KB, snap);

    Rcpp::NumericVector growth(N_POP), probs(NTR);
    for (int i = 0; i < N_POP; ++i) growth[i] = snap.growth[i];
    for (int k = 0; k < NTR; ++k)   probs[k]  = snap.probs[k];

    return Rcpp::List::create(
        Rcpp::_["psi"] = snap.psi, Rcpp::_["growth"] = growth,
        Rcpp::_["dr"] = snap.dr, Rcpp::_["dB"] = snap.dB, Rcpp::_["dA"] = snap.dA,
        Rcpp::_["probs"] = probs,
        Rcpp::_["from_1based"] = from, Rcpp::_["to_1based"] = to);
}

// ---- debug_step1_fired ------------------------------------------------------
//
// Evaluates snapshot_kernel() once for the supplied state and applies the
// Monte Carlo rule once, with a generator seeded from seed. The rule and the
// draw order are those of simulate_hr(), so for the same inputs and seed the
// result is the set of transitions simulate_hr() fires in its first step.
//
// Arguments
//   As for debug_step1(), plus seed, which is treated as in simulate_hr().
//   to does not enter the calculation.
//
// Preconditions
//   As for simulate_hr(); none is checked here.
//
// Return value
//   A logical vector of length NTR in transition-table order: TRUE where the
//   transition fired.
// [[Rcpp::export]]
Rcpp::LogicalVector debug_step1_fired(
    Rcpp::NumericVector init, double r0, double A0, double B0,
    Rcpp::NumericVector v, Rcpp::NumericVector vm, Rcpp::NumericVector kH,
    Rcpp::NumericVector M, Rcpp::NumericVector Pb,
    Rcpp::IntegerVector from, Rcpp::IntegerVector to, Rcpp::NumericVector rate_vec,
    double dt, double kmon, double eR, double d0, double KB,
    double seed)
{
    std::vector<double> pop(N_POP), Vv(N_POP), Vm(N_POP), Kh(N_POP), Mm(N_POP), Pp(N_POP);
    for (int i = 0; i < N_POP; ++i) {
        pop[i] = init[i]; Vv[i] = v[i]; Vm[i] = vm[i]; Kh[i] = kH[i]; Mm[i] = M[i]; Pp[i] = Pb[i];
    }
    const int NTR = rate_vec.size();
    std::vector<int> f(NTR);
    std::vector<double> rate(NTR);
    for (int k = 0; k < NTR; ++k) { f[k] = from[k] - 1; rate[k] = rate_vec[k]; }

    Snapshot snap; snap.growth.assign(N_POP, 0.0); snap.probs.assign(NTR, 0.0);
    snapshot_kernel(pop, r0, A0, B0, Vv, Vm, Kh, Mm, Pp, f, rate, dt, kmon, eR, d0, KB, snap);

    Xoshiro256pp rng; rng.seed((uint64_t) seed);
    Rcpp::LogicalVector fired(NTR);
    for (int k = 0; k < NTR; ++k) {
        const double pk = snap.probs[k];
        // No draw is made when pk <= 0 or pk >= 1, as in simulate_hr().
        fired[k] = (pk <= 0.0) ? false : (pk >= 1.0) ? true : (rng.unif() < pk);
    }
    return fired;
}
