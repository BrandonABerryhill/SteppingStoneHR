# gate1b_dump_step1.R
#
# Part of the tests of the hrsim package, the simulation engine for Berryhill
# and Andersson (2026), "Unstable Heteroresistance as a Stepping Stone to
# Stable Resistance via Compensatory Mutations".
#
# Writes tests/fixtures/step1_internals.rds, the reference values read by
# gate1b_check.R. The expressions below are those of the first step of the R
# reference engine (run_model() in 03_engine_R.R), evaluated once for the
# Scenario 1 condition at its initial state. No random number is drawn, so no
# seed is involved. Running this script overwrites the fixture. Run from the
# package root.

ref_dir <- getOption("hrsim.ref", "tests/ref")
source(file.path(ref_dir, "load_canonical.R"))

DT  <- p$DT
pop <- as.numeric(y0$pop); r <- y0$r; A <- y0$A; B <- y0$B
v <- p$v; vm <- p$vm; kH <- p$kH; M <- p$M; Pb <- p$P
base     <- c(gamma = p$gamma, mu = p$mu, omega = p$omega)
rate_vec <- as.numeric(transitions$mult * base[transitions$type])
from <- transitions$from; to <- transitions$to
kmon <- p$kmon; eR <- p$e

psi    <- r / (r + kmon)
AM     <- (A / M)^kH
HH     <- ((v - vm) * AM) / (AM - vm / v)
Q      <- (v - HH) * psi
growth <- DT * pop * Q
dr     <- -eR * psi * sum(v  * pop)
dB     <-  eR * psi * sum(Pb * pop)
dA     <- -A * p$d - p$KB * B * A
# Probability that each transition fires in the step. A run with the Monte
# Carlo layer disabled never uses these values, which is why they are tested
# here.
probs  <- rate_vec * DT * pop[from] * psi

step1 <- list(psi = psi, growth = growth, dr = dr, dB = dB, dA = dA,
              probs = probs, from = from, to = to, rate_vec = rate_vec,
              pop_from = pop[from])
dir.create("tests/fixtures", recursive = TRUE, showWarnings = FALSE)
saveRDS(step1, "tests/fixtures/step1_internals.rds")

cat(sprintf("step1 dumped. psi=%.9f growth[1]=%.6e growth[3]=%.4f dr=%.7f dB=%.7f dA=%g probs[1]=%.6e max(probs)=%.5f\n",
            psi, growth[1], growth[3], dr, dB, dA, probs[1], max(probs)))
