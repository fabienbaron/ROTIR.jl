# ROTIR — outstanding fixes

Compiled 2026-10-03. Everything here is backed by a measurement or a verified code reference;
where something is inferred rather than observed it says so.

Ordered by **how the failure presents**, not by effort: a wrong number that looks right is
worse than a crash, and a crash is worse than an ugly message.

---

## Context: what was fixed and verified on 2026-10-03

No action needed on these; listed so the items below are not re-derived.

- `_fit_nested` silently dropped `threaded`; `fit_sphere_ld`, `fit_ellipsoid_ld` and
  `fit_orbit` each dropped it on the way in. All now accept and forward it, and since the
  `:t3` pool landed the default is `true` again on all four.
- **`:t3` IS NOW RE-ENTRANT** — the mutable scratch moved out of the globally cached
  `Type3Plan` into a per-plan pool (`t3_checkout!`/`t3_checkin!` in a `try`/`finally`), and the
  four `Threads.@threads :static` schedules became `:dynamic` over `p.nchunks` chunks. This
  closed P1 items 1, 2 and 5, P2 item 9 and P3 entirely; see the design record below.
- `src/fit_pigeons.jl` had an unterminated string — the file did not parse, so
  `ROTIRPigeonsExt` never loaded. Fixed and committed.
- `component_vis` ported to OITOOLS 0.15's model API (the standalone `visibility_*` functions
  were removed upstream in 0.13.2). Pinned by a regression test against the originals.
- `rpole = 0` crashed Makie's tick formatter through an all-NaN temperature map; the blocking
  validator, `_map_range`, `map_colors` and `show_binary3d!` were all hardened.
- GUI parameter labels were elided at `dp(118)`; now `dp(168)`.
- χ² observable weights `(V², T3amp, T3φ)` thread through the primitive, its rrule, the three
  log-posterior builders and both samplers.
- The GUI's Pigeons runners share `PIGEONS_GUI_OPTS` (stabilized variational,
  `first_tuning_round = 2`). Measured on α Cen A: 2 round trips → 62 with a variational
  reference, 40 000× smaller barrier at the reference end.

Suite state: main **2892/2892** (OITOOLS 0.15), GUI **1350/1350**, GUI samplers **73/0/0**,
nested-threading **11/11**.

---

## P1 — silent wrongness

A user gets a number that is wrong and looks fine. The big one here is **done**; what is left
is item 4, which is what would make the *next* one of these audible instead of silent.

### 1-2. `:t3` re-entrancy — **DONE 2026-10-03**, and not with a guard

The scratch pool superseded the proposed throw-on-concurrent-entry guard: a guard would have
made the failure audible, the pool makes it not happen. Kept here because the measurement is
the reason the design changed.

**What the bug was.** `type3_plan_for` returns a **process-wide cached** `Type3Plan`
(`TYPE3_PLANS_LOCK` guarded construction, not use), and the plan carried `grids`, `scratch`,
`centred` and `realbuf`, which `type3_quads!` opens by zeroing. Two tasks on one plan, 2000
calls each, three trials: 2000 throws, 1999 correct, and **1 silently wrong** each time, worst
relative error 1.30–4.34. The window was structural — the zeroing and `_t3_finish!`'s
reduction/FFT run serially, outside any threaded region, so Julia's own `:static` precondition
passed and the second entrant walked onto the first's scratch. A finite, wrong visibility,
which an objective with `catch → return 1e30` *accepts*.

**What shipped.** `T3Scratch{T}` holding the four mutable buffers; `Type3Plan` keeps only
read-only plan data plus `nchunks`, `pool` and `poollock`. `type3_cvis`/`type3_cvis_adj!`
borrow one set for the whole call in a `try`/`finally`; the three exported entry points
(`type3_quads!`, `type3_quads_adj!`, `type3_points!`) keep their old signatures as fixed-arity
wrappers that borrow for the duration. `T3_POOL_MAX = 32` bounds steady-state memory (a set is
0.09 MiB at the `nf = 24` a lam And mesh plans).

**Verified.** Bit-identical to the pre-pool kernel — byte-for-byte, not to a tolerance — on
forward, adjoint and the point-set route at HEALPix 3 and 5, by running the same script against
both files. Clean over 2000 concurrent `fused_cvis` calls and 960 concurrent gradients where
the old code corrupted about 1 in 2000. Pinned by a new testset in `test/test_type3_nufft.jl`
("re-entrancy: concurrent transforms on one cached plan"), which asserts bit-identity against
the serial answer, that the pool does not leak a set when the transform throws, and the cap.

Two beliefs stated earlier in this file's history were wrong and shaped the design, so:

- **`:static` was not required for `grids[tid]`.** There is no `Threads.threadid()` anywhere in
  the file. `tid` is the loop index, and both the grid and the `lo`/`hi` quad block are keyed on
  it, so each iteration owns its grid whichever OS thread runs it.
- **The `:static` throw was not about shared scratch.** `jl_in_threaded_region` is a *global*
  counter, so two concurrent transforms on **different** plans, sharing nothing, also threw.
  The scheduler change was needed independently of the scratch fix, and they had to land
  together: shared scratch with `:dynamic` gave **0/8 correct, worst relative error 1.03e+28**.

**Cost.** Checkout/checkin is 374 ns against a 206 µs call (0.18 %), with no measurable
contention (16 tasks on one pool, 89 ns/pair aggregate). End to end through `fused_cvis` the
single-call cost is **unchanged within noise** at 16 threads (0.108/0.248/0.806 ms at HEALPix
3/4/5 against 0.105/0.232/0.799 before) — the earlier prototype's 5–17 % single-caller *win*
was measured at the plan level and does not survive end to end, so the pool should be taken as
free, not as faster. What it buys is concurrency: 8 simultaneous transforms complete **1.96×**
faster than the same 8 in a loop, on top of the ~4.8× the inner threading already gives.

**Rejected alternatives, with the measurement:**

- **Fresh scratch per call** — correct and trivial, but **+24 to +28 %** at HEALPix 2–4 (+5–7 %
  at 5–6): the scratch is mesh-independent, so a flat cost hurts most where `:t3` is cheapest.
  Sustained 27.6 µs per set, 84 % of it GC.
- **Lock the plan** — bit-identical, ~3 lines, and it does silence the `:static` error, but
  throughput *degrades*: 3834 → 3819 → 3405 → **3084** calls/s for 1, 2, 4, 8 callers. The right
  emergency patch, the wrong destination.
- **Task-local storage** — 358 ns lookup (no win) and the lifetime is wrong both ways:
  `@spawn`-per-call samplers get a fresh task every time, degenerating to per-call allocation;
  long-lived workers accumulate one set per plan touched with nothing to free them.
- **Per-thread scratch keyed beyond `threadid`** — that *is* the pool with a fixed backing
  store, and `threadid()` is unsound under task migration on Julia ≥ 1.7 anyway.

**Incidental, worth keeping:** `:dynamic` is free-to-faster, never slower (−11.7 %, −16.5 %,
−7.3 %, −6.6 %, +2.7 % at HEALPix 2–6) because `:static` pins iteration *i* to thread *i*, so
one contended thread stalls the whole barrier — exactly what a multi-chain sampler creates. And
**no `FFTW.UNALIGNED` is needed**: `plan_fft!` records the buffer's alignment class, which
looked like a blocker for per-call buffers, but 400/400 fresh 72×72 `ComplexF64` arrays measured
class 0, so the cached plan applies to any per-call `Array`. There is an `FFTW.alignment_of`
assertion at scratch construction to keep that honest; `UNALIGNED` would cost +6 µs per
transform.

### 3. `bootstrap_parametric` at its default `parallel = :outer` — **mostly closed**

The failure mode was: OITOOLS runs replicates under `Threads.@spawn`, `run_chunk` swallows the
throw, every replicate NaNs out, and you get `"fewer than 2 usable replicates"` or a σ computed
from whichever handful happened not to collide. With `:t3` re-entrant, a `:t3`-route `fitfun` no
longer throws there.

What remains is the swallowing itself, which is upstream and applies to **any** `fitfun` that is
not re-entrant — including one a user writes. Worth a line in `bootstrap_parametric`'s docstring
saying `parallel = :outer` requires a re-entrant `fitfun` and will quietly drop replicates
otherwise. Not worth code.

### 4. `shell_fit`'s bare `catch → return 1e30`

This is what turned a hard failure into a plausible fit with an evidence: 40 000 of 40 000
likelihood calls hit the sentinel, and the panel reported `log(Z) = -5e29` (finite),
`N_eff = 1`, the flat prior as the posterior, and a median 1.402 against the true 1.37717.

A derivative-free objective genuinely needs to survive geometry failures, so this is a design
call, not a mechanical fix. Minimum viable: refuse a run in which *every* evaluation failed.
Better: console-log the first exception, then continue.

### 5. `type3_quads_adj!` — **DONE**, same change

It shared `centred`/`scratch`/`realbuf` and had the same loop; it takes the scratch from the
same pool, is bit-identical, and the new testset covers it (20 rounds × 8 concurrent callers,
all bit-identical to serial). It has no production caller yet — `compute_adjoint_cvis!` has no
`:t3` branch on purpose, see `src/fused_polyft.jl:348` — but it is exported and precompiled.

---

## P2 — loud failure, or claims in our own source that are wrong

### 6. The timing comment at `src/gui/shell.jl` — **DONE**

It claimed serial `:t3` was 0.20 ms per evaluation at HEALPix 3; that was the 16-thread figure
and the serial one is 0.94 ms, so the comment argued for serialising on a number that only held
if it did not serialise. The whole block is rewritten: the branch now passes `threaded = true`
and the comment records what the bug was and how the pool closed it.

### 7. `kernel = nothing` should be spelled `:scalar`

`shell.jl:3631` does `kernel === nothing ? body : with_polyft_backend(…)`, so `nothing` means
no scope and the worker runs the process default `:t3` — while the console says *"this fit runs
the reference kernel instead"*.

**There is no live bug**: `_fit_kernel` returns `:turbo` only for
`KERNEL_TURBO_METHODS = (:gradient, :hmc, :pigeons)`, and those three route through
`interferometric_chi2` / `shape_chi2_fg!`, which have no `:t3` branch — verified empirically
by probing `length(TYPE3_PLANS)` (0 plans cached, bit-identical to forced `:scalar`). Pigeons in
a PackageCompiler bundle is therefore safe.

But the message's accuracy rests entirely on that tuple's contents. Adding `:neldermead`,
`:bobyqa` or `:nautilus` to it — plausible, since `:turbo` is the better Float32 kernel at
HEALPix 3–4 — would make the console line false *and* the fit non-re-entrant. Spelling the
fallback `:scalar` makes it true by construction. An assertion in `_fit_kernel` is the
alternative.

### 8. The GUI thread is not quiesced during a fit — **downgraded to cosmetic**

`start_job!` refuses a second job, so job-vs-job is impossible. But only **6** file-IO buttons
carry `enabled: !win.jobRunning` in `src/gui/qml/Main.qml`; the parameter fields, epoch stepper,
tab bar and kernel combo stay live. `shell_set_param` clears `sh.chi2key[]`, and the next
refresh runs `epoch_chi2` → `model_state` → `observables` → `fused_cvis` → `:t3` **on the GUI
thread** while the worker is also on `:t3`.

That used to hand the worker a 1e30 sentinel at that trial point. With the pool the two simply
both run, correctly. What is left is that a user can edit the parameters underneath a running
fit and be confused by the result, which is a UI question, not a correctness one. Extending
`enabled: !win.jobRunning` to the parameter form and epoch stepper is still the fix.

### 9. `spheroid_chi2_allepochs_f` threads over epochs into `:t3` — **DONE** by the pool

`src/oichi2_spheroid.jl:409` — `Threads.@threads for i=1:nepochs` → `spheroid_chi2_f` →
`observables`, which takes the matrix-free `:t3` branch whenever `setup_oi!` has not been run.
Exported, and latent only because both in-repo callers happen to run `setup_oi!` first, which
the function cannot know. Concurrent `:t3` from an epoch loop is now correct, so nothing is
needed. Same for item 12 below (`binary_cvis`'s matrix-free branch).

---

## P3 — re-entrancy: **CLOSED**

Everything that was in this section shipped on 2026-10-03. The design record, the rejected
alternatives and the numbers are in P1 item 1-2 above, which is now the single place `:t3`
re-entrancy is documented. One consequence worth stating on its own:

### 10. The three `threaded = false` defaults are back to `true`

`src/parametric_fit.jl` (`fit_sphere_ld`, `fit_ellipsoid_ld`), `src/orbit_fit.jl` (`fit_orbit`)
and the GUI's `_run_fit` nested-sampler branch each carried `threaded = false` solely because of
this bug. All four now default to `true`, and the stale `@threads :static` rationale in their
comments is replaced by what actually holds. `grep ':static' src ext` returns nothing.

### 11. `binary_cvis`'s matrix-free branch

`src/oichi2_binary.jl:168-170`, public, on the `:t3` route with the same `isempty(polyft)`
auto-switch. No threaded caller in-repo, and now safe if one appears.

---

## P4 — hygiene

### 12. CI has no test workflow

`.github/workflows/` contains only `docs.yml`. **`test/runtests.jl` has never run on push or
PR** — 2892 assertions that execute only locally. The main suite needs no display and would run
on `ubuntu-latest` as-is. This is the single biggest CI gap.

### 13. Documenter's `checkdocs = :none` hides a backlog

Undocumented or unreachable exports surface only when something `@ref`s them, at which point
the build fails hard — this broke the docs build three times in one session. Known:
`fused_cvis`, `binary_cvis`, `poly_to_cvis`, `limb_mu` are exported but in no `@docs` block
(`api/fused_polyft.md` documents by markdown table, which is the structural cause);
`provider_intensity` is exported with no docstring at all, and it is the single call the whole
provider interface exists to offer. Consider `checkdocs = :exports`, expecting a backlog.

### 14. The parse-scan test is syntax-only

`test/runtests.jl`'s "every shipped source file parses" uses `Meta.parseall`, which validates
syntax but does **not** lower. It would have caught the `fit_pigeons.jl` unterminated string
that motivated it, but it misses `const` on a local, undefined variables and bad signatures —
it did in fact miss exactly that in `test_nested_threading.jl`. Say so in its docstring so it is
not read as a general validity check.

### 15. Uncommitted work and the AppImage

The AppImage is built (562 600 440 B) and smoke-tested but **not deployed**, and corresponds to
no commit. The deploy target is
`joy:/nfs/morgan/users/baron/public_html/ROTIR-x86_64.AppImage` — a public URL; upload to a
`.new` name, verify SHA-256, then move into place.

### 16. ζ And

The ellipsoid fit railed against `f_bounds`' upper wall at 0.3997 (16–84 interval
[0.3993, 0.3999]) — and 0.40 is **past the physical maximum**, since a rigidly rotating star at
breakup has R_eq/R_pol = 3/2, i.e. flattening 0.333. It is absorbing the ~20 % spot signal, not
measuring a shape. Δlog Z = +15 834 should not be read as "ellipsoid preferred".

Rerun with `f_bounds = (0.0, 0.333)` to confirm it rails at the physical limit too. The real
answer is the joint model + imaging fit, where the spots have a surface to live on.

Established for that work: D = 2.5603 mas, α = 0.2004 from the symmetric fit, with
`corr(D, α) = +0.74`; honest errors are the per-epoch scatter (±0.0469 mas, ±0.0742), not the
posterior widths, because χ²ᵣ = 9.17 says the quoted error bars are wrong.

### 17. Chromatic imaging and a model-atmosphere intensity

Measured 2026-10-04 while reproducing Roettenbacher et al.'s Figure 10. Three numbers decide
the design, and none of them was obvious beforehand.

**(a) Limb darkening is MORE chromatic than temperature is, by 4x.** Kurucz ATLAS9 at
4500 K / logg 2.5, across MIRC-X's six H-band channels:

| lambda (um) | 1.5154 | 1.5510 | 1.5869 | 1.6228 | 1.6587 | 1.6947 |
|---|---|---|---|---|---|---|
| alpha | 0.2585 | 0.2466 | 0.2273 | 0.2249 | 0.2189 | 0.2170 |

alpha falls **16 %** across the band, which moves the brightness at mu = 0.3 by **+5.1 %**.
A 273 K spot's contrast moves **+1.24 %** over the same span. So a chromatic fit run with a
fixed scalar alpha would dump the limb-darkening chromaticity straight into the temperature
map and the artefact would be four times the signal. Chromatic work therefore REQUIRES a
lambda-dependent limb darkening; it is not an optional refinement.

**(b) But the power law is the right functional form.** V^2 of a limb-darkened disc computed
from the full Kurucz mu-profile against the best-fit mu^alpha, over zeta And's baselines
(x = pi*theta*B/lambda up to 8.3, past the third lobe): worst |dV2| = 0.00031 at x = 3.0,
against a median inflated sigma(V2) of 0.0012 — **0.26 sigma at worst per point**. The power
law only fails at mu < 0.05, which carries almost no projected area. What the grid adds is
alpha(lambda) and alpha(T), NOT a different law.

**(c) A Kurucz T -> I conversion is cheaper than Planck.** 768 tessels: `:linear` 0.13 us,
`:planck` 2.07 us, a precomputed uniform-T table lookup **0.98 us**. No `expm1` per tessel.
Cost is a non-issue against `fused_cvis`.

**(d) The chromatic signal in zeta And's own data is marginal.** Computing the full chromatic
model against a grey one on the real uv coverage, at a fixed map: Dchi2 = 33.8 on published
errors (12.4 on their 1.65-inflated ones), **91 % of it in T3phi**. With the over-contrasted
published map it is 152.6 / 56.1. Treat these as UPPER BOUNDS — they compare two models at a
fixed map, and a real fit absorbs much of the difference by adjusting the grey map.

#### Decisions

**No new `ldtype` number.** Three reasons. (i) `ldtype` selects a parametric law whose
coefficients live in the schema as `ld1..ld4`; a grid has no coefficients, so
`ld_coefficients_used` would have nothing to return and the GUI dropdown would need a file
picker. (ii) `ldtype = 0` already means "none — the intensity provider supplies the mu
dependence", with `owns_mu(provider)` as the double-counting guard; a second spelling would
let `ldtype = 5` + a provider double-count, needing a new `validate_star_params` rule to
police a problem we invented. (iii) measurement (b) says the grid's value is alpha(lambda),
which fits inside `ldtype = 3` by setting `ld1` per channel.

So instead: **`kurucz_ld_alpha(grid, Teff, logg, lambda; mu_min = 0.05) -> alpha`**, an
area-weighted power-law fit returning the number to put in `ld1` for that channel.

**`:kurucz` as an intensity law, built from the existing grid.** `intensity(x, model, band)`
takes `model::Symbol`; relax that annotation to accept a law object as well, and every call
site works untouched — including the 8 `xw` sites and `spheroid_chi2_fg`'s `intensity_model`
/`band` plumbing. The law:

    struct GridIntensity{T}
        Tgrid::Vector{T}   # UNIFORM, so lookup is multiply-and-floor, not a search
        I::Vector{T}       # I(T) at (logg, mu = 1, lambda), normalised
        slope::Vector{T}   # dI/dT per interval — the rrule reads this
        logg::T; lambda::T # provenance, for the guard
    end

Build it FROM the `RectGrid4`/`TabulatedProvider`, do not reimplement the interpolation —
one grid, one source of truth (the same reason there is no `ParametricLDProvider`). The
existing `intensity` rrule is diagonal in `x` and a table lookup's derivative is also
diagonal, so no new rule SHAPE is needed. Make the type grid-agnostic so a Korg/MARCS grid
serves a cool dwarf identically; `:kurucz` is only how you obtain it.

**Guard the provenance.** A table built at the wrong lambda or logg gives a plausible wrong
answer that survives every relative test — the same failure mode as the Kurucz
ratio-vs-intensity trap. The chi2 path must check `data`'s band against the table's lambda and
refuse a mismatch, in the spirit of the existing "`:planck` requires `band`" error.

**Scope it to scalar logg, and say so.** The table is 1-D in T. alpha moves 0.004 across
logg 2.0-3.0 and a 6 %-distorted giant spans ~0.05 dex, so this is fine for zeta And. It is
NOT fine for beta Cas, whose surface spans 0.47 dex; that case is the provider path.

#### Order

The three steps as proposed, with one gate inserted and one GUI step folded in:

1. **`kurucz_ld_alpha`** — small, immediately useful, and what the Fisher gate in step 3 needs.
2. **`GridIntensity` + relaxing `intensity`'s `model::Symbol`**, with an FD gradient test
   beside the existing `:planck` one. Then the GUI step below, which is cheap and is what makes
   steps 1-2 usable by anyone but a script.
3. **A Fisher gate — INSERTED, not in the original plan, and it gates step 4.** sigma(T_spot)
   on zeta And's actual coverage with the grid supplying alpha(lambda), following
   `demos/betcas_fisher_predicted_ld.jl`. Measurement (d) above is an upper bound taken at a
   fixed map; this says how much survives marginalising over a grey map. Expectation: H alone
   is marginal; J+K gives ~4.5x the lever, since dlnI/dlnT runs 2.84 (1.25 um) to 1.93 (2.2 um).
   If the gate says marginal, step 4 is a multi-band proposal argument rather than code to
   write now.
4. **The channel loop** in `spheroid_crit_allepochs_fg` (`src/oichi2_spheroid.jl:460`):
   per-channel `polyft`, per-channel flux normalisation, and
   `dchi2/dT_i = sum_{e,c} (dchi2_{e,c}/dI_{i,c}) * (dI/dT)`. The per-term pieces exist; the
   accumulation is new. `readoifits(polychromatic = true)` already works on these files and
   gives each channel its own lambda-scaled uv.

#### A layout case that does not exist yet

`parametric_layout` (`src/bootstrap.jl:395`) drops `ld1`/`ld2` only when `ldtype == 0`.
Predicted limb darkening at `ldtype == 3` — alpha supplied per channel by `kurucz_ld_alpha`
rather than fitted — is a NEW case: the law is still the power law, so `ldtype` stays 3, but
`ld1` must leave theta or the fit will move a number the grid already decided. Add a
`ld_predicted::Bool` flag, and make `ld_coefficients_used` (or its callers) agree, because the
GUI reads it to decide which rows are editable.

#### GUI

Surveyed 2026-10-04; these are the exact seams.

**(1) The law selector is a one-line gate.** `src/gui/shell.jl:1912` reads
`m in (:linear, :planck) || return "unknown intensity model"`, with `sh.intensity_model[]` a
`Ref{Symbol}` (`:91`). `:kurucz` is a third entry there plus a grid to point at.

**(2) There is no file dialog anywhere in the GUI.** So the grid arrives as a PATH in a text
field, pointing at the derived FITS grid (0.3 MB for one band) rather than the 67 MB ATLAS9
pack — `read_kurucz_intensity` already persists to FITS, and `demos/spica_kurucz_grid.jl`
establishes `ARGS[1]` / `KURUCZ_PACK` / `demos/data/` as the convention. Beside the field,
show the grid's own coverage (Teff, logg, lambda ranges) and a live verdict on whether the
current surface falls inside it. A silently extrapolated grid is the failure mode worth
spending a label on.

**(3) `tpole` must join the `:planck` condition, and this one will bite.** `shell.jl:3394`
frees `tpole` only under `:planck`, because under `:linear` a uniform scale on the map divides
out of a normalised visibility. A grid intensity is nonlinear in T for the same reason Planck
is, so `:kurucz` belongs in that condition. Add it to `:1912` and forget `:3394` and `tpole`
stays silently fixed — the fit then cannot set its own flux scale and will absorb the error
elsewhere. Same for `:1919`'s `band <= 0` guard: a grid needs a band too.

**(4) The display already switches units.** `shell.jl:2011` returns `"I (arb.)"` under
`:planck` and `"T (K)"` otherwise, and `:1953` maps the map through `intensity` for rendering.
`:kurucz` reuses both. Under chromatic imaging the reconstructed quantity really is T, so the
label is right without further work.

**(5) Limb darkening needs a "predict from grid" checkbox, not a new dropdown entry.** The LD
rows' editability is driven by `ld_coefficients_used` (`:638`, `:1307`, `:1432`). With the box
ticked, `ld1` becomes read-only and shows the per-channel alpha the grid gives — a short table
or six values, since the whole point is that alpha moves 16 % across the band. Keep `ldtype`
at 3: this is a predicted coefficient, not a different law. See the layout case above.

**(6) Channel awareness is a real change, and the GUI currently throws the axis away.**
`shell.jl:357` exports `readoifits_multiepochs(files; ...)[1, :]` — row 1 of an
`(nwav, nepochs)` array, i.e. the first channel only. Chromatic mode means dropping that
index, and the panel should then state the channel count and their lambda range, because the
forward model becomes nchannel times more expensive. Six channels is 6x the fit time and the
user should learn that before pressing Fit, not after.

**(7) Make the chromaticity visible.** One small plot — alpha vs lambda from the grid, and the
spot contrast vs lambda — because the governing fact of this whole item is that the first is
four times the second, and a GUI that hides it invites exactly the mistake of fitting a
chromatic map with a fixed alpha.

**(8) Guards belong in `validate_star_params`.** It doubles as the GUI's build gate, and a
check moved out of it stops being enforced there (this has happened once already). So
"provider implies `ldtype = 0`", "a grid intensity needs a band", and "the surface lies inside
the grid" all go there rather than in GUI-side code.

**One trap to record:** the lambda^5 Planck prefactor that `src/intensity.jl` drops is
HARMLESS for per-channel normalised visibilities — it is common to every tessel at one lambda
and cancels in `F / dot(polyflux, xw)`. It must be restored before anyone fits an absolute
spectrum. Do not "fix" it blindly.

**Also worth noting**, separately: zeta And's 300:1 companion is relatively brighter at longer
lambda, so the literal SPARCO case — one image plus a parametric companion with its own
spectral index — applies here too and would reuse `src/oichi2_binary.jl`. Simpler job,
different payoff.

---

### 18. An optional calibration scale `V0`

Roettenbacher et al.'s Table 6 fits a "visibility at the origin", V0 = 1.02 +- 0.02. ROTIR has
no equivalent: `compute_polyflux_and_cvis!` divides by `dot(polyflux, xw)`, so the model
visibility is exactly 1 at zero baseline by construction. Freeing a scale on zeta And recovers
**1.020** (uniform surface) and **1.022** (their map) — their published number, from a
parameter we were not told existed. If that 2 % is real calibration error rather than a fitting
artefact, then every ROTIR fit is currently absorbing it into geometry and limb darkening.

**This is a plan, not a task list. Nothing below is implemented.**

#### D1 — where it multiplies: the COMPLEX visibility

`cvis -> V0 * cvis`, once, where cvis is formed. Then V^2 scales by V0^2, T3amp by V0^3,
T3phi is unchanged (V0 real positive), visamp by V0, visphi unchanged. This is the only choice
self-consistent across all seven observables; scaling V^2 alone would contradict T3amp.

#### D2 — scalar or per-dataset, and why V0 must NOT go in `star_params`

Physically each night has its own calibration, so V0 is per-DATASET. Their fit uses one global
value across three epochs. Support scalar first with the plumbing shaped for a vector.

The consequence is structural: `create_star_multiepochs` takes ONE `star_params` and makes N
stars, so a per-epoch V0 cannot live in the schema. Keeping it out also keeps it away from
`validate_star_params`, `derived_quantities` and `surface_geometry_io`, where a non-physical
nuisance parameter does not belong. **Decision: V0 is a calibration parameter in the layout,
never a schema field.**

#### D3 — how it enters theta

`ParametricLayout` (`src/bootstrap.jl:353`) maps each theta slot to a `star_params` FIELD, and
`layout_merge` (`:422`) merges every field into the NamedTuple. A calibration entry has no
field. Add a kind marker (or a sentinel symbol) plus `ncalib`, make `layout_merge` skip
calibration slots, and add `layout_calib(L, theta) -> Vector`. `layout_theta` (`:414`) needs
the mirror. Both are exported and pinned by tests but have **no in-tree callers**, so the blast
radius is the tests, not the samplers.

`_PARAMETRIC_ENTRIES` (`:365`) gains a `V0` entry behind a `v0_free` flag in
`parametric_layout` (`:395`), bounds `[0.8, 1.2]`.

#### D4 — prior, and why it is not optional

A free V0 with no prior will eat real signal, especially where the data do not reach the first
null. Use the existing `logprior` hook with N(1, 0.02) rather than hard-wiring one — but the
DEFAULT when `v0_free = true` should be that prior, not a flat box, and the docstring should
say why.

#### D5 — identifiability

V0 is degenerate with a uniform halo, which a tessellated surface cannot represent, so it IS a
new direction here. It is partially degenerate with theta_LD when the data stop short of the
first null. **Gate before trusting any fitted V0: a Fisher corr(V0, rpole) on the target's own
coverage.** zeta And reaches x = 8.3 so the two separate well there; a short-baseline dataset
is a different story.

#### D6 — every function that must change

Besides `parametric_chi2` (`src/parametric_fit.jl:171`):

| what | where | change |
|---|---|---|
| `interferometric_chi2` | `src/parametric_gradient.jl:328` | `cvis = (V0/flux) .* F` |
| its rrule | `:352` | `adj_F = adj_cvis .* (V0/flux)`; new component `dchi2/dV0 = sum_k Re(adj_cvis[k] * F[k]/flux)`; the return tuple gains a slot |
| `build_parametric_logpi` | `:473` | `slot("V0")`, 0 meaning "held at 1" |
| `ParametricLayout` et al. | `src/bootstrap.jl:353,365,395,414,422` | D3 |
| `parametric_posterior_spec` | `src/parametric_fit.jl:50` | inherits once the layout carries V0 |
| `fit_parametric` | `ext/ROTIRZygoteExt.jl:45` | the hard theta-length assert |
| `bootstrap_parametric` | `ext/ROTIRZygoteExt.jl:150` | thread the flag — it silently dropped `provider`/`layout` once before |
| `_fit_hmc`, `_fit_pigeons`, `fit_parametric_ultranest` | `src/fit_hmc.jl:69`, `src/fit_pigeons.jl:123`, `src/ultranest.jl:45` | inherit via the spec; labels/reporting need the name |
| `fit_sphere_ld`, `fit_ellipsoid_ld` | `src/parametric_fit.jl:300,497` | their OWN theta vectors, not `ParametricLayout` |
| GUI | `src/gui/shell.jl:3185,3193,3576` | below |

**The reporting trap, and it is the dangerous one.** `chi2s` (`src/oichi2_spheroid.jl:335`),
`chi2_breakdown` (`:306`), `observables` (`:268`) and the residual plots all rebuild the model
WITHOUT V0. Left alone they would print a different chi2 from the one the optimiser minimised,
which looks like a bug in the fit rather than in the report. Either thread V0 through all of
them or refuse to accept a fitted V0 in them — do not leave the discrepancy.

#### D7 — GUI

Because V0 is deliberately not a schema field, the schema-driven parameter form will NOT show
it — and that is correct, because it belongs with the DATA, not the model. So it needs its own
control, and the right home is the fit panel beside the weights, not the parameter table.

`_fit_layout` (`shell.jl:3185`) gains `v0_free = :V0 in names`. `_theta_slot` (`:3193`)
resolves on `L.fields`, so it works unchanged given a sentinel symbol. `shell.jl:3576` calls
`parametric_chi2` directly for its own readout and must be passed the same V0 the fit used, or
the panel will disagree with the fit — the GUI-side instance of the reporting trap in D6.

Three presentation points:

* **Show it as a percentage, not a bare number.** "calibration +2.0 %" is readable; "V0 =
  1.0200" invites someone to think it is a property of the star. The control is a checkbox
  ("fit a visibility calibration") plus the fitted value and its error, greyed when unticked.
* **Show the prior, because a free V0 with none eats signal (D4).** The panel should state the
  prior it is using — "N(1, 2 %)" — rather than leaving the user to discover that the default
  is not a flat box.
* **Warn when it is unidentifiable, rather than silently absorbing geometry.** D5's degeneracy
  with `theta_LD` is real when the data stop short of the first null. The cheap proxy the GUI
  can compute without a Fisher matrix is the maximum `pi*theta*B/lambda` the loaded data reach;
  below ~4 (the first null of a uniform disc) say so next to the checkbox. zeta And reaches
  8.3, so it is fine there and the warning stays quiet.

Serialisation: `shell.jl:1891` already writes a tab-separated settings line carrying the
intensity flag, model and band. V0's state belongs in the same line, and the script-export
path (`:354`) has to emit the matching keyword, or a GUI fit and its exported script will
quietly differ.

#### D8 — minimisation

V0 is just another bounded parameter to OptimPack/NLopt. For NUTS, `_box_transform` maps
`[0.8, 1.2]` through a scaled logistic, which is the case it already handles. Scale is benign
(V0 ~ 1 beside rpole ~ 1 mas and tpole ~ 5000 K, a range the fitters already span).

#### D9 — the test that makes this safe

A null-space gate in the style of the `tpole` and `d` ones in
`test/test_parametric_gradient.jl:83-90`: with `weights = (0, 0, 1)`, i.e. **T3phi only**,
`dlogpi/dV0` must be EXACTLY 0.0, because closure phase is invariant under a real positive
scale. Non-zero means V0 has been wired into the phase path by mistake. Then the ordinary FD
check with `central_fdm(5, 1)` for the amplitude weights, and a round trip: scale synthetic
data by a known V0 and recover it.

#### D10 — what is explicitly OUT of scope for a first pass

The imaging path. `spheroid_chi2_fg` (`src/oichi2_spheroid.jl:378`),
`fused_spheroid_chi2_fg` and `binary_chi2_fg` carry three hand-written adjoints built on
linearity in `x`. V0 is a separate scalar so it is a NEW gradient component rather than a
change to the `x` adjoint — tractable, but three more adjoints to validate. Do the parametric
path first and let it prove the sign conventions.

#### Suggested order

1. `interferometric_chi2` + rrule + the FD and null-space tests. Nothing else can be trusted
   until the gradient is right.
2. `parametric_chi2` (value path) and the layout.
3. `fit_parametric` end to end on zeta And; confirm it returns 1.02.
4. Reporting (D6) — before the GUI, so the GUI cannot display an inconsistent chi2.
5. GUI.
6. Imaging, if the parametric result justifies it.

---

## What is left, in order

**P1 item 4** — `shell_fit`'s bare `catch → return 1e30`. The pool removed the exception this
was swallowing, but it will swallow the next one just as quietly. Minimum viable: refuse a run
in which *every* evaluation failed; better, console-log the first exception and continue.

**P2 item 7** — spell the Pigeons fallback kernel `:scalar` rather than `nothing`, so the
console line stays true if someone adds a threading kernel to that tuple.

**P4 item 12** — CI runs the docs build and nothing else. The suite is the thing most likely to
catch a regression and it does not run on a push.

Then the science: **item 16**, ζ And's joint model + imaging fit, which is what the ellipsoid
railing at the physical flattening limit is asking for.

**Item 18** (`V0`) is the smallest of the new work and the one with the widest reach — if the
2 % is real, every existing fit has been absorbing it. Do its gradient and null-space test
first; everything else in that item is plumbing.

**Item 17** (chromatic) is gated on a Fisher calculation, not on code. Run the gate before
writing the channel loop.
