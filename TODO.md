# ROTIR — outstanding fixes

Compiled 2026-10-03, audited 2026-10-04. Everything here is backed by a measurement or a
verified code reference; where something is inferred rather than observed it says so.

Ordered by **how the failure presents**, not by effort: a wrong number that looks right is worse
than a crash, and a crash is worse than an ugly message.

Items 1-2, 5, 6, 9, 10, 11 and 14 were closed and have been removed; the numbering is kept so
older references still resolve. They were the `:t3` re-entrancy work and its consequences — the
design record lives in `src/type3_nufft.jl`'s own comments and in git history, not here.

---

## P1 — silent wrongness

A user gets a number that is wrong and looks fine.

### 4. `shell_fit`'s bare `catch → return 1e30`

This is what turned a hard failure into a plausible fit with an evidence: 40 000 of 40 000
likelihood calls hit the sentinel, and the panel reported `log(Z) = -5e29` (finite),
`N_eff = 1`, the flat prior as the posterior, and a median 1.402 against the true 1.37717.

A derivative-free objective genuinely needs to survive geometry failures, so this is a design
call, not a mechanical fix. Minimum viable: refuse a run in which *every* evaluation failed.
Better: console-log the first exception, then continue.

---

## P2 — loud failure, or claims in our own source that are wrong

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

---

## P4 — hygiene

### 3. `bootstrap_parametric`'s `parallel = :outer` is undocumented

The failure mode it had — OITOOLS runs replicates under `Threads.@spawn`, `run_chunk` swallows
the throw, every replicate NaNs out — is gone for `:t3` routes now that the pool landed. But the
swallowing is upstream and applies to **any** `fitfun` that is not re-entrant, including one a
user writes.

One sentence in the docstring (`src/bootstrap.jl:228`): `parallel = :outer` requires a
re-entrant `fitfun` and will quietly drop replicates otherwise. Not worth code.

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

### 15. The AppImage is built but not deployed

562 600 440 B, smoke-tested, and corresponds to no commit. Deploy target is
`joy:/nfs/morgan/users/baron/public_html/ROTIR-x86_64.AppImage` — a public URL, so upload to a
`.new` name, verify SHA-256, then move into place.

(The "uncommitted work" half of this item is closed: the working tree is clean.)

---

## Science and new capability

### 16. ζ And — **rewritten 2026-10-04; the shape question is settled**

The previous version asked to rerun the ellipsoid fit with `f_bounds = (0.0, 0.333)` because it
railed at 0.3997. **That run can no longer be made and no longer needs to be**:
`fit_ellipsoid_ld` refuses `f_bounds` (`src/parametric_fit.jl:430-432`) — the oblate
parameterization it belonged to was replaced by three free semi-axes, which is the right model
for a tidally distorted star and has no wall to rail against.

What is now established, from the data and independently from Roettenbacher et al.'s own
Figure 10 (`demos/zetand_fig10_match.jl`, `zetand_fig10_family.jl`, `zetand_alpha_audit.jl`):

| | their Figure 10 | ROTIR, their family | ROTIR, triaxial |
|---|---|---|---|
| r_major (mas) | 1.3264 | 1.3271 | 1.3275 |
| r_minor (mas) | 1.2479 | 1.2489 | ry 1.2663, rz 1.2120 |
| major/minor | 1.0629 | 1.0626 | — |
| inclination | 69.44° | 75.6° | 64.7° |
| PA | 124.30° | 120.8° | 122.8° |
| alpha | 0.269 (fixed) | 0.154 | 0.139 |
| chi2/n | — | 4.523 | 3.626 |

Two independent routes — interferometry, and a rasterised PDF page — agree on the ellipsoid to
0.05 % in the major axis and 0.08 % in the minor. The figure is **prolate** (tidal), not oblate:
its silhouette is phase-dependent, which an oblate spheroid cannot be. Their θ_LD = 2.495 mas is
the **minor** axis.

Also settled: their map is an intensity map wearing a temperature label (I ∝ T gives r = 0.9909
against 0.9152 for Planck), so their spot deficits are overstated by 2.31× — the real range is
4098–4463 K, not 3742–4589. Every convention needed to use their maps is in
`demos/zetand_surfing_conventions.jl`.

**The limb-darkening disagreement is not a convention error, and it is mostly ours.** alpha
climbs 0.139 → 0.154 → 0.17 → 0.18 → 0.19 as their modelling choices are adopted one at a time,
against 0.22–0.26 from a Kurucz atmosphere and their fixed 0.269. Their value sits where the
physics puts it; a free alpha from an incomplete surface model is the number to distrust, which
is the argument they give for fixing it. The residual 0.19 vs 0.22 is the size of the misfit
still present (chi2/n 2.27 even with their own map).

#### What is actually left on ζ And

1. **The joint model + imaging fit** — still the real answer, still not done. Machinery is
   `joint_reconstruct_oi` (`src/shape_gradient.jl`), which has three documented defects.
   Starting point: the triaxial column above, with T3φ weighted in — it is the spot signal and
   a uniform ellipsoid cannot use it.
2. **Independent λ selection.** The MacKay evidence port gave γ ≈ 1.9e6 against a hard bound of
   768 (CG hitting non-positive curvature) and was never committed; it needs a Gauss-Newton PSD
   operator.
3. **Rerun HMC with σ ≈ 0.04**, not 0.1 — her map's log-contrast sd is 0.04, so the prior was
   2.5× too loose.
4. **V0 (item 18)** was measured here at 1.020–1.022 from a scratch script; making it a fitted
   parameter is what turns that into a result.

---

---

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

---

## What is left, in order

Audited 2026-10-04.

**P1 item 4** — `shell_fit`'s bare `catch → return 1e30` (`src/gui/shell.jl:3341`). The pool
removed the exception this was swallowing, but it will swallow the next one just as quietly.
Minimum viable: refuse a run in which *every* evaluation failed; better, console-log the first
exception and continue. This is the one that makes the next silent failure audible.

**P2 item 7** — spell the Pigeons fallback kernel `:scalar` rather than `nothing`
(`src/gui/shell.jl:3631`), so the console line stays true if someone adds a threading kernel to
`KERNEL_TURBO_METHODS`.

**P4 item 12** — `.github/workflows/` still contains only `docs.yml`. The suite is the thing
most likely to catch a regression and it does not run on a push. Biggest CI gap.

**P4 item 13** — `checkdocs = :none` (`docs/make.jl:10`) still hides the backlog.

**P4 item 3** (one sentence) — the docstring line saying `parallel = :outer` requires a
re-entrant `fitfun` was never added to `bootstrap_parametric`.

**P4 item 15** — the AppImage deploy to `joy:…/ROTIR-x86_64.AppImage` is still unverified.

Then the science: **item 16**, ζ And's joint model + imaging fit. The shape and the conventions
are now settled, so it is no longer blocked on understanding the published model.

**Item 18** (`V0`) is the smallest of the new work and the one with the widest reach — if the
2 % is real, every existing fit has been absorbing it. Do its gradient and null-space test
first; everything else in that item is plumbing.

**Item 17** (chromatic) is gated on a Fisher calculation, not on code. Run the gate before
writing the channel loop.
