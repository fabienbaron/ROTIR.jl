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
