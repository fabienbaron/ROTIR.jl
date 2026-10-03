# ROTIR — outstanding fixes

Compiled 2026-10-03. Everything here is backed by a measurement or a verified code reference;
where something is inferred rather than observed it says so.

Ordered by **how the failure presents**, not by effort: a wrong number that looks right is
worse than a crash, and a crash is worse than an ugly message.

---

## Context: what was fixed and verified on 2026-10-03

No action needed on these; listed so the items below are not re-derived.

- `_fit_nested` silently dropped `threaded`; `fit_sphere_ld`, `fit_ellipsoid_ld` and
  `fit_orbit` each dropped it on the way in. All now accept and forward it, default `false`.
- The GUI's `_run_fit` nested-sampler branch passes `threaded = false`.
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

A user gets a number that is wrong and looks fine. Fix these first; they turn every other
item on this list from "silently wrong" into "loudly wrong".

### 1. `:t3` re-entrancy guard

**The bug.** `src/type3_nufft.jl` is the process-wide default forward kernel and is not
re-entrant on two independent counts:

- `type3_plan_for` (~`:784`) returns a **globally cached** `Type3Plan`; `TYPE3_PLANS_LOCK`
  guards construction, not use. Concurrent calls on the same data share one plan object.
- `Type3Plan` (`:91`) carries mutable scratch — `grids`, `scratch`, `centred`, `realbuf` —
  and `type3_quads!` (~`:455`) opens by zeroing it.

(Note: the `@threads :static` schedule at `:409`, `:458`, `:533`, `:612` is *not* load-bearing —
see P3. It is a separate, global constraint, not the reason the scratch is unsafe.)

**Measured, and worse than a throw.** Two tasks on one plan, 2000 calls each, three trials:
2000 throws, 1999 correct, and **1 silently wrong** each time, worst relative error 1.30–4.34.
The window is structural: `type3_quads!` zeroes serially and `_t3_finish!` reduces/FFTs
serially before its `:static` gather, and in those stretches no threaded region is open, so
Julia's own precondition passes and the second entrant walks onto the first's scratch. A
finite, wrong visibility — which an objective with `catch → return 1e30` *accepts*.

**Fix.** `jl_in_threaded_region` check + atomic CAS, on `type3_quads!`, `type3_quads_adj!`
and `type3_points!` — the three functions that touch the scratch. **Process-wide, not
per-plan**: Julia's check is global, and `type3_plan_for` caps the cache with a wholesale
`empty!` at 32 entries, so two tasks can legitimately hold different plan objects for the same
key and a per-plan flag waves both through (demonstrated).

**Cost.** 4.3 ns per call against a 0.37 ms cheapest real call — below the noise floor of any
A/B on this kernel (the identical no-op control scattered ±4.5 %).

Do **not** make it serialise by default: that converts a caller-side design error into an
unexplained wall-clock cliff. Offer serialisation as an opt-in; throw by default.

### 2. Make the guard a typed exception

`T3NotReentrant <: Exception`, not a bare `ErrorException`. Three places swallow exceptions
today and would swallow the guard too:

- OITOOLS `bootstrap_driver` → `run_chunk` — comments the row to NaN
- `src/gui/shell.jl:3576-3580` — `catch → return 1e30`
- `model_state` — `catch err → console!("χ² failed: …")`

Without the type, item 3 stays silent even with the guard in place.

### 3. `bootstrap_parametric` at its default `parallel = :outer`

With any `:t3`-route `fitfun`, OITOOLS runs replicates under `Threads.@spawn` and `run_chunk`
swallows the throw. Result: every replicate NaNs out, and you get either
`"fewer than 2 usable replicates"` or a σ computed from whichever handful happened not to
collide. **This is the default setting and the failure is quiet.**

`parallel = :inner` is serial and safe. The convenience overload
`bootstrap_parametric(DATA, TESSELS, TEPOCHS, BASE; …)` is also safe — its `fitfun` is
`fit_parametric`, which routes through `interferometric_chi2`, a different kernel.

### 4. `shell_fit`'s bare `catch → return 1e30`

This is what turned a hard failure into a plausible fit with an evidence: 40 000 of 40 000
likelihood calls hit the sentinel, and the panel reported `log(Z) = -5e29` (finite),
`N_eff = 1`, the flat prior as the posterior, and a median 1.402 against the true 1.37717.

A derivative-free objective genuinely needs to survive geometry failures, so this is a design
call, not a mechanical fix. Minimum viable: refuse a run in which *every* evaluation failed.
Better: console-log the first exception, then continue.

### 5. `type3_quads_adj!`

Same `:static` pattern over the same shared scratch (~`:612`). Never tested. Covered by item 1
if the guard is placed there, which it should be.

---

## P2 — loud failure, or claims in our own source that are wrong

### 6. The timing comment at `src/gui/shell.jl:4152-4153` is wrong

It reads *"Serial `:t3` is also the faster of the two repairs: 0.20 ms per evaluation at
HEALPix 3"*. **0.20 ms is the 16-thread figure; serial is 0.94 ms.** Measured speedup from
`:t3`'s internal threading, serial ÷ threaded at equal thread count:

| nside | tessels | 4 thr | 8 thr | 16 thr |
|---|---|---|---|---|
| 2 | 192 | 2.80× | 3.78× | 3.88× |
| 3 | 768 | 3.14× | 5.02× | 5.07× |
| 5 | 12 288 | 3.85× | 6.86× | 7.50× |
| 6 | 49 152 | 3.67× | 6.57× | 8.54× |

It earns its keep at every mesh size tested. The `threaded = false` fix is still correct — it
serialises the *sampler's outer loop* while the kernel keeps its inner threading, which is
exactly the right split — but the stated justification must go. It has already propagated once:
a later investigation cited this comment as evidence.

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

### 8. The GUI thread is not quiesced during a fit

`start_job!` refuses a second job, so job-vs-job is impossible. But only **6** file-IO buttons
carry `enabled: !win.jobRunning` in `src/gui/qml/Main.qml`; the parameter fields, epoch stepper,
tab bar and kernel combo stay live. `shell_set_param` clears `sh.chi2key[]`, and the next
refresh runs `epoch_chi2` → `model_state` → `observables` → `fused_cvis` → `:t3` **on the GUI
thread** while the worker is also on `:t3` (for `:neldermead`/`:bobyqa`/`:nautilus` and any
binary fit). The worker's `catch → 1e30` then silently takes a sentinel at that trial point.

Fix: extend `enabled: !win.jobRunning` to the parameter form and epoch stepper, or early-bail
in `model_state` when `sh.job !== nothing`. Source-verified; not reproduced in a live GUI.

### 9. `spheroid_chi2_allepochs_f` threads over epochs into `:t3`

`src/oichi2_spheroid.jl:409` — `Threads.@threads for i=1:nepochs` → `spheroid_chi2_f` →
`observables`, which takes the matrix-free `:t3` branch whenever `setup_oi!` has not been run.
Exported. Latent only because both in-repo callers happen to run `setup_oi!` first, which the
function cannot know. Fix: route the epoch loop serially when `isempty(stars[i].polyft)` — the
dense route keeps its threading and loses nothing.

---

## P3 — real re-entrancy

**This is measured, prototyped, bit-identical and FASTER than what ships. It is not a
trade-off.** Two beliefs stated earlier in this file's history were wrong and are corrected
here, because both shaped the design:

- **`:static` is NOT required for `grids[tid]`.** There is no `Threads.threadid()` anywhere in
  `src/type3_nufft.jl` (verified: `grep threadid` returns nothing). In
  `for tid in 1:Threads.nthreads()`, `tid` is the **loop index**, and both `p.grids[tid]` and the
  `lo`/`hi` quad block are keyed on it, so each iteration owns its grid whichever OS thread runs
  it. `:dynamic` is bit-identical. `:static` buys only first-touch locality, which measures as
  nothing (nf = 72 at every HEALPix level 2–6; grids are 41 KiB and live in L2).
- **The `:static` throw is not about shared scratch.** `jl_in_threaded_region` is a *global*
  counter, so two concurrent `:t3` transforms on **different plans**, sharing nothing, also
  throw. The scheduler change is therefore needed independently of the scratch fix.

What does hold is the warning that they must land together: shared scratch with `:dynamic`,
8 concurrent callers → **0/8 correct, worst relative error 1.03e+28**.

### 10. Scratch pool on the plan + `:dynamic` — the recommended fix

Keep `tbase`/`twx`/`twy`/`phihat`/`kcoef`/`psicoef`/`fft!` shared read-only in the cached plan;
move `grids`/`scratch`/`centred`/`realbuf` into a pool field, checked out per call in a
`try`/`finally`. Four `:static` → `:dynamic`. No public signature changes.

Per-call cost, single caller, min of 11 interleaved rounds:

| nside | npix | shipped `:static` | pooled `:dynamic` |
|---|---|---|---|
| 2 | 108 | 0.158 ms | **0.137 (−13.1 %)** |
| 3 | 422 | 0.247 | **0.206 (−16.6 %)** |
| 4 | 1676 | 0.437 | **0.395 (−9.6 %)** |
| 5 | 6714 | 0.843 | **0.798 (−5.3 %)** |
| 6 | 26822 | 2.587 | **2.621 (+1.3 %)** |

Checkout/checkin is 374 ns on a 206 µs call (0.18 %) with no measurable contention — 16 tasks on
one pool run at 89 ns/pair. The pool grows to the number of *simultaneous* callers (5 sets for 8
tasks), 853 KiB each, mesh-independent. It lives on the plan so it dies with it when
`TYPE3_PLANS` hits its 32-entry cap; a free-standing pool keyed on the plan would leak there.

Concurrency restored **with scaling**, nside 3: 3547 → 6393 → 7888 → **9982** calls/s for 1, 2,
4, 8 callers (2.8×), and **10 819** with `nchunks = 4` instead of 16 — a knob the design opens
up, since with private scratch the chunk count need no longer equal `nthreads()`.

Bit-identical to the shipped kernel — not 1e-12, byte-for-byte — across fresh/pooled/held/locked
scratch × `:static`/`:dynamic`/`:spawn`/serial, at nside 3 and 5, forward **and adjoint**.

### 11. The adjoint takes the same change

`type3_quads_adj!` shares `centred`/`scratch`/`realbuf` and has the same loop. Prototyped:
bit-identical, 8/8 concurrent callers correct, and **10.7 % faster** at nside 4 (0.645 → 0.576
ms). The shipped adjoint throws under 8 callers.

### 12. Rejected alternatives, with the reason

- **Fresh scratch per call** — simplest, correct, but **+24 to +28 %** at HEALPix 2–4 (+5–7 % at
  5–6). The scratch is 853 KiB and mesh-independent, so the cost is flat while the kernel is not,
  hurting most where `:t3` is cheapest. Sustained cost is 27.6 µs per set, 84 % of it GC. Fallback
  only if the pool is judged too much machinery.
- **Lock the plan** — bit-identical, ~3 lines, and it legitimately silences the `:static` error.
  But throughput *degrades*: 3834 → 3819 → 3405 → **3084** calls/s for 1, 2, 4, 8 callers, i.e.
  8 concurrent chains run 20 % slower in aggregate than one. Right emergency patch, wrong
  destination.
- **Task-local storage** — lookup is 358 ns (no win over the pool) and the lifetime is wrong both
  ways: `@spawn`-per-call samplers get a fresh task each time, degenerating to per-call
  allocation; long-lived workers accumulate 853 KiB per plan touched with nothing to free them.
- **Per-thread scratch keyed beyond `threadid`** — that *is* the pool with a fixed backing store,
  and `threadid()` is unsound under task migration on Julia ≥ 1.7 anyway.

### 13. Incidental, worth keeping

- **`:dynamic` is free-to-faster, never slower**: −11.7 %, −16.5 %, −7.3 %, −6.6 %, +2.7 % at
  nside 2–6. `:static` pins iteration *i* to thread *i*, so one contended thread stalls the whole
  barrier — exactly the condition a multi-chain sampler creates.
- **No `FFTW.UNALIGNED` needed.** `plan_fft!` records the buffer's alignment class, which looked
  like a blocker for per-call buffers; measured class 0 for **400/400** fresh 72×72 `ComplexF64`
  arrays, so the cached plan applies to any per-call Array. Add a one-line `FFTW.alignment_of`
  assertion at scratch construction. `UNALIGNED` would cost +6 µs per transform.

### 14. What this unblocks

The three shipped `threaded = false` defaults exist solely because of this bug —
`src/parametric_fit.jl:273` and `:377`, `src/gui/shell.jl:4154`. The pool is what lets them go
back to `true`, worth **2.8–3.1× at 8 chains**, against the 0.80× the lock would give.

### 12. `binary_cvis`'s matrix-free branch

`src/oichi2_binary.jl:168-170`, public, on the `:t3` route with the same `isempty(polyft)`
auto-switch. No threaded caller in-repo today. Covered by the guard; listed so it is not
forgotten if someone parallelises over epochs.

---

## P4 — hygiene

### 13. CI has no test workflow

`.github/workflows/` contains only `docs.yml`. **`test/runtests.jl` has never run on push or
PR** — 2892 assertions that execute only locally. The main suite needs no display and would run
on `ubuntu-latest` as-is. This is the single biggest CI gap.

### 14. Documenter's `checkdocs = :none` hides a backlog

Undocumented or unreachable exports surface only when something `@ref`s them, at which point
the build fails hard — this broke the docs build three times in one session. Known:
`fused_cvis`, `binary_cvis`, `poly_to_cvis`, `limb_mu` are exported but in no `@docs` block
(`api/fused_polyft.md` documents by markdown table, which is the structural cause);
`provider_intensity` is exported with no docstring at all, and it is the single call the whole
provider interface exists to offer. Consider `checkdocs = :exports`, expecting a backlog.

### 15. The parse-scan test is syntax-only

`test/runtests.jl`'s "every shipped source file parses" uses `Meta.parseall`, which validates
syntax but does **not** lower. It would have caught the `fit_pigeons.jl` unterminated string
that motivated it, but it misses `const` on a local, undefined variables and bad signatures —
it did in fact miss exactly that in `test_nested_threading.jl`. Say so in its docstring so it is
not read as a general validity check.

### 16. Uncommitted work and the AppImage

The AppImage is built (562 600 440 B) and smoke-tested but **not deployed**, and corresponds to
no commit. The deploy target is
`joy:/nfs/morgan/users/baron/public_html/ROTIR-x86_64.AppImage` — a public URL; upload to a
`.new` name, verify SHA-256, then move into place.

### 17. ζ And

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

## If you only do three

**1 + 2 together** (the guard, typed), then **4**. Those close the corruption window and make
every remaining failure audible, for 4.3 ns.

Then **P3 item 10** — and note it is not a trade-off to weigh: the pooled `:dynamic` kernel is
bit-identical, 5–17 % *faster* single-threaded at HEALPix 2–5, and 2.8–3.1× faster at 8
concurrent callers. It supersedes the guard rather than competing with it, and it is what lets
the three `threaded = false` defaults go back to `true`.
