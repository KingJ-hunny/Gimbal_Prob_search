# Gimbal Probabilistic Search — Track & FOU projection (Stage 1)

MATLAB/Octave code for the **satellite laser-ranging gimbal search** problem.
This first stage builds the pieces that everything else will stand on:

1. a **mean-orbit track generator** (LEO / MEO / GEO / high-eccentricity, plus a
   realistic TLE→SGP4 mode),
2. an **along-track-elongated orbit covariance** generator, and
3. a **projection test** (`main_projection_test.m`) that shows how that 3‑D
   covariance, carried through `ECI → ECEF → ENU`, becomes an **elliptical
   angular Field of Uncertainty (FOU)** on the sky, and how its size and
   orientation change along a pass.

The physics follows two references supplied for the project:

* **Yoon, *AE550 HW1 Solution*** — the topocentric az/alt pointing pipeline
  (`LOS in ECI → gimbal/ENU frame → az, alt`).
* **Kim, Yun & Yoon, *Expanding-amplitude Lissajous scanning…*, Optics
  Communications 620 (2026) 133652** — the elliptical FOU model: the
  orbit-prediction error is anisotropic and in-track dominated, and when
  projected onto the plane of sky it subtends an angle `≈ σ_pos / range`,
  giving an ellipse whose major axis is the in-track direction.

---

## Files

| File | Role |
|------|------|
| `generateTrack.m`      | **Function 1.** Mean-orbit track. `mode='tle'` (SGP4) or `mode='kepler'` (two-body + optional J2), any regime. |
| `generateCovariance.m` | **Function 2.** In-track-elongated RIC covariance, rotated to ECI at every sample. |
| `projectFoU.m`         | Projects the ECI covariance onto the az/alt tangent plane → angular FOU (ellipse axes, ratio, tilt). |
| `main_projection_test.m` | Driver/demo. Korea ground station, the 3 supplied TLEs, best-pass search, plots. |
| `parseTLE.m`           | TLE two-line parser → the struct `sgp4.m` expects (+ epoch JD). |
| `sgp4.m`               | Provided SGP4 propagator (unchanged). |
| `findBestPass.m`       | Finds the highest-elevation pass over the station within a search window. |
| `computeAzAlt.m`       | Topocentric az/alt (Yoon pipeline) + `ECI→ENU` rotation for covariance projection. |
| `geodetic2ecef.m`, `eci2ecefMatrix.m`, `gmstRad.m`, `jday.m` | Frame/time helpers. |
| `plotSatelliteFOU.m`, `plotSummary.m` | Figures. |

Run it:

```matlab
main_projection_test        % writes ./figures/*.png and prints a summary
```

Tested end-to-end in GNU Octave 8.4 (headless, gnuplot). All code is written to
run unmodified in MATLAB as well.

---

## What the projection test shows

Ground station: Daejeon / KAIST (36.371 °N, 127.362 °E), fixed in ECEF.
Orbit-error 1‑σ in RIC (default, ≈ Kim et al. Table 2, 24–48 h TLE age):
`σ_R = 0.205`, `σ_I = 0.808`, `σ_C = 0.168` km → `k_orb = 0.208`.

For each of the three supplied satellites the code finds the best pass and
reports how the **3σ angular FOU** evolves:

```
Starlink 35978  range  1754 → 488 → 1736 km   a: 634 → 4962 → 645 urad   k: 0.45 → 0.21 → 0.45
Starlink 31742  range  1738 → 491 → 1649 km   a: 641 → 4940 → 691 urad   k: 0.45 → 0.21 → 0.45
COSMOS 1052     range  3610 → 1473 → 3663 km  a: 421 → 1646 → 423 urad   k: 0.33 → 0.21 → 0.33
```
(values at rise → peak → set)

Three effects, all visible in `figures/`:

* **1/range law.** The FOU grows as the pass climbs (range shrinks) and shrinks
  again on the way down. At culmination the projected major axis matches
  `3 σ_I / range` to < 0.1 %. COSMOS is farther, so its angular FOU is smaller
  for the same physical error — all three fall on one `a` vs `range` hyperbola.
* **Orientation follows the in-track direction.** The ellipse rotates on the sky
  as the ground track curves overhead (see the coloured-by-elevation overlay and
  the tilt-vs-elevation panel).
* **Anisotropy is geometry-dependent.** `k = b/a` → `k_orb ≈ 0.208` near zenith
  (in-track projects fully onto the sky) but is rounder (≈ 0.45) at low
  elevation, where part of the in-track error lies along the line of sight and
  does not move the pointing.

---

## Conventions & modelling choices (so they are easy to change)

* **Frames.** SGP4 output is treated as ECI (TEME); `ECI→ECEF` is a single GMST
  rotation about the pole (IAU‑1982 GMST), neglecting precession/nutation/polar
  motion — consistent with the Yoon derivation and adequate for pointing-scale
  geometry. Swap `eci2ecefMatrix.m` for a full reduction if sub-arcsecond ECEF
  accuracy is ever needed.
* **RIC frame.** `R = r/|r|`, `C = (r×v)/|r×v|`, `I = C×R` (in-track ≈ velocity).
  Rotation to ECI is `Q = [R̂ Î Ĉ]`, `Σ_ECI = Q Σ_RIC Qᵀ`.
* **FOU projection.** Keep only the two components ⟂ to the line of sight
  (radial/range error does not move az/alt), then divide by `range²` for the
  angular covariance. Optional isotropic angular term `σ_iso` (attitude /
  boresight / thermoelastic) can be added in quadrature via `projOpt.sigmaIsoRad`
  (default 0, i.e. pure orbit term, so `k` reflects orbit anisotropy alone).
* **Covariance model.** Static RIC covariance over a single pass (radial→in-track
  drift from differential period is neglected over a few minutes). Pass a full
  `Σ_RIC` (3×3) to `generateCovariance` if you want correlations.
* **Beam / PRF / kinematics** are **not** in this stage yet — this is the FOU
  geometry only. They enter the scan-path design stage.

### Switching orbit regimes / to a new TLE

* New TLE: replace entries in the `sats` cell in `main_projection_test.m`
  (two standard 69-column lines each).
* MEO / GEO / HEO analytic study: set `USE_TLE = false`; edit `keplerSats`
  (classical elements). `generateTrack` labels the regime automatically and can
  add J2 secular rates (`cfg.useJ2`).

---

## Stage 2 — conventional Archimedean-spiral scan time

`scanArchimedeanSpiral.m` + `main_conventional_scan.m` answer: **using the
conventional method (a circular Archimedean spiral centred on the mean
trajectory), how long does it take to scan the whole FOU** under the real
gimbal limits?

Hardware (editable in the main): gimbal `Vaz=20`, `Vel=10` deg/s, `Aaz≥5`,
`Ael≥2` deg/s²; pulse rate `PRF=2 kHz`; beam width selectable from
`{5,10,20,50,100,200}` arcsec.

Model: the spiral reaches `r = a` (FOU semi-major, so the whole ellipse is
enclosed); arm pitch = along-arm pulse spacing `Δs = beamwidth·(1−α)`. The
completion time is a **time-optimal traversal of the fixed spiral path** under
three limits — PRF no-gap (`v_sky ≤ Δs·PRF`), gimbal rate, and gimbal
acceleration — with the azimuth axis stretched by `1/cos(el)` near zenith
(keyhole). Evaluated at culmination (largest FOU = worst case) and swept along
the pass.

Headline results (at culmination, α = 0):

```
                a[urad]  pass[s]   5"     10"    20"    50"   100"   200"   binding
Starlink 35978    4962     480    541*   271    135     54     27     13    az accel (keyhole)
Starlink 31742    4940     468    444    222    111     44     22     11    az accel (keyhole)
COSMOS 1052       1645    1068     88     44     22      9      4      2     az accel (keyhole)
   (* 5" beam on Starlink 35978 exceeds its own pass duration -> cannot finish in one pass)
```

Three findings, all in `figures/scan_*`:

* **Scan time ∝ 1/beamwidth²** (finer beam ⇒ tighter arms *and* a lower
  PRF-limited speed): a 5″ beam is ~40× slower than 200″.
* **The azimuth keyhole dominates near zenith.** For every case the binding
  constraint is azimuth acceleration: at el≈84° the az axis is stretched ×9.6,
  so the effective on-sky az accel is only `Aaz·cos(el) ≈ 0.5 deg/s²`, and the
  spiral's twice-per-turn az reversals throttle the whole scan (`T_full` is
  ~10× the velocity/PRF-only bound). Scan time climbs steeply with elevation.
* **The circular spiral massively over-scans the thin ellipse** (see
  `scan_spiral_demo.png`): a disk of radius `a` for an FOU whose minor axis is
  `k≈0.2` of that — the exact inefficiency the elliptical-FOU scan patterns
  (Kim et al.) were designed to remove, and the motivation for the smart method.

Run it:

```matlab
main_conventional_scan      % prints the tables, writes figures/scan_*.png
```

## Stage 3 — spiral-vs-FOU race (time-varying FOU)

`spiralExitTime.m` + `main_spiral1_test.m` relax the frozen-FOU assumption of
Stage 2. Acquisition begins at first visibility (**elevation 10°, rising**) and
an Archimedean spiral is started from the FOU centre, winding **outward** at the
gimbal/PRF-limited rate. At every epoch the spiral radius `r_spiral(t)` is
compared with the current 3σ FOU semi-axes `a(t)`, `b(t)`.

Key physical point (corrects the "FOU shrinks" intuition): from 10° on the way
up the angular FOU **grows** (1/range) to a maximum at culmination, then shrinks.
So an outward spiral and the FOU can cross **more than once**. Two times are
reported per beam width:

* `t_first` — spiral radius first reaches `a(t)` (first encloses the small
  initial FOU), and
* `t_persist` — the epoch after which `r_spiral(t) ≥ a(t)` for the **rest of the
  pass** (the FOU never overtakes the spiral again). This is the meaningful
  "FOU permanently covered" time.

Result (start elev 10°, α = 0):

```
                       t_first[s]   t_persist[s]   elev@persist
Starlink 35978  5"        17           259            63°   <- FOU outruns the slow spiral to culmination
                10"        8             8            11°
               20-200"    1             1           ~10°
COSMOS 1052    5"         7             7            10°
```

Reading `spiral1_race.png`: the FOU `a(t)` is the red hump (683 → 4962 → 645
µrad). Beams ≥ 10″ expand fast enough to clear the hump in the first seconds and
stay outside for the whole pass. The **5″ spiral is too slow**: it clears the
tiny initial FOU at 17 s, but the ballooning FOU climbs back **above** it from
~30 s to 259 s — i.e. for ~4 minutes around culmination the FOU rim sits
*outside* the 5″ spiral (uncovered) — and only on the descent, as the FOU
shrinks, does the spiral permanently enclose it (dot at 259 s, elev 63°).

So with the finest beam a single outward spiral cannot keep up with the FOU
growth near zenith; a coarser beam, or a scan that adapts to the growing FOU, is
required. This is the frozen-FOU caveat of Stage 2 made quantitative, and points
straight at the time-coupled / FOU-adaptive scan design.

Model note: `spiralExitTime.m` integrates `r_spiral(t)` with quasi-steady speed
caps (PRF no-gap, axis velocity, and centripetal acceleration, all with the
`1/cos(el)` azimuth keyhole). It is a geometric race (spiral radius vs FOU
radius), independent of the Stage-2 time-optimal traversal, and slightly
conservative.

Run it:

```matlab
main_spiral1_test           % prints the tables, writes figures/spiral1_*.png
```

## Stage 4 — Monte-Carlo acquisition of moving points (scan failure)

`spiralAcqSim.m` + `eciToAzElVec.m` + `main_points_acq_test.m` replace the
analytic FOU with **N = 3000 candidate points** drawn from the in-track
covariance. Each point is a fixed RIC offset (a hypothesis) propagated through
`ECI→ENU`, so the cloud translates, rotates and breathes exactly like the
analytic FOU — but now every point has its own moving sky track.

Above a detectable elevation (default 20°) the gimbal scans the mean path with
an **overlapping Archimedean spiral**: point `k` at `θ_k = √(4πk)`,
`r_k = p·√(k/π)`, pitch `p = α·Rbeam` (α<2 ⇒ footprint overlap — the requested
`θ_k ∝ √k` construction). Point `k` is visited at time `t_k` from traversing the
spiral under the gimbal/PRF limits (velocity, PRF no-gap, centripetal
acceleration, all with the `1/cos(el)` keyhole). A candidate is **acquired** the
first time a scan point lands within `Rbeam` of the candidate's position **at that
point's visit time** — so a point is missed if the moving cloud is never
coincident with the beam while the beam is there.

Result (N=3000, detElev 20°, α=1):

```
                     5"      10"    20"   50"  100"  200"
Starlink 35978     92.8%   99.2%  100%  100%  100%  100%
Starlink 31742     92.5%   99.2%  100%  100%  100%  100%
COSMOS 1052        100%    100%   100%  100%  100%  100%
```

The scan **fails** (acq < 100%) for the fine beams on the fast Starlink passes:

* `acq_sky_snapshot.png` — the missed points are the **extreme cross-elevation
  tips**, i.e. the far along-track (timing) outliers of the elongated FOU. Near
  culmination these reach the largest angular excursion and move fastest, and
  the slow fine-beam spiral never reaches them in time.
* `acq_cumulative.png` — the failure mechanism in time: the 5″ curve climbs to
  ~87 % in the first ~50 s (low elevation, compact FOU), then **plateaus through
  the whole high-elevation phase** — the ballooning FOU has outrun the spiral, so
  no new points are acquired around culmination — and only resumes on the descent
  as the FOU shrinks, ending at ~93 % when the pass runs out.
* COSMOS (longer pass, smaller angular FOU) keeps up at every beam.

So a single conventional overlapping spiral, with a fine beam and the real gimbal
limits, cannot guarantee acquisition of a **moving** FOU on a fast overhead pass.
The result is sensitive to the detectable-elevation floor (a lower floor gives the
spiral more low-elevation time to sweep the compact FOU early). This is the
quantitative "scan can fail as the FOU moves" demonstration and the motivation for
a motion-aware / FOU-adaptive scan.

Run it:

```matlab
main_points_acq_test        % prints acquisition tables, writes figures/acq_*.png
```

## Stage 4b — constant-linear-velocity spiral, overlap factor, N = 40000

The spiral scan-point allocation was set to the constant-linear-velocity
Archimedean rule the user specified:

```
r(theta) = theta * R_s/(2*pi),   R_s = 2*alpha*Rbeam,   theta_k = sqrt(4*pi*k)
=> r_k = R_s*sqrt(k/pi)   (one turn advances r by exactly R_s; arc spacing = R_s)
```

`alpha` is the overlap factor: `alpha = 1/sqrt(2)` → `R_s = sqrt(2)·Rbeam` (the
critical square-packing spacing, full coverage); `alpha = 1` → `R_s = 2·Rbeam`
(footprints tangent, leaving interstitial gaps). `main_points_acq_test.m` now
runs **N = 40000** points, both `alpha`, all beam widths, and reports overall
**and rim** (outer 10 %) acquisition.

```
Starlink 35978        alpha = 1/sqrt(2)          alpha = 1
   beam"          overall     rim          overall     rim
     5             96.1%     63.7%           81.4%     83.3%
    10             99.8%     98.5%           80.3%     85.8%
    20             99.7%    100.0%           78.6%     83.0%
    50             98.3%    100.0%           72.5%     80.9%
   100             94.2%    100.0%           58.9%     78.8%
   200             84.3%     98.9%           40.1%     66.0%
```

This resolves the earlier "why is it >90 %?" puzzle — that run used a much
tighter spacing. With the correct `R_s = 2·alpha·Rbeam` the conventional scan
misses a great deal, and the pattern is instructive:

* **`alpha = 1` (tangent) fails badly — 40–80 %.** The missed points sit in the
  DENSE CENTRE of the cloud (`acq_alpha_snapshot.png`, right): a tangent spiral
  leaves interstitial gaps everywhere, and near the centre the arms are angularly
  sparse, so the most probable region is the worst covered. rim > overall
  throughout (the Archimedean spiral under-samples its own centre).
* **`alpha = 1/sqrt(2)` is near-complete only in a sweet spot (~10–20").** It
  fails at both ends for different reasons: the **fine beam** (5") reaches only
  ~64 % of the RIM — too slow to sweep out to the tips within the pass (the
  kinematic/time limit of Stages 2–4); the **coarse beam** (100–200") drops to
  ~84 % overall because the discrete dwell points get sparse near the centre.
* The dense-centre / sparse-arm effect means the **overall** rate can sit *below*
  the **rim** rate — opposite to the naive expectation.

Caveat (model choice): acquisition is tested at discrete dwell points. A truly
continuous sweep would close the along-path and centre gaps for
`alpha = 1/sqrt(2)` (whose radial pitch already overlaps), so the coarse-beam
drop is partly a dwell-point artifact; the `alpha = 1` radial gaps and the
fine-beam rim-timing losses are real either way.

## Stage 4c — uniform coverage, detectable elevation, and PRF

Two additions to the point-acquisition study.

**Uniform distribution (coverage measure).** `main_points_acq_test.m` now also
draws the N = 40000 points **uniformly inside the 3σ RIC ellipsoid** alongside
the Gaussian set. With equal weight per hypothesis the acquisition rate is the
**covered fraction of the 3σ region** — a pure coverage measure with none of the
Gaussian centre-weighting. It also sweeps the **detectable elevation** at which
the scan starts (10° and 20°). `acq_dist_det_vs_beam.png`: rows = gaussian /
uniform, colour = α, line style = detElev.

```
Starlink 35978, alpha = 1/sqrt(2)      gaussian            uniform (coverage)
   beam"                              det10   det20       det10   det20
     5                               100.0%   95.9%      100.0%   87.4%
    10                                99.8%   99.8%       99.9%  100.0%
    50                                96.2%   98.3%       98.7%   99.4%
   200                                79.6%   84.3%       82.8%   91.7%
```

* **Uniform reads lower than Gaussian at the edges** (e.g. 5", det20: 95.9 % vs
  87.4 %) — the coverage measure exposes the unswept rim that the probability
  weight hides.
* **A lower detectable elevation helps the fine beams**: starting at 10° gives
  the spiral more low-elevation time to sweep the compact FOU before it balloons,
  so 5" rises from 87 % (det20) to 100 % (det10). Coarse beams, limited by
  spatial gaps rather than time, barely move with detElev.

**Beam width vs PRF (`main_beam_frequency_test.m`).** With α = 1/sqrt(2) fixed,
one detectable elevation, and the uniform cloud, it sweeps PRF ∈ {100, 500,
1000, 2000, 3000} Hz across beam widths to explain *why* the error grows with
beam width. PRF enters only through the speed cap `v ≤ R_s·PRF` (one pulse per
point spacing), so it changes scan *speed*, not spatial spacing. The result
decomposes the beam-width error into two regimes (`beam_freq_coverage.png`,
`beam_freq_vs_prf.png`):

```
Starlink 35978 coverage %      100Hz  500Hz  1kHz  2kHz  3kHz
   5"   (fine)                  61.0   87.2  87.4  87.4  87.4   PRF-sensitive, then saturates
  10"                           93.5  100.0 100.0 100.0 100.0
  20"+  (coarse)               ~const across PRF (99.9 / 99.4 / 97.8 / 91.7)
```

* **Fine beams are time-limited** — the spiral is slow and needs a high pulse
  rate; coverage climbs steeply from 100→500 Hz, then **saturates** (at 5″ it
  plateaus at ~87 %, now capped by the gimbal velocity/acceleration keyhole, not
  PRF).
* **Coarse beams are spatially limited** — every PRF curve lands on the same
  point; the drop toward 200″ is the discrete-dwell interstitial gap, which no
  amount of PRF fixes.
* At the operating 2 kHz we are already in the saturated regime; raising PRF
  further buys nothing, and the coarse-beam loss is fundamental to the spacing.

## Stage 4d — low vs high culmination (is it the high elevation?)

`main_low_elev_test.m` tests whether the near-zenith culmination is what causes
the fine-beam misses, by comparing two passes over Korea under identical settings
(α = 1/sqrt(2), PRF = 2 kHz, N = 40000, uniform + Gaussian):

* a **synthetic LEO tuned to culminate ~40°** — same Starlink-like `a, i` (so the
  FOU physics is identical), only the geometry differs. Found by tuning RAAN:
  `a = 6861.2 km, e = 0.001, i = 43°, RAAN = 284.5°, argp = 0, M0 = 0` →
  best pass peaks **40.7°** (range 724 km, a_max 3350 µrad, keyhole 1/cos = 1.3);
* **Starlink 35978** — best pass peaks **83.9°** (range 488 km, a_max 4962 µrad,
  keyhole 1/cos = 9.5).

The result depends on where the scan *starts* (`low_vs_high_elev.png`):

```
uniform coverage %      start elev 10 deg        start elev 20 deg
  beam"               40 deg    84 deg          40 deg    84 deg
    5                 100.0     100.0            94.8      87.4
   10                  99.9      99.9           100.0     100.0
  200                  85.4      82.8            93.8      91.7
```

* **Start at first visibility (10°):** the 40° and 84° curves nearly coincide —
  the near-zenith keyhole does **not** cause a fine-beam failure, because the
  spiral sweeps the compact low-elevation FOU and catches almost everything
  before reaching the keyhole region.
* **Start higher (20°):** now the high pass IS worse at the finest beam — 5″
  drops to **87.4 %** at 84° vs **94.8 %** at 40° — the near-zenith keyhole plus
  the larger 1/range FOU cost ~7 pp. So the "high-elevation 찐빠" is real, but
  only when you cannot use the early low-elevation window.
* **Coarse beams (200″)** are roughly elevation-independent (both ~83–94 %): that
  loss is the discrete-dwell spatial gap, not the elevation.

Bottom line: the high culmination hurts the finest beam **only if the scan starts
high**; starting at first visibility removes it. The coarse-beam gap is not an
elevation effect at all.

## Stage 6 — conventional spiral, control-frequency sweep, harvest probability

`main_sept_spiral_test.m` (+ `buildSpiralPoints.m`, `harvestScan.m`,
`precomputeOffsets.m`) is a particle-based conceptual study of a **conventional
Archimedean spiral with overlap = 1** (`alpha = 1`, `R_s = 2·Rbeam`, footprints
exactly tangent → interstitial gaps remain). The spiral is drawn for the
**largest (zenith) FOU** out to its major axis (K points). Scan point k is
visited at `t_k = t_start + (k-1)/f` — pure **beam-control-frequency** timing,
no gimbal dynamics — and f is swept over {100, 250, 500, 1000, 2000} Hz.
Harvest = fraction of N = 40000 Gaussian particles covered.

Starlink 35978, beam 50″, peak 83.9°, `a_zenith = 4962 µrad`, `R_s = 242 µrad`,
`K = 1316`:

```
(A) STATIONARY (frozen zenith FOU):  harvest = 77.7%  (frequency-independent)
    scan time K/f:  100Hz->13.2s  250->5.3  500->2.6  1k->1.3  2k->0.66s

(B) ORBITING (same scan time K/f, spiral tracks the mean orbit):
     f[Hz]   harvest@10deg   harvest@zenith
     100         64.9%          77.7%
     ...          (flat)         (flat)
     2000        64.9%          77.8%
```

Findings:

* **The conventional overlap-1 spiral misses ~22% even on a FIXED FOU** (77.7%):
  the tangent spacing leaves interstitial gaps, plus a few % of Gaussian mass
  lies beyond the 3σ major axis the spiral stops at.
* **Scan time = K/f exactly** (0.66–13 s), so a higher control frequency only
  makes the same scan faster — it does not change the stationary harvest.
* **Harvest is essentially frequency-independent, and start-at-zenith ≈
  stationary** — with control-frequency timing even the slowest scan (13 s) is
  fast compared with how quickly the FOU evolves (~minutes), so the FOU's motion
  over one scan is negligible (its rotation near zenith even helps slightly, by
  letting particles drift out of the fixed gaps).
* **Start-at-10° is lower (64.9%)**, but this is a *size-mismatch*, not a motion
  effect: the small low-elevation FOU sits in the sparse inner turns of the
  zenith-sized spiral (where the α=1 gaps are widest), so it is harvested worse.

Caveat: because timing here is `t_k = k/f` (the requested control-frequency
model), the beam is assumed to slew as fast as f demands, ignoring the gimbal
vmax/amax. Under the *real* kinematic scan speed the scan would take 10–100×
longer and the FOU motion WOULD bite — that is the regime of Stages 3–4.

## Stage 5 — control: FOU frozen onto the trajectory (motion vs size)

`main_points_acq_fixed_test.m` isolates the cause of the Stage-4 failure. It
freezes the FOU rigidly onto the mean trajectory — every candidate keeps a
**constant** angular offset from the mean pointing (no 1/range breathing, no
rotation, no differential motion), the cloud only translating with the mean —
and runs the same overlapping spiral. Two freeze sizes bracket it: at the small
scan-start projection and at the large culmination projection.

Result (Starlink 35978, N=3000, detElev 20°, α=1):

```
 beam   moving   fixed@culm   fixed@start
  5"    92.8%      92.0%        100.0%
 10"    99.2%      99.5%        100.0%
 >=20"  100%       100%         100%
```

The finding refines Stage 4: **the failure is a size-vs-scan-budget effect, not
the motion itself.**

* Freeze the FOU at its **small** start size → **100 % at every beam**: a compact
  static FOU is fully covered even by the 5″ beam.
* Freeze it at its **large** culmination size → **~92 % at 5″**, essentially
  identical to the moving case. A large region simply cannot be swept by a fine
  beam within the pass, whether it moves or not.
* Moving vs fixed@culmination differ by <1 pp (moving is even marginally better,
  because early on the real FOU is compact and the spiral catches some tips while
  they are close in).

So what breaks the fine-beam scan is that the FOU **grows** toward zenith (1/range)
until it exceeds the beam's coverage budget — the rigid translation of the cloud
is almost irrelevant. This points the fix at the coverage budget (coarser beam,
overlap, or a size-adaptive scan), not merely at "tracking the motion."

Run it:

```matlab
main_points_acq_fixed_test
```

## Where this sits in the larger problem

The projected FOU here is exactly the **moving, rotating, elongated uncertainty
ribbon** the gimbal must sweep. Its parametrisation by the along-track (timing)
axis — the ellipse major axis — is what makes the search-path design tractable:
the candidates line up along the sky track, so the scan reduces to a smooth
lead/lag schedule over that axis plus a thin cross-track dither. Stage 2 will add
the beam footprint, the PRF, and the gimbal velocity/acceleration limits on top
of the geometry established here.
