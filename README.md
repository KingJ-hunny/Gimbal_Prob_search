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

## Where this sits in the larger problem

The projected FOU here is exactly the **moving, rotating, elongated uncertainty
ribbon** the gimbal must sweep. Its parametrisation by the along-track (timing)
axis — the ellipse major axis — is what makes the search-path design tractable:
the candidates line up along the sky track, so the scan reduces to a smooth
lead/lag schedule over that axis plus a thin cross-track dither. Stage 2 will add
the beam footprint, the PRF, and the gimbal velocity/acceleration limits on top
of the geometry established here.
