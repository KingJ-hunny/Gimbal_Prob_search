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

## Where this sits in the larger problem

The projected FOU here is exactly the **moving, rotating, elongated uncertainty
ribbon** the gimbal must sweep. Its parametrisation by the along-track (timing)
axis — the ellipse major axis — is what makes the search-path design tractable:
the candidates line up along the sky track, so the scan reduces to a smooth
lead/lag schedule over that axis plus a thin cross-track dither. Stage 2 will add
the beam footprint, the PRF, and the gimbal velocity/acceleration limits on top
of the geometry established here.
