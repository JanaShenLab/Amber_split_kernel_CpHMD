# CpHMD and AREX input reference

Input options for continuous constant-pH molecular dynamics (CpHMD) and native
asynchronous pH replica exchange (AREX) in Amber24 `pmemd`, accompanying:

**Dayhoff, G. W. and Shen, J. (2026).** [*High-Performance All-Atom
Continuous Constant-pH Molecular Dynamics and Native Asynchronous Replica
Exchange in Amber*](https://doi.org/10.64898/2026.11.11.766022).
*bioRxiv*.

## Input files

| File | Namelist | Purpose |
|---|---|---|
| `mdin` | `&cntrl` | CpHMD and AREX run controls. |
| `phmdin` | `&phmdin` | Titration dynamics, residue selection, and output. |
| `phmdparm` | `&phmdparm` | Titration parameters for residue types. |
| CpHMD restart | `&phmdrestrt` | Titration coordinates and velocities. |

Pass the inputs with `-i`, `-phmdin`, and `-phmdparm`. The lambda output file
is selected with `-phmdout`. CpHMD restart input/output uses the associated
PHMD restart command-line options.

**Status key:** **Required** = needed to use the named feature;
**Conditional** = needed only with the specified option;
**Optional** = has a default and may be omitted. **Unset** means no usable
value is supplied by default.

## `mdin` (`&cntrl`)

### CpHMD controls

| Keyword | Status | Default | Description |
|---|---|---|---|
| `iphmd` | Required for CpHMD | `0` | CpHMD mode: `1` (GB), `3` (PME); `0` disables CpHMD. |
| `solvph` | Optional | `7.0` | Solution pH for a fixed-pH run; AREX uses its active ladder value. |

Example: PME-CpHMD at pH 7.0 (include in a complete MD input):

```fortran
&cntrl
  iphmd = 3,
  solvph = 7.0,
  ! Other MD settings for the chosen simulation protocol
/
```

### Native AREX controls

| Keyword | Status | Default | Description |
|---|---|---|---|
| `arex` | Optional | `0` | `1` enables AREX; `0` disables it. |
| `arex_type` | Required for AREX | `''` | Exchange type; currently `'phrex'`. |
| `arex_ladder` | Required for AREX | Unset | Ordered pH conditions; at least two contiguous entries (maximum 100). |
| `arex_segment_steps` | Required for AREX | `0` | MD steps per segment; must be positive. |
| `arex_warmup_segments` | Optional | `0` | Additional warmup segments per pH slot before exchange; not included in `nstlim`. |
| `arex_chkfrq` | Optional | `0` | Completed AREX segments between checkpoints; `0` = final only. |
| `arex_seed` | Optional | `0` | AREX random seed: `0` uses default RNG state, positive fixes the seed, negative uses a clock seed. |
| `arex_num_workers` | Optional | `0` | Number of AREX worker processes. |
| `arex_emulate_restart` | Optional | `0` | Set to `1` to emulate stop/restart behavior across MD segments. |
| `nstlim` | Conditional | `1` | Production MD steps per pH slot; excludes warmup and must be divisible by `arex_segment_steps`. |
| `dt` | Optional | `0.001` ps | MD timestep; converts the step count to simulated time. |
| `ioutfm` | Conditional | `1` with BINTRAJ; otherwise `0` | Set to `1` (NetCDF) when AREX writes coordinate, velocity, or force trajectories. |

AREX is available in the MPI build. Each pH condition owns its output files;
file names contain the associated pH label (for example, `ph4_0`).

Example AREX settings (include in a complete MD input):

```fortran
&cntrl
  iphmd = 3,
  arex = 1,
  arex_type = 'phrex',
  arex_ladder = 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0,
  arex_segment_steps = 1000,
  arex_warmup_segments = 6,
  nstlim = 20000,
  dt = 0.002,
  ioutfm = 1,
/
```

Array-valued namelist options accept **comma-separated lists** or
**individual indexed assignments**. For example, the `arex_ladder` line
above is equivalent to:

```fortran
  arex_ladder(1) = 4.0,
  arex_ladder(2) = 4.5,
  arex_ladder(3) = 5.0,
  arex_ladder(4) = 5.5,
  arex_ladder(5) = 6.0,
  arex_ladder(6) = 6.5,
  arex_ladder(7) = 7.0,
```

Use either form; list entries are assigned in order starting at index 1.

### AREX simulation length

`arex_segment_steps` defines the number of MD steps in one AREX segment.
`nstlim` sets the **production steps per pH slot**. Warmup is additional:
`arex_warmup_segments` does **not** subtract from `nstlim`.

For each pH slot:

```text
Production segments = nstlim / arex_segment_steps
Warmup steps        = arex_warmup_segments * arex_segment_steps
Total steps         = nstlim + warmup steps
Simulated time (ps) = total steps * dt
```

`nstlim` must be an exact multiple of `arex_segment_steps`.

Using the example above (`dt = 0.002` ps, or 2 fs):

| Quantity | Calculation | Result per pH slot |
|---|---|---|
| Production | `20000 / 1000` segments | 20 segments (40 ps) |
| Warmup | `6 * 1000` steps | 6 segments (12 ps) |
| Total | `20000 + 6000` steps | **26,000 steps (52 ps)** |

With seven pH values, the aggregate simulated time across all slots is
`7 * 52 = 364 ps`, equivalent to `7 * 26000 = 182000` MD steps.
This sum is **not wall-clock execution time**: elapsed runtime depends on
hardware, worker count, and exchange overhead.

### AREX positional-restraint schedule

The optional scheduler adjusts positional-restraint force constants at AREX
slot-assignment boundaries. Its progress variable is the **number of completed
MD segments for that pH slot**, beginning at zero; it is not the number of
segments completed by a particular worker.

| Keyword | Status | Default | Description |
|---|---|---|---|
| `arex_restraint_npoints` | Optional | `0` | Number of schedule points; `0` disables scheduling (maximum 100). |
| `arex_restraint_segment(i)` | Conditional | `-1` | Zero-based segment boundary for schedule point `i`. |
| `arex_restraint_force(i)` | Conditional | `0.0` | Absolute positional-restraint force constant beginning at boundary `i`. |
| `ntr` | Conditional | `0` | Must be `1` when scheduling restraints. |
| `restraintmask` | Conditional | `''` | Select atoms for positional restraints using an AMBER mask. |
| `restraint_wt` | Conditional | `0.0` | Initial, uniform positional-restraint force constant. |
| `-ref` | Conditional | — | Reference coordinates for positional restraints (command-line argument). |

Schedule rules:

- Requires `arex = 1`, `arex_warmup_segments > 0`, and `ntr = 1`.
- The first boundary must be `0`; subsequent boundaries must be strictly
  increasing and lie within `0` to `arex_warmup_segments`.
- Each force constant must be finite and nonnegative. It applies from its
  boundary up to the next one; the last value continues through production.
- The force constant applies **uniformly** to the restrained atoms. An AMBER
  `restraintmask` or GROUP-based selection is supported; the initial weights
  must be uniform.

Example: decrease the restraint force over six warmup segments, then leave
restraints off during production. The schedule arrays use comma-separated
lists here:

```fortran
  ntr = 1,
  restraintmask = '@CA',
  restraint_wt = 10.0,
  arex_warmup_segments = 6,
  arex_restraint_npoints = 3,
  arex_restraint_segment = 0, 3, 6,
  arex_restraint_force = 10.0, 5.0, 0.0,
```

The same schedule can instead use indexed assignments:

```fortran
  arex_restraint_segment(1) = 0,
  arex_restraint_segment(2) = 3,
  arex_restraint_segment(3) = 6,
  arex_restraint_force(1) = 10.0,
  arex_restraint_force(2) = 5.0,
  arex_restraint_force(3) = 0.0,
```

This gives a force constant of **10.0** for completed segment counts `0–2`,
**5.0** for `3–5`, and **0.0** starting at count `6`. Supply reference
coordinates with `-ref`. Without an enabled schedule, ordinary PMEMD
positional-restraint behavior is unchanged.

## `phmdin` (`&phmdin`)

Example titration selection and output:

```fortran
&phmdin
  masktitrres(1) = 'ASP',
  masktitrres(2) = 'GLU',
  masktitrres(3) = 'HIS',
  prlam = .true.,
  nprint_phmd = 500,
  lambda_format = 'binary',
/
```

The residue selectors may also be written as a list:

```fortran
  masktitrres = 'ASP', 'GLU', 'HIS',
```

### Residue selection and lambda output

| Keyword | Status | Default | Description |
|---|---|---|---|
| `masktitrres` | Required for titration | `''` | Select titratable residues (syntax below). |
| `prlam` | Optional | `.false.` | Enable ordinary lambda trajectory output. |
| `lambda_format` | Optional | `'binary'` | `'binary'` (`.lmb`) or `'ascii'` (`.lambda`). |
| `nprint_phmd` | Optional | `500` | Lambda output interval, in MD steps. |
| `prnlev` | Optional | `6` | Output level; lambda trajectory records require `prnlev > 2`. |
| `qphmdstart` | Optional | `.true.` | Start new CpHMD state; `.false.` reads a CpHMD restart. |

| Selector | Meaning |
|---|---|
| `'ASP'` | Select all ASP residues. |
| `'ASP:42'` | Select only topology residue 42 (if ASP). |
| `'ASP:!42'` | Select all ASP residues except topology residue 42. |

Multiple include or exclude selectors may be used for the same residue type.
Do not mix **all**, **include**, and **exclude** modes for a single type.
Selected residues must have parameters in `phmdparm`.

### Titration dynamics

| Keyword | Status | Default | Description |
|---|---|---|---|
| `qmass_phmd` | Optional | `10.0` | Mass of titration-coordinate particles. |
| `temp_phmd` | Optional | `300.0` | Titration-coordinate temperature (K). |
| `phbeta` | Optional | `5.0` | Langevin friction for titration-coordinate dynamics. |
| `iphfrq` | Optional | `1` | Titration-coordinate update interval (steps). |
| `phmdcut` | Optional | `1000.0` | GB-CpHMD interaction cutoff (Å). |
| `add_energies` | Optional | `.true.` | Include titration kinetic energy in reported kinetic energy. |
| `use_forces` | Optional | `.true.` | CpHMD force-control flag. |

### Additional output and specialized modes

| Keyword | Status | Default | Description |
|---|---|---|---|
| `prderiv` | Optional | `.false.` | Write theta and derivative pairs in ASCII. |
| `prdudls` | Optional | `.false.` | Write lambda and dU/dλ pairs in ASCII. |
| `prke` | Optional | `.false.` | Write titration kinetic energies in ASCII. |
| `qphmd_ti` | Optional | `.false.` | Enable CpHMD thermodynamic-integration mode. |
| `phtest` | Optional | `0` | `1` enables fixed-titration-coordinate testing. |
| `phtest_theta` | Conditional | Unset | Initial theta values (radians), required for `qphmd_ti`. |
| `qbarread` | Optional | `.false.` | Use barrier values from `phmdin` instead of `phmdparm`. |
| `barr` | Conditional | `0.0` | Per-coordinate barrier parameters for `qbarread`. |
| `link_res_list` | Conditional | Unset | Residue indices for linked titration groups. |
| `outu` | Optional | `mdout` unit | Output unit for CpHMD messages. |
| `masktitrrestypes` | Deprecated | Ignored | Types are inferred from `masktitrres`. |
| `nsolute` | Legacy | Not initialized | Accepted input retained for compatibility. |

If several output modes are enabled, priority is `prdudls`, `prderiv`,
`prlam`, then `prke`. Binary output applies only to ordinary lambda records.

With `qphmd_ti = .true.`, `phtest_theta` is required; the code enables
derivative output, disables `prlam` and `prdudls`, sets `phtest = 1`, and
starts fresh CpHMD state rather than reading a PHMD restart.

## Lambda trajectory formats

Set `prlam = .true.` to write lambda trajectories. Binary is the default.

| Format | Extension | Representation |
|---|---|---|
| `'binary'` | `.lmb` | Uncompressed one-byte lambda values; implicit timesteps. |
| `'ascii'` | `.lambda` | Original text format with two decimal places. |

In `.lmb`, each lambda is stored as an integer `0–100` and reconstructed by
dividing by `100`. The `LMB1` header includes the number of titration
coordinates (`ntitr`), first step (`step0`), output stride (`stride`), and
the original `ititr`, `ires`, `itauto`, and `ParK` metadata records. Each
following frame contains exactly `ntitr` bytes. No compression is used.

With `-phmdout nvt.lambda`, binary output is named `nvt.lmb`.
Convert between formats with `lambda_convert.py`:

```bash
python3 guys_lambda_file_converter.py nvt.lambda  # ASCII -> binary
python3 guys_lambda_file_converter.py nvt.lmb     # Binary -> ASCII
```

## Citation

Dayhoff, G. W., and Shen, J. (2026). *High-Performance All-Atom Continuous
Constant-pH Molecular Dynamics and Native Asynchronous Replica Exchange in
Amber*. *bioRxiv*. https://doi.org/10.64898/2026.11.11.766022

