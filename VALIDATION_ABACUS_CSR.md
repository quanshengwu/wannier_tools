# ABACUS CSR adapter: local validation

Validation date: 2026-11-03. Upstream baseline:
`362171ce35450118828c6960a882e3a6d08c73f5`.
This documents a proposed, opt-in, bulk-band-only input adapter, not the
production research implementation or a merged upstream feature.

## Self-contained regression

| Environment | Passed / total | Compared legacy outputs |
|---|---:|---|
| Windows, GNU Fortran 8.1, sequential | 166 / 166 | No numerical difference |
| Ubuntu 22.04 on WSL2, GNU Fortran 11.4, sequential | 164 / 164 | No numerical difference |
| Same Linux build, OpenMPI 4.1.2, 2 ranks | 164 / 164 | No numerical difference |
| Same Linux build, OpenMPI 4.1.2, 4 ranks | 164 / 164 | No numerical difference |

These are repeated configurations of the same tests, not independent materials.
Windows includes two optional comparisons with the frozen research executable.
Linux does not count those two tests. Builds enable runtime checks and use
single-thread BLAS. MPI is real multi-process execution on one host; neither
multi-node networking nor performance scaling is tested.

Rejected-input tests require a nonzero exit, the expected adapter diagnostic,
no normal-completion marker and no Fortran runtime crash. A launcher failure
or timeout does not count as a successful rejection. The seven legacy groups
cover orthogonal bands, ABACUS-labelled old input, an unused WCC parameter,
nonorthogonal bands, Gap, Wilson loop and surface spectral outputs.

## Public real ABACUS matrices

Both original HR/SR files are preserved byte-for-byte. At original SCF k points,
the test independently constructs H(k), S(k) and uses SciPy's generalized
Hermitian eigensolver. It compares absolute energies with E_FERMI=0, with no
fitted shift, symmetrization, matrix repair or band selection to hide errors.

| Sample | Basis | Unique k points | Max error vs ABACUS (eV) | Max error vs SciPy (eV) |
|---|---:|---:|---:|---:|
| Al, non-SOC | 21 | 47 | 4.2749e-5 | <5.0e-10 |
| Fe, SOC | 108 | 35 | 7.6799e-5 | <4.3e-9 |

Sequential and MPI 1/2/4-rank energies have zero difference at output precision
for both samples. Each run also checks completion, finite numbers and shapes.
ABACUS reference thresholds were set to 1e-4 eV (Al) and 2e-4 eV (Fe), before
comparing WT output, to accommodate the printed reference precision. All WT
bands are compared with SciPy at 1e-6 eV; the ABACUS comparison uses every band
available in the original reference output (which may contain fewer bands).

- Al: [deepmodeling/dftio](https://github.com/deepmodeling/dftio), commit
  `c9d128f24a74ef2911e1a28f5640357488beb196`, `test/data/abacus_scf/`.
- Fe: [mailhexu/TB2J_examples](https://github.com/mailhexu/TB2J_examples), commit
  `97624aa3756f63a49f6752a3e900cafc7df31d13`, `Abacus/Fe_SOC/DFT/`.

Public fixture manifests record Git blob IDs and full SHA256 hashes. Both
repositories include MIT licenses; the reproduction download retains them.
Al is metallic and Fe is magnetic: these are energy/input tests, not evidence
of an insulating occupied manifold, time reversal, Wilson loops or Z2.

## Reproduction

Official Graphene and Bi2Se3 Wannier examples from the same pinned upstream
commit also have identical candidate/upstream band energies and projection
weights, in sequential and 2-rank MPI tests. The Graphene input is unchanged.
For Bi2Se3, only the extra DOS and discrete-point control flags are disabled;
the original Hamiltonian, path and physical parameters are unchanged. The
original upstream DOS run exceeded the local 180-second test budget and is
not counted as validated. Its files are retained in the local audit.

Use `tests/abacus_csr/build_local.py` and `run_tests.py` as documented in
`README_ABACUS_CSR.md`, with matching upstream and candidate builds. Repeat
the Linux suite without a launcher, then with `mpirun -np 2` and `mpirun -np 4`.
Use a separate output directory for each run.

The accompanying small `public_validation_bundle.zip` contains independent
reference scripts and pinned download manifests, not large matrix files.
On a normal Python environment with NumPy and SciPy, run:

```sh
python3 run_public_checks.py --seq /path/to/sequential/wt.x \
  --mpi /path/to/mpi/wt.x --mpirun mpirun --output /path/to/new/results
```

The script downloads only pinned public data and validates SHA256 before use.
It does not execute repository scripts, rebuild ABACUS, access a cluster, or
change original data. A cached fixture directory can be supplied explicitly.

## Remaining boundaries

- Only the legacy `STEP:` layout is supported, not the new ionic-step layout.
  See the official [CSR migration guide](https://abacus-develop.readthedocs.io/en/latest/advanced/interface/migration-guide-csr-format.html).
- The adapter stores matrices densely and supports bulk-band energies only.
- Other compiler families, multi-node MPI and broad material coverage remain
  untested. Passing selected legacy tests is not a proof that every old feature
  is unaffected.
- Independent maintainer review and agreement on API/format/documentation scope
  are still required. No PR has been submitted by this validation workflow.
