# Experimental ABACUS CSR bulk-band adapter

This is a deliberately limited PR candidate, not a replacement for the broader
research version. It adds explicitly selected ABACUS HR/SR input for bulk band
energies to the existing nonorthogonal WannierTools implementation. It does not
add Wilson loops, Z2, Gap, surface states, or response calculations for this input
format. Those calculations remain available through the original input routes.

## Selection and scope

```fortran
&TB_FILE
 Abacus_CSR = T
 Package = 'ABACUS'
 Orthogonal_Basis = F
 Hrfile = 'HR.csr'
 Overlapfile = 'SR.csr'
 Is_HrFile = T
 Is_Sparse = F
 Is_Sparse_Hr = F
/
&CONTROL
 BulkBand_calc = T
/
&SYSTEM
 SOC = 0
 NumOccupied = 1
 E_FERMI = 0.15
/
```

Provide the normal lattice, atoms, projectors, parameters, and KPATH_BULK cards.
The included two-orbital example is a **synthetic parser/numerical fixture**, not
an ABACUS calculation or evidence about a real material.

- `Abacus_CSR` defaults to false. `Package='ABACUS'` alone does not enable it;
  existing input files keep the original reader and algorithms.
- Only `BulkBand_calc` and its `BulkBand_line_calc` alias are supported. Other
  CONTROL flags are rejected for this adapter, rather than silently ignored.
- Use `SOC=0` for real CSR entries and `SOC=1` for complex `(real,imag)` entries.
  Projector counts and SOC determine the basis dimension, which must match both
  CSR headers. Occupied-band indices must be within that dimension.
- `Nk1` must be at least 2 for the band path. Invalid range selections cause a
  nonzero exit for this adapter; the legacy error behavior is unchanged.
- Additional magnetic/Zeeman and electric-field modifications are not supported
  here; spin polarization already present in the input matrices is retained.
- The sparse **file format** is converted to upstream dense H(R)/S(R) storage.
  The flags named `Is_Sparse`/`Is_Sparse_Hr` refer to other WT representations;
  they must remain false. No sparse-memory or production-performance claim is made.

## Data contract

HR values are in Ry and are converted once to Hartree; SR is dimensionless.
CSR indices are zero-based. Each file must contain exactly one STEP record.
The supported header is the legacy `STEP:`, `Matrix Dimension of H(R):`,
`Matrix number of H(R):` layout (with `S(R)` in the overlap file).
The newer `--- Ionic Step ... ---` layout is not supported by this adapter.
Identify the format from the actual header, not the filename or version alone.
HR/SR must have the same STEP, dimension, and R-vector set; SR blocks may be
reordered. Zero-nnz blocks have no value/index/pointer arrays. R=0 must exist.
Duplicate blocks/columns, malformed pointers, non-finite values, missing files,
and extra trailing snapshots are rejected. Inputs are never symmetrized or repaired.

Each nonempty block uses one record for values, one for column indices, and one
for row pointers, with exactly the declared number of whitespace-separated fields.
Records may be long; spaces within complex `(real,imag)` pairs are allowed.
Fortran list-directed shortcuts such as `/`, null fields, repetitions, and inline
comments are not part of this CSR contract and are rejected, not filled with zero.

The existing reduced-coordinate Fourier transforms and symmetric
orthogonalization are retained. H(k)/S(k) must be Hermitian to relative tolerance
`1e-10 * max(1, max(abs(matrix)))`; overlap eigenvalues must be finite and positive.
Positive definiteness alone is **not** a conditioning or basis-quality guarantee.
Non-finite orthogonalized matrices or band energies cause an explicit nonzero
exit before band files are written; no eigenvalue clipping or regularization is used.
The requested Fermi energy is subtracted once by the existing nonorthogonal
solver. `bulkek.dat` contains only k distance and energy; `bulkek.dat-matlab`
contains all band energies. No orbital projection weights are exported for CSR.
PROJECTORS must use WannierTools orbital spellings, not ABACUS labels verbatim;
for example `dz^2` becomes `dz2`, and `fz^3` becomes `FZ3`.
This spelling conversion does not reorder or modify the input matrices.
No spin/symmetry expectation values are validated for this adapter.

## Build and self-contained tests

Use the normal WT build, or the out-of-tree test helper with GNU Fortran and a
compatible LP64 BLAS/LAPACK library. The helper reads the existing source list;
it does not modify source files or install dependencies.

```sh
python3 tests/abacus_csr/build_local.py --source . --output ../build-csr \
  --fc gfortran --libs="-llapack -lblas"
python3 tests/abacus_csr/run_tests.py --candidate ../build-csr/wt.x \
  --output ../test-csr
```

Build a clean upstream checkout with the same toolchain and pass
`--baseline ../build-upstream/wt.x` to additionally compare legacy orthogonal and
nonorthogonal bands, an unchanged ABACUS-labelled legacy input, an unused WCC
parameter in a band-only job, Gap, Wilson-loop, and surface outputs.

An optional `--research PATH` compares two synthetic band-energy cases with the
frozen research version. This is not a substitute for real-material regression.

For MPI validation build **both** executables with `--fc mpifort --mpi`, then use
`--launcher "mpirun -np 2"`. All supplied executables must be MPI builds when a
launcher is used. Repeat with one rank and compare against sequential results.
The suite has timeouts and checks exit status, expected diagnostic text, normal
completion, output sizes, and finite numbers; nonzero exit alone is not a pass.

## Review status

Local validation is recorded outside this source tree: Windows sequential
(GNU Fortran 8.1), Linux sequential (GNU Fortran 11.4), and real single-host
OpenMPI 4.1.2 runs. The self-contained suite covers analytic energies, rejected
inputs, MPI abort behavior, and seven legacy regression groups. Public Al
non-SOC and Fe SOC matrices are also compared with original ABACUS energies
and an independent generalized eigensolver, using sequential and 1/2/4-rank
MPI runs. Fe is magnetic; its test checks band energies, not time reversal.

See `VALIDATION_ABACUS_CSR.md` for the pinned sources, numerical results and
reproduction instructions. Multi-node communication, other compiler families,
broader material coverage, and independent maintainer review remain outside
this validation. This work does not establish Wilson-loop convergence or
reliable topology for any material. No PR has been submitted.

The scientific WCC, surface, Gap, and response sources are retained from upstream;
the legacy call paths are unchanged when the new flag is false. Numerical
regression demonstrates only the cases actually tested, not universal immunity
from regressions.
