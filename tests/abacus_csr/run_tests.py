"""Self-contained synthetic tests. No ABACUS calculation or private data needed."""
import argparse
import cmath
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import time

EV_PER_HARTREE = 27.211385  # Same conversion constant as the pinned WT baseline.
EF = 0.15
RVECS = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]


def model(soc=False, orthogonal=False):
    v = 0.02 + (0.03j if soc else 0j)
    h = [[[0.2, v], [v.conjugate(), 0.4]],
         [[0.01, 0.004], [0.0, -0.02]],
         [[0.01, 0.0], [0.004, -0.02]]]
    s = [[[1., 0.1], [0.1, 1.2]], [[0.02, 0.], [0., 0.01]], [[0.02, 0.], [0., 0.01]]]
    if orthogonal:
        s = [[[1., 0.], [0., 1.]], [[0., 0.], [0., 0.]], [[0., 0.], [0., 0.]]]
    return h, s


def csr(matrices, kind, soc=False, order=None):
    lines = ['STEP: 0', f'Matrix Dimension of {kind}(R): 2', f'Matrix number of {kind}(R): 3']
    for b in (range(3) if order is None else order):
        vals, cols, pointers = [], [], [0]
        for row in matrices[b]:
            for j, value in enumerate(row):
                if value != 0:
                    vals.append(complex(value))
                    cols.append(j)
            pointers.append(len(vals))
        lines.append(' '.join(map(str, (*RVECS[b], len(vals)))))
        if vals:
            if soc:
                lines.append(' '.join(f'({x.real:.17g},{x.imag:.17g})' for x in vals))
            else:
                lines.append(' '.join(f'{x.real:.17g}' for x in vals))
            lines.append(' '.join(map(str, cols)))
            lines.append(' '.join(map(str, pointers)))
    return '\n'.join(lines) + '\n'


def write_legacy(directory, h, s, nonorthogonal):
    for filename, matrices, factor in [('hr.dat', h, 0.5 * EV_PER_HARTREE), ('sr.dat', s, 1.)]:
        sparse = nonorthogonal
        lines = ['synthetic regression fixture']
        if sparse:
            lines += ['12', '2', '3']
        else:
            lines += ['2', '3', '1 1 1']
        for r, matrix in zip(RVECS, matrices):
            for j in range(2):
                for i in range(2):
                    v = complex(matrix[i][j]) * factor
                    lines.append(' '.join(map(str, (*r, i + 1, j + 1))) + f' {v.real:.17g} {v.imag:.17g}')
        (directory / filename).write_text('\n'.join(lines) + '\n', encoding='ascii')


def wt_input(soc=False, adapter=True, nonorth=True, package=None, control='BulkBand_calc=T', extra=''):
    package = package or ('ABACUS' if adapter else 'VASP')
    tb = ("Abacus_CSR=T, Hrfile='HR.csr', Overlapfile='SR.csr'," if adapter else
          "Hrfile='hr.dat', Overlapfile='sr.dat',")
    sparse = 'F' if adapter or not nonorth else 'T'
    return f'''&TB_FILE
 {tb}
 Package='{package}', Orthogonal_Basis={"F" if nonorth else "T"}, Is_Sparse_Hr={sparse}
/
&CONTROL
 {control}
/
&SYSTEM
 SOC={1 if soc else 0}, NumOccupied=1, E_FERMI={EF}, Nslab=3
/
&PARAMETERS
 Nk1=9, Nk2=5, Nk3=3, Np=1, OmegaNum=11, OmegaMin=-1, OmegaMax=6, Eta_Arc=0.01
 {extra}
/
LATTICE
Angstrom
1 0 0
0 1 0
0 0 1
ATOM_POSITIONS
1
Direct
X 0 0 0
PROJECTORS
{1 if soc else 2}
X {"s" if soc else "s pz"}
SURFACE
1 0 0
0 1 0
KPATH_BULK
1
G 0 0 0 X 0.5 0 0
KPLANE_BULK
0 0 0
1 0 0
0 0.5 0
KCUBE_BULK
-0.5 -0.5 -0.5
1 0 0
0 1 0
0 0 1
KPATH_SLAB
1
G 0 0 X 0.5 0
KPLANE_SLAB
-0.5 -0.5
1 0
0 1
'''


def fixture(directory, soc=False, adapter=True, nonorth=True, **kwargs):
    directory.mkdir(parents=True, exist_ok=False)
    h, s = model(soc, not nonorth)
    (directory / 'HR.csr').write_text(csr(h, 'H', soc), encoding='ascii')
    (directory / 'SR.csr').write_text(csr(s, 'S', soc), encoding='ascii')
    write_legacy(directory, h, s, nonorth)
    (directory / 'wt.in').write_text(wt_input(soc, adapter, nonorth, **kwargs), encoding='ascii')
    return h, s


def read_numbers(path):
    rows = []
    for line in path.read_text(encoding='utf-8', errors='replace').splitlines():
        if not line.strip() or line.lstrip().startswith(('#', '%')):
            continue
        if path.name == 'GapCube.dat' and line.lstrip().startswith('kx (1/A)'):
            continue
        row = [float(x.replace('D', 'E')) for x in line.split()]
        if not all(math.isfinite(x) for x in row):
            raise AssertionError(f'Non-finite data: {path.name}')
        rows.append(row)
    if not rows:
        raise AssertionError(f'Empty data: {path.name}')
    return rows


def analytic_bands(h, s, k):
    def fourier(matrices, factor):
        return [[sum(factor * m[i][j] * cmath.exp(2j * math.pi * k * r[0])
                     for r, m in zip(RVECS, matrices)) for j in range(2)] for i in range(2)]
    hh, ss = fourier(h, .5), fourier(s, 1.)
    a = (ss[0][0] * ss[1][1] - ss[0][1] * ss[1][0]).real
    b = (hh[0][0] * ss[1][1] + hh[1][1] * ss[0][0] - hh[0][1] * ss[1][0] - hh[1][0] * ss[0][1]).real
    c = (hh[0][0] * hh[1][1] - hh[0][1] * hh[1][0]).real
    d = math.sqrt(b*b - 4*a*c)
    return [(b-d)/(2*a)*EV_PER_HARTREE-EF, (b+d)/(2*a)*EV_PER_HARTREE-EF]


def controls_from_source(source):
    text = (source / 'src/module.f90').read_text(encoding='utf-8')
    block = re.search(r'(?is)namelist\s*/\s*Control\s*/(.*?)(?:\n\s*\n)', text).group(1)
    return [n.strip() for n in re.sub(r'[&\r\n]', '', block).split(',') if n.strip()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--research', type=Path, help='Optional frozen ABACUS direct-reader executable')
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--dll-dir', type=Path)
    parser.add_argument('--launcher', default='', help='Optional MPI launch command, e.g. mpirun -np 2')
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    args.candidate = args.candidate.resolve()
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    if args.dll_dir:
        env['PATH'] = str(args.dll_dir.resolve()) + os.pathsep + env.get('PATH', '')
    results = []
    started = time.time()

    def launch(exe, directory, expect=None):
        command = shlex.split(args.launcher) + [str(exe)]
        run = subprocess.run(command, cwd=directory, capture_output=True, text=True,
                             encoding='utf-8', errors='replace', timeout=120, env=env)
        (directory / 'process.stdout').write_text(run.stdout, encoding='utf-8')
        (directory / 'process.stderr').write_text(run.stderr, encoding='utf-8')
        log = (directory / 'WT.out').read_text(encoding='utf-8', errors='replace') if (directory / 'WT.out').exists() else ''
        text = log + run.stdout + run.stderr
        (directory / 'execution.json').write_text(json.dumps({'command': command, 'returncode': run.returncode}, indent=2))
        if expect:
            assert run.returncode != 0, 'Error returned zero exit status'
            assert 'ERROR: ABACUS CSR:' in text and expect in text, f'Expected adapter diagnostic: {expect}'
            assert 'whole program' not in log, 'Error reached normal completion'
            assert 'Fortran runtime error' not in text, 'Runtime crash instead of validation error'
        else:
            assert run.returncode == 0, f'Nonzero exit: {run.returncode}'
            assert 'whole program' in log, 'Missing normal completion marker'
            assert 'ERROR' not in text, 'Unexpected error diagnostic'
        return directory

    def record(name, test):
        try:
            details = test() or {}
            results.append({'name': name, 'pass': True, **details})
        except Exception as exc:
            results.append({'name': name, 'pass': False, 'error': str(exc)})
        print(json.dumps(results[-1]), flush=True)
        payload = {'scope': 'synthetic numerical and compatibility tests, not material topology validation',
                   'candidate_sha256': hashlib.sha256(args.candidate.read_bytes()).hexdigest(),
                   'mpi_launcher': args.launcher or None, 'results': results,
                   'passed': sum(x['pass'] for x in results), 'total': len(results),
                   'all_pass': all(x['pass'] for x in results), 'seconds': time.time()-started}
        (args.output / 'results.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')

    def valid_test(name, soc=False, reverse_sr=False, zero_block=False, layout=None):
        d = args.output / name
        h, s = fixture(d, soc=soc)
        if zero_block:
            h[1] = h[2] = [[0., 0.], [0., 0.]]
            s[1] = s[2] = [[0., 0.], [0., 0.]]
            (d / 'HR.csr').write_text(csr(h, 'H', soc))
        if reverse_sr or zero_block:
            (d / 'SR.csr').write_text(csr(s, 'S', soc, [2, 0, 1] if reverse_sr else None))
        for filename in ('HR.csr', 'SR.csr'):
            p = d / filename
            text = p.read_text(encoding='ascii')
            if layout == 'long_record':
                lines = text.splitlines()
                lines[4] = ' ' * 9000 + lines[4]
                text = '\n'.join(lines) + '\n'
            elif layout == 'whitespace':
                lines = text.splitlines()
                lines[3:] = [line.replace(' ', '\t') if not soc else line.replace(',', ',\t ')
                             for line in lines[3:]]
                text = '\n'.join(lines) + '\n'
            elif layout == 'no_final_newline':
                text = text.rstrip('\n')
            if layout:
                p.write_text(text, encoding='ascii')
        if layout == 'selected_range':
            p = d / 'wt.in'
            p.write_text(p.read_text()+'SELECTED_OCCUPIED_BANDS\n1-2\n', encoding='ascii')
        launch(args.candidate, d)
        rows = read_numbers(d / 'bulkek.dat-matlab')
        assert len(rows) == 9 and all(len(r) == 3 for r in rows), 'Wrong band output dimensions'
        assert all(len(r) == 2 for r in read_numbers(d/'bulkek.dat')), 'Unexpected orbital-weight output'
        error = max(abs(row[j+1] - analytic_bands(h, s, .5*i/8)[j])
                    for i, row in enumerate(rows) for j in range(2))
        assert error < 2e-8, f'Analytic generalized eigenvalue error {error}'
        return {'max_error_eV': error}

    record('csr_real_bands', lambda: valid_test('csr_real_bands'))
    if args.smoke:
        raise SystemExit(0 if results[-1]['pass'] else 1)
    record('csr_complex_soc_bands', lambda: valid_test('csr_complex_soc_bands', soc=True))
    record('sr_reordered_blocks', lambda: valid_test('sr_reordered_blocks', reverse_sr=True))
    record('zero_nnz_blocks', lambda: valid_test('zero_nnz_blocks', zero_block=True))
    for name, soc, layout in [
        ('csr_long_record', False, 'long_record'),
        ('csr_tabs', False, 'whitespace'),
        ('csr_complex_spaces', True, 'whitespace'),
        ('csr_no_final_newline', False, 'no_final_newline'),
        ('csr_selected_range', False, 'selected_range'),
    ]:
        record(name, lambda n=name, s=soc, f=layout: valid_test(n, soc=s, layout=f))

    def invalid(name, filename, transform, expected, soc=False):
        d = args.output / name
        fixture(d, soc=soc)
        p = d / filename
        text = p.read_text(encoding='ascii')
        p.write_text(transform(text), encoding='ascii')
        launch(args.candidate, d, expected)

    cases = [
        ('dimension_mismatch', 'HR.csr', lambda t: t.replace('Dimension of H(R): 2', 'Dimension of H(R): 3'), 'dimension differs'),
        ('negative_nnz', 'HR.csr', lambda t: t.replace('0 0 0 4', '0 0 0 -1', 1), 'negative nonzero count'),
        ('column_out_of_range', 'HR.csr', lambda t: t.replace('0 1 0 1', '0 2 0 1', 1), 'column index out of range'),
        ('row_out_of_range', 'HR.csr', lambda t: t.replace('0 2 4', '0 5 4', 1), 'row pointer out of range'),
        ('row_nonmonotonic', 'HR.csr', lambda t: t.replace('0 2 4', '0 3 2', 1), 'row-pointer endpoints'),
        ('duplicate_column', 'HR.csr', lambda t: t.replace('0 1 0 1', '0 0 0 1', 1), 'duplicate CSR column'),
        ('nonfinite_values', 'HR.csr', lambda t: t.replace('0.20000000000000001', 'NaN', 1), 'non-finite CSR value'),
        ('truncated_file', 'HR.csr', lambda t: t[:t.index('0 0 0 4')], 'invalid R-block header'),
        ('trailing_snapshot', 'HR.csr', lambda t: t+'STEP: 1\n', 'unexpected trailing data'),
        ('step_mismatch', 'SR.csr', lambda t: t.replace('STEP: 0', 'STEP: 1'), 'STEP values differ'),
        ('sr_dimension', 'SR.csr', lambda t: t.replace('Dimension of S(R): 2', 'Dimension of S(R): 3'), 'dimension differs'),
        ('sr_missing_r', 'SR.csr', lambda t: t.replace('1 0 0 2', '2 0 0 2', 1), 'absent from HR'),
        ('hr_duplicate_r', 'HR.csr', lambda t: t.replace('-1 0 0 3', '1 0 0 3'), 'duplicate HR R vector'),
        ('sr_duplicate_r', 'SR.csr', lambda t: t.replace('-1 0 0 2', '1 0 0 2'), 'duplicate SR R vector'),
        ('missing_hr', 'wt.in', lambda t: t.replace("Hrfile='HR.csr'", "Hrfile='missing.csr'"), 'cannot open missing.csr'),
        ('missing_sr', 'wt.in', lambda t: t.replace("Overlapfile='SR.csr'", "Overlapfile='missing.csr'"), 'cannot open missing.csr'),
        ('zero_occupied', 'wt.in', lambda t: t.replace('NumOccupied=1', 'NumOccupied=0'), 'NumOccupied must be positive'),
        ('excess_occupied', 'wt.in', lambda t: t.replace('NumOccupied=1', 'NumOccupied=3'), 'NumOccupied exceeds basis dimension'),
        ('negative_occupied', 'wt.in', lambda t: t.replace('NumOccupied=1', 'NumOccupied=-1'), 'NumOccupied must be positive'),
        ('bad_selected_band', 'wt.in', lambda t: t+'SELECTED_OCCUPIED_BANDS\n3\n', 'selection is out of range'),
        ('wrong_package', 'wt.in', lambda t: t.replace("Package='ABACUS'", "Package='VASP'"), 'requires Package=ABACUS'),
        ('orthogonal_flag', 'wt.in', lambda t: t.replace('Orthogonal_Basis=F', 'Orthogonal_Basis=T'), 'requires Orthogonal_Basis=F'),
        ('sparse_flag', 'wt.in', lambda t: t.replace('Is_Sparse_Hr=F', 'Is_Sparse_Hr=T'), 'Is_Sparse=Is_Sparse_Hr=F'),
        ('invalid_tb', 'wt.in', lambda t: t.replace("Hrfile='HR.csr'", "unknown_keyword=1, Hrfile='HR.csr'"), 'invalid TB_FILE'),
        ('invalid_parameters', 'wt.in', lambda t: t.replace('Nk1=9', 'unknown_keyword=1, Nk1=9'), 'invalid PARAMETERS'),
        ('electric_field', 'wt.in', lambda t: t.replace('Nslab=3', 'Nslab=3, Electric_field_in_eVpA=1'), 'electric-field modification'),
        ('nonfinite_fermi', 'wt.in', lambda t: t.replace('E_FERMI=0.15', 'E_FERMI=NaN'), 'E_FERMI must be finite'),
        ('negative_overlap', 'SR.csr', lambda t: t.replace('\n1 0.10000000000000001', '\n-1 0.10000000000000001', 1), 'not positive definite'),
        ('nonhermitian_h', 'HR.csr', lambda t: t.replace('0.02 0.02', '0.03 0.02', 1), 'H(k) is not Hermitian'),
        ('zero_band_points', 'wt.in', lambda t: t.replace('Nk1=9', 'Nk1=0'), 'Nk1 must be at least 2'),
        ('one_band_point', 'wt.in', lambda t: t.replace('Nk1=9', 'Nk1=1'), 'Nk1 must be at least 2'),
        ('selected_invalid_token', 'wt.in', lambda t: t+'SELECTED_OCCUPIED_BANDS\nx\n', 'SelectedOccupiedBands'),
        ('selected_empty', 'wt.in', lambda t: t+'SELECTED_OCCUPIED_BANDS\n\n', 'SelectedOccupiedBands'),
        ('selected_eof', 'wt.in', lambda t: t+'SELECTED_OCCUPIED_BANDS\n', 'SELECTED_OCCUPIED_BANDS'),
        ('selected_duplicate', 'wt.in', lambda t: t+'SELECTED_OCCUPIED_BANDS\n1,1\n', 'duplicate values'),
        ('selected_bad_range', 'wt.in', lambda t: t+'SELECTED_OCCUPIED_BANDS\n1-x\n', 'SelectedOccupiedBands'),
        ('selected_missing_endpoint', 'wt.in', lambda t: t+'SELECTED_OCCUPIED_BANDS\n1-\n', 'missing range endpoint'),
        ('selected_invalid_range_prefix', 'wt.in', lambda t: t+'SELECTED_OCCUPIED_BANDS\n1-x2\n', 'incorrect range'),
        ('near_singular_overlap', 'SR.csr', lambda t: csr([
            [[1e-310, 0.], [0., 1.]], [[0., 0.], [0., 0.]], [[0., 0.], [0., 0.]]], 'S'),
         'non-finite orthogonalized H(k)'),
    ]
    for name, filename, transform, expected in cases:
        record(name, lambda n=name, f=filename, t=transform, e=expected: invalid(n, f, t, e))

    def replace_record(text, index, replacement):
        lines = text.splitlines()
        lines[index] = replacement
        return '\n'.join(lines) + '\n'

    malformed_records = [
        ('step_slash', 0, 'STEP: /', 'header: STEP:'),
        ('step_extra', 0, 'STEP: 0 1', 'header: STEP:'),
        ('dimension_null', 1, 'Matrix Dimension of H(R): ,', 'header:'),
        ('block_missing', 3, '0 0 0', 'R-block header'),
        ('block_slash', 3, '0 0 0 /', 'R-block header'),
        ('values_slash', 4, '/', 'CSR value array'),
        ('values_partial', 4, '0.2 /', 'CSR value array'),
        ('values_missing', 4, '0.2 0.02 0.02', 'CSR value array'),
        ('values_extra', 4, '0.2 0.02 0.02 0.4 999', 'CSR value array'),
        ('values_null', 4, '0.2,,0.02,0.4', 'CSR value array'),
        ('values_repeat', 4, '4*', 'CSR value array'),
        ('values_comment', 4, '0.2 ! missing values', 'CSR value array'),
        ('columns_slash', 5, '/', 'CSR column array'),
        ('columns_missing', 5, '0 1 0', 'CSR column array'),
        ('columns_extra', 5, '0 1 0 1 0', 'CSR column array'),
        ('pointers_slash', 6, '/', 'CSR row-pointer array'),
        ('pointers_null', 6, '0,,4', 'CSR row-pointer array'),
        ('pointers_extra', 6, '0 2 4 4', 'CSR row-pointer array'),
    ]
    for name, index, replacement, expected in malformed_records:
        record(name, lambda n=name, i=index, r=replacement, e=expected: invalid(
            n, 'HR.csr', lambda t: replace_record(t, i, r), e))
    for name, replacement in [
        ('complex_slash', '/'),
        ('complex_partial', '(0.2,0) /'),
        ('complex_null_real', '(,0) (0.02,0.03) (0.02,-0.03) (0.4,0)'),
        ('complex_null_imag', '(0.2,) (0.02,0.03) (0.02,-0.03) (0.4,0)'),
        ('complex_extra', '(0.2,0) (0.02,0.03) (0.02,-0.03) (0.4,0) (9,0)'),
        ('complex_unclosed', '(0.2,0 (0.02,0.03) (0.02,-0.03) (0.4,0)'),
    ]:
        record(name, lambda n=name, r=replacement: invalid(
            n, 'HR.csr', lambda t: replace_record(t, 4, r), 'CSR value array', soc=True))

    if args.research:
        def research_compare(soc):
            d = args.output / ('research_soc' if soc else 'research_real')
            fixture(d, soc=soc)
            p = d / 'wt.in'
            p.write_text(p.read_text().replace('Abacus_CSR=T, ', ''), encoding='ascii')
            launch(args.research.resolve(), d)
            reference = args.output / ('csr_complex_soc_bands' if soc else 'csr_real_bands')
            a, b = read_numbers(d/'bulkek.dat-matlab'), read_numbers(reference/'bulkek.dat-matlab')
            assert len(a) == len(b) and all(len(x) == len(y) for x,y in zip(a,b))
            error = max(abs(x-y) for ar,br in zip(a,b) for x,y in zip(ar,br))
            assert error < 2e-8, f'Research energy mismatch: {error}'
            return {'max_error_eV': error, 'scope': 'synthetic energy outputs only'}
        record('research_real_energy', lambda: research_compare(False))
        record('research_soc_energy', lambda: research_compare(True))

    for flag in controls_from_source(args.source):
        if flag.lower() in ('bulkband_calc', 'bulkband_line_calc'):
            continue
        name = 'reject_' + flag.lower()
        record(name, lambda n=name, f=flag: invalid(n, 'wt.in',
               lambda t: t.replace('BulkBand_calc=T', f'BulkBand_calc=T, {f}=T'), 'unsupported CONTROL flag: '+f))

    if args.baseline:
        args.baseline = args.baseline.resolve()
        legacy_cases = [
            ('legacy_orthogonal_bands', {}, ['bulkek.dat']),
            ('legacy_package_abacus', {'package': 'ABACUS'}, ['bulkek.dat']),
            ('legacy_unused_wcc_parameter', {'extra': 'wcc_calc_tol=0.6'}, ['bulkek.dat']),
            ('legacy_nonorthogonal_bands', {'nonorth': True}, ['bulkek.dat']),
            ('legacy_gap', {'control': 'BulkGap_cube_calc=T', 'extra': 'Gap_threshold=100'}, ['GapCube.dat']),
            ('legacy_wilson_loop', {'control': 'Wilsonloop_calc=T'}, ['wcc.dat']),
            ('legacy_surface', {'control': 'SlabSS_calc=T'}, ['dos.dat_l', 'dos.dat_r', 'dos.dat_bulk']),
        ]

        def legacy(name, kw, files):
            directories = []
            for label, exe in [('upstream', args.baseline), ('candidate', args.candidate)]:
                d = args.output / (name + '_' + label)
                params = {'adapter': False, 'nonorth': False, **kw}
                fixture(d, **params)
                launch(exe, d)
                directories.append(d)
            maximum = 0.
            for file in files:
                a, b = [read_numbers(d/file) for d in directories]
                if file == 'GapCube.dat':
                    assert len(a) == 9*5*3, 'Gap grid is incomplete'
                assert len(a) == len(b) and all(len(x) == len(y) for x,y in zip(a,b)), 'Output dimensions changed'
                maximum = max(maximum, max(abs(x-y) for ar,br in zip(a,b) for x,y in zip(ar,br)))
                assert a == b, f'Legacy numerical output changed: {file}'
            return {'compared_files': files, 'max_absolute_difference': maximum}

        for name, kw, files in legacy_cases:
            record(name, lambda n=name,k=kw,f=files: legacy(n,k,f))
    print(json.dumps({'passed': sum(r['pass'] for r in results), 'total': len(results)}))
    raise SystemExit(0 if all(r['pass'] for r in results) else 1)


if __name__ == '__main__':
    main()
