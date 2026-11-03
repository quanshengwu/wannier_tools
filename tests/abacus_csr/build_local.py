"""Build unmodified source files out of tree with a caller-supplied toolchain."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import time


def build(source, output, compiler, libraries, mpi=False):
    source, output = Path(source).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    production_makefile = source / "build/Makefile.seq"
    if production_makefile.is_file():
        makefile = production_makefile.read_text(encoding="utf-8")
        block = makefile.split("SOURCES =", 1)[1].split("OBJECTS", 1)[0]
        objects = re.findall(r"([\w]+)\.f90\b", block)
    else:
        makefile = (source / "src/Makefile.gfortran").read_text(encoding="utf-8")
        block = makefile.split("OBJ =", 1)[1].split("# compiler", 1)[0]
        objects = re.findall(r"([\w]+)\.o\b", block)
    assert objects and len(objects) == len(set(objects)), "Invalid source list"
    flags = ["-cpp", "-O0", "-g", "-fcheck=all", "-ffree-line-length-none"]
    version = subprocess.check_output([compiler, "--version"], text=True).splitlines()[0]
    if mpi:
        flags += ["-DMPI"]
        # Existing upstream mpif.h callers mix argument types across calls.
        major = re.search(r"\b(\d+)\.\d+\.\d+", version)
        if major and int(major.group(1)) >= 10:
            flags += ["-fallow-argument-mismatch"]
    records = []
    for name in objects:
        path = source / "src" / (name + ".f90")
        command = [compiler, *flags, "-I" + str(source / "src"), "-c", str(path), "-o", name + ".o"]
        start = time.time()
        run = subprocess.run(command, cwd=output, text=True, capture_output=True, timeout=180)
        (output / (name + ".log")).write_text(run.stdout + run.stderr, encoding="utf-8")
        records.append({"source": name + ".f90", "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "command": command, "returncode": run.returncode, "seconds": time.time() - start})
        (output / "build.json").write_text(json.dumps({"compiler": version, "records": records}, indent=2))
        if run.returncode:
            raise RuntimeError(f"Compile failed: {name}; see {output / (name + '.log')}")
    executable = output / ("wt.exe" if os.name == "nt" else "wt.x")
    command = [compiler, *[n + ".o" for n in objects], *libraries, "-o", str(executable)]
    run = subprocess.run(command, cwd=output, capture_output=True, text=True, timeout=180)
    (output / "link.log").write_text(run.stdout + run.stderr, encoding="utf-8")
    if run.returncode:
        raise RuntimeError(f"Link failed; see {output / 'link.log'}")
    result = {"compiler": version, "records": records, "link_command": command,
              "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest()}
    (output / "build.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({"compiled_files": len(objects), "executable": str(executable)}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fc", default=shutil.which("gfortran"))
    parser.add_argument("--libs", default="-llapack -lblas")
    parser.add_argument("--mpi", action="store_true")
    args = parser.parse_args()
    if not args.fc:
        parser.error("Fortran compiler not found; use --fc")
    build(args.source, args.output, args.fc, shlex.split(args.libs), args.mpi)
