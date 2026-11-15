"""Re-running a Brinkman boundary-mesh case must reuse its indicator cache.

The first run computes the indicator from an STL mesh and writes it to the
cache. The second run finds the cache and reads the indicator back instead.
Both runs must succeed and give the same indicator, compared through the
Brinkman output of each run.

The cache holds the indicator in the pressure slot of an fld file. Caches
written by earlier versions of Neko hold it as the first passive scalar, so
the cache is also re-read after relabelling it to that layout.

A Brinkman file object without a field name must read the indicator back from
the Brinkman output of a first run, which holds it as the first passive scalar.

Finally, a Brinkman file object asking for a field that its fld file does not
contain must stop with a clear error.
"""

import json
import subprocess
from pathlib import Path

import json5
import numpy as np
import pytest

import conftest
from testlib import (
    configure_nprocs,
    get_genmeshbox,
    get_neko,
    get_neko_dir,
    run_neko,
)

NPROCS = 2
# Reading the cache maps between identical spaces, so the indicators of the
# two runs differ by round-off only.
TOLERANCE = {"dp": 1.0e-12, "sp": 1.0e-5}
# Position of the rdcode, which lists the fields of an fld file, in the
# 132-byte header.
RDCODE = slice(83, 93)


def _generate_mesh(directory):
    """Generate a 32-element box mesh enclosing the bunny."""
    result = subprocess.run(
        [
            str(Path(get_genmeshbox()).resolve()),
            "0", "4", "0", "1", "0", "1", "8", "2", "2",
            ".false.", ".false.", ".false.",
        ],
        cwd=directory,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout
    return directory / "box.nmsh"


def _case(mesh, cache, output_dir):
    """Turn the immersed bunny example into a short case with a cache."""
    example = Path(get_neko_dir(), "examples", "immersed_bunny").resolve()
    with open(example / "immersed_bunny.case", "r") as f:
        config = json5.load(f)

    case = config["case"]
    case["mesh_file"] = str(mesh)
    case["output_directory"] = str(output_dir)
    case["output_at_end"] = False
    case["output_boundary"] = False
    case["numerics"]["polynomial_order"] = 5
    case["time"]["end_time"] = 2 * case["time"]["timestep"]
    case["fluid"]["output_control"] = "never"

    brinkman = case["fluid"]["source_terms"][0]
    brinkman["output"] = {
        "enable": True,
        "format": "nek5000",
        "precision": "dp",
        "path": str(output_dir),
    }
    bunny = brinkman["objects"][0]
    bunny["name"] = str(example / "bunny.stl")
    bunny["cache"] = True
    bunny["cache_file"] = str(cache)
    return config


def _run(launcher_script, nprocs, config, run_dir, log_file):
    """Write the case file into its own directory and run Neko on it."""
    run_dir.mkdir()
    case_file = run_dir / "brinkman_cache.case"
    case_file.write_text(json.dumps(config, indent=4))
    return run_neko(launcher_script, nprocs, str(case_file), get_neko(),
                    log_file)


def _rdcode(fld_file):
    """Return the rdcode of an fld file."""
    with open(fld_file, "rb") as f:
        return f.read(132)[RDCODE].decode("ascii").strip()


def _relabel(fld_file, rdcode):
    """Overwrite the rdcode of an fld file, leaving its data untouched."""
    with open(fld_file, "r+b") as f:
        f.seek(RDCODE.start)
        f.write(rdcode.ljust(RDCODE.stop - RDCODE.start).encode("ascii"))


def _read_indicator(fld_file):
    """Read the first field after the coordinates, ordered by element id.

    In the Brinkman output this is the indicator, whatever its slot.
    """
    with open(fld_file, "rb") as f:
        header = f.read(132).decode("ascii")
        f.read(4)  # Endianness test pattern
        word_size, lx, ly, lz, n_elements = map(int, header.split()[1:6])
        element_ids = np.fromfile(f, dtype="<i4", count=n_elements)
        data = np.fromfile(f, dtype="<f8" if word_size == 8 else "<f4")

    assert header[RDCODE].startswith("X"), f"No coordinates in {fld_file}"
    n = lx * ly * lz * n_elements
    values = data[3 * n:4 * n].reshape(n_elements, -1)
    return values[np.argsort(element_ids)]


@pytest.mark.parametrize("layout", ["pressure", "scalar"])
def test_brinkman_cache(layout, launcher_script, log_file, tmp_path):
    """Write the cache in a first run and read it in a second run."""
    nprocs = configure_nprocs(NPROCS)
    mesh = _generate_mesh(tmp_path)
    cache = tmp_path / "cache_"
    cache_data = tmp_path / "cache_0.f00000"
    write_dir = tmp_path / "write"
    read_dir = tmp_path / "read"
    write_log = log_file.replace(".log", "_write.log")
    read_log = log_file.replace(".log", "_read.log")

    result = _run(launcher_script, nprocs, _case(mesh, cache, write_dir),
                  write_dir, write_log)
    assert (
        result.returncode == 0
    ), f"Cache-writing run failed with exit code {result.returncode}"
    assert cache_data.is_file(), "The first run did not write the cache"
    written = _rdcode(cache_data)
    if layout == "scalar":
        _relabel(cache_data, "XS01")

    result = _run(launcher_script, nprocs, _case(mesh, cache, read_dir),
                  read_dir, read_log)
    assert (
        result.returncode == 0
    ), f"Cache-reading run failed with exit code {result.returncode}"
    with open(read_log, "r") as f:
        assert (
            f"Reading fld file {cache_data}" in f.read()
        ), "The second run did not read the cache"
    assert written == "XP", f"Cache written as {written}, not XP"

    computed = _read_indicator(write_dir / "brinkman0.f00000")
    cached = _read_indicator(read_dir / "brinkman0.f00000")
    # The bunny must cover some points, or two empty indicators would match
    assert computed.max() > 0.5 > computed.min()
    assert np.abs(cached - computed).max() <= TOLERANCE[conftest.RP]


def test_brinkman_file_default_field(launcher_script, log_file, tmp_path):
    """Read the Brinkman output of a first run without naming the field."""
    nprocs = configure_nprocs(NPROCS)
    mesh = _generate_mesh(tmp_path)
    cache = tmp_path / "cache_"
    write_dir = tmp_path / "write"
    read_dir = tmp_path / "read"
    write_log = log_file.replace(".log", "_write.log")
    read_log = log_file.replace(".log", "_read.log")

    result = _run(launcher_script, nprocs, _case(mesh, cache, write_dir),
                  write_dir, write_log)
    assert (
        result.returncode == 0
    ), f"Indicator-writing run failed with exit code {result.returncode}"
    output = write_dir / "brinkman0.f00000"
    written = _rdcode(output)
    assert written == "XS02", f"Brinkman output written as {written}, not XS02"

    config = _case(mesh, cache, read_dir)
    config["case"]["fluid"]["source_terms"][0]["objects"] = [
        {"type": "file", "file_name": str(write_dir / "brinkman0.fld")}
    ]
    result = _run(launcher_script, nprocs, config, read_dir, read_log)
    assert (
        result.returncode == 0
    ), f"File-reading run failed with exit code {result.returncode}"

    computed = _read_indicator(output)
    imported = _read_indicator(read_dir / "brinkman0.f00000")
    assert computed.max() > 0.5 > computed.min()
    assert np.abs(imported - computed).max() <= TOLERANCE[conftest.RP]


def test_brinkman_file_missing_field(launcher_script, log_file, tmp_path):
    """Request the pressure field from a file that only holds a scalar."""
    nprocs = configure_nprocs(NPROCS)
    mesh = _generate_mesh(tmp_path)
    cache = tmp_path / "cache_"
    cache_data = tmp_path / "cache_0.f00000"
    write_dir = tmp_path / "write"
    read_dir = tmp_path / "read"
    write_log = log_file.replace(".log", "_write.log")
    read_log = log_file.replace(".log", "_read.log")

    result = _run(launcher_script, nprocs, _case(mesh, cache, write_dir),
                  write_dir, write_log)
    assert (
        result.returncode == 0
    ), f"Cache-writing run failed with exit code {result.returncode}"
    assert cache_data.is_file(), "The first run did not write the cache"
    _relabel(cache_data, "XS01")

    config = _case(mesh, cache, read_dir)
    config["case"]["fluid"]["source_terms"][0]["objects"] = [
        {
            "type": "file",
            "file_name": str(tmp_path / "cache_0.fld"),
            "field_name": "p",
        }
    ]
    result = _run(launcher_script, nprocs, config, read_dir, read_log)
    assert result.returncode != 0, "Importing a missing field did not fail"
    with open(read_log, "r") as f:
        assert (
            "requested for import is not in the fld file" in f.read()
        ), "Importing a missing field did not raise the expected error"
