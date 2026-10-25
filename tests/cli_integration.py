"""Run with Python 3: python tests/cli_integration.py build/placer[.exe]."""
import pathlib
import subprocess
import sys
import tempfile

exe = pathlib.Path(sys.argv[1]).resolve()
root = pathlib.Path(__file__).resolve().parents[1]
benchmark = root / "benchmarks/sample/sample.aux"

def run(*args, ok=True):
    result = subprocess.run([str(exe), str(benchmark), *args], cwd=root,
                            capture_output=True, text=True, timeout=90)
    assert (result.returncode == 0) == ok, result.stdout + result.stderr
    return result.stdout

with tempfile.TemporaryDirectory() as temp:
    directory = pathlib.Path(temp)
    output = directory / "initial.pl"
    text = run("--method", "compare", "--no-sa", "--output", str(output))
    assert text.count("INITIAL PLACEMENT:") == 3 and "LOCAL SA AFTER:" not in text
    assert text.count("Boundary/Fixed Violations: 0") == 3
    assert "COMPARE SUMMARY" in text
    assert all(label in text for label in ("method", "HPWL", "overflow%", "runtime", "status"))
    assert all(method in text for method in ("random", "cluster", "quadratic"))
    comparison = {}
    for method in ("random", "cluster", "quadratic"):
        comparison[method] = (directory / f"initial.{method}.pl").read_bytes()
        single = directory / f"single.{method}.pl"
        text = run("--method", method, "--output", str(single))
        assert text.count("INITIAL PLACEMENT:") == 1 and "LOCAL SA AFTER:" not in text
        assert single.read_bytes() == comparison[method], "Cross-method state contamination"
    text = run("--method", "compare", "--sa", "--output", str(output))
    assert text.count("LOCAL SA AFTER:") == 3
    assert text.count("Boundary/Fixed Violations: 0") == 6
    for method in comparison:
        assert (directory / f"initial.{method}.pl").read_bytes() == comparison[method]
    text = run("--method", "compare", "--sa", "--check-only")
    assert "INITIAL PLACEMENT:" not in text and "LOCAL SA AFTER:" not in text
    assert "LOCAL SA AFTER:" not in run("--method", "compare", "--sa", "--no-sa")
    run("--method", "unknown", ok=False)
    run("--clusters", "-1", ok=False)
    run("--output", str(directory / "missing" / "failed.pl"), ok=False)
print("CLI integration tests passed: independent methods, SA, check-only, exports, errors")
