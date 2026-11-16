"""The shape of the package: error mapping, type stubs, docstrings, `help()` and the README's examples."""

from __future__ import annotations

import ast
import contextlib
import importlib
import importlib.util
import inspect
import io
import math
import pydoc
import re
import subprocess
import sys
import types
from pathlib import Path

import pytest

import srkit
from helpers import INF, ORIGIN, X_AXIS
from srkit import Segment, Time384, Worldline, causal, kinematics, optics, plan

CRATE = Path(__file__).resolve().parents[1]
PACKAGE = Path(srkit.__file__).resolve().parent
MODULES = ["srkit", *(f"srkit.{name}" for name in ("causal", "kepler", "kinematics", "optics", "plan", "rocket", "units"))]


# ------------------------------------------------------------------------------------------------ errors


def test_core_errors_map_to_python_exceptions():
    """Each variant of the core's `Error`, reached through the public API."""
    # InvalidInput -> ValueError, with the core's own message
    with pytest.raises(ValueError, match="invalid or non-finite physical input"):
        srkit.Trip(-1.0, 1.0)
    # OutOfRange -> ValueError
    with pytest.raises(ValueError, match="outside the supported domain"):
        Time384(1).duration_since(Time384(2))
    # Discontinuous -> ValueError
    apart = Worldline([Segment.rest(0.0, 1.0, ORIGIN), Segment.rest(1.0, 2.0, X_AXIS)])
    with pytest.raises(ValueError, match="do not join continuously"):
        apart.check_continuity()
    # Overflow -> OverflowError
    with pytest.raises(OverflowError, match="exceeds the representation"):
        Time384.MAX + 1
    # PrecisionLoss -> srkit.PrecisionLoss, an ArithmeticError
    with pytest.raises(srkit.PrecisionLoss, match="promised precision"):
        Time384(10**30).local_seconds_since(Time384(0))
    # no intersection -> srkit.NoIntersection
    with pytest.raises(srkit.NoIntersection, match="does not reach"):
        Segment.rest(0.0, 1.0, ORIGIN).recv(0.0, (5.0, 0.0, 0.0), 0.0, 1.0)
    with pytest.raises(srkit.NoIntersection):
        Worldline([Segment.hyperbolic(0.0, INF, ORIGIN, X_AXIS, 1.0)]).seen_from(-5.0, (-3.0, 0.0, 0.0))


def test_the_exception_classes():
    assert srkit.NoIntersection.__mro__[1] is Exception
    assert srkit.PrecisionLoss.__mro__[1] is ArithmeticError
    assert srkit.NoIntersection.__module__ == srkit.PrecisionLoss.__module__ == "srkit"
    assert srkit.NoIntersection.__name__ == "NoIntersection"
    # neither is caught by a handler for bad input
    for exc in (srkit.NoIntersection, srkit.PrecisionLoss):
        assert not issubclass(exc, (ValueError, TypeError, OverflowError))
        with pytest.raises(exc, match="by hand"):
            raise exc("by hand")


def test_validation_added_by_the_bindings_is_value_error():
    """Where the core's constructors take anything, the bindings refuse what is not physics."""
    for bad in [
        lambda: Segment.inertial(0.0, 1.0, ORIGIN, (2.0, 0.0, 0.0)),
        lambda: srkit.CircularOrbit(1.0, 1.5),
        lambda: plan.line(plan.Start(), (1.0, 0.0, 0.0), -1.0),
        lambda: srkit.rocket.mass_ratio(1.0, 2.0),
        lambda: optics.playback((1.0, 0.0, 0.0), X_AXIS, ORIGIN),
        lambda: srkit.units.Units(0.0),
        lambda: Worldline([]),
    ]:
        with pytest.raises(ValueError) as caught:
            bad()
        assert str(caught.value)  # every refusal says why


def test_wrong_types_are_type_errors():
    for bad in [
        lambda: Segment.rest(0.0, 1.0, None),
        lambda: Segment.rest(None, 1.0, ORIGIN),
        lambda: srkit.Trip("far", 1.0),
        lambda: srkit.reach(0.0, 1.0),
        lambda: kinematics.gamma(0.5),
        lambda: causal.light_ticks((1.0, 2.0, 3.0)),
        lambda: Time384(1.5),
        lambda: plan.line(None, ORIGIN, 1.0),
    ]:
        with pytest.raises(TypeError):
            bad()


# ------------------------------------------------------------------------------------------------- stubs


def stub_tree(module_name: str) -> ast.Module:
    relative = "__init__" if module_name == "srkit" else module_name.split(".", 1)[1]
    path = PACKAGE / f"{relative}.pyi"
    assert path.is_file(), f"no stub file for {module_name}"
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def is_public(name: str) -> bool:
    return not name.startswith("_") or name == "__version__"


def stub_names(tree: ast.Module) -> set[str]:
    """The public names a stub file defines or re-exports."""
    names: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.ClassDef, ast.FunctionDef)):
            names.add(node.name)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
        elif isinstance(node, ast.Assign):
            names.update(t.id for t in node.targets if isinstance(t, ast.Name))
        elif isinstance(node, (ast.ImportFrom, ast.Import)):
            # `import x as x` is the stub convention for a re-export
            names.update(a.asname for a in node.names if a.asname and a.asname == a.name.split(".")[-1])
    return {n for n in names if is_public(n)}


def stub_all(tree: ast.Module) -> list[str]:
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "__all__" for t in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError("the stub has no __all__")


@pytest.mark.parametrize("module_name", MODULES)
def test_all_matches_the_stub_file(module_name):
    module = importlib.import_module(module_name)
    tree = stub_tree(module_name)
    assert len(module.__all__) == len(set(module.__all__))
    assert stub_names(tree) == set(module.__all__)
    assert sorted(stub_all(tree)) == sorted(module.__all__)
    # and __all__ is the whole public surface: nothing public is left out of it
    public = {n for n in vars(module) if not n.startswith("_")} | ({"__version__"} if module_name == "srkit" else set())
    assert public == set(module.__all__)
    for name in module.__all__:
        assert hasattr(module, name)


def test_the_submodules_import_as_modules():
    for module_name in MODULES[1:]:
        module = importlib.import_module(module_name)
        short = module_name.split(".")[1]
        assert module is getattr(srkit, short) and module.__name__ == module_name
        assert isinstance(module, types.ModuleType)
    from srkit.causal import light_ticks  # noqa: PLC0415
    from srkit.units import G0_PER_YEAR, Units  # noqa: PLC0415

    assert light_ticks is causal.light_ticks and Units is srkit.units.Units
    assert G0_PER_YEAR == pytest.approx(1.0323, abs=1e-4)
    assert (PACKAGE / "py.typed").is_file()


def stub_classes(tree: ast.Module) -> dict[str, ast.ClassDef]:
    return {n.name: n for n in tree.body if isinstance(n, ast.ClassDef) and is_public(n.name)}


def class_members(node: ast.ClassDef) -> set[str]:
    members: set[str] = set()
    for item in node.body:
        if isinstance(item, ast.FunctionDef):
            members.add(item.name)
        elif isinstance(item, ast.AnnAssign) and isinstance(item.target, ast.Name):
            members.add(item.target.id)
    return {m for m in members if not m.startswith("_")}


def parameters(node: ast.FunctionDef) -> list[tuple[str, bool]]:
    """`(name, has a default)` for each parameter of a stub function, without self or cls."""
    args = node.args
    positional = [*args.posonlyargs, *args.args]
    defaults = [False] * (len(positional) - len(args.defaults)) + [True] * len(args.defaults)
    out = list(zip((a.arg for a in positional), defaults))
    out += [(a.arg, d is not None) for a, d in zip(args.kwonlyargs, args.kw_defaults)]
    return [(name, has_default) for name, has_default in out if name not in ("self", "cls")]


def runtime_parameters(obj) -> list[tuple[str, bool]]:
    signature = inspect.signature(obj)
    return [
        (p.name, p.default is not inspect.Parameter.empty)
        for p in signature.parameters.values()
        if p.name not in ("self", "cls")
    ]


def is_property(node: ast.FunctionDef) -> bool:
    return any(isinstance(d, ast.Name) and d.id == "property" for d in node.decorator_list)


@pytest.mark.parametrize("module_name", MODULES)
def test_stubs_list_the_members_and_signatures_the_module_has(module_name):
    module = importlib.import_module(module_name)
    tree = stub_tree(module_name)
    checked = 0
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and is_public(node.name):
            assert parameters(node) == runtime_parameters(getattr(module, node.name)), node.name
            checked += 1
    for name, node in stub_classes(tree).items():
        cls = getattr(module, name)
        if issubclass(cls, BaseException):
            continue
        assert class_members(node) == {m for m in dir(cls) if not m.startswith("_")}, name
        for item in node.body:
            if not isinstance(item, ast.FunctionDef):
                continue
            if item.name == "__new__":
                assert parameters(item) == runtime_parameters(cls), f"{name}()"
                checked += 1
            elif item.name.startswith("__"):
                assert hasattr(cls, item.name), f"{name}.{item.name}"
            elif is_property(item):
                attribute = inspect.getattr_static(cls, item.name)
                assert inspect.isgetsetdescriptor(attribute) or isinstance(attribute, property), f"{name}.{item.name}"
            else:
                assert parameters(item) == runtime_parameters(getattr(cls, item.name)), f"{name}.{item.name}"
                checked += 1
    assert checked > 0 or module_name == "srkit.units"


def test_mypy_stubtest_agrees(tmp_path):
    """mypy's own comparison of the stubs with the module, where mypy is installed."""
    pytest.importorskip("mypy")
    # PyO3 fills the reflected slots of `-` and `%` on the time types; they only return NotImplemented.
    allowlist = tmp_path / "allowlist.txt"
    allowlist.write_text("srkit\\.Time(384|512)\\.__r(sub|mod)__\n", encoding="utf-8")
    result = subprocess.run(
        [sys.executable, "-m", "mypy.stubtest", "srkit", "--allowlist", str(allowlist)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


# -------------------------------------------------------------------------------------------- docstrings


def documented(obj) -> bool:
    doc = inspect.getdoc(obj)
    return bool(doc and doc.strip())


@pytest.mark.parametrize("module_name", MODULES)
def test_everything_public_has_a_docstring(module_name):
    module = importlib.import_module(module_name)
    assert documented(module), module_name
    missing = []
    for name in module.__all__:
        obj = getattr(module, name)
        if isinstance(obj, (int, float, str)):
            continue  # constants are documented in their module's docstring and in the stubs
        if not documented(obj):
            missing.append(name)
        if inspect.isclass(obj) and not issubclass(obj, BaseException):
            for member in dir(obj):
                if member.startswith("_"):
                    continue
                attribute = inspect.getattr_static(obj, member)
                if isinstance(attribute, (int, float)) or type(attribute) is obj:
                    continue  # class constants: BITS, ZERO, MAX, Units.YEARS
                if not documented(getattr(obj, member, attribute)) and not documented(attribute):
                    missing.append(f"{name}.{member}")
    assert not missing, f"no docstring: {missing}"


def test_the_unit_constants_are_listed_in_the_module_docstring():
    units = srkit.units
    for name in units.__all__:
        if name != "Units":
            assert isinstance(getattr(units, name), float) and name in units.__doc__


def test_help_renders():
    text = pydoc.render_doc(srkit, renderer=pydoc.plaintext)
    for expected in [
        "special relativity in 3+1 Minkowski spacetime",
        "class Segment",
        "class Worldline",
        "class Time512",
        "at_many(self, ts)",
        "recv(self, /, event_t, event_x, lo=-inf, hi=inf)",
        "Trip(d, a, cap=inf, stops=0, round_trip=False)",
        "reach(phi, a, budget, heading)",
        "NoIntersection",
    ]:
        assert expected in text, expected
    for module_name in MODULES[1:]:
        module = importlib.import_module(module_name)
        rendered = pydoc.render_doc(module, renderer=pydoc.plaintext)
        for name in module.__all__:
            assert name in rendered, f"{module_name}.{name}"
    plan_text = pydoc.render_doc(plan, renderer=pydoc.plaintext)
    assert "circ_piece(start, orbit, dtau=inf)" in plan_text and "phi_cap=inf" in plan_text


# ------------------------------------------------------------------------------------------- packaging


def test_version_and_packaging():
    cargo = (CRATE / "Cargo.toml").read_text(encoding="utf-8")
    version = re.search(r'^version = "([^"]+)"', cargo, re.M).group(1)
    assert srkit.__version__ == version
    tomllib = pytest.importorskip("tomllib")
    project = tomllib.loads((CRATE / "pyproject.toml").read_text(encoding="utf-8"))
    assert "Private :: Do Not Upload" in project["project"]["classifiers"]
    # the licence is not decided: none is declared, here or in Cargo.toml
    assert "license" not in project["project"] and "license-files" not in project["project"]
    assert not any(c.startswith("License ::") for c in project["project"]["classifiers"])
    assert not re.search(r"^license", cargo, re.M)
    assert project["tool"]["maturin"]["module-name"] == "srkit._srkit"
    assert 'features = ["abi3-py39"]' in cargo  # one wheel for CPython 3.9 and later
    assert "numpy" not in " ".join(project["project"].get("dependencies", []))  # NumPy stays optional


def readme_blocks() -> list[tuple[str, str]]:
    text = (CRATE / "README.md").read_text(encoding="utf-8")
    return [(m.group(1), m.group(2)) for m in re.finditer(r"^```(\w+)\n(.*?)^```$", text, re.S | re.M)]


def test_the_readme_examples_run_and_print_what_it_says():
    blocks = readme_blocks()
    ran = compared = 0
    for i, (language, body) in enumerate(blocks):
        if language != "python":
            continue
        if "import numpy" in body and importlib.util.find_spec("numpy") is None:
            ran += 1  # the NumPy example needs NumPy; the others must still run
            continue
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            exec(compile(body, f"README.md block {i}", "exec"), {"__name__": "__readme__"})
        ran += 1
        if i + 1 < len(blocks) and blocks[i + 1][0] == "text":
            assert out.getvalue() == blocks[i + 1][1]
            compared += 1
    assert ran == 3 and compared == 2


def test_the_readme_example_is_the_physics():
    """The numbers the README prints, from the closed forms."""
    g, d = srkit.units.G0_PER_YEAR, 4.2465
    fb = srkit.FlipBurn(d, g)
    assert f"{fb.t:.3f} {fb.tau:.3f} {math.tanh(fb.phi_peak):.4f}" == "5.873 3.542 0.9497"
    # seen from home at t = 3: the emission is on the boost, where t + x = 3
    trip = srkit.Trip(d, g)
    state = trip.state(trip.tau_of_v(3.0))
    assert f"{state.t:.3f} {state.x:.3f} {math.tanh(state.phi):.3f} {state.tau:.3f}" == "1.866 1.134 0.888 1.366"
