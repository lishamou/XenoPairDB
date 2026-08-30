"""
__init__.py — Test runner infrastructure for the registry-driven audit.

Each gate module exports a `run(canonical_dir, facts)` function that returns a
list of TestResult. The audit orchestrator imports registered gates and aggregates
the results into the consistency_report.html.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Any

# Status constants
PASS = "PASS"
FAIL = "FAIL"
WARN = "WARN"


@dataclass
class TestResult:
    id: str
    name: str
    status: str
    message: str = ""
    actual: str = ""
    expected: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "status": self.status,
            "message": self.message,
            "actual": self.actual,
            "expected": self.expected,
        }


# Gate registry
_GATES: list[tuple[str, Callable, int]] = []


def register_gate(name: str, weight: int = 10):
    """Decorator that registers a gate's run function."""
    def decorator(fn):
        _GATES.append((name, fn, weight))
        return fn
    return decorator


def all_gates() -> list[tuple[str, Callable, int]]:
    return list(_GATES)
