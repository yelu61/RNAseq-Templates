"""Unit tests for the bulk RNA-seq Skill backend router."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SCRIPT = (
    Path(__file__).parents[1]
    / "skills"
    / "bulk-rnaseq-analysis"
    / "scripts"
    / "backend_router.py"
)
SPEC = importlib.util.spec_from_file_location("backend_router", SCRIPT)
assert SPEC and SPEC.loader
ROUTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ROUTER)


class BackendRouterTests(unittest.TestCase):
    def test_local_deg_routes_to_rnaseq_templates(self) -> None:
        result = ROUTER.route("local", "raw-counts", {"deg", "enrichment"})
        self.assertEqual(result["backends"], ["rnaseq-templates"])
        self.assertFalse(result["warnings"])

    def test_tcga_source_always_uses_tcga_toolkit(self) -> None:
        result = ROUTER.route("tcga", "raw-counts", {"deg", "wgcna"})
        self.assertEqual(result["backends"], ["tcga-toolkit"])

    def test_specialized_cancer_task_uses_tcga_toolkit(self) -> None:
        result = ROUTER.route("local", "maf", {"mutation", "tmb"})
        self.assertEqual(result["backends"], ["tcga-toolkit"])

    def test_discovery_and_validation_request_uses_both(self) -> None:
        result = ROUTER.route(
            "local", "raw-counts", {"deg", "external-validation"}
        )
        self.assertEqual(result["backends"], ["rnaseq-templates", "tcga-toolkit"])

    def test_geo_survival_remains_generic(self) -> None:
        result = ROUTER.route("geo", "normalized", {"survival"})
        self.assertEqual(result["backends"], ["rnaseq-templates"])

    def test_invalid_scales_emit_warnings(self) -> None:
        result = ROUTER.route("local", "vst", {"deg", "tme"})
        self.assertEqual(len(result["warnings"]), 2)

    def test_unknown_request_requires_clarification(self) -> None:
        result = ROUTER.route("unknown", "unknown", set())
        self.assertTrue(result["needs_clarification"])
        self.assertIsNone(result["primary_backend"])

    def test_tcga_toolkit_env_accepts_toolkit_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            repo = Path(temp_dir) / "TCGA"
            toolkit = repo / "tcga_toolkit"
            for sentinel in ROUTER.TCGA_SENTINELS:
                path = repo / sentinel
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("test\n", encoding="utf-8")
            with mock.patch.dict(
                "os.environ", {"TCGA_TOOLKIT_ROOT": str(toolkit)}, clear=False
            ):
                result = ROUTER.discover_one(
                    "TCGA",
                    "TCGA_TOOLKIT_ROOT",
                    ROUTER.TCGA_SENTINELS,
                    Path(temp_dir),
                )
            self.assertTrue(result["found"])
            self.assertEqual(result["root"], str(repo.resolve()))


if __name__ == "__main__":
    unittest.main()
