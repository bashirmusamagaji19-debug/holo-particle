import csv
import json
import tempfile
import typing
import unittest
from pathlib import Path

import numpy as np

from generate_hologram_shape_classifier import (
    Particle,
    PolyhedronShape,
    _adaptive_voxel_um,
    generate_hologram,
)


class ShapeClassifierParameterTests(unittest.TestCase):
    def test_custom_polyhedron_parameters_are_recorded_without_preset_override(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            generate_hologram(
                N=96,
                n_circular=0,
                n_irregular=1,
                radius_um_min=5.0,
                radius_um_max=7.0,
                poly_radial_jitter=0.18,
                poly_spike_fraction=0.05,
                poly_indent_fraction=0.08,
                poly_vertex_min=9,
                poly_vertex_max=11,
                fragmentation="custom",
                output_dir=tmpdir,
                prefix="case_",
                seed=123,
            )

            metadata = json.loads((Path(tmpdir) / "case_metadata.json").read_text(encoding="utf-8"))

        self.assertEqual(metadata["fragmentation"], "custom")
        self.assertEqual(metadata["poly_radial_jitter"], 0.18)
        self.assertEqual(metadata["poly_spike_fraction"], 0.05)
        self.assertEqual(metadata["poly_indent_fraction"], 0.08)
        self.assertEqual(metadata["poly_vertex_min"], 9)
        self.assertEqual(metadata["poly_vertex_max"], 11)

    def test_polyhedron_fragment_parameters_are_recorded_in_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            generate_hologram(
                N=96,
                n_circular=0,
                n_irregular=1,
                radius_um_min=5.0,
                radius_um_max=7.0,
                fragmentation="severe",
                output_dir=tmpdir,
                prefix="case_",
                seed=7,
            )

            metadata = json.loads((Path(tmpdir) / "case_metadata.json").read_text(encoding="utf-8"))

        self.assertEqual(metadata["poly_num_fragments"], 8)
        self.assertEqual(metadata["poly_fragment_overlap"], 0.35)

    def test_polyhedron_voxel_grid_contains_severe_fragment_union(self):
        for seed in range(1, 6):
            rng = np.random.default_rng(seed)
            shape = PolyhedronShape(
                target_radius_um=50.0,
                num_vertices=45,
                aspect_ratio=(1.0, 1.0, 1.0),
                radial_jitter=0.55,
                spike_fraction=0.20,
                indent_fraction=0.30,
                num_fragments=8,
                fragment_overlap=0.35,
            )
            shape.initialize(rng, voxel_um=_adaptive_voxel_um(50.0))

            thickness = shape._thickness_map
            touches_edge = (
                thickness[0, :].any()
                or thickness[-1, :].any()
                or thickness[:, 0].any()
                or thickness[:, -1].any()
            )
            self.assertFalse(touches_edge, f"seed {seed} was clipped at the voxel boundary")

    def test_polyhedron_fragment_count_supports_finer_severe_shapes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            generate_hologram(
                N=96,
                n_circular=0,
                n_irregular=1,
                radius_um_min=5.0,
                radius_um_max=7.0,
                fragmentation="severe",
                output_dir=tmpdir,
                prefix="case_",
                seed=7,
            )

            with open(Path(tmpdir) / "case_ground_truth.csv", newline="", encoding="utf-8") as f:
                row = next(csv.DictReader(f))

        num_fragments = int(row["num_fragments"])
        self.assertGreaterEqual(num_fragments, 7)
        self.assertLessEqual(num_fragments, 8)

    def test_particle_type_hints_resolve(self):
        hints = typing.get_type_hints(Particle)

        self.assertIn("shape", hints)


if __name__ == "__main__":
    unittest.main()
