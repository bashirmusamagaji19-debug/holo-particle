import json
import tempfile
import unittest
from pathlib import Path

from generate_hologram_shape_classifier import generate_hologram


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


if __name__ == "__main__":
    unittest.main()
