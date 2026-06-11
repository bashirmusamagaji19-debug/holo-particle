import os
import tempfile
import unittest

from generate_hologram_matlab_validation_new import generate_hologram


class GenerateHologramOutputTests(unittest.TestCase):
    def test_no_noise_creates_numbered_folder_with_three_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "0001")

            generate_hologram(
                N=32,
                num_particles=2,
                z_step=40e-6,
                snr_db=0,
                output_dir=output_dir,
                prefix="0001_",
            )

            self.assertTrue(os.path.isdir(output_dir))
            self.assertTrue(os.path.isfile(os.path.join(output_dir, "0001_hologram.bmp")))
            self.assertTrue(os.path.isfile(os.path.join(output_dir, "0001_object_reference.bmp")))
            self.assertTrue(os.path.isfile(os.path.join(output_dir, "0001_particle_ground_truth.csv")))
            self.assertFalse(os.path.exists(os.path.join(output_dir, "0001_hologram_with_noise.bmp")))
            self.assertFalse(os.path.exists(os.path.join(output_dir, "0001_background_noise.bmp")))

    def test_noise_creates_extra_noisy_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "0001")

            generate_hologram(
                N=32,
                num_particles=2,
                z_step=40e-6,
                snr_db=20,
                output_dir=output_dir,
                prefix="0001_",
            )

            self.assertTrue(os.path.isfile(os.path.join(output_dir, "0001_hologram.bmp")))
            self.assertTrue(os.path.isfile(os.path.join(output_dir, "0001_object_reference.bmp")))
            self.assertTrue(os.path.isfile(os.path.join(output_dir, "0001_particle_ground_truth.csv")))
            self.assertTrue(os.path.isfile(os.path.join(output_dir, "0001_hologram_with_noise.bmp")))
            self.assertTrue(os.path.isfile(os.path.join(output_dir, "0001_background_noise.bmp")))


if __name__ == "__main__":
    unittest.main()
