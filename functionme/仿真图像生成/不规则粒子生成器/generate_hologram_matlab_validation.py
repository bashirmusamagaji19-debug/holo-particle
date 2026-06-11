#!/usr/bin/env python3
"""
Generate an inline digital hologram for MATLAB reconstruction validation.

Outputs (default):
  - hologram.bmp      : 8-bit grayscale hologram
  - object.bmp        : 8-bit binary object-plane map
  - ground_truth.csv  : particle metadata for error analysis in MATLAB
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy.fft import fft2, fftshift, ifft2, ifftshift


@dataclass
class Particle:
    particle_id: int
    x_pixel: float  # 0-based column index
    y_pixel: float  # 0-based row index
    radius_m: float
    z_position_m: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate BMP/JPG hologram + CSV ground truth for MATLAB validation."
    )
    parser.add_argument(
        "--img_size", "--N", "--size", dest="img_size", type=int, default=256, help="Image size NxN."
    )
    parser.add_argument(
        "--lambda_", "--wavelength", "--lambda", dest="lambda_", type=float, default=660e-9, help="Wavelength (m)."
    )
    parser.add_argument(
        "--pixel_size", "--dx", "--pixel-size", dest="pixel_size", type=float, default=2.2e-6, help="Pixel size dx=dy (m)."
    )
    parser.add_argument(
        "--particle_diam_range",
        "--diameter-range",
        dest="particle_diam_range",
        nargs=2,
        type=float,
        default=[20e-6, 30e-6],
        metavar=("D_MIN", "D_MAX"),
        help="Particle diameter range (m): [d_min, d_max].",
    )
    parser.add_argument(
        "--z_range",
        dest="z_range",
        nargs=2,
        type=float,
        default=[1800e-6, 4800e-6],
        metavar=("Z_MIN", "Z_MAX"),
        help="Axial range (m): [z_min, z_max].",
    )
    parser.add_argument("--z_step", type=float, default=100e-6, help="Axial layer spacing (m).")
    parser.add_argument("--num_particles", "--M", type=int, default=None, help="Fixed particle count per hologram.")
    parser.add_argument(
        "--num_particles_range",
        nargs=2,
        type=int,
        default=[30, 100],
        metavar=("N_MIN", "N_MAX"),
        help="Particle count range. Used when --num_particles is not provided.",
    )
    parser.add_argument("--R", type=float, default=1.0, help="Reference wave amplitude.")
    parser.add_argument(
        "--dc-mode",
        choices=["none", "subtract_mean", "subtract_1"],
        default="none",
        help="Optional DC operation before normalization.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Output directory.")
    parser.add_argument(
        "--holo-out", "--hologram-name", dest="hologram_name", type=str, default="hologram.bmp", help="Hologram filename."
    )
    parser.add_argument(
        "--obj-out", "--object-name", dest="object_name", type=str, default="object.bmp", help="Object filename."
    )
    parser.add_argument(
        "--gt-out", "--csv-name", dest="csv_name", type=str, default="ground_truth.csv", help="Ground truth CSV filename."
    )
    return parser.parse_args()


def draw_disk(mask: np.ndarray, cx: float, cy: float, radius_px: float, fill_value: float) -> None:
    ny, nx = mask.shape
    x = np.arange(nx, dtype=np.float64)
    y = np.arange(ny, dtype=np.float64)
    xx, yy = np.meshgrid(x, y, indexing="xy")
    disk = (xx - cx) ** 2 + (yy - cy) ** 2 <= radius_px ** 2
    mask[disk] = fill_value


def place_particles(
    n: int,
    num_particles: int,
    radius_min_px: float,
    radius_max_px: float,
    pixel_size: float,
    z_levels: np.ndarray,
    rng: np.random.Generator,
    max_attempts: int = 20000,
) -> list[Particle]:
    particles: list[Particle] = []
    attempts = 0
    while len(particles) < num_particles and attempts < max_attempts:
        attempts += 1
        r_px = float(rng.uniform(radius_min_px, radius_max_px))
        margin = r_px + 2.0
        if margin >= (n / 2.0):
            raise ValueError("Particle radius too large for image size.")

        x = float(rng.uniform(margin, n - 1 - margin))
        y = float(rng.uniform(margin, n - 1 - margin))

        no_overlap = True
        for p in particles:
            p_r_px = p.radius_m / pixel_size
            dist = np.hypot(x - p.x_pixel, y - p.y_pixel)
            if dist < (r_px + p_r_px + 2.0):
                no_overlap = False
                break

        if no_overlap:
            particles.append(
                Particle(
                    particle_id=len(particles) + 1,
                    x_pixel=x,
                    y_pixel=y,
                    radius_m=r_px * pixel_size,
                    z_position_m=float(rng.choice(z_levels)),
                )
            )

    if len(particles) < num_particles:
        raise RuntimeError(
            f"Could only place {len(particles)} particles out of {num_particles}. "
            "Try reducing particle count/size or increasing image size."
        )
    return particles


def build_object_plane(n: int, particles: list[Particle], pixel_size: float) -> np.ndarray:
    # Amplitude object for MATLAB validation: background=0 (dark), particles=1 (bright).
    obj = np.zeros((n, n), dtype=np.float64)
    for p in particles:
        radius_px = p.radius_m / pixel_size
        draw_disk(obj, p.x_pixel, p.y_pixel, radius_px, fill_value=1.0)
    return obj


def build_depth_planes(n: int, particles: list[Particle], pixel_size: float) -> dict[float, np.ndarray]:
    planes: dict[float, np.ndarray] = {}
    for p in particles:
        if p.z_position_m not in planes:
            planes[p.z_position_m] = np.zeros((n, n), dtype=np.float64)
        radius_px = p.radius_m / pixel_size
        draw_disk(planes[p.z_position_m], p.x_pixel, p.y_pixel, radius_px, fill_value=1.0)
    return planes


def angular_spectrum_propagation(obj: np.ndarray, wavelength: float, pixel_size: float, z: float) -> np.ndarray:
    n_rows, n_cols = obj.shape
    k = 2.0 * np.pi / wavelength

    dfx = 1.0 / (n_cols * pixel_size)
    dfy = 1.0 / (n_rows * pixel_size)
    fx = (np.arange(n_cols) - n_cols / 2.0) * dfx
    fy = (np.arange(n_rows) - n_rows / 2.0) * dfy
    fxx, fyy = np.meshgrid(fx, fy, indexing="xy")

    root_arg = 1.0 - (wavelength * fxx) ** 2 - (wavelength * fyy) ** 2
    h = np.zeros_like(root_arg, dtype=np.complex128)
    valid = root_arg >= 0.0
    h[valid] = np.exp(1j * k * z * np.sqrt(root_arg[valid]))

    obj_f = fftshift(fft2(ifftshift(obj)))
    u_f = obj_f * h
    u = fftshift(ifft2(ifftshift(u_f)))
    return u


def quantize_to_uint8(intensity: np.ndarray, dc_mode: str) -> np.ndarray:
    i_proc = intensity.astype(np.float64)
    if dc_mode == "subtract_mean":
        i_proc = i_proc - np.mean(i_proc)
    elif dc_mode == "subtract_1":
        i_proc = i_proc - 1.0

    i_min = np.min(i_proc)
    i_max = np.max(i_proc)
    if np.isclose(i_max, i_min):
        i_norm = np.zeros_like(i_proc)
    else:
        i_norm = (i_proc - i_min) / (i_max - i_min)

    return np.uint8(np.round(i_norm * 255.0))


def save_gray_image(arr_uint8: np.ndarray, path: Path) -> None:
    Image.fromarray(arr_uint8, mode="L").save(path)


def save_ground_truth_csv(particles: list[Particle], path: Path) -> None:
    df = pd.DataFrame(
        {
            "Particle_ID": [p.particle_id for p in particles],
            "x_pixel": [p.x_pixel for p in particles],  # 0-based, x=column
            "y_pixel": [p.y_pixel for p in particles],  # 0-based, y=row
            "radius_m": [p.radius_m for p in particles],
            "z_position_m": [p.z_position_m for p in particles],
        }
    )
    df.to_csv(path, index=False)


if __name__ == "__main__":
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    diameter_min, diameter_max = args.particle_diam_range
    z_min, z_max = args.z_range
    if z_min > z_max:
        raise ValueError("Invalid z_range: z_min must be <= z_max.")
    if args.z_step <= 0:
        raise ValueError("z_step must be positive.")

    z_levels = np.arange(z_min, z_max + 0.5 * args.z_step, args.z_step, dtype=np.float64)
    if z_levels.size == 0:
        raise ValueError("No axial layers generated. Check z_range and z_step.")

    if args.num_particles is None:
        n_min, n_max = args.num_particles_range
        if n_min <= 0 or n_max < n_min:
            raise ValueError("Invalid num_particles_range.")
        num_particles = int(rng.integers(n_min, n_max + 1))
    else:
        if args.num_particles <= 0:
            raise ValueError("num_particles must be positive.")
        num_particles = args.num_particles

    radius_min_px = (diameter_min * 0.5) / args.pixel_size
    radius_max_px = (diameter_max * 0.5) / args.pixel_size
    if radius_min_px <= 0 or radius_max_px <= 0 or radius_max_px < radius_min_px:
        raise ValueError("Invalid diameter range.")

    particles = place_particles(
        n=args.img_size,
        num_particles=num_particles,
        radius_min_px=radius_min_px,
        radius_max_px=radius_max_px,
        pixel_size=args.pixel_size,
        z_levels=z_levels,
        rng=rng,
    )

    obj = build_object_plane(args.img_size, particles, args.pixel_size)
    depth_planes = build_depth_planes(args.img_size, particles, args.pixel_size)
    u_obj = np.zeros((args.img_size, args.img_size), dtype=np.complex128)
    for z in sorted(depth_planes.keys()):
        u_obj += angular_spectrum_propagation(depth_planes[z], args.lambda_, args.pixel_size, z)

    reference = complex(args.R, 0.0)
    intensity = np.abs(reference + u_obj) ** 2
    hologram_uint8 = quantize_to_uint8(intensity, dc_mode=args.dc_mode)

    # Binary object map for visual comparison.
    object_uint8 = np.uint8(np.round(obj * 255.0))

    save_gray_image(hologram_uint8, args.output_dir / args.hologram_name)
    save_gray_image(object_uint8, args.output_dir / args.object_name)
    save_ground_truth_csv(particles, args.output_dir / args.csv_name)

    print("Done.")
    print(f"Hologram: {args.output_dir / args.hologram_name}")
    print(f"Object  : {args.output_dir / args.object_name}")
    print(f"CSV     : {args.output_dir / args.csv_name}")
    print("CSV coordinates are 0-based pixel indices: x=column, y=row.")
