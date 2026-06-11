#!/usr/bin/env python3
"""
generate_hologram_shape_classifier.py
======================================
Generate synthetic inline digital holograms for particle shape classification.

Supports three particle types:
  - sphere     : circular particles (standard spheres)
  - aggregate  : irregular particles (aggregated sphere clusters, like PTFE powder)
  - polyhedron : irregular particles (convex polyhedra, like ice crystal fragments)

Physics model — Complex transmittance thin screen + angular spectrum propagation:
  t(x,y) = A(x,y) * exp(i * k * dn * h(x,y))

  h(x,y) = particle thickness along optical axis at lateral position (x,y)
  A(x,y) = amplitude attenuation (0~1, with smooth edge transition)
  dn      = n_particle - n_medium (refractive index contrast)
  k       = 2*pi / lambda

  The scattered field (t - 1) at each particle plane is propagated to the sensor
  via angular spectrum method. Fields from all particles are superposed, then
  interfered with the unit-amplitude reference wave to produce the inline hologram:
      I = |1 - sum(U_scatter_i)|^2

Outputs per hologram (into numbered folder):
  {prefix}hologram.bmp          : 8-bit grayscale hologram
  {prefix}hologram.mat          : float64 exact hologram (for MATLAB verification)
  {prefix}object_reference.bmp  : ground-truth object-plane amplitude projection
  {prefix}ground_truth.csv      : particle metadata

Usage:
  python generate_hologram_shape_classifier.py --output_dir ./output/0001 --seed 42

  # Batch generation:
  for i in {1..50}; do
      python generate_hologram_shape_classifier.py --output_dir ./output/$(printf "%04d" $i) --seed $i
  done
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from PIL import Image
from scipy.fft import fft2, fftshift, ifft2, ifftshift
from scipy.ndimage import gaussian_filter
from scipy.io import savemat
from scipy.spatial import ConvexHull

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SPEED_OF_LIGHT = 299792458.0  # m/s


# ===================================================================
# Shape generators
# ===================================================================

@dataclass
class SphereShape:
    """Circular particle — standard homogeneous sphere."""
    radius_um: float

    def get_thickness_map(self, x_um: np.ndarray, y_um: np.ndarray) -> np.ndarray:
        """h(x,y) = 2 * sqrt(R^2 - (x-x0)^2 - (y-y0)^2) for points inside the sphere."""
        r_sq = x_um**2 + y_um**2
        inside = r_sq < self.radius_um**2
        h = np.zeros_like(x_um)
        h[inside] = 2.0 * np.sqrt(np.maximum(0, self.radius_um**2 - r_sq[inside]))
        return h

    def compute_volume(self) -> float:
        return (4.0 / 3.0) * np.pi * self.radius_um**3

    @property
    def shape_type(self) -> str:
        return "sphere"

    @property
    def shape_params(self) -> dict:
        return {"radius_um": self.radius_um}


def _generate_aggregate_sub_spheres(
    target_radius_um: float,
    num_sub_spheres: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate an irregular cluster of overlapping sub-spheres.

    Sub-spheres are biased toward the surface (not the center) to avoid
    the cluster looking like a smooth sphere. Fewer, larger sub-spheres
    with looser connectivity create more irregular silhouettes.

    Returns (centers_um, radii_um) — 3D coordinates relative to centroid.
    """
    centers = np.zeros((num_sub_spheres, 3))
    radii = np.zeros(num_sub_spheres)

    # Seed at a random off-center position
    r0 = target_radius_um * rng.uniform(0.3, 0.7)
    phi = rng.uniform(0, 2 * np.pi)
    theta = np.arccos(rng.uniform(-1, 1))
    centers[0] = [r0 * np.sin(theta) * np.cos(phi),
                  r0 * np.sin(theta) * np.sin(phi),
                  r0 * np.cos(theta)]
    radii[0] = rng.uniform(0.30, 0.55) * target_radius_um

    for i in range(1, num_sub_spheres):
        placed = False
        for _ in range(200):
            # Push sub-spheres toward the outer shell (0.5–1.0 of target_radius)
            # r^3 distribution gives uniform-in-volume; raise to <1 for surface bias
            phi = rng.uniform(0, 2 * np.pi)
            cos_theta = rng.uniform(-1, 1)
            theta = np.arccos(cos_theta)
            radial_frac = rng.uniform(0.0, 1.0) ** (1.0 / 4.0)  # bias toward surface
            r = target_radius_um * (0.3 + 0.7 * radial_frac)

            cx = r * np.sin(theta) * np.cos(phi)
            cy = r * np.sin(theta) * np.sin(phi)
            cz = r * np.cos(theta)

            # Larger, more varied sub-spheres
            sub_r = rng.uniform(0.25, 0.60) * target_radius_um

            # Looser connectivity: only need to touch (not deeply overlap)
            dists = np.sqrt(
                (centers[:i, 0] - cx) ** 2
                + (centers[:i, 1] - cy) ** 2
                + (centers[:i, 2] - cz) ** 2
            )
            if np.any(dists < (radii[:i] + sub_r) * 1.05):
                centers[i] = [cx, cy, cz]
                radii[i] = sub_r
                placed = True
                break

        if not placed:
            parent = rng.integers(0, i)
            phi = rng.uniform(0, 2 * np.pi)
            cos_theta = rng.uniform(-1, 1)
            theta = np.arccos(cos_theta)
            sub_r = rng.uniform(0.25, 0.55) * target_radius_um
            # Place at surface of parent (minimal overlap)
            direction = np.array(
                [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
            )
            direction = direction / np.linalg.norm(direction)
            centers[i] = centers[parent] + direction * (radii[parent] + sub_r) * rng.uniform(0.85, 1.0)
            radii[i] = sub_r

    # Recenter so centroid is at origin
    centroid = np.mean(centers, axis=0)
    centers -= centroid

    return centers, radii


def _adaptive_voxel_um(target_radius_um: float) -> float:
    """Choose voxel size based on particle radius to keep meshgrid manageable.

    At 3.45 um/pixel with 3x oversample, the effective query resolution is ~1.15 um.
    We target ~1-3 um voxels, scaling with radius to cap n_voxels at ~60-80 per dim.
    """
    return max(1.5, target_radius_um / 50.0)


def _voxelize_sphere_cluster(
    centers_um: np.ndarray,
    radii_um: np.ndarray,
    voxel_um: float,
    padding_um: float = 5.0,
) -> tuple[np.ndarray, float, np.ndarray]:
    """Voxelize a cluster of spheres via slice-by-slice 2D processing.

    Avoids 3D meshgrid — only a single 2D (X,Y) meshgrid is allocated.
    Memory: O(N^2) instead of O(N^3).

    Returns:
        thickness_map_um: 2D array (ny, nx) — thickness along z at each (x,y)
        volume_um3: float — total occupied volume
        voxel_centers: 1D array of voxel center coordinates in um
    """
    extent = np.max(np.abs(centers_um)) + np.max(radii_um) + padding_um
    n_voxels = int(np.ceil(2 * extent / voxel_um))
    if n_voxels % 2 == 0:
        n_voxels += 1  # odd so origin is a voxel center

    voxel_centers = (np.arange(n_voxels) - n_voxels // 2) * voxel_um

    # Single 2D meshgrid reused for every z-slice
    X2D, Y2D = np.meshgrid(voxel_centers, voxel_centers, indexing="ij")

    thickness_map_um = np.zeros((n_voxels, n_voxels), dtype=np.float64)
    total_occupied = 0

    for z in voxel_centers:
        occupied_2d = np.zeros((n_voxels, n_voxels), dtype=bool)
        for c, r in zip(centers_um, radii_um):
            dz_sq = (z - c[2]) ** 2
            r_xy_sq = r**2 - dz_sq
            if r_xy_sq <= 0:
                continue
            dist_sq = (X2D - c[0]) ** 2 + (Y2D - c[1]) ** 2
            occupied_2d |= dist_sq < r_xy_sq
        thickness_map_um += occupied_2d.astype(np.float64)
        total_occupied += int(occupied_2d.sum())

    thickness_map_um *= voxel_um
    volume_um3 = total_occupied * voxel_um**3

    return thickness_map_um, volume_um3, voxel_centers


@dataclass
class AggregateShape:
    """Irregular particle — cluster of overlapping sub-spheres (like PTFE powder)."""
    target_radius_um: float
    num_sub_spheres: int
    _sub_centers: np.ndarray = field(default=None, repr=False)
    _sub_radii: np.ndarray = field(default=None, repr=False)
    _volume_um3: float = field(default=None, repr=False)
    _thickness_map: np.ndarray = field(default=None, repr=False)
    _voxel_centers: np.ndarray = field(default=None, repr=False)

    def __post_init__(self):
        # Will be initialized lazily
        pass

    def initialize(self, rng: np.random.Generator, voxel_um: float = 0.5):
        self._sub_centers, self._sub_radii = _generate_aggregate_sub_spheres(
            self.target_radius_um, self.num_sub_spheres, rng
        )
        self._thickness_map, self._volume_um3, self._voxel_centers = (
            _voxelize_sphere_cluster(self._sub_centers, self._sub_radii, voxel_um)
        )

    def get_thickness_map(self, x_um: np.ndarray, y_um: np.ndarray) -> np.ndarray:
        """Interpolate precomputed voxel thickness onto arbitrary query grid."""
        from scipy.interpolate import RegularGridInterpolator

        if self._thickness_map is None:
            raise RuntimeError("AggregateShape not initialized — call initialize(rng)")

        interp = RegularGridInterpolator(
            (self._voxel_centers, self._voxel_centers),
            self._thickness_map,
            bounds_error=False,
            fill_value=0.0,
            method="linear",
        )
        pts = np.stack([x_um.ravel(), y_um.ravel()], axis=-1)
        h = interp(pts).reshape(x_um.shape)
        return np.maximum(h, 0.0)

    def compute_volume(self) -> float:
        if self._volume_um3 is None:
            raise RuntimeError("AggregateShape not initialized")
        return self._volume_um3

    @property
    def shape_type(self) -> str:
        return "aggregate"

    @property
    def shape_params(self) -> dict:
        return {
            "target_radius_um": self.target_radius_um,
            "num_sub_spheres": self.num_sub_spheres,
        }


@dataclass
class PolyhedronShape:
    """Irregular particle — convex polyhedron (like ice crystal fragment).

    A random 3D convex hull is generated from scattered points, optionally
    stretched along one axis to create plate-like or columnar habits.
    """
    target_radius_um: float
    num_vertices: int
    aspect_ratio: tuple[float, float, float] = (1.0, 1.0, 1.0)
    _hull_vertices: np.ndarray = field(default=None, repr=False)
    _volume_um3: float = field(default=None, repr=False)
    _thickness_map: np.ndarray = field(default=None, repr=False)
    _voxel_centers: np.ndarray = field(default=None, repr=False)

    def initialize(self, rng: np.random.Generator, voxel_um: float = 0.5):
        # Generate random points from a stretched 3D normal distribution.
        # This produces a cloud with natural anisotropy, not a spherical shell.
        # When fed into ConvexHull, the result is a faceted polyhedron with
        # sharp edges and flat faces — ideal for ice-crystal-like fragments.
        ax, ay, az = self.aspect_ratio
        n_total = self.num_vertices + self.num_vertices // 2  # fewer interior points

        # Surface-biased points: use random radial scaling so some points
        # are near the surface and some are deeper, creating flat faces
        points = rng.normal(0, 1, (n_total, 3))

        # Radial rescaling: push points toward the surface to form sharp faces
        norms = np.linalg.norm(points, axis=1, keepdims=True)
        norms = np.where(norms < 1e-9, 1e-9, norms)
        # Mix of surface and interior: 70% pushed to surface, 30% scattered
        surface_mask = rng.uniform(0, 1, size=(n_total, 1)) < 0.7
        r_uniform = rng.uniform(0.6, 1.0, size=(n_total, 1))
        scaling = np.where(surface_mask, 1.0, r_uniform)
        points = points / norms * self.target_radius_um * scaling

        # Apply aspect ratio stretching
        points[:, 0] *= ax
        points[:, 1] *= ay
        points[:, 2] *= az

        try:
            hull = ConvexHull(points)
            self._hull_vertices = points[hull.vertices]
            self._volume_um3 = hull.volume
        except Exception:
            self._hull_vertices = points
            self._volume_um3 = (4.0 / 3.0) * np.pi * self.target_radius_um**3

        # Voxelize the convex hull
        self._thickness_map, self._volume_um3, self._voxel_centers = (
            _voxelize_convex_hull(self._hull_vertices, voxel_um)
        )

    def get_thickness_map(self, x_um: np.ndarray, y_um: np.ndarray) -> np.ndarray:
        from scipy.interpolate import RegularGridInterpolator

        if self._thickness_map is None:
            raise RuntimeError("PolyhedronShape not initialized — call initialize(rng)")

        interp = RegularGridInterpolator(
            (self._voxel_centers, self._voxel_centers),
            self._thickness_map,
            bounds_error=False,
            fill_value=0.0,
            method="linear",
        )
        pts = np.stack([x_um.ravel(), y_um.ravel()], axis=-1)
        h = interp(pts).reshape(x_um.shape)
        return np.maximum(h, 0.0)

    def compute_volume(self) -> float:
        if self._volume_um3 is None:
            raise RuntimeError("PolyhedronShape not initialized")
        return self._volume_um3

    @property
    def shape_type(self) -> str:
        return "polyhedron"

    @property
    def shape_params(self) -> dict:
        return {
            "target_radius_um": self.target_radius_um,
            "num_vertices": self.num_vertices,
            "aspect_ratio": list(self.aspect_ratio),
        }



def _voxelize_convex_hull(
    hull_vertices: np.ndarray,
    voxel_um: float,
    padding_um: float = 5.0,
) -> tuple[np.ndarray, float, np.ndarray]:
    """Voxelize a convex hull via slice-by-slice half-space testing.

    Avoids 3D meshgrid — for each z-slice, tests all (x,y) against the
    hull's face equations using vectorized 2D operations.
    Memory: O(N^2) instead of O(N^3).

    Returns:
        thickness_map_um: 2D array — thickness along z at each (x,y)
        volume_um3: float
        voxel_centers: 1D array
    """
    extent = np.max(np.abs(hull_vertices)) + padding_um
    n_voxels = int(np.ceil(2 * extent / voxel_um))
    if n_voxels % 2 == 0:
        n_voxels += 1

    voxel_centers = (np.arange(n_voxels) - n_voxels // 2) * voxel_um

    # Precompute half-space equations
    hull = ConvexHull(hull_vertices)
    eq = hull.equations  # (F, 4): A*x + B*y + C*z + D <= 0 inside

    # Single 2D meshgrid reused for every z-slice
    X2D, Y2D = np.meshgrid(voxel_centers, voxel_centers, indexing="ij")

    thickness_map_um = np.zeros((n_voxels, n_voxels), dtype=np.float64)
    total_occupied = 0

    for z in voxel_centers:
        occupied_2d = np.ones((n_voxels, n_voxels), dtype=bool)
        for face_eq in eq:
            A, B, C, D = face_eq
            lhs = A * X2D + B * Y2D + C * z + D
            occupied_2d &= lhs <= 1e-9
        thickness_map_um += occupied_2d.astype(np.float64)
        total_occupied += int(occupied_2d.sum())

    thickness_map_um *= voxel_um
    volume_um3 = total_occupied * voxel_um**3

    return thickness_map_um, volume_um3, voxel_centers


# ===================================================================
# Particle placement
# ===================================================================

@dataclass
class Particle:
    particle_id: int
    shape_type: str  # "sphere" | "aggregate" | "polyhedron"
    x_pixel: float  # 0-based column index in image
    y_pixel: float  # 0-based row index in image
    z_position_m: float
    shape: SphereShape | AggregateShape | Polyhedron
    volume_um3: float = 0.0
    equivalent_diameter_um: float = 0.0


def _compute_equivalent_diameter(volume_um3: float) -> float:
    """Equivalent spherical diameter from volume."""
    return 2.0 * (3.0 * volume_um3 / (4.0 * np.pi)) ** (1.0 / 3.0)


def place_particles(
    N: int,
    pixel_size_m: float,
    n_circular: int,
    n_irregular: int,
    radius_um_range: tuple[float, float],
    z_range_m: tuple[float, float],
    rng: np.random.Generator,
    irregular_type: str = "aggregate",
    max_attempts: int = 5000,
) -> list[Particle]:
    """Place particles on the image plane without overlap.

    Irregular particles are split equally between aggregate and polyhedron types.
    """
    total = n_circular + n_irregular
    pixel_size_um = pixel_size_m * 1e6
    z_min, z_max = z_range_m
    r_min_um, r_max_um = radius_um_range

    particles: list[Particle] = []
    placed_positions_um: list[tuple[float, float, float]] = []  # (x, y, radius)

    # Build type list
    type_list = []
    for _ in range(n_circular):
        type_list.append("sphere")
    half_irregular = n_irregular // 2
    for _ in range(half_irregular):
        type_list.append("aggregate")
    for _ in range(n_irregular - half_irregular):
        type_list.append("polyhedron")
    rng.shuffle(type_list)

    particle_id = 0
    for shape_type in type_list:
        placed = False
        for _ in range(max_attempts):
            radius_um = rng.uniform(r_min_um, r_max_um)
            radius_px = radius_um / pixel_size_um
            margin = radius_px + 5.0

            if margin >= N / 2.0:
                continue

            x_px = rng.uniform(margin, N - 1 - margin)
            y_px = rng.uniform(margin, N - 1 - margin)
            x_um = x_px * pixel_size_um
            y_um = y_px * pixel_size_um

            # Check overlap in image plane
            overlap = False
            for px, py, pr in placed_positions_um:
                dist = np.hypot(x_um - px, y_um - py)
                if dist < (radius_um + pr + 2.0 * pixel_size_um):
                    overlap = True
                    break

            if not overlap:
                z_m = rng.uniform(z_min, z_max)

                # Create shape object
                if shape_type == "sphere":
                    shape = SphereShape(radius_um=radius_um)
                elif shape_type == "aggregate":
                    n_sub = rng.integers(5, 16)
                    shape = AggregateShape(target_radius_um=radius_um, num_sub_spheres=n_sub)
                    shape.initialize(rng, voxel_um=_adaptive_voxel_um(radius_um))
                else:  # polyhedron
                    n_vertices = rng.integers(6, 14)
                    # Stronger, more varied aspect ratios for plate/column/irregular habits
                    aspect = (
                        rng.uniform(0.25, 1.0),
                        rng.uniform(0.25, 1.0),
                        rng.uniform(0.15, 1.0),
                    )
                    shape = PolyhedronShape(
                        target_radius_um=radius_um,
                        num_vertices=n_vertices,
                        aspect_ratio=aspect,
                    )
                    shape.initialize(rng, voxel_um=_adaptive_voxel_um(radius_um))

                volume = shape.compute_volume()
                eq_diam = _compute_equivalent_diameter(volume)

                particle_id += 1
                particle = Particle(
                    particle_id=particle_id,
                    shape_type=shape_type,
                    x_pixel=x_px,
                    y_pixel=y_px,
                    z_position_m=z_m,
                    shape=shape,
                    volume_um3=volume,
                    equivalent_diameter_um=eq_diam,
                )
                particles.append(particle)
                placed_positions_um.append((x_um, y_um, radius_um))
                placed = True
                break

        if not placed:
            print(f"  [警告] 无法放置粒子 #{particle_id+1} ({shape_type})，跳过")

    # Renumber sequentially
    for i, p in enumerate(particles):
        p.particle_id = i + 1

    return particles


# ===================================================================
# Complex transmittance computation
# ===================================================================

def compute_complex_transmittance(
    shape: SphereShape | AggregateShape | Polyhedron,
    x_pixel: float,
    y_pixel: float,
    N: int,
    pixel_size_um: float,
    wavelength_m: float,
    dn: float,
    attenuation: float,
    edge_sigma_px: float,
    roughness: float,
    rng: np.random.Generator,
    oversample: int = 3,
) -> np.ndarray:
    """Compute complex transmittance t(x,y) for a single particle.

    t(x,y) = A(x,y) * exp(i * k * dn * h(x,y))

    Supersampling: h(x,y) and A(x,y) are computed at sub-pixel resolution
    (oversample * pixel_size_um) and then averaged down to the image grid.
    This prevents aliasing artifacts (e.g. spheres looking like rounded squares)
    and preserves fine shape features.

    Returns:
        t: complex128 array of shape (N, N), t=1 outside particle region
    """
    k = 2.0 * np.pi / wavelength_m

    # Determine bounding box in image pixels (add margin for blur + oversample)
    radius_um = _estimate_radius_um(shape)
    radius_px = radius_um / pixel_size_um + edge_sigma_px * 2 + 6
    x1 = max(0, int(np.floor(x_pixel - radius_px)))
    x2 = min(N, int(np.ceil(x_pixel + radius_px)) + 1)
    y1 = max(0, int(np.floor(y_pixel - radius_px)))
    y2 = min(N, int(np.ceil(y_pixel + radius_px)) + 1)

    nx_img = x2 - x1
    ny_img = y2 - y1

    # ---- Fine grid (oversampled) ----
    fine_pixel_um = pixel_size_um / oversample
    nx_fine = nx_img * oversample
    ny_fine = ny_img * oversample

    # Fine-grid pixel centers in image coordinates
    x_fine_px = np.linspace(x1 + 0.5 / oversample, x2 - 0.5 / oversample, nx_fine)
    y_fine_px = np.linspace(y1 + 0.5 / oversample, y2 - 0.5 / oversample, ny_fine)
    xx_fine, yy_fine = np.meshgrid(x_fine_px, y_fine_px)

    # Fine-grid coordinates in um relative to particle center
    x_um = (xx_fine - x_pixel) * pixel_size_um
    y_um = (yy_fine - y_pixel) * pixel_size_um

    # Get thickness map at fine resolution
    h_fine = shape.get_thickness_map(x_um, y_um)

    # Binary mask from thickness
    mask_fine = (h_fine > 1e-9).astype(np.float64)

    # ---- Downsample to image resolution ----
    # Reshape and average over oversample x oversample blocks
    h_img = h_fine.reshape(ny_img, oversample, nx_img, oversample).mean(axis=(1, 3))
    mask_img = mask_fine.reshape(ny_img, oversample, nx_img, oversample).mean(axis=(1, 3))

    # Smooth edge transition at image resolution
    mask_smooth = gaussian_filter(mask_img, sigma=edge_sigma_px)

    # Amplitude: attenuation inside particle, 1 outside
    A = 1.0 - mask_smooth * (1.0 - attenuation)

    # Add surface roughness for irregular particles
    if roughness > 0 and shape.shape_type != "sphere":
        roughness_map = rng.normal(0, roughness, A.shape)
        roughness_map = gaussian_filter(roughness_map, sigma=1.0)
        A = A + roughness_map * mask_smooth
        A = np.clip(A, 0.05, 1.0)

    # Complex transmittance — h_img is in um, convert to meters
    phase = k * dn * (h_img * 1e-6)
    t_local = A * np.exp(1j * phase)

    # Place onto full grid (default t=1 everywhere)
    t_full = np.ones((N, N), dtype=np.complex128)
    t_full[y1:y2, x1:x2] = t_local

    return t_full


def _estimate_radius_um(shape) -> float:
    """Estimate the lateral extent of a shape in um."""
    if isinstance(shape, SphereShape):
        return shape.radius_um
    elif isinstance(shape, AggregateShape):
        return shape.target_radius_um * 1.15
    elif isinstance(shape, PolyhedronShape):
        return shape.target_radius_um * max(shape.aspect_ratio[:2]) * 1.2
    return 50.0


# ===================================================================
# Angular spectrum propagation
# ===================================================================

def angular_spectrum_propagation(
    field: np.ndarray,
    wavelength_m: float,
    pixel_size_m: float,
    z_m: float,
) -> np.ndarray:
    """Propagate a complex field by distance z using the angular spectrum method.

    Args:
        field: complex field at source plane, shape (Ny, Nx)
        wavelength_m: wavelength in meters
        pixel_size_m: pixel size in meters (assumes square pixels)
        z_m: propagation distance in meters (positive = forward)

    Returns:
        complex field at destination plane
    """
    ny, nx = field.shape
    k = 2.0 * np.pi / wavelength_m

    # Frequency grid
    dfx = 1.0 / (nx * pixel_size_m)
    dfy = 1.0 / (ny * pixel_size_m)
    fx = (np.arange(nx) - nx / 2.0) * dfx
    fy = (np.arange(ny) - ny / 2.0) * dfy
    fxx, fyy = np.meshgrid(fx, fy)

    # Transfer function (evanescent waves are filtered)
    root_arg = 1.0 - (wavelength_m * fxx) ** 2 - (wavelength_m * fyy) ** 2
    transfer = np.zeros_like(root_arg, dtype=np.complex128)
    valid = root_arg >= 0.0
    transfer[valid] = np.exp(1j * k * z_m * np.sqrt(root_arg[valid]))

    field_f = fftshift(fft2(ifftshift(field)))
    propagated_f = field_f * transfer
    propagated = fftshift(ifft2(ifftshift(propagated_f)))

    return propagated


# ===================================================================
# Hologram synthesis
# ===================================================================

def _to_uint8(image: np.ndarray) -> np.ndarray:
    """Quantize a float image to uint8 with percentile clipping (same as original)."""
    clipped = image.copy()
    low = np.percentile(clipped, 0.1)
    high = np.percentile(clipped, 99.9)
    clipped = np.clip(clipped, low, high)

    imin = clipped.min()
    imax = clipped.max()
    if imax == imin:
        return np.zeros(clipped.shape, dtype=np.uint8)

    normalized = (clipped - imin) / (imax - imin)
    return np.uint8(np.round(normalized * 255.0))


def synthesize_hologram(
    particles: list[Particle],
    N: int,
    wavelength_m: float,
    pixel_size_m: float,
    dn: float,
    attenuation: float,
    edge_sigma_px: float,
    roughness: float,
    snr_db: Optional[float],
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Synthesize an inline hologram from a list of particles.

    Returns:
        hologram_float: float64 hologram (N×N), values roughly [0, ~4]
        object_projection: float64 amplitude projection (N×N), [0, 1]
        hologram_noisy: float64 hologram with noise added (if SNR provided)
    """
    pixel_size_um = pixel_size_m * 1e6

    # Accumulate scattered field
    u_total = np.zeros((N, N), dtype=np.complex128)

    # Also build object projection for ground truth visualization
    object_proj = np.zeros((N, N), dtype=np.float64)

    for p in particles:
        # Complex transmittance
        t = compute_complex_transmittance(
            shape=p.shape,
            x_pixel=p.x_pixel,
            y_pixel=p.y_pixel,
            N=N,
            pixel_size_um=pixel_size_um,
            wavelength_m=wavelength_m,
            dn=dn,
            attenuation=attenuation,
            edge_sigma_px=edge_sigma_px,
            roughness=roughness,
            rng=rng,
        )

        # Scattered field component
        scatter = t - 1.0

        # Propagate to sensor
        u_scatter = angular_spectrum_propagation(scatter, wavelength_m, pixel_size_m, p.z_position_m)
        u_total += u_scatter

        # Object projection: amplitude of the particle mask
        mask = np.abs(scatter) > 1e-9
        object_proj = np.maximum(object_proj, mask.astype(np.float64))

    # Reference wave
    reference = 1.0 + 0.0j

    # Inline hologram
    hologram_float = np.abs(reference - u_total) ** 2

    # Add noise if requested
    hologram_noisy = hologram_float.copy()
    if snr_db is not None and snr_db > 0:
        signal_var = np.var(hologram_float)
        noise_var = signal_var / (10 ** (snr_db / 10.0)) if signal_var > 0 else 0.0
        noise_std = np.sqrt(noise_var)
        hologram_noisy = hologram_float + rng.normal(0, noise_std, hologram_float.shape)

    return hologram_float, object_proj, hologram_noisy


# ===================================================================
# Output
# ===================================================================

def save_outputs(
    output_dir: str,
    prefix: str,
    hologram_float: np.ndarray,
    hologram_noisy: np.ndarray,
    object_proj: np.ndarray,
    particles: list[Particle],
    add_noise: bool,
    metadata: dict,
) -> None:
    """Save all output files: BMP, .mat, and CSV."""
    os.makedirs(output_dir, exist_ok=True)

    holo_bmp_path = os.path.join(output_dir, f"{prefix}hologram.bmp")
    holo_mat_path = os.path.join(output_dir, f"{prefix}hologram.mat")
    obj_bmp_path = os.path.join(output_dir, f"{prefix}object_reference.bmp")
    csv_path = os.path.join(output_dir, f"{prefix}ground_truth.csv")
    meta_path = os.path.join(output_dir, f"{prefix}metadata.json")

    # 8-bit BMP (from noisy version if noise is on, otherwise clean)
    holo_uint8 = _to_uint8(hologram_noisy if add_noise else hologram_float)
    Image.fromarray(holo_uint8, mode="L").save(holo_bmp_path)
    print(f"  [OK] 全息图 BMP: {holo_bmp_path}")

    # .mat file with exact float values
    savemat(
        holo_mat_path,
        {
            "hologram": hologram_float.astype(np.float64),
            "hologram_noisy": hologram_noisy.astype(np.float64)
            if add_noise
            else hologram_float.astype(np.float64),
            "object_projection": object_proj.astype(np.float64),
            "metadata": metadata,
        },
    )
    print(f"  [OK] 全息图 MAT: {holo_mat_path}")

    # Object projection BMP
    obj_uint8 = _to_uint8(object_proj)
    Image.fromarray(obj_uint8, mode="L").save(obj_bmp_path)
    print(f"  [OK] 物面投影 BMP: {obj_bmp_path}")

    # Ground truth CSV
    rows = []
    for p in particles:
        rows.append(
            {
                "particle_id": p.particle_id,
                "shape_type": p.shape_type,
                "x_pixel": round(p.x_pixel, 3),
                "y_pixel": round(p.y_pixel, 3),
                "x_pixel_matlab": round(p.x_pixel + 1, 3),
                "y_pixel_matlab": round(p.y_pixel + 1, 3),
                "z_position_m": p.z_position_m,
                "z_position_mm": p.z_position_m * 1e3,
                "equivalent_diameter_um": round(p.equivalent_diameter_um, 3),
                "volume_um3": round(p.volume_um3, 3),
                **p.shape.shape_params,
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(csv_path, index=False)
    print(f"  [OK] 真值 CSV: {csv_path}  ({len(rows)} 粒子)")

    # Metadata JSON
    import json

    meta_serializable = {}
    for k, v in metadata.items():
        if isinstance(v, np.ndarray):
            meta_serializable[k] = v.tolist()
        elif isinstance(v, (np.integer,)):
            meta_serializable[k] = int(v)
        elif isinstance(v, (np.floating,)):
            meta_serializable[k] = float(v)
        else:
            meta_serializable[k] = v
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta_serializable, f, indent=2, ensure_ascii=False)
    print(f"  [OK] 元数据 JSON: {meta_path}")


# ===================================================================
# Main generation function
# ===================================================================

def generate_hologram(
    N: int = 1024,
    wavelength_m: float = 638e-9,
    pixel_size_m: float = 3.45e-6,
    n_circular: int = 10,
    n_irregular: int = 10,
    radius_um_min: float = 50.0,
    radius_um_max: float = 100.0,
    z_min_m: float = 20e-3,
    z_max_m: float = 35e-3,
    n_particle: float = 1.5,
    n_medium: float = 1.0,
    attenuation: float = 0.6,
    edge_sigma_px: float = 0.8,
    roughness: float = 0.03,
    snr_db: float = 0.0,
    output_dir: str = "./output/0001",
    prefix: str = "",
    seed: int = 42,
) -> str:
    """Generate one synthetic inline hologram for shape classification.

    Parameters
    ----------
    N : int
        Image size (N×N pixels).
    wavelength_m : float
        Laser wavelength in meters (default 638 nm).
    pixel_size_m : float
        Pixel pitch in meters (default 2.2 µm).
    n_circular : int
        Number of circular (sphere) particles.
    n_irregular : int
        Number of irregular particles (split between aggregate & polyhedron).
    radius_um_min, radius_um_max : float
        Particle radius range in µm.
    z_min_m, z_max_m : float
        Depth range in meters (distance from sensor).
    n_particle : float
        Particle refractive index.
    n_medium : float
        Surrounding medium refractive index.
    attenuation : float
        Amplitude attenuation inside particle (0=opaque, 1=transparent).
    edge_sigma_px : float
        Gaussian blur sigma at particle edges (pixels).
    roughness : float
        Surface roughness std for irregular particles.
    snr_db : float
        SNR in dB for additive Gaussian noise. 0 or negative = no noise.
    output_dir : str
        Output directory path.
    prefix : str
        Filename prefix (e.g. "0001_").
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    output_dir : str
        The output directory path.
    """
    rng = np.random.default_rng(seed)

    dn = n_particle - n_medium
    pixel_size_um = pixel_size_m * 1e6
    radius_um_range = (radius_um_min, radius_um_max)
    z_range = (z_min_m, z_max_m)

    add_noise = snr_db > 0

    print("=" * 60)
    print("  数字全息仿真 -- 粒子形状分类数据生成器")
    print("=" * 60)
    print(f"  图像尺寸:  {N}x{N} px")
    print(f"  波长:      {wavelength_m * 1e9:.1f} nm")
    print(f"  像元尺寸:  {pixel_size_um:.1f} um")
    print(f"  折射率:    n_particle={n_particle}, n_medium={n_medium} (Δn={dn:.2f})")
    print(f"  粒径范围:  {radius_um_min:.0f}-{radius_um_max:.0f} um (半径)")
    print(f"  深度范围:  {z_min_m * 1e3:.0f}-{z_max_m * 1e3:.0f} mm")
    print(f"  圆形粒子:  {n_circular} 个")
    print(f"  不规则粒子: {n_irregular} 个 (聚合球体 + 凸多面体)")
    print(f"  振幅衰减:  {attenuation:.2f}")
    print(f"  噪声 SNR:  {snr_db} dB" if add_noise else "  噪声:      关闭")
    print(f"  随机种子:  {seed}")
    print(f"  输出目录:  {output_dir}")
    print("-" * 60)

    # Place particles
    print("[1/4] 放置粒子...")
    particles = place_particles(
        N=N,
        pixel_size_m=pixel_size_m,
        n_circular=n_circular,
        n_irregular=n_irregular,
        radius_um_range=radius_um_range,
        z_range_m=z_range,
        rng=rng,
    )
    n_sphere = sum(1 for p in particles if p.shape_type == "sphere")
    n_agg = sum(1 for p in particles if p.shape_type == "aggregate")
    n_poly = sum(1 for p in particles if p.shape_type == "polyhedron")
    print(f"  已放置 {len(particles)} 个粒子: {n_sphere} 球形, {n_agg} 聚合体, {n_poly} 多面体")

    # Synthesize hologram
    print("[2/4] 合成全息图 (复透射率薄屏 + 角谱传播)...")
    hologram_float, object_proj, hologram_noisy = synthesize_hologram(
        particles=particles,
        N=N,
        wavelength_m=wavelength_m,
        pixel_size_m=pixel_size_m,
        dn=dn,
        attenuation=attenuation,
        edge_sigma_px=edge_sigma_px,
        roughness=roughness,
        snr_db=snr_db if add_noise else None,
        rng=rng,
    )
    print(f"  全息图强度范围: [{hologram_float.min():.4f}, {hologram_float.max():.4f}]")

    # Metadata
    print("[3/4] 收集元数据...")
    metadata = {
        "N": N,
        "wavelength_m": wavelength_m,
        "pixel_size_m": pixel_size_m,
        "pixel_size_um": pixel_size_um,
        "n_particle": n_particle,
        "n_medium": n_medium,
        "dn": dn,
        "radius_um_range": list(radius_um_range),
        "z_range_m": [z_min_m, z_max_m],
        "n_circular": n_circular,
        "n_irregular": n_irregular,
        "n_placed": len(particles),
        "n_sphere": n_sphere,
        "n_aggregate": n_agg,
        "n_polyhedron": n_poly,
        "attenuation": attenuation,
        "edge_sigma_px": edge_sigma_px,
        "roughness": roughness,
        "snr_db": snr_db,
        "seed": seed,
    }

    # Save
    print("[4/4] 保存输出文件...")
    save_outputs(
        output_dir=output_dir,
        prefix=prefix,
        hologram_float=hologram_float,
        hologram_noisy=hologram_noisy if add_noise else hologram_float,
        object_proj=object_proj,
        particles=particles,
        add_noise=add_noise,
        metadata=metadata,
    )

    print("=" * 60)
    print("  生成完成！")
    print("=" * 60)
    return output_dir


# ===================================================================
# CLI
# ===================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate synthetic holograms for particle shape classification.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single hologram
  python generate_hologram_shape_classifier.py -o ./output/0001 -s 42

  # Batch: generate 50 holograms
  for i in $(seq 1 50); do
      python generate_hologram_shape_classifier.py -o ./output/$(printf "%%04d" $i) -s $i
  done

  # With noise and custom particle counts
  python generate_hologram_shape_classifier.py -o ./test -N 1024 --circular 15 --irregular 5 --snr 25
        """,
    )
    p.add_argument("-o", "--output-dir", default="./output/0001", help="Output directory.")
    p.add_argument("-p", "--prefix", default="", help="Filename prefix.")
    p.add_argument("-s", "--seed", type=int, default=42, help="Random seed.")
    p.add_argument("-N", "--img-size", type=int, default=1024, help="Image size NxN.")
    p.add_argument("--wavelength", type=float, default=638e-9, help="Wavelength (m).")
    p.add_argument("--pixel-size", type=float, default=3.45e-6, help="Pixel pitch (m).")
    p.add_argument("--circular", type=int, default=10, help="Number of circular particles.")
    p.add_argument("--irregular", type=int, default=10, help="Number of irregular particles.")
    p.add_argument("--radius-min", type=float, default=50.0, help="Min particle radius (um).")
    p.add_argument("--radius-max", type=float, default=100.0, help="Max particle radius (um).")
    p.add_argument("--z-min", type=float, default=20.0, help="Min depth (mm).")
    p.add_argument("--z-max", type=float, default=35.0, help="Max depth (mm).")
    p.add_argument("--n-particle", type=float, default=1.5, help="Particle refractive index.")
    p.add_argument("--n-medium", type=float, default=1.0, help="Medium refractive index.")
    p.add_argument("--attenuation", type=float, default=0.6, help="Amplitude attenuation in particle.")
    p.add_argument("--edge-sigma", type=float, default=0.8, help="Edge blur sigma (pixels).")
    p.add_argument("--roughness", type=float, default=0.03, help="Surface roughness std (irregular only).")
    p.add_argument("--snr", type=float, default=0.0, help="Noise SNR (dB). 0=no noise.")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    generate_hologram(
        N=args.img_size,
        wavelength_m=args.wavelength,
        pixel_size_m=args.pixel_size,
        n_circular=args.circular,
        n_irregular=args.irregular,
        radius_um_min=args.radius_min,
        radius_um_max=args.radius_max,
        z_min_m=args.z_min * 1e-3,
        z_max_m=args.z_max * 1e-3,
        n_particle=args.n_particle,
        n_medium=args.n_medium,
        attenuation=args.attenuation,
        edge_sigma_px=args.edge_sigma,
        roughness=args.roughness,
        snr_db=args.snr,
        output_dir=args.output_dir,
        prefix=args.prefix,
        seed=args.seed,
    )
