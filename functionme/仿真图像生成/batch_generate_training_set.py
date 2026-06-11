#!/usr/bin/env python3
"""
batch_generate_training_set.py
===============================
Batch generate holograms for particle shape classification training.

Generates N holograms, each in its own numbered folder under output_dir.
Each hologram contains a mix of circular (sphere) and irregular (aggregate +
polyhedron) particles with known ground truth.

Usage:
  python batch_generate_training_set.py --num 50 --output_dir ./training_set
  python batch_generate_training_set.py --num 10 --output_dir ./val_set --circular 8 --irregular 12
"""

from __future__ import annotations

import argparse
import os
import time
import json

from generate_hologram_shape_classifier import generate_hologram


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Batch generate holograms for shape classification training.",
    )
    p.add_argument("--num", type=int, default=50, help="Number of holograms to generate.")
    p.add_argument("--output-dir", "-o", default="./training_set", help="Root output directory.")
    p.add_argument("--N", type=int, default=1024, help="Image size.")
    p.add_argument("--circular", type=int, default=10, help="Circular particles per hologram.")
    p.add_argument("--irregular", type=int, default=10, help="Irregular particles per hologram.")
    p.add_argument("--radius-min", type=float, default=50.0, help="Min particle radius (um).")
    p.add_argument("--radius-max", type=float, default=100.0, help="Max particle radius (um).")
    p.add_argument("--z-min", type=float, default=20.0, help="Min depth (mm).")
    p.add_argument("--z-max", type=float, default=35.0, help="Max depth (mm).")
    p.add_argument("--snr", type=float, default=0.0, help="Noise SNR (dB), 0=no noise.")
    p.add_argument("--seed-start", type=int, default=1, help="Starting seed (incremented per hologram).")
    p.add_argument("--wavelength", type=float, default=638e-9, help="Wavelength (m).")
    p.add_argument("--pixel-size", type=float, default=3.45e-6, help="Pixel size (m).")
    p.add_argument("--n-particle", type=float, default=1.38, help="Particle refractive index.")
    p.add_argument("--attenuation", type=float, default=0.30, help="Amplitude attenuation.")
    p.add_argument("--roughness", type=float, default=0.0, help="Surface roughness std (0=off).")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Save batch config
    config = vars(args).copy()
    config_path = os.path.join(args.output_dir, "batch_config.json")
    with open(config_path, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2, ensure_ascii=False)

    print("=" * 60)
    print(f"  批量生成训练集: {args.num} 张全息图")
    print(f"  输出目录: {args.output_dir}")
    print(f"  每张: {args.circular} 圆形 + {args.irregular} 不规则")
    print(f"  总粒子数: {args.num * (args.circular + args.irregular)}")
    print("=" * 60)

    total_start = time.time()
    success = 0

    for i in range(args.num):
        idx = i + 1
        seed = args.seed_start + i
        folder = os.path.join(args.output_dir, f"{idx:04d}")
        prefix = f"{idx:04d}_"

        print(f"\n[{idx}/{args.num}] seed={seed} ...")

        try:
            t0 = time.time()
            generate_hologram(
                N=args.N,
                wavelength_m=args.wavelength,
                pixel_size_m=args.pixel_size,
                n_circular=args.circular,
                n_irregular=args.irregular,
                radius_um_min=args.radius_min,
                radius_um_max=args.radius_max,
                z_min_m=args.z_min * 1e-3,
                z_max_m=args.z_max * 1e-3,
                n_particle=args.n_particle,
                attenuation=args.attenuation,
                roughness=args.roughness,
                snr_db=args.snr,
                output_dir=folder,
                prefix=prefix,
                seed=seed,
            )
            elapsed = time.time() - t0
            print(f"  -> 完成, 耗时 {elapsed:.1f}s")
            success += 1
        except Exception as e:
            print(f"  -> 失败: {e}")
            import traceback
            traceback.print_exc()

    total_elapsed = time.time() - total_start
    print("\n" + "=" * 60)
    print(f"  批量生成完成: {success}/{args.num} 成功")
    print(f"  总耗时: {total_elapsed:.1f}s ({total_elapsed/60:.1f} min)")
    print("=" * 60)
