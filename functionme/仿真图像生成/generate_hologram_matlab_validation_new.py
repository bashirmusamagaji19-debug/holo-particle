import os

import numpy as np
import pandas as pd
from PIL import Image
from scipy.fft import fft2, fftshift, ifft2, ifftshift


def _to_uint8(image):
    clipped = image.copy()
    low = np.percentile(clipped, 0.1)
    high = np.percentile(clipped, 99.9)
    clipped[clipped < low] = low
    clipped[clipped > high] = high

    image_min = clipped.min()
    image_max = clipped.max()
    if image_max == image_min:
        return np.zeros(clipped.shape, dtype=np.uint8)

    normalized = (clipped - image_min) / (image_max - image_min)
    return np.uint8(normalized * 255)


def generate_hologram(
    wavelength=632.8e-9,
    dx=2.2e-6,
    dy=2.2e-6,
    N=256,
    diameter_min=10e-6,
    diameter_max=20e-6,
    num_particles=15,
    z_min=2000e-6,
    z_max=8000e-6,
    z_step=40e-6,
    snr_db=20.0,
    output_dir="./",
    prefix="",
):
    print("开始生成用于 Matlab 验证的数字全息仿真图像...")

    os.makedirs(output_dir, exist_ok=True)

    holo_filename = os.path.normpath(os.path.join(output_dir, f"{prefix}hologram.bmp"))
    noisy_holo_filename = os.path.normpath(os.path.join(output_dir, f"{prefix}hologram_with_noise.bmp"))
    bg_noise_filename = os.path.normpath(os.path.join(output_dir, f"{prefix}background_noise.bmp"))
    obj_filename = os.path.normpath(os.path.join(output_dir, f"{prefix}object_reference.bmp"))
    csv_filename = os.path.normpath(os.path.join(output_dir, f"{prefix}particle_ground_truth.csv"))

    print(f"[*] 参数: 波长={wavelength * 1e9:.1f}nm, 像素尺寸={dx * 1e6:.1f}um, 图片尺寸={N}x{N}")
    print(f"[*] z分布: {z_min * 1e3:.1f}mm - {z_max * 1e3:.1f}mm, 粒子数: {num_particles}")

    add_noise = False
    try:
        if snr_db is not None and float(snr_db) > 0:
            add_noise = True
    except (TypeError, ValueError):
        pass

    if add_noise:
        print(f"[*] 噪声: 开启, SNR={snr_db}dB")
    else:
        print("[*] 噪声: 关闭 (无噪声输入或SNR<=0)")

    u_total = np.zeros((N, N), dtype=complex)
    obj_visual = np.zeros((N, N), dtype=float)
    particle_data = []

    k = 2 * np.pi / wavelength

    dfx = 1.0 / (N * dx)
    dfy = 1.0 / (N * dy)
    fx = (np.arange(N) - N / 2) * dfx
    fy = (np.arange(N) - N / 2) * dfy
    fx_grid, fy_grid = np.meshgrid(fx, fy)
    sqrt_arg = 1 - (wavelength * fx_grid) ** 2 - (wavelength * fy_grid) ** 2
    valid_mask = sqrt_arg >= 0

    for i in range(num_particles):
        r = np.random.uniform(diameter_min / 2, diameter_max / 2)
        r_pixel = r / dx

        margin = int(np.ceil(r_pixel)) + 5
        xo = np.random.uniform(margin, N - margin)
        yo = np.random.uniform(margin, N - margin)

        z_random = np.random.uniform(z_min, z_max)
        zi = np.round(z_random / z_step) * z_step

        particle_data.append(
            {
                "Particle_ID": i + 1,
                "x_pixel": xo,
                "y_pixel": yo,
                "radius_m": r,
                "diameter_um": r * 2 * 1e6,
                "z_position_m": zi,
            }
        )

        obj_plane = np.zeros((N, N), dtype=float)
        y_grid, x_grid = np.ogrid[0:N, 0:N]
        dist_sq = (x_grid - xo) ** 2 + (y_grid - yo) ** 2
        obj_plane[dist_sq <= r_pixel**2] = 1.0

        obj_visual = np.maximum(obj_visual, obj_plane)

        obj = fftshift(fft2(ifftshift(obj_plane)))
        transfer = np.zeros((N, N), dtype=complex)
        transfer[valid_mask] = np.exp(1j * k * zi * np.sqrt(sqrt_arg[valid_mask]))

        propagated = obj * transfer
        u_plane = fftshift(ifft2(ifftshift(propagated)))
        u_total += u_plane

    pd.DataFrame(particle_data).to_csv(csv_filename, index=False)
    print(f"[*] 成功保存真值坐标至:   {csv_filename}")

    reference = 1.0
    clean_hologram = np.abs(reference - u_total) ** 2
    noisy_hologram = clean_hologram.copy()
    background_noise = None

    if add_noise:
        signal_var = np.var(clean_hologram)
        noise_var = signal_var / (10 ** (snr_db / 10.0)) if signal_var > 0 else 0.0
        noise_std = np.sqrt(noise_var)

        noisy_hologram = noisy_hologram + np.random.normal(0, noise_std, clean_hologram.shape)
        background_noise = np.full(clean_hologram.shape, np.abs(reference) ** 2, dtype=float)
        background_noise = background_noise + np.random.normal(0, noise_std, clean_hologram.shape)

    clean_hologram_uint8 = _to_uint8(clean_hologram)
    obj_uint8 = np.uint8(obj_visual * 255)

    Image.fromarray(clean_hologram_uint8).save(holo_filename)
    Image.fromarray(obj_uint8).save(obj_filename)

    if add_noise:
        Image.fromarray(_to_uint8(noisy_hologram)).save(noisy_holo_filename)
        Image.fromarray(_to_uint8(background_noise)).save(bg_noise_filename)

    print(f"[*] 成功保存8位数字全息图: {holo_filename}")
    print(f"[*] 成功保存二值物体截面图: {obj_filename}")
    if add_noise:
        print(f"[*] 成功保存带噪声全息图: {noisy_holo_filename}")
        print(f"[*] 成功保存背景噪声图: {bg_noise_filename}")

    print("\n================== 仿真完成 ==================")
    print("您现在可以在 Matlab 中进行验证测试。建议方法:")
    print(" 1. I = imread('hologram.bmp');")
    print(" 2. I = double(I) / 255.0;")
    print(" 3. I = I - mean(I(:));  % 在 Matlab 内减去直流")
    print(" 4. 真值 CSV 中的 (x, y) 坐标需要 +1 才是 Matlab 中的真实像素坐标。")

    return output_dir


if __name__ == "__main__":
    generate_hologram()
