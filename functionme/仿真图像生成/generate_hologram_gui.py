#!/usr/bin/env python3
"""
generate_hologram_gui.py
=========================
GUI for generating digital holograms.

Two modes:
  - 标准仿真: Original MATLAB-validation hologram generator
               (binary mask + Born approximation, small particles 10-20um)
  - 形状分类: Shape-classification hologram generator
               (complex transmittance thin screen, 50-100um particles,
                sphere / aggregate / polyhedron types, 1024x1024)

Usage:
  python generate_hologram_gui.py
"""

import os
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from generate_hologram_matlab_validation_new import generate_hologram as gen_standard
from generate_hologram_shape_classifier import generate_hologram as gen_shape


class HologramGeneratorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("数字全息仿真图像生成器")
        self.root.geometry("560x820")
        self.root.resizable(True, True)
        self.root.minsize(520, 700)

        style = ttk.Style()
        style.theme_use("clam")

        # ---- Top: mode selector ----
        top_bar = ttk.Frame(self.root, padding="10 10 10 5")
        top_bar.pack(fill=tk.X)

        ttk.Label(top_bar, text="生成模式:", font=("Arial", 10, "bold")).pack(side=tk.LEFT)
        self.mode_var = tk.StringVar(value="标准仿真")
        mode_combo = ttk.Combobox(
            top_bar,
            textvariable=self.mode_var,
            values=["标准仿真", "形状分类"],
            state="readonly",
            width=16,
        )
        mode_combo.pack(side=tk.LEFT, padx=10)
        mode_combo.bind("<<ComboboxSelected>>", self.on_mode_changed)

        # ---- Notebook for parameter pages ----
        self.notebook = ttk.Notebook(self.root, padding="10 5 10 10")
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Page 1: standard mode
        self.page_standard = ttk.Frame(self.notebook)
        self.notebook.add(self.page_standard, text="标准仿真")

        # Page 2: shape classification mode
        self.page_shape = ttk.Frame(self.notebook)
        self.notebook.add(self.page_shape, text="形状分类")

        # ---- Shared: output + action bar ----
        bottom_bar = ttk.Frame(self.root, padding="10 5 10 10")
        bottom_bar.pack(fill=tk.X, side=tk.BOTTOM)

        out_frame = ttk.Frame(bottom_bar)
        out_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Label(out_frame, text="保存至:").pack(side=tk.LEFT)
        self.outdir_var = tk.StringVar(value=os.path.abspath("./"))
        ttk.Entry(out_frame, textvariable=self.outdir_var, width=40).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(out_frame, text="浏览...", command=self.browse_dir).pack(side=tk.LEFT)

        self.generate_btn = ttk.Button(bottom_bar, text="生成全息图及真值", command=self.start_generation)
        self.generate_btn.pack(pady=8, ipadx=10, ipady=5)

        self.status_var = tk.StringVar(value="就绪")
        ttk.Label(bottom_bar, textvariable=self.status_var, foreground="gray").pack()

        # ---- Build parameter pages ----
        self._build_standard_page()
        self._build_shape_page()

        # Show default page
        self.on_mode_changed()

    # ==================================================================
    # Standard simulation page (unchanged from original)
    # ==================================================================

    def _build_standard_page(self):
        parent = self.page_standard
        canvas = tk.Canvas(parent, highlightthickness=0)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scroll_frame = ttk.Frame(canvas)

        scroll_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all")),
        )
        canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        canvas.bind_all("<MouseWheel>", _on_mousewheel)

        row = 0

        # System params
        ttk.Label(scroll_frame, text="系统参数", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.wavelength_var = tk.DoubleVar(value=632.8)
        self.dx_var = tk.DoubleVar(value=3.45)
        self.dy_var = tk.DoubleVar(value=3.45)
        self.n_var = tk.IntVar(value=1024)

        self._add_entry(scroll_frame, "波长 (nm):", self.wavelength_var, row)
        row += 1
        self._add_entry(scroll_frame, "像素尺寸 dx (um):", self.dx_var, row)
        row += 1
        self._add_entry(scroll_frame, "像素尺寸 dy (um):", self.dy_var, row)
        row += 1
        self._add_entry(scroll_frame, "图像尺寸 N:", self.n_var, row)
        row += 1

        self._add_sep(scroll_frame, row)
        row += 1

        # Particle params
        ttk.Label(scroll_frame, text="粒子参数（圆形）", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.dmin_var = tk.DoubleVar(value=10.0)
        self.dmax_var = tk.DoubleVar(value=20.0)
        self.num_var = tk.IntVar(value=15)

        self._add_entry(scroll_frame, "最小直径 (um):", self.dmin_var, row)
        row += 1
        self._add_entry(scroll_frame, "最大直径 (um):", self.dmax_var, row)
        row += 1
        self._add_entry(scroll_frame, "粒子数量:", self.num_var, row)
        row += 1

        self._add_sep(scroll_frame, row)
        row += 1

        # Z-axis
        ttk.Label(scroll_frame, text="轴向分布 (Z轴)", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.zmin_var = tk.DoubleVar(value=20.0)
        self.zmax_var = tk.DoubleVar(value=35.0)
        self.zstep_var = tk.DoubleVar(value=40.0)

        self._add_entry(scroll_frame, "最小距离 Z_min (mm):", self.zmin_var, row)
        row += 1
        self._add_entry(scroll_frame, "最大距离 Z_max (mm):", self.zmax_var, row)
        row += 1
        self._add_entry(scroll_frame, "轴向层间距 (um):", self.zstep_var, row)
        row += 1

        self._add_sep(scroll_frame, row)
        row += 1

        # Noise
        ttk.Label(scroll_frame, text="噪声设置", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1
        ttk.Label(scroll_frame, text="(留空或输入 0 表示无噪声)").grid(row=row, column=1, sticky=tk.W, padx=5)
        row += 1

        self.snr_var = tk.DoubleVar(value=20.0)
        self._add_entry(scroll_frame, "信噪比 SNR (dB):", self.snr_var, row)
        row += 1

        # Spacer
        ttk.Label(scroll_frame, text="").grid(row=row, column=0)
        row += 1
        ttk.Label(scroll_frame, text="提示: 标准仿真使用二值掩膜 + Born 近似", foreground="gray").grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(10, 0)
        )
        row += 1
        ttk.Label(scroll_frame, text="适用于 10-20um 级小粒子，输出 BMP + CSV", foreground="gray").grid(
            row=row, column=0, columnspan=3, sticky=tk.W
        )

    # ==================================================================
    # Shape classification page
    # ==================================================================

    def _build_shape_page(self):
        parent = self.page_shape
        canvas = tk.Canvas(parent, highlightthickness=0)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scroll_frame = ttk.Frame(canvas)

        scroll_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all")),
        )
        canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        canvas.bind_all("<MouseWheel>", _on_mousewheel)

        row = 0

        # ---- System params ----
        ttk.Label(scroll_frame, text="系统参数", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.s_wavelength_var = tk.DoubleVar(value=638.0)
        self.s_pixel_var = tk.DoubleVar(value=3.45)
        self.s_N_var = tk.IntVar(value=1024)

        self._add_entry(scroll_frame, "波长 (nm):", self.s_wavelength_var, row)
        row += 1
        self._add_entry(scroll_frame, "像素尺寸 (um):", self.s_pixel_var, row)
        row += 1
        self._add_entry(scroll_frame, "图像尺寸 N:", self.s_N_var, row)
        row += 1

        self._add_sep(scroll_frame, row)
        row += 1

        # ---- Particle count ----
        ttk.Label(scroll_frame, text="粒子数量与类型", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.s_circular_var = tk.IntVar(value=10)
        self.s_irregular_var = tk.IntVar(value=10)

        self._add_entry(scroll_frame, "圆形粒子数:", self.s_circular_var, row)
        ttk.Label(scroll_frame, text="(sphere, 标准球体)", foreground="gray").grid(
            row=row, column=2, sticky=tk.W, padx=5
        )
        row += 1
        self._add_entry(scroll_frame, "不规则粒子数:", self.s_irregular_var, row)
        ttk.Label(scroll_frame, text="(聚合球体 + 凸多面体)", foreground="gray").grid(
            row=row, column=2, sticky=tk.W, padx=5
        )
        row += 1

        self._add_sep(scroll_frame, row)
        row += 1

        # ---- Particle size ----
        ttk.Label(scroll_frame, text="粒子尺寸（半径）", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.s_radius_min_var = tk.DoubleVar(value=50.0)
        self.s_radius_max_var = tk.DoubleVar(value=100.0)

        self._add_entry(scroll_frame, "最小半径 (um):", self.s_radius_min_var, row)
        row += 1
        self._add_entry(scroll_frame, "最大半径 (um):", self.s_radius_max_var, row)
        row += 1

        self._add_sep(scroll_frame, row)
        row += 1

        # ---- Optical properties ----
        ttk.Label(scroll_frame, text="光学参数", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.s_n_particle_var = tk.DoubleVar(value=1.5)
        self.s_n_medium_var = tk.DoubleVar(value=1.0)
        self.s_attenuation_var = tk.DoubleVar(value=0.6)

        self._add_entry(scroll_frame, "粒子折射率 n:", self.s_n_particle_var, row)
        row += 1
        self._add_entry(scroll_frame, "介质折射率 n0:", self.s_n_medium_var, row)
        row += 1
        self._add_entry(scroll_frame, "振幅衰减 (0~1):", self.s_attenuation_var, row)
        ttk.Label(scroll_frame, text="(0=不透明, 1=全透明)", foreground="gray").grid(
            row=row, column=2, sticky=tk.W, padx=5
        )
        row += 1

        self._add_sep(scroll_frame, row)
        row += 1

        # ---- Z range ----
        ttk.Label(scroll_frame, text="轴向范围 (Z轴)", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.s_zmin_var = tk.DoubleVar(value=20.0)
        self.s_zmax_var = tk.DoubleVar(value=35.0)

        self._add_entry(scroll_frame, "最小距离 (mm):", self.s_zmin_var, row)
        row += 1
        self._add_entry(scroll_frame, "最大距离 (mm):", self.s_zmax_var, row)
        row += 1

        self._add_sep(scroll_frame, row)
        row += 1

        # ---- Edge & roughness ----
        ttk.Label(scroll_frame, text="表面与边缘参数", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.s_edge_sigma_var = tk.DoubleVar(value=0.8)
        self.s_roughness_var = tk.DoubleVar(value=0.03)

        self._add_entry(scroll_frame, "边缘模糊 sigma (px):", self.s_edge_sigma_var, row)
        ttk.Label(scroll_frame, text="(粒子边缘高斯过渡宽度)", foreground="gray").grid(
            row=row, column=2, sticky=tk.W, padx=5
        )
        row += 1
        self._add_entry(scroll_frame, "表面粗糙度:", self.s_roughness_var, row)
        ttk.Label(scroll_frame, text="(不规则粒子 A(x,y) 扰动 std)", foreground="gray").grid(
            row=row, column=2, sticky=tk.W, padx=5
        )
        row += 1

        self._add_sep(scroll_frame, row)
        row += 1

        # ---- Polyhedron fragmentation ----
        ttk.Label(scroll_frame, text="多面体破碎度", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.s_frag_var = tk.StringVar(value="custom")
        frag_frame = ttk.Frame(scroll_frame)
        frag_frame.grid(row=row, column=0, columnspan=3, sticky=tk.W, pady=2)
        ttk.Label(frag_frame, text="预设:").pack(side=tk.LEFT)
        frag_combo = ttk.Combobox(frag_frame, textvariable=self.s_frag_var,
                                   values=["custom", "mild", "medium", "severe"],
                                   state="readonly", width=10)
        frag_combo.pack(side=tk.LEFT, padx=5)
        ttk.Label(frag_frame, text="(custom=使用下方参数 / mild / medium / severe)", foreground="gray").pack(side=tk.LEFT)
        row += 1

        ttk.Label(scroll_frame, text="(选择 mild/medium/severe 后忽略下方精细参数)", foreground="gray").grid(
            row=row, column=0, columnspan=3, sticky=tk.W, padx=5
        )
        row += 1

        self.s_poly_jitter_var = tk.DoubleVar(value=0.40)
        self.s_poly_spike_var = tk.DoubleVar(value=0.12)
        self.s_poly_indent_var = tk.DoubleVar(value=0.20)
        self.s_poly_vmin_var = tk.IntVar(value=15)
        self.s_poly_vmax_var = tk.IntVar(value=35)

        self._add_entry(scroll_frame, "径向抖动 (jitter):", self.s_poly_jitter_var, row)
        row += 1
        self._add_entry(scroll_frame, "突刺比例 (spike):", self.s_poly_spike_var, row)
        row += 1
        self._add_entry(scroll_frame, "内缩比例 (indent):", self.s_poly_indent_var, row)
        row += 1
        self._add_entry(scroll_frame, "顶点数范围:", self.s_poly_vmin_var, row)
        ttk.Label(scroll_frame, text=f"(最小, 最大={self.s_poly_vmax_var.get()})", foreground="gray").grid(
            row=row, column=2, sticky=tk.W, padx=5
        )
        row += 1
        self._add_entry(scroll_frame, "顶点数最大:", self.s_poly_vmax_var, row)
        row += 1

        self.s_num_frags_var = tk.IntVar(value=3)
        self.s_frag_overlap_var = tk.DoubleVar(value=0.55)

        self._add_entry(scroll_frame, "碎片数量:", self.s_num_frags_var, row)
        ttk.Label(scroll_frame, text="(2-8, 多碎片=非凸凹口多)", foreground="gray").grid(
            row=row, column=2, sticky=tk.W, padx=5
        )
        row += 1
        self._add_entry(scroll_frame, "碎片重叠度:", self.s_frag_overlap_var, row)
        ttk.Label(scroll_frame, text="(~0.3=松散, ~0.7=紧密)", foreground="gray").grid(
            row=row, column=2, sticky=tk.W, padx=5
        )
        row += 1

        self._add_sep(scroll_frame, row)
        row += 1

        # ---- Noise ----
        ttk.Label(scroll_frame, text="噪声设置", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.s_snr_var = tk.DoubleVar(value=0.0)
        self._add_entry(scroll_frame, "信噪比 SNR (dB):", self.s_snr_var, row)
        ttk.Label(scroll_frame, text="(0 = 无噪声)", foreground="gray").grid(
            row=row, column=2, sticky=tk.W, padx=5
        )
        row += 1

        self._add_sep(scroll_frame, row)
        row += 1

        # ---- Seed ----
        ttk.Label(scroll_frame, text="随机种子", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.s_seed_var = tk.IntVar(value=42)
        self._add_entry(scroll_frame, "随机种子:", self.s_seed_var, row)
        ttk.Label(scroll_frame, text="(同种子 = 可复现)", foreground="gray").grid(
            row=row, column=2, sticky=tk.W, padx=5
        )
        row += 1

        # Spacer + hints
        ttk.Label(scroll_frame, text="").grid(row=row, column=0)
        row += 1
        hints = [
            "物理模型: 复透射率薄屏 t=A*exp(ik*Dn*h) + 角谱传播",
            "输出: BMP + .mat (float64) + CSV (含 particle_id/type/volume)",
            "不规则粒子: 一半聚合球体 (PTFE-like), 一半凸多面体 (ice-crystal-like)",
            "圆形/不规则走同一物理管线, 差异仅来自形态本身",
        ]
        for hint in hints:
            ttk.Label(scroll_frame, text=f"  > {hint}", foreground="gray").grid(
                row=row, column=0, columnspan=3, sticky=tk.W
            )
            row += 1

    # ==================================================================
    # Helpers
    # ==================================================================

    @staticmethod
    def _add_entry(parent, label_text, variable, row):
        ttk.Label(parent, text=label_text).grid(row=row, column=0, sticky=tk.W, pady=2)
        ttk.Entry(parent, textvariable=variable, width=16).grid(row=row, column=1, sticky=tk.W, padx=5, pady=2)

    @staticmethod
    def _add_sep(parent, row):
        ttk.Separator(parent, orient="horizontal").grid(
            row=row, column=0, columnspan=3, sticky="ew", pady=8
        )

    def on_mode_changed(self, event=None):
        """Switch the visible notebook tab based on mode selection."""
        mode = self.mode_var.get()
        if mode == "标准仿真":
            self.notebook.select(self.page_standard)
        else:
            self.notebook.select(self.page_shape)

    def browse_dir(self):
        dir_name = filedialog.askdirectory(initialdir=self.outdir_var.get(), title="选择输出目录")
        if dir_name:
            self.outdir_var.set(dir_name)

    def get_next_prefix(self, out_dir):
        import re

        max_num = 0
        if os.path.exists(out_dir):
            file_pattern = re.compile(r"^(\d+)_hologram\.bmp$")
            dir_pattern = re.compile(r"^(\d+)$")
            for entry in os.listdir(out_dir):
                dir_match = dir_pattern.match(entry)
                if dir_match:
                    max_num = max(max_num, int(dir_match.group(1)))
                    continue
                file_match = file_pattern.match(entry)
                if file_match:
                    max_num = max(max_num, int(file_match.group(1)))
        return f"{max_num + 1:04d}_"

    def start_generation(self):
        self.generate_btn.config(state=tk.DISABLED)
        self.status_var.set("正在生成，请稍候...")
        self.root.update()

        worker = threading.Thread(target=self.run_generation)
        worker.daemon = True
        worker.start()

    def run_generation(self):
        try:
            mode = self.mode_var.get()
            out_dir = self.outdir_var.get()
            prefix = self.get_next_prefix(out_dir)
            run_dir = os.path.join(out_dir, prefix.rstrip("_"))

            if mode == "标准仿真":
                saved_dir = self._run_standard(run_dir, prefix)
            else:
                saved_dir = self._run_shape(run_dir, prefix)

            self.root.after(0, self.generation_complete, True, (prefix, saved_dir))
        except Exception as exc:
            import traceback

            traceback.print_exc()
            self.root.after(0, self.generation_complete, False, str(exc))

    def _run_standard(self, run_dir, prefix):
        wl = self.wavelength_var.get() * 1e-9
        dx = self.dx_var.get() * 1e-6
        dy = self.dy_var.get() * 1e-6
        N = self.n_var.get()

        dmin = self.dmin_var.get() * 1e-6
        dmax = self.dmax_var.get() * 1e-6
        num = self.num_var.get()

        zmin = self.zmin_var.get() * 1e-3
        zmax = self.zmax_var.get() * 1e-3
        zstep = self.zstep_var.get() * 1e-6

        snr = self.snr_var.get()

        return gen_standard(
            wavelength=wl,
            dx=dx,
            dy=dy,
            N=N,
            diameter_min=dmin,
            diameter_max=dmax,
            num_particles=num,
            z_min=zmin,
            z_max=zmax,
            z_step=zstep,
            snr_db=snr,
            output_dir=run_dir,
            prefix=prefix,
        )

    def _run_shape(self, run_dir, prefix):
        return gen_shape(
            N=self.s_N_var.get(),
            wavelength_m=self.s_wavelength_var.get() * 1e-9,
            pixel_size_m=self.s_pixel_var.get() * 1e-6,
            n_circular=self.s_circular_var.get(),
            n_irregular=self.s_irregular_var.get(),
            radius_um_min=self.s_radius_min_var.get(),
            radius_um_max=self.s_radius_max_var.get(),
            z_min_m=self.s_zmin_var.get() * 1e-3,
            z_max_m=self.s_zmax_var.get() * 1e-3,
            n_particle=self.s_n_particle_var.get(),
            n_medium=self.s_n_medium_var.get(),
            attenuation=self.s_attenuation_var.get(),
            edge_sigma_px=self.s_edge_sigma_var.get(),
            roughness=self.s_roughness_var.get(),
            snr_db=self.s_snr_var.get(),
            fragmentation=self.s_frag_var.get(),
            poly_radial_jitter=self.s_poly_jitter_var.get(),
            poly_spike_fraction=self.s_poly_spike_var.get(),
            poly_indent_fraction=self.s_poly_indent_var.get(),
            poly_vertex_min=self.s_poly_vmin_var.get(),
            poly_vertex_max=self.s_poly_vmax_var.get(),
            poly_num_fragments=self.s_num_frags_var.get(),
            poly_fragment_overlap=self.s_frag_overlap_var.get(),
            output_dir=run_dir,
            prefix=prefix,
            seed=self.s_seed_var.get(),
        )

    def generation_complete(self, success, msg):
        self.generate_btn.config(state=tk.NORMAL)
        if success:
            prefix, saved_dir = msg
            self.status_var.set("生成完成")
            messagebox.showinfo(
                "成功",
                f"全息图及真值已成功生成并保存至:\n{saved_dir}\n\n(自动分配前缀编号: {prefix})",
            )
        else:
            self.status_var.set("生成失败")
            messagebox.showerror("错误", f"生成过程中出现错误:\n{msg}")


if __name__ == "__main__":
    root = tk.Tk()
    app = HologramGeneratorGUI(root)

    root.update_idletasks()
    width = root.winfo_width()
    height = root.winfo_height()
    x = (root.winfo_screenwidth() // 2) - (width // 2)
    y = (root.winfo_screenheight() // 2) - (height // 2)
    root.geometry(f"+{x}+{y}")

    root.mainloop()
