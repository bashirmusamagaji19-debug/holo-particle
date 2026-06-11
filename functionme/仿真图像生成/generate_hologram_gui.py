import os
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from generate_hologram_matlab_validation_new import generate_hologram


class HologramGeneratorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("数字全息仿真图像生成器")
        self.root.geometry("500x650")
        self.root.resizable(False, False)

        style = ttk.Style()
        style.theme_use("clam")

        main_frame = ttk.Frame(self.root, padding="20 20 20 20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        self.wavelength_var = tk.DoubleVar(value=632.8)
        self.dx_var = tk.DoubleVar(value=2.2)
        self.dy_var = tk.DoubleVar(value=2.2)
        self.n_var = tk.IntVar(value=256)

        self.dmin_var = tk.DoubleVar(value=10.0)
        self.dmax_var = tk.DoubleVar(value=20.0)
        self.num_var = tk.IntVar(value=15)

        self.zmin_var = tk.DoubleVar(value=2.0)
        self.zmax_var = tk.DoubleVar(value=8.0)
        self.zstep_var = tk.DoubleVar(value=40.0)

        self.snr_var = tk.DoubleVar(value=20.0)
        self.outdir_var = tk.StringVar(value=os.path.abspath("./"))

        self.create_widgets(main_frame)

    def create_widgets(self, parent):
        row = 0

        ttk.Label(parent, text="系统参数", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.add_entry(parent, "波长 (nm):", self.wavelength_var, row)
        row += 1
        self.add_entry(parent, "像素尺寸 dx (μm):", self.dx_var, row)
        row += 1
        self.add_entry(parent, "像素尺寸 dy (μm):", self.dy_var, row)
        row += 1
        self.add_entry(parent, "图像尺寸 N:", self.n_var, row)
        row += 1

        ttk.Separator(parent, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="粒子参数", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.add_entry(parent, "最小直径 (μm):", self.dmin_var, row)
        row += 1
        self.add_entry(parent, "最大直径 (μm):", self.dmax_var, row)
        row += 1
        self.add_entry(parent, "粒子数量:", self.num_var, row)
        row += 1

        ttk.Separator(parent, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="轴向分布 (Z轴)", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        self.add_entry(parent, "最小距离 Z_min (mm):", self.zmin_var, row)
        row += 1
        self.add_entry(parent, "最大距离 Z_max (mm):", self.zmax_var, row)
        row += 1
        self.add_entry(parent, "轴向层间距 (μm):", self.zstep_var, row)
        row += 1

        ttk.Separator(parent, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="噪声设置", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        ttk.Label(parent, text="(留空或输入 0 表示无噪声)").grid(row=row, column=1, sticky=tk.W, padx=5)
        row += 1
        self.add_entry(parent, "信噪比 SNR (dB):", self.snr_var, row)
        row += 1

        ttk.Separator(parent, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="ew", pady=10)
        row += 1

        ttk.Label(parent, text="导出设置", font=("Arial", 10, "bold")).grid(
            row=row, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        row += 1

        ttk.Label(parent, text="保存至:").grid(row=row, column=0, sticky=tk.W)
        ttk.Entry(parent, textvariable=self.outdir_var, state="readonly", width=30).grid(row=row, column=1, padx=5)
        ttk.Button(parent, text="浏览...", command=self.browse_dir).grid(row=row, column=2, sticky=tk.E)
        row += 1

        self.generate_btn = ttk.Button(parent, text="生成全息图及真值", command=self.start_generation)
        self.generate_btn.grid(row=row, column=0, columnspan=3, pady=20, ipadx=10, ipady=5)

        self.status_var = tk.StringVar(value="就绪")
        ttk.Label(parent, textvariable=self.status_var, foreground="gray").grid(
            row=row + 1, column=0, columnspan=3, sticky=tk.W
        )

    def add_entry(self, parent, label_text, variable, row):
        ttk.Label(parent, text=label_text).grid(row=row, column=0, sticky=tk.W, pady=2)
        ttk.Entry(parent, textvariable=variable, width=15).grid(row=row, column=1, sticky=tk.W, padx=5, pady=2)

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
            wl = self.wavelength_var.get() * 1e-9
            dx = self.dx_var.get() * 1e-6
            dy = self.dy_var.get() * 1e-6
            n = self.n_var.get()

            dmin = self.dmin_var.get() * 1e-6
            dmax = self.dmax_var.get() * 1e-6
            num = self.num_var.get()

            zmin = self.zmin_var.get() * 1e-3
            zmax = self.zmax_var.get() * 1e-3
            zstep = self.zstep_var.get() * 1e-6

            snr = self.snr_var.get()
            out_dir = self.outdir_var.get()

            prefix = self.get_next_prefix(out_dir)
            run_dir = os.path.join(out_dir, prefix.rstrip("_"))

            saved_dir = generate_hologram(
                wavelength=wl,
                dx=dx,
                dy=dy,
                N=n,
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

            self.root.after(0, self.generation_complete, True, (prefix, saved_dir))
        except Exception as exc:
            self.root.after(0, self.generation_complete, False, str(exc))

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
