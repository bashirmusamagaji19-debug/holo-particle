import numpy as np
from scipy.signal import butter, filtfilt
import matplotlib.pyplot as plt

# 设置 matplotlib 正常显示中文和负号
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']  # 兼容 Win 和 Mac
plt.rcParams['axes.unicode_minus'] = False


class SapphireHeterodyneInterferometer:
    """
    倒装LED芯片蓝宝石衬底外差干涉厚度测量仿真模型
    包含：信号生成、热光漂移模拟、偏振串扰注入、DLIA相位解调与可视化
    """

    def __init__(self):
        # 1. 物理参数与材料参数初始化
        self.wavelength = 632.8e-9  # 氦氖激光器真空波长 (m)
        self.n_sapphire_o = 1.768  # 蓝宝石寻常折射率 no (标称值 @632.8nm)
        self.dn_dT = 13.0e-6  # 蓝宝石热光系数 (K^-1)
        self.alpha_cte = 5.8e-6  # 蓝宝石c轴热膨胀系数 (K^-1)
        self.nominal_thickness = 100.0e-6  # 衬底标称几何厚度 100微米 (m)

        # 2. 系统信号与电子学参数
        self.f_beat = 2.0e6  # 外差拍频频率 2.0 MHz
        self.fs = 20.0e6  # ADC高频采样率 20.0 MHz
        self.t_duration = 10e-6  # 单次测量采样时间窗 10 us
        self.t = np.arange(0, self.t_duration, 1 / self.fs)

        # 3. DLIA 滤波器设计 (4阶Butterworth低通滤波器)
        nyquist = 0.5 * self.fs
        cutoff = 0.05 * self.f_beat  # 截止频率设为拍频的5% (100 kHz)
        self.b, self.a = butter(4, cutoff / nyquist, btype='low')

    def calculate_physical_state(self, delta_T):
        """ 依据温度波动 delta_T 计算真实的物理厚度与折射率 """
        # 热膨胀导致的真实厚度改变
        true_d = self.nominal_thickness * (1 + self.alpha_cte * delta_T)
        # 热光效应导致的折射率改变
        effective_n = self.n_sapphire_o + self.dn_dT * delta_T
        return true_d, effective_n

    def generate_heterodyne_signal(self, delta_T, crosstalk_ratio, snr_db):
        """
        生成混合了真实相位、偏振串扰、以及加性高斯白噪声的测量信号
        """
        true_d, effective_n = self.calculate_physical_state(delta_T)

        # 真实的光程差(OPD)导致的目标相位
        # 垂直反射式架构：OPD = 2 * n * d
        ideal_phase = (4 * np.pi * effective_n * true_d) / self.wavelength

        # 将相位包裹到 [-pi, pi] 以模拟干涉仪周期性特性
        ideal_phase_wrapped = (ideal_phase + np.pi) % (2 * np.pi) - np.pi

        # 光束主振幅
        A_main, B_main = 1.0, 1.0
        # 由偏振不理想产生的串扰寄生振幅
        a_cross = A_main * crosstalk_ratio
        b_cross = B_main * crosstalk_ratio

        # 信号合成（根据偏振串扰物理方程）
        # 项1：理想干涉信号
        sig_ideal = A_main * B_main * np.cos(2 * np.pi * self.f_beat * self.t + ideal_phase_wrapped)
        # 项2：混叠产生的零相位基底波动
        sig_crosstalk_dc = (A_main * b_cross + B_main * a_cross) * np.cos(2 * np.pi * self.f_beat * self.t)
        # 项3：混叠产生的负相位二次交调项
        sig_crosstalk_neg = a_cross * b_cross * np.cos(2 * np.pi * self.f_beat * self.t - ideal_phase_wrapped)

        clean_signal = sig_ideal + sig_crosstalk_dc + sig_crosstalk_neg

        # 根据给定的信噪比注入高斯白噪声 (AWGN)
        signal_power = np.mean(clean_signal ** 2)
        noise_power = signal_power / (10 ** (snr_db / 10))
        noise = np.random.normal(0, np.sqrt(noise_power), len(self.t))

        return clean_signal + noise, true_d, ideal_phase_wrapped

    def digital_lock_in_amplifier(self, signal):
        """
        数字锁相放大器算法：正交混频与低通滤波，提取相位
        """
        # 生成内部正交本地振荡信号
        lo_i = 2 * np.cos(2 * np.pi * self.f_beat * self.t)
        lo_q = 2 * np.sin(2 * np.pi * self.f_beat * self.t)

        # 正交混频
        mixed_i = signal * lo_i
        mixed_q = signal * lo_q

        # 零相移数字低通滤波
        filt_i = filtfilt(self.b, self.a, mixed_i)
        filt_q = filtfilt(self.b, self.a, mixed_q)

        # 规避滤波器边缘效应，取信号后半段稳定均值
        stable_idx = int(len(self.t) * 0.5)
        I_mean = np.mean(filt_i[stable_idx:])
        Q_mean = np.mean(filt_q[stable_idx:])

        # 相位重建
        extracted_phase = np.arctan2(Q_mean, I_mean)
        return extracted_phase

    def run_monte_carlo_validation(self, iterations=1000):
        """
        蒙特卡洛误差分析：随机化注入环境变量和硬件缺陷，输出不确定度
        """
        measurement_errors_um = []  # 初始化为空列表

        # 获取系统在理想标称条件下的基准参考相位
        _, _, base_phase = self.generate_heterodyne_signal(delta_T=0, crosstalk_ratio=0, snr_db=100)

        for _ in range(iterations):
            # 随机工况抽样
            d_T = np.random.uniform(-2.0, 2.0)  # 温度扰动：±2 K
            cross_ratio = np.random.uniform(0.005, 0.05)  # 偏振串扰：0.5% 至 5%
            snr = np.random.uniform(30, 50)  # 信噪比：30 dB 至 50 dB

            # 生成受干扰的实时信号
            noisy_signal, true_d, _ = self.generate_heterodyne_signal(d_T, cross_ratio, snr)

            # DLIA 相位解调
            meas_phase = self.digital_lock_in_amplifier(noisy_signal)

            # 计算测量相位相对于基准相位的偏移，处理包裹边界
            phase_diff = meas_phase - base_phase
            phase_diff = (phase_diff + np.pi) % (2 * np.pi) - np.pi

            # 采用标准常温折射率解算“表观厚度增量”
            calc_d_delta = (phase_diff * self.wavelength) / (4 * np.pi * self.n_sapphire_o)
            calculated_thickness = self.nominal_thickness + calc_d_delta

            # 计算绝对测量误差 = 表观厚度 - 真实的物理厚度
            abs_error_m = calculated_thickness - true_d
            measurement_errors_um.append(abs_error_m * 1e6)  # 转化为微米

        # 提取统计特征
        mean_error = np.mean(measurement_errors_um)
        std_uncertainty = np.std(measurement_errors_um)
        max_peak_error = np.max(np.abs(measurement_errors_um))

        return mean_error, std_uncertainty, max_peak_error, measurement_errors_um

    def plot_signal_and_demodulation(self, delta_T=0, crosstalk_ratio=0.02, snr_db=40):
        """
        可视化单次测量的时域波形与 DLIA 解调过程（纯中文界面）
        """
        # 生成对比信号：理想基准 vs 带噪受扰信号
        ideal_sig, _, _ = self.generate_heterodyne_signal(delta_T, 0, 100)
        noisy_sig, _, _ = self.generate_heterodyne_signal(delta_T, crosstalk_ratio, snr_db)

        # 提取 DLIA 的中间变量 (I/Q 信号) 用于展示
        lo_i = 2 * np.cos(2 * np.pi * self.f_beat * self.t)
        lo_q = 2 * np.sin(2 * np.pi * self.f_beat * self.t)
        filt_i = filtfilt(self.b, self.a, noisy_sig * lo_i)
        filt_q = filtfilt(self.b, self.a, noisy_sig * lo_q)
        stable_idx = int(len(self.t) * 0.5)

        plt.figure(figsize=(10, 8))

        # 子图 1：高频时域干涉信号
        plt.subplot(2, 1, 1)
        plt.plot(self.t * 1e6, ideal_sig, label='理想基准信号', color='black', alpha=0.6, linestyle='--')
        plt.plot(self.t * 1e6, noisy_sig, label=f'含噪受扰信号 (SNR={snr_db}dB)', color='C0', alpha=0.8)
        plt.title('外差干涉信号 (时域)')
        plt.xlabel('时间 (μs)')
        plt.ylabel('幅值')
        plt.legend(loc='upper right')
        plt.grid(True, linestyle=':', alpha=0.7)

        # 子图 2：DLIA 低通滤波后的正交 I/Q 信号
        plt.subplot(2, 1, 2)
        plt.plot(self.t * 1e6, filt_i, label='滤波后的 I 路 (同相)', color='C1', linewidth=2)
        plt.plot(self.t * 1e6, filt_q, label='滤波后的 Q 路 (正交)', color='C2', linewidth=2)

        # 标注有效采样的稳定区间
        plt.axvspan(self.t[stable_idx] * 1e6, self.t[-1] * 1e6, color='gray', alpha=0.2, label='稳定平均窗口')
        plt.title('DLIA 解调：I 路与 Q 路信号')
        plt.xlabel('时间 (μs)')
        plt.ylabel('幅值')
        plt.legend(loc='lower right')
        plt.grid(True, linestyle=':', alpha=0.7)

        plt.tight_layout()
        plt.show()

    def plot_error_distribution(self, errors_um):
        """
        绘制蒙特卡洛误差分布直方图（纯中文界面）
        """
        plt.figure(figsize=(8, 5))

        # 绘制直方图
        counts, bins, patches = plt.hist(errors_um, bins=50, color='royalblue', edgecolor='black', alpha=0.7)

        # 标记均值和 1σ 范围
        mean_err = np.mean(errors_um)
        std_err = np.std(errors_um)

        plt.axvline(mean_err, color='red', linestyle='dashed', linewidth=2, label=f'均值: {mean_err:.4f} μm')
        plt.axvline(mean_err + std_err, color='orange', linestyle='dotted', linewidth=2,
                    label=f'+1σ: {mean_err + std_err:.4f} μm')
        plt.axvline(mean_err - std_err, color='orange', linestyle='dotted', linewidth=2,
                    label=f'-1σ: {mean_err - std_err:.4f} μm')

        plt.title('蒙特卡洛测量误差分布')
        plt.xlabel('绝对厚度误差 (μm)')
        plt.ylabel('频数')
        plt.legend()
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    # 实例化仿真器
    simulator = SapphireHeterodyneInterferometer()

    # 1. 选一个典型工况展示波形和解调过程
    print("正在生成单次测量的时域波形图...")
    simulator.plot_signal_and_demodulation(delta_T=1.0, crosstalk_ratio=0.03, snr_db=35)

    # 2. 执行蒙特卡洛验证
    N_iters = 5000
    print(f"\n正在执行 {N_iters} 次蒙特卡洛仿真计算，请稍候...")
    mean_err, std_dev, max_err, errors_array = simulator.run_monte_carlo_validation(iterations=N_iters)

    print(f"\n--- 蓝宝石衬底厚度外差干涉仿真结果 ({N_iters} 次迭代) ---")
    print(f"标称厚度: 100.0 μm")
    print(f"误差均值: {mean_err:>.6f} μm")
    print(f"标准不确定度 (1σ): {std_dev:>.6f} μm")
    print(f"最大峰值误差: {max_err:>.6f} μm")

    if max_err < 0.1:
        print("结论: 验证通过，系统完全满足 0.1 μm 的不确定度测量要求。")
    else:
        print("结论: 验证失败，误差超出限度。")

    # 3. 绘制误差分布图
    print("正在生成误差分布直方图...")
    simulator.plot_error_distribution(errors_array)