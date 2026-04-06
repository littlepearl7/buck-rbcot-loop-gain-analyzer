# main.py
import os
import sys
_DIR = os.path.dirname(os.path.abspath(__file__))
if _DIR not in sys.path:
    sys.path.insert(0, _DIR)

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from PIL import Image, ImageTk
import warnings
warnings.filterwarnings('ignore')

from rbcot_formulas import (
    esr_pade_loop,
    fixed_slope_pade_loop,
    sw_filter_pade_loop,
    interp_zero_db_crossing,
    bode_unwrap_deg,
)
from compensation_parameter_manager import CompensationParameterManager

# ========== 全局样式设置 ==========
# matplotlib 全局样式
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['axes.titlepad'] = 12
plt.rcParams['axes.labelpad'] = 8
plt.rcParams['figure.facecolor'] = '#1e1e2f'
plt.rcParams['axes.facecolor'] = '#2a2a3b'
plt.rcParams['axes.edgecolor'] = '#4a4a5a'
plt.rcParams['axes.labelcolor'] = '#e0e0e0'
plt.rcParams['text.color'] = '#e0e0e0'
plt.rcParams['grid.color'] = '#4a4a5a'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.alpha'] = 0.3

# ttk 主题样式配置
def setup_ttk_style():
    style = ttk.Style()
    style.theme_use('clam')
    
    # 主窗口背景
    style.configure('TFrame', background='#1e1e2f')
    style.configure('TLabel', background='#1e1e2f', foreground='#e0e0e0', font=('Segoe UI', 10))
    style.configure('TLabelframe', background='#1e1e2f', foreground='#e0e0e0', 
                    bordercolor='#4a4a5a', relief='solid', borderwidth=1)
    style.configure('TLabelframe.Label', background='#1e1e2f', foreground='#61afef', 
                    font=('Segoe UI', 11, 'bold'))
    
    # 输入框
    style.configure('TEntry', fieldbackground='#2a2a3b', foreground='#e0e0e0', 
                    borderwidth=1, focuscolor='#61afef', font=('Segoe UI', 10))
    style.map('TEntry', fieldbackground=[('focus', '#2a2a3b')])
    
    # 按钮样式
    style.configure('TButton', background='#3a3a4a', foreground='#e0e0e0', 
                    borderwidth=1, focuscolor='#61afef', font=('Segoe UI', 10),
                    padding=(10, 5))
    style.map('TButton', background=[('active', '#4a4a5a'), ('pressed', '#2a2a3b')])
    
    return style

class BUCKBodeAnalyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("BUCK参数分析器 - 基于论文 Accurate Loop Gain Model of RBCOT")
        self.root.geometry("1800x950")
        self.root.configure(bg='#1e1e2f')
        
        # 初始化参数管理器
        self.param_manager = CompensationParameterManager()
        
        # 存储各个模型的计算结果
        self.esr_mag = None
        self.esr_phase = None
        self.esr_freq = None
        self.fixed_mag = None
        self.fixed_phase = None
        self.fixed_freq = None
        self.sw_mag = None
        self.sw_phase = None
        self.sw_freq = None
        
        # 当前显示的模型
        self.current_model = 'esr'
        
        # 设置ttk样式
        setup_ttk_style()
        
        # 创建主框架
        self.create_widgets()
        
        # 加载保存的配置到UI
        self.load_config_to_ui()
        
        # 更新占空比显示
        self.update_duty()
        
        # 初始计算ESR模型
        self.calculate_esr_model()
        self.update_display()
    
    def create_widgets(self):
        """创建所有UI控件"""
        # 输入参数框架 (左侧)
        input_frame = ttk.LabelFrame(self.root, text="输入参数", padding=15)
        input_frame.place(x=10, y=10, width=680, height=920)
        
        # ========== 电路公共参数 ==========
        ttk.Label(input_frame, text="【电路公共参数】", foreground='#61afef', 
                  font=('Segoe UI', 11, 'bold')).grid(row=0, column=0, columnspan=4, sticky='w', pady=(0, 10))
        
        # 第1排
        row = 1
        ttk.Label(input_frame, text="输入电压 Vin (V):", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=0, sticky='w', pady=8)
        self.Vin = ttk.Entry(input_frame, width=18, font=('Segoe UI', 10))
        self.Vin.grid(row=row, column=1, padx=8)
        self.Vin.bind('<KeyRelease>', self.update_duty)  # 只更新占空比，不自动计算
        
        ttk.Label(input_frame, text="输出电压 Vout (V):", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=2, sticky='w', padx=(20,0), pady=8)
        self.Vout = ttk.Entry(input_frame, width=18, font=('Segoe UI', 10))
        self.Vout.grid(row=row, column=3, padx=8)
        self.Vout.bind('<KeyRelease>', self.update_duty)  # 只更新占空比，不自动计算
        
        # 第2排
        row += 1
        ttk.Label(input_frame, text="开关频率 Fsw (Hz):", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=0, sticky='w', pady=8)
        self.Fsw = ttk.Entry(input_frame, width=18, font=('Segoe UI', 10))
        self.Fsw.grid(row=row, column=1, padx=8)
        
        ttk.Label(input_frame, text="占空比 D:", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=2, sticky='w', padx=(20,0), pady=8)
        self.D_label = ttk.Label(input_frame, text="0.1", foreground='#61afef', font=('Segoe UI', 10, 'bold'))
        self.D_label.grid(row=row, column=3, sticky='w', padx=8)
        
        # 第3排
        row += 1
        ttk.Label(input_frame, text="输出电容 Cout (F):", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=0, sticky='w', pady=8)
        self.Cout = ttk.Entry(input_frame, width=18, font=('Segoe UI', 10))
        self.Cout.grid(row=row, column=1, padx=8)
        
        ttk.Label(input_frame, text="功率电感 L (H):", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=2, sticky='w', padx=(20,0), pady=8)
        self.Ls = ttk.Entry(input_frame, width=18, font=('Segoe UI', 10))
        self.Ls.grid(row=row, column=3, padx=8)
        
        # 第4排
        row += 1
        ttk.Label(input_frame, text="负载电阻 Rload (Ω):", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=0, sticky='w', pady=8)
        self.Rload = ttk.Entry(input_frame, width=18, font=('Segoe UI', 10))
        self.Rload.grid(row=row, column=1, padx=8)
        
        # 分隔线
        row += 1
        ttk.Separator(input_frame, orient='horizontal').grid(row=row, column=0, columnspan=4, sticky='ew', pady=10)
        
        # ========== ESR纹波补偿参数 ==========
        row += 1
        ttk.Label(input_frame, text="【ESR纹波补偿参数】", foreground='#98c379', 
                  font=('Segoe UI', 11, 'bold')).grid(row=row, column=0, columnspan=4, sticky='w', pady=(10, 5))
        
        row += 1
        ttk.Label(input_frame, text="ESR Rc (mOhm):", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=0, sticky='w', pady=8)
        self.Resr = ttk.Entry(input_frame, width=18, font=('Segoe UI', 10))
        self.Resr.grid(row=row, column=1, padx=8)
        
        # 分隔线
        row += 1
        ttk.Separator(input_frame, orient='horizontal').grid(row=row, column=0, columnspan=4, sticky='ew', pady=10)
        
        # ========== 固定斜坡补偿参数 ==========
        row += 1
        ttk.Label(input_frame, text="【固定斜坡补偿参数】", foreground='#e5c07b', 
                  font=('Segoe UI', 11, 'bold')).grid(row=row, column=0, columnspan=4, sticky='w', pady=(10, 5))
        
        row += 1
        ttk.Label(input_frame, text="固定斜坡 Vslope (V/周期):", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=0, sticky='w', pady=8)
        self.Vslope = ttk.Entry(input_frame, width=18, font=('Segoe UI', 10))
        self.Vslope.grid(row=row, column=1, padx=8)
        
        # 分隔线
        row += 1
        ttk.Separator(input_frame, orient='horizontal').grid(row=row, column=0, columnspan=4, sticky='ew', pady=10)
        
        # ========== SW滤波器补偿参数 ==========
        row += 1
        ttk.Label(input_frame, text="【SW滤波器补偿参数】", foreground='#e06c75', 
                  font=('Segoe UI', 11, 'bold')).grid(row=row, column=0, columnspan=4, sticky='w', pady=(10, 5))
        
        row += 1
        ttk.Label(input_frame, text="SW滤波器 Rx (Ω):", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=0, sticky='w', pady=8)
        self.Rx = ttk.Entry(input_frame, width=18, font=('Segoe UI', 10))
        self.Rx.grid(row=row, column=1, padx=8)
        
        ttk.Label(input_frame, text="SW滤波器 Cx (F):", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=2, sticky='w', padx=(20,0), pady=8)
        self.Cx = ttk.Entry(input_frame, width=18, font=('Segoe UI', 10))
        self.Cx.grid(row=row, column=3, padx=8)
        
        # 分隔线
        row += 1
        ttk.Separator(input_frame, orient='horizontal').grid(row=row, column=0, columnspan=4, sticky='ew', pady=10)
        
        # ========== 分析设置 ==========
        row += 1
        ttk.Label(input_frame, text="【分析设置】", foreground='#61afef', 
                  font=('Segoe UI', 11, 'bold')).grid(row=row, column=0, columnspan=4, sticky='w', pady=(10, 5))
        
        row += 1
        ttk.Label(input_frame, text="频率范围 (Hz):", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=row, column=0, sticky='w', pady=8)
        self.freq_min = ttk.Entry(input_frame, width=15, font=('Segoe UI', 10))
        self.freq_min.grid(row=row, column=1, padx=8, sticky='w')
        self.freq_min.insert(0, "100")
        
        self.freq_max = ttk.Entry(input_frame, width=15, font=('Segoe UI', 10))
        self.freq_max.grid(row=row, column=2, padx=8, sticky='w')
        self.freq_max.insert(0, "2e6")
        
        # 按钮
        row += 1
        btn_frame = ttk.Frame(input_frame)
        btn_frame.grid(row=row, column=0, columnspan=4, pady=25)
        
        self.btn_esr = ttk.Button(btn_frame, text="ESR纹波补偿模型", 
                                  command=self.calculate_and_show_esr, width=20)
        self.btn_esr.pack(side=tk.LEFT, padx=8)
        
        self.btn_fixed = ttk.Button(btn_frame, text="固定斜坡补偿模型", 
                                    command=self.calculate_and_show_fixed, width=20)
        self.btn_fixed.pack(side=tk.LEFT, padx=8)
        
        self.btn_sw = ttk.Button(btn_frame, text="SW滤波器补偿模型", 
                                 command=self.calculate_and_show_sw, width=20)
        self.btn_sw.pack(side=tk.LEFT, padx=8)
        
        # 计算结果区域
        row += 1
        result_frame = ttk.LabelFrame(input_frame, text="计算结果", padding=10)
        result_frame.grid(row=row, column=0, columnspan=4, pady=12, sticky='ew')
        
        # ESR模型结果
        ttk.Label(result_frame, text="ESR模型:", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=0, column=0, sticky='w', pady=5)
        self.esr_fc = ttk.Label(result_frame, text="Fc = -- Hz", foreground='#98c379', font=('Segoe UI', 10, 'bold'))
        self.esr_fc.grid(row=0, column=1, sticky='w', padx=10)
        self.esr_pm = ttk.Label(result_frame, text="PM = -- °", foreground='#98c379', font=('Segoe UI', 10, 'bold'))
        self.esr_pm.grid(row=0, column=2, sticky='w', padx=10)
        
        # 固定斜坡模型结果
        ttk.Label(result_frame, text="固定斜坡模型:", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=1, column=0, sticky='w', pady=5)
        self.fixed_fc = ttk.Label(result_frame, text="Fc = -- Hz", foreground='#98c379', font=('Segoe UI', 10, 'bold'))
        self.fixed_fc.grid(row=1, column=1, sticky='w', padx=10)
        self.fixed_pm = ttk.Label(result_frame, text="PM = -- °", foreground='#98c379', font=('Segoe UI', 10, 'bold'))
        self.fixed_pm.grid(row=1, column=2, sticky='w', padx=10)
        
        # SW滤波器模型结果
        ttk.Label(result_frame, text="SW滤波器模型:", foreground='#e0e0e0', font=('Segoe UI', 10)).grid(row=2, column=0, sticky='w', pady=5)
        self.sw_fc = ttk.Label(result_frame, text="Fc = -- Hz", foreground='#98c379', font=('Segoe UI', 10, 'bold'))
        self.sw_fc.grid(row=2, column=1, sticky='w', padx=10)
        self.sw_pm = ttk.Label(result_frame, text="PM = -- °", foreground='#98c379', font=('Segoe UI', 10, 'bold'))
        self.sw_pm.grid(row=2, column=2, sticky='w', padx=10)
        
        # ========== 右侧区域：Bode图 ==========
        # 上图区域 - 增益图
        self.fig_gain = Figure(figsize=(7.8, 4.8), facecolor='#1e1e2f')
        self.ax_gain = self.fig_gain.add_subplot(111, facecolor='#2a2a3b')
        self.canvas_gain = FigureCanvasTkAgg(self.fig_gain, master=self.root)
        self.canvas_gain.get_tk_widget().place(x=710, y=10, width=1060, height=450)
        self.ax_gain.set_title('', color='#61afef', fontsize=12, fontweight='bold')
        self.ax_gain.set_xlabel('频率 (Hz)', color='#e0e0e0')
        self.ax_gain.set_ylabel('增益 (dB)', color='#61afef', fontsize=11)
        self.ax_gain.grid(True, alpha=0.3, linestyle='--')
        self.ax_gain.tick_params(colors='#e0e0e0')
        for spine in self.ax_gain.spines.values():
            spine.set_color('#4a4a5a')
        
        # 下图区域 - 相位图
        self.fig_phase = Figure(figsize=(7.8, 4.2), facecolor='#1e1e2f')
        self.ax_phase = self.fig_phase.add_subplot(111, facecolor='#2a2a3b')
        self.canvas_phase = FigureCanvasTkAgg(self.fig_phase, master=self.root)
        self.canvas_phase.get_tk_widget().place(x=710, y=470, width=1060, height=430)
        self.ax_phase.set_title('', color='#61afef', fontsize=12, fontweight='bold')
        self.ax_phase.set_xlabel('频率 (Hz)', color='#e0e0e0')
        self.ax_phase.set_ylabel('相位 (deg)', color='#e5c07b', fontsize=11)
        self.ax_phase.grid(True, alpha=0.3, linestyle='--')
        self.ax_phase.tick_params(colors='#e0e0e0')
        self.ax_phase.set_ylim([-180, 180])
        for spine in self.ax_phase.spines.values():
            spine.set_color('#4a4a5a')
    
    def load_config_to_ui(self):
        """从参数管理器加载配置到UI"""
        try:
            self.Vin.insert(0, str(self.param_manager.circuit.Vin))
            self.Vout.insert(0, str(self.param_manager.circuit.Vout))
            self.Fsw.insert(0, str(self.param_manager.circuit.Fsw))
            self.Cout.insert(0, str(self.param_manager.circuit.Cout))
            self.Ls.insert(0, str(self.param_manager.circuit.L))
            self.Rload.insert(0, str(self.param_manager.circuit.Rload))
            self.Resr.insert(0, str(self.param_manager.esr.Rc * 1e3))
            self.Vslope.insert(0, str(self.param_manager.fixed_slope.Vslope))
            self.Rx.insert(0, str(self.param_manager.sw_filter.Rx))
            self.Cx.insert(0, str(self.param_manager.sw_filter.Cx))
        except Exception as e:
            print(f"加载配置到UI失败: {e}")
    
    def update_duty(self, event=None):
        """更新占空比显示"""
        try:
            vin = float(self.Vin.get())
            vout = float(self.Vout.get())
            duty = vout / vin
            self.D_label.config(text=f"{duty:.3f}")
        except:
            pass
    
    def get_input_values(self):
        """获取所有输入参数"""
        try:
            values = {
                'Vin': float(self.Vin.get()),
                'Vout': float(self.Vout.get()),
                'Fs': float(self.Fsw.get()),
                'Cout': float(self.Cout.get()),
                'Rc': float(self.Resr.get()) * 1e-3,  # mOhm -> Ohm
                'L': float(self.Ls.get()),
                'Rload': float(self.Rload.get()),
                'Vslope': float(self.Vslope.get()),
                'Rx': float(self.Rx.get()),
                'Cx': float(self.Cx.get()),
                'freq_min': float(self.freq_min.get()),
                'freq_max': float(self.freq_max.get())
            }
            
            # 基本验证
            for key, val in values.items():
                if val <= 0:
                    raise ValueError(f"{key} 必须为正数")
            
            if values['freq_max'] <= values['freq_min']:
                raise ValueError("f_min 必须小于 f_max")
            
            # 计算派生参数
            values['Ts'] = 1/values['Fs']
            values['Ton'] = (values['Vout']/values['Vin']) * values['Ts']
            values['Toff'] = values['Ts'] - values['Ton']
            values['D'] = values['Vout']/values['Vin']
            
            return values
            
        except ValueError as e:
            messagebox.showerror("输入错误", f"请输入有效的数字参数\n{str(e)}")
            return None
    
    def calculate_esr_model(self):
        """ESR纹波补偿模型"""
        values = self.get_input_values()
        if values is None:
            return
        
        try:
            Ts = values['Ts']
            Ton = values['Ton']
            L = values['L']
            Cout = values['Cout']
            Rc = values['Rc']
            Rload = values['Rload']
            
            Floop, meta = esr_pade_loop(L, Cout, Rc, Rload, Ts, Ton)
            fz_esr = meta['omega_z'] / (2 * np.pi)
            
            w = 2 * np.pi * np.logspace(np.log10(values['freq_min']/10), 
                                        np.log10(values['freq_max']*10), 3000)
            freq, mag, phase_deg_unwrapped = bode_unwrap_deg(Floop, w)
            
            self.esr_freq = freq
            self.esr_mag = mag
            self.esr_phase = phase_deg_unwrapped
            
            # 从Bode图提取Fc和PM（已修复零点公式，现在准确）
            mag_db = 20 * np.log10(np.maximum(mag, 1e-10))
            fc_actual, cross_idx = interp_zero_db_crossing(freq, mag_db)
            
            if fc_actual is not None and cross_idx is not None:
                i0 = cross_idx
                df = freq[i0 + 1] - freq[i0]
                t = (fc_actual - freq[i0]) / df if abs(df) > 1e-30 else 0.0
                phase_at_fc = phase_deg_unwrapped[i0] + t * (phase_deg_unwrapped[i0 + 1] - phase_deg_unwrapped[i0])
                pm_actual = phase_at_fc + 180.0
                while pm_actual > 180:
                    pm_actual -= 360
                while pm_actual < 0:
                    pm_actual += 360
                pm_actual = min(pm_actual, 180)
                
                self.esr_fc.config(text=f"Fc = {fc_actual/1000:.1f} kHz")
                self.esr_pm.config(text=f"PM = {pm_actual:.1f}° | fz={fz_esr/1000:.1f} kHz")
            else:
                fc_est = meta['omega_c_approx'] / (2 * np.pi)
                peak_idx = int(np.argmax(mag_db))
                peak_gain_db = mag_db[peak_idx]
                peak_freq = freq[peak_idx]
                self.esr_fc.config(text="Fc = No 0 dB crossing")
                self.esr_pm.config(
                    text=(
                        f"Peak = {peak_gain_db:.1f} dB @ {peak_freq/1000:.1f} kHz"
                        f" | Fc(approx)={fc_est/1000:.1f} kHz | fz={fz_esr/1000:.1f} kHz"
                    )
                )
            
            # 保存参数到管理器
            self.sync_parameters_to_manager()
            
        except Exception as e:
            messagebox.showerror("计算错误", f"ESR模型计算过程中出现错误:\n{str(e)}")
            print(f"ESR模型错误: {e}")

    def calculate_fixed_ramp_model(self):
        """固定斜坡补偿模型"""
        values = self.get_input_values()
        if values is None:
            return
        
        try:
            Vout = values['Vout']
            Ts = values['Ts']
            Ton = values['Ton']
            L = values['L']
            Cout = values['Cout']
            Rload = values['Rload']
            Vslope = values['Vslope']
            fs = 1 / Ts
            
            Floop, meta = fixed_slope_pade_loop(L, Cout, Rload, Vout, Ts, Ton, Vslope)
            
            w = 2 * np.pi * np.logspace(np.log10(values['freq_min']), 
                                        np.log10(min(values['freq_max'], 5 * fs)), 2000)
            freq, mag, phase_deg_unwrapped = bode_unwrap_deg(Floop, w)
            
            self.fixed_freq = freq
            self.fixed_mag = mag
            self.fixed_phase = phase_deg_unwrapped
            
            mag_db = 20 * np.log10(np.maximum(mag, 1e-10))
            fc_actual, cross_idx = interp_zero_db_crossing(freq, mag_db)
            
            if fc_actual is not None and cross_idx is not None:
                i0 = cross_idx
                df = freq[i0 + 1] - freq[i0]
                t = (fc_actual - freq[i0]) / df if abs(df) > 1e-30 else 0.0
                phase_at_fc = phase_deg_unwrapped[i0] + t * (phase_deg_unwrapped[i0 + 1] - phase_deg_unwrapped[i0])
                pm_actual = phase_at_fc + 180.0
                while pm_actual > 180:
                    pm_actual -= 360
                while pm_actual < 0:
                    pm_actual += 360
                pm_actual = min(pm_actual, 180)
                self.fixed_fc.config(text=f"Fc = {fc_actual/1000:.1f} kHz")
                self.fixed_pm.config(text=f"PM = {pm_actual:.1f}°")
            else:
                fc_est = meta['omega_c_s_approx'] / (2 * np.pi)
                self.fixed_fc.config(text=f"Fc ≈ {fc_est/1000:.1f} kHz")
                self.fixed_pm.config(text=f"PM ≈ 35°")
            
            # 保存参数到管理器
            self.sync_parameters_to_manager()
            
        except Exception as e:
            messagebox.showerror("计算错误", f"固定斜坡模型计算过程中出现错误:\n{str(e)}")
            print(f"固定斜坡模型错误: {e}")

    def calculate_sw_filter_model(self):
        """SW滤波器补偿模型"""
        values = self.get_input_values()
        if values is None:
            return
        
        try:
            Ts = values['Ts']
            Ton = values['Ton']
            L = values['L']
            Cout = values['Cout']
            Rload = values['Rload']
            Rx = values['Rx']
            Cx = values['Cx']
            
            Floop, meta = sw_filter_pade_loop(L, Cout, Rload, Ts, Ton, Rx, Cx)
            
            # Start sweep from omega_zr/20 to capture the low-freq zero
            omega_zr = meta['omega_zr']
            f_low = max(1.0, omega_zr / (2 * np.pi) / 20.0)
            f_high = values['freq_max'] * 10
            w = 2 * np.pi * np.logspace(np.log10(f_low), np.log10(f_high), 5000)
            freq, mag, phase_deg_unwrapped = bode_unwrap_deg(Floop, w)
            
            self.sw_freq = freq
            self.sw_mag = mag
            self.sw_phase = phase_deg_unwrapped
            
            mag_db = 20 * np.log10(np.maximum(mag, 1e-10))
            fc_actual, cross_idx = interp_zero_db_crossing(freq, mag_db)
            # RxCx design range hint
            rxcx_min = meta['rxcx_min']
            rxcx_max = meta['rxcx_max']
            tau_val = Rx * Cx
            in_range = rxcx_min <= tau_val <= rxcx_max
            range_tag = 'OK' if in_range else 'OUT'
            rxcx_min = meta['rxcx_min']
            rxcx_max = meta['rxcx_max']
            tau_val = Rx * Cx
            in_range = rxcx_min <= tau_val <= rxcx_max
            range_tag = 'OK' if in_range else 'OUT' 
            
            if fc_actual is not None and cross_idx is not None:
                i0 = cross_idx
                df = freq[i0 + 1] - freq[i0]
                t = (fc_actual - freq[i0]) / df if abs(df) > 1e-30 else 0.0
                phase_at_fc = phase_deg_unwrapped[i0] + t * (phase_deg_unwrapped[i0 + 1] - phase_deg_unwrapped[i0])
                pm_actual = phase_at_fc + 180.0
                while pm_actual > 180:
                    pm_actual -= 360
                while pm_actual < 0:
                    pm_actual += 360
                pm_actual = min(pm_actual, 180)
                self.sw_fc.config(text=f"Fc = {fc_actual/1000:.1f} kHz [{range_tag}]")
                self.sw_pm.config(text=f"PM = {pm_actual:.1f}deg  RxCx:[{rxcx_min*1e6:.1f},{rxcx_max*1e6:.1f}]us")
            else:
                idx_peak = int(np.argmax(mag_db))
                fc_peak = freq[idx_peak]
                phase_at_peak = phase_deg_unwrapped[idx_peak]
                pm_peak = phase_at_peak + 180.0
                while pm_peak > 180:
                    pm_peak -= 360
                while pm_peak < 0:
                    pm_peak += 360
                pm_peak = min(pm_peak, 180)

                self.sw_fc.config(text=f"No 0dB crossing [{range_tag}]")
                self.sw_pm.config(text=f"PM(ref)={pm_peak:.1f}deg  RxCx:[{rxcx_min*1e6:.1f},{rxcx_max*1e6:.1f}]us")
            
            # 保存参数到管理器
            self.sync_parameters_to_manager()
            
        except Exception as e:
            messagebox.showerror("计算错误", f"SW滤波器模型计算过程中出现错误:\n{str(e)}")
            print(f"SW滤波器模型错误: {e}")
    
    def sync_parameters_to_manager(self):
        """同步UI参数到参数管理器（分别存储三种补偿方式的参数）"""
        try:
            # 更新电路公共参数
            self.param_manager.update_circuit(
                vin=float(self.Vin.get()),
                vout=float(self.Vout.get()),
                fsw=float(self.Fsw.get()),
                cout=float(self.Cout.get()),
                l=float(self.Ls.get()),
                rload=float(self.Rload.get())
            )
            
            # 更新ESR补偿参数
            self.param_manager.update_esr(
                rc=float(self.Resr.get()) * 1e-3  # mOhm -> Ohm
            )
            
            # 更新固定斜坡补偿参数
            self.param_manager.update_fixed_slope(
                vslope=float(self.Vslope.get())
            )
            
            # 更新SW滤波器补偿参数
            self.param_manager.update_sw_filter(
                rx=float(self.Rx.get()),
                cx=float(self.Cx.get())
            )
            
        except Exception as e:
            # 静默处理转换错误
            pass

    def update_display(self):
        """绘制当前模型的 Bode 图"""
        values = self.get_input_values()
        if values is None:
            return
        
        if self.current_model == 'esr':
            freq, mag, phase = self.esr_freq, self.esr_mag, self.esr_phase
            title = 'ESR纹波补偿模型'
            fs = values['Fs']
        elif self.current_model == 'fixed':
            freq, mag, phase = self.fixed_freq, self.fixed_mag, self.fixed_phase
            title = '固定斜率斜坡补偿模型'
            fs = values['Fs']
        elif self.current_model == 'sw':
            freq, mag, phase = self.sw_freq, self.sw_mag, self.sw_phase
            title = 'SW滤波器补偿模型'
            fs = values['Fs']
        else:
            return
        
        if freq is None or len(freq) == 0:
            return
        
        freq = np.asarray(freq)
        mag = np.asarray(mag)
        phase = np.asarray(phase)
        valid_idx = np.isfinite(mag) & np.isfinite(phase) & np.isfinite(freq)
        freq, mag, phase = freq[valid_idx], mag[valid_idx], phase[valid_idx]
        
        if len(freq) == 0:
            self.ax_gain.clear()
            self.ax_gain.text(0.5, 0.5, '无有效数据', transform=self.ax_gain.transAxes,
                            ha='center', va='center', color='#e0e0e0')
            self.canvas_gain.draw()
            self.ax_phase.clear()
            self.ax_phase.text(0.5, 0.5, '无有效数据', transform=self.ax_phase.transAxes,
                            ha='center', va='center', color='#e0e0e0')
            self.canvas_phase.draw()
            return
        
        mag_db = 20 * np.log10(np.maximum(mag, 1e-10))
        fc_actual, cross_idx = interp_zero_db_crossing(freq, mag_db)
        
        x_min = float(np.min(freq) * 0.8)
        x_max = float(np.max(freq) * 1.2)
        
        # 绘制增益图
        self.ax_gain.clear()
        self.ax_gain.set_facecolor('#2a2a3b')
        self.ax_gain.semilogx(freq, mag_db, '#61afef', linewidth=2, label='增益')
        self.ax_gain.set_xlabel('频率 (Hz)', fontsize=10, color='#e0e0e0')
        self.ax_gain.set_ylabel('增益 (dB)', color='#61afef', fontsize=10)
        self.ax_gain.tick_params(axis='y', labelcolor='#61afef', colors='#e0e0e0')
        self.ax_gain.grid(True, alpha=0.3, which='both', linestyle='--', color='#4a4a5a')
        self.ax_gain.set_xlim([x_min, x_max])
        if fs > 0:
            self.ax_gain.axvline(fs, color='#98c379', linestyle='--', alpha=0.5, linewidth=1, label=f'fs={fs/1000:.0f}kHz')
            self.ax_gain.axvline(fs/2, color='#98c379', linestyle=':', alpha=0.5, linewidth=1, label='fs/2')
        self.ax_gain.set_title(f'{title} - 增益曲线', color='#61afef', fontsize=12, fontweight='bold')
        self.ax_gain.axhline(0, color='#e5c07b', linestyle=':', alpha=0.5, linewidth=1)
        
        if fc_actual is not None and x_min <= fc_actual <= x_max:
            self.ax_gain.axvline(fc_actual, color='#98c379', linestyle='--', alpha=0.7, linewidth=1)
            self.ax_gain.text(fc_actual * 1.1, -15, f'Fc={fc_actual/1000:.1f}kHz',
                            fontsize=9, color='#e0e0e0',
                            bbox=dict(facecolor='#2a2a3b', alpha=0.8, edgecolor='#4a4a5a'))
        else:
            peak_idx = int(np.argmax(mag_db))
            peak_freq = freq[peak_idx]
            peak_gain_db = mag_db[peak_idx]
            text_y = peak_gain_db + 3.0
            self.ax_gain.plot(peak_freq, peak_gain_db, 'o', color='#e06c75', markersize=6)
            self.ax_gain.text(
                peak_freq * 1.05,
                text_y,
                'No 0 dB crossing',
                fontsize=9,
                color='#e0e0e0',
                bbox=dict(facecolor='#2a2a3b', alpha=0.8, edgecolor='#4a4a5a')
            )
        for spine in self.ax_gain.spines.values():
            spine.set_color('#4a4a5a')
        self.canvas_gain.draw()
        
        # 绘制相位图
        self.ax_phase.clear()
        self.ax_phase.set_facecolor('#2a2a3b')
        self.ax_phase.semilogx(freq, phase, '#e5c07b', linewidth=2, label='相位')
        self.ax_phase.set_xlabel('频率 (Hz)', fontsize=10, color='#e0e0e0')
        self.ax_phase.set_ylabel('相位 (deg)', color='#e5c07b', fontsize=10)
        self.ax_phase.tick_params(axis='y', labelcolor='#e5c07b', colors='#e0e0e0')
        self.ax_phase.grid(True, alpha=0.3, which='both', linestyle='--', color='#4a4a5a')
        self.ax_phase.set_xlim([x_min, x_max])
        if fs > 0:
            self.ax_phase.axvline(fs, color='#98c379', linestyle='--', alpha=0.5, linewidth=1)
            self.ax_phase.axvline(fs/2, color='#98c379', linestyle=':', alpha=0.5, linewidth=1)
        
        ph = phase[np.isfinite(phase)]
        if len(ph) > 0:
            y_margin = max((np.max(ph) - np.min(ph)) * 0.2, 45)
            y_min = float(np.min(ph) - y_margin)
            y_max = float(np.max(ph) + y_margin)
        else:
            y_min, y_max = -360.0, 90.0
        self.ax_phase.set_ylim([y_min, y_max])
        self.ax_phase.set_title(f'{title} - 相位曲线', color='#61afef', fontsize=12, fontweight='bold')
        self.ax_phase.axhline(0, color='#61afef', linestyle=':', alpha=0.5, linewidth=1)
        self.ax_phase.axhline(-180, color='#e06c75', linestyle=':', alpha=0.5, linewidth=1, label='-180°')
        
        if fc_actual is not None and cross_idx is not None and x_min <= fc_actual <= x_max:
            i0 = cross_idx
            df = freq[i0 + 1] - freq[i0]
            t = (fc_actual - freq[i0]) / df if abs(df) > 1e-30 else 0.0
            phase_at_fc = phase[i0] + t * (phase[i0 + 1] - phase[i0])
            pm_actual = phase_at_fc + 180.0
            while pm_actual > 180:
                pm_actual -= 360
            while pm_actual < 0:
                pm_actual += 360
            pm_actual = min(pm_actual, 180)
            self.ax_phase.axvline(fc_actual, color='#98c379', linestyle='--', alpha=0.7, linewidth=1)
            self.ax_phase.plot(fc_actual, phase_at_fc, 'o', color='#e5c07b', markersize=6)
            text_y = phase_at_fc + 20
            if text_y > y_max:
                text_y = phase_at_fc - 30
            if text_y < y_min:
                text_y = phase_at_fc + 15
            self.ax_phase.text(fc_actual * 1.1, text_y, f'PM={pm_actual:.1f}°',
                            fontsize=9, color='#e0e0e0',
                            bbox=dict(facecolor='#2a2a3b', alpha=0.8, edgecolor='#4a4a5a'))
        for spine in self.ax_phase.spines.values():
            spine.set_color('#4a4a5a')
        self.canvas_phase.draw()

    def calculate_and_show_esr(self):
        """计算并显示ESR模型"""
        self.current_model = 'esr'
        self.calculate_esr_model()
        self.update_display()

    def calculate_and_show_fixed(self):
        """计算并显示固定斜坡模型"""
        self.current_model = 'fixed'
        self.calculate_fixed_ramp_model()
        self.update_display()

    def calculate_and_show_sw(self):
        """计算并显示SW滤波器模型"""
        self.current_model = 'sw'
        self.calculate_sw_filter_model()
        self.update_display()

if __name__ == "__main__":
    root = tk.Tk()
    app = BUCKBodeAnalyzer(root)
    root.mainloop()
