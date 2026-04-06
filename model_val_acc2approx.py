import tkinter as tk
from tkinter import messagebox

import numpy as np

from model_val_acc2imp import BUCKBodeAnalyzer
from rbcot_formulas import (
    esr_approx_bode,
    interp_zero_db_crossing,
)


class BUCKApproxBodeAnalyzer(BUCKBodeAnalyzer):
    def __init__(self, root):
        super().__init__(root)
        self.root.title(
            "BUCK Bode Analyzer - RBCOT Approximate ESR Loop Gain Model"
        )
        self.btn_esr.config(text="ESR纹波补偿模型")
        self.esr_fc.config(text="Fc = -- Hz (approx)")
        self.esr_pm.config(text="PM = -- deg")
        self.calculate_esr_model()
        self.update_display()

    def calculate_esr_model(self):
        """Use the paper's simplified ESR loop gain model for Bode plotting."""
        values = self.get_input_values()
        if values is None:
            return

        try:
            Ts = values["Ts"]
            Ton = values["Ton"]
            L = values["L"]
            Cout = values["Cout"]
            Rc = values["Rc"]
            Rload = values["Rload"]

            w = 2 * np.pi * np.logspace(
                np.log10(values["freq_min"] / 10.0),
                np.log10(values["freq_max"] * 10.0),
                3000,
            )
            freq, mag, phase_deg_unwrapped, meta = esr_approx_bode(
                L, Cout, Rc, Rload, Ts, Ton, w
            )
            fz_esr = meta["omega_z"] / (2 * np.pi)

            self.esr_freq = freq
            self.esr_mag = mag
            self.esr_phase = phase_deg_unwrapped

            mag_db = 20 * np.log10(np.maximum(mag, 1e-10))
            fc_actual, cross_idx = interp_zero_db_crossing(freq, mag_db)

            fc_approx = meta["omega_c_approx"] / (2 * np.pi)
            pm_approx = meta["PM_approx_rad"] * 180.0 / np.pi
            pm_approx = min(max(pm_approx, 0.0), 180.0)

            if fc_actual is not None and cross_idx is not None:
                i0 = cross_idx
                df = freq[i0 + 1] - freq[i0]
                t = (fc_actual - freq[i0]) / df if abs(df) > 1e-30 else 0.0
                phase_at_fc = phase_deg_unwrapped[i0] + t * (
                    phase_deg_unwrapped[i0 + 1] - phase_deg_unwrapped[i0]
                )
                pm_actual = phase_at_fc + 180.0
                while pm_actual > 180:
                    pm_actual -= 360
                while pm_actual < 0:
                    pm_actual += 360
                pm_actual = min(pm_actual, 180)

                self.esr_fc.config(
                    text=f"Fc = {fc_actual/1000:.1f} kHz | Eq.(13)={fc_approx/1000:.1f} kHz"
                )
                self.esr_pm.config(
                    text=f"PM = {pm_actual:.1f} deg | Eq.(14)={pm_approx:.1f} deg | fz={fz_esr/1000:.1f} kHz"
                )
            else:
                peak_idx = int(np.argmax(mag_db))
                peak_gain_db = mag_db[peak_idx]
                peak_freq = freq[peak_idx]
                self.esr_fc.config(text="Fc = No 0 dB crossing")
                self.esr_pm.config(
                    text=(
                        f"Peak = {peak_gain_db:.1f} dB @ {peak_freq/1000:.1f} kHz"
                        f" | Eq.(13)={fc_approx/1000:.1f} kHz | Eq.(14)={pm_approx:.1f} deg"
                    )
                )

            self.sync_parameters_to_manager()

        except Exception as e:
            messagebox.showerror(
                "Calculation Error",
                f"Approximate ESR model calculation failed:\n{str(e)}",
            )
            print(f"Approximate ESR model error: {e}")


if __name__ == "__main__":
    root = tk.Tk()
    app = BUCKApproxBodeAnalyzer(root)
    root.mainloop()
