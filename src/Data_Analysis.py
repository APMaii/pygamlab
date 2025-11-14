import math
import statistics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from matplotlib.ticker import AutoMinorLocator
import seaborn as sns 
from scipy import stats
from scipy.integrate import solve_bvp
from scipy.signal import savgol_filter, find_peaks , peak_widths
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from scipy.integrate import simpson  # Corrected import




def DSC(data, application="plot", prominence=0.5, distance=5, sample_mass=1.0, heating_rate=1., orientation=None):
    """
    Perform Differential Scanning Calorimetry (DSC) data processing, analysis, and visualization.

    This function allows for automated DSC curve plotting, peak detection, 
    transition temperature determination (Tg, Tm, Tc), enthalpy (ΔH) estimation, 
    and kinetic analysis from experimental DSC datasets. The analysis can be 
    adapted for both exothermic-up and endothermic-up instrument conventions.

    Parameters
    ----------
    data : pandas.DataFrame
        Input DataFrame containing DSC measurement data. It must include one of:
        - Columns ["t", "Value"] for time-based measurements, or
        - Columns ["Temperature", "Value"] for temperature-based measurements.
    
    application : str, optional, default="plot"
        The type of analysis or operation to perform. 
        Supported options include:
        - `"plot"` : Plot the raw DSC curve.
        - `"peak_detection"` : Detect and label endothermic and exothermic peaks.
        - `"Tg"` : Estimate the glass transition temperature (Tg).
        - `"Tm"` : Determine the melting temperature (Tm).
        - `"Tc"` : Determine the crystallization temperature (Tc).
        - `"dH"` : Compute enthalpy changes (ΔH) for detected events.
        - `"kinetics"` : Estimate reaction onset, peak, endset, and corresponding ΔH.

    prominence : float, optional, default=0.5
        Minimum prominence of peaks for detection. Higher values filter out smaller peaks.
        Passed to `scipy.signal.find_peaks`.

    distance : int, optional, default=5
        Minimum number of data points between detected peaks. 
        Helps to separate closely spaced transitions.

    sample_mass : float, optional, default=1.0
        Sample mass in milligrams (mg). Used to normalize enthalpy (ΔH) values.

    heating_rate : float, optional, default=1.0
        Heating or cooling rate in °C/min. Used to normalize ΔH for temperature-based data.

    orientation : str or None, optional, default=None
        Defines the thermal orientation of the DSC instrument:
        - `"exo_up"` : Exothermic events produce positive peaks.
        - `"endo_up"` : Endothermic events produce positive peaks.
        If None, the user is prompted interactively to choose.

    Returns
    -------
    varies depending on `application`
        - `"plot"` : None
        - `"peak_detection"` : dict
            Contains coordinates of detected endothermic and exothermic peaks:
                {
                    "endothermic": [(x1, y1), (x2, y2), ...],
                    "exothermic": [(x1, y1), (x2, y2), ...]
                }
        - `"Tg"`, `"Tm"`, `"Tc"` : float
            The estimated transition temperature value in the same units as the x-axis.
        - `"dH"` : list of tuples
            Each tuple contains (Temperature, Signal, ΔH) for detected events.
        - `"kinetics"` : list of dict
            Each dictionary contains:
                {
                    "Onset": float,
                    "Peak": float,
                    "End": float,
                    "ΔH (J/g)": float
                }

    Raises
    ------
    ValueError
        If the required data columns are missing or if `application` 
        is not one of the supported analysis modes.

    Notes
    -----
    - The function automatically handles both time-based (`t`) and temperature-based (`Temperature`) DSC data.
    - The `orientation` parameter affects sign convention in peak detection and ΔH calculation.
      For example, `exo_up` instruments produce positive exothermic peaks, 
      while `endo_up` instruments produce negative ones.
    - The area under peaks (ΔH) is numerically integrated using the trapezoidal rule.

    Examples
    --------
    >>> import pandas as pd
    >>> data = pd.read_csv("sample_dsc.csv")
    >>> DSC(data, application="plot")
    # Displays the DSC curve.

    >>> results = DSC(data, application="peak_detection", orientation="exo_up")
    >>> results["exothermic"]
    [(134.2, -0.023), (276.4, -0.018)]

    >>> Tg = DSC(data, application="Tg", orientation="exo_up")
    Estimated Glass Transition Temperature (Tg): 65.12 °C

    >>> dH_values = DSC(data, application="dH", sample_mass=5.0, heating_rate=10.0, orientation="endo_up")
    Enthalpy Changes (ΔH):
    Peak at 135.50 °C, ΔH ≈ 25.432 J/g
    """
        
    # Determine X and Y columns
    if "t" in data.columns and "Value" in data.columns:
        x_col = "t"
        y_col = "Value"
    elif "Temperature" in data.columns and "Value" in data.columns:
        x_col = "Temperature"
        y_col = "Value"
    else:
        raise ValueError("Data must contain either 't' or 'Temperature' and 'Value' columns.")
    
    x = data[x_col].values
    y = data[y_col].values
    y_plot = y.copy()  # always plot raw data
    
    # Orientation handling
    if orientation is None:
        ans = input("Is your data exo up or endo up? (type 'exo' or 'endo'): ").strip().lower()
        orientation = "exo_up" if ans.startswith("exo") else "endo_up"
    
    exo_is_negative = True if orientation == "endo_up" else False
    # If endo_up → exo = negative peaks
    # If exo_up → exo = positive peaks
    
    if application == "plot":
        plt.figure(figsize=(8,5))
        plt.plot(x, y_plot, label="DSC Curve")
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.title("DSC Curve")
        plt.grid(True)
        plt.legend()
        plt.show()
    
    elif application == "peak_detection":
        if exo_is_negative:
            peaks, _ = find_peaks(y, prominence=prominence, distance=distance)       # endo
            troughs, _ = find_peaks(-y, prominence=prominence, distance=distance)    # exo
        else:
            peaks, _ = find_peaks(-y, prominence=prominence, distance=distance)      # endo
            troughs, _ = find_peaks(y, prominence=prominence, distance=distance)     # exo
        
        plt.figure(figsize=(8,5))
        plt.plot(x, y_plot, label="DSC Curve")
        plt.plot(x[peaks], y_plot[peaks], "ro", label="Endothermic Peaks")
        plt.plot(x[troughs], y_plot[troughs], "bo", label="Exothermic Peaks")
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.title("DSC Curve with Peak Detection")
        plt.grid(True)
        plt.legend()
        plt.show()
        
        return {
            "endothermic": [(x[i], y[i]) for i in peaks],
            "exothermic": [(x[i], y[i]) for i in troughs]
        }
    
    elif application == "Tg":
        if x_col == "t":
            print("Tg calculation requires Temperature as x-axis, not Time.")
            return None
        dy_dx = np.gradient(y, x)
        Tg_index = np.argmax(np.abs(dy_dx))
        Tg_value = x[Tg_index]
        
        plt.figure(figsize=(8,5))
        plt.plot(x, y_plot, label="DSC Curve")
        plt.axvline(Tg_value, color='green', linestyle='--', label=f"Tg ≈ {Tg_value:.2f}")
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.title("DSC Curve with Tg")
        plt.grid(True)
        plt.legend()
        plt.show()
        
        print(f"Estimated Glass Transition Temperature (Tg): {Tg_value:.2f} {x_col}")
        return Tg_value
    
    elif application == "Tm":
        if x_col == "t":
            print("Tm calculation requires Temperature as x-axis, not Time.")
            return None
        
        if exo_is_negative:
            peaks, _ = find_peaks(y, prominence=prominence, distance=distance)  # endo melting
        else:
            peaks, _ = find_peaks(-y, prominence=prominence, distance=distance) # endo melting
        
        if len(peaks) == 0:
            print("No melting peak detected.")
            return None
        max_peak_idx = peaks[np.argmax(np.abs(y[peaks]))]
        Tm_value = x[max_peak_idx]
        
        plt.figure(figsize=(8,5))
        plt.plot(x, y_plot, label="DSC Curve")
        plt.plot(x[max_peak_idx], y_plot[max_peak_idx], "ro", label=f"Tm ≈ {Tm_value:.2f}")
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.title("DSC Curve with Melting Point")
        plt.grid(True)
        plt.legend()
        plt.show()
        
        print(f"Estimated Melting Point (Tm): {Tm_value:.2f} {x_col}")
        return Tm_value
    
    elif application == "Tc":
        if x_col == "t":
            print("Tc calculation requires Temperature as x-axis, not Time.")
            return None
        
        if exo_is_negative:
            troughs, _ = find_peaks(-y, prominence=prominence, distance=distance)  # exo crystallization
        else:
            troughs, _ = find_peaks(y, prominence=prominence, distance=distance)   # exo crystallization
        
        if len(troughs) == 0:
            print("No crystallization peak detected.")
            return None
        max_trough_idx = troughs[np.argmax(np.abs(y[troughs]))]
        Tc_value = x[max_trough_idx]
        
        plt.figure(figsize=(8,5))
        plt.plot(x, y_plot, label="DSC Curve")
        plt.plot(x[max_trough_idx], y_plot[max_trough_idx], "bo", label=f"Tc ≈ {Tc_value:.2f}")
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.title("DSC Curve with Crystallization Temperature")
        plt.grid(True)
        plt.legend()
        plt.show()
        
        print(f"Estimated Crystallization Temperature (Tc): {Tc_value:.2f} {x_col}")
        return Tc_value
    
    elif application == "dH":
        if exo_is_negative:
            peaks, _ = find_peaks(y, prominence=prominence, distance=distance)
            troughs, _ = find_peaks(-y, prominence=prominence, distance=distance)
        else:
            peaks, _ = find_peaks(-y, prominence=prominence, distance=distance)
            troughs, _ = find_peaks(y, prominence=prominence, distance=distance)
        
        all_events = np.sort(np.concatenate((peaks, troughs)))
        if len(all_events) == 0:
            print("No events found for enthalpy calculation.")
            return None
        
        results = []
        plt.figure(figsize=(8,5))
        plt.plot(x, y_plot, label="DSC Curve")
        
        for i in all_events:
            left = max(0, i-20)
            right = min(len(x)-1, i+20)
            dx = np.diff(x[left:right+1])
            avg_y = (y[left:right] + y[left+1:right+1]) / 2
            area = np.sum(avg_y * dx)
            
            if x_col == "t":
                dH = area / sample_mass
            else:
                dH = area / (heating_rate * sample_mass)
            
            # Apply orientation sign convention
            if orientation == "exo_up":
                dH = -dH
            
            results.append((x[i], y[i], dH))
            plt.fill_between(x[left:right+1], y_plot[left:right+1], alpha=0.3, label=f"ΔH at {x[i]:.2f}")
        
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.title("DSC Curve with Enthalpy Areas")
        plt.grid(True)
        plt.legend()
        plt.show()
        
        print("Enthalpy Changes (ΔH):")
        for r in results:
            print(f"Peak at {r[0]:.2f} {x_col}, ΔH ≈ {r[2]:.3f} J/g")
        
        return results
    
    elif application == "kinetics":
        peaks, _ = find_peaks(np.abs(y), prominence=prominence, distance=distance)
        if len(peaks) == 0:
            print("No reaction events detected.")
            return None
        
        results = []
        plt.figure(figsize=(8,5))
        plt.plot(x, y_plot, label="DSC Curve")
        
        for i in peaks:
            left = max(0, i-30)
            right = min(len(x)-1, i+30)
            onset = x[left]
            endset = x[right]
            peak_temp = x[i]
            
            dx = np.diff(x[left:right+1])
            avg_y = (y[left:right] + y[left+1:right+1]) / 2
            area = np.sum(avg_y * dx)
            
            if x_col == "t":
                dH = area / sample_mass
            else:
                dH = area / (heating_rate * sample_mass)
            
            # Orientation adjustment
            if orientation == "exo_up":
                dH = -dH
            
            results.append({
                "Onset": onset,
                "Peak": peak_temp,
                "End": endset,
                "ΔH (J/g)": dH
            })
            
            plt.fill_between(x[left:right+1], y_plot[left:right+1], alpha=0.3, label=f"Reaction at {peak_temp:.2f}")
            plt.axvline(onset, color="green", linestyle="--")
            plt.axvline(endset, color="red", linestyle="--")
            plt.plot(peak_temp, y_plot[i], "ro")
        
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.title("DSC Reaction Kinetics Analysis")
        plt.legend()
        plt.grid(True)
        plt.show()
        
        print("Reaction Kinetics Results:")
        for r in results:
            print(f"Onset: {r['Onset']:.2f} {x_col}, Peak: {r['Peak']:.2f}, End: {r['End']:.2f}, ΔH ≈ {r['ΔH (J/g)']:.3f} J/g")
        
        return results
    
    else:
        raise ValueError("Application must be one of: 'plot','peak_detection','Tg','Tm','Tc','dH','kinetics'.")



def EELS_Analysis(data, application):
    """
    Perform quantitative and visual analysis of Electron Energy Loss Spectroscopy (EELS) data.

    This function allows detailed inspection of EELS spectra across different energy-loss regions,
    including Zero-Loss Peak (ZLP), low-loss, and core-loss regions. It supports both raw plotting 
    and automated analysis such as peak detection, band gap estimation, plasmon peak identification, 
    and fine structure analysis (ELNES/EXELFS).

    Parameters
    ----------
    data : list of tuples/lists or pandas.DataFrame
        Input EELS data.  
        - If a list, each element should be `(energy_loss, intensity)`.  
        - If a DataFrame, it must contain columns:
          - `"energy_loss"` : float — Energy loss values in eV.
          - `"Intensity"`   : float — Measured intensity (arbitrary units).

    application : str
        Specifies the type of analysis to perform. Options include:
        
        - `"plot"` :
            Simply plot the EELS spectrum for visual inspection.
        
        - `"ZLP"` :
            Analyze the **Zero-Loss Peak (ZLP)** region near 0 eV.  
            Automatically detects the main elastic scattering peak and estimates:
              - Peak position (energy in eV)
              - Peak height (intensity)
              - Full Width at Half Maximum (FWHM) if determinable.  
            The results are printed and visualized with the smoothed curve and annotations.
        
        - `"low_loss"` :
            Analyze the **Low-Loss Region (−5 to 50 eV)** including pre-zero baseline.  
            Performs:
              - Baseline smoothing and visualization
              - Detection of **plasmon peaks** (typically <25 eV)
              - Estimation of **optical band gap (Eg)** via derivative onset method.  
            Prints and plots plasmon peaks and band gap position.
        
        - `"core_loss"` :
            Analyze the **Core-Loss (High-Loss) Region (>50 eV)**.  
            Performs:
              - Edge onset detection using signal derivative
              - Step height estimation at the absorption edge
              - Identification of fine structure features:
                  * ELNES (Energy-Loss Near Edge Structure) within ~30 eV above onset  
                  * EXELFS (Extended Energy-Loss Fine Structure) oscillations beyond onset  
            Results include detected edges, peaks, and oscillations with visualized spectrum.

    Returns
    -------
    None
        The function primarily displays plots and prints analysis results to the console.
        Key detected parameters (peak positions, FWHM, etc.) are reported in the output text.

    Notes
    -----
    - Smoothing is performed using a Savitzky-Golay filter (`scipy.signal.savgol_filter`)
      with a default window length of 11 and polynomial order of 3.
    - Peak detection uses `scipy.signal.find_peaks` with adaptive height thresholds.
    - Energy regions are automatically segmented as:
        * ZLP: around 0 eV
        * Low-loss: −5–50 eV
        * Core-loss: >50 eV
    - The function assumes intensity units are arbitrary and energy loss is in electronvolts (eV).

    Examples
    --------
    >>> import pandas as pd
    >>> data = pd.DataFrame({
    ...     "energy_loss": np.linspace(-10, 200, 500),
    ...     "Intensity": np.random.random(500) * np.exp(-np.linspace(-10, 200, 500)/100)
    ... })
    >>> EELS_Analysis(data, "plot")
    # Displays the EELS spectrum

    >>> EELS_Analysis(data, "ZLP")
    # Detects and plots Zero-Loss Peak with FWHM estimation

    >>> EELS_Analysis(data, "low_loss")
    # Identifies plasmon peaks and estimates band gap

    >>> EELS_Analysis(data, "core_loss")
    # Detects absorption edge and ELNES/EXELFS features
    """

    # Handle both DataFrame and list input
    if isinstance(data, pd.DataFrame):
        energy_loss = data["energy_loss"].values
        intensity = data["Intensity"].values
    else:
        energy_loss = np.array([point[0] for point in data])
        intensity = np.array([point[1] for point in data])

    # ====== PLOT APPLICATION ======
    if application == "plot":
        plt.figure(figsize=(8, 5))
        plt.plot(energy_loss, intensity, color="blue", linewidth=1.5)
        plt.title("EELS Spectrum", fontsize=14)
        plt.xlabel("Energy Loss (ΔE) [eV]", fontsize=12)
        plt.ylabel("Intensity (a.u.)", fontsize=12)
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()

    # ====== ZERO LOSS PEAK ======
    elif application == "ZLP":
        smoothed_intensity = savgol_filter(intensity, window_length=11, polyorder=3)
        peaks, _ = find_peaks(smoothed_intensity, height=np.max(smoothed_intensity)*0.5)
        if len(peaks) == 0:
            print("⚠️ No Zero-Loss Peak detected.")
            return
        
        peak_idx = peaks[np.argmin(np.abs(energy_loss[peaks]))]
        peak_position = energy_loss[peak_idx]
        peak_height = smoothed_intensity[peak_idx]

        half_max = peak_height / 2
        left_idx = np.where(smoothed_intensity[:peak_idx] <= half_max)[0]
        right_idx = np.where(smoothed_intensity[peak_idx:] <= half_max)[0]
        if len(left_idx) > 0 and len(right_idx) > 0:
            fwhm = energy_loss[peak_idx + right_idx[0]] - energy_loss[left_idx[-1]]
        else:
            fwhm = None

        print("📊 Zero-Loss Peak (ZLP) Analysis:")
        print(f"   Peak Position: {peak_position:.3f} eV")
        print(f"   Peak Height  : {peak_height:.3f} a.u.")
        if fwhm:
            print(f"   Peak Width (FWHM): {fwhm:.3f} eV")
        else:
            print("   Peak Width (FWHM): Could not be determined")

        plt.figure(figsize=(8, 5))
        plt.plot(energy_loss, intensity, label="Raw Data", color="lightgray")
        plt.plot(energy_loss, smoothed_intensity, label="Smoothed", color="blue")
        plt.axvline(peak_position, color="red", linestyle="--", label="ZLP Peak")
        if fwhm:
            plt.axhline(half_max, color="green", linestyle="--", label="Half Maximum")
        plt.title("Zero-Loss Peak (ZLP) Analysis", fontsize=14)
        plt.xlabel("Energy Loss (ΔE) [eV]", fontsize=12)
        plt.ylabel("Intensity (a.u.)", fontsize=12)
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()

    # ====== LOW LOSS REGION ======
    elif application == "low_loss":
    # Keep a window that includes negative values for proper baseline
        mask_full = (energy_loss >= -5) & (energy_loss <= 50)   # allow pre-zero baseline
        E_full = energy_loss[mask_full]
        I_full = intensity[mask_full]

        if len(E_full) == 0:
            print("⚠️ No data in -5–50 eV range.")
            return

        smoothed_full = savgol_filter(I_full, window_length=11, polyorder=3)

        # Extract 0–50 eV part (for plasmon & bandgap plotting)
        mask_pos = (E_full >= 0) & (E_full <= 50)
        E = E_full[mask_pos]
        smoothed_I = smoothed_full[mask_pos]
        I = I_full[mask_pos]

        # ---- Detect plasmon peaks (<20 eV typically) ----
        plasmon_mask = (E >= 5) & (E <= 25)   # restrict to plasmon range
        E_plasmon = E[plasmon_mask]
        I_plasmon = smoothed_I[plasmon_mask]

        peaks, props = find_peaks(I_plasmon, 
                                height=np.max(I_plasmon)*0.08,  # lower threshold
                                distance=5)

        plasmon_info = []
        for p in peaks:
            peak_pos = E_plasmon[p]
            peak_height = I_plasmon[p]
            half_max = peak_height / 2
            left_idx = np.where(I_plasmon[:p] <= half_max)[0]
            right_idx = np.where(I_plasmon[p:] <= half_max)[0]
            if len(left_idx) > 0 and len(right_idx) > 0:
                width = E_plasmon[p + right_idx[0]] - E_plasmon[left_idx[-1]]
            else:
                width = None
            plasmon_info.append((peak_pos, peak_height, width))

        # ---- Band gap estimation ----
        # baseline = min around -5 to 0 eV
        baseline_region = (E_full >= -5) & (E_full <= 0)
        baseline = np.mean(smoothed_full[baseline_region])
        dI = np.gradient(smoothed_I, E)   # derivative
        dI_threshold = np.max(dI) * 0.1  

        valid_idx = np.where(E > 2)[0]
        onset_candidates = valid_idx[dI[valid_idx] > dI_threshold]

        if len(onset_candidates) > 0:
            band_gap = E[onset_candidates[0]]
        else:
            band_gap = None

        # ---- Print results ----
        print("📊 Low-Loss Region Analysis (-5–50 eV):")
        if band_gap:
            print(f"   Estimated Band Gap Eg: {band_gap:.3f} eV")
        else:
            print("   Band Gap: Could not be determined")

        if plasmon_info:
            for i, (pos, h, w) in enumerate(plasmon_info, 1):
                print(f"   Plasmon Peak {i}: {pos:.3f} eV, Height={h:.3f}, Width={w if w else 'N/A'} eV")
        else:
            print("   No plasmon peaks detected.")

        # ---- Plot ----
        # ---- Plot ----
        plt.figure(figsize=(8, 5))
        plt.plot(E_full, I_full, label="Raw Data (-5–50 eV)", color="lightgray")
        plt.plot(E_full, smoothed_full, label="Smoothed", color="blue")
        if band_gap:
            plt.axvline(band_gap, color="green", linestyle="--", label=f"Band Gap ≈ {band_gap:.2f} eV")
        for i, (pos, h, w) in enumerate(plasmon_info, 1):
            plt.axvline(pos, color="red", linestyle="--", label=f"Plasmon {i}: {pos:.2f} eV")
        plt.title("Low-Loss Region Analysis (with pre-zero baseline)", fontsize=14)
        plt.xlabel("Energy Loss (ΔE) [eV]", fontsize=12)
        plt.ylabel("Intensity (a.u.)", fontsize=12)
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()

    # ====== CORE LOSS REGION ======
    elif application == "core_loss":
        mask = (energy_loss > 50)
        E = energy_loss[mask]
        I = intensity[mask]
        if len(E) == 0:
            print("⚠️ No data in >50 eV range.")
            return

        smoothed_I = savgol_filter(I, window_length=11, polyorder=3)

        # Detect edges via derivative
        derivative = np.gradient(smoothed_I, E)
        edge_idx = np.where(derivative > np.max(derivative)*0.5)[0]
        if len(edge_idx) > 0:
            edge_onset = E[edge_idx[0]]
            step_height = smoothed_I[edge_idx[-1]] - smoothed_I[edge_idx[0]]
        else:
            edge_onset, step_height = None, None

        # ELNES peaks: local maxima near edge (within 30 eV above onset)
        elnes_info = []
        if edge_onset:
            near_edge_mask = (E >= edge_onset) & (E <= edge_onset + 30)
            E_near, I_near = E[near_edge_mask], smoothed_I[near_edge_mask]
            peaks, _ = find_peaks(I_near, height=np.max(I_near)*0.2)
            for p in peaks:
                elnes_info.append((E_near[p], I_near[p]))

        # EXELFS: oscillations further above the edge (>30 eV after onset)
        exelfs_info = []
        if edge_onset:
            exelfs_mask = (E > edge_onset + 30)
            E_far, I_far = E[exelfs_mask], smoothed_I[exelfs_mask]
            peaks, _ = find_peaks(I_far, height=np.mean(I_far))
            for p in peaks:
                exelfs_info.append((E_far[p], I_far[p]))

        # Print results
        print("📊 Core-Loss Region Analysis (>50 eV):")
        if edge_onset:
            print(f"   Edge Onset Energy: {edge_onset:.3f} eV")
            print(f"   Step Height: {step_height:.3f} a.u.")
        else:
            print("   No clear edge detected.")
        if elnes_info:
            for i, (pos, h) in enumerate(elnes_info, 1):
                print(f"   ELNES Peak {i}: {pos:.3f} eV, Intensity={h:.3f}")
        if exelfs_info:
            print(f"   Detected {len(exelfs_info)} EXELFS oscillation peaks.")

        # Plot
        plt.figure(figsize=(8, 5))
        plt.plot(E, I, label="Raw Data", color="lightgray")
        plt.plot(E, smoothed_I, label="Smoothed", color="blue")
        if edge_onset:
            plt.axvline(edge_onset, color="red", linestyle="--", label=f"Edge Onset ≈ {edge_onset:.2f} eV")
        for i, (pos, h) in enumerate(elnes_info, 1):
            plt.plot(pos, h, "go", label=f"ELNES {i}" if i==1 else "")
        for i, (pos, h) in enumerate(exelfs_info, 1):
            plt.plot(pos, h, "mo", label="EXELFS" if i==1 else "")
        plt.title("Core-Loss Region Analysis (>50 eV)", fontsize=14)
        plt.xlabel("Energy Loss (ΔE) [eV]", fontsize=12)
        plt.ylabel("Intensity (a.u.)", fontsize=12)
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()




def Raman_Analysis(data, application):
    """
    Perform quantitative and visual analysis of Raman spectroscopy data.

    This function provides flexible tools for visualizing and analyzing
    Raman spectra. It supports basic spectrum plotting and automated
    peak detection for identifying characteristic Raman bands.

    Parameters
    ----------
    data : list of tuples or list of lists
        Raman spectrum data, where each element corresponds to one measurement point:
            (wavenumber, intensity)
        - wavenumber : float
            Raman shift in inverse centimeters (cm⁻¹)
        - intensity : float
            Measured Raman intensity in arbitrary units (a.u.)

        Example:
        >>> data = [(100, 0.1), (150, 0.5), (200, 1.2)]

    application : str
        Defines the type of analysis to perform. Supported options:

        - `"plot"` :
            Plot the Raman spectrum with labeled axes and gridlines for quick visual inspection.
        
        - `"peak_detect"` :
            Automatically detect and highlight prominent peaks in the Raman spectrum.  
            Peak detection is performed using `scipy.signal.find_peaks` with:
              * Minimum peak height = 10% of maximum intensity
              * Minimum distance between peaks = 5 data points  
            The detected peaks are printed (wavenumber and intensity) and plotted with red markers.

    Raises
    ------
    ValueError
        If the data format is invalid or the specified application is not supported.

    Returns
    -------
    None
        The function generates plots and prints peak data to the console when applicable.
        No explicit return value.

    Notes
    -----
    - The function assumes Raman shift values are given in cm⁻¹ and intensity in arbitrary units.
    - The x-axis is plotted as Raman shift (increasing rightward). Uncomment the `invert_xaxis()`
      line to follow the traditional Raman plotting convention (decreasing Raman shift).
    - Peak detection parameters (height and distance) can be fine-tuned based on spectral resolution.

    Examples
    --------
    >>> import numpy as np
    >>> # Generate synthetic Raman data
    >>> wavenumbers = np.linspace(100, 2000, 500)
    >>> intensities = np.exp(-((wavenumbers - 1350)/40)**2) + 0.5*np.exp(-((wavenumbers - 1580)/30)**2)
    >>> data = list(zip(wavenumbers, intensities))

    >>> Raman_Analysis(data, "plot")
    # Displays the Raman spectrum

    >>> Raman_Analysis(data, "peak_detect")
    # Detects and highlights Raman peaks in the spectrum
    """
    
    # --- Convert data into two lists ---
    try:
        wavenumbers = [point[0] for point in data]
        intensities = [point[1] for point in data]
    except Exception as e:
        raise ValueError("Data format must be list of (wavenumber, intensity) pairs.") from e
    
    wavenumbers = np.array(wavenumbers)
    intensities = np.array(intensities)
    
    # --- Application handling ---
    if application == "plot":
        plt.figure(figsize=(8,5))
        plt.plot(wavenumbers, intensities, color="blue", linewidth=1.5)
        plt.xlabel("Raman Shift (cm⁻¹)")
        plt.ylabel("Intensity (a.u.)")
        plt.title("Raman Spectrum")
        #plt.gca().invert_xaxis()  # Raman convention
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()
    
    elif application == "peak_detect":
        # Detect peaks (threshold: 10% of max, min distance=5 points)
        peaks, properties = find_peaks(intensities, height=np.max(intensities)*0.1, distance=5)
        
        plt.figure(figsize=(8,5))
        plt.plot(wavenumbers, intensities, color="blue", linewidth=1.5)
        plt.plot(wavenumbers[peaks], intensities[peaks], "ro", label="Detected Peaks")
        plt.xlabel("Raman Shift (cm⁻¹)")
        plt.ylabel("Intensity (a.u.)")
        plt.title("Raman Spectrum with Peak Detection")
        #plt.gca().invert_xaxis()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.show()
        
        # Print detected peaks
        print("Detected peak positions (cm⁻¹):")
        for wn, inten in zip(wavenumbers[peaks], intensities[peaks]):
            print(f"Peak at {wn:.2f} cm⁻¹, Intensity = {inten:.4f}")
    
    
    else:
        raise ValueError(f"Application '{application}' not implemented yet.")




def TGA(data, application):
    """
    Perform multi-mode Thermogravimetric Analysis (TGA) for material characterization.

    This function enables comprehensive TGA data analysis for studying
    thermal stability, composition, surface modification, and reaction kinetics.
    It supports visualization, derivative thermogravimetry (DTG), decomposition
    step identification, and moisture or solvent content determination.

    Parameters
    ----------
    data : pandas.DataFrame
        Experimental TGA dataset with the following required columns:
        - 'Temp' : float
            Temperature in degrees Celsius (°C)
        - 'Mass' : float
            Corresponding sample mass in percentage (%)

        Example:
        >>> data = pd.DataFrame({
        ...     "Temp": [25, 100, 200, 300],
        ...     "Mass": [100, 99.5, 80.2, 10.5]
        ... })

    application : str
        Defines the type of analysis to perform. Supported options include:

        - `"plot"` :
            Plot the raw TGA curve (Mass vs. Temperature).

        - `"peaks"` :
            Compute and display the derivative thermogravimetry (DTG) curve
            and identify key decomposition peaks using `scipy.signal.find_peaks`.

        - `"stability"` :
            Estimate the onset temperature of thermal degradation by
            tangent extrapolation from the baseline region.

        - `"moisture"` :
            Calculate moisture or solvent content based on mass loss before
            the first decomposition event (typically below 150 °C).

        - `"functionalization"` :
            Identify surface functionalization or modification steps by
            detecting multiple degradation peaks above 150 °C.

        - `"composition"` :
            Estimate polymer and filler content from the initial and final
            mass values (residue analysis).

        - `"DTG"` :
            Compute and plot the first derivative of the TGA curve (dM/dT)
            for insight into reaction rate behavior.

        - `"decomposition_steps"` :
            Identify and quantify major decomposition events (DTG peaks),
            returning their temperatures and mass values.

        - `"kinetics"` :
            Evaluate relative reaction rates and identify the fastest
            decomposition step (maximum |dM/dT| above 150 °C).

    Raises
    ------
    TypeError
        If the input is not a pandas DataFrame.
    ValueError
        If required columns ('Temp', 'Mass') are missing or
        if the specified application is not supported.

    Returns
    -------
    object
        Depends on the application:

        - `"plot"` :
            Displays the TGA curve; returns `None`.
        - `"peaks"` :
            DataFrame containing detected DTG peak temperatures and intensities.
        - `"stability"` :
            Dictionary with onset temperature and mass at onset.
        - `"moisture"` :
            Dictionary with moisture content, cutoff temperature, and mass loss.
        - `"functionalization"` :
            DataFrame listing detected modification steps.
        - `"composition"` :
            Dictionary with polymer and filler content percentages.
        - `"DTG"` :
            DataFrame of temperatures and corresponding dM/dT values.
        - `"decomposition_steps"` :
            DataFrame of decomposition step information.
        - `"kinetics"` :
            Dictionary with step-wise reaction rate data and the fastest decomposition step.

    Notes
    -----
    - TGA data should be preprocessed to ensure monotonic temperature increase.
    - The function uses numerical differentiation (`np.gradient`) for DTG calculations.
    - Peak prominence thresholds can be adjusted to improve detection sensitivity.
    - Onset temperatures are approximate and depend on the slope estimation method.

    Examples
    --------
    >>> import pandas as pd, numpy as np
    >>> T = np.linspace(25, 800, 300)
    >>> M = 100 - 0.05*(T - 25) + 10*np.exp(-((T-400)/50)**2)
    >>> data = pd.DataFrame({"Temp": T, "Mass": M})

    >>> TGA(data, "plot")
    # Displays the TGA curve.

    >>> peaks_info = TGA(data, "peaks")
    >>> print(peaks_info.head())

    >>> stability = TGA(data, "stability")
    >>> print(stability)
    """

    if not isinstance(data, pd.DataFrame):
        raise TypeError("Input data must be a pandas DataFrame.")

    if "Temp" not in data.columns or "Mass" not in data.columns:
        raise ValueError("DataFrame must contain 'Temp' and 'Mass' columns.")

    T = data["Temp"].values
    M = data["Mass"].values

    # --- Plot Application ---
    if application == "plot":
        plt.figure(figsize=(8, 5))
        plt.plot(T, M, label="TGA Curve", color="blue")
        plt.xlabel("Temperature (°C)")
        plt.ylabel("Mass (%)")
        plt.title("Thermogravimetric Analysis (TGA)")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.show()
        return None

    # --- Peak Detection Application ---
    elif application == "peaks":
        dM_dT = np.gradient(M, T)
        dM_dT_pos = -dM_dT
        peaks, _ = find_peaks(dM_dT_pos, prominence=0.01)

        peaks_info = pd.DataFrame({
            "Temperature": T[peaks],
            "dM/dT": dM_dT_pos[peaks],
            "Mass": M[peaks]
        })

        plt.figure(figsize=(10, 6))
        plt.subplot(2, 1, 1)
        plt.plot(T, M, label="TGA Curve", color="blue")
        plt.xlabel("Temperature (°C)")
        plt.ylabel("Mass (%)")
        plt.title("TGA Curve")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()

        plt.subplot(2, 1, 2)
        plt.plot(T, dM_dT_pos, label="DTG Curve (dM/dT)", color="red")
        plt.plot(T[peaks], dM_dT_pos[peaks], "ko", label="Peaks")
        for p in peaks:
            plt.text(T[p], dM_dT_pos[p], f"T={T[p]:.1f}", fontsize=8, ha="center", va="bottom")
        plt.xlabel("Temperature (°C)")
        plt.ylabel("dM/dT")
        plt.title("DTG Curve with Peak Detection")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.tight_layout()
        plt.show()

        return peaks_info

    # --- Stability Application ---
    elif application == "stability":
        dM_dT = np.gradient(M, T)

        # Only consider data above 150 °C
        mask = T >= 150
        T_sub = T[mask]
        M_sub = M[mask]
        dM_dT_sub = dM_dT[mask]

        onset_temp = None
        onset_mass = None
        tangent_line_T = None
        tangent_line_M = None

        if len(T_sub) > 0:
            # Find the index of the maximum negative derivative (steepest drop)
            min_deriv_idx = np.argmin(dM_dT_sub)
            
            # Get the corresponding temperature and mass
            max_degrad_temp = T_sub[min_deriv_idx]
            max_degrad_mass = M_sub[min_deriv_idx]
            
            # Get the slope at this point (the maximum negative derivative)
            slope = dM_dT_sub[min_deriv_idx]

            # Define the initial mass baseline for extrapolation
            # This is often the average mass in a stable region before degradation
            # We'll use the average of the first few points after 150 C
            baseline_mass = np.mean(M_sub[:10])

            # Calculate the onset temperature using the tangent line equation:
            # M_tangent = slope * T_tangent + intercept
            # intercept = max_degrad_mass - slope * max_degrad_temp
            # At onset, M_tangent = baseline_mass
            # baseline_mass = slope * onset_temp + intercept
            # onset_temp = (baseline_mass - intercept) / slope
            onset_temp = max_degrad_temp - (max_degrad_mass - baseline_mass) / slope
            
            # Now find the corresponding mass on the original curve
            # We need to find the index of the temperature closest to onset_temp
            onset_idx = np.abs(T - onset_temp).argmin()
            onset_mass = M[onset_idx]

            # Generate points for plotting the tangent line
            tangent_line_T = np.array([onset_temp, max_degrad_temp])
            tangent_line_M = slope * (tangent_line_T - max_degrad_temp) + max_degrad_mass

        # --- Plot ---
        plt.figure(figsize=(10, 6))
        plt.plot(T, M, label="TGA Curve", color="blue")
        
        if onset_temp is not None:
            plt.axvline(onset_temp, color="red", linestyle="--", label=f"Onset ~ {onset_temp:.1f} °C")
            plt.scatter(onset_temp, onset_mass, color="red", zorder=5) # zorder ensures it's on top
            
            # Plot the tangent line used for calculation
            plt.plot(tangent_line_T, tangent_line_M, color="green", linestyle="-", label="Tangent Line")
            
        plt.xlabel("Temperature (°C)")
        plt.ylabel("Mass (%)")
        plt.title("Thermal Stability (Onset of Degradation)")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.show()

        return {"Onset_Temperature": onset_temp, "Mass_at_Onset": onset_mass}



    # --- Moisture Content Application ---
    elif application == "moisture":
        dM_dT = np.gradient(M, T)

        # Dynamic cutoff detection (first significant drop)
        threshold = -0.01 * np.max(np.abs(dM_dT))
        idx = np.where(dM_dT < threshold)[0]

        if len(idx) > 0:
            cutoff_temp = T[idx[0]]
        else:
            cutoff_temp = 150  # fallback cutoff

        # Mass values at start and at cutoff
        initial_mass = M[0]
        final_mass_cutoff = M[T <= cutoff_temp][-1]

        # Moisture/solvent loss up to cutoff
        moisture_loss_cutoff = initial_mass - final_mass_cutoff
        moisture_percent_cutoff = (moisture_loss_cutoff / initial_mass) * 100

        # --- Additional calculation: Mass loss up to 150 °C (fixed moisture definition) ---
        if np.any(T <= 150):
            final_mass_150 = M[T <= 150][-1]
            moisture_loss_150 = initial_mass - final_mass_150
            moisture_percent_150 = (moisture_loss_150 / initial_mass) * 100
        else:
            final_mass_150 = None
            moisture_loss_150 = None
            moisture_percent_150 = None

        # Plot
        plt.figure(figsize=(8, 5))
        plt.plot(T, M, label="TGA Curve", color="blue")
        plt.axvline(cutoff_temp, color="green", linestyle="--", label=f"Cutoff ~ {cutoff_temp:.1f} °C")
        plt.scatter([T[0], cutoff_temp], [initial_mass, final_mass_cutoff], color="red", label="Dynamic Moisture Loss")
        plt.axvline(150, color="orange", linestyle="--", label="150 °C Reference")
        if final_mass_150 is not None:
            plt.scatter([150], [final_mass_150], color="purple", label="Mass @ 150 °C")
        plt.xlabel("Temperature (°C)")
        plt.ylabel("Mass (%)")
        plt.title("Moisture/Solvent Content Analysis")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.show()

        return {
            "Initial_Mass": initial_mass,
            "Final_Mass_at_cutoff": final_mass_cutoff,
            "Moisture_Loss_Cutoff(%)": moisture_percent_cutoff,
            "Cutoff_Temperature": cutoff_temp,
            "Final_Mass_at_150C": final_mass_150,
            "Moisture_Loss_150C(%)": moisture_percent_150
        }


    # --- Functionalization / Surface Modification Application ---
    elif application == "functionalization":
        dM_dT = np.gradient(M, T)
        peaks, _ = find_peaks(-dM_dT, prominence=0.01)

        # Filter out peaks below 150 °C (moisture/solvent)
        peaks = [p for p in peaks if T[p] >= 150]

        results = []
        for i, p in enumerate(peaks):
            if i < len(peaks) - 1:
                mass_loss = M[p] - M[peaks[i+1]]
            else:
                mass_loss = M[p] - M[-1]
            results.append({
                "Step": i+1,
                "Temperature": T[p],
                "Mass_at_Peak": M[p],
                "Estimated_Mass_Loss": abs(mass_loss)
            })

        results_df = pd.DataFrame(results)

        plt.figure(figsize=(8, 5))
        plt.plot(T, M, label="TGA Curve", color="blue")
        plt.plot(T[peaks], M[peaks], "ro", label="Functionalization Steps")
        for r in results:
            plt.text(r["Temperature"], r["Mass_at_Peak"],
                    f"Step {r['Step']}\n{r['Temperature']:.1f} °C",
                    fontsize=8, ha="center", va="bottom")
        plt.xlabel("Temperature (°C)")
        plt.ylabel("Mass (%)")
        plt.title("Functionalization / Surface Modification Detection")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.show()

        return results_df


    # --- Composition Application ---
    elif application == "composition":
        initial_mass = M[0]
        final_mass = M[-1]  # residue at high T
        polymer_loss = initial_mass - final_mass

        polymer_percent = (polymer_loss / initial_mass) * 100
        filler_percent = (final_mass / initial_mass) * 100

        plt.figure(figsize=(8, 5))
        plt.plot(T, M, label="TGA Curve", color="blue")
        plt.axhline(final_mass, color="orange", linestyle="--",
                    label=f"Residue ~ {final_mass:.1f}%")
        plt.xlabel("Temperature (°C)")
        plt.ylabel("Mass (%)")
        plt.title("Composition Estimation (Polymer vs Filler)")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.show()

        return {
            "Initial_Mass": initial_mass,
            "Final_Residue": final_mass,
            "Polymer_Content(%)": polymer_percent,
            "Filler_Content(%)": filler_percent
        }
    elif application == "DTG":
        dM_dT = np.gradient(M, T)

        plt.figure(figsize=(8, 5))
        plt.plot(T, dM_dT, label="DTG Curve (dM/dT)", color="red")
        plt.xlabel("Temperature (°C)")
        plt.ylabel("dM/dT (%/°C)")
        plt.title("Derivative Thermogravimetry (DTG)")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.show()

        return pd.DataFrame({"Temperature": T, "dM/dT": dM_dT})
    
    elif application == "decomposition_steps":
        dM_dT = np.gradient(M, T)
    
    # Detect peaks in DTG (negative peaks for weight loss)
        peaks, properties = find_peaks(-dM_dT, prominence=0.01)
        
        # Prepare DataFrame with decomposition info
        decomposition_info = pd.DataFrame({
            "Step": np.arange(1, len(peaks)+1),
            "Temperature": T[peaks],
            "Mass_at_Peak": M[peaks],
            "dM/dT": dM_dT[peaks]
        })

        # Plot TGA with DTG peaks highlighted
        plt.figure(figsize=(8,5))
        plt.plot(T, M, label="TGA Curve", color="blue")
        plt.plot(T[peaks], M[peaks], "ro", label="Decomposition Steps (DTG Peaks)")
        for r in decomposition_info.itertuples():
            plt.text(r.Temperature, r.Mass_at_Peak,
                    f"Step {r.Step}\n{r.Temperature:.1f} °C",
                    fontsize=8, ha="center", va="bottom")
        plt.xlabel("Temperature (°C)")
        plt.ylabel("Mass (%)")
        plt.title("Decomposition Steps Identified via DTG")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.show()

        return decomposition_info

    elif application == "kinetics":
        dM_dT = np.gradient(M, T)

        # Detect peaks in DTG
        peaks, properties = find_peaks(-dM_dT, prominence=0.01)

        # Keep only peaks above 150 °C
        peaks = [p for p in peaks if T[p] >= 150]

        # Quantify reaction rates (absolute value of dM/dT at peaks)
        kinetics_info = pd.DataFrame({
            "Step": np.arange(1, len(peaks)+1),
            "Temperature": [T[p] for p in peaks],
            "Mass_at_Peak": [M[p] for p in peaks],
            "Reaction_Rate_dM_dT": [-dM_dT[p] for p in peaks]  # positive values
        })

        # Identify fastest decomposition (highest reaction rate)
        if not kinetics_info.empty:
            fastest_idx = kinetics_info["Reaction_Rate_dM_dT"].idxmax()
            fastest_step = kinetics_info.loc[fastest_idx].to_dict()
        else:
            fastest_step = None

        # Plot DTG with reaction rates
        plt.figure(figsize=(8,5))
        plt.plot(T, dM_dT, label="DTG Curve (dM/dT)", color="red")
        plt.plot(T[peaks], dM_dT[peaks], "ko", label="Reaction Peaks")

        for r in kinetics_info.itertuples():
            plt.text(r.Temperature, dM_dT[T == r.Temperature][0], 
                    f"Step {r.Step}\nRate={r.Reaction_Rate_dM_dT:.3f}",
                    fontsize=8, ha="center", va="bottom")

        if fastest_step is not None:
            plt.axvline(fastest_step["Temperature"], color="blue", linestyle="--",
                        label=f"Fastest Decomp: Step {fastest_step['Step']}")

        plt.xlabel("Temperature (°C)")
        plt.ylabel("dM/dT (%/°C)")
        plt.title("Reaction Rates / Thermal Degradation Kinetics (Above 150 °C)")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.show()

        return {
            "Kinetics_Info": kinetics_info,
            "Fastest_Decomposition": fastest_step
        }


    # --- Catalytic Effect Application ---
    '''elif application == "catalysis":
        """
        Studying catalytic effects of nanoparticles on degradation.
        Requires DTG analysis to identify accelerated degradation steps.
        """
        # Compute DTG
        dM_dT = np.gradient(M, T)

        # Detect DTG peaks
        peaks, properties = find_peaks(-dM_dT, prominence=0.01)

        # Prepare DataFrame with catalytic info
        catalysis_info = pd.DataFrame({
            "Step": np.arange(1, len(peaks)+1),
            "Temperature": T[peaks],
            "Mass_at_Peak": M[peaks],
            "Rate_dM_dT": -dM_dT[peaks]  # positive rate
        })

        # Plot DTG curve highlighting catalytic effects
        plt.figure(figsize=(8,5))
        plt.plot(T, dM_dT, label="DTG Curve (dM/dT)", color="red")
        plt.plot(T[peaks], dM_dT[peaks], "ko", label="Catalytic Degradation Steps")
        for r in catalysis_info.itertuples():
            plt.text(r.Temperature, r.Rate_dM_dT,
                    f"Step {r.Step}\nRate={r.Rate_dM_dT:.3f}",
                    fontsize=8, ha="center", va="bottom")
        plt.xlabel("Temperature (°C)")
        plt.ylabel("dM/dT (%/°C)")
        plt.title("Catalytic Effects of Nanoparticles on Degradation")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.show()

        return catalysis_info
    else:
        raise ValueError(f"Application '{application}' not recognized.")'''

    
    """elif application == "oxidation":
        if data_inert is None:
            raise ValueError("For 'oxidation', you must provide 'data_inert' (TGA under inert gas).")
        if "Temp" not in data_inert.columns or "Mass" not in data_inert.columns:
            raise ValueError("'data_inert' must contain 'Temp' and 'Mass' columns.")

        T_inert = data_inert["Temp"].values
        M_inert = data_inert["Mass"].values

        # Interpolate if temperature points do not match
        M_inert_interp = np.interp(T, T_inert, M_inert)

        # Mass difference
        mass_diff = M - M_inert_interp

        plt.figure(figsize=(8,5))
        plt.plot(T, M, label="Oxygen (Oxidative)", color="red")
        plt.plot(T, M_inert_interp, label="Inert (N2/Ar)", color="blue")
        plt.fill_between(T, M_inert_interp, M, color="orange", alpha=0.3, label="Oxidation Loss")
        plt.xlabel("Temperature (°C)")
        plt.ylabel("Mass (%)")
        plt.title("Oxidation Resistance of Nanomaterial")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.legend()
        plt.show()

        # Onset of oxidation: first point where M < M_inert - threshold
        threshold = 0.5  # % mass loss to define onset
        idx = np.where(M < M_inert_interp - threshold)[0]
        if len(idx) > 0:
            onset_temp = T[idx[0]]
        else:
            onset_temp = None

        return {
            "Onset_of_Oxidation_Temperature": onset_temp,
            "Mass_Difference_Oxygen_vs_Inert": mass_diff
        }

    else:
        raise ValueError(f"Application '{application}' not recognized.")"""
    # --- DTG Application ---
    



def UV_Visible_Analysis(data, application, **kwargs):
    """
    Perform multi-mode UV–Visible Spectroscopy analysis for optical and electronic characterization.

    This function provides tools for analyzing UV–Vis absorbance spectra, including
    visualization, Beer–Lambert law concentration estimation, peak identification,
    Landau maximum detection, and Tauc plot-based band gap estimation.

    Parameters
    ----------
    data : pandas.DataFrame or dict
        Experimental UV–Vis dataset containing the columns:
        - 'Wavelength' : float
            Wavelength values in nanometers (nm)
        - 'Absorbance' : float
            Measured absorbance at each wavelength

        Example:
        >>> data = pd.DataFrame({
        ...     "Wavelength": [200, 250, 300, 350],
        ...     "Absorbance": [0.2, 0.8, 1.1, 0.4]
        ... })

    application : str
        Defines the analysis mode. Supported applications:

        - `"plot"` :
            Plot the UV–Vis spectrum (Absorbance vs. Wavelength).

        - `"beer_lambert"` :
            Apply Beer–Lambert law to calculate molar concentration:
            A = ε × l × c, where:
            ε = molar extinction coefficient,
            l = optical path length (cm),
            c = concentration (M).

            Required keyword arguments:
            - `molar_extinction_coefficient` : float
            - `path_length` : float, optional (default=1.0)

        - `"peak_detection"` or `"identify_peaks"` :
            Detect spectral peaks using `scipy.signal.find_peaks`.
            Optional keyword arguments:
            - `height` : float, threshold for peak height.
            - `distance` : int, minimum number of points between peaks.

        - `"band_gap"` :
            Generate a Tauc plot for optical band gap determination.
            Uses the relation (αhν)^n vs. hν, where n = 0.5 for direct
            and n = 2 for indirect transitions.

            Keyword arguments:
            - `n` : float, exponent type (default=0.5)

        - `"landau_max"` :
            Identify the wavelength corresponding to maximum absorbance
            (Landau maximum). If Beer–Lambert parameters are provided,
            the function estimates the sample concentration at that point.

            Optional keyword arguments:
            - `molar_extinction_coefficient` : float
            - `path_length` : float, optional (default=1.0)

    Keyword Arguments
    -----------------
    molar_extinction_coefficient : float, optional
        Required for Beer–Lambert law or Landau Max concentration estimation.
    path_length : float, default=1.0
        Optical path length of the cuvette (in cm).
    height : float, optional
        Minimum absorbance for peak detection.
    distance : int, optional
        Minimum distance between adjacent detected peaks.
    n : float, default=0.5
        Exponent in the Tauc plot for direct/indirect band gap transitions.

    Returns
    -------
    object
        Depends on the analysis mode:

        - `"plot"` :
            Displays spectrum; returns `None`.
        - `"beer_lambert"` :
            pandas.DataFrame with calculated concentration values.
        - `"peak_detection"` / `"identify_peaks"` :
            pandas.DataFrame listing detected peak wavelengths and absorbances.
        - `"band_gap"` :
            pandas.DataFrame with photon energy and Tauc Y-values.
        - `"landau_max"` :
            dict with wavelength, absorbance, and (if applicable) concentration.

    Raises
    ------
    ValueError
        If input format or application type is invalid.
    KeyError
        If required columns ('Wavelength', 'Absorbance') are missing.

    Notes
    -----
    - Band gap energy (Eg) is estimated by extrapolating the linear portion
      of the Tauc plot to the energy axis.
    - The Landau maximum provides insights into π–π* or n–π* transitions.
    - Beer–Lambert analysis assumes linearity in the absorbance–concentration range.
    - Wavelengths must be sorted in ascending order for accurate results.

    Examples
    --------
    >>> data = pd.DataFrame({
    ...     "Wavelength": np.linspace(200, 800, 300),
    ...     "Absorbance": np.exp(-((np.linspace(200, 800, 300) - 400) / 50)**2)
    ... })
    >>> UV_Visible_Analysis(data, "plot")
    # Displays the UV–Vis spectrum.

    >>> UV_Visible_Analysis(data, "peak_detection", height=0.2)
    # Detects and highlights spectral peaks.

    >>> UV_Visible_Analysis(data, "beer_lambert",
    ...     molar_extinction_coefficient=15000, path_length=1.0)
    # Computes sample concentration using Beer–Lambert law.

    >>> UV_Visible_Analysis(data, "band_gap", n=0.5)
    # Displays the Tauc plot for band gap estimation.

    >>> UV_Visible_Analysis(data, "landau_max",
    ...     molar_extinction_coefficient=20000, path_length=1.0)
    # Identifies Landau maximum and estimates concentration.
    """
    
    # Convert to DataFrame
    if isinstance(data, dict):
        df = pd.DataFrame(data)
    elif isinstance(data, pd.DataFrame):
        df = data.copy()
    else:
        raise ValueError("Data must be a pandas DataFrame or dict with keys 'Wavelength' and 'Absorbance'")
    
    if 'Wavelength' not in df.columns or 'Absorbance' not in df.columns:
        raise ValueError("Data must have 'Wavelength' and 'Absorbance' columns")
    
    # ---------------- PLOTTING ----------------
    if application == "plot":
        plt.figure(figsize=(8,5))
        plt.plot(df['Wavelength'], df['Absorbance'], color='blue', lw=2)
        plt.title("UV-Vis Spectrum")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Absorbance")
        plt.grid(True)
        plt.show()
    
    # ---------------- BEER-LAMBERT ----------------
    elif application == "beer_lambert":
        epsilon = kwargs.get('molar_extinction_coefficient')
        path_length = kwargs.get('path_length', 1.0)
        if epsilon is None:
            raise ValueError("For Beer–Lambert calculation, provide 'molar_extinction_coefficient'.")
        
        df['Concentration (M)'] = df['Absorbance'] / (epsilon * path_length)
        print("Calculated concentrations (M):")
        print(df[['Wavelength', 'Absorbance', 'Concentration (M)']])
        
        plt.figure(figsize=(8,5))
        plt.plot(df['Concentration (M)'], df['Absorbance'], "bo-", lw=2)
        plt.title("Absorbance vs Concentration")
        plt.xlabel("Concentration (M)")
        plt.ylabel("Absorbance")
        plt.grid(True)
        plt.show()
        return df
    
    # ---------------- PEAK DETECTION / IDENTIFY ----------------
    elif application in ["peak_detection", "identify_peaks"]:
        height = kwargs.get('height', None)
        distance = kwargs.get('distance', None)
        
        peaks, _ = find_peaks(df['Absorbance'], height=height, distance=distance)
        peak_wavelengths = df['Wavelength'].iloc[peaks].values
        peak_absorbances = df['Absorbance'].iloc[peaks].values
        
        print("Detected peaks (Wavelength, Absorbance):")
        for wl, ab in zip(peak_wavelengths, peak_absorbances):
            print(f"{wl} nm : {ab}")
        
        plt.figure(figsize=(8,5))
        plt.plot(df['Wavelength'], df['Absorbance'], color='blue', lw=2)
        plt.plot(peak_wavelengths, peak_absorbances, "rx", markersize=8, label="Detected Peaks")
        plt.title("UV-Vis Spectrum with Detected Peaks")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Absorbance")
        plt.legend()
        plt.grid(True)
        plt.show()
        return pd.DataFrame({"Wavelength": peak_wavelengths, "Absorbance": peak_absorbances})
    
    # ---------------- BAND GAP ----------------
    elif application == "band_gap":
        n = kwargs.get('n', 0.5)  # 0.5 for direct, 2 for indirect
        # Convert wavelength (nm) to photon energy (eV)
        hv = 1240 / df['Wavelength']  # eV
        # Tauc plot y-axis
        tauc_y = (df['Absorbance'] * hv) ** n
        
        plt.figure(figsize=(8,5))
        plt.plot(hv, tauc_y, "bo-", lw=2)
        plt.title("Tauc Plot for Band Gap Determination")
        plt.xlabel("Photon Energy (eV)")
        plt.ylabel("(Absorbance × hν)^n")
        plt.grid(True)
        plt.show()
        
        print("Use the linear region in the plot to extrapolate x-axis to determine band gap energy Eg (eV).")
        return pd.DataFrame({"Photon Energy (eV)": hv, "Tauc Y": tauc_y})
    
    # ---------------- LANDAU MAX ----------------
    elif application == "landau_max":
        # Find max absorbance
        idx_max = df['Absorbance'].idxmax()
        max_wavelength = df.loc[idx_max, 'Wavelength']
        max_absorbance = df.loc[idx_max, 'Absorbance']
        
        print(f"Landau Max: {max_wavelength} nm with Absorbance {max_absorbance}")
        
        epsilon = kwargs.get('molar_extinction_coefficient')
        path_length = kwargs.get('path_length', 1.0)
        concentration = None
        if epsilon is not None:
            concentration = max_absorbance / (epsilon * path_length)
            print(f"Estimated concentration at Landau Max using Beer–Lambert Law: {concentration} M")
        
        plt.figure(figsize=(8,5))
        plt.plot(df['Wavelength'], df['Absorbance'], color='blue', lw=2)
        plt.axvline(max_wavelength, color='red', linestyle='--', label=f"Landau Max: {max_wavelength} nm")
        plt.scatter([max_wavelength], [max_absorbance], color='red', zorder=5)
        plt.title("UV-Vis Spectrum with Landau Max")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Absorbance")
        plt.legend()
        plt.grid(True)
        plt.show()
        
        return {"Wavelength (nm)": max_wavelength, "Absorbance": max_absorbance, "Concentration (M)": concentration}
    
    else:
        raise ValueError(f"Application '{application}' not implemented yet.")



def CV(data, application):
    """
    Perform Cyclic Voltammetry (CV) data analysis for electrochemical characterization.

    This function provides core analytical and visualization tools for cyclic
    voltammetry experiments, including voltammogram plotting, oxidation/reduction
    peak detection, and peak shape analysis for assessing reversibility of redox
    processes.

    Parameters
    ----------
    data : list of tuples, list of lists, or pandas.DataFrame
        Experimental CV dataset containing the following columns or structure:
        - 'E' : float
            Applied potential (V vs. reference electrode)
        - 'I' : float
            Measured current (A)

        Example:
        >>> data = pd.DataFrame({
        ...     "E": [-0.5, -0.3, 0.0, 0.3, 0.5],
        ...     "I": [-0.0001, 0.0003, 0.0012, 0.0005, -0.0002]
        ... })

    application : str
        Defines the analysis type. Supported options include:

        - `"plot"` :
            Display the cyclic voltammogram (current vs. potential).

        - `"peaks"` :
            Detect and highlight oxidation and reduction peaks using
            `scipy.signal.find_peaks` with a default prominence of 0.001 A.
            The function identifies the most intense oxidation peak and
            up to two reduction peaks.

        - `"shape"` :
            Analyze the shape and symmetry of oxidation/reduction peaks to
            determine the reversibility of the redox process.
            It computes:
                - E_pa : anodic (oxidation) peak potential (V)
                - E_pc : cathodic (reduction) peak potential (V)
                - ΔE_p : peak separation (V)
                - |I_pc/I_pa| : peak current ratio

            Based on electrochemical theory:
            - Reversible systems exhibit ΔE_p ≈ 59 mV/n (for one-electron transfer)
              and |I_pc/I_pa| ≈ 1.
            - Quasi-reversible systems show moderate deviations.
            - Irreversible systems display large separations and asymmetric peaks.

    Returns
    -------
    None
        The function primarily displays visualizations and prints analysis
        results directly to the console.

    Raises
    ------
    TypeError
        If the input data format is invalid.
    ValueError
        If the specified application is not supported.

    Notes
    -----
    - Ensure that potentials (E) are in ascending or cyclic order for
      accurate peak detection.
    - Peak prominence and smoothing parameters can be tuned for noisy data.
    - The reversibility classification is heuristic and assumes one-electron
      transfer unless otherwise known.

    Examples
    --------
    >>> data = pd.DataFrame({
    ...     "E": np.linspace(-0.5, 0.5, 200),
    ...     "I": 0.001 * np.sin(4 * np.pi * np.linspace(-0.5, 0.5, 200))
    ... })
    >>> CV(data, "plot")
    # Displays the cyclic voltammogram.

    >>> CV(data, "peaks")
    # Detects and highlights oxidation/reduction peaks.

    >>> CV(data, "shape")
    # Computes ΔEp and |Ipc/Ipa| to infer redox reversibility.
    """
    
    # Convert list input to DataFrame for consistency
    if isinstance(data, list):
        data = pd.DataFrame(data, columns=["E", "I"])
    elif not isinstance(data, pd.DataFrame):
        raise TypeError("Data must be a list of (E,I) pairs or a pandas DataFrame.")
    
    E = data["E"].values
    I = data["I"].values
    
    # -------------------------------
    # Application 1: Simple CV plot
    # -------------------------------
    if application == "plot":
        plt.figure(figsize=(6,5))
        plt.plot(E, I, color="blue", lw=1.5)
        plt.xlabel("Potential (V vs Ref)")
        plt.ylabel("Current (A)")
        plt.title("Cyclic Voltammogram")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()
    
    # -------------------------------
    # Application 2: Peak detection
    # -------------------------------
    elif application == "peaks":
        # Oxidation peaks = local maxima
        ox_peaks, _ = find_peaks(I, prominence=0.001)
        # Reduction peaks = local minima
        red_peaks, _ = find_peaks(-I, prominence=0.001)

        # --- Select strongest peaks ---
        # Keep the highest oxidation peak (max current)
        if len(ox_peaks) > 0:
            ox_peaks = [ox_peaks[np.argmax(I[ox_peaks])]]

        # Keep the two most negative reduction peaks (lowest currents)
        if len(red_peaks) > 2:
            red_peaks = red_peaks[np.argsort(I[red_peaks])[:2]]

        # --- Plot ---
        plt.figure(figsize=(6,5))
        plt.plot(E, I, color="black", lw=1.5, label="CV curve")
        plt.scatter(E[ox_peaks], I[ox_peaks], color="red", s=70, label="Oxidation peak")
        plt.scatter(E[red_peaks], I[red_peaks], color="blue", s=70, label="Reduction peaks")

        plt.xlabel("Potential (V vs Ref)")
        plt.ylabel("Current (A)")
        plt.title("Cyclic Voltammetry with Selected Peaks")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()

        # --- Print results ---
        print("Oxidation peak (Potential, Current):")
        for e, i in zip(E[ox_peaks], I[ox_peaks]):
            print(f"  {e:.3f} V, {i:.5f} A")

        print("\nReduction peaks (Potential, Current):")
        for e, i in zip(E[red_peaks], I[red_peaks]):
            print(f"  {e:.3f} V, {i:.5f} A")
    
    # -------------------------------
    # Application 3: Shape analysis
    # -------------------------------
    elif application == "shape":
        ox_peaks, _ = find_peaks(I, prominence=0.001)
        red_peaks, _ = find_peaks(-I, prominence=0.001)
        
        plt.figure(figsize=(7,6))
        plt.plot(E, I, color="blue", lw=1.5, label="CV curve")
        plt.scatter(E[ox_peaks], I[ox_peaks], color="red", s=70, label="Oxidation peaks")
        plt.scatter(E[red_peaks], I[red_peaks], color="green", s=70, label="Reduction peaks")
        
        for e, i in zip(E[ox_peaks], I[ox_peaks]):
            plt.annotate(f"{e:.2f}V", (e, i), textcoords="offset points", xytext=(0,10), ha="center", color="red")
        for e, i in zip(E[red_peaks], I[red_peaks]):
            plt.annotate(f"{e:.2f}V", (e, i), textcoords="offset points", xytext=(0,-15), ha="center", color="green")
        
        plt.xlabel("Potential (V vs Ref)")
        plt.ylabel("Current (A)")
        plt.title("Cyclic Voltammetry: Shape of Peaks (Reversibility)")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()
        
        print("\n--- Shape / Reversibility Analysis ---")
        for idx in range(min(len(ox_peaks), len(red_peaks))):
            E_pa, I_pa = E[ox_peaks[idx]], I[ox_peaks[idx]]
            E_pc, I_pc = E[red_peaks[idx]], I[red_peaks[idx]]
            
            delta_Ep = abs(E_pc - E_pa)
            current_ratio = abs(I_pc) / abs(I_pa) if I_pa != 0 else np.nan
            
            print(f"\nPair {idx+1}:")
            print(f"  E_pa = {E_pa:.3f} V, I_pa = {I_pa:.5f} A")
            print(f"  E_pc = {E_pc:.3f} V, I_pc = {I_pc:.5f} A")
            print(f"  ΔEp = {delta_Ep*1000:.1f} mV")
            print(f"  |Ipc/Ipa| = {current_ratio:.2f}")
            
            # Classification
            if abs(delta_Ep - 0.0592) < 0.02 and abs(current_ratio - 1) < 0.2:
                print("  → Likely REVERSIBLE (mirror-like peaks, fast electron transfer)")
            elif delta_Ep < 0.2 and 0.5 < current_ratio < 1.5:
                print("  → Likely QUASI-REVERSIBLE (moderate electron transfer)")
            else:
                print("  → Likely IRREVERSIBLE (asymmetric peaks, slow transfer or side reactions)")
    
    else:
        raise ValueError(f"Unknown application: {application}")




def SAXS_Analysis(data, application):
    """
    Perform Small-Angle X-ray Scattering (SAXS) data analysis for nanostructural characterization.

    This function provides key analytical tools to extract structural information
    from SAXS profiles, including visualization, peak position detection, intensity
    integration, peak width (FWHM) determination, and Guinier (radius of gyration) analysis.

    Parameters
    ----------
    data : pandas.DataFrame
        Experimental SAXS dataset containing:
        - 'q' : float
            Scattering vector magnitude (1/nm)
        - 'I' : float
            Scattered intensity I(q)

        Example:
        >>> data = pd.DataFrame({
        ...     "q": [0.01, 0.02, 0.03, 0.04],
        ...     "I": [300, 800, 400, 200]
        ... })

    application : str
        Defines the type of analysis to perform. Supported options include:

        - `"plot"` :
            Plot I(q) vs. q to visualize the SAXS curve and scattering profile.

        - `"peak_position"` :
            Identify the q position of the main scattering peak and calculate
            the corresponding real-space characteristic spacing:
            d = 2π / q_peak

        - `"peak_intensity"` :
            Quantify the intensity and integrated area under the most intense
            scattering peak using numerical integration (`numpy.trapz`).

        - `"peak_width"` :
            Compute the full width at half maximum (FWHM) of the main scattering
            peak, which provides information about domain size and order distribution.

        - `"rog"` :
            Perform Guinier analysis (low-q region) by linear fitting of
            ln I(q) vs. q² to estimate the radius of gyration (Rg) and I(0):
                ln I(q) = ln I(0) − (Rg² * q²) / 3

    Returns
    -------
    tuple or None
        Depending on the analysis:
        - `"peak_position"` → (q_peak, d_spacing)
        - `"peak_intensity"` → (I_peak, area)
        - `"peak_width"` → FWHM
        - `"rog"` → (Rg, I0)
        - `"plot"` → None

    Raises
    ------
    TypeError
        If the input data is not a pandas DataFrame or lacks required columns.
    ValueError
        If the specified application is not supported or no peaks are detected.

    Notes
    -----
    - q is defined as q = (4π/λ) sin(θ), where θ is half the scattering angle.
    - The characteristic spacing (d) corresponds to periodicity or average
      interparticle distance.
    - FWHM can be used to estimate crystalline order (via Scherrer-like relations).
    - The Guinier approximation is valid only for q·Rg < 1.3.

    Examples
    --------
    >>> SAXS_Analysis(data, "plot")
    # Displays the SAXS intensity profile.

    >>> SAXS_Analysis(data, "peak_position")
    # Prints and plots q_peak and corresponding d-spacing.

    >>> SAXS_Analysis(data, "peak_intensity")
    # Calculates peak height and integrated scattering area.

    >>> SAXS_Analysis(data, "peak_width")
    # Determines full width at half maximum (FWHM) in q-space.

    >>> SAXS_Analysis(data, "rog")
    # Performs Guinier analysis to estimate radius of gyration (Rg).
    """

    if application == "plot":
        plt.figure(figsize=(6,5))
        plt.plot(data['q'], data['I'], 'o-', markersize=4, label="SAXS Data")
        plt.xlabel("q (1/nm)")
        plt.ylabel("Intensity I(q)")
        plt.title("SAXS Curve")
        plt.legend()
        plt.grid(True, which="both", ls="--", lw=0.5)
        #plt.xlim(0, 0.08)
        #plt.ylim(0, 4000)
        plt.show()

    elif application == "peak_position":
        idx_max = data['I'].idxmax()
        q_peak = data.loc[idx_max, 'q']
        I_peak = data.loc[idx_max, 'I']
        d_spacing = 2 * np.pi / q_peak

        print(f"q_peak = {q_peak:.4f} 1/nm (I = {I_peak})")
        print(f"Characteristic spacing d = {d_spacing:.2f} nm")

        plt.figure(figsize=(6,5))
        plt.plot(data['q'], data['I'], 'o-', markersize=4, label="SAXS Data")
        plt.axvline(q_peak, color='r', linestyle='--', label=f"q_peak = {q_peak:.4f}")
        plt.text(q_peak, I_peak, f"d = {d_spacing:.2f} nm", 
                 rotation=90, va='bottom', ha='right', color='r')
        plt.xlabel("q (1/nm)")
        plt.ylabel("Intensity I(q)")
        plt.title("SAXS Peak Position and Spacing")
        plt.legend()
        plt.grid(True, which="both", ls="--", lw=0.5)
        #plt.xlim(0, 0.08)
        #plt.ylim(0, 4000)
        plt.show()

        return q_peak, d_spacing

    elif application == "peak_intensity":
        from scipy.signal import find_peaks, peak_widths

        # Find peaks
        peaks, _ = find_peaks(data['I'])
        if len(peaks) == 0:
            print("⚠️ No peaks found.")
            return None, None

        # Select the highest peak
        idx_max = peaks[np.argmax(data['I'].iloc[peaks])]
        q_peak = data.loc[idx_max, 'q']
        I_peak = data.loc[idx_max, 'I']

        # Get peak width at half max
        results_half = peak_widths(data['I'], [idx_max], rel_height=0.5)
        left_idx = int(results_half[2][0])
        right_idx = int(results_half[3][0])

        # Slice the region under the peak
        q_region = data['q'].iloc[left_idx:right_idx+1]
        I_region = data['I'].iloc[left_idx:right_idx+1]

        # Integrate under the peak
        area = np.trapz(I_region, q_region)

        # Print results
        print(f"Peak height (I_peak) = {I_peak:.2f}")
        print(f"Integrated area under peak = {area:.2f}")

        # Plot with shaded peak region
        plt.figure(figsize=(6,5))
        plt.plot(data['q'], data['I'], 'o-', markersize=4, label="SAXS Data")
        plt.fill_between(q_region, I_region, alpha=0.3, color='orange', label="Peak area")
        plt.axvline(q_peak, color='r', linestyle='--', label=f"q_peak = {q_peak:.4f}")
        plt.scatter(q_peak, I_peak, color='red', zorder=5)
        plt.text(q_peak, I_peak, f"Peak = {I_peak:.1f}", va='bottom', ha='left', color='r')

        plt.xlabel("q (1/nm)")
        plt.ylabel("Intensity I(q)")
        plt.title("SAXS Peak Intensity")
        plt.legend()
        plt.grid(True, which="both", ls="--", lw=0.5)
        plt.show()

        return I_peak, area

    elif application == "peak_width":
        from scipy.signal import find_peaks, peak_widths
        # --- Find all peaks ---
        peaks, _ = find_peaks(data['I'])
        if len(peaks) == 0:
            print("⚠️ No peaks found.")
            return None

        # --- Select the highest peak ---
        idx_max = peaks[np.argmax(data['I'].iloc[peaks])]
        q_peak = data.loc[idx_max, 'q']
        I_peak = data.loc[idx_max, 'I']

        # --- Calculate FWHM using scipy ---
        results_half = peak_widths(data['I'], [idx_max], rel_height=0.5)

        # Width in index space → convert to q units
        FWHM = results_half[0][0] * (data['q'].iloc[1] - data['q'].iloc[0])  

        # Left & right positions at half max
        q_left = data['q'].iloc[int(results_half[2][0])]
        q_right = data['q'].iloc[int(results_half[3][0])]

        # Print results
        print(f"q_peak = {q_peak:.4f} 1/nm")
        print(f"I_peak = {I_peak:.2f}")
        print(f"FWHM = {FWHM:.4f} 1/nm")

        # --- Plot ---
        plt.figure(figsize=(6,5))
        plt.plot(data['q'], data['I'], 'o-', markersize=4, label="SAXS Data")
        plt.axhline(I_peak/2, color='g', linestyle='--', label="Half max")
        plt.axvline(q_left, color='r', linestyle='--')
        plt.axvline(q_right, color='r', linestyle='--')
        plt.scatter([q_peak], [I_peak], color='red', zorder=5, label="q_peak")
        plt.text(q_peak, I_peak, f"FWHM = {FWHM:.4f}", va='bottom', ha='center', color='r')

        plt.xlabel("q (1/nm)")
        plt.ylabel("Intensity I(q)")
        plt.title("SAXS Peak Width (FWHM)")
        plt.legend()
        plt.grid(True, which="both", ls="--", lw=0.5)
        plt.show()

        return FWHM

    elif application == "rog":
        # Use low-q data (heuristic: q < 0.03 1/nm)
        low_q = data[data['q'] < 0.03].copy()
        q2 = low_q['q']**2
        lnI = np.log(low_q['I'])

        # Linear regression: lnI = intercept + slope * q^2
        coeffs = np.polyfit(q2, lnI, 1)
        slope, intercept = coeffs[0], coeffs[1]

        Rg = np.sqrt(-3 * slope)
        I0 = np.exp(intercept)

        print(f"Slope = {slope:.4f}")
        print(f"Intercept = {intercept:.4f}")
        print(f"Radius of gyration Rg = {Rg:.2f} nm")
        print(f"I(0) = {I0:.2f}")

        # Plot Guinier plot
        plt.figure(figsize=(6,5))
        plt.plot(q2, lnI, 'o', label="Data (ln I vs q²)")
        plt.plot(q2, slope*q2 + intercept, 'r-', label="Guinier Fit")
        plt.xlabel("q² (1/nm²)")
        plt.ylabel("ln I(q)")
        plt.title("Guinier Analysis")
        plt.legend()
        plt.grid(True, ls="--", lw=0.5)
        plt.show()

        return Rg, I0

    else:
        print(f"Application '{application}' not implemented yet.")




def WAXS_Analysis(data, application, **kwargs):
    """
    Perform Wide-Angle X-ray Scattering (WAXS) data analysis for crystallographic and
    nanostructural characterization.

    This function analyzes WAXS diffraction patterns to determine structural
    information such as peak positions, d-spacings, peak widths, crystallite size,
    degree of crystallinity, and peak shape classification.

    Parameters
    ----------
    data : pandas.DataFrame or array-like
        Experimental WAXS dataset containing two columns:
        - 'q' : float
            Scattering vector (Å⁻¹) or 2θ values (degrees)
        - 'I' : float
            Scattering intensity (a.u.)
        
        Example:
        >>> data = pd.DataFrame({
        ...     "q": [0.5, 1.0, 1.5, 2.0],
        ...     "I": [200, 600, 300, 100]
        ... })

    application : str
        Defines the type of analysis to perform. Supported options include:
        
        - `"plot"` :
            Plot the WAXS pattern (Intensity vs q or 2θ).

        - `"peak_position"` :
            Detect the most intense diffraction peaks, compute their
            corresponding d-spacings using:
                d = 2π / q
            Returns a table of q values, d-spacings, and intensities.

        - `"peak_intensity"` :
            Determine the intensity and integrated area under the
            strongest diffraction peaks, useful for semi-quantitative
            crystallinity assessment.

        - `"peak_width"` :
            Compute full width at half maximum (FWHM) of main peaks and
            estimate crystallite size using the Scherrer equation:
                L = Kλ / (β cosθ)
            Also estimates overall percent crystallinity from integrated peak areas.

        - `"peak_shape"` :
            Classify peak sharpness based on FWHM(2θ) and estimate
            the crystallinity percentage. Sharp peaks imply high
            crystallinity, broad peaks indicate amorphous domains.

    Optional Keyword Arguments
    ---------------------------
    threshold : float, optional
        Minimum relative intensity (fraction of max) to detect peaks.
        Default = 0.1 (10% of max intensity).

    top_n : int, optional
        Number of top peaks to consider. Default = 3.

    wavelength : float, optional
        X-ray wavelength in Ångströms (required for `"peak_width"` and `"peak_shape"`).

    K : float, optional
        Scherrer constant, typically between 0.89–0.94. Default = 0.9.

    width_threshold : float, optional
        Threshold in degrees for classifying peak shapes. Default = 2.0° (2θ).

    Returns
    -------
    pandas.DataFrame or tuple
        Depending on the analysis type:
        - `"peak_position"` → DataFrame of q, d-spacing, and intensity
        - `"peak_intensity"` → DataFrame of peak positions and intensities
        - `"peak_width"` → (DataFrame of peak properties, crystallinity_percent)
        - `"peak_shape"` → (DataFrame of peak classification, crystallinity_percent)
        - `"plot"` → None

    Raises
    ------
    ValueError
        If an unsupported application is specified or if wavelength is missing
        for analyses that require it.
    TypeError
        If the input data format is invalid.

    Notes
    -----
    - q and 2θ are related by: q = (4π / λ) sin(θ)
    - d-spacing provides interplanar distances according to Bragg’s law.
    - Crystallite size estimation assumes negligible strain and instrumental broadening.
    - The degree of crystallinity is estimated from the ratio of crystalline
      (peak) area to total scattered intensity.

    Examples
    --------
    >>> WAXS_Analysis(data, "plot")
    # Displays the WAXS pattern.

    >>> WAXS_Analysis(data, "peak_position")
    # Returns major peaks and corresponding d-spacings.

    >>> WAXS_Analysis(data, "peak_intensity")
    # Calculates integrated areas of main peaks.

    >>> WAXS_Analysis(data, "peak_width", wavelength=1.54)
    # Estimates FWHM, crystallite size, and crystallinity.

    >>> WAXS_Analysis(data, "peak_shape", wavelength=1.54)
    # Classifies peaks as sharp/broad and returns crystallinity percent.
    """
    
    # Convert array-like to DataFrame
    if not isinstance(data, pd.DataFrame):
        data = pd.DataFrame(data, columns=["q", "I"])
    
    x = data["q"].values
    y = data["I"].values
    
    # --------------------- PLOT ---------------------
    if application == "plot":
        plt.figure(figsize=(8,5))
        plt.plot(x, y, color="blue", linewidth=1.5)
        plt.xlabel("Scattering vector q (Å⁻¹) or 2θ (degrees)")
        plt.ylabel("Intensity (a.u.)")
        plt.title("WAXS Pattern")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()
    
    # --------------------- PEAK POSITION ---------------------
    elif application == "peak_position":
        q_vals = data["q"].values     # q in Å⁻¹
        intensity = data["I"].values

        # Peak detection
        threshold = kwargs.get("threshold", 0.1)  # keep peaks above 10% of max
        peaks, props = find_peaks(intensity, height=np.max(intensity) * threshold)

        # Keep only strongest peaks
        top_n = kwargs.get("top_n", 3)
        sorted_idx = np.argsort(props["peak_heights"])[::-1][:top_n]
        peaks = peaks[sorted_idx]

        # d-spacing using q (Å⁻¹)
        d_spacings = (2 * np.pi) / q_vals[peaks]

        # Results table
        results = pd.DataFrame({
            "q (Å⁻¹)": q_vals[peaks],
            "d_spacing (Å)": d_spacings,
            "Intensity": intensity[peaks]
        }).sort_values("q (Å⁻¹)").reset_index(drop=True)

        print("Main Peaks and d-spacings:")
        print(results)

        # Plot
        plt.figure(figsize=(8, 5))
        plt.plot(q_vals, intensity, label="WAXS Pattern")
        plt.plot(q_vals[peaks], intensity[peaks], "rx", label="Main Peaks")
        plt.xlabel("q (Å⁻¹)")
        plt.ylabel("Intensity (a.u.)")
        plt.title("WAXS Peaks (q-based)")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()

        return results
    
    # --------------------- PEAK INTENSITY ---------------------
    elif application == "peak_intensity":
        top_n = kwargs.get("top_n", 3)   # default keep 3 strongest peaks
        threshold = kwargs.get("threshold", 0.05)  # 5% of max intensity
        
        # Peak detection
        peaks, props = find_peaks(y, height=np.max(y) * threshold)
        
        # Keep top_n strongest peaks
        sorted_idx = np.argsort(props["peak_heights"])[::-1][:top_n]
        peaks = peaks[sorted_idx]
        
        peak_positions = x[peaks]
        peak_intensities = y[peaks]
        
        # Degree of crystallinity
        Ac = simpson(peak_intensities, peak_positions)
        Aa = simpson(y, x) - Ac
        #Xc = (Ac / (Ac + Aa)) * 100
        
        results = pd.DataFrame({
            "Peak_2θ(deg)": peak_positions,
            "Intensity(a.u.)": peak_intensities
        }).sort_values("Peak_2θ(deg)").reset_index(drop=True)
        
        print("Main Peak Intensities:")
        print(results)
        #print(f"\nEstimated Degree of Crystallinity (Xc): {Xc:.2f}%")
        
        plt.figure(figsize=(8,5))
        plt.plot(x, y, label="WAXS Pattern")
        plt.plot(peak_positions, peak_intensities, "rx", label="Main Peaks")
        plt.xlabel("q")
        plt.ylabel("Intensity (a.u.)")
        plt.title("WAXS Peaks and Intensities")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()
        
        return results, #Xc

# --------------------- PEAK WIDTH / CRYSTALLITE SIZE ---------------------
   # --------------------- PEAK WIDTH + CRYSTALLINITY ---------------------
    elif application == "peak_width":
        if "wavelength" not in kwargs:
            raise ValueError("Please provide X-ray wavelength in Angstroms using wavelength=<value>.")
        wavelength = kwargs["wavelength"]
        K = kwargs.get("K", 0.9)
        top_n = kwargs.get("top_n", 3)
        threshold = kwargs.get("threshold", 0.05)

        # Find peaks
        peaks, props = find_peaks(y, height=np.max(y) * threshold)
        sorted_idx = np.argsort(props["peak_heights"])[::-1][:top_n]
        peaks = peaks[sorted_idx]
        q_peaks = x[peaks]   # x is q (Å⁻¹)

        # FWHM calculation
        results_half = peak_widths(y, peaks, rel_height=0.5)
        FWHM_points = results_half[0]
        dq = x[1] - x[0]   # step size in q
        FWHM_q = FWHM_points * dq   # convert width in index → width in q

        # Convert q to theta (radians)
        theta_rad = np.arcsin((q_peaks * wavelength) / (4 * np.pi))

        # Convert FWHM in q to FWHM in radians of 2θ
        # Approximation: d(2θ)/dq = λ / (2 cosθ)
        FWHM_rad = (wavelength / (2 * np.cos(theta_rad))) * FWHM_q

        # Scherrer equation
        L = K * wavelength / (FWHM_rad * np.cos(theta_rad))

        # Peak areas
        peak_areas = []
        for i, p in enumerate(peaks):
            left, right = int(results_half[2][i]), int(results_half[3][i])
            peak_area = np.trapz(y[left:right], x[left:right])
            peak_areas.append(peak_area)

        total_area = np.trapz(y, x)
        crystallinity_percent = (np.sum(peak_areas) / total_area) * 100

        results = pd.DataFrame({
            "q (Å⁻¹)": q_peaks,
            "Theta (deg)": np.degrees(theta_rad),
            "FWHM(q)": FWHM_q,
            "Crystallite_Size(Å)": L,
            "Peak_Area": peak_areas
        }).sort_values("q (Å⁻¹)").reset_index(drop=True)

        print("Main Peak Widths, Crystallite Sizes, and Areas:")
        print(results)
        print(f"\nEstimated Percent Crystallinity = {crystallinity_percent:.2f}%")

        # Plot
        plt.figure(figsize=(8,5))
        plt.plot(x, y, label="WAXS Pattern")
        plt.plot(q_peaks, y[peaks], "rx", label="Main Peaks")
        for i, p in enumerate(peaks):
            plt.fill_between(x[int(results_half[2][i]):int(results_half[3][i])],
                            y[int(results_half[2][i]):int(results_half[3][i])],
                            alpha=0.3, label=f"Peak {i+1} Area")
        plt.xlabel("q (Å⁻¹)")
        plt.ylabel("Intensity (a.u.)")
        plt.title("WAXS Peaks with FWHM & Crystallinity")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()

        return results, crystallinity_percent


       # --------------------- PEAK SHAPE / BACKGROUND ---------------------
    elif application == "peak_shape":
        if "wavelength" not in kwargs:
            raise ValueError("Please provide X-ray wavelength in Angstroms using wavelength=<value>.")
        wavelength = kwargs["wavelength"]
        
        top_n = kwargs.get("top_n", 3)
        threshold = kwargs.get("threshold", 0.05)
        
        # Peak detection
        peaks, props = find_peaks(y, height=np.max(y) * threshold)
        sorted_idx = np.argsort(props["peak_heights"])[::-1][:top_n]
        peaks = peaks[sorted_idx]
        q_peaks = x[peaks]

        # FWHM in q-units
        results_half = peak_widths(y, peaks, rel_height=0.5)
        dq = x[1] - x[0]
        FWHM_q = results_half[0] * dq

        # Convert q → theta
        theta_rad = np.arcsin((q_peaks * wavelength) / (4 * np.pi))

        # Convert FWHM(q) → FWHM in 2θ (radians), then to degrees
        FWHM_rad = (wavelength / (2 * np.cos(theta_rad))) * FWHM_q
        FWHM_deg = np.degrees(FWHM_rad)

        # Peak areas & crystallinity
        peak_areas = []
        for i, p in enumerate(peaks):
            left, right = int(results_half[2][i]), int(results_half[3][i])
            peak_area = np.trapz(y[left:right], x[left:right])
            peak_areas.append(peak_area)
        total_area = np.trapz(y, x)
        crystallinity_percent = (np.sum(peak_areas) / total_area) * 100

        # Classification
        width_threshold = kwargs.get("width_threshold", 2.0)  # default 2° 2θ
        shape = [
            "Sharp (Highly Crystalline)" if w < width_threshold else
            "Broad (Mostly Amorphous)"
            for w in FWHM_deg
        ]
        
        results = pd.DataFrame({
            "q (Å⁻¹)": q_peaks,
            "Theta (deg)": np.degrees(theta_rad),
            "FWHM(2θ deg)": FWHM_deg,
            "Shape": shape,
            "Peak_Area": peak_areas
        }).sort_values("q (Å⁻¹)").reset_index(drop=True)
        
        print("Main Peak Shapes (θ-based):")
        print(results)
        print(f"\nEstimated Percent Crystallinity = {crystallinity_percent:.2f}%")
        
        # Plot
        plt.figure(figsize=(8,5))
        plt.plot(x, y, label="WAXS Pattern")
        plt.plot(q_peaks, y[peaks], "rx", label="Main Peaks")
        for i, p in enumerate(peaks):
            plt.fill_between(x[int(results_half[2][i]):int(results_half[3][i])],
                            y[int(results_half[2][i]):int(results_half[3][i])],
                            alpha=0.3, label=f"Peak {i+1} Area")
        plt.xlabel("q (Å⁻¹)")
        plt.ylabel("Intensity (a.u.)")
        plt.title("WAXS Peak Shapes (θ-based)")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.show()
        
        return results, crystallinity_percent


    
    else:
        raise ValueError(f"Application '{application}' not supported yet.")







def NMR_Analysis(df, application,peak_regions=None,peak_info=None):
    """
    Analyze and visualize ¹H NMR spectra for different applications.

    This function provides multiple modes:
    1. Plotting the raw NMR spectrum (`'plot'`).
    2. Plotting the spectrum with integrated peak steps (`'plot_with_integrals'`).
    3. Estimating mole fractions of compounds in a mixture (`'mixture_composition'`).
    4. Calculating percentage impurity of a compound (`'calculate_impurity'`).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing NMR data with columns:
        - 'ppm' : chemical shift values (x-axis)
        - 'Spectrum' : intensity values (y-axis)
    application : str
        Mode of operation. Options:
        - 'plot' : generates a professional NMR spectrum plot.
        - 'plot_with_integrals' : generates a plot with integral steps (requires `peak_regions`).
        - 'mixture_composition' : calculates mole fractions of compounds (requires `peak_info`).
        - 'calculate_impurity' : calculates impurity percentage (requires `peak_info` with main and impurity info).
    peak_regions : dict, optional
        Dictionary specifying integration regions for peaks (required for `'plot_with_integrals'`).
        Format: {region_name: (start_ppm, end_ppm)}
    peak_info : dict, optional
        Dictionary with compound information for mixture analysis or impurity calculation.
        For `'mixture_composition'`:
            {compound_name: {'region': (start_ppm, end_ppm), 'protons': int}}
        For `'calculate_impurity'`:
            {
                'main_compound': {'region': (start, end), 'protons': int},
                'impurity': {'region': (start, end), 'protons': int}
            }

    Returns
    -------
    None
        The function either displays plots or prints calculated results.

    Examples
    --------
    # 1. Simple plot of NMR spectrum
    >>> NMR_Analysis(df, application='plot')

    # 2. Plot spectrum with integrals
    >>> peak_regions = {'peak1': (7.0, 7.5), 'peak2': (3.5, 4.0)}
    >>> NMR_Analysis(df, application='plot_with_integrals', peak_regions=peak_regions)

    # 3. Mixture composition analysis
    >>> peak_info = {
    ...     'CompoundA': {'region': (7.0, 7.5), 'protons': 5},
    ...     'CompoundB': {'region': (3.5, 4.0), 'protons': 3}
    ... }
    >>> NMR_Analysis(df, application='mixture_composition', peak_info=peak_info)
    --- Mixture Composition ---
    Mole Fraction of CompoundA: 0.62
    Mole Fraction of CompoundB: 0.38

    # 4. Impurity calculation
    >>> peak_info = {
    ...     'main_compound': {'region': (7.0, 7.5), 'protons': 5},
    ...     'impurity': {'region': (3.5, 4.0), 'protons': 1}
    ... }
    >>> NMR_Analysis(df, application='calculate_impurity', peak_info=peak_info)
    --- Impurity Analysis ---
    Main Compound Integral per Proton: 0.1234
    Impurity Integral per Proton: 0.0123
    Estimated Impurity: 9.09%
    """
    if application == 'plot':
        # Create the plot
        plt.figure(figsize=(12, 6))
        plt.plot(df['ppm'], df['Spectrum'], color='darkblue', linewidth=1)
        
        # Customize the plot for a professional look
        plt.title('¹H NMR Spectrum', fontsize=16, fontweight='bold')
        plt.xlabel('Chemical Shift (ppm)', fontsize=12)
        plt.ylabel('Intensity', fontsize=12)
        plt.gca().invert_xaxis()  # Invert the x-axis, which is standard for NMR spectra
        plt.grid(True, linestyle='--', alpha=0.6)
        
        # Add a light gray background to the plot area
        plt.gca().set_facecolor('#f7f7f7')
        
        # Remove the top and right spines for a cleaner look
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.show()

    elif application == 'plot_with_integrals':
        if not peak_regions:
            print("Error: 'peak_regions' dictionary must be provided for 'plot_with_integrals' application.")
            return

        integrals = {}
        for region, (start, end) in peak_regions.items():
            peak_data = df[(df['ppm'] >= start) & (df['ppm'] <= end)]
            if not peak_data.empty:
                # Use the new simpson function
                area = simpson(peak_data['Spectrum'], x=peak_data['ppm'])
                integrals[region] = abs(area)

        if not integrals:
            print("No peaks found in the defined regions.")
            return

        min_integral = min(integrals.values())
        proton_ratios = {region: round(area / min_integral) for region, area in integrals.items()}

        # Start plotting with integral steps
        fig, ax1 = plt.subplots(figsize=(12, 6))
        
        # Plot the NMR spectrum
        ax1.plot(df['ppm'], df['Spectrum'], color='darkblue', linewidth=1, label='NMR Spectrum (scaled)')
        ax1.set_xlabel('Chemical Shift (ppm)', fontsize=12)
        ax1.set_ylabel('Intensity (a.u.)', fontsize=12)
        ax1.invert_xaxis()
        ax1.grid(True, linestyle='--', alpha=0.6)
        ax1.set_facecolor('#f7f7f7')
        
        # Create a second y-axis for the integral
        ax2 = ax1.twinx()
        
        # Build the integral step plot
        integral_y = np.zeros(len(df))
        
        # Use a list to store the total integral for each step
        total_integral_steps = [0]
        
        for region, (start, end) in peak_regions.items():
            ratio = proton_ratios[region]
            
            # Find the index of the start of the current peak region
            start_idx = df['ppm'].sub(end).abs().idxmin()
            
            # Add the new integral value to the step list
            total_integral_steps.append(total_integral_steps[-1] + ratio)
            
            # Update the integral step for the points after the peak
            integral_y[start_idx:] = total_integral_steps[-1]
            
            # Add annotation for the integral
            mid_ppm = (start + end) / 2
            ax2.annotate(f"$\\Delta$Integral = {ratio} H", 
                         xy=(mid_ppm, total_integral_steps[-2] + ratio / 2),
                         xycoords='data',
                         textcoords='data',
                         ha='center',
                         fontsize=10,
                         color='red',
                         fontweight='bold')
            
        # Plot the integral steps
        ax2.plot(df['ppm'], integral_y, drawstyle='steps-post', color='red', linewidth=2, label='Integral (stepped)')
        ax2.set_ylabel('Integral (H units)', color='red', fontsize=12)
        ax2.tick_params(axis='y', colors='red')
        
        # Set titles and legends
        plt.title('Annotated Simulated ¹H NMR Spectrum - Integrals', fontsize=16, fontweight='bold')
        fig.legend(loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)
        
        plt.tight_layout()
        plt.show()

    elif application == 'mixture_composition':
        if not peak_info:
            print("Error: 'peak_info' dictionary with compound and proton information is required.")
            return

        integrals_per_proton = {}
        for compound, info in peak_info.items():
            region = info['region']
            protons = info['protons']
            
            peak_data = df[(df['ppm'] >= region[0]) & (df['ppm'] <= region[1])]
            
            if not peak_data.empty:
                integral = simpson(peak_data['Spectrum'], x=peak_data['ppm'])
                integral_per_proton = abs(integral) / protons
                integrals_per_proton[compound] = integral_per_proton
            else:
                print(f"Warning: No data found for {compound} in the specified region {region}.")
                integrals_per_proton[compound] = 0

        # Calculate total integral per proton to find mole fractions
        total_integral_per_proton = sum(integrals_per_proton.values())
        
        if total_integral_per_proton == 0:
            print("Could not determine composition. Check peak regions and data.")
            return

        mole_fractions = {
            compound: (value / total_integral_per_proton)
            for compound, value in integrals_per_proton.items()
        }

        print("\n--- Mixture Composition ---")
        for compound, mole_fraction in mole_fractions.items():
            print(f"Mole Fraction of {compound}: {mole_fraction:.2f}")




    elif application == 'calculate_impurity':
        if not peak_info or 'main_compound' not in peak_info or 'impurity' not in peak_info:
            print("Error: 'peak_info' must contain keys for 'main_compound' and 'impurity'.")
            return
        
        main_info = peak_info['main_compound']
        impurity_info = peak_info['impurity']
        
        # Calculate integral per proton for the main compound
        main_peak_data = df[(df['ppm'] >= main_info['region'][0]) & (df['ppm'] <= main_info['region'][1])]
        if not main_peak_data.empty:
            main_integral = simpson(main_peak_data['Spectrum'], x=main_peak_data['ppm'])
            main_integral_per_proton = abs(main_integral) / main_info['protons']
        else:
            print("Error: Main compound peak not found.")
            return
        
        # Calculate integral per proton for the impurity
        impurity_peak_data = df[(df['ppm'] >= impurity_info['region'][0]) & (df['ppm'] <= impurity_info['region'][1])]
        if not impurity_peak_data.empty:
            impurity_integral = simpson(impurity_peak_data['Spectrum'], x=impurity_peak_data['ppm'])
            impurity_integral_per_proton = abs(impurity_integral) / impurity_info['protons']
        else:
            print("Error: Impurity peak not found.")
            return
            
        # Calculate total integral (on a per-proton basis)
        total_integral_per_proton = main_integral_per_proton + impurity_integral_per_proton
        
        if total_integral_per_proton == 0:
            print("Could not calculate impurity percentage. Check peak regions and data.")
            return
        
        # Calculate percentage impurity
        percent_impurity = (impurity_integral_per_proton / total_integral_per_proton) * 100
        
        print("\n--- Impurity Analysis ---")
        print(f"Main Compound Integral per Proton: {main_integral_per_proton:.4f}")
        print(f"Impurity Integral per Proton: {impurity_integral_per_proton:.4f}")
        print(f"Estimated Impurity: {percent_impurity:.2f}%")

    
    else:
        print(f"Application '{application}' is not supported.")
    




def BET_Analysis(df, application, mass_of_sample=None, cross_sectional_area=None, T=None, Pa=None, total_surface_area=None,pore_volume=None):
    """
    Perform BET (Brunauer–Emmett–Teller) analysis on adsorption data, including surface area 
    determination, pore volume, and pore radius calculations.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing adsorption data with columns:
        - 'Relative Pressure (P/P0)' : relative pressure of adsorbate
        - 'Adsorbed Volume (cm3/g STP)' : adsorbed gas volume
    application : str
        Mode of operation. Options:
        - 'plot_isotherm' : plots the adsorption isotherm.
        - 'calculate_surface_area' : plots BET plot and calculates the specific surface area.
        - 'pore_volume_calculation' : calculates the total pore volume.
        - 'pore_radius_calculations' : calculates average pore radius.
    mass_of_sample : float, optional
        Mass of the sample in grams. Required for 'calculate_surface_area'.
    cross_sectional_area : float, optional
        Cross-sectional area of adsorbate molecule (m^2). Required for 'calculate_surface_area'.
    T : float, optional
        Ambient temperature in Kelvin. Required for 'pore_volume_calculation'.
    Pa : float, optional
        Ambient pressure in Pa. Required for 'pore_volume_calculation'.
    total_surface_area : float, optional
        Total surface area (St in m^2) for pore radius calculation. Required for 'pore_radius_calculations'.
    pore_volume : float, optional
        Total pore volume (V_liq) in m^3/g. Can be used instead of recalculating from data.

    Returns
    -------
    dict or None
        Depending on the application, returns a dictionary with calculated values:
        - 'calculate_surface_area' : {'slope': m, 'intercept': b, 'v_m': vm, 'constant': c, 'sbet': SBET}
        - 'pore_volume_calculation' : {'pore_volume': V_liq}
        - 'pore_radius_calculations' : {'pore_radius_nm': r_p}
        Returns None for simple plots or if calculations fail.

    Examples
    --------
    # 1. Plot adsorption isotherm
    >>> BET_Analysis(df, application='plot_isotherm')

    # 2. Calculate BET surface area
    >>> BET_Analysis(df, application='calculate_surface_area', mass_of_sample=0.05, cross_sectional_area=0.162e-18)
    --- BET Surface Area Calculation ---
    Slope (m): 10.1234
    Y-intercept (b): 2.3456
    Monolayer Adsorbed Volume (vm): 0.1234 cm^3/g STP
    BET Constant (c): 5.32
    Specific Surface Area (SBET): 45.67 m^2/g

    # 3. Calculate pore volume
    >>> BET_Analysis(df, application='pore_volume_calculation', T=77, Pa=101325)
    --- Pore Volume Calculation ---
    Volume of gas adsorbed (V_ads): 150.0 cm^3/g STP
    Total Pore Volume (V_liq): 0.000150 m^3/g

    # 4. Calculate average pore radius
    >>> BET_Analysis(df, application='pore_radius_calculations', total_surface_area=45.67, pore_volume=0.000150)
    --- Pore Radius Calculation ---
    Total Pore Volume (V_liq): 0.000150 m^3/g
    Total Surface Area (S): 45.67 m^2
    Average Pore Radius (r_p): 6.57 nm
    """

    # Molar volume of liquid adsorbate (Vm) for N2 at 77K
    Vm = 34.65  # cm^3/mol
    
    if application == 'plot_isotherm':
        # Simple plot of Adsorption Isotherm
        plt.figure(figsize=(10, 6))
        plt.plot(df['Relative Pressure (P/P0)'], df['Adsorbed Volume (cm3/g STP)'], 'o-', color='blue')
        plt.xlabel('Relative Pressure ($P/P_0$)', fontsize=12)
        plt.ylabel('Adsorbed Volume ($cm^3/g$ STP)', fontsize=12)
        plt.title('Adsorption Isotherm', fontsize=16, fontweight='bold')
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.show()

    elif application == 'calculate_surface_area':
        if mass_of_sample is None or cross_sectional_area is None:
            print("Error: 'mass_of_sample' and 'cross_sectional_area' must be provided for this application.")
            return

        # Prepare data for the BET plot
        p = df['Relative Pressure (P/P0)']
        v = df['Adsorbed Volume (cm3/g STP)']
        
        # Calculate the y-axis values for the BET plot
        # y = 1 / (v * ((p0/p) - 1)) = (1 / (v * ((1/p) - 1)))
        bet_y = 1 / (v * ((1/p) - 1))
        
        # Plot the BET plot
        plt.figure(figsize=(10, 6))
        plt.plot(p, bet_y, 'o', color='darkgreen')
        
        # Fit a linear regression line to the data
        # We need to filter for the linear region (typically P/P0 from ~0.05 to ~0.35)
        linear_region = df[(df['Relative Pressure (P/P0)'] >= 0.05) & (df['Relative Pressure (P/P0)'] <= 0.35)]
        p_linear = linear_region['Relative Pressure (P/P0)']
        v_linear = linear_region['Adsorbed Volume (cm3/g STP)']
        
        if not p_linear.empty:
            bet_y_linear = 1 / (v_linear * ((1/p_linear) - 1))
            m, b = np.polyfit(p_linear, bet_y_linear, 1)
            
            # Plot the linear fit
            x_fit = np.linspace(p_linear.min(), p_linear.max(), 100)
            y_fit = m * x_fit + b
            plt.plot(x_fit, y_fit, '-', color='red', label=f'Linear Fit: y = {m:.2f}x + {b:.2f}')
            
            # Label the plot
            plt.xlabel('Relative Pressure ($P/P_0$)', fontsize=12)
            plt.ylabel('$\\frac{1}{v\\left(\\frac{P_0}{P}-1\\right)}$', fontsize=16)
            plt.title('BET Linear Plot', fontsize=16, fontweight='bold')
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.legend()
            plt.tight_layout()
            plt.show()

            # Calculate vm and c from the slope (m) and y-intercept (b)
            # m = (c-1) / (vm * c)
            # b = 1 / (vm * c)
            vm = 1 / (m + b)
            c = (m / b) + 1
            
            # Calculate total surface area (St)
            # St = (vm * N * s) / V  where V is molar volume at STP
            N_avogadro = 6.022e23  # molecules/mol
            V_molar_STP = 22414 # cm^3/mol
            St = (vm * N_avogadro * cross_sectional_area) / V_molar_STP
            
            # Calculate specific surface area (SBET)
            # SBET = St / a
            SBET = St / mass_of_sample
            
            print("\n--- BET Surface Area Calculation ---")
            print(f"Slope (m): {m:.4f}")
            print(f"Y-intercept (b): {b:.4f}")
            print(f"Monolayer Adsorbed Volume (vm): {vm:.4f} cm^3/g STP")
            print(f"BET Constant (c): {c:.4f}")
            print(f"Specific Surface Area (SBET): {SBET:.2f} m^2/g")
            results={'slope': m , 'intercept' : b, 'v_m': vm , 'constant' : c , 'sbet':SBET}
            return results
        else:
            print("Error: No data found in the linear region (0.05 < P/P0 < 0.35).")
            return None

    elif application == 'pore_volume_calculation':
        if T is None or Pa is None:
            print("Error: 'T' (temperature in K) and 'Pa' (ambient pressure in Pa) are required.")
            return

        # Find the adsorbed volume at the highest relative pressure (closest to unity)
        high_pressure_data = df.iloc[-1]
        V_ads = high_pressure_data['Adsorbed Volume (cm3/g STP)']
        
        # Molar volume of liquid N2 at 77K is approximately 34.65 cm^3/mol
        # Using the ideal gas law to convert V_ads (STP) to moles
        R = 8.314 # J/(mol·K)
        # Convert V_ads from cm^3 to m^3 for unit consistency (1 m^3 = 1,000,000 cm^3)
        V_ads_m3 = V_ads / 1e6
        
        # Molar volume of gas at STP (273.15K, 101325Pa)
        V_gas_STP = 22.414e-3 # m^3/mol
        n = V_ads_m3 / V_gas_STP # moles
        
        # Calculate V_liq using the formula V_liq = n * Vm
        V_liq = n * (Vm / 1e6) # m^3/g
        
        print("\n--- Pore Volume Calculation ---")
        print(f"Volume of gas adsorbed (V_ads): {V_ads:.2f} cm^3/g STP")
        print(f"Total Pore Volume (V_liq): {V_liq:.6f} m^3/g")
        results={'pore_volume': V_liq}
        return results
        
    elif application == 'pore_radius_calculations':
        if total_surface_area is None:
            print("Error: 'total_surface_area' (St) is required to calculate pore radius. Please run 'calculate_surface_area' first.")
            return

        # A more direct approach is to have the user pass the V_liq
        # For this function, let's assume V_liq is passed or calculated within.
        if 'pore_volume' not in df.columns and pore_volume is None:
            print("Error: Pore volume calculation must be performed first.")
            return

        #V_liq_final = df['pore_volume'].iloc[-1]

        V_liq_final=pore_volume 
        
        # rp = 2 * V_liq / S
        r_p = (2 * V_liq_final) / total_surface_area
        
        print("\n--- Pore Radius Calculation ---")
        print(f"Total Pore Volume (V_liq): {V_liq_final:.6f} m^3/g")
        print(f"Total Surface Area (S): {total_surface_area:.2f} m^2")
        print(f"Average Pore Radius (r_p): {r_p * 1e9:.2f} nm")

        results={'pore_radius_nm': r_p * 1e9}
        return results

    else:
        print(f"Application '{application}' is not supported.")








# ---- Reader (same as before) ----
def read_msa(filename):
    energy = []
    counts = []
    with open(filename, "r") as f:
        data_section = False
        for line in f:
            line = line.strip()
            if line.startswith("#SPECTRUM"):
                data_section = True
                continue
            if data_section and line and not line.startswith("#"):
                parts = line.replace(",", " ").split()
                if len(parts) >= 2:
                    energy.append(float(parts[0]))
                    counts.append(float(parts[1]))
    return np.array(energy), np.array(counts)


# Characteristic X-ray line energies (keV) – simplified reference values
EDS_element_lines = {
    "C": 0.277,   # Kα
    "N": 0.392,   # Kα
    "O": 0.525,   # Kα
    "F": 0.677,   # Kα
    "Na": 1.041,  # Kα
    "Mg": 1.254,  # Kα
    "Al": 1.486,  # Kα
    "Si": 1.740,  # Kα
    "P": 2.013,   # Kα
    "S": 2.307,   # Kα
    "Cl": 2.622,  # Kα
    "K": 3.312,   # Kα
    "Ca": 3.690,  # Kα
    "Ti": 4.508,  # Kα
    "V": 4.952,   # Kα
    "Cr": 5.414,  # Kα
    "Mn": 5.899,  # Kα
    "Fe": 6.404,  # Kα
    "Co": 6.930,  # Kα
    "Ni": 7.471,  # Kα
    "Cu": 8.041,  # Kα
    "Zn": 8.638,  # Kα
    "Zr": 15.775, # Kα
    "Ag": 22.162, # Kα
    "Sn": 25.271, # Kα
    "Au": 9.713,  # Mα (often observed in SEM)
    "Pb": 10.55,  # Mα
}

# Atomic weights (IUPAC 2019, rounded)
EDS_atomic_weights = {
    "C": 12.01,
    "N": 14.01,
    "O": 16.00,
    "F": 19.00,
    "Na": 22.99,
    "Mg": 24.31,
    "Al": 26.98,
    "Si": 28.09,
    "P": 30.97,
    "S": 32.06,
    "Cl": 35.45,
    "K": 39.10,
    "Ca": 40.08,
    "Ti": 47.87,
    "V": 50.94,
    "Cr": 52.00,
    "Mn": 54.94,
    "Fe": 55.85,
    "Co": 58.93,
    "Ni": 58.69,
    "Cu": 63.55,
    "Zn": 65.38,
    "Zr": 91.22,
    "Ag": 107.87,
    "Sn": 118.71,
    "Au": 196.97,
    "Pb": 207.2,
}


def EDS_Analysis(file_path, application, elements=["C","O","Fe"]):
    """
    Perform analysis on Energy Dispersive X-ray Spectroscopy (EDS) data.

    This function can plot the EDS spectrum, return raw data, quantify elemental composition,
    or detect peaks in the spectrum.

    Parameters
    ----------
    file_path : str
        Path to the EDS data file in `.msa` format.
    application : str
        Mode of operation:
        - 'plot' : Plot the EDS spectrum.
        - 'data' : Return raw energy and counts arrays.
        - 'quantify' : Estimate elemental weight and atomic percentages.
        - 'find_peak' : Detect peaks in the spectrum and plot them.
    elements : list of str, optional
        List of elements to quantify when application='quantify'. Default is ["C","O","Fe"].

    Returns
    -------
    varies
        - If application='data' : tuple (energy, counts) as numpy arrays.
        - If application='quantify' : dict with keys:
            - 'weight_percent' : dict of elements and their weight percentages
            - 'atomic_percent' : dict of elements and their atomic percentages
        - If application='find_peak' : list of tuples [(energy_keV, counts), ...] for detected peaks.
        - If application='plot' : None (displays plot only).

    Raises
    ------
    ValueError
        If the 'application' argument is not one of 'plot', 'data', 'quantify', or 'find_peak'.

    Examples
    --------
    # 1. Plot the EDS spectrum
    >>> EDS_Analysis("sample.msa", application='plot')

    # 2. Get raw energy and counts data
    >>> energy, counts = EDS_Analysis("sample.msa", application='data')

    # 3. Quantify elemental composition
    >>> results = EDS_Analysis("sample.msa", application='quantify', elements=["C","O","Fe"])
    >>> results['weight_percent']
    {'C': 12.3, 'O': 30.1, 'Fe': 57.6}
    >>> results['atomic_percent']
    {'C': 35.2, 'O': 40.8, 'Fe': 24.0}

    # 4. Find and plot peaks
    >>> peaks = EDS_Analysis("sample.msa", application='find_peak')
    >>> peaks
    [(0.28, 100), (0.53, 250), (6.40, 1200)]
    """

    energy, counts = read_msa(file_path)

    if application == 'plot':
        plt.figure(figsize=(8,5))
        plt.plot(energy, counts, color="red", linewidth=1)
        plt.xlabel("Energy (keV)")
        plt.ylabel("Counts (Intensity)")
        plt.title("EDS Spectrum")
        plt.grid(True)
        plt.show()

    elif application == 'data':
        return energy, counts

    elif application == 'quantify':
        intensities = {}
        for el in elements:
            peak_e = EDS_element_lines[el]
            mask = (energy > peak_e - 0.05) & (energy < peak_e + 0.05)
            intensities[el] =  counts[mask].sum()

        
        # --- Placeholder ZAF factors ---
        zaf_factors = {el: 1.0 for el in elements}
        corrected = {el: intensities[el] * zaf_factors[el] for el in elements}
        
        # --- Normalize to weight % ---
        total = sum(corrected.values()) if sum(corrected.values()) > 0 else 1
        weight_percent = {el: (val/total)*100 for el, val in corrected.items()}
        
        # --- Convert to atomic % ---
        mols = {el: weight_percent[el]/EDS_atomic_weights[el] for el in elements}
        total_mol = sum(mols.values()) if sum(mols.values()) > 0 else 1
        atomic_percent = {el: (val/total_mol)*100 for el, val in mols.items()}
        wt=weight_percent
        at=atomic_percent
        return {"weight_percent": wt, "atomic_percent": at}

    elif application == 'find_peak':
        # Use scipy to find peaks
        peaks, _ = find_peaks(counts, height=np.max(counts)*0.05)  # >5% of max
        peak_positions = energy[peaks]
        peak_heights = counts[peaks]

        # Plot peaks
        plt.figure(figsize=(8,5))
        plt.plot(energy, counts, color="blue")
        plt.plot(peak_positions, peak_heights, "rx", label="Detected Peaks")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Counts (Intensity)")
        plt.title("EDS Spectrum with Peaks")
        plt.legend()
        plt.grid(True)
        plt.show()

        # Return peak list
        return list(zip(peak_positions, peak_heights))

    else:
        raise ValueError("Invalid application. Use 'plot', 'data', 'quantify', or 'find_peak'.")











# Reference binding energies for elements (expand as needed)
ELEMENT_BINDING_ENERGIES = {
    'C 1s': 285.0,
    'O 1s': 532.0,
    'Ti 2p': 458.5,
    'Mo 3d': 232.0,
    'S 2p': 164.0,
}

def XPS_Analysis(df, application='plot', sensitivity_factors=None, tolerance=1.5,
                 peak_prominence=None, peak_distance=None, smoothing_window=11, smoothing_poly=3):
    """
    Perform X-ray Photoelectron Spectroscopy (XPS) data analysis.

    This function allows for plotting the XPS spectrum, returning raw data, performing
    surface composition analysis based on sensitivity factors, and detecting peaks with optional smoothing.

    Parameters
    ----------
    df : pd.DataFrame
        XPS data containing columns 'eV' (binding energy) and 'Counts / s' (intensity).
    application : str, optional
        Mode of operation (default='plot'):
        - 'plot' : Plot the XPS spectrum.
        - 'data' : Return raw energy and counts arrays.
        - 'composition' : Estimate atomic composition using peak areas and sensitivity factors.
        - 'peak_detection' : Detect peaks, optionally smooth the spectrum, and plot.
    sensitivity_factors : dict, optional
        Element-specific sensitivity factors required for 'composition' application.
        Example: {'C': 1.0, 'O': 2.93, 'Fe': 3.5}
    tolerance : float, optional
        Binding energy tolerance in eV for peak assignment (default=1.5 eV).
    peak_prominence : float, optional
        Minimum prominence of peaks for detection (used in 'composition' and 'peak_detection').
    peak_distance : int, optional
        Minimum distance between peaks in number of points (used in 'composition' and 'peak_detection').
    smoothing_window : int, optional
        Window length for Savitzky-Golay smoothing (must be odd, default=11).
    smoothing_poly : int, optional
        Polynomial order for Savitzky-Golay smoothing (default=3).

    Returns
    -------
    varies
        - If application='plot' : None (displays plot only)
        - If application='data' : tuple (energy, counts) as numpy arrays
        - If application='composition' : dict of atomic percentages {element: atomic %}
        - If application='peak_detection' : list of dicts with peak information, e.g.
          [{'energy': eV, 'counts': intensity, 'smoothed_counts': value, 
            'width': FWHM, 'start_energy': eV_start, 'end_energy': eV_end}, ...]

    Raises
    ------
    ValueError
        - If 'df' does not contain required columns
        - If 'application' is invalid
        - If sensitivity_factors are not provided for 'composition'

    Examples
    --------
    # 1. Plot XPS spectrum
    >>> XPS_Analysis(df, application='plot')

    # 2. Get raw data
    >>> energy, counts = XPS_Analysis(df, application='data')

    # 3. Compute atomic composition
    >>> sensitivity_factors = {'C': 1.0, 'O': 2.93, 'Fe': 3.5}
    >>> composition = XPS_Analysis(df, application='composition', sensitivity_factors=sensitivity_factors)
    >>> composition
    {'C': 45.3, 'O': 32.1, 'Fe': 22.6}

    # 4. Detect peaks and plot
    >>> peaks_info = XPS_Analysis(df, application='peak_detection', peak_prominence=50, smoothing_window=11)
    >>> peaks_info[0]
    {'energy': 284.8, 'counts': 1200, 'smoothed_counts': 1185, 'width': 1.2, 'start_energy': 284.0, 'end_energy': 285.6}
    """
    if 'eV' not in df.columns or 'Counts / s' not in df.columns:
        raise ValueError("DataFrame must contain 'eV' and 'Counts / s' columns")
    
    energy = df['eV'].values
    counts = df['Counts / s'].values
    
    if application == 'plot':
        plt.figure(figsize=(8,5))
        plt.plot(energy, counts, color='blue', linewidth=1)
        plt.xlabel('Binding Energy (eV)')
        plt.ylabel('Counts / s')
        plt.title('XPS Spectrum')
        plt.grid(True)
        plt.show()
        return
    
    elif application == 'data':
        return energy, counts
    
    elif application == 'composition':
        if sensitivity_factors is None:
            raise ValueError("Sensitivity factors must be provided")
        
        # Find peaks
        peaks, _ = find_peaks(counts, prominence=peak_prominence, distance=peak_distance)
        peak_energies = energy[peaks]
        peak_counts = counts[peaks]
        
        # Assign peaks to elements
        assigned_peaks = {}
        for elem, ref_energy in ELEMENT_BINDING_ENERGIES.items():
            matched = [(e, c) for e, c in zip(peak_energies, peak_counts)
                       if abs(e - ref_energy) <= tolerance]
            if matched:
                areas = []
                for e_peak, _ in matched:
                    idx = (energy >= e_peak - tolerance) & (energy <= e_peak + tolerance)
                    area = np.trapz(counts[idx], energy[idx])
                    areas.append(area)
                assigned_peaks[elem] = sum(areas) / sensitivity_factors.get(elem, 1.0)
        
        total_corrected = sum(assigned_peaks.values())
        atomic_percentages = {elem: (area/total_corrected)*100 for elem, area in assigned_peaks.items()}
        return atomic_percentages
    
    elif application == 'peak_detection':
        # Smooth spectrum
        smoothed_counts = savgol_filter(counts, window_length=smoothing_window, polyorder=smoothing_poly)
        
        # Detect peaks
        peaks, _ = find_peaks(smoothed_counts, prominence=peak_prominence, distance=peak_distance)
        widths_result = peak_widths(smoothed_counts, peaks, rel_height=0.5)
        
        # Collect peak info
        peaks_info = []
        for i, p in enumerate(peaks):
            peak_dict = {
                'energy': energy[p],
                'counts': counts[p],
                'smoothed_counts': smoothed_counts[p],
                'width': widths_result[0][i],
                'start_energy': energy[int(widths_result[2][i])],
                'end_energy': energy[int(widths_result[3][i])]
            }
            peaks_info.append(peak_dict)
        
        # Plot spectrum with detected peaks
        plt.figure(figsize=(8,5))
        plt.plot(energy, counts, color='blue', linewidth=1, label='Raw Spectrum')
        plt.plot(energy, smoothed_counts, color='orange', linewidth=1, label='Smoothed Spectrum')
        plt.plot(energy[peaks], smoothed_counts[peaks], 'rx', label='Detected Peaks')
        plt.xlabel('Binding Energy (eV)')
        plt.ylabel('Counts / s')
        plt.title('XPS Spectrum with Detected Peaks')
        plt.legend()
        plt.grid(True)
        plt.show()
        
        return peaks_info
    
    else:
        raise ValueError("Invalid application. Use 'plot', 'data', 'composition', or 'peak_detection'.")









# Physical constants
h = 6.62607015e-34  # Planck constant (J·s)
c = 3e8             # Speed of light (m/s)
e = 1.602176634e-19 # Electron charge (J/eV)

def Photoluminescence_analysis(data_frame, application="plot"):
    """
    Perform photoluminescence (PL) data analysis and visualization.

    This function analyzes a PL spectrum, identifies the main emission peak, 
    calculates bandgap energy, estimates FWHM, and provides various plots.

    Parameters
    ----------
    data_frame : pd.DataFrame
        DataFrame containing PL spectrum data with columns:
        - 'wavelength' : wavelength in nanometers (nm)
        - 'intensity' : emission intensity (arbitrary units)
    application : str, optional
        Specifies the type of analysis or visualization (default='plot'):
        - 'plot' : Plot the full PL spectrum.
        - 'peak_position' : Identify and return the wavelength of the main peak.
        - 'peak_intensity' : Identify and return the intensity of the main peak.
        - 'bandgap_energy' : Calculate bandgap energy (eV) from the peak wavelength.
        - 'fwhm' : Calculate and return the full width at half maximum (FWHM) in nm.

    Returns
    -------
    varies
        - 'plot' : None (displays a plot)
        - 'peak_position' : float, wavelength of main peak in nm
        - 'peak_intensity' : float, intensity of main peak
        - 'bandgap_energy' : float, bandgap energy in eV
        - 'fwhm' : float, full width at half maximum in nm
        - {} : empty dictionary if no peak is detected or invalid application

    Raises
    ------
    ValueError
        - If the DataFrame does not contain required columns.
        - If an invalid application string is provided.

    Notes
    -----
    - Bandgap energy is calculated using Eg = h*c / λ, where:
        h : Planck constant (J·s)
        c : speed of light (m/s)
        λ : peak wavelength (m)
        e : elementary charge (C)
    - FWHM is estimated using linear interpolation and root-finding.

    Examples
    --------
    # 1. Plot PL spectrum
    >>> Photoluminescence_analysis(df, application="plot")

    # 2. Get peak wavelength
    >>> peak_wl = Photoluminescence_analysis(df, application="peak_position")
    >>> print(f"Peak wavelength: {peak_wl:.2f} nm")

    # 3. Get peak intensity
    >>> peak_int = Photoluminescence_analysis(df, application="peak_intensity")
    >>> print(f"Peak intensity: {peak_int:.3f}")

    # 4. Calculate bandgap energy
    >>> Eg = Photoluminescence_analysis(df, application="bandgap_energy")
    >>> print(f"Bandgap: {Eg:.3f} eV")

    # 5. Calculate FWHM
    >>> fwhm = Photoluminescence_analysis(df, application="fwhm")
    >>> print(f"FWHM: {fwhm:.2f} nm")
    """
    wavelength = data_frame['wavelength'].values
    intensity = data_frame['intensity'].values

    # Find main peak
    peaks, _ = find_peaks(intensity)
    if len(peaks) == 0:
        print("No clear peak detected.")
        return {}
    main_peak_idx = peaks[np.argmax(intensity[peaks])]
    peak_wavelength = wavelength[main_peak_idx]
    peak_intensity = intensity[main_peak_idx]

    results = {
        "peak_wavelength_nm": float(peak_wavelength),
        "peak_intensity": float(peak_intensity),
        "bandgap_energy_eV": None,
        "fwhm_nm": None
    }

    plt.figure(figsize=(7,5))
    plt.plot(wavelength, intensity, label="PL Spectrum", color="blue", linewidth=2)

    if application == "plot":
        plt.title("Photoluminescence Spectrum", fontsize=14, fontweight="bold")
        plt.xlabel("Wavelength (nm)", fontsize=12)
        plt.ylabel("Intensity (a.u.)", fontsize=12)
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.tight_layout()
        plt.show()
        

    elif application == "peak_position":
        plt.axvline(peak_wavelength, color="red", linestyle="--", 
                    label=f"Peak = {peak_wavelength:.2f} nm")
        plt.title("PL Spectrum - Peak Position")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Intensity (a.u.)")
        plt.legend()
        return results['peak_wavelength_nm']

    elif application == "peak_intensity":
        plt.axhline(peak_intensity, color="green", linestyle="--", 
                    label=f"Peak Intensity = {peak_intensity:.3f}")
        plt.scatter(peak_wavelength, peak_intensity, color="red", zorder=5)
        plt.title("PL Spectrum - Peak Intensity")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Intensity (a.u.)")
        plt.legend()
        plt.tight_layout()
        plt.show()
        return results['peak_intensity']

    elif application == "bandgap_energy":
        Eg = (h * c) / (peak_wavelength * 1e-9) / e
        results["bandgap_energy_eV"] = float(Eg)
        plt.axvline(peak_wavelength, color="purple", linestyle="--", 
                    label=f"Eg ≈ {Eg:.3f} eV")
        plt.title("PL Spectrum - Bandgap Energy")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Intensity (a.u.)")
        plt.legend()
        plt.tight_layout()
        plt.show()
        return results['bandgap_energy_eV']

    elif application == "fwhm":
        half_max = peak_intensity / 2
        f_interp = interp1d(wavelength, intensity - half_max, kind='linear', fill_value="extrapolate")
        try:
            left = brentq(f_interp, wavelength[0], peak_wavelength)
            right = brentq(f_interp, peak_wavelength, wavelength[-1])
            fwhm = right - left
            results["fwhm_nm"] = float(fwhm)
            plt.axhline(half_max, color="orange", linestyle="--", 
                        label=f"FWHM = {fwhm:.2f} nm")
            plt.axvline(left, color="gray", linestyle="--")
            plt.axvline(right, color="gray", linestyle="--")
        except:
            plt.text(0.5,0.5,"FWHM not found", transform=plt.gca().transAxes)
        plt.title("PL Spectrum - FWHM")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Intensity (a.u.)")
        plt.legend()
        return results['fwhm_nm']


    else:
        print("Invalid application. Choose from: 'plot', 'peak_position', 'peak_intensity', 'bandgap_energy', 'fwhm'")
        return {}








def Dynamic_Light_Scattering_Analysis(df, application=None):
    """
    Analyze and visualize Dynamic Light Scattering (DLS) data.

    This function provides professional plotting of DLS data and
    extraction of key metrics such as the particle size corresponding
    to the maximum intensity.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing DLS data. Expected columns include:
        - 'Size (nm)' : Particle size in nanometers
        - 'Intensity (%)' : Corresponding intensity in percentage
        - 'Lag time (µs)' : Lag time for autocorrelation measurements
        - 'Autocorrelation' : Autocorrelation function values
    application : str, optional
        Type of analysis to perform:
        - 'plot' : Generate professional plots based on available columns.
            - If 'Size (nm)' and 'Intensity (%)' exist, plots Intensity vs Size.
            - If 'Lag time (µs)' and 'Autocorrelation' exist, plots Autocorrelation vs Lag time.
        - 'max_intensity' : Returns the particle size corresponding to maximum intensity.

    Returns
    -------
    dict or None
        - If `application='max_intensity'`:
            Dictionary with keys:
            - "Peak Size (nm)" : particle size at maximum intensity
            - "Peak Intensity (%)" : intensity at that size
        - If `application='plot'` or None, returns None and displays plots.

    Raises
    ------
    ValueError
        - If required columns are missing for the selected application.
        - If `application` is invalid (not 'plot' or 'max_intensity').

    Examples
    --------
    # 1. Plot DLS Intensity vs Size
    >>> Dynamic_Light_Scattering_Analysis(df, application="plot")

    # 2. Plot Autocorrelation vs Lag time
    >>> Dynamic_Light_Scattering_Analysis(df_with_autocorr, application="plot")

    # 3. Get particle size at maximum intensity
    >>> result = Dynamic_Light_Scattering_Analysis(df, application="max_intensity")
    >>> print(result)
    {'Peak Size (nm)': 120.5, 'Peak Intensity (%)': 85.2}
    """

    if application == "plot":
        
        # Set professional style
        sns.set(style="whitegrid", context="talk", palette="deep")

        # Plot Intensity vs Size
        if 'Intensity (%)' in df.columns and 'Size (nm)' in df.columns:
            plt.figure(figsize=(10, 6))
            plt.plot(df['Size (nm)'], df['Intensity (%)'], marker='o', linestyle='-', linewidth=2, markersize=6, color='tab:blue')
            plt.xlabel('Size (nm)', fontsize=14)
            plt.ylabel('Intensity (%)', fontsize=14)
            plt.title('Dynamic Light Scattering: Intensity vs Size', fontsize=16, fontweight='bold')
            plt.grid(True, which='both', linestyle='--', linewidth=0.5)
            plt.minorticks_on()
            plt.tight_layout()
            plt.show()

        # Plot Autocorrelation vs Lag time
        elif 'Lag time (µs)' in df.columns and 'Autocorrelation' in df.columns:
            plt.figure(figsize=(10, 6))
            plt.plot(df['Lag time (µs)'], df['Autocorrelation'], marker='s', linestyle='-', linewidth=2, markersize=6, color='tab:orange')
            plt.xlabel('Lag time (µs)', fontsize=14)
            plt.ylabel('Autocorrelation', fontsize=14)
            plt.title('Dynamic Light Scattering: Autocorrelation vs Lag time', fontsize=16, fontweight='bold')
            plt.grid(True, which='both', linestyle='--', linewidth=0.5)
            plt.minorticks_on()
            plt.tight_layout()
            plt.show()

        else:
            print("DataFrame columns not recognized for plotting.")

    elif application == "max_intensity":
        if 'Intensity (%)' in df.columns and 'Size (nm)' in df.columns:
            max_idx = df['Intensity (%)'].idxmax()
            peak_size = df.loc[max_idx, 'Size (nm)']
            peak_intensity = df.loc[max_idx, 'Intensity (%)']
            return {"Peak Size (nm)": peak_size, "Peak Intensity (%)": peak_intensity}
        else:
            print("DataFrame must contain 'Size (nm)' and 'Intensity (%)' columns for max_intensity.")

    else:
        print("Invalid application. Use 'plot' or 'max_intensity'.")







def Auger_Electron_Spectroscopy_analysis(df, application=None, sensitivity_factors=None):
    """
    Analyze and visualize Auger Electron Spectroscopy (AES) data.

    This function provides options to plot AES spectra, detect peak positions, 
    and estimate atomic percentages using sensitivity factors.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing AES data. Must include columns:
        - 'Energy (eV)' : Electron energy values in eV
        - 'Intensity (Counts)' : Corresponding measured intensity
    application : str, optional
        Type of analysis to perform:
        - 'plot' : Generates a professional plot of Intensity vs Energy.
        - 'peak_position' : Detects peaks and returns their energy positions and intensities.
        - 'atomic' : Calculates atomic percentages based on provided sensitivity factors.
    sensitivity_factors : dict, optional
        Dictionary mapping element symbols to their sensitivity factors.
        Example: {'C': 0.25, 'O': 0.66, 'Fe': 2.5}.
        Required if `application='atomic'`.

    Returns
    -------
    dict or list or None
        - If `application='plot'` : None (displays plot)
        - If `application='peak_position'` : dict with keys:
            - "Peak Positions (eV)" : numpy array of peak energies
            - "Peak Intensities (Counts)" : numpy array of peak intensities
        - If `application='atomic'` : list of dicts for each element, e.g.:
            [{"Element": "C", "Atomic %": 25.4}, {"Element": "O", "Atomic %": 74.6}]

    Raises
    ------
    ValueError
        If `application='atomic'` and `sensitivity_factors` is not provided.

    Examples
    --------
    # 1. Plot AES spectrum
    >>> Auger_Electron_Spectroscopy_analysis(df, application='plot')

    # 2. Detect peak positions
    >>> peaks = Auger_Electron_Spectroscopy_analysis(df, application='peak_position')
    >>> print(peaks)
    {'Peak Positions (eV)': array([280, 530]), 'Peak Intensities (Counts)': array([150, 200])}

    # 3. Estimate atomic composition
    >>> sensitivity = {'C': 0.25, 'O': 0.66, 'Fe': 2.5}
    >>> composition = Auger_Electron_Spectroscopy_analysis(df, application='atomic', sensitivity_factors=sensitivity)
    >>> print(composition)
    [{'Element': 'C', 'Atomic %': 30.5}, {'Element': 'O', 'Atomic %': 69.5}]
    """
    if application == "plot":
        sns.set(style="whitegrid", context="talk", palette="deep")
        plt.figure(figsize=(10, 6))
        plt.plot(df['Energy (eV)'], df['Intensity (Counts)'], color="tab:blue", linewidth=1.8)
        plt.xlabel("Energy (eV)", fontsize=14)
        plt.ylabel("Intensity (Counts)", fontsize=14)
        plt.title("AES Spectrum: Energy vs Intensity", fontsize=16, fontweight="bold")
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.tight_layout()
        plt.show()

    elif application == "peak_position":
        # Detect peaks
        peaks, _ = find_peaks(df['Intensity (Counts)'], height=np.mean(df['Intensity (Counts)']))
        peak_positions = df['Energy (eV)'].iloc[peaks].values
        peak_intensities = df['Intensity (Counts)'].iloc[peaks].values
        return {
            "Peak Positions (eV)": peak_positions,
            "Peak Intensities (Counts)": peak_intensities
        }

    elif application == "atomic":
        if sensitivity_factors is None:
            raise ValueError("Please provide sensitivity_factors dictionary for atomic analysis.")

        # Detect peaks (areas approximation by trapezoidal integration around peaks)
        peaks, _ = find_peaks(df['Intensity (Counts)'], height=np.mean(df['Intensity (Counts)']))
        peak_positions = df['Energy (eV)'].iloc[peaks].values
        peak_intensities = df['Intensity (Counts)'].iloc[peaks].values

        results = []
        corrected_areas = {}
        total = 0.0

        # Loop through elements in sensitivity_factors and assign peaks (user maps manually in real AES)
        for element, S in sensitivity_factors.items():
            # Approximate peak area (here we just use intensity, could integrate around peak in real analysis)
            if len(peak_intensities) > 0:
                A = max(peak_intensities)  # simplified peak area ~ max intensity
                corrected_area = A / S
                corrected_areas[element] = corrected_area
                total += corrected_area

        # Calculate atomic percentages
        for element, corrected_area in corrected_areas.items():
            Ci = (corrected_area / total) * 100 if total > 0 else 0
            results.append({"Element": element, "Atomic %": Ci})

        return results

    else:
        print("Invalid application. Use 'plot', 'peak_position', or 'atomic'.")





#================================================================
#================================================================
#================================================================
#================================================================
#================================================================
#================================================================
#================================================================
#================================================================





def Tensile_Analysis(dataframe, gauge_length=1, width=1, thickness=1,
                           application='plot-force', save=False,):
    """
    Parameters:
    - dataframe: raw data from Excel (Force vs Displacement)
    - gauge_length: Initial length of the sample in mm
    - width: Width of the sample in mm
    - thickness: Thickness of the sample in mm
    - application: 'plot-force' or 'plot-stress'
    - save: True to save the plot
    - show_peaks: True to annotate peaks (e.g. UTS)
    - fname: Filename to save if save=True
    """
    dataframe.drop(labels='1 _ 1',axis=1,inplace=True)
    dataframe.drop(labels='Unnamed: 3',axis=1,inplace=True)
    dataframe.drop(index=0,inplace=True)
    dataframe.drop(index=1,inplace=True)
    dataframe.reset_index(inplace=True,drop=True)


    d2=np.array(dataframe)
    d2=d2.astype(float)
    
    force = d2[:, 0]       # in N
    displacement = d2[:, 1]  # in mm

    # Cross-sectional Area (mm²)
    area = width * thickness

    # Compute strain and stress
    strain = displacement / int(gauge_length)
    stress = force / int(area)  # in MPa

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.rcParams["figure.dpi"] = 600

    if application == 'plot-force':
        plt.plot(displacement, force, label='Force vs Displacement', c='blue')
        plt.xlabel('Displacement (mm)')
        plt.ylabel('Force (N)')
        plt.title('Tensile Test - Force vs Displacement', size=16)
        plt.show()
        if save==True:
        
            plt.savefig('stress_strain',dpi=600,format='eps')
        
        
    elif application == 'plot-stress':
        plt.plot(strain, stress, label='Stress vs Strain', c='green')
        plt.xlabel('Strain')
        plt.ylabel('Stress (MPa)')
        plt.title('Tensile Test - Stress vs Strain', size=16)
        plt.show()
        if save==True:
        
            plt.savefig('stress_strain',dpi=600,format='eps')


    elif application == 'UTS' :
        uts = np.max(stress)
        print(f"Ultimate Tensile Strength (UTS): {uts:.2f} MPa")
        return uts
        
    elif application == 'Young Modulus':
        linear_region = int(len(strain) * 0.1)
        E = np.polyfit(strain[:linear_region], stress[:linear_region], 1)[0]

        print(f" Young’s Modulus (E): {E:.2f} MPa")
        return E

        
    elif application == 'Fracture Stress':
        
        stress_at_break = stress[-1]
        print(f"Fracture Stress: {stress_at_break:.2f} MPa")

        return stress_at_break
        
    elif application == 'Strain at break':
       
        strain_at_break = strain[-1]
        print(f"Strain at Break: {strain_at_break:.4f}")
        return strain_at_break
    
    


def FtirAnalysis(dataframe, application, prominence=0.5, distance=10, save=False):
    """
    Parameters:
    - dataframe: pandas.DataFrame
        Raw FTIR data (expects one column with tab-separated values 'X Y').
    - application: str
        One of ['plot', 'peak'].
        'plot' will generate an FTIR plot.
        'peak' will detect and return peak positions and properties.
    - prominence: float, default=0.5
        Required prominence of peaks (used in peak detection).
    - distance: int, default=10
        Minimum horizontal distance (in number of samples) between peaks.
    - save: bool, default=False
        If True, save the generated plot.
    """

    xValues = []
    yValues = []
    for i in range(len(dataframe)):
        x, y = dataframe['X\tY'][i].split()
        xValues.append(float(x))
        yValues.append(float(y))

    if application == 'plot':
        plt.figure(figsize=(10, 6))
        plt.rcParams["figure.dpi"] = 600
        plt.plot(xValues, yValues, c='k')
        plt.title('FTIR Result', size=20)
        plt.xlabel('Wavenumber (cm⁻¹)')
        plt.ylabel('Transmittance (a.u.)')
        plt.xlim(4000, 400)
        plt.ylim(28, 40)

        ax = plt.gca()
        ax.tick_params(axis='both', which='both', direction='in', bottom=True, top=False, left=True, right=False)
        plt.show()

        if save:
            plt.savefig('ftir', dpi=600, format='eps')

    elif application == 'peak':
        peaks, properties = find_peaks(yValues, prominence=prominence, distance=distance)
        return peaks, properties



   
    
def FTIR(data1,application,prominence=0.5, distance=10,save=False):
    '''
    OLD Version V1.00


    '''
    
    '''
    xx=[]
    yy=[]
    for i in range(0,len(data1)):
        b=data1['X\tY'][i].split()
        xx.append(float(b[0]))
        yy.append(float(b[1]))
        
        
    if application=='plot':
    
        plt.figure(figsize=(10, 6))  # Set the figure size (width: 10 inches, height: 6 inches)
        plt.rcParams["figure.dpi"] = 600 
        plt.plot(xx,yy,c='k')
        plt.title('FTIR Result',size=20)
        plt.xlabel('Wavenumber (cm-1)')
        plt.ylabel('Transmitance(a.u.')
        plt.xlim(4000,400)
        plt.ylim(28,40)
        #plt.invert_xaxis()
        ax=plt.gca()
        ax.tick_params(axis='both',which='both',direction='in',bottom=True,top=False,left=True,right=False)
        plt.show()
        
        if save==True:
            plt.savefig('ftir',dpi=600,format='eps')
    
    elif application=='peak':
        peaks, properties = find_peaks(yy, prominence=prominence, distance=distance)
        
        return peaks, properties
    '''

    print('Warning: This function is deprecated. Please use "FtirAnalysis" instead.')
    return None
        
    


def XrdZno(dataframe, application):
    """
    Parameters:
    - dataframe: pandas.DataFrame
        Data containing XRD data. Expected columns: ['Angle', 'Det1Disc1'].
    - application: str
        One of ['plot', 'FWHM', 'Scherrer'].
        'plot'      → Draw the XRD pattern.
        'FWHM'      → Calculate Full Width at Half Maximum.
        'Scherrer'  → Calculate crystallite size using Scherrer equation.

    Returns:
    - float or None
        Returns FWHM (float) if application='FWHM'.
        Returns crystallite size (float) if application='Scherrer'.
        Returns None if application='plot'.
    """
    angles = np.array(dataframe['Angle'])
    intensities = np.array(dataframe['Det1Disc1'])

    if application == 'plot':
        plt.plot(angles, intensities, c='red')
        plt.title('XRD Pattern')
        plt.xlabel('2θ (degrees)')
        plt.ylabel('Intensity')
        plt.show()

    elif application in ['FWHM', 'Scherrer']:
        maxIntensity = np.max(intensities)
        halfMax = maxIntensity / 2

        indices = [i for i, val in enumerate(intensities) if val >= halfMax]

        if len(indices) > 0:
            leftIndex = np.min(indices)
            rightIndex = np.max(indices)
            fwhm = angles[rightIndex] - angles[leftIndex]

            if application == 'FWHM':
                return fwhm

            elif application == 'Scherrer':
                mean2theta = angles[indices].mean()
                theta = mean2theta / 2

                fwhmRad = np.deg2rad(fwhm)
                thetaRad = np.deg2rad(theta)

                crystalSize = (0.9 * 1.5406) / (fwhmRad * np.cos(thetaRad))
                return crystalSize


    
def XRD_ZnO(XRD,application):
    '''
    

    Parameters
    ----------
    XRD : DataFrame
        Data containing XRD data.
    application : str
        Type of application 'plot','FWHM','Scherrer'.
        plot:To draw the figure.
        FWHM:Width at Half Maximum.
        Scherrer:To calculate the crystallite size.

    Returns
    FWHM,Scherrer
    -------
    None.

    '''
    '''
    Angles=np.array(XRD['Angle'])
    Intensities=np.array(XRD['Det1Disc1'])
    if  application=='plot':
        plt.plot(Angles,Intensities,c='red')
        plt.title('XRD Pattern')
        plt.xlabel('2theta (degrees)')
        plt.ylabel('Intensity')
        plt.show()
    elif application in ['FWHM', 'Scherrer']:
        max_intensity = np.max(Intensities)
        half_max = max_intensity / 2
        indices = []
        half_max = max_intensity / 2
        for i in range(len(Intensities)):
           if Intensities[i] >= half_max:
               indices.append(i)

        
        if len(indices) > 0:
            left_index = np.min(indices)
            right_index = np.max(indices)
       
    
            FWHM = Angles[right_index] - Angles[left_index]
            if application == 'FWHM':
                return FWHM
           
            elif application =='Scherrer':
                mean_2theta = Angles[indices].mean()


                theta = mean_2theta / 2
                FWHM_rad = ((3.14/180)*FWHM)
                theta_rad = ((3.14/180)*theta)  
                crystal_size = (0.9 * 1.5406) / (FWHM_rad * np.cos(theta_rad))
                
                return crystal_size
    '''
    print('Warning: This function is deprecated. Please use "XrdZno" instead.')
    return None


    
    
def PressureVolumeIdealGases(dataframe, application):
    """
    Parameters:
    - dataframe: pandas.DataFrame
        Must contain 'pressure' and 'volume' columns.
    - application: str
        One of ['plot', 'min pressure', 'max pressure', 'min volume', 
                'max volume', 'average pressure', 'average volume', 'temperature'].

    Returns:
    - float, pandas.Series, or None
        Depending on the selected application.
    """
    if application == 'plot':
        pressure = dataframe['pressure']
        volume = dataframe['volume']
        plt.plot(volume, pressure)
        plt.title('Volume-Pressure Chart')
        plt.xlabel('Volume')
        plt.ylabel('Pressure')
        plt.show()

    elif application == 'min pressure':
        return dataframe['pressure'].min()

    elif application == 'max pressure':
        return dataframe['pressure'].max()

    elif application == 'min volume':
        return dataframe['volume'].min()

    elif application == 'max volume':
        return dataframe['volume'].max()

    elif application == 'average pressure':
        return dataframe['pressure'].mean()

    elif application == 'average volume':
        return dataframe['volume'].mean()

    elif application == 'temperature':
        n = 1
        R = 0.821
        return (dataframe['pressure'] * dataframe['volume']) / (n * R)

    else:
        print("Invalid application selected.")


    
    
    
    
def EnergieAnalysis(dataframe, application):
    """
    Parameters:
    - dataframe: pandas.DataFrame
        Must contain motor energy data with columns 
        ['Angle[°]', 'Energie', 'Power[mW]', 'Time for a Cycle'].
    - application: str
        One of ['draw', 'calculate'].
        'draw'      → Plot energy vs angle.
        'calculate' → Calculate total consumption energy in Ws.

    Returns:
    - float or None
        Energy consumption in Ws if application='calculate'.
        None if application='draw'.
    """
    if application == 'plot':
        angle = dataframe['Angle[°]']
        energy = dataframe['Energie']
        plt.plot(angle, energy, color='green')
        plt.title('Energy of OTT Motor (185000 Cycles)')
        plt.xlabel('Angle [°]')
        plt.ylabel('Consumption Energy')
        plt.show()

    elif application == 'calculate':
        dataframe = dataframe[['Angle[°]', 'Power[mW]', 'Time for a Cycle', 'Energie']]
        summ = dataframe['Energie'].sum()
        summ = (summ * 2) / 1000  # Convert to Ws
        return summ







def Stress_Strain1(df,operation,L0=90,D0=9):
    '''
    
    
    This function gets data and an operation .
    It plots Stress-Strain curve if the oepration is plot 
    and finds the UTS value (which is the ultimate tensile strength) otherwise.
    ------------------------------
    Parameters
    ----------
    df : DataFrame
       It has 2 columns: DL(which is length in mm) & F (which is the force in N).
    operation :
       It tells the function to whether PLOT the curve or find the UTS valu. 
       
     L0: initial length of the sample
     D0: initial diameter of the sample

    Returns
    -------
    The Stress-Strain curve or the amount of UTS
    
    '''
    
    A0 = math.pi / 4 * (D0 ** 2)
    df['e'] = df['DL'] / L0
    df['S'] = df['F'] / A0
    if operation == 'PLOT':
        plt.scatter(df['e'], df['S'])
        plt.xlabel('e')
        plt.ylabel('S')
        plt.title('S vs e Plot')
        plt.grid(True)
        plt.show()
    elif operation == 'UTS':
        return df['S'].max()
    else:
        print("Please enter proper operation")
        return







def Stress_Strain2(input_file,which,count):
    '''
    This function claculates the stress and strain
    Parameters from load and elongation data
    ----------
    input_file : .csv format
        the file must be inserted in csv.
    whcih : str
        please say which work we do ( plot or calculate?).
    count: int
        please enter the yarn count in Tex
    remember: gauge length has been set in 250 mm
    '''

    #convert the file
    mydf=pd.read_csv(input_file)

    if which=='plot':
       
        stress=mydf['Load']/count
        strain=mydf['Extension']/250
        plt.plot(stress,strain)
        plt.title('stress-strain curve')
        plt.xlabel('stress')
        plt.ylabel('strain')
        plt.show()
    
    
    if which=='max stress':
        stress_max=mydf['stress'].max()
        return stress_max

    if which=='max strain':
        strain_max=mydf['strain'].max()
        return strain_max
    
    
def Stress_Strain3(input_data, action):
    stress = input_data['Stress (MPa)']
    strain = input_data['Strain (%)']
    
    if action == 'plot':
        # Plotting data
        fig, ax = plt.subplots(figsize=(8, 6))
        plt.plot(strain, stress, linewidth=2, color='royalblue', marker='o', markersize=5, label='Stress-Strain Curve')
        plt.title('Stress-Strain Curve', fontsize=16)
        plt.xlabel('Strain (%)', fontsize=14)
        plt.ylabel('Stress (MPa)', fontsize=14)
        plt.xlim([0, strain.max()])
        plt.ylim([0, stress.max()])
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend()
   
    elif action == 'max stress':
        # Calculation of the maximum stress
        stress_max = stress.max()
        return stress_max
    
    elif action == 'young modulus':
        # Calculation of Young's Modulus
        slope_intercept = np.polyfit(strain, stress, 1)
        return slope_intercept[0]

def Stress_Strain4(file_path, D0, L0):
    '''
    This function uses the data file
    that contains length and force, calculates the engineering, true
    and yielding stress and strain and also draws a graph of these.
    
    Parameters:
    D0(mm): First Qatar to calculate stress
    L0(mm): First Length to canculate strain
    F(N): The force applied to the object during this test
    DL(mm): Length changes
    
    Returns:
    Depending on the operation selected,
    it returns calculated values, plots,
    advanced analysis, or saves results.
    '''
    try:
        data = pd.read_excel(file_path)
        
    except FileNotFoundError:
        print("File not found. Please check the file path.")
        return

    A0 = math.pi * (D0/2)**2

    data['stress'] = data['F (N)'] / A0
    data['strain'] = (data['DL (mm)'] - L0) / L0

    data['true_stress'] = data['F (N)'] / A0
    data['true_strain'] = np.log(1 + data['strain'])

    yield_point = data.loc[data['stress'].idxmax()]
    permanent_strain = data['strain'].iloc[-1]

    plt.figure(figsize=(12, 8))
    plt.plot(data['strain'], data['stress'], label='Engineering Stress-Strain', marker='o', color='b', linestyle='-')
    plt.plot(data['true_strain'], data['true_stress'], label='True Stress-Strain', marker='x', color='r', linestyle='--')
    plt.scatter(yield_point['strain'], yield_point['stress'], color='g', label='Yield Point')
    plt.annotate(f"Yield Point: Strain={yield_point['strain']:.2f}, Stress={yield_point['stress']:.2f}", (yield_point['strain'], yield_point['stress']), textcoords="offset points", xytext=(0,10), ha='center')
    plt.xlabel('Strain')
    plt.ylabel('Stress (MPa)')
    plt.title('Stress-Strain Curve')
    plt.grid(True)
    plt.legend()
    plt.show()

    print("Columns in the data:")
    print(data.columns)
    
    print("\nFirst few rows of the data:")
    print(data.head())

    print("\nYield Point Information:")
    print(yield_point)
    print("Permanent Strain:", permanent_strain)




def Stress_Strain5(input_data, action):
    stress = input_data['Stress (MPa)']
    strain = input_data['Strain (%)']
    
    if action == 'plot':
        # Plotting data
        fig, ax = plt.subplots(figsize=(8, 6))
        plt.plot(strain, stress, linewidth=2, color='royalblue', marker='o', markersize=5, label='Stress-Strain Curve')
        plt.title('Stress-Strain Curve', fontsize=16)
        plt.xlabel('Strain (%)', fontsize=14)
        plt.ylabel('Stress (MPa)', fontsize=14)
        plt.xlim([0, strain.max()])
        plt.ylim([0, stress.max()])
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend()
   
    elif action == 'max stress':
        # Calculation of the maximum stress
        stress_max = stress.max()
        return stress_max
    
    elif action == 'young modulus':
        # Calculation of Young's Modulus
        slope_intercept = np.polyfit(strain, stress, 1)
        return slope_intercept[0]





def Stress_Strain6(data,application):
    '''
    this function converts F and dD to Stress and Strain by thickness(1.55mm), width(3.2mm) and parallel length(35mm).

    Parameters
    ----------
    data : DataFrame
        this DataFrame contains F(N) and dD(mm) received from the tensil test machine.
    application : str
        application determines the expected output of Stress_Strain function.

    Returns
    -------
    int, float or plot
        return may be elongation at break, strength or a plot.

    '''
    
    stress=np.array([data['F']/(1.55*3.2)])
    strain=np.array([(data['dD']/35)*100])
    if application.upper()=='ELONGATION AT BREAK':
        elongation_at_break=np.max(strain)
        print(elongation_at_break,'%')
        return elongation_at_break
    elif application.upper()=='STRENGTH':
        strength=np.max(stress)
        print(strength,'N/mm2')
        return strength
    elif application.upper()=='PLOT':
        myfont_title={'family':'sans-serif',
                      'color':'black',
                      'size':20}
        myfont_lables={'family':'Tahoma',
                       'color':'green',
                       'size':16}
        plt.plot(strain,stress,ls='--',c='g',linewidth=10)
        plt.title('Stress-Strain',fontdict=myfont_title)
        plt.xlabel('Strain(%)',fontdict=myfont_lables)
        plt.ylabel('Stress(N/mm2)',fontdict=myfont_lables)
        plt.show()



def AerospaceAnalysis(dataframe, application):
    """

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Must contain two columns: ['Newton', 'Area'].
        Values should be in Newtons (N) and square meters (m²).
    application : str
        One of ['plot', 'maxPressure'].
        'plot'        → Plot Newton vs Area.
        'maxPressure' → Return maximum pressure value.

    Returns:
    -------
    float or None
        Maximum pressure if application='maxPressure'.
        None if application='plot'.
    """
    # Ensure proper DataFrame format
    df = dataframe.copy()
    df['Pressure'] = df['Newton'] / df['Area']

    if application == 'plot':
        plt.plot(df['Area'], df['Newton'])
        plt.xlabel('Area (m²)')
        plt.ylabel('Force (N)')
        plt.title('Force vs Area')
        plt.show()

    elif application == 'maxPressure':
        return df['Pressure'].max()



def XRD_Analysis(file,which,peak=0):
    '''
    

    Parameters
    ----------
    file : str
        the variable in which you saved the .cvs file path         
    which : str
        which operation you want to perform on the file      
    peak : float, optional
        2θ for the peak you want to analyse. The default is 0.     

    Returns
    -------
    fwhm : float
        value of FWHM for the peak you specified.

    '''
    '''
    
    df=pd.read_csv(file)
    npar=pd.DataFrame.to_numpy(df)

    if which=='plot':
        angle=df['angle']
        intensity=df['intensity']
        plt.plot(angle,intensity,color='k')
        font_title={'family':'serif','color':'blue','size':20}
        plt.title('XRD pattern',fontdict=font_title)
        font_label={'family':'times new roman','color':'black','size':15}
        plt.xlabel('angle (2θ)',fontdict=font_label)
        plt.ylabel('intensity (a.u.)',fontdict=font_label)
        plt.grid(axis='x',which='both')
        plt.xticks(np.arange(0,max(angle),5))
        plt.xlim([np.min(npar,axis=0)[0], np.max(npar,axis=0)[0]])
        plt.yticks([])
        plt.ylim([0, 1.1*np.max(npar,axis=0)[1]])
        plt.tick_params(axis='x',direction='in')
        plt.show()
        return None
    elif which=='fwhm':
        diff=int((npar[1,0]-npar[0,0])*1000)/2000
        for i in range(int(len(npar)/2)+1):
            if -diff<npar[i,0]-peak<diff:
                pl=i
                ph=i
                p=i
                break
        while pl>0:
            if ((npar[pl,1]-npar[pl-1,1])/(npar[pl-1,1]-npar[pl-2,1]))>1.04 and (npar[pl-1,1]-np.min(npar,axis=0)[1])/(np.max(npar,axis=0)[1]-np.min(npar,axis=0)[1])<0.4:
                in_low_1=npar[pl-1,1]
                break
            pl=pl-1
        while ph>0:
            if ((npar[ph+2,1]-npar[ph+1,1])/(npar[ph+1,1]-npar[ph,1]))<0.96 and (npar[ph+1,1]-np.min(npar,axis=0)[1])/(np.max(npar,axis=0)[1]-np.min(npar,axis=0)[1])<0.4:
                in_low_2=npar[ph+1,1]
                break
            ph=ph+1
        in_low=(in_low_1+in_low_2)/2
        h=npar[p,1]-in_low
        hm=in_low+h/2
        diff_in=[]
        hm_i=[]
        for l in range(len(npar)-1):
            diff_in.append((npar[l+1,1]-npar[l,1])/2)
        for j in range(2):
            for k in range(int(len(npar)/2)+1):
                c=((-1)**j)*k
                if abs(npar[p+c,1]-hm)<abs(max(diff_in)):
                    hm_i.append(p+c)
                    break
        fwhm=npar[hm_i[0],0]-npar[hm_i[1],0]
        return fwhm
    else:
        print('The which argument not valid')
        return None
    '''
    print('Warning: This function is deprecated. Please use "XrdAnalysis" instead.')
    return None




def XrdAnalysis(df: pd.DataFrame, which: str, peak: float = 0):
    """
    Perform XRD (X-ray Diffraction) analysis on a given DataFrame containing 'angle' and 'intensity'.

    Parameters
    ----------
    df : pd.DataFrame
        A pandas DataFrame with at least two columns: 'angle' and 'intensity'.
    which : str
        Operation to perform on the DataFrame. Options:
        - 'plot' : Plots the XRD pattern.
        - 'fwhm' : Calculates the Full Width at Half Maximum (FWHM) for a given peak.
    peak : float, optional
        The 2θ angle of the peak to analyze. Default is 0.

    Returns
    -------
    fwhm : float or None
        - If which == 'fwhm', returns the FWHM value of the specified peak.
        - If which == 'plot', returns None.

    Example
    -------
    >>> data = pd.DataFrame({'angle': [20, 21, 22, 23, 24],
    ...                      'intensity': [5, 20, 50, 20, 5]})
    >>> xrdAnalysis(data, which='plot')   # Plots the XRD pattern
    >>> xrdAnalysis(data, which='fwhm', peak=22)
    0.5
    """
    npar = df.to_numpy()

    if which == 'plot':
        angle = df['angle']
        intensity = df['intensity']
        plt.plot(angle, intensity, color='k')
        font_title = {'family': 'serif', 'color': 'blue', 'size': 20}
        plt.title('XRD pattern', fontdict=font_title)
        font_label = {'family': 'times new roman', 'color': 'black', 'size': 15}
        plt.xlabel('angle (2θ)', fontdict=font_label)
        plt.ylabel('intensity (a.u.)', fontdict=font_label)
        plt.grid(axis='x', which='both')
        plt.xticks(np.arange(0, max(angle), 5))
        plt.xlim([np.min(npar, axis=0)[0], np.max(npar, axis=0)[0]])
        plt.yticks([])
        plt.ylim([0, 1.1 * np.max(npar, axis=0)[1]])
        plt.tick_params(axis='x', direction='in')
        plt.show()
        return None

    elif which == 'fwhm':
        diff = int((npar[1, 0] - npar[0, 0]) * 1000) / 2000
        for i in range(int(len(npar) / 2) + 1):
            if -diff < npar[i, 0] - peak < diff:
                pl = i
                ph = i
                p = i
                break
        while pl > 0:
            if ((npar[pl, 1] - npar[pl - 1, 1]) / (npar[pl - 1, 1] - npar[pl - 2, 1])) > 1.04 and \
               (npar[pl - 1, 1] - np.min(npar, axis=0)[1]) / (np.max(npar, axis=0)[1] - np.min(npar, axis=0)[1]) < 0.4:
                in_low_1 = npar[pl - 1, 1]
                break
            pl -= 1
        while ph > 0:
            if ((npar[ph + 2, 1] - npar[ph + 1, 1]) / (npar[ph + 1, 1] - npar[ph, 1])) < 0.96 and \
               (npar[ph + 1, 1] - np.min(npar, axis=0)[1]) / (np.max(npar, axis=0)[1] - np.min(npar, axis=0)[1]) < 0.4:
                in_low_2 = npar[ph + 1, 1]
                break
            ph += 1
        in_low = (in_low_1 + in_low_2) / 2
        h = npar[p, 1] - in_low
        hm = in_low + h / 2
        diff_in = []
        hm_i = []
        for l in range(len(npar) - 1):
            diff_in.append((npar[l + 1, 1] - npar[l, 1]) / 2)
        for j in range(2):
            for k in range(int(len(npar) / 2) + 1):
                c = ((-1) ** j) * k
                if abs(npar[p + c, 1] - hm) < abs(max(diff_in)):
                    hm_i.append(p + c)
                    break
        fwhm = npar[hm_i[0], 0] - npar[hm_i[1], 0]
        return fwhm

    else:
        raise ValueError("Invalid argument for 'which'. Use 'plot' or 'fwhm'.")


def LN_S_E(df, operation):
    """
    This function analyzes the elastic part of a true stress-strain curve.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain 2 columns:
        - 'DL' : elongation (length change in mm)
        - 'F'  : force in Newtons

    operation : str
        - 'PLOT' : plots the elastic region of the true stress-strain curve
        - 'YOUNG_MODULUS' : calculates and returns Young's Modulus (E)

    Returns
    -------
    None if operation='PLOT'
    float if operation='YOUNG_MODULUS'
    """

    # ---- Material Test Constants ----
    L0 = 40  # initial gauge length (mm)
    D0 = 9   # initial diameter (mm)
    A0 = math.pi / 4 * (D0 ** 2)  # initial cross-sectional area (mm^2)

    # ---- Engineering Strain & Stress ----
    df['e'] = df['DL'] / L0
    df['S'] = df['F'] / A0

    # ---- True Strain & True Stress ----
    df['eps'] = np.log(1 + df['e'])      # true strain
    df['sig'] = df['S'] * (1 + df['e'])  # true stress

    # ---- Select Elastic Region (strain 0.04–0.08) ----
    mask = (df['eps'] >= 0.04) & (df['eps'] <= 0.08)
    elastic_df = df.loc[mask, ['eps', 'sig']].copy()

    # ---- Log Transform ----
    elastic_df['ln_eps'] = np.log(elastic_df['eps'])
    elastic_df['ln_sig'] = np.log(elastic_df['sig'])

    # ---- Perform Requested Operation ----
    if operation.upper() == 'PLOT':
        plt.scatter(elastic_df['ln_eps'], elastic_df['ln_sig'], color="blue", label="Elastic Region")
        plt.xlabel("ln(eps)")
        plt.ylabel("ln(sig)")
        plt.title("Elastic Part of Stress-Strain Curve")
        plt.legend()
        plt.grid(True)
        plt.show()

    elif operation.upper() == 'YOUNG_MODULUS':
        # Fit linear regression to ln(sig) ~ ln(eps)
        X = elastic_df['ln_eps'].values.reshape(-1, 1)
        y = elastic_df['ln_sig'].values.reshape(-1, 1)

        model = LinearRegression()
        model.fit(X, y)

        intercept = model.intercept_[0]  # log(E)
        return math.exp(intercept)       # Young's modulus in same units

    else:
        raise ValueError("Invalid operation. Use 'PLOT' or 'YOUNG_MODULUS'.")



def old_LN_S_E(df, operation):
    """
    This function analyzes the elastic part of a true stress-strain curve.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain 2 columns:
        - 'DL' : elongation (length change in mm)
        - 'F'  : force in Newtons

    operation : str
        - 'PLOT' : plots the elastic region of the true stress-strain curve
        - 'YOUNG_MODULUS' : calculates and returns Young's Modulus (E)

    Returns
    -------
    None if operation='PLOT'
    float if operation='YOUNG_MODULUS'
    """
    '''

    L0 = 40
    D0 = 9
    A0 = math.pi / 4 * (D0 ** 2)
    df['e'] = df['DL'] / L0
    df['S'] = df['F'] / A0
    df['eps'] = np.log(1 + df['e'])
    df['sig'] = df['S'] * (1 + df['e'])
    filtered_index = df[(df['eps'] >= 0.04) & (df['eps'] <= 0.08)].index
    "the elastic part of the curve is where the true strain(eps) is from 0.04 to 0.08"

    df['selected_ln_eps'] = np.nan

    df.loc[filtered_index, 'selected_ln_eps'] = df.loc[filtered_index, 'eps']
    df['selected_ln_eps'] = np.where(~df['selected_ln_eps'].isna(), np.log(df['selected_ln_eps']), df['selected_ln_eps'])
    df['selected_ln_sig'] = np.nan

    df.loc[filtered_index, 'selected_ln_sig'] = df.loc[filtered_index, 'sig']
    df['selected_ln_sig'] = np.where(~df['selected_ln_sig'].isna(), np.log(df['selected_ln_sig']), df['selected_ln_sig'])


    if operation == 'PLOT':
        plt.scatter(df['selected_ln_eps'].dropna(), df['selected_ln_sig'].dropna())
        plt.xlabel('ln_eps')
        plt.ylabel('ln_sig')
        plt.title('ln(sig) vs ln(eps) Plot')
        plt.grid(True)
        plt.show()
    elif operation == 'YOUNG_MODULUS':
        X = df['selected_ln_eps'].dropna().values.reshape(-1, 1)  # Independent variable
        y = df['selected_ln_sig'].dropna().values.reshape(-1, 1)  # Dependent variable
        model = LinearRegression()
        model.fit(X, y)
        intercept = model.intercept_[0]
        return math.exp(intercept)
        
    else: 
        print("Please enter proper operation")
        return
    '''
    print('Warning: This function is deprecated. Please use "LN_S_E" instead.')
    return None





def Oxygen_HeatCapacity_Analysis(df):
    """
    Calculate enthalpy and entropy of oxygen from heat capacity data
    and plot Cp, enthalpy, and entropy versus temperature.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing at least two columns:
        - 'T': Temperature values
        - 'Cp': Heat capacity at constant pressure

    Returns
    -------
    pandas.DataFrame
        Original DataFrame with added 'Enthalpy' and 'Entropy' columns.
    Also shows plots of:
        - Heat capacity vs temperature
        - Enthalpy and entropy vs temperature

    Example
    -------
    >>> df = pd.DataFrame({'T':[100,200,300],'Cp':[0.9,1.1,1.3]})
    >>> Oxygen_HeatCapacity_Analysis(df)
    """
    # Ensure temperature is sorted
    df = df.sort_values('T').reset_index(drop=True)

    # Calculate enthalpy as cumulative sum of Cp*dT (approximation)
    dT = df['T'].diff().fillna(0)  # Temperature differences
    df['Enthalpy'] = (df['Cp'] * dT).cumsum()

    # Calculate entropy as cumulative sum of Cp/T*dT (numerical approximation)
    df['Entropy'] = ((df['Cp'] / df['T']) * dT).cumsum()

    # Plotting
    fig, axs = plt.subplots(1, 2, figsize=(16, 6))
    
    # Heat capacity vs temperature
    axs[0].plot(df['T'], df['Cp'], color='blue')
    axs[0].set_xlabel('Temperature (T)')
    axs[0].set_ylabel('Heat Capacity (Cp)')
    axs[0].set_title('Heat Capacity vs Temperature')

    # Enthalpy and entropy vs temperature
    axs[1].plot(df['T'], df['Enthalpy'], label='Enthalpy', color='red')
    axs[1].plot(df['T'], df['Entropy'], label='Entropy', color='green')
    axs[1].set_xlabel('Temperature (T)')
    axs[1].set_ylabel('Enthalpy and Entropy')
    axs[1].set_title('Enthalpy and Entropy vs Temperature')
    axs[1].legend()

    plt.tight_layout()
    plt.show()

    return df







def Compression_TestAnalysis(df, operator, sample_name, density=0):
    """
    Analyze compression test data: plot stress-strain curve or calculate maximum strength.

    Parameters
    ----------
    df : pandas.DataFrame
        Compression test data containing at least two columns:
        - 'e': strain
        - 'S (Mpa)': stress in MPa
    operator : str
        Action to perform on data:
        - 'plot': plots stress-strain diagram
        - 'S_max': returns maximum stress
        - 'S_max/Density': returns specific maximum stress (requires density != 0)
    sample_name : str
        Name of the sample (used for plot label)
    density : float, optional
        Density of the sample (needed for 'S_max/Density'). Default is 0.

    Returns
    -------
    float or None
        Maximum stress if operator is 'S_max'.
        Specific maximum stress if operator is 'S_max/Density'.
        None if operator is 'plot'.

    Example
    -------
    >>> df = pd.DataFrame({'e':[0,0.01,0.02],'S (Mpa)':[10,20,15]})
    >>> Compression_TestAnalysis(df, 'S_max', 'Sample1')
    20
    >>> Compression_TestAnalysis(df, 'plot', 'Sample1')
    """
    # Check required columns exist
    if not {'e', 'S (Mpa)'}.issubset(df.columns):
        raise ValueError("DataFrame must contain 'e' and 'S (Mpa)' columns.")

    e = df['e']
    S = df['S (Mpa)']

    if operator == "S_max":
        return S.max()

    elif operator == "S_max/Density":
        if density == 0:
            raise ValueError("Density must be non-zero for 'S_max/Density' calculation.")
        return S.max() / density

    elif operator == "plot":
        font_label = {'family': 'Times New Roman', 'color': 'black', 'size': 15}
        font_legend = {'family': 'Times New Roman', 'size': 13}

        plt.plot(e, S, label=sample_name, linewidth=3)
        plt.xlabel("Strain (e)", fontdict=font_label, labelpad=5)
        plt.ylabel("Stress (S MPa)", fontdict=font_label, labelpad=5)

        plt.autoscale()
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
        plt.legend(frameon=False, prop=font_legend)

        plt.tick_params(axis='both', width=2)
        plt.tick_params(axis='both', which='minor', width=1)
        for spine in plt.gca().spines.values():
            spine.set_linewidth(2)
        plt.tick_params(axis='both', labelsize=11)

        plt.show()
        return None

    else:
        raise ValueError("Operator must be one of 'plot', 'S_max', 'S_max/Density'.")




def DMTA_TestAnalysis(df, operator, sample_name):
    """
    Analyze DMTA test data: find maxima or plot storage modulus, loss modulus, and tanδ.

    Parameters
    ----------
    df : pandas.DataFrame
        DMTA test data containing at least these columns:
        - 'Frequency (Hz)'
        - "E'-Storage Modulus (Mpa)"
        - Column 13 (loss modulus) or specify proper column
        - 'Tanδ'
    operator : str
        Action to perform on data:
        - 'storage_max': returns maximum storage modulus
        - 'loss_max': returns maximum loss modulus
        - 'tan_max': returns maximum Tanδ
        - 'plot_storage', 'plot_loss', 'plot_tan': plots corresponding data
    sample_name : str
        Name of the sample (used for plot label)

    Returns
    -------
    float or None
        Maximum value for storage, loss, or Tanδ if requested.
        None if plotting.

    Example
    -------
    >>> df = pd.DataFrame({
    ... 'Frequency (Hz)':[1,10,100],
    ... "E'-Storage Modulus (Mpa)":[100,150,200],
    ... df.columns[13]:[10,20,30],  # Loss modulus column
    ... 'Tanδ':[0.1,0.15,0.2]})
    >>> DMTA_TestAnalysis(df, 'storage_max', 'Sample1')
    200
    """
    # Check required columns exist
    required_cols = ["Frequency (Hz)", "E'-Storage Modulus (Mpa)", "Tanδ"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"DataFrame must contain '{col}' column.")
    # Assuming 14th column is loss modulus
    if df.shape[1] < 14:
        raise ValueError("DataFrame must have at least 14 columns to access Loss Modulus.")
    
    frequency = df["Frequency (Hz)"]
    storage_modulus = df["E'-Storage Modulus (Mpa)"]
    loss_modulus = df.iloc[:, 13].copy()
    tan_delta = df["Tanδ"]

    def plot_data(x, y, y_label):
        """Helper function to reduce repeated plotting code."""
        font_label = {'family': 'Times New Roman', 'color': 'black', 'size': 15}
        font_legend = {'family': 'Times New Roman', 'size': 13}

        plt.plot(x, y, label=sample_name, linewidth=3)
        plt.xlabel("Frequency (Hz)", fontdict=font_label, labelpad=5)
        plt.ylabel(y_label, fontdict=font_label, labelpad=5)

        plt.autoscale()
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
        plt.legend(frameon=False, prop=font_legend)

        plt.tick_params(axis='both', width=2)
        plt.tick_params(axis='both', which='minor', width=1)
        for spine in plt.gca().spines.values():
            spine.set_linewidth(2)
        plt.tick_params(axis='both', labelsize=11)
        plt.show()

    # Operator handling
    if operator == "storage_max":
        return storage_modulus.max()
    elif operator == "loss_max":
        return loss_modulus.max()
    elif operator == "tan_max":
        return tan_delta.max()
    elif operator == "plot_storage":
        plot_data(frequency, storage_modulus, "E'-Storage Modulus (Mpa)")
    elif operator == "plot_loss":
        plot_data(frequency, loss_modulus, "E''-Loss Modulus (Mpa)")
    elif operator == "plot_tan":
        plot_data(frequency, tan_delta, "Tanδ")
    else:
        raise ValueError("Operator must be one of 'storage_max','loss_max','tan_max','plot_storage','plot_loss','plot_tan'.")





def Find_MaxVerticalVelocity(df):
    """
    Find the maximum vertical flow velocity and its location from simulation data,
    and plot velocity versus position.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing at least two columns:
        - 'x(m)': position in meters
        - 'u(m/s)': vertical velocity in m/s

    Returns
    -------
    tuple
        (maximum velocity, location of maximum velocity)

    Example
    -------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'x(m)':[0,0.5,1.0],'u(m/s)':[0.1,0.3,0.2]})
    >>> Find_MaxVerticalVelocity(df)
    The maximum value of Flow Velocity for this problem is: 0.3
    Also this maximum value occurs in this location: 0.5
    (0.3, 0.5)
    """
    # Check required columns exist
    if not {'x(m)', 'u(m/s)'}.issubset(df.columns):
        raise ValueError("DataFrame must contain 'x(m)' and 'u(m/s)' columns.")
    
    x = np.array(df['x(m)'])
    u = np.array(df['u(m/s)'])

    # Find maximum velocity and its location
    u_max = np.max(u)
    index_max = np.argmax(u)
    loc_max = x[index_max]

    print(f'The maximum value of Flow Velocity for this problem is: {u_max}')
    print(f'Also this maximum value occurs in this location: {loc_max}')

    # Plot velocity vs position
    plt.scatter(x, u)
    plt.title('Flow Velocity', c='blue', family='Times New Roman', size=20, pad=20)
    plt.xlabel('x (m)', size=15, c='green')
    plt.ylabel('u velocity (m/s)', size=15, c='green')
    plt.show()

    return u_max, loc_max

        
def SolidificationStart(df, temp_sol):
    """
    Determine if solidification has started based on temperature profile,
    and plot temperature along the centerline.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing at least two columns:
        - 'x(m)': position in meters
        - 'T(K)': temperature in Kelvin
    temp_sol : float
        Solidus temperature of the material in Kelvin.

    Returns
    -------
    bool
        True if solidification has started (temperature <= solidus temperature),
        False otherwise.

    Example
    -------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'x(m)':[0,0.5,1],'T(K)':[1600,1550,1500]})
    >>> SolidificationStart(df, 1520)
    The solidification process has started.
    True
    """
    # Check required columns exist
    if not {'x(m)', 'T(K)'}.issubset(df.columns):
        raise ValueError("DataFrame must contain 'x(m)' and 'T(K)' columns.")

    x = np.array(df['x(m)'])
    y = np.array(df['T(K)'])

    # Plot temperature profile
    plt.scatter(x, y)
    plt.title('Temperature Profile', c='blue', family='Times New Roman', size=20, pad=20)
    plt.xlabel('x (m)', size=15, c='green')
    plt.ylabel('Temperature (K)', size=15, c='green')
    plt.show()

    temp_min = np.min(y)

    if temp_min > temp_sol:
        print('The solidification process has not started yet.')
        return False
    else:
        print('The solidification process has started.')
        return True




def Water_Hardness(df):
    """
    Evaluate water hardness based on metal content and pyrogenic compounds,
    filter out unsuitable water, calculate hardness (ppm), and plot results.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing at least the following columns:
        - 'name': sample name
        - 'Cu', 'Ni', 'Zn', 'pyro', 'Cya', 'Mg', 'Ca'

    Returns
    -------
    tuple
        - Filtered DataFrame with suitable water samples
        - List of DataFrames containing names of unavailable water samples
        - Displays a bar plot of water hardness (ppm) vs sample names

    Example
    -------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ... 'name':['W1','W2','W3'],
    ... 'Cu':[10,25,5],'Ni':[5,3,15],'Zn':[5,8,12],
    ... 'pyro':[50,120,90],'Cya':[1,3,0.5],'Mg':[10,15,5],'Ca':[20,25,15]})
    >>> Water_Hardness(df)
    """
    # Check required columns exist
    required_cols = ['name','Cu','Ni','Zn','pyro','Cya','Mg','Ca']
    if not set(required_cols).issubset(df.columns):
        raise ValueError(f"DataFrame must contain columns: {required_cols}")

    unavailable_water = []

    # Thresholds for filtering
    thresholds = {'Cu':20, 'Ni':10, 'Zn':10, 'pyro':100, 'Cya':2}

    filtered_df = df.copy()
    for col, thresh in thresholds.items():
        drop_index = filtered_df[filtered_df[col] > thresh].index
        for i in drop_index:
            unavailable_water.append(filtered_df.loc[i, ['name']])
        filtered_df = filtered_df.drop(drop_index)

    # Calculate hardness in ppm: Ca*2.5 + Mg*4.12
    ca_ppm = filtered_df['Ca'] * 2.5
    mg_ppm = filtered_df['Mg'] * 4.12
    hardness_ppm = ca_ppm + mg_ppm

    # Plot hardness
    plt.figure(figsize=(8,5))
    plt.bar(filtered_df['name'], hardness_ppm, color='skyblue')
    plt.xlabel('Sample Name')
    plt.ylabel('Hardness (ppm)')
    plt.title('Water Hardness')
    plt.show()

    return filtered_df, unavailable_water





def WearRate_Calculation(df, S, F, work='wear rate'):
    """
    Calculate wear rate of samples based on weight loss during a wear test.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing two columns:
        - 'weight before test': sample weight before the test
        - 'weight after test': sample weight after the test
    S : float
        Sliding distance during the test (in meters)
    F : float
        Normal force applied during the test (in Newtons)
    work : str, optional
        Type of calculation, default is 'wear rate'

    Returns
    -------
    float
        Wear rate (WR) in units of mass/(force*distance)

    Example
    -------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ... 'weight before test':[5.0,4.8,5.2],
    ... 'weight after test':[4.9,4.7,5.1]})
    >>> WearRate_Calculation(df, S=100, F=50)
    0.002
    """
    # Check required columns exist
    required_cols = ['weight before test', 'weight after test']
    if not set(required_cols).issubset(df.columns):
        raise ValueError(f"DataFrame must contain columns: {required_cols}")

    wb = np.array(df['weight before test'])
    wa = np.array(df['weight after test'])

    # Weight loss for each sample
    wl = wb - wa

    # Average weight loss
    m = wl.mean()

    if work.lower() == 'wear rate':
        WR = m / (S * F)  # Wear rate = average weight loss / (sliding distance * force)
        return WR
    else:
        raise ValueError("Unsupported work type. Only 'wear rate' is supported.")




def WearBar_Plot(df_list, S=300, F=5, work='bar'):
    """
    Calculate wear rate for multiple samples and plot as a bar chart.

    Parameters
    ----------
    df_list : list of pandas.DataFrame
        Each DataFrame must contain columns:
        - 'weight before test'
        - 'weight after test'
    S : float, optional
        Sliding distance in meters (default 300)
    F : float, optional
        Normal force in Newtons (default 5)
    work : str, optional
        Currently only 'bar' supported (default 'bar')

    Returns
    -------
    None
        Displays a bar plot of wear rates for the samples.

    Example
    -------
    >>> df1 = pd.DataFrame({'weight before test':[5.0],'weight after test':[4.9]})
    >>> df2 = pd.DataFrame({'weight before test':[4.8],'weight after test':[4.7]})
    >>> WearBar_Plot([df1, df2])
    """
    if work.lower() != 'bar':
        raise ValueError("Only 'bar' work type is supported.")

    if not isinstance(df_list, list) or len(df_list) == 0:
        raise ValueError("df_list must be a non-empty list of DataFrames.")

    wear_rates = []
    for i, df in enumerate(df_list):
        # Calculate wear rate for each sample using WearRate_Calculation()
        wr = WearRate_Calculation(df, S=S, F=F)
        wear_rates.append(wr)

    # Define sample labels (A, B, C, ...)
    sample_labels = [chr(65 + i) for i in range(len(wear_rates))]

    # Plot bar chart
    plt.figure(figsize=(8,5))
    plt.bar(sample_labels, wear_rates, color='green', width=0.5)
    plt.title('Wear Rate for Samples', color='green', size=14)
    plt.ylabel('Wear Rate (mg/N.m)', size=12, color='black')
    plt.xlabel('Sample', size=12, color='black')
    plt.show()









def PolarizationAnalysis(df, work):
    """
    Analyze polarization data: plot polarization curve or calculate corrosion potential.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing at least two columns:
        - 'Current density': current density in A/cm2
        - 'Potential': potential in V vs Ag/AgCl
    work : str
        Action to perform:
        - 'plot': plots the polarization curve (log(current) vs potential)
        - 'corrosion potential': returns the potential corresponding to the minimum current density

    Returns
    -------
    float or None
        Corrosion potential in volts if work='corrosion potential', None if plotting.

    Example
    -------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'Current density':[1e-6,1e-5,1e-4],'Potential':[0.1,0.2,0.3]})
    >>> PolarizationAnalysis(df, 'plot')  # Displays the plot
    >>> PolarizationAnalysis(df, 'corrosion potential')
    0.1
    """
    # Check required columns exist
    required_cols = ['Current density', 'Potential']
    if not set(required_cols).issubset(df.columns):
        raise ValueError(f"DataFrame must contain columns: {required_cols}")

    current_density = np.array(df['Current density'])
    potential = np.array(df['Potential'])

    # Take natural logarithm of current density
    log_cd = np.log(current_density)

    if work.lower() == 'plot':
        font_title = {'family': 'serif', 'color': 'b', 'size': 14}
        font_label = {'family': 'serif', 'color': 'k', 'size': 12}

        plt.plot(log_cd, potential, c='r', linewidth=3)
        plt.title('Polarization Curve', fontdict=font_title, pad=15)
        plt.xlabel('Log i (A/cm²)', fontdict=font_label)
        plt.ylabel('E (V vs Ag/AgCl)', fontdict=font_label)
        plt.show()
        return None

    elif work.lower() == 'corrosion potential':
        # Corrosion potential corresponds to minimum current density
        min_index = np.argmin(current_density)
        corrosion_potential = potential[min_index]
        return corrosion_potential

    else:
        raise ValueError("Unsupported work type. Use 'plot' or 'corrosion potential'.")







def StatisticalAnalysis(df, operation):
    """
    Perform statistical analysis or plots on a DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame with numeric features.
    operation : str
        Operation to perform:
        - 'statistics' : prints min, max, median, quantiles, IQR, and z-score for each numeric feature
        - 'histogram'  : plots histograms for numeric features
        - 'correlation': plots correlation heatmap
        - 'pairplot'   : plots pairplot with regression lines

    Returns
    -------
    None
        Prints statistics or displays plots.

    Example
    -------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'A':[1,2,3,4],'B':[4,3,2,1]})
    >>> StatisticalAnalysis(df, 'statistics')
    >>> StatisticalAnalysis(df, 'histogram')
    >>> StatisticalAnalysis(df, 'correlation')
    >>> StatisticalAnalysis(df, 'pairplot')
    """
    if df.empty:
        raise ValueError("Input DataFrame is empty.")

    font_title = {'family': 'serif', 'color': 'black', 'size': 15}
    font_x = {'family': 'serif', 'color': 'black', 'size': 10}
    font_y = {'family': 'serif', 'color': 'black', 'size': 10}

    numeric_cols = df.select_dtypes(include=[np.number]).columns

    if operation.lower() == 'statistics':
        for c in numeric_cols:
            Q1 = df[c].quantile(0.25)
            Q3 = df[c].quantile(0.75)
            IQR = Q3 - Q1
            min_value = df[c].min()
            max_value = df[c].max()
            median_value = df[c].median()
            z_score = np.abs(stats.zscore(df[c]))
            print(f'Feature: {c}\n Min: {min_value}, Max: {max_value}, Median: {median_value}, '
                  f'Q1: {Q1}, Q3: {Q3}, IQR: {IQR}, Z-score: {z_score}\n')

    elif operation.lower() == 'histogram':
        for c in numeric_cols:
            plt.hist(df[c], label=c, color='green', edgecolor='black', linewidth=0.5)
            plt.legend()
            plt.title('Distribution', fontdict=font_title)
            plt.xlabel('Bins', fontdict=font_x)
            plt.ylabel('Frequency', fontdict=font_y)
            plt.show()

    elif operation.lower() == 'correlation':
        plt.figure(figsize=(18, 12), dpi=200)
        sns.heatmap(df.corr(), xticklabels=numeric_cols, yticklabels=numeric_cols, center=0,
                    annot=True, cmap='coolwarm')
        plt.title('Correlation', fontdict=font_title)
        plt.show()

    elif operation.lower() == 'pairplot':
        sns.pairplot(df[numeric_cols], markers='*', kind='reg')
        plt.suptitle('Relation Between Columns', fontsize=16)
        plt.show()

    else:
        raise ValueError("Unsupported operation. Choose from 'statistics', 'histogram', 'correlation', or 'pairplot'.")

         


 
def ParticleSizeAnalysis(df, operation):
    """
    Analyze particle size distribution: calculate average size or plot size distribution.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing at least two columns:
        - 'size': particle sizes (nm)
        - 'distribution': intensity (%) corresponding to each size
    operation : str
        Action to perform:
        - 'calculate': calculate and return the average particle size
        - 'plot'     : plot the particle size distribution curve

    Returns
    -------
    float or None
        Average particle size if operation='calculate', None if plotting.

    Example
    -------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'size':[10,20,30],'distribution':[30,50,20]})
    >>> ParticleSizeAnalysis(df, 'calculate')
    20
    >>> ParticleSizeAnalysis(df, 'plot')  # Displays the plot
    """
    # Check required columns exist
    required_cols = ['size', 'distribution']
    if not set(required_cols).issubset(df.columns):
        raise ValueError(f"DataFrame must contain columns: {required_cols}")

    if operation.lower() == 'calculate':
        average_size = statistics.mean(df['size'])
        print('Average particle size is:', average_size)
        return average_size

    elif operation.lower() == 'plot':
        x = df['size']
        y = df['distribution']
        font_x = {'family':'Times New Roman', 'color':'k', 'size':15}
        font_y = {'family':'Times New Roman', 'color':'k', 'size':15}

        plt.plot(x, y, marker='o', color='blue')
        plt.xlabel('Size (nm)', fontdict=font_x)
        plt.ylabel('Intensity (%)', fontdict=font_y)
        plt.title('Particle Size Distribution', fontsize=16)
        plt.grid(True)
        plt.show()

    else:
        raise ValueError("Unsupported operation. Choose 'calculate' or 'plot'.")



    





def LoadPositionAnalysis(df, operation, area, length):
    """
    Analyze Load-Position data: generate curves, calculate stress-strain, normalized stress-strain, or energy absorption density.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing two columns:
        - 'Load (kN)': load values
        - 'Position (mm)': position values
    operation : str
        Operation to perform:
        - 'LPC' or 'Load-Position Curve'
        - 'SSCal' or 'Stress-Strain Calculation'
        - 'SSC' or 'Stress-Strain Curve'
        - 'NSSCal' or 'Normal Stress-Strain Calculation'
        - 'NSSC' or 'Normal Stress-Strain Curve'
        - 'EADCal' or 'EAD Calculation'
    area : float
        Cross-sectional area (mm²) for stress calculation
    length : float
        Gauge length (mm) for strain calculation

    Returns
    -------
    np.ndarray or float or None
        Depends on operation:
        - Stress-Strain arrays for 'SSCal' and 'NSSCal'
        - Energy absorption density for 'EADCal'
        - None for plotting operations

    Example
    -------
    >>> df = pd.DataFrame({'Load (kN)':[1,2,3],'Position (mm)':[0,1,2]})
    >>> LoadPositionAnalysis(df, 'LPC', 100, 50)  # Plot load-position curve
    >>> LoadPositionAnalysis(df, 'SSCal', 100, 50)  # Returns Stress-Strain array
    >>> LoadPositionAnalysis(df, 'EADCal', 100, 50)  # Returns energy absorption density
    """
    required_cols = ['Load (kN)', 'Position (mm)']
    if not set(required_cols).issubset(df.columns):
        raise ValueError(f"DataFrame must contain columns: {required_cols}")

    Load = np.array(df['Load (kN)'])
    Position = np.array(df['Position (mm)'])

    # Compute Stress and Strain for operations that need them
    Strain = Position / length
    Stress = (Load * 1000) / area  # Convert kN to N

    title_font = {'family': 'Times New Roman', 'color': 'black', 'size': 14}
    label_font = {'family': 'Times New Roman', 'color': 'black', 'size': 12}

    op = operation.lower()
    
    if op in ['lpc', 'load-position curve']:
        plt.plot(Position, Load, c='teal', lw=1)
        plt.title('Load-Position', fontdict=title_font, loc='center', pad=10)
        plt.xlabel('Position (mm)', fontdict=label_font, labelpad=5)
        plt.ylabel('Load (kN)', fontdict=label_font, labelpad=5)
        plt.xlim(0, np.max(Position))
        plt.ylim(0, np.max(Load))
        plt.grid(linewidth=0.5, color='grey', alpha=0.4)
        plt.show()
        return None

    elif op in ['sscal', 'stress-strain calculation']:
        return np.column_stack((Strain, Stress))

    elif op in ['ssc', 'stress-strain curve']:
        plt.plot(Strain, Stress, c='teal', lw=1)
        plt.title('Stress-Strain', fontdict=title_font, loc='center', pad=10)
        plt.xlabel('Strain (-)', fontdict=label_font, labelpad=5)
        plt.ylabel('Stress (MPa)', fontdict=label_font, labelpad=5)
        plt.grid(linewidth=0.5, color='grey', alpha=0.4)
        plt.show()
        return None

    elif op in ['nsscal', 'normal stress-strain calculation']:
        N_Strain = Strain / np.max(Strain)
        N_Stress = Stress / np.max(Stress)
        return np.column_stack((N_Strain, N_Stress))

    elif op in ['nssc', 'normal stress-strain curve']:
        N_Strain = Strain / np.max(Strain)
        N_Stress = Stress / np.max(Stress)
        plt.plot(N_Strain, N_Stress, c='teal', lw=1)
        plt.title('Normalized Stress-Strain', fontdict=title_font, loc='center', pad=10)
        plt.xlabel('Normalized Strain (-)', fontdict=label_font, labelpad=5)



def SI_Calculation(df, P, PC, Density=1):
    """
    Calculate Separation Index (SI) and plot Flux, Rejection, and SI charts.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing columns: 'Mem Code', 'Flux', 'Rejection'.
    P : float
        Pressure (bar)
    PC : float
        Pollutant concentration in Feed (g/L)
    Density : float, optional
        Feed Density (g/cm3), default is 1.

    Returns
    -------
    SI : numpy.ndarray
        Array of Separation Index values for each membrane.

    Example
    -------
    >>> df = pd.DataFrame({
    ...     'Mem Code': ['M1', 'M2', 'M3'],
    ...     'Flux': [90, 150, 250],
    ...     'Rejection': [0.5, 0.7, 0.8]
    ... })
    >>> SI = SI_Calculation(df, P=5, PC=50)
    """
    J = df['Flux'].values
    R = df['Rejection'].values
    Mem_Code = df['Mem Code'].values

    # Separation Index calculation
    SI = (Density - (1 - R) * PC / 1000) * (J / (P * ((1 - R) ** 0.41)))

    font = {'family': 'serif', 'color': 'k', 'size': 20}

    def plot_bar(values, title, ylabel, thresholds):
        """
        Internal function to plot bar chart with hatch patterns based on thresholds
        """
        hatches = np.array([
            '.' if val < thresholds[0] else 'o' if val < thresholds[1] else 'O' 
            for val in values
        ])
        fig, ax = plt.subplots()
        bar_labels = ['.:Low', 'o:Medium', 'O:High']
        ax.bar(Mem_Code, values, color='w', edgecolor='c', hatch=hatches,
               linewidth=1, yerr=0.01 if ylabel=='Rejection' else 10,
               ecolor='c', width=0.85, label=bar_labels)
        plt.title(title, fontdict=font)
        plt.xlabel('Membrane code', fontdict=font)
        plt.ylabel(ylabel, fontdict=font)
        ax.legend(title=f'{ylabel} Range')
        plt.show()

    # Plot charts
    plot_bar(J, 'Flux Chart', 'Flux', [100, 200])
    plot_bar(R, 'Rejection Chart', 'Rejection', [0.6, 0.75])
    plot_bar(SI, 'SI Chart', 'SI', [250, 500])

    return SI



def SICalculation(f_loc,P,PC,Density=1):

    '''

    This function is used for Separation Index Calculation
    

    P : Pressure (bar)

    Density : Feed Density(g/cm3)

    PC :  Pollutant concentration in Feed (g/L)

    Returns Separation Index and Flux & Rejection & Rejection Charts

    '''

    '''

    Data=pd.read_excel(f_loc)

    Data.columns

    J=Data['Flux']

    R=Data['Rejection']



    SI=(Density-(1-R)*PC/1000)*(J/(P*((1-R)**0.41)))

    Mem_Code=np.array(Data['Mem Code'])

    Flux=np.array(Data['Flux'])

    Rejection=np.array(Data['Rejection'])



    font={'family':'serif','color':'k','size':'20'}

    

    c=np.array([])

    for i in range (0,len(Flux)):

        if Flux[i]<100:

            a=np.array(['.'])

            c=np.concatenate((c,a))

        elif Flux[i]<200:

            a=np.array(['o'])

            c=np.concatenate((c,a))

        else:

            a=np.array(['O'])

            c=np.concatenate((c,a))

    fig, ax = plt.subplots()

    

    bar_labels=['.:Low Flux','o:Medium Flux','O:High Flux']

    for i in range(0,len(Flux)-3):

        m=['_.']

        bar_labels=bar_labels+m

        

    ax.bar(Mem_Code,Flux,color='w',edgecolor='c',hatch=c,linewidth=1,yerr=10,ecolor='c',width=0.85,label=bar_labels)   

    plt.title('Flux Chart',fontdict=font)

    plt.xlabel('Membrane code',fontdict=font)

    plt.ylabel('Flux',fontdict=font)

    ax.legend(title='Flux Range')

    plt.show()

    

    d=np.array([])

    for i in range (0,len(Rejection)):

        if Rejection[i]<0.6:

            a=np.array(['.'])

            d=np.concatenate((d,a))

        elif Rejection[i]<0.75:

            a=np.array(['o'])

            d=np.concatenate((d,a))

        else:

            a=np.array(['O'])

            d=np.concatenate((d,a))

    fig, ax = plt.subplots()

    

    bar_labels=['.:Low Rejection','o:Medium Rejection','O:High Rejection']

    for i in range(0,len(Rejection)-3):

        m=['_.']

        bar_labels=bar_labels+m

        

    ax.bar(Mem_Code,Rejection,color='w',edgecolor='c',hatch=d,linewidth=1,yerr=0.01,ecolor='c',width=0.85,label=bar_labels)   

    plt.title('Rejection Chart',fontdict=font)

    plt.xlabel('Membrane code',fontdict=font)

    plt.ylabel('Rejection',fontdict=font)

    ax.legend(title='Rejection Range')

    plt.show()

    

    f=np.array([])

    for i in range (0,len(SI)):

        if SI[i]<250:

            a=np.array(['.'])

            f=np.concatenate((f,a))

        elif SI[i]<500:

            a=np.array(['o'])

            f=np.concatenate((f,a))

        else:

            a=np.array(['O'])

            f=np.concatenate((f,a))

    fig, ax = plt.subplots()

    

    bar_labels=['.:Low SI','o:Medium SI','O:High SI']

    for i in range(0,len(SI)-3):

        m=['_.']

        bar_labels=bar_labels+m

        

    ax.bar(Mem_Code,SI,color='w',edgecolor='c',hatch=f,linewidth=1,yerr=10,ecolor='c',width=0.85,label=bar_labels)   

    plt.title('SI Chart',fontdict=font)

    plt.xlabel('Membrane code',fontdict=font)

    plt.ylabel('SI',fontdict=font)

    ax.legend(title='SI Range')

    plt.show()    

    

    return SI
    '''
    print('Warning: This function is deprecated. Please use "SI_Calculation" instead.')
    return None








def Porosity(df, Density=1):
    """
    Calculate porosity of membranes and plot a Porosity Chart.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing columns: 'membrane', 'Ww', 'Wd', 'V'
        where Ww = weight of wet sample, Wd = weight of dry sample,
        V = sample volume.
    Density : float, optional
        Water density (g/cm3). Default is 1.

    Returns
    -------
    Porosity : numpy.ndarray
        Array of porosity values for each membrane.

    Example
    -------
    >>> df = pd.DataFrame({
    ...     'membrane': ['M1', 'M2', 'M3'],
    ...     'Ww': [2.5, 3.0, 2.8],
    ...     'Wd': [2.0, 2.4, 2.3],
    ...     'V': [1.0, 1.2, 1.1]
    ... })
    >>> porosity_values = Porosity(df)
    """
    Ww = df['Ww'].values
    Wd = df['Wd'].values
    V = df['V'].values
    membrane = df['membrane'].values

    # Porosity calculation
    Porosity = (Ww - Wd) / (Density * V)

    font = {'family': 'serif', 'color': 'k', 'size': 20}

    # Assign hatches based on porosity value
    hatches = np.array([
        '.' if val < 0.9 else 'o' if val < 1 else 'O'
        for val in Porosity
    ])

    fig, ax = plt.subplots()
    bar_labels = ['.:Low Porosity', 'o:Medium Porosity', 'O:High Porosity']

    ax.bar(membrane, Porosity, color='w', edgecolor='c', hatch=hatches,
           linewidth=1, yerr=0.05, ecolor='c', width=0.85, label=bar_labels)
    plt.title('Porosity Chart', fontdict=font)
    plt.xlabel('Membrane', fontdict=font)
    plt.ylabel('Porosity', fontdict=font)
    ax.legend(title='Porosity Range')
    plt.show()

    return Porosity




def Tortuosity(df, Density=1):
    """
    Calculate the pore tortuosity of membranes and plot a Tortuosity Chart.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing columns: 'membrane', 'Ww', 'Wd', 'V'
        where Ww = weight of wet sample, Wd = weight of dry sample, V = sample volume.
    Density : float, optional
        Water density (g/cm3). Default is 1.

    Returns
    -------
    Tortuosity : numpy.ndarray
        Array of tortuosity values for each membrane.

    Example
    -------
    >>> df = pd.DataFrame({
    ...     'membrane': ['M1', 'M2', 'M3'],
    ...     'Ww': [2.5, 3.0, 2.8],
    ...     'Wd': [2.0, 2.4, 2.3],
    ...     'V': [1.0, 1.2, 1.1]
    ... })
    >>> tort_values = Tortuosity(df)
    """
    # Calculate porosity first
    Ww = df['Ww'].values
    Wd = df['Wd'].values
    V = df['V'].values
    Porosity = (Ww - Wd) / (Density * V)

    # Tortuosity calculation
    Tortuosity = ((2 - Porosity)**2) / Porosity
    membrane = df['membrane'].values

    font = {'family': 'serif', 'color': 'k', 'size': 20}

    # Assign hatches based on tortuosity value
    hatches = np.array([
        '.' if val < 0.75 else 'o' if val < 1.25 else 'O'
        for val in Tortuosity
    ])

    fig, ax = plt.subplots()
    bar_labels = ['.:Low Tortuosity', 'o:Medium Tortuosity', 'O:High Tortuosity']

    ax.bar(membrane, Tortuosity, color='w', edgecolor='c', hatch=hatches,
           linewidth=1, yerr=0.05, ecolor='c', width=0.85, label=bar_labels)
    plt.title('Tortuosity Chart', fontdict=font)
    plt.xlabel('Membrane', fontdict=font)
    plt.ylabel('Tortuosity', fontdict=font)
    ax.legend(title='Tortuosity Range')
    plt.show()

    return Tortuosity








def Pore_Size(df, A, P, Vis=8.9e-4, Density=1):
    """
    Calculate the pore size of membranes and plot a Pore Size Chart.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing columns: 'membrane', 'Ww', 'Wd', 'V', 'q', 'l'
        Ww = weight of wet sample (g), Wd = weight of dry sample (g), V = sample volume (cm3)
        q = flow rate (m3/s), l = membrane thickness (m)
    A : float
        Effective surface area of the membrane (m2)
    P : float
        Operational pressure (Pa)
    Vis : float, optional
        Water viscosity (Pa.s). Default is 8.9e-4
    Density : float, optional
        Water density (g/cm3). Default is 1

    Returns
    -------
    Pore_Size : numpy.ndarray
        Array of pore size values in nm.

    Example
    -------
    >>> df = pd.DataFrame({
    ...     'membrane': ['M1', 'M2', 'M3'],
    ...     'Ww': [2.5, 3.0, 2.8],
    ...     'Wd': [2.0, 2.4, 2.3],
    ...     'V': [1.0, 1.2, 1.1],
    ...     'q': [1e-6, 1.2e-6, 0.9e-6],
    ...     'l': [1e-3, 1.1e-3, 0.9e-3]
    ... })
    >>> pore_sizes = Pore_Size(df, A=0.01, P=2e5)
    """

    # Calculate Porosity
    Ww = df['Ww'].values
    Wd = df['Wd'].values
    V = df['V'].values
    Porosity = (Ww - Wd) / (Density * V)

    q = df['q'].values
    l = df['l'].values
    Pore_Size = ((2.9 - 1.75*Porosity) * (8*Vis*q*l/1000) / (Porosity*A*P))**0.5 * 1e9
    membrane = df['membrane'].values

    font = {'family':'serif','color':'k','size':20}

    # Assign hatches based on pore size
    hatches = np.array(['.' if val < 4.5 else 'o' if val < 5.5 else 'O' for val in Pore_Size])

    fig, ax = plt.subplots()
    bar_labels = ['.:Low Pore Size', 'o:Medium Pore Size', 'O:High Pore Size']
    ax.bar(membrane, Pore_Size, color='w', edgecolor='c', hatch=hatches,
           linewidth=1, yerr=0.05, ecolor='c', width=0.85, label=bar_labels)

    plt.title('Pore Size Chart', fontdict=font)
    plt.xlabel('Membrane', fontdict=font)
    plt.ylabel('Pore Size (nm)', fontdict=font)
    ax.legend(title='Pore Size Range')
    plt.show()

    return Pore_Size




def Signal_To_Noise_Ratio(data, application):
    """
    Calculate and optionally plot signal, noise, or SNR from experimental data.

    Parameters
    ----------
    data : DataFrame
        Experimental data with columns:
            1- 'location': measurement locations
            2- 'Signal Strength': signal power in dBm
            3- 'Noise Power': noise power in dBm
    application : str
        One of the following:
            'plot signal' - plots the signal column
            'plot noise' - plots the noise column
            'plot snr'   - plots the signal-to-noise ratio

    Returns
    -------
    mx : float
        Maximum signal-to-noise ratio in dB
    """
    location = np.array(data['location'])
    signal = np.array(data['Signal Strength'])
    noise = np.array(data['Noise Power'])
    
    snr = signal - noise  # Compute signal-to-noise ratio

    # Plot based on application
    app = str(application).lower()
    if app == 'plot signal':
        plt.plot(location, signal, color='blue', marker='o')
        plt.title('Signal Power at Each Location')
        plt.xlabel('Location')
        plt.ylabel('Signal Power (dBm)')
        plt.grid(True)
        plt.show()
    elif app == 'plot noise':
        plt.plot(location, noise, color='red', marker='x')
        plt.title('Noise Power at Each Location')
        plt.xlabel('Location')
        plt.ylabel('Noise Power (dBm)')
        plt.grid(True)
        plt.show()
    elif app == 'plot snr':
        plt.plot(location, snr, color='green', marker='s')
        plt.title('Signal-to-Noise Ratio at Each Location')
        plt.xlabel('Location')
        plt.ylabel('SNR (dB)')
        plt.grid(True)
        plt.show()
    
    mx = snr.max()
    return mx




def Polarization_Control(data: pd.DataFrame, application: str):
    """
    Analyze polymerization process data and either visualize trends or return key values.

    Parameters
    ----------
    data : pd.DataFrame
        A DataFrame containing the following required columns:
        - 'time' (float or int): Time in seconds
        - 'temp' (float): Temperature in °C
        - 'pressure' (float): Pressure in Pa
        - 'percent' (float): Reaction progress percentage (0–100)

    application : str
        Selects the analysis/plotting mode. Options:
        - 'temp_time'       : Plot Temperature vs Time
        - 'pressure_time'   : Plot Pressure vs Time
        - 'percent_time'    : Plot Reaction Percent vs Time
        - '100% reaction'   : Return (temperature, pressure) when polymerization reaches 100%
        - 'Max_pressure'    : Return maximum process pressure
        - 'Max_temp'        : Return maximum process temperature

    Returns
    -------
    tuple | float | None
        - (temp, pressure) if application is '100% reaction'
        - max pressure (float) if application is 'Max_pressure'
        - max temperature (float) if application is 'Max_temp'
        - None if plotting is performed

    Examples
    --------
    >>> df = pd.DataFrame({
    ...     'time': [0, 10, 20, 30],
    ...     'temp': [25, 50, 75, 100],
    ...     'pressure': [1, 2, 3, 4],
    ...     'percent': [0, 30, 70, 100]
    ... })
    >>> Polarization_Control(df, 'temp_time')  # Plots Temperature vs Time
    >>> Polarization_Control(df, 'Max_temp')
    100
    >>> Polarization_Control(df, '100% reaction')
    (100, 4)
    """
    
    # Extract arrays
    time = np.array(data['time'])
    temp = np.array(data['temp'])
    pressure = np.array(data['pressure'])   # FIX: corrected spelling ("pessure" → "pressure")
    reaction_percent = np.array(data['percent'])
    
    if application == 'temp_time':
        plt.plot(time, temp, c='g', linewidth=1.5)
        plt.title('Temperature variation', fontdict={'family':'serif','color':'black','size':16})
        plt.xlabel('time(s)', fontdict={'family':'serif','color':'black','size':16})
        plt.ylabel('Temperature (°C)', fontdict={'family':'serif','color':'black','size':16})
        plt.show()

    elif application == 'pressure_time':
        plt.plot(time, pressure, c='r', linewidth=1.5)
        plt.title('Pressure variation', fontdict={'family':'serif','color':'black','size':16})
        plt.xlabel('time(s)', fontdict={'family':'serif','color':'black','size':16})
        plt.ylabel('Pressure (Pa)', fontdict={'family':'serif','color':'black','size':16})
        plt.show()

    elif application == 'percent_time':  # FIX: lowercase name for consistency
        plt.plot(time, reaction_percent, c='b', linewidth=1.5)
        plt.title('Reaction progress', fontdict={'family':'serif','color':'black','size':16})
        plt.xlabel('time(s)', fontdict={'family':'serif','color':'black','size':16})
        plt.ylabel('Reaction Percent (%)', fontdict={'family':'serif','color':'black','size':16})
        plt.show()



def Desulfurization_Rate(data, application):
    """
    Analyze desulfurization rate with and without ultrasonic assistance.

    Parameters
    ----------
    data : pandas.DataFrame
        A dataframe containing the following columns:
        - 'Time': Measurement times
        - 'Desulfurization_With_Ultrasonic': Removal efficiency with ultrasonic
        - 'Desulfurization_Without_Ultrasonic': Removal efficiency without ultrasonic

    application : str
        Choose one of the following options:
        - "plot": plots the desulfurization with and without ultrasonic
        - "Max_Removal_With_Ultrasonic": returns maximum removal efficiency with ultrasonic
        - "Max_Removal_Without_Ultrasonic": returns maximum removal efficiency without ultrasonic

    Returns
    -------
    float or None
        - Returns the maximum value (float) if application is "Max_Removal_With_Ultrasonic"
          or "Max_Removal_Without_Ultrasonic".
        - Returns None if application is "plot".

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     "Time": [0, 10, 20, 30],
    ...     "Desulfurization_With_Ultrasonic": [5, 20, 45, 60],
    ...     "Desulfurization_Without_Ultrasonic": [3, 15, 35, 50]
    ... })
    >>> Desulfurization_Rate(df, "Max_Removal_With_Ultrasonic")
    60
    >>> Desulfurization_Rate(df, "plot")
    # Displays plot
    """

    # Extract columns as numpy arrays
    x = np.array(data['Time'])
    y1 = np.array(data['Desulfurization_With_Ultrasonic'])
    y2 = np.array(data['Desulfurization_Without_Ultrasonic'])
    
    if application == 'plot':
        # Plot both cases
        plt.plot(x, y1, marker='*', mec='r', mfc='y', ms=10, ls='-', linewidth=2,
                 color='r', label='With Ultrasonic')
        plt.plot(x, y2, marker='o', mec='g', mfc='y', ms=10, ls='--', linewidth=2,
                 color='g', label='Without Ultrasonic')

        myfont = {'family': 'serif', 'color': 'red', 'size': 15}
        
        plt.xlabel('Time')
        plt.ylabel('Desulfurization (%)')
        plt.title('Sulfur Removal Plot', fontdict=myfont)
        plt.legend()
        plt.grid(True)
        plt.show()
        return None

    elif application == 'Max_Removal_With_Ultrasonic':
        return y1.max()

    elif application == 'Max_Removal_Without_Ultrasonic':
        return y2.max()

    else:
        raise ValueError("Invalid application. Choose from: 'plot', "
                         "'Max_Removal_With_Ultrasonic', "
                         "'Max_Removal_Without_Ultrasonic'")







def Imerssion_Test(data, application):
    """
    Analyze immersion test data for weight gain/loss over time.

    Parameters
    ----------
    data : pandas.DataFrame
        A DataFrame containing immersion test results with the following columns:
        - 'time'
        - 'Mg'
        - 'Mg_H'
        - 'Mg_Pl'
        - 'Mg_HPl'
    application : str
        The analysis to perform:
            - 'plot' : Plot the changes of weight (%) vs. time (days).
            - 'More_Bioactive' : Return the sample with the highest weight gain (more bioactive).
            - 'Less_Bioactive' : Return the sample with the greatest weight loss (less bioactive).

    Returns
    -------
    None or str
        - If application == 'plot', displays a plot and returns None.
        - If application == 'More_Bioactive', returns the name of the most bioactive sample.
        - If application == 'Less_Bioactive', returns the name of the least bioactive sample.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'time': [1, 2, 3],
    ...     'Mg': [0.1, 0.2, 0.3],
    ...     'Mg_H': [0.2, 0.3, 0.5],
    ...     'Mg_Pl': [0.15, 0.25, 0.35],
    ...     'Mg_HPl': [0.25, 0.35, 0.45]
    ... })
    >>> Imerssion_Test(df, 'More_Bioactive')
    'Mg_HPl'
    >>> Imerssion_Test(df, 'Less_Bioactive')
    'Mg'
    >>> Imerssion_Test(df, 'plot')  # plots the data
    """
    x = np.array(data['time'])
    y1 = np.array(data['Mg_HPl'])
    y2 = np.array(data['Mg_H'])
    y3 = np.array(data['Mg_Pl'])
    y4 = np.array(data['Mg'])

    if application.lower() == 'plot':
        # Plot all curves
        plt.plot(x, y1, marker='o', label='Mg_HPl')
        plt.plot(x, y2, marker='*', label='Mg_H')
        plt.plot(x, y3, marker='^', label='Mg_Pl')
        plt.plot(x, y4, marker='+', label='Mg')
        plt.title('The graph of changes in the weight of the samples in the SBF solution', c='r')
        plt.xlabel('Immersion Time (days)', c='g')
        plt.ylabel('Weight Gain (%)', c='g')
        plt.legend()
        plt.show()

    elif application.lower() == 'more_bioactive':
        # Find sample with highest maximum weight gain
        max_col = data[['Mg', 'Mg_H', 'Mg_Pl', 'Mg_HPl']].max().idxmax()
        return max_col

    elif application.lower() == 'less_bioactive':
        # Find sample with lowest minimum weight (largest weight loss)
        min_col = data[['Mg', 'Mg_H', 'Mg_Pl', 'Mg_HPl']].min().idxmin()
        return min_col






def Reaction_Conversion_Analysis(data: pd.DataFrame, app: str):
    """
    Analyze and visualize conversion data from a chemical reaction experiment.

    Parameters
    ----------
    data : pandas.DataFrame
        A DataFrame that must contain the following columns:
        - 'time' : time in seconds
        - 'temp' : temperature in Celsius
        - 'pressure' : pressure in bar
        - 'conv' : conversion percentage
        
    app : str
        Determines the action:
        - "PLOT_TEMP"        → plots Temperature vs. Time
        - "PLOT_PRESSURE"    → plots Pressure vs. Time
        - "PLOT_CONVERSION"  → plots Conversion vs. Time
        - "MAXIMUM_CONVERSION" → returns index and values at maximum conversion

    Returns
    -------
    result : int or None
        - If `app="MAXIMUM_CONVERSION"`, returns the index of maximum conversion.
        - Otherwise, returns None (just shows plots).

    Raises
    ------
    TypeError
        If `app` is not one of the accepted values.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     "time": [0, 1, 2, 3],
    ...     "temp": [300, 310, 315, 320],
    ...     "pressure": [1, 1.2, 1.3, 1.4],
    ...     "conv": [10, 20, 30, 50]
    ... })
    >>> Conversion_Analysis(df, "PLOT_TEMP")  # plots Temperature vs Time
    >>> Conversion_Analysis(df, "MAXIMUM_CONVERSION")
    maximum of temperature is  320
    maximum of conversion is  50
    The temperature in maximum conversion is  320 and the pressure is  1.4
    3
    """
    # Ensure required columns exist
    required_cols = {"time", "temp", "pressure", "conv"}
    if not required_cols.issubset(data.columns):
        raise ValueError(f"Data must contain columns: {required_cols}")

    time = np.array(data['time'])
    temp = np.array(data['temp'])
    pressure = np.array(data['pressure'])
    conv = np.array(data['conv'])

    if app.upper() == 'PLOT_TEMP':
        plt.plot(time, temp, color='black')
        plt.title('Temperature over Time')
        plt.xlabel('Time (s)')
        plt.ylabel('Temperature (°C)')
        plt.grid()
        plt.show()

    elif app.upper() == 'PLOT_PRESSURE':
        plt.plot(time, pressure, color='red')
        plt.title('Pressure over Time')
        plt.xlabel('Time (s)')
        plt.ylabel('Pressure (bar)')
        plt.grid()
        plt.show()

    elif app.upper() == 'PLOT_CONVERSION':
        plt.plot(time, conv, color='blue')
        plt.title('Conversion over Time')
        plt.xlabel('Time (s)')
        plt.ylabel('Conversion (%)')
        plt.grid()
        plt.show()

    elif app.upper() == 'MAXIMUM_CONVERSION':
        max_conv = conv.max()
        idx_max_conv = np.argmax(conv)
        print('maximum of temperature is', temp.max())
        print('maximum of conversion is', max_conv)
        print(
            'The temperature at maximum conversion is', temp[idx_max_conv],
            'and the pressure is', pressure[idx_max_conv]
        )
        return idx_max_conv

    else:
        raise TypeError('The data or application argument is not entered correctly.')





def Import_Data(File_Directory=None):
    
    
    if File_Directory is None:
        raise TypeError('Please enter the file directory to open!')
        
        
        
    try:
        data=pd.read_csv(File_Directory)
        
        
#for examle -->data=pd.read_csv('C:\\Users\\Parsis.Co\\Desktop\\CV.csv')
         
         
#1--->File not find  (by 'FileNotFoundError')       
    except FileNotFoundError:
        raise FileNotFoundError('Unfortunately, the desired file <',File_Directory,'>is not available ')
        
        
        
#2---> File is  empty (by 'pd.errors.EmptyDataError')
    except pd.errors.EmptyDataError:
         raise ValueError('The file: ',data, ' is empty. please a correct file.')
         
   
#3--->Format is wrong  (by 'pd.errors.ParserError')     
    except pd.errors.ParserError:
             raise ValueError('The file format is not valid, please import a <csv> format file')
             
             
 #4--->remove empty cells     
    if data.isnull().values.any():
        print('Empty cell founded and removed ')
        data.dropna(inplace=True)
        
        
 #5--->turn object to numeric form for both columns    
    for x in data['P']:
        if data['P'].dtype=='object':
            print('object element founded in potential column and converted to numeric form')
            data['P']=pd.to_numeric(data['P'])
            
    for y in data['C']:
        if data['C'].dtype=='object':
            print('object element founded in current density column and converted to numeric')
            data['P']=pd.to_numeric(data['P'])
            
            
 #6--->remove duplicated data         
    if data.duplicated().any():
        print('Duolicated elemets in rows founded and removed ')
        data=data.drop_duplicates()
        
        
    return data
        
            




def Fatigue_Test_Analysis(data, application):
    """
    Analyze fatigue test data and provide multiple metrics and plots.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain columns:
        - 'stress_amplitude' : Stress amplitude in MPa
        - 'number_of_cycles' : Number of cycles to failure (N)
    
    application : str
        Determines the operation:
        - "plot"                  : S-N plot (Stress vs. Number of Cycles)
        - "max stress amplitude"  : Maximum stress amplitude
        - "fatigue strength"      : Mean stress amplitude
        - "fatigue life"          : Mean number of cycles
        - "stress in one cycle"   : Basquin's equation for stress at N=1
        - "Sa"                    : Stress amplitude at N=1
        - "fatigue limit"         : Cycle where stress becomes constant
        - "std stress"            : Standard deviation of stress
        - "std cycles"            : Standard deviation of cycles

    Returns
    -------
    value : float or array-like or None
        Result depending on the chosen application.
    """

    stress_amplitude = np.array(data["stress_amplitude"])
    number_of_cycles = np.array(data["number_of_cycles"])

    if application.lower() == "plot":
        title_font = {"color": "black", 'family': 'Merriweather', 'size': 20}
        xy_label_font = {"color": "Magenta", 'family': 'Merriweather', 'size': 12}

        plt.plot(number_of_cycles, stress_amplitude, marker='o', c='c', mec='k', mfc='r', 
                 label='S-N Curve', linewidth=2)
        plt.title('S-N Curve (Fatigue Test)', fontdict=title_font, pad=13)
        plt.xscale('log')
        plt.xlabel('Number of Cycles to Failure (N)', fontdict=xy_label_font)
        plt.ylabel('Stress Amplitude (MPa)', fontdict=xy_label_font)
        plt.grid(True)
        plt.legend()
        plt.show()
        return None

    if application.lower() == 'max stress amplitude':
        return np.max(stress_amplitude)

    if application.lower() == 'fatigue strength':
        return np.mean(stress_amplitude)

    if application.lower() == 'fatigue life':
        return np.mean(number_of_cycles)

    if application.lower() == 'stress in one cycle':
        # Basquin's law: Sa = S_f * N^b, solve for stress at N=1
        log_N = np.log10(number_of_cycles)
        log_S = np.log10(stress_amplitude)
        b, log_Sf = np.polyfit(log_N, log_S, 1)  # Linear fit: log(S) = b*log(N) + log(Sf)
        stress_one_cycle = 10**(log_Sf)  # Stress at N=1
        return stress_one_cycle

    if application.lower() == 'sa':
        # Stress amplitude at N=1, if exact not present, take nearest
        idx = (np.abs(number_of_cycles - 1)).argmin()
        return stress_amplitude[idx]

    if application.lower() == 'fatigue limit':
        fatigue_limit = None
        for i in range(1, len(stress_amplitude)):
            if stress_amplitude[i] == stress_amplitude[i-1]:
                fatigue_limit = number_of_cycles[i]
                break
        return fatigue_limit

    if application.lower() == 'std stress':
        return np.std(stress_amplitude)

    if application.lower() == 'std cycles':
        return np.std(number_of_cycles)

    raise ValueError("Invalid application. Choose a valid option.")
    






