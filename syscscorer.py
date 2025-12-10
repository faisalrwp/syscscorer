#!/usr/bin/env python3
"""
CyReScoF GUI - compute compromise vector C
Requirements:
  pip install numpy ttkbootstrap
  pip install pandas   # optional but recommended

Save as cyrescof_gui_c.py and run: python cyrescof_gui_c.py
"""
import os
import json
import traceback
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# try pandas, else fallback to csv module
try:
    import pandas as pd
    PANDAS = True
except Exception:
    pd = None
    PANDAS = False

import numpy as np

# ttkbootstrap for theming
try:
    import ttkbootstrap as tb
    from ttkbootstrap.constants import *
except Exception as e:
    raise RuntimeError("ttkbootstrap is required. Install with: pip install ttkbootstrap") from e

SETTINGS_FILE = "settings.json"
DEFAULT_THEME = "cosmo"

# -------------------------
# Settings handling
# -------------------------
def load_settings(path=SETTINGS_FILE):
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as fh:
                s = json.load(fh)
                if "theme" not in s:
                    s["theme"] = DEFAULT_THEME
                return s
        except Exception:
            return {"theme": DEFAULT_THEME}
    else:
        s = {"theme": DEFAULT_THEME}
        try:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(s, fh, indent=2)
        except Exception:
            pass
        return s

def save_settings(settings, path=SETTINGS_FILE):
    try:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(settings, fh, indent=2)
    except Exception:
        pass

# -------------------------
# CSV readers (pandas or csv)
# -------------------------
def read_matrix_csv(path):
    """
    Read dependency matrix CSV.
    Accepts labeled form: first column row labels + header of column names,
    or numeric square with header names.
    Returns: A (numpy.ndarray), names (list of strings)
    """
    if PANDAS:
        df = pd.read_csv(path, header=0)
        # detect labeled form: header length = n+1
        if df.shape[0] == df.shape[1] - 1:
            # use first column as index
            first_col = df.columns[0]
            df_indexed = df.set_index(first_col)
            if df_indexed.shape[0] != df_indexed.shape[1]:
                raise ValueError("Matrix CSV labeled form is not square.")
            names = list(df_indexed.index.astype(str))
            A = df_indexed.values.astype(float)
            return A, names
        elif df.shape[0] == df.shape[1]:
            names = list(df.columns.astype(str))
            A = df.values.astype(float)
            return A, names
        else:
            # try plain numeric CSV (no header)
            df2 = pd.read_csv(path, header=None)
            if df2.shape[0] == df2.shape[1]:
                n = df2.shape[0]
                names = [f"C{i+1}" for i in range(n)]
                return df2.values.astype(float), names
            raise ValueError("Matrix CSV shape not square. Provide n x n matrix.")
    else:
        import csv
        with open(path, newline='', encoding='utf-8') as fh:
            reader = csv.reader(fh)
            rows = list(reader)
        if not rows:
            raise ValueError("Empty matrix CSV.")
        header = rows[0]
        # labeled form?
        if len(header) >= 2 and len(rows) == len(header):
            names = []
            mat = []
            for r in rows[1:]:
                names.append(r[0])
                mat.append([float(x) for x in r[1:]])
            A = np.array(mat, dtype=float)
            if A.shape[0] != A.shape[1]:
                raise ValueError("Labeled matrix CSV not square.")
            return A, names
        else:
            # numeric with header columns
            try:
                names = header
                mat = []
                for r in rows[1:]:
                    mat.append([float(x) for x in r])
                A = np.array(mat, dtype=float)
                if A.shape[0] != A.shape[1] or A.shape[1] != len(names):
                    raise ValueError("Matrix CSV numeric form malformed.")
                return A, names
            except Exception as e:
                raise ValueError("Cannot parse matrix CSV: " + str(e)) from e

def read_scores_csv(path):
    """
    Read scores CSV. Detects component and score columns heuristically.
    Returns Svec (numpy array), names (list)
    """
    if PANDAS:
        df = pd.read_csv(path, header=0)
        cols_lower = [c.lower() for c in df.columns]
        comp_col = None
        score_col = None
        for orig, low in zip(df.columns, cols_lower):
            if low in ("component","name","system","cs","component_name"):
                comp_col = orig
            if low in ("s_s","s","score","resilience","resilience_score"):
                score_col = orig
        if comp_col is None:
            comp_col = df.columns[0]
        if score_col is None:
            if df.shape[1] >= 2:
                score_col = df.columns[1]
            else:
                raise ValueError("Scores CSV missing a score column.")
        names = [str(x) for x in df[comp_col].astype(str).tolist()]
        Svec = df[score_col].astype(float).values
        return Svec.astype(float), names
    else:
        import csv
        with open(path, newline='', encoding='utf-8') as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)
            if not rows:
                raise ValueError("Empty scores CSV.")
            header = reader.fieldnames
            comp_col = None
            score_col = None
            for h in header:
                if h.lower() in ("component","name","system"):
                    comp_col = h
                if h.lower() in ("s_s","s","score","resilience"):
                    score_col = h
            if comp_col is None:
                comp_col = header[0]
            if score_col is None:
                score_col = header[1] if len(header) > 1 else header[0]
            names = [r[comp_col] for r in rows]
            Svec = np.array([float(r[score_col]) for r in rows], dtype=float)
            return Svec, names

# -------------------------
# Compute C
# -------------------------
def compute_compromise_vector(A, Svec, beta, regularize_eps=1e-9):
    """
    Solve (I - beta*A) C = (1 - Svec)
    Returns Cvec (numpy array) and diagnostics (dict)
    """
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        raise ValueError("A must be square.")
    if len(Svec) != n:
        raise ValueError("Length of S vector must match matrix dimension.")
    I = np.eye(n)
    M = I - beta * A
    b = 1.0 - Svec

    diag = {}
    # spectral radius
    try:
        eigs = np.linalg.eigvals(beta * A)
        diag["spectral_radius_betaA"] = float(np.max(np.abs(eigs)))
    except Exception:
        diag["spectral_radius_betaA"] = None
    # determinant and condition
    try:
        diag["det_M"] = float(np.round(np.linalg.det(M), 8))
    except Exception:
        diag["det_M"] = None
    try:
        diag["cond_M"] = float(np.linalg.cond(M))
    except Exception:
        diag["cond_M"] = None

    used_reg = False
    try:
        C = np.linalg.solve(M, b)
    except np.linalg.LinAlgError:
        # regularize
        M_reg = M + regularize_eps * I
        try:
            C = np.linalg.solve(M_reg, b)
            used_reg = True
        except np.linalg.LinAlgError:
            C, *_ = np.linalg.lstsq(M_reg, b, rcond=None)
            used_reg = True

    diag["used_regularization"] = used_reg
    return np.array(C), diag

# -------------------------
# GUI App
# -------------------------
class CyReScoFGui:
    def __init__(self, root):
        # load settings
        self.settings = load_settings()
        theme = self.settings.get("theme", DEFAULT_THEME)

        # create window with theme
        self.root = tb.Window(themename=theme)
        self.root.title("CyReScoF - Compromise Vector (C) Calculator")
        self.root.geometry("980x640")

        # style object (access available themes)
        self.style = self.root.style

        # state
        self.matrix_path = tk.StringVar()
        self.scores_path = tk.StringVar()
        self.beta_var = tk.DoubleVar(value=0.65)
        self.last_names = []
        self.last_S = None
        self.last_C = None
        self.last_diag = {}

        # Notebook
        nb = ttk.Notebook(self.root)
        nb.pack(fill="both", expand=True, padx=10, pady=8)

        # Compute tab
        tab_compute = ttk.Frame(nb)
        nb.add(tab_compute, text="Compute")

        # Settings tab
        tab_settings = ttk.Frame(nb)
        nb.add(tab_settings, text="Settings")

        # Sensitivity tab
        tab_sensitivity = ttk.Frame(nb)
        nb.add(tab_sensitivity, text="Sensitivity")
        self.build_sensitivity_tab(tab_sensitivity)

        # Sensitivity Graph tab
        tab_graph = ttk.Frame(nb)
        nb.add(tab_graph, text="Sensitivity Analysis Graph")
        self.build_sensitivity_graph_tab(tab_graph)


        # --- Compute tab layout ---
        top_frame = ttk.Frame(tab_compute, padding=8)
        top_frame.pack(fill="x")

        # Matrix input
        ttk.Label(top_frame, text="Dependency matrix (A) CSV:").grid(row=0, column=0, sticky="w")
        ttk.Entry(top_frame, textvariable=self.matrix_path, width=70).grid(row=0, column=1, padx=6)
        tb.Button(top_frame, text="Browse", bootstyle="info-outline", command=self.browse_matrix).grid(row=0, column=2, padx=4)

        # Scores input
        ttk.Label(top_frame, text="Component scores (S) CSV:").grid(row=1, column=0, sticky="w", pady=6)
        ttk.Entry(top_frame, textvariable=self.scores_path, width=70).grid(row=1, column=1, padx=6)
        tb.Button(top_frame, text="Browse", bootstyle="info-outline", command=self.browse_scores).grid(row=1, column=2, padx=4)

        # Beta and Compute
        ttk.Label(top_frame, text="Propagation factor β:").grid(row=2, column=0, sticky="w")
        ttk.Entry(top_frame, textvariable=self.beta_var, width=12).grid(row=2, column=1, sticky="w")
        tb.Button(top_frame, text="Compute C", bootstyle="success", command=self.on_compute).grid(row=2, column=2, padx=4, pady=6)

        tb.Button(top_frame, text="Save Results...", bootstyle="primary", command=self.save_results).grid(row=0, column=3, padx=6)
        tb.Button(top_frame, text="Create Templates", bootstyle="outline-secondary", command=self.create_templates).grid(row=1, column=3, padx=6)

        # Diagnostics area
        diag_frame = ttk.Labelframe(tab_compute, text="Diagnostics", padding=8)
        diag_frame.pack(fill="x", padx=8, pady=6)
        self.diag_text = tk.Text(diag_frame, height=5, wrap="word")
        self.diag_text.pack(fill="x")

        # Results table
        res_frame = ttk.Labelframe(tab_compute, text="Results (C)", padding=8)
        res_frame.pack(fill="both", expand=True, padx=8, pady=6)

        cols = ("component", "S_s", "C")
        self.tree = ttk.Treeview(res_frame, columns=cols, show="headings")
        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, anchor="center")
        vsb = ttk.Scrollbar(res_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        self.tree.pack(fill="both", expand=True)

        # status bar
        self.status_var = tk.StringVar(value="Ready")
        status = ttk.Label(self.root, textvariable=self.status_var, relief="sunken", anchor="w")
        status.pack(fill="x", side="bottom")

        # --- Settings tab layout ---
        settings_frame = ttk.Frame(tab_settings, padding=12)
        settings_frame.pack(fill="both", expand=True)

        ttk.Label(settings_frame, text="Select Theme:").grid(row=0, column=0, sticky="w")
        themes = list(self.style.theme_names())
        self.theme_var = tk.StringVar(value=theme)
        self.theme_cb = ttk.Combobox(settings_frame, values=themes, textvariable=self.theme_var, state="readonly", width=30)
        self.theme_cb.grid(row=0, column=1, padx=6, pady=6)
        tb.Button(settings_frame, text="Apply Theme", bootstyle="secondary", command=self.on_apply_theme).grid(row=0, column=2, padx=6)

        ttk.Label(settings_frame, text="Current theme:").grid(row=1, column=0, sticky="w", pady=6)
        self.current_theme_label = ttk.Label(settings_frame, text=theme)
        self.current_theme_label.grid(row=1, column=1, sticky="w", pady=6)

        ttk.Label(settings_frame, text="Theme preview controls:").grid(row=2, column=0, sticky="w", pady=8)
        # sample controls with different bootstyles to show theme colors
        preview_frame = ttk.Frame(settings_frame)
        preview_frame.grid(row=3, column=0, columnspan=3, sticky="w", pady=6)
        tb.Button(preview_frame, text="Primary", bootstyle="primary").grid(row=0, column=0, padx=6)
        tb.Button(preview_frame, text="Info", bootstyle="info").grid(row=0, column=1, padx=6)
        tb.Button(preview_frame, text="Success", bootstyle="success").grid(row=0, column=2, padx=6)
        tb.Button(preview_frame, text="Warning", bootstyle="warning").grid(row=0, column=3, padx=6)
        tb.Button(preview_frame, text="Danger", bootstyle="danger").grid(row=0, column=4, padx=6)

        ttk.Label(settings_frame, text="(Theme selection is saved to settings.json)").grid(row=4, column=0, columnspan=3, sticky="w", pady=12)

            # ---------------------------------------------------------
    # Save Graph as Image
    # ---------------------------------------------------------
    def save_graph_image(self):
        if self.sensitivity_results is None:
            messagebox.showwarning("No Data", "Run sensitivity analysis and plot a graph first.")
            return

        p = filedialog.asksaveasfilename(
            title="Save Graph Image",
            defaultextension=".png",
            filetypes=[
                ("PNG Image", "*.png"),
                ("JPEG Image", "*.jpg"),
                ("PDF Document", "*.pdf"),
                ("SVG Vector Image", "*.svg")
            ]
        )

        if not p:
            return

        try:
            self.fig.savefig(p, dpi=300, bbox_inches="tight")
            messagebox.showinfo("Saved", f"Graph saved successfully:\n{p}")
        except Exception as e:
            messagebox.showerror("Save Error", f"Failed to save image:\n{e}")


    # ---------------------------------------------------------
    # Sensitivity Analysis Tab
    # ---------------------------------------------------------
    def build_sensitivity_tab(self, tab):
        frame = ttk.Frame(tab, padding=10)
        frame.pack(fill="both", expand=True)

        # β-range inputs
        ttk.Label(frame, text="β minimum:").grid(row=0, column=0, sticky="w")
        self.beta_min_var = tk.DoubleVar(value=0.10)
        ttk.Entry(frame, textvariable=self.beta_min_var, width=10).grid(row=0, column=1, padx=6)

        ttk.Label(frame, text="β maximum:").grid(row=1, column=0, sticky="w")
        self.beta_max_var = tk.DoubleVar(value=0.95)
        ttk.Entry(frame, textvariable=self.beta_max_var, width=10).grid(row=1, column=1, padx=6)

        ttk.Label(frame, text="β step size:").grid(row=2, column=0, sticky="w")
        self.beta_step_var = tk.DoubleVar(value=0.05)
        ttk.Entry(frame, textvariable=self.beta_step_var, width=10).grid(row=2, column=1, padx=6)

        # Run button
        tb.Button(frame, text="Run Sensitivity Analysis",
                  bootstyle="success",
                  command=self.run_sensitivity_analysis).grid(row=3, column=0, columnspan=2, pady=10)

        # Export button
        tb.Button(frame, text="Export to CSV",
                  bootstyle="info",
                  command=self.export_sensitivity_csv).grid(row=3, column=2, padx=10)

        # Results table
        table_frame = ttk.Labelframe(frame, text="Sensitivity Results", padding=8)
        table_frame.grid(row=4, column=0, columnspan=4, pady=10, sticky="nsew")

        self.sens_tree = ttk.Treeview(table_frame, columns=["beta"], show="headings")
        self.sens_tree.heading("beta", text="β")
        self.sens_tree.column("beta", width=70)

        vsb = ttk.Scrollbar(table_frame, orient="vertical", command=self.sens_tree.yview)
        self.sens_tree.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        self.sens_tree.pack(fill="both", expand=True)

        frame.rowconfigure(4, weight=1)
        frame.columnconfigure(3, weight=1)

        # Storage for graph tab
        self.sensitivity_results = None


    # ---------------------------------------------------------
    # Execute Sensitivity Analysis
    # ---------------------------------------------------------
    def run_sensitivity_analysis(self):
        if self.last_S is None or self.last_names is None:
            messagebox.showwarning("Missing Data",
                                   "Load A & S and run Compute before running sensitivity analysis.")
            return

        # read beta inputs
        try:
            bmin = float(self.beta_min_var.get())
            bmax = float(self.beta_max_var.get())
            bstep = float(self.beta_step_var.get())
            if not (0 < bmin < 1 and 0 < bmax < 1 and 0 < bstep and bmin < bmax):
                raise ValueError()
        except Exception:
            messagebox.showerror("Input Error", "Invalid β-range or step size.")
            return

        # load A and S
        A, _ = read_matrix_csv(self.matrix_path.get())
        Svec = self.last_S

        # generate beta sequence
        betas = np.arange(bmin, bmax + bstep, bstep)

        # configure columns
        cols = ["beta"] + [f"C_{nm}" for nm in self.last_names]
        self.sens_tree.configure(columns=cols)
        for c in cols:
            self.sens_tree.heading(c, text=c)
            self.sens_tree.column(c, width=90, anchor="center")

        # clear table
        for item in self.sens_tree.get_children():
            self.sens_tree.delete(item)

        results = []

        # compute rows
        for beta in betas:
            try:
                Cvec, _ = compute_compromise_vector(A, Svec, beta)
                row = [beta] + list(Cvec)
                results.append(row)

                display_row = [f"{beta:.4f}"] + [f"{c:.6f}" for c in Cvec]
                self.sens_tree.insert("", "end", values=display_row)
            except Exception:
                row = [beta] + ["ERROR"] * len(self.last_names)
                results.append(row)
                self.sens_tree.insert("", "end",
                                      values=[f"{beta:.4f}"] + ["ERROR"] * len(self.last_names))

        # save for graph tab
        self.sensitivity_results = {
            "betas": betas,
            "C": np.array([r[1:] for r in results], dtype=object),
            "names": self.last_names
        }

        messagebox.showinfo("Done", "Sensitivity analysis completed.")


    # ---------------------------------------------------------
    # Export Sensitivity Results to CSV
    # ---------------------------------------------------------
    def export_sensitivity_csv(self):
        if self.sensitivity_results is None:
            messagebox.showwarning("No Data", "Run sensitivity analysis first.")
            return

        p = filedialog.asksaveasfilename(defaultextension=".csv",
                                         filetypes=[("CSV Files", "*.csv")])
        if not p:
            return

        import csv
        betas = self.sensitivity_results["betas"]
        C = self.sensitivity_results["C"]
        names = self.sensitivity_results["names"]

        with open(p, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["beta"] + names)
            for beta, Crow in zip(betas, C):
                w.writerow([beta] + list(Crow))

        messagebox.showinfo("Exported", f"Sensitivity results exported to:\n{p}")

    # ---------------------------------------------------------
    # Sensitivity Graph Tab
    # ---------------------------------------------------------
    def build_sensitivity_graph_tab(self, tab):
        frame = ttk.Frame(tab, padding=10)
        frame.pack(fill="both", expand=True)

        # Component selector
        ttk.Label(frame, text="Select Component(s):").grid(row=0, column=0, sticky="w")
        self.graph_listbox = tk.Listbox(frame, selectmode="extended", height=6)
        self.graph_listbox.grid(row=1, column=0, sticky="nsw", pady=5)

        tb.Button(frame, text="Plot",
                  bootstyle="primary",
                  command=self.plot_sensitivity_graph).grid(row=2, column=0, pady=10)

        tb.Button(frame, text="Save Graph as Image",
          bootstyle="info",
          command=self.save_graph_image).grid(row=3, column=0, pady=5)

        # Figure area
        graph_frame = ttk.Frame(frame)
        graph_frame.grid(row=0, column=1, rowspan=3, sticky="nsew", padx=10)

        self.fig = Figure(figsize=(6, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=graph_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        frame.columnconfigure(1, weight=1)
        frame.rowconfigure(1, weight=1)


    # ---------------------------------------------------------
    # Plot Sensitivity Graph
    # ---------------------------------------------------------
    def plot_sensitivity_graph(self):
        if self.sensitivity_results is None:
            messagebox.showwarning("No Data", "Run sensitivity analysis first.")
            return

        betas = self.sensitivity_results["betas"]
        C = self.sensitivity_results["C"]
        names = self.sensitivity_results["names"]

        # fill listbox if empty
        if self.graph_listbox.size() == 0:
            for nm in names:
                self.graph_listbox.insert("end", nm)

        selection = self.graph_listbox.curselection()
        if not selection:
            messagebox.showwarning("No Selection", "Select at least one component to plot.")
            return

        self.ax.clear()

        for idx in selection:
            comp_name = names[idx]
            y = C[:, idx].astype(float)
            self.ax.plot(betas, y, marker="o", label=comp_name)

        self.ax.set_title("Sensitivity Analysis: C vs β")
        self.ax.set_xlabel("β")
        self.ax.set_ylabel("C value")
        self.ax.grid(True)
        self.ax.legend()

        self.canvas.draw()


    # -------------------------
    # Callbacks
    # -------------------------
    def browse_matrix(self):
        p = filedialog.askopenfilename(title="Select dependency matrix CSV", filetypes=[("CSV files","*.csv"),("All files","*.*")])
        if p:
            self.matrix_path.set(p)

    def browse_scores(self):
        p = filedialog.askopenfilename(title="Select scores CSV", filetypes=[("CSV files","*.csv"),("All files","*.*")])
        if p:
            self.scores_path.set(p)

    def on_apply_theme(self):
        new_theme = self.theme_var.get()
        try:
            self.root.style.theme_use(new_theme)
            self.settings["theme"] = new_theme
            save_settings(self.settings)
            self.current_theme_label.config(text=new_theme)
            messagebox.showinfo("Theme Applied", f"Theme '{new_theme}' applied and saved.")
        except Exception as e:
            messagebox.showerror("Theme Error", f"Failed to apply theme '{new_theme}':\n{e}")

    def create_templates(self):
        # prompt folder
        folder = filedialog.askdirectory(title="Choose folder to save templates")
        if not folder:
            return
        names = ["CS1_IoMT","CS2_Edge","CS3_EHR","CS4_PACS","CS5_Ambulance","CS6_HCC","CS7_MDMS"]
        # matrix
        mat_path = os.path.join(folder, "dependency_matrix_template.csv")
        with open(mat_path, "w", newline='', encoding='utf-8') as fh:
            writer = csv.writer(fh)
            writer.writerow(["component"] + names)
            for n in names:
                writer.writerow([n] + ["0.00"] * len(names))
        # scores
        scores_path = os.path.join(folder, "component_scores_template.csv")
        with open(scores_path, "w", newline='', encoding='utf-8') as fh:
            writer = csv.writer(fh)
            writer.writerow(["component","S_s"])
            for n in names:
                writer.writerow([n, "0.85"])
        messagebox.showinfo("Templates Created", f"Templates saved:\n{mat_path}\n{scores_path}")

    def on_compute(self):
        mat_path = self.matrix_path.get().strip()
        scores_path = self.scores_path.get().strip()
        if not mat_path or not os.path.exists(mat_path):
            messagebox.showwarning("Missing Matrix", "Select a valid dependency matrix CSV first.")
            return
        if not scores_path or not os.path.exists(scores_path):
            messagebox.showwarning("Missing Scores", "Select a valid component scores CSV first.")
            return
        try:
            A, names_matrix = read_matrix_csv(mat_path)
        except Exception as e:
            messagebox.showerror("Matrix Error", f"Failed to read matrix:\n{e}")
            return
        try:
            Svec_raw, names_scores = read_scores_csv(scores_path)
        except Exception as e:
            messagebox.showerror("Scores Error", f"Failed to read scores:\n{e}")
            return

        # align by name if possible
        try:
            Svec, names_used = self.align_scores(names_matrix, names_scores, Svec_raw)
        except Exception as e:
            # Ask user if they want to try index alignment
            resp = messagebox.askyesno("Alignment error", f"Failed to align names: {e}\nTry aligning by index order instead?")
            if not resp:
                return
            if len(Svec_raw) != A.shape[0]:
                messagebox.showerror("Dimension mismatch", f"Scores length ({len(Svec_raw)}) != matrix size ({A.shape[0]}). Cannot compute.")
                return
            Svec = Svec_raw
            names_used = names_matrix

        # validate beta
        try:
            beta = float(self.beta_var.get())
            if not (0.0 < beta < 1.0):
                raise ValueError("beta must be in (0,1)")
        except Exception as e:
            messagebox.showerror("Beta error", f"Invalid beta: {e}")
            return

        self.status_var.set("Computing...")
        try:
            Cvec, diag = compute_compromise_vector(A, Svec, beta)
        except Exception as e:
            messagebox.showerror("Compute error", f"Computation failed:\n{e}\n{traceback.format_exc()}")
            self.status_var.set("Ready")
            return

        # populate results table
        for item in self.tree.get_children():
            self.tree.delete(item)
        for nm, s, c in zip(names_used, Svec, Cvec):
            self.tree.insert("", "end", values=(nm, f"{s:.6f}", f"{c:.6f}"))

        # diagnostics
        self.diag_text.delete("1.0", tk.END)
        lines = [
            f"spectral_radius_betaA: {diag.get('spectral_radius_betaA')}",
            f"det_M: {diag.get('det_M')}",
            f"cond_M: {diag.get('cond_M')}",
            f"used_regularization: {diag.get('used_regularization')}"
        ]
        self.diag_text.insert(tk.END, "\n".join(lines))

        self.last_names = names_used
        self.last_S = Svec
        self.last_C = Cvec
        self.last_diag = diag

        self.status_var.set("Computation complete")
        messagebox.showinfo("Done", "Computation finished. See results and diagnostics.")

    def align_scores(self, names_matrix, names_scores, Svec_raw):
        """
        Align scores to matrix names case-insensitively.
        Raises if some names missing.
        """
        map_scores = {n.lower(): i for i, n in enumerate(names_scores)}
        aligned = []
        for nm in names_matrix:
            key = nm.lower()
            if key not in map_scores:
                raise ValueError(f"Component '{nm}' not found in scores CSV.")
            idx = map_scores[key]
            aligned.append(Svec_raw[idx])
        return np.array(aligned, dtype=float), names_matrix

    def save_results(self):
        if self.last_C is None or not self.last_names:
            messagebox.showwarning("No results", "No computed results to save. Run Compute first.")
            return
        p = filedialog.asksaveasfilename(title="Save results CSV", defaultextension=".csv", filetypes=[("CSV files","*.csv")])
        if not p:
            return
        try:
            import csv as _csv
            with open(p, "w", newline='', encoding='utf-8') as fh:
                writer = _csv.writer(fh)
                writer.writerow(["component","S_s","C"])
                for nm, s, c in zip(self.last_names, self.last_S, self.last_C):
                    writer.writerow([nm, f"{float(s):.6f}", f"{float(c):.6f}"])
            messagebox.showinfo("Saved", f"Results saved to:\n{p}")
        except Exception as e:
            messagebox.showerror("Save error", f"Failed to save:\n{e}")

# -------------------------
# Run
# -------------------------
def main():
    settings = load_settings()
    theme = settings.get("theme", DEFAULT_THEME)
    # Start the app
    app = CyReScoFGui(None)
    # override root created inside class (tb.Window created there); start mainloop
    app.root.mainloop()

if __name__ == "__main__":
    main()
