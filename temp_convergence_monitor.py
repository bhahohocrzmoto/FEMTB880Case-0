# -*- coding: utf-8 -*-
# ============================================================
# temp_convergence_monitor.py
# ------------------------------------------------------------
# Logs per-iteration temperatures to CSV and (optionally) shows
# a live convergence plot in a WinForms window.
# ============================================================

import os
import sys
import csv

_USE_DOTNET_IO = False
File = None
Directory = None

try:
    import clr
    clr.AddReference("System")
    from System.IO import File as _File, Directory as _Directory
    File = _File
    Directory = _Directory
    _USE_DOTNET_IO = True
except Exception:
    _USE_DOTNET_IO = False


def _file_exists(path):
    try:
        if _USE_DOTNET_IO and File is not None:
            return File.Exists(path)
        return os.path.exists(path)
    except Exception:
        return False


def _dir_exists(path):
    try:
        if _USE_DOTNET_IO and Directory is not None:
            return Directory.Exists(path)
        return os.path.isdir(path)
    except Exception:
        return False


def _ensure_dir(path):
    if _dir_exists(path):
        return
    if _USE_DOTNET_IO and Directory is not None:
        Directory.CreateDirectory(path)
    else:
        os.makedirs(path)


def _open_csv_text(path, mode):
    # ANSYS Mechanical commonly hosts IronPython 2.7, whose built-in open()
    # does not accept the Python 3-only newline keyword argument. The csv
    # module therefore needs a runtime-specific open mode.
    if sys.version_info[0] >= 3:
        return open(path, mode, newline="")
    return open(path, mode + "b")


_CHARTING_AVAILABLE = False
Form = None
Application = None
DockStyle = None
Chart = None
ChartArea = None
Series = None
SeriesChartType = None
AxisType = None
AxisEnabled = None
ChartImageFormat = None

try:
    import clr
    clr.AddReference("System.Windows.Forms")
    clr.AddReference("System.Windows.Forms.DataVisualization")

    from System.Windows.Forms import Form as _Form, Application as _Application, DockStyle as _DockStyle
    from System.Windows.Forms.DataVisualization.Charting import (
        Chart as _Chart,
        ChartArea as _ChartArea,
        Series as _Series,
        SeriesChartType as _SeriesChartType,
        AxisType as _AxisType,
        AxisEnabled as _AxisEnabled,
        ChartImageFormat as _ChartImageFormat,
    )

    Form = _Form
    Application = _Application
    DockStyle = _DockStyle
    Chart = _Chart
    ChartArea = _ChartArea
    Series = _Series
    SeriesChartType = _SeriesChartType
    AxisType = _AxisType
    AxisEnabled = _AxisEnabled
    ChartImageFormat = _ChartImageFormat
    _CHARTING_AVAILABLE = True
except Exception:
    _CHARTING_AVAILABLE = False


class TempConvergenceMonitor(object):
    def __init__(self, script_dir, cable_ids,
                 csv_name="temp_history.csv",
                 clear_csv=True,
                 live_plot=True,
                 window_title="Core temperature convergence"):
        self.script_dir = script_dir
        self.cable_ids = list(cable_ids)
        self.csv_path = os.path.join(script_dir, csv_name)
        self.live_plot_requested = bool(live_plot)
        self.window_title = window_title

        self._header = ["iter", "stage", "max_dT_C"]
        for cid in self.cable_ids:
            self._header.append("Tcore_{0}_C".format(cid))

        _ensure_dir(self.script_dir)

        if clear_csv and _file_exists(self.csv_path):
            if _USE_DOTNET_IO and File is not None:
                File.Delete(self.csv_path)
            else:
                os.remove(self.csv_path)

        self._write_header_once()

        self._form = None
        self._chart = None
        self._series_by_cable = {}
        self._series_dT = None

        if self.live_plot_requested and _CHARTING_AVAILABLE:
            self._init_live_plot()

    def _write_header_once(self):
        if _file_exists(self.csv_path):
            return
        with _open_csv_text(self.csv_path, "w") as f:
            w = csv.writer(f)
            w.writerow(self._header)

    def log(self, it, stage, T_core_dict, max_dT_C=None):
        row = [int(it), str(stage), "" if (max_dT_C is None) else float(max_dT_C)]
        for cid in self.cable_ids:
            row.append(float(T_core_dict.get(cid, 0.0)))

        with _open_csv_text(self.csv_path, "a") as f:
            w = csv.writer(f)
            w.writerow(row)

        if self._chart is not None:
            self._update_plot(int(it), T_core_dict, max_dT_C)

    def _init_live_plot(self):
        try:
            self._form = Form()
            self._form.Text = self.window_title
            self._form.Width = 900
            self._form.Height = 550

            self._chart = Chart()
            self._chart.Dock = DockStyle.Fill

            area = ChartArea("Main")
            area.AxisX.Title = "Iteration"
            area.AxisY.Title = "Core temperature [C]"
            area.AxisY2.Title = "max_dT [C]"
            area.AxisY2.Enabled = getattr(AxisEnabled, "True")
            self._chart.ChartAreas.Add(area)

            for cid in self.cable_ids:
                s = Series("Tcore_{0}".format(cid))
                s.ChartType = SeriesChartType.Line
                s.ChartArea = "Main"
                s.YAxisType = AxisType.Primary
                self._chart.Series.Add(s)
                self._series_by_cable[cid] = s

            sdt = Series("max_dT")
            sdt.ChartType = SeriesChartType.Line
            sdt.ChartArea = "Main"
            sdt.YAxisType = AxisType.Secondary
            self._chart.Series.Add(sdt)
            self._series_dT = sdt

            self._form.Controls.Add(self._chart)
            self._form.Show()

            # If DoEvents is used, avoid interacting with Mechanical UI during solve.
            Application.DoEvents()
        except Exception as exc:
            print("Warning: live plot disabled: {0}".format(exc))
            self._form = None
            self._chart = None
            self._series_by_cable = {}
            self._series_dT = None

    def _update_plot(self, it, T_core_dict, max_dT_C):
        try:
            for cid, s in self._series_by_cable.items():
                if cid in T_core_dict:
                    s.Points.AddXY(it, float(T_core_dict[cid]))
            if (self._series_dT is not None) and (max_dT_C is not None):
                self._series_dT.Points.AddXY(it, float(max_dT_C))
            Application.DoEvents()
        except Exception as exc:
            print("Warning: live plot update failed: {0}".format(exc))

    def save_png(self, png_name="temp_history.png"):
        if self._chart is None:
            return None
        out_path = os.path.join(self.script_dir, png_name)
        try:
            self._chart.SaveImage(out_path, ChartImageFormat.Png)
            return out_path
        except Exception as exc:
            print("Warning: could not save convergence chart: {0}".format(exc))
            return None
