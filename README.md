# Swath Accuracy Plotter

A PyQt6-based application for visualizing and analyzing multibeam echosounder swath accuracy data. Provides tools for plotting swath data, reference surfaces, uncertainty analysis, and exporting results.

## Description

The Swath Accuracy Plotter is a comprehensive toolkit for assessing multibeam echosounder performance. It enables users to:

- Load and analyze crossline data from multibeam systems
- Compare crossline soundings against reference surfaces
- Visualize accuracy metrics (bias, standard deviation) as a function of beam angle
- Apply various filters to reference surfaces and crossline data
- Export results as plots and GeoTIFF files
- Save and load analysis sessions

## Features

- **Interactive Plotting**: Visualization of swath accuracy data with customizable point sizes and opacities
- **Reference Surface Analysis**: Compare crossline data against reference surfaces with filtering options
- **Multiple Plot Types**: 
  - Accuracy plots (bias and standard deviation vs. beam angle)
  - Reference surface visualizations (depth, density, slope, uncertainty)
  - Crossline coverage plots
  - Tide plots
- **Advanced Filtering**: Filter data by depth, angle, backscatter, slope, density, and uncertainty
- **Session Management**: Save and load complete analysis sessions
- **Export Capabilities**: Export plots as PNG/JPG/TIF and reference surfaces as GeoTIFF
- **Unique Line Colors**: Color-code points by crossline file for easy identification
- **IHO Standards**: Visualize Special Order, Order 1a, 1b, 2, and 3 accuracy limits
- **Dark Mode GUI**: Fusion style with dark palette so the interface looks consistent on both light and dark system themes

## Requirements

- Python 3.7 or higher
- PyQt6
- matplotlib
- numpy
- scipy
- pyproj
- utm
- rasterio (for GeoTIFF export)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/seamapper/SwathAccuracy.git
cd SwathAccuracy
```

2. Install dependencies:
```bash
pip install PyQt6 matplotlib numpy scipy pyproj utm rasterio
```

## Usage

### Run from source

```bash
python swath_accuracy_plotter.py
```

(Use your project’s Python/venv if you use one, e.g. `path\to\.venv\Scripts\python.exe swath_accuracy_plotter.py`.)

### Build a standalone executable

The executable name is **Swath_Accuracy_Plotter_v** plus the version from the code (e.g. `Swath_Accuracy_Plotter_v2026.01.exe`).

- **Option 1 – Batch file (Windows)**  
  Double‑click **`build_exe.bat`** in the project root. It uses the Python in `..\.venv\Scripts\python.exe`; edit the bat file if your venv is elsewhere.

- **Option 2 – Python script**  
  From the project root:
  ```bash
  python build_exe.py
  ```
  (Use the same Python that has PyInstaller and the project dependencies installed.)

The built exe is written to the **`dist`** folder. The build uses PyInstaller with `--onefile`, `--noconsole`, and the app icon from `media/mac.ico` if present.

### Basic Workflow

1. **Load Reference Surface**: Use "Add Ref. Surface" to load a reference surface file (XYZ format in UTM projection)
2. **Load Crosslines**: Use "Add Crosslines" to load crossline data files (.all, .kmall, or ASCII.txt formats)
3. **Configure Filters**: Set up filters in the "Filter" tab to refine your analysis
4. **Calculate Accuracy**: Click "Calc Accuracy" to process the data
5. **Customize Plots**: Adjust plot settings in the "Plot" tab
6. **Export Results**: Use "Save All Plots" or "Export All to GeoTIFF" to save your results

### Supported File Formats

- **Reference Surfaces**: XYZ files (UTM projection, depth in meters, positive up)
- **Density Surfaces**: XYD files (sounding density per cell)
- **Crosslines**: 
  - Kongsberg .all files
  - Kongsberg .kmall files
  - Qimera ASCII.txt files
- **Tide Files**: .tid text files (space-delimited: YYYY/MM/DD hh:mm:ss amplitude)

## Project Structure

```
SwathAccuracy/
├── swath_accuracy_plotter.py  # Main application entry point
├── build_exe.py               # Build script: reads __version__, runs PyInstaller
├── build_exe.bat              # Windows batch launcher for build_exe.py
├── libs/                      # Core library modules
│   ├── swath_accuracy_lib.py  # Main plotting and analysis functions
│   ├── file_fun.py            # File I/O operations
│   ├── gui_widgets.py         # Custom GUI widgets
│   ├── parseEM.py             # EM file parsing
│   ├── readEM.py              # EM file reading
│   ├── kmall.py               # KMALL file support
│   └── swath_fun.py           # Swath utility functions
├── media/                     # Application media files
│   ├── CCOM_MAC.png           # Logo
│   └── mac.ico                # Application icon (for built exe)
├── LICENSE                     # BSD 3-Clause License
└── README.md                   # This file
```

## Version History

- **2026.01**: Dark mode GUI (Fusion style + dark palette); versioned exe build (build_exe.py / build_exe.bat)
- **2025.8**: Fix for save_all_plots import
- **2025.7**: Fix for save_all_plots function
- **2025.6**: Point size and opacity controls for accuracy and coverage plots, unique line colors feature
- **2025.5**: File management and export all to GeoTIFF functionality
- **2025.3**: GUI improvements and GeoTIFF export button
- **2025.2**: Data management improvements and GUI changes
- **2025.1**: Major rewrite with shaded relief, special order, and IHO order 1a-3 support

## License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details.

## Copyright

Copyright (c) 2025 Center for Coastal and Ocean Mapping / Joint Hydrographic Center, University of New Hampshire

## Authors

- kjerram
- pjohnson

## Acknowledgments

Developed at the Center for Coastal and Ocean Mapping / Joint Hydrographic Center, University of New Hampshire.

## Support

For issues, questions, or contributions, please use the GitHub Issues page.




