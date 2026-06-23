# -*- coding: utf-8 -*-
"""
Build script for Swath Accuracy Plotter.
Reads __version__ from swath_accuracy_plotter.py and runs PyInstaller with
exe name: Swath_Accuracy_Plotter_v<version>
"""
import ast
import subprocess
import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
main_script = os.path.join(script_dir, "swath_accuracy_plotter.py")

with open(main_script, "r", encoding="utf-8") as f:
    source = f.read()

version = None
tree = ast.parse(source, filename=main_script)
for node in tree.body:
    if isinstance(node, ast.Assign):
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id == "__version__":
                if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                    version = node.value.value
                elif isinstance(node.value, ast.Str):  # compatibility
                    version = node.value.s
                break
    if version is not None:
        break

if version is None:
    sys.exit("Could not find string __version__ assignment in swath_accuracy_plotter.py")

exe_name = "Swath_Accuracy_Plotter_v" + version
icon_path = os.path.join(script_dir, "media", "mac.ico")
preferred_python = os.path.join(
    os.path.expanduser("~"),
    "AppData", "Local", "miniforge3", "python.exe"
)
python_exe = preferred_python if os.path.isfile(preferred_python) else sys.executable

cmd = [
    python_exe, "-m", "PyInstaller",
    "--onefile",
    "--noconsole",
    "--name=" + exe_name,
    "swath_accuracy_plotter.py",
    "--clean",
]
if os.path.isfile(icon_path):
    cmd.insert(-2, "--icon=" + icon_path)

print("Using Python:", python_exe)
print("Detected version:", version)
print("Building exe:", exe_name + ".exe")
os.chdir(script_dir)
subprocess.run(cmd, check=True)
