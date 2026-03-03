# -*- coding: utf-8 -*-
"""
Build script for Swath Accuracy Plotter.
Reads __version__ from swath_accuracy_plotter.py and runs PyInstaller with
exe name: Swath_Accuracy_Plotter_v<version>
"""
import re
import subprocess
import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
main_script = os.path.join(script_dir, "swath_accuracy_plotter.py")

with open(main_script, "r", encoding="utf-8") as f:
    for line in f:
        m = re.match(r'^__version__\s*=\s*["\']([^"\']+)["\']', line)
        if m:
            version = m.group(1)
            break
    else:
        sys.exit("Could not find __version__ in swath_accuracy_plotter.py")

exe_name = "Swath_Accuracy_Plotter_v" + version
icon_path = os.path.join(script_dir, "media", "mac.ico")

cmd = [
    sys.executable, "-m", "PyInstaller",
    "--onefile",
    "--noconsole",
    "--name=" + exe_name,
    "swath_accuracy_plotter.py",
    "--clean",
]
if os.path.isfile(icon_path):
    cmd.insert(-2, "--icon=" + icon_path)

print("Building exe:", exe_name + ".exe")
os.chdir(script_dir)
subprocess.run(cmd)
