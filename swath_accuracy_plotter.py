# -*- coding: utf-8 -*-
"""
===============================================================================
Multibeam Echosounder Assessment Toolkit: Swath Accuracy Plotter
===============================================================================

Description:
    A PyQt6-based application for visualizing and analyzing multibeam 
    echosounder swath accuracy data. Provides tools for plotting swath data,
    reference surfaces, uncertainty analysis, and exporting results.

Author(s):
    kjerram, pjohnson

Created:
    Thu Apr 11 14:45:21 2019

Version:
    2025.5

Dependencies:
    - PyQt6
    - matplotlib
    - numpy
    - scipy
    - pyproj
    - utm

License:
    BSD 3-Clause License
    
    Copyright (c) 2019, kjerram, pjohnson
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
    1. Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.
    
    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.
    
    3. Neither the name of the copyright holders nor the names of its
       contributors may be used to endorse or promote products derived from
       this software without specific prior written permission.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.

===============================================================================
"""
# __version__ = "20240217"  # testing Qimera ASCII import
# __version__ = "0.1.2"  # new version with position time series duplicate filtering per IB Nuyina EM712 example
# __version__ = "0.1.3"  # added EM2042 model and update pyproj 3.6.1
# __version__ = "2025.1"  # rewrite with Cursor AI, added shaded relief, special order, and order 1a-3
# __version__ = "2025.2"  # changes in data management, gui changes, and swath pkl file management
# __version__ = "2025.3"  # changes in gui, added file management, and added export all to geotiff button
# __version__ = "2025.5"  # changes in gui, added file management, and added export all to geotiff button
__version__ = "2025.6"  # changes in gui, added point size and opacity for accuracy and coverage plots


import sys
import os

# Add the script's directory to the Python path to ensure imports work
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

# Set matplotlib backend before any other matplotlib imports
import matplotlib
matplotlib.use('QtAgg')

import numpy as np

from PyQt6 import QtWidgets, QtGui
from PyQt6.QtGui import QDoubleValidator, QIntValidator
from PyQt6.QtCore import Qt, QSize

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.colors import LightSource
from libs.gui_widgets import *
from libs.swath_accuracy_lib import *


class MainWindow(QtWidgets.QMainWindow):

    media_path = os.path.join(os.path.dirname(__file__), "media")

    def __init__(self, parent=None):
        print("MainWindow __init__ called")
        super(MainWindow, self).__init__()

        print("Setting up main window...")
        # set up main window
        self.mainWidget = QtWidgets.QWidget(self)
        self.setCentralWidget(self.mainWidget)
        self.setMinimumWidth(1000)
        self.setMinimumHeight(900)
        self.setWindowTitle('Swath Accuracy Plotter AI v.%s' % __version__)
        self.setWindowIcon(QtGui.QIcon(os.path.join(self.media_path, "icon.png")))

        print("Setting up Windows taskbar icon...")
        if os.name == 'nt':  # necessary to explicitly set taskbar icon
            import ctypes
            current_app_id = 'MAC.AccuracyPlotter.' + __version__  # arbitrary string
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(current_app_id)

        print("Calling setup()...")
        setup(self)  # initialize variables and plotter params

        print("Setting up layouts...")
        # set up three layouts of main window
        self.set_left_layout()
        print("Left layout done")
        self.set_center_layout()
        print("Center layout done")
        self.set_right_layout()
        print("Right layout done")
        self.set_main_layout()
        print("Main layout done")
        
        print("Initializing axes...")
        init_all_axes(self)
        print("Axes initialized")
        
        # init_swath_ax(self)
        # init_surf_ax(self)
        # init_tide_ax(self)
        print("Updating buttons...")
        update_buttons(self)
        print("Buttons updated")
        
        # Initialize default filter values
        print("Initializing default filter values...")
        self.initialize_default_filters()
        print("Default filters initialized")
        
        # Initialize default plot settings
        print("Initializing default plot settings...")
        self.initialize_default_plot_settings()
        print("Default plot settings initialized")

        # set up button controls for specific action other than refresh_plot
        # self.add_file_btn.clicked.connect(lambda: self.add_files('Kongsberg(*.all *.kmall)'))
        # self.add_file_btn.clicked.connect(lambda: add_acc_files(self, 'Kongsberg (*.all *.kmall)'))
        self.add_file_btn.clicked.connect(lambda: add_acc_files(self, 'Kongsberg or Qimera ASCII (*.all *.kmall *ASCII.txt)'))

        self.get_indir_btn.clicked.connect(lambda: add_acc_files(self, ['.all', '.kmall'], input_dir='',
                                                                 include_subdir=self.include_subdir_chk.isChecked()))
        self.get_outdir_btn.clicked.connect(lambda: get_output_dir(self))
        self.add_ref_surf_btn.clicked.connect(lambda: add_ref_file(self, 'Reference surface XYZ (*.xyz)'))
        self.add_dens_surf_btn.clicked.connect(lambda: add_dens_file(self, 'Density surface XYD (*.xyd)'))
        self.export_geotiff_btn.clicked.connect(lambda: export_all_to_geotiff(self))
        # self.add_tide_btn.clicked.connect(lambda: add_tide_file(self, 'Tide file (*.tid)'))
        self.add_tide_btn.clicked.connect(lambda: add_tide_file(self, 'Tide file (*.tid *.txt)'))
        self.rmv_file_btn.clicked.connect(lambda: remove_acc_files(self))
        # self.rmv_file_btn.clicked.connect(lambda: remove_files(self))
        self.clr_file_btn.clicked.connect(lambda: clear_files(self))
        self.show_path_chk.stateChanged.connect(lambda: show_file_paths(self))
        self.calc_accuracy_btn.clicked.connect(lambda: calc_accuracy(self))
        self.save_all_plots_btn.clicked.connect(self.handle_save_all_plots)
        # self.ref_proj_cbox.activated.connect(lambda: parse_ref_depth(self))
        self.ref_proj_cbox.activated.connect(lambda: update_ref_utm_zone(self))
        self.slope_win_cbox.activated.connect(lambda: update_ref_slope(self))
        self.slope_win_cbox.activated.connect(lambda: calc_accuracy(self, recalc_dz_only=True))
        # self.tide_unit_cbox.activated.connect(lambda: parse_tide(self, force_refresh=True))
        self.tide_unit_cbox.activated.connect(lambda: process_tide(self, unit_set_by_user=True))
        self.depth_mode_cbox.activated.connect(lambda: calc_accuracy(self, recalc_dz_only=True))
        self.waterline_tb.returnPressed.connect(lambda: update_buttons(self, recalc_acc=True))
        
        # connect filter management buttons
        self.save_filters_btn.clicked.connect(lambda: save_current_filters(self))
        self.load_last_filters_btn.clicked.connect(lambda: load_last_filters(self))
        self.default_filters_btn.clicked.connect(lambda: load_default_filters(self))
        
        # connect plot settings management buttons
        self.save_plot_settings_btn.clicked.connect(lambda: save_current_plot_settings(self))
        self.load_last_plot_settings_btn.clicked.connect(lambda: load_last_plot_settings(self))
        self.default_plot_settings_btn.clicked.connect(lambda: load_default_plot_settings(self))
        

        


        # set up event actions that call refresh_plot with appropriate lists of which plots to refresh

        # testing maps with dicts of more options... seems to refresh only on last item
        # gb_map = {self.custom_info_gb: {'refresh': ['acc'], 'active_tab': 1},
        #           self.plot_lim_gb: {'refresh': ['acc'], 'active_tab': 1},
        #           self.angle_gb: {'refresh': ['acc'], 'active_tab': 1},
        #           self.depth_xline_gb: {'refresh': ['acc'], 'active_tab': 1},
        #           self.depth_ref_gb: {'refresh': ['ref'], 'active_tab': 2},
        #           self.bs_gb: {'refresh': ['acc'], 'active_tab': 1},
        #           self.slope_gb: {'refresh': ['ref'], 'active_tab': 2},
        #           self.density_gb: {'refresh': ['ref'], 'active_tab': 2}}
                  # self.depth_gb: {'refresh': ['ref'], 'active_tab': 2}}

        # gb_map = {0: {'widget': self.custom_info_gb, 'refresh': ['acc'], 'active_tab': 1},
        #           1: {'widget': self.plot_lim_gb, 'refresh': ['acc'], 'active_tab': 1},
        #           2: {'widget': self.angle_gb, 'refresh': ['acc'], 'active_tab': 1},
        #           3: {'widget': self.depth_xline_gb, 'refresh': ['acc'], 'active_tab': 1},
        #           4: {'widget': self.depth_ref_gb, 'refresh': ['ref'], 'active_tab': 2},
        #           5: {'widget': self.bs_gb, 'refresh': ['acc'], 'active_tab': 1},
        #           6: {'widget': self.slope_gb, 'refresh': ['ref'], 'active_tab': 2},
        #           7: {'widget': self.density_gb, 'refresh': ['ref'], 'active_tab': 2}}
        #           # self.depth_gb: {'refresh': ['ref'], 'active_tab': 2}}

        gb_acc_map = [self.custom_info_gb]
                      # self.plot_lim_gb,
                      # self.angle_xline_gb,
                      # self.depth_xline_gb,
                      # self.bs_xline_gb,
                      # self.dz_gb]

        gb_ref_map = [self.slope_gb,
                      self.density_gb,
                      self.depth_ref_gb,
                      self.uncertainty_gb]

        gb_all_map = [self.pt_count_gb,
                      self.bin_decimation_gb,
                      self.plot_lim_gb]  # refresh acc and ref plot if depth_gb is activated

        cbox_map = [self.model_cbox,
                    self.pt_size_cbox,
                    self.pt_size_cov_cbox,
                    self.ref_cbox,
                    self.depth_mode_cbox]  # only one reference in cbox so far

        chk_acc_map = [self.show_acc_proc_text_chk,
                       # self.grid_lines_toggle_chk,
                       self.IHO_lines_toggle_chk,
                       self.set_zero_mean_chk,
                       self.unique_line_pt_colors_chk]

        chk_ref_map = [self.update_ref_plots_chk,
                       self.show_xline_cov_chk,
                       self.show_ref_proc_text_chk]

        chk_all_map = [self.grid_lines_toggle_chk,
                       self.show_u_plot_chk,
                       self.show_model_chk,
                       self.show_ship_chk,
                       self.show_cruise_chk,
                       self.show_shaded_relief_chk,
                       self.show_special_order_chk,
                       self.show_order_1a_chk,
                       self.show_order_1b_chk,
                       self.show_order_2_chk,
                       self.show_order_3_chk]

        tb_all_map = [self.ship_tb,
                      self.cruise_tb,
                      self.max_count_tb,
                      self.dec_fac_tb,
                      self.max_points_per_bin_tb,
                      self.pt_alpha_acc_tb,
                      self.pt_alpha_cov_tb]
                      # self.waterline_tb]

        tb_acc_map = [self.max_beam_angle_tb,
                      self.angle_spacing_tb,
                      self.max_bias_tb,
                      self.max_std_tb]

        tb_ref_map = [self.min_depth_ref_tb,
                      self.max_depth_ref_tb,
                      self.max_slope_tb,
                      self.min_dens_tb,
                      self.max_u_tb,
                      self.axis_margin_tb]

        tb_recalc_bins_map = [self.min_depth_xline_tb,
                              self.max_depth_xline_tb,
                              self.min_angle_xline_tb,
                              self.max_angle_xline_tb,
                              self.max_bs_xline_tb,
                              self.min_bs_xline_tb,
                              self.max_dz_tb,
                              self.max_dz_wd_tb,
                              self.min_bin_count_tb,
                              self.mean_bias_angle_lim_tb]

        tb_recalc_dz_map = [self.min_depth_ref_tb,
                            self.max_depth_ref_tb,
                            self.max_slope_tb,
                            self.min_dens_tb,
                            self.max_u_tb]

        gb_recalc_bins_map = [self.depth_xline_gb,
                              self.angle_xline_gb,
                              self.bs_xline_gb,
                              self.dz_abs_gb,
                              self.dz_pct_gb,
                              self.depth_mode_gb,
                              self.bin_count_gb,
                              self.flatten_mean_gb]

        gb_recalc_dz_map = [self.depth_ref_gb,
                            self.slope_gb,
                            self.uncertainty_gb,
                            self.density_gb]

        # for i, gb in gb_map.items():
        #     print('running through gb_map with gb=', gb)
        #     print('settings[refresh]=', gb['refresh'], 'and active_tab = ', gb['active_tab'])
        # # groupboxes tend to not have objectnames, so use generic sender string
        #     gb['widget'].clicked.connect(lambda: refresh_plot(self,
        #                                             refresh_list=gb['refresh'],
        #                                             sender=str(gb['widget'].objectName()).split('.')[-1],
        #                                             set_active_tab=gb['active_tab']))

        for gb in gb_acc_map:
            # groupboxes tend to not have objectnames, so use generic sender string
            gb.clicked.connect(lambda:
                               refresh_plot(self, refresh_list=['acc'], sender='GROUPBOX_CHK', set_active_tab=0))

        for gb in gb_ref_map:
            # groupboxes tend to not have objectnames, so use generic sender string
            gb.clicked.connect(lambda:
                               refresh_plot(self, refresh_list=['ref'], sender='GROUPBOX_CHK', set_active_tab=1))

            gb.clicked.connect(lambda: calc_accuracy(self, recalc_dz_only=True))  # recalculate xline dz after ref filt
            # Add plot_ref_surf call for crossline coverage on reference surface tabs
            gb.clicked.connect(lambda: plot_ref_surf(self))

        for gb in gb_all_map:
            # groupboxes tend to not have objectnames, so use generic sender string
            gb.clicked.connect(lambda:
                               refresh_plot(self, refresh_list=['acc', 'ref'], sender='GROUPBOX_CHK', set_active_tab=0))
            # Add plot_ref_surf call for crossline coverage on reference surface tabs
            gb.clicked.connect(lambda: plot_ref_surf(self))

        for gb in gb_recalc_bins_map:
            gb.clicked.connect(lambda: calc_accuracy(self, recalc_bins_only=True))

        for cbox in cbox_map:
            # lambda needs _ for cbox
            cbox.activated.connect(lambda _, sender=cbox.objectName():
                                   refresh_plot(self, sender=sender))  # refresh_list not specified, will update all
            # Add plot_ref_surf call for crossline coverage on reference surface tabs
            cbox.activated.connect(lambda _, sender=cbox.objectName():
                                   plot_ref_surf(self))

        # for cbox in cbox_recalc_map:
        #     lambda needs _ for cbox
            # cbox.activated.connect(lambda _, sender=cbox.objectName():
            #                        refresh_plot(self, sender=sender))  # refresh_list not specified, will update all

        for chk in chk_acc_map:
            # lambda needs _ for chk
            if chk == self.set_zero_mean_chk:
                # set_zero_mean_chk affects data calculation, so trigger recalculation
                chk.stateChanged.connect(lambda _, sender=chk.objectName():
                                         calc_accuracy(self, recalc_bins_only=True))
            else:
                # Other accuracy checkboxes are just display options
                chk.stateChanged.connect(lambda _, sender=chk.objectName():
                                         refresh_plot(self, refresh_list=['acc'], sender=sender, set_active_tab=0))

        for chk in chk_ref_map:
            # lambda needs _ for chk
            chk.stateChanged.connect(lambda _, sender=chk.objectName():
                                     refresh_plot(self, refresh_list=['ref'], sender=sender, set_active_tab=1))
            # Add plot_ref_surf call for crossline coverage on reference surface tabs
            chk.stateChanged.connect(lambda _, sender=chk.objectName():
                                     plot_ref_surf(self))

        for chk in chk_all_map:
            # lambda needs _ for chk
            chk.stateChanged.connect(lambda _, sender=chk.objectName():
                                     refresh_plot(self, sender=sender))  # refresh_list not specified, will update all
            # Add plot_ref_surf call for crossline coverage on reference surface tabs
            chk.stateChanged.connect(lambda _, sender=chk.objectName():
                                     plot_ref_surf(self))
        

        for tb in tb_acc_map:
            # lambda seems to not need _ for tb
            tb.returnPressed.connect(lambda sender=tb.objectName():
                                     calc_accuracy(self, recalc_bins_only=True))

        for tb in tb_ref_map:
            # lambda seems to not need _ for tb
            tb.returnPressed.connect(lambda sender=tb.objectName():
                                     refresh_plot(self, refresh_list=['ref'], sender=sender, set_active_tab=1))
            # Add plot_ref_surf call for crossline coverage on reference surface tabs
            tb.returnPressed.connect(lambda sender=tb.objectName():
                                     plot_ref_surf(self))


        for tb in tb_all_map:
            # lambda seems to not need _ for tb
            tb.returnPressed.connect(lambda sender=tb.objectName():
                                     refresh_plot(self, sender=sender))
            # Add plot_ref_surf call for crossline coverage on reference surface tabs
            tb.returnPressed.connect(lambda sender=tb.objectName():
                                     plot_ref_surf(self))

        for tb in tb_recalc_bins_map:
            # lambda seems to not need _ for tb
            tb.returnPressed.connect(lambda sender=tb.objectName():
                                     calc_accuracy(self, recalc_bins_only=True))

        for tb in tb_recalc_dz_map:
            tb.returnPressed.connect(lambda sender=tb.objectName():
                                     calc_accuracy(self, recalc_dz_only=True))

    def set_left_layout(self):
        btnh = 20  # height of file control button
        btnw = 140  # width of file control button

        # add reference surface import options (processed elsewhere, XYZ in meters positive up, UTM 1-60N through 1-60S)
        self.add_ref_surf_btn = PushButton('Add Ref. Surface', btnw, btnh, 'add_ref_surf_btn',
                                           'Add a reference surface in UTM projection with northing, easting, and '
                                           'depth in meters, with depth positive up / negative down.')
        self.add_dens_surf_btn = PushButton('Add Dens. Surface', btnw, btnh, 'add_ref_surf_btn',
                                           'Add a sounding density surface corresponding to the reference surface, if '
                                           'needed (e.g., Qimera .xyz exports do not include sounding density. '
                                           'A density layer can be exported from the same surface as .xyz, changed to '
                                           '.xyd for clarity, and imported here to support filtering by sounding count')

        proj_list = [str(i) + 'N' for i in range(1, 61)]  # list of all UTM zones, 1-60N and 1-60S
        proj_list.extend([str(i) + 'S' for i in range(1, 61)])
        EPSG_list = [str(i) for i in range(32601, 32661)]  # list of EPSG codes for WGS84 UTM1-60N
        EPSG_list.extend([str(i) for i in range(32701, 32761)])  # add EPSG codes for WGS84 UTM1-60S
        self.proj_dict = dict(zip(proj_list, EPSG_list))  # save for lookup during xline UTM zone conversion with pyproj
        self.ref_proj_cbox = ComboBox(proj_list, 80, 20, 'ref_proj_cbox', 'Select the reference surface UTM projection')
        ref_cbox_layout = BoxLayout([Label('Proj.:', 50, 20, 'ref_proj_lbl', (Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)),
                                     self.ref_proj_cbox], 'h')
        # ref_cbox_layout.addStretch()
        
        # Add Export All to GeoTIFF button
        self.export_geotiff_btn = PushButton('Export All to GeoTIFF', btnw, btnh, 'export_geotiff_btn',
                                            'Export Final Surface, Depth, Uncertainty, Density, and Slope as GeoTIFF files')
        
        ref_btn_layout = BoxLayout([self.add_ref_surf_btn, self.add_dens_surf_btn, ref_cbox_layout, self.export_geotiff_btn], 'v')
        ref_utm_gb = GroupBox('Reference Surface', ref_btn_layout, False, False, 'ref_surf_gb')

        # add crossline file control buttons and file list
        self.add_file_btn = PushButton('Add Crosslines', btnw, btnh, 'add_xlines_btn', 'Add crossline files')
        self.get_indir_btn = PushButton('Add Directory', btnw, btnh, 'get_indir_btn', 'Add a directory')
        self.include_subdir_chk = CheckBox('Incl. subfolders', False, 'include_subdir_chk',
                                           'Include subdirectories when adding a directory')
        self.show_path_chk = CheckBox('Show file paths', False, 'show_paths_chk', 'Show file paths')
        self.get_outdir_btn = PushButton('Select Output Dir.', btnw, btnh, 'get_outdir_btn',
                                         'Select the output directory (see current directory below)')
        source_btn_layout = BoxLayout([self.add_file_btn, self.get_indir_btn, self.get_outdir_btn, 
                                       self.include_subdir_chk, self.show_path_chk], 'v')
        source_btn_gb = GroupBox('Crosslines', source_btn_layout, False, False, 'source_btn_gb')

        # add file management buttons
        self.rmv_file_btn = PushButton('Remove Selected', btnw, btnh, 'rmv_file_btn', 'Remove selected files')
        self.clr_file_btn = PushButton('Remove All Files', btnw, btnh, 'clr_file_btn', 'Remove all files')
        file_mgmt_layout = BoxLayout([self.rmv_file_btn, self.clr_file_btn], 'v')
        file_mgmt_gb = GroupBox('File Management', file_mgmt_layout, False, False, 'file_mgmt_gb')

        # add tide file control buttons
        self.add_tide_btn = PushButton('Add Tide', btnw, btnh, 'add_tide_btn',
                                       'Add a tide (.tid) text file to apply to the accuracy crosslines.\n\n'
                                       'Each line of the tide file is space-delimited with'
                                       '[YYYY/MM/DD hh:mm:ss amplitude] in meters, positive up.\n\n'
                                       'The time zone is assumed to match that used in the accuracy crossline files '
                                       'and the vertical datum is assumed to match that used during processing of the '
                                       'reference surface.')

        tide_btn_gb = GroupBox('Tide', BoxLayout([self.add_tide_btn], 'v'), False, False, 'tide_btn_gb')

        self.calc_accuracy_btn = PushButton('Calc Accuracy', btnw, btnh, 'calc_accuracy_btn',
                                            'Calculate accuracy from loaded files')
        self.save_all_plots_btn = PushButton('Save All Plots', btnw, btnh, 'save_all_plots_btn', 'Save all plot tabs')
        
        plot_btn_gb = GroupBox('Plot Data', BoxLayout([self.calc_accuracy_btn, self.save_all_plots_btn], 'v'),
                               False, False, 'plot_btn_gb')
        # plot_btn_gb.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.MinimumExpanding)
        

        
        file_btn_layout = BoxLayout([ref_utm_gb, source_btn_gb, tide_btn_gb, file_mgmt_gb, plot_btn_gb], 'v', add_stretch=True)
        self.file_list = FileList()  # add file list with extended selection and icon size = (0,0) to avoid indent
        self.file_list.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.file_list.setMaximumWidth(350)  # Limit left panel width
        # file_layout = BoxLayout([self.file_list, file_btn_layout], 'h')
        # Compose a vertical layout for buttons and CCOM_MAC logo
        file_btn_and_logo_layout = BoxLayout([file_btn_layout], 'v')
        ccom_mac_logo_path = os.path.join(self.media_path, 'CCOM_MAC.png')
        logo_row = QtWidgets.QHBoxLayout()
        # Add CCOM_MAC logo
        if os.path.exists(ccom_mac_logo_path):
            logo_label = QtWidgets.QLabel()
            logo_pixmap = QtGui.QPixmap(ccom_mac_logo_path)
            logo_pixmap = logo_pixmap.scaled(150, 75, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation)
            logo_label.setPixmap(logo_pixmap)
            logo_label.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignBottom)
            logo_row.addWidget(logo_label)
        else:
            print(f"Warning: Logo file not found at {ccom_mac_logo_path}")
        # Add the logo row to the vertical layout, right-aligned
        file_btn_and_logo_layout.addLayout(logo_row)
        # Now use this vertical layout in the Sources group box
        file_gb = GroupBox('Sources', BoxLayout([self.file_list, file_btn_and_logo_layout], 'h'), False, False, 'file_gb')
        file_gb.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.MinimumExpanding)

        # add activity log widget
        self.log = TextEdit("background-color: lightgray", True, 'log')
        self.log.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.Minimum)
        # self.log.setMaximumWidth(500)

        # self.log.setMaximumSize(350,50)
        # self.log.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)

        update_log(self, '*** New swath accuracy processing log ***')
        log_gb = GroupBox('Activity Log', BoxLayout([self.log], 'v'), False, False, 'log_gb')

        # add progress bar for total file list and layout
        self.current_file_lbl = Label('Current File:')
        self.current_outdir_lbl = Label('Current output directory:\n' + self.output_dir)
        calc_pb_lbl = Label('Total Progress:')
        self.calc_pb = QtWidgets.QProgressBar()
        self.calc_pb.setGeometry(0, 0, 100, 30)
        self.calc_pb.setMaximum(100)  # this will update with number of files
        self.calc_pb.setValue(0)
        self.calc_pb.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.Minimum)
        # self.calc_pb.setMaximumSize(300,30)
        # self.calc_pb.setMaximumWidth(500)
        calc_pb_layout = BoxLayout([calc_pb_lbl, self.calc_pb], 'h')
        self.prog_layout = BoxLayout([self.current_file_lbl, self.current_outdir_lbl], 'v')
        self.prog_layout.addLayout(calc_pb_layout)

        # set the left panel layout with file controls on top and log on bottom
        self.left_layout = BoxLayout([file_gb, log_gb, self.prog_layout], 'v')

    def set_center_layout(self):  # set center layout with swath coverage plot
        # add figure instance and layout for swath accuracy plots
        self.swath_canvas_height = 10
        self.swath_canvas_width = 10
        self.swath_figure = Figure(figsize=(self.swath_canvas_width, self.swath_canvas_height))
        self.swath_canvas = FigureCanvas(self.swath_figure)  # canvas widget that displays the figure
        self.swath_canvas.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.MinimumExpanding)
        self.swath_toolbar = NavigationToolbar(self.swath_canvas, self) # swath plot toolbar
        self.x_max = 0.0
        self.y_max = 0.0
        self.swath_layout = BoxLayout([self.swath_toolbar, self.swath_canvas], 'v')

        # add figure instance and layout for reference surface plots
        self.surf_canvas_height = 10
        self.surf_canvas_width = 10
        self.surf_figure = Figure(figsize=(self.surf_canvas_width, self.surf_canvas_height))
        self.surf_canvas = FigureCanvas(self.surf_figure)
        self.surf_canvas.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.MinimumExpanding)
        self.surf_toolbar = NavigationToolbar(self.surf_canvas, self)
        self.x_max_surf = 0.0
        self.y_max_surf = 0.0
        self.surf_layout = BoxLayout([self.surf_toolbar, self.surf_canvas], 'v')

        # add figure instance and layout for large final masked reference surface
        # self.surf_canvas_height = 10
        # self.surf_canvas_width = 10
        self.surf_final_figure = Figure(figsize=(self.surf_canvas_width, self.surf_canvas_height))
        self.surf_final_canvas = FigureCanvas(self.surf_final_figure)
        self.surf_final_canvas.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.MinimumExpanding)
        self.surf_final_toolbar = NavigationToolbar(self.surf_final_canvas, self)
        # self.x_max_surf = 0.0
        # self.y_max_surf = 0.0
        self.surf_final_layout = BoxLayout([self.surf_final_toolbar, self.surf_final_canvas], 'v')

        # add figure instance and layout for tide plot
        self.tide_canvas_height = 10
        self.tide_canvas_width = 10
        self.tide_figure = Figure(figsize=(self.tide_canvas_width, self.tide_canvas_height))
        self.tide_canvas = FigureCanvas(self.tide_figure)
        self.tide_canvas.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.MinimumExpanding)
        self.tide_toolbar = NavigationToolbar(self.tide_canvas, self)
        self.x_max_tide = 0.0
        self.y_max_tide = 0.0
        self.tide_layout = BoxLayout([self.tide_toolbar, self.tide_canvas], 'v')

        # set up tabs
        self.plot_tabs = QtWidgets.QTabWidget()
        self.plot_tabs.setStyleSheet("background-color: none")
        self.plot_tabs.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum, QtWidgets.QSizePolicy.Policy.Maximum)

        # set up tab 1: accuracy results
        self.plot_tab1 = QtWidgets.QWidget()
        self.plot_tab1.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum, QtWidgets.QSizePolicy.Policy.Maximum)
        self.plot_tab1.layout = self.swath_layout
        self.plot_tab1.setLayout(self.plot_tab1.layout)

        # set up tab 2: reference surface
        self.plot_tab2 = QtWidgets.QWidget()
        try:
            self.plot_tab2.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum, QtWidgets.QSizePolicy.Policy.Maximum)
        except AttributeError:
            self.plot_tab2.setSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        self.plot_tab2.layout = self.surf_layout
        self.plot_tab2.setLayout(self.plot_tab2.layout)

        # set up tab 3: final masked reference surface
        self.plot_tab3 = QtWidgets.QWidget()
        try:
            self.plot_tab3.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum, QtWidgets.QSizePolicy.Policy.Maximum)
        except AttributeError:
            self.plot_tab3.setSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        self.plot_tab3.layout = self.surf_final_layout
        self.plot_tab3.setLayout(self.plot_tab3.layout)

        # NEW: set up tab 4: Ref Density Final
        self.plot_tab_density_final = QtWidgets.QWidget()
        self.density_final_figure = Figure(figsize=(self.surf_canvas_width, self.surf_canvas_height))
        self.density_final_canvas = FigureCanvas(self.density_final_figure)
        self.density_final_canvas.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.MinimumExpanding)
        self.density_final_toolbar = NavigationToolbar(self.density_final_canvas, self)
        self.density_final_layout = BoxLayout([self.density_final_toolbar, self.density_final_canvas], 'v')
        self.plot_tab_density_final.setLayout(self.density_final_layout)

        # NEW: set up tab 5: Ref Slope Final
        self.plot_tab_slope_final = QtWidgets.QWidget()
        self.slope_final_figure = Figure(figsize=(self.surf_canvas_width, self.surf_canvas_height))
        self.slope_final_canvas = FigureCanvas(self.slope_final_figure)
        self.slope_final_canvas.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.MinimumExpanding)
        self.slope_final_toolbar = NavigationToolbar(self.slope_final_canvas, self)
        self.slope_final_layout = BoxLayout([self.slope_final_toolbar, self.slope_final_canvas], 'v')
        self.plot_tab_slope_final.setLayout(self.slope_final_layout)

        # NEW: set up tab 4: Depth
        self.plot_tab_depth = QtWidgets.QWidget()
        self.depth_figure = Figure(figsize=(self.surf_canvas_width, self.surf_canvas_height))
        self.depth_canvas = FigureCanvas(self.depth_figure)
        self.depth_canvas.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.MinimumExpanding)
        self.depth_toolbar = NavigationToolbar(self.depth_canvas, self)
        self.depth_layout = BoxLayout([self.depth_toolbar, self.depth_canvas], 'v')
        self.plot_tab_depth.setLayout(self.depth_layout)

        # NEW: set up tab 6: Ref Uncertainty
        self.plot_tab_uncertainty_final = QtWidgets.QWidget()
        self.uncertainty_final_figure = Figure(figsize=(self.surf_canvas_width, self.surf_canvas_height))
        self.uncertainty_final_canvas = FigureCanvas(self.uncertainty_final_figure)
        self.uncertainty_final_canvas.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.MinimumExpanding)
        self.uncertainty_final_toolbar = NavigationToolbar(self.uncertainty_final_canvas, self)
        self.uncertainty_final_layout = BoxLayout([self.uncertainty_final_toolbar, self.uncertainty_final_canvas], 'v')
        self.plot_tab_uncertainty_final.setLayout(self.uncertainty_final_layout)

        # set up tab 4: crossline tide (now will be tab 9)
        self.plot_tab4 = QtWidgets.QWidget()
        try:
            self.plot_tab4.setSizePolicy(QtWidgets.QSizePolicy.Policy.Maximum, QtWidgets.QSizePolicy.Policy.Maximum)
        except AttributeError:
            self.plot_tab4.setSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        self.plot_tab4.layout = self.tide_layout
        self.plot_tab4.setLayout(self.plot_tab4.layout)

        # NEW: set up tab 8: Soundings per bin
        self.plot_tab_soundings = QtWidgets.QWidget()
        self.soundings_figure = Figure(figsize=(self.surf_canvas_width, self.surf_canvas_height))
        self.soundings_canvas = FigureCanvas(self.soundings_figure)
        self.soundings_canvas.setSizePolicy(QtWidgets.QSizePolicy.Policy.MinimumExpanding, QtWidgets.QSizePolicy.Policy.MinimumExpanding)
        self.soundings_toolbar = NavigationToolbar(self.soundings_canvas, self)
        self.soundings_layout = BoxLayout([self.soundings_toolbar, self.soundings_canvas], 'v')
        self.plot_tab_soundings.setLayout(self.soundings_layout)



        # add tabs to tab layout
        self.plot_tabs.addTab(self.plot_tab1, 'Accuracy')
        self.plot_tabs.addTab(self.plot_tab2, 'Surface Filters')
        self.plot_tabs.addTab(self.plot_tab3, 'Final Surface')
        self.plot_tabs.addTab(self.plot_tab_depth, 'Depth')
        self.plot_tabs.addTab(self.plot_tab_uncertainty_final, 'Uncertainty')
        self.plot_tabs.addTab(self.plot_tab_density_final, 'Density')
        self.plot_tabs.addTab(self.plot_tab_slope_final, 'Slope')
        self.plot_tabs.addTab(self.plot_tab_soundings, 'Soundings')
        self.plot_tabs.addTab(self.plot_tab4, 'Tide')

        self.center_layout = BoxLayout([self.plot_tabs], 'v')
        # self.center_layout.addStretch()

    def initialize_default_filters(self):
        """Initialize all filter controls with their default values"""
        try:
            # Reference surface defaults
            self.min_depth_ref_tb.setText('0')
            self.max_depth_ref_tb.setText('10000')
            self.max_slope_tb.setText('5')
            self.min_dens_tb.setText('10')
            self.max_u_tb.setText('10')
            self.slope_win_cbox.setCurrentIndex(1)  # 3x3
            
            # Crossline defaults
            self.min_depth_xline_tb.setText('0')
            self.max_depth_xline_tb.setText('10000')
            self.min_angle_xline_tb.setText('-75')
            self.max_angle_xline_tb.setText('75')
            self.min_bs_xline_tb.setText('-50')
            self.max_bs_xline_tb.setText('0')
            self.max_dz_tb.setText('10')
            self.max_dz_wd_tb.setText('5')
            self.min_bin_count_tb.setText('10')
            self.depth_mode_cbox.setCurrentIndex(0)  # None
            
            # Point count defaults
            self.max_count_tb.setText('50000')
            self.dec_fac_tb.setText('1')
            
            # Disable all filters by default (except point count limit)
            self.depth_ref_gb.setChecked(False)
            self.slope_gb.setChecked(False)
            self.density_gb.setChecked(False)
            self.uncertainty_gb.setChecked(False)
            self.depth_xline_gb.setChecked(False)
            self.angle_xline_gb.setChecked(False)
            self.bs_xline_gb.setChecked(False)
            if hasattr(self, 'dz_abs_gb'):
                self.dz_abs_gb.setChecked(False)
            if hasattr(self, 'dz_pct_gb'):
                self.dz_pct_gb.setChecked(False)
            self.depth_mode_gb.setChecked(False)
            self.bin_count_gb.setChecked(True)
            # Keep point count limit disabled by default, use bin-based decimation instead
            self.pt_count_gb.setChecked(False)
            # Enable bin-based decimation by default for uniform coverage
            self.bin_decimation_gb.setChecked(True)
            
        except Exception as e:
            print(f"Error initializing default filters: {str(e)}")

    def initialize_default_plot_settings(self):
        """Initialize all plot settings with their default values"""
        try:
            # Custom info defaults
            self.custom_info_gb.setChecked(False)
            self.model_cbox.setCurrentIndex(0)  # First model in list
            self.ship_tb.setText('R/V Unsinkable II')
            self.cruise_tb.setText('A 3-hour tour')
            self.show_model_chk.setChecked(True)
            self.show_ship_chk.setChecked(True)
            self.show_cruise_chk.setChecked(True)
            
            # Data reference defaults
            self.ref_cbox.setCurrentIndex(0)  # Waterline
            self.tide_unit_cbox.setCurrentIndex(0)  # Meter
            self.waterline_tb.setText('0.00')
            
            # Point style defaults
            self.pt_size_cbox.setCurrentIndex(1)  # Point size 1
            self.pt_size_cov_cbox.setCurrentIndex(5)  # Point size 5 for coverage
            self.pt_alpha_acc_tb.setText('100')  # default opacity for accuracy
            self.pt_alpha_cov_tb.setText('6')  # default opacity for coverage
            
            # Plot limits defaults
            self.plot_lim_gb.setChecked(False)
            self.max_beam_angle_tb.setText(str(self.x_max_default))
            self.angle_spacing_tb.setText(str(self.x_spacing_default))
            self.max_bias_tb.setText(str(self.y_max_bias_default))
            self.max_std_tb.setText(str(self.y_max_std_default))
            self.axis_margin_tb.setText(str(self.axis_margin_default))
            self.tide_range_tb.setText('12')  # Default to 12 hours
            
            # Toggle options defaults
            self.show_acc_proc_text_chk.setChecked(False)
            self.show_ref_proc_text_chk.setChecked(False)
            self.grid_lines_toggle_chk.setChecked(True)
            self.IHO_lines_toggle_chk.setChecked(True)
            self.update_ref_plots_chk.setChecked(True)
            self.show_xline_cov_chk.setChecked(True)
            self.show_u_plot_chk.setChecked(True)
            self.show_shaded_relief_chk.setChecked(True)
            self.show_special_order_chk.setChecked(False)
            self.show_order_1a_chk.setChecked(False)
            self.show_order_1b_chk.setChecked(False)
            self.show_order_2_chk.setChecked(False)
            self.show_order_3_chk.setChecked(False)
            if hasattr(self, 'unique_line_pt_colors_chk'):
                self.unique_line_pt_colors_chk.setChecked(False)
            
            # Flatten swath defaults
            self.flatten_mean_gb.setChecked(False)
            self.set_zero_mean_chk.setChecked(False)
            self.mean_bias_angle_lim_tb.setText('40')
            
        except Exception as e:
            print(f"Error initializing default plot settings: {str(e)}")

    def set_right_layout(self):
        # set right layout with swath plot controls

        # TAB 1: PLOT OPTIONS
        # add text boxes for system, ship, cruise
        model_tb_lbl = Label('Model:', width=100, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.model_cbox = ComboBox(self.model_list, 100, 20, 'model_cbox', 'Select the model')
        self.show_model_chk = CheckBox('', True, 'show_model_chk', 'Show model in plot title')
        model_info_layout_left = BoxLayout([model_tb_lbl, self.model_cbox], 'h', alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        model_info_layout = BoxLayout([model_info_layout_left, self.show_model_chk], 'h', alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))

        ship_tb_lbl = Label('Ship Name:', width=100, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.ship_tb = LineEdit('R/V Unsinkable II', 100, 20, 'ship_tb', 'Enter the ship name')
        self.show_ship_chk = CheckBox('', True, 'show_ship_chk', 'Show ship name in plot title')
        ship_info_layout_left = BoxLayout([ship_tb_lbl, self.ship_tb], 'h', alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        ship_info_layout = BoxLayout([ship_info_layout_left, self.show_ship_chk], 'h', alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))

        cruise_tb_lbl = Label('Description:', width=100, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.cruise_tb = LineEdit('A 3-hour tour', 100, 20, 'cruise_tb', 'Enter the cruise name')
        self.show_cruise_chk = CheckBox('', True, 'show_cruise_chk', 'Show cruise in plot title')
        cruise_info_layout_left = BoxLayout([cruise_tb_lbl, self.cruise_tb], 'h', alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        cruise_info_layout = BoxLayout([cruise_info_layout_left, self.show_cruise_chk], 'h', alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))

        self.custom_info_gb = GroupBox('Use custom system information',
                                       BoxLayout([model_info_layout, ship_info_layout, cruise_info_layout], 'v'),
                                       True, False, 'custom_info_gb')
        self.custom_info_gb.setToolTip('Add system/cruise info; system info parsed from the file is used if available')

        # add depth reference options and groupbox
        self.ref_cbox = ComboBox(self.data_ref_list, 100, 20, 'ref_cbox',
                                 'Select the reference for plotting depth and acrosstrack distance\n\n'
                                 'As parsed, .all depths are referenced to the TX array and .kmall depths are '
                                 'referenced to the mapping system origin in SIS\n\n'
                                 'Waterline reference is appropriate for normal surface vessel data; '
                                 'other options are available for special cases (e.g., underwater vehicles or '
                                 'troubleshooting installation offset discrepancies)\n\n'
                                 'Overview of adjustments:\n\nWaterline: change reference to the waterline '
                                 '(.all: shift Y and Z ref from TX array to origin, then Z ref to waterline; '
                                 '.kmall: shift Z ref from origin to waterline)\n\n'
                                 'Origin: change reference to the mapping system origin '
                                 '(.all: shift Y and Z ref from TX array to origin; .kmall: no change)\n\n'
                                 'TX Array: change reference to the TX array reference point '
                                 '(.all: no change; .kmall: shift Y and Z ref from origin to TX array)\n\n'
                                 'Raw: use the native depths and acrosstrack distances parsed from the file '
                                 '(.all: referenced to TX array; .kmall: referenced to mapping system origin)')

        data_ref_lbl = Label('Reference data to:', 100, 20, 'data_ref_lbl', (Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter))
        data_ref_layout = BoxLayout([data_ref_lbl, self.ref_cbox], 'h')

        # add tide units
        self.tide_unit_cbox = ComboBox([unit for unit in self.tide_unit_dict.keys()], 100, 20, 'tide_unit_cbox',
                                 'Select the tide amplitude units')
        tide_unit_lbl = Label('Tide file units:', 100, 20, 'tide_unit_lbl', (Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter))
        tide_unit_layout = BoxLayout([tide_unit_lbl, self.tide_unit_cbox], 'h')


        # add waterline offset (SIS format, with WL in meters Z+ down from origin)
        waterline_lbl = Label('Adjust waterline (m, pos. down):', width=110, alignment=(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter))
        self.waterline_tb = LineEdit('0.00', 50, 20, 'waterline_tb',
                                     'Adjust the SIS waterline (WL) value applied.  WL is the Z offset from the system '
                                     'origin to waterline in meters, positive DOWN (Kongsberg convention).\n\n'
                                     'The value entered here will be added to the WL parsed from the crossline.  For '
                                     'instance, if a crossline is acquired with WL = -3.00 but WL is actually -3.50 m, '
                                     '(WL is 0.50 m higher above the origin, or more negative, than the configured '
                                     'WL), then a WL adjustment of -0.50 can be entered to effectively shift the '
                                     'crossline WL from -3.00 to -3.50 m for processing purposes.  This will cause the '
                                     'crossline data to appear deeper by 0.50 m compared to default as acquired.\n\n'
                                     'NOTE: Waterline adjustment applies to crosslines only.  No change is made to the '
                                     'reference surface.')
        self.waterline_tb.setValidator(QDoubleValidator(-1*np.inf, np.inf, 2))
        waterline_layout = BoxLayout([waterline_lbl, self.waterline_tb], 'h')

        # ref_unit_layout = BoxLayout([data_ref_layout, tide_unit_layout], 'v')
        ref_unit_layout = BoxLayout([data_ref_layout, tide_unit_layout, waterline_layout], 'v')

        # self.data_ref_gb = GroupBox('Data reference', data_ref_layout, False, False, 'data_ref_gb')
        self.data_ref_gb = GroupBox('Data reference and units', ref_unit_layout, False, False, 'data_ref_gb')

        # add point size and opacity comboboxes
        # pt_size_lbl = Label('Point size:', width=60, alignment=(Qt.AlignRight | Qt.AlignVCenter))
        # self.pt_size_cbox = ComboBox([str(pt) for pt in range(11)], 45, 20, 'pt_size_cbox', 'Select point size')
        # self.pt_size_cbox.setCurrentIndex(5)
        # pt_size_layout = BoxLayout([pt_size_lbl, self.pt_size_cbox], 'h', add_stretch=True)

        # Accuracy groupbox
        acc_pt_size_lbl = Label('Point size:', width=60, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.pt_size_cbox = ComboBox([str(pt) for pt in range(11)], 45, 20, 'pt_size_cbox',
                                     'Select point size for soundings in the accuracy plot')
        self.pt_size_cbox.setCurrentIndex(1)
        acc_pt_size_layout = BoxLayout([acc_pt_size_lbl, self.pt_size_cbox], 'h')
        
        acc_pt_alpha_lbl = Label('Opacity (%):', width=60, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.pt_alpha_acc_tb = LineEdit('100', 45, 20, 'pt_alpha_acc_tb', 'Set opacity for accuracy plot (0-100)')
        self.pt_alpha_acc_tb.setValidator(QDoubleValidator(0, 100, 1))
        acc_pt_alpha_layout = BoxLayout([acc_pt_alpha_lbl, self.pt_alpha_acc_tb], 'h', add_stretch=True)
        
        # Add checkbox for unique line point colors
        self.unique_line_pt_colors_chk = CheckBox('Unique Line Colors', False, 'unique_line_pt_colors_chk',
                                                  'Color points in the accuracy plot differently based on which crossline file they come from')
        
        acc_groupbox_layout = BoxLayout([acc_pt_size_layout, acc_pt_alpha_layout, self.unique_line_pt_colors_chk], 'v')
        acc_groupbox = GroupBox('Accuracy', acc_groupbox_layout, False, False, 'acc_groupbox')
        
        # Surface Coverage groupbox
        cov_pt_size_lbl = Label('Point size:', width=60, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.pt_size_cov_cbox = ComboBox([str(pt) for pt in range(11)], 45, 20, 'pt_size_cov_cbox',
                                         'Select point size for soundings in the coverage plot (if shown)')
        self.pt_size_cov_cbox.setCurrentIndex(5)
        cov_pt_size_layout = BoxLayout([cov_pt_size_lbl, self.pt_size_cov_cbox], 'h')
        
        cov_pt_alpha_lbl = Label('Opacity (%):', width=60, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.pt_alpha_cov_tb = LineEdit('6', 45, 20, 'pt_alpha_cov_tb', 'Set opacity for coverage plot (0-100)')
        self.pt_alpha_cov_tb.setValidator(QDoubleValidator(0, 100, 1))
        cov_pt_alpha_layout = BoxLayout([cov_pt_alpha_lbl, self.pt_alpha_cov_tb], 'h', add_stretch=True)
        
        cov_groupbox_layout = BoxLayout([cov_pt_size_layout, cov_pt_alpha_layout], 'v')
        cov_groupbox = GroupBox('Surface Coverage', cov_groupbox_layout, False, False, 'cov_groupbox')
        
        # set final point parameter layout with two groupboxes side by side
        pt_param_layout = BoxLayout([acc_groupbox, cov_groupbox], 'h')
        pt_param_gb = GroupBox('Point Styles', pt_param_layout, False, False, 'pt_param_gb')

        # add custom plot axis limits
        max_beam_angle_lbl = Label('Max Angle (deg):', width=110, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.max_beam_angle_tb = LineEdit('', 50, 20, 'max_beam_angle_tb', 'Set the maximum plot angle (X axis)')
        self.max_beam_angle_tb.setValidator(QDoubleValidator(0, 90, 2))
        max_beam_angle_layout = BoxLayout([max_beam_angle_lbl, self.max_beam_angle_tb], 'h')

        angle_spacing_lbl = Label('Label (deg):', width=110, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.angle_spacing_tb = LineEdit('', 50, 20, 'angle_spacing_tb', 'Set the angle tick spacing')
        self.angle_spacing_tb.setValidator(QDoubleValidator(0, 90, 2))
        angle_spacing_layout = BoxLayout([angle_spacing_lbl, self.angle_spacing_tb], 'h')
        
        # Combine Max Angle and Label on the same line
        angle_combined_layout = BoxLayout([max_beam_angle_layout, angle_spacing_layout], 'h', add_stretch=True)

        max_bias_lbl = Label('Max Bias ' + self.unit_mode + ':', width=110,
                             alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.max_bias_tb = LineEdit('', 50, 20, 'max_bias_tb', 'Set the maximum plot bias (Y axis)')
        self.max_bias_tb.setValidator(QDoubleValidator(0, 100, 2))
        max_bias_layout = BoxLayout([max_bias_lbl, self.max_bias_tb], 'h')

        max_std_lbl = Label('Max SDEV ' + self.unit_mode + ':', width=110,
                            alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.max_std_tb = LineEdit('', 50, 20, 'max_std_tb', 'Set the maximum plot standard deviation (Y axis)')
        self.max_std_tb.setValidator(QDoubleValidator(0, 100, 2))
        max_std_layout = BoxLayout([max_std_lbl, self.max_std_tb], 'h')
        
        # Combine Max Bias and Max STDEV on the same line
        bias_std_combined_layout = BoxLayout([max_bias_layout, max_std_layout], 'h', add_stretch=True)

        axis_margin_lbl = Label('Plot Buffer (%)', width=110, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.axis_margin_tb = LineEdit('', 50, 20, 'axis_margin_tb', 'Set the reference plot axis margins (%)')
        self.axis_margin_tb.setValidator(QDoubleValidator(0, 100, 2))
        axis_margin_layout = BoxLayout([axis_margin_lbl, self.axis_margin_tb], 'h')

        # Add tide range parameter
        tide_range_lbl = Label('Tide Range (hrs)', width=110, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.tide_range_tb = LineEdit('12', 50, 20, 'tide_range_tb', 'Set the ± time range in hours for tide plotting')
        self.tide_range_tb.setValidator(QDoubleValidator(0.1, 72, 1))  # Allow 0.1 to 72 hours
        tide_range_layout = BoxLayout([tide_range_lbl, self.tide_range_tb], 'h')
        
        # Combine Plot Buffer and Tide Range on the same line
        buffer_tide_combined_layout = BoxLayout([axis_margin_layout, tide_range_layout], 'h', add_stretch=True)

        # autoscale_lbl = Label('Autoscale')

        # plot_lim_layout = BoxLayout([max_beam_angle_layout, angle_spacing_layout, max_bias_layout, max_std_layout], 'v')
        plot_lim_layout = BoxLayout([angle_combined_layout, bias_std_combined_layout,
                                     buffer_tide_combined_layout], 'v')

        self.plot_lim_gb = GroupBox('Use custom plot limits', plot_lim_layout, True, False, 'plot_lim_gb')

        # add check boxes with other options
        self.show_acc_proc_text_chk = CheckBox('Show crossline proc. params.', False, 'show_acc_proc_text_chk',
                                               'Show text box with crossline processing/filtering information')
        self.show_ref_proc_text_chk = CheckBox('Show reference proc. params.', False, 'show_ref_proc_text_chk',
                                               'Show text box with reference surface processing/filtering information')
        self.grid_lines_toggle_chk = CheckBox('Show grid lines', True, 'show_grid_lines_chk', 'Show grid lines')
        self.IHO_lines_toggle_chk = CheckBox('Show IHO lines', True, 'show_IHO_lines_chk', 'Show IHO lines')
        self.update_ref_plots_chk = CheckBox('Show filtered reference data', True, 'show_filtered_ref_chk',
                                             'Update reference surface subplots with depth, density, and slope filters'
                                             'applied (uncheck to show the reference surface data as parsed).\n\n'
                                             'This option affects plotting only; all reference surface filters are '
                                             'applied prior to crossline analysis.')
        self.show_xline_cov_chk = CheckBox('Show crossline coverage', True, 'show_xline_cov_chk',
                                           'Show crossline soundings (unfiltered) on reference grids')
        self.show_u_plot_chk = CheckBox('Show uncertainty plot if parsed', True, 'show_u_plot_chk',
                                        'Plot the reference surface uncertainty if parsed from .xyz file (zeros if not'
                                        'parsed.\nThis will replace the subplot for the "final" masked surface.')

        self.show_shaded_relief_chk = CheckBox('Show shaded relief', True, 'show_shaded_relief_chk',
                                               'Add shaded relief visualization to reference surface plots using '
                                               'multidirectional hillshade technique for enhanced terrain visualization.')

        # Bathymetric order checkboxes
        self.show_special_order_chk = CheckBox('Show Special Order (±0.25+0.0075×depth m)', False, 'show_special_order_chk',
                                               'Show Special Order accuracy limits (±0.25+0.0075×depth m) on accuracy plot')
        self.show_order_1a_chk = CheckBox('Show Order 1a (±0.5+0.013×depth m)', False, 'show_order_1a_chk',
                                           'Show Order 1a accuracy limits (±0.5+0.013×depth m) on accuracy plot')
        self.show_order_1b_chk = CheckBox('Show Order 1b (±0.5+0.013×depth m)', False, 'show_order_1b_chk',
                                           'Show Order 1b accuracy limits (±0.5+0.013×depth m) on accuracy plot')
        self.show_order_2_chk = CheckBox('Show Order 2 (±1.0+0.023×depth m)', False, 'show_order_2_chk',
                                         'Show Order 2 accuracy limits (±1.0+0.023×depth m) on accuracy plot')
        self.show_order_3_chk = CheckBox('Show Order 3 (±2.0+0.05×depth m)', False, 'show_order_3_chk',
                                         'Show Order 3 accuracy limits (±2.0+0.05×depth m) on accuracy plot')

        toggle_chk_layout = BoxLayout([self.show_acc_proc_text_chk,
                                       self.show_ref_proc_text_chk,
                                       self.grid_lines_toggle_chk,
                                       self.update_ref_plots_chk,
                                       self.show_xline_cov_chk,
                                       self.show_u_plot_chk,
                                       self.show_shaded_relief_chk,
                                       self.show_special_order_chk,
                                       self.show_order_1a_chk,
                                       self.show_order_1b_chk,
                                       self.show_order_2_chk,
                                       self.show_order_3_chk], 'v')  # self.IHO_lines_toggle_chk], 'v')

        toggle_gb = QtWidgets.QGroupBox('Plot options')
        toggle_gb.setLayout(toggle_chk_layout)

        # options to flatten the mean bias curve to reduce impacts of refraction on the visualization of sounding dist.
        self.set_zero_mean_chk = CheckBox('Force zero mean (remove all biases)', False, 'zero_mean_plot_chk',
                                        'Force the mean to zero for each angle bin; this is to be used only for '
                                        'visualizing the distribution of soundings without other biases (e.g., '
                                        'refraction issues) and does not represent the observed swath performance.')

        mean_bias_angle_lim_lbl = Label('Angle limit for bias calc. (deg)',
                                        width=110, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))

        self.mean_bias_angle_lim_tb = LineEdit('40', 45, 20, 'mean_bias_angle_lim_tb',
                                               'Set the angle limit (+/- deg to each side) for the desired portion of '
                                               'the swath to use for the mean bias calculations; this is useful for '
                                               'reducing the impacts of significant refraction (e.g., outer swath) on '
                                               'visualization of the sounding distribition, thereby highlighting other '
                                               'biases (e.g., waterline errors) that may have been masked in part by '
                                               'the refraction issues.')

        self.mean_bias_angle_lim_tb.setValidator(QDoubleValidator(1, np.inf, 2))
        mean_bias_angle_lim_layout = BoxLayout([mean_bias_angle_lim_lbl, self.mean_bias_angle_lim_tb], 'h')

        flatten_mean_chk_layout = BoxLayout([mean_bias_angle_lim_layout, self.set_zero_mean_chk], 'v')
        self.flatten_mean_gb = GroupBox('Flatten swath', flatten_mean_chk_layout, True, False, 'flatten_mean_gb')
        # flatten_mean_gb.setLayout(flatten_mean_chk_layout)


        # TAB 2: FILTER OPTIONS: REFERENCE SURFACE
        # add custom depth limits for ref surf
        min_depth_ref_lbl = Label('Min:', alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        max_depth_ref_lbl = Label('Max:', alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))

        self.min_depth_ref_tb = LineEdit('0', 40, 20, 'min_depth_ref_tb', 'Min depth of the reference surface data')
        self.max_depth_ref_tb = LineEdit('10000', 40, 20, 'max_depth_ref_tb', 'Max depth of the reference surface data')
        self.min_depth_ref_tb.setValidator(QDoubleValidator(0, float(self.max_depth_ref_tb.text()), 2))
        self.max_depth_ref_tb.setValidator(QDoubleValidator(float(self.min_depth_ref_tb.text()), np.inf, 2))
        depth_ref_layout = BoxLayout([min_depth_ref_lbl, self.min_depth_ref_tb,
                                      max_depth_ref_lbl, self.max_depth_ref_tb], 'h', add_stretch=True)
        self.depth_ref_gb = GroupBox('Depth (m, ref. surf.)', depth_ref_layout, True, False, 'depth_ref_gb')
        self.depth_ref_gb.setToolTip('Hide reference surface data by depth (m, positive down).\n\n'
                                     'Acceptable min/max fall within [0 inf].')

        # add slope calc and filtering options
        slope_win_lbl = Label('Window (cells):', width=50, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.slope_win_cbox = ComboBox(['1x1', '3x3', '5x5', '7x7', '9x9'], 60, 20, 'slope_window_cbox',
                                       'Select the ref. surface moving average window size for slope calculation')
        max_slope_lbl = Label('Max (deg):', width=50, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.max_slope_tb = LineEdit('5', 40, 20, 'max_slope_tb',
                                     'Set the maximum reference surface slope allowed for crossline analysis (maximum '
                                     'slope for each ref. surface node is estimated from maximum N-S and E-W gradients '
                                     'of the reference depth grid after any averaging using the window size selected)')
        self.max_slope_tb.setValidator(QDoubleValidator(0, np.inf, 2))

        slope_win_layout = BoxLayout([slope_win_lbl, self.slope_win_cbox], 'h', add_stretch=True)
        slope_max_layout = BoxLayout([max_slope_lbl, self.max_slope_tb], 'h', add_stretch=True)
        slope_layout = BoxLayout([slope_win_layout, slope_max_layout], 'v')
        self.slope_gb = GroupBox('Slope', slope_layout, True, False, 'slope_win_gb')

        # ref surf sounding density filtering
        min_dens_lbl = Label('Min:', width=50, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.min_dens_tb = LineEdit('10', 40, 20, 'min_dens_tb',
                                    'Set the min. reference surface sounding density allowed for crossline analysis)')
        self.min_dens_tb.setValidator(QDoubleValidator(0, np.inf, 2))

        density_layout = BoxLayout([min_dens_lbl, self.min_dens_tb], 'h', add_stretch=True)
        self.density_gb = GroupBox('Density (soundings/cell)', density_layout, True, False, 'density_gb')

        # ref surf uncertainty filtering
        max_u_lbl = Label('Max:', width=50, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.max_u_tb = LineEdit('10', 40, 20, 'max_u_tb',
                                 'Set the max. reference surface uncertainty allowed for crossline analysis)')
        self.max_u_tb.setValidator(QDoubleValidator(0, np.inf, 2))

        uncertainty_layout = BoxLayout([max_u_lbl, self.max_u_tb], 'h', add_stretch=True)
        self.uncertainty_gb = GroupBox('Uncertainty (m)', uncertainty_layout, True, False, 'uncertainty_gb')

        # set up layout and groupbox for tabs
        tab2_ref_filter_layout = BoxLayout([self.depth_ref_gb, self.uncertainty_gb,
                                            self.density_gb, self.slope_gb], 'v')
        tab2_ref_filter_gb = GroupBox('Reference surface', tab2_ref_filter_layout, False, False, 'tab2_ref_filter_gb')

        # TAB 2: FILTER OPTIONS: ACCURACY CROSSLINES
        # add custom depth limits for crosslines
        min_depth_xline_lbl = Label('Min:', alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        max_depth_xline_lbl = Label('Max:', alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.min_depth_xline_tb = LineEdit('0', 40, 20, 'min_depth_tb', 'Min depth of the crossline data')
        self.max_depth_xline_tb = LineEdit('10000', 40, 20, 'max_depth_tb', 'Max depth of the crossline data')
        self.min_depth_xline_tb.setValidator(QDoubleValidator(0, float(self.max_depth_xline_tb.text()), 2))
        self.max_depth_xline_tb.setValidator(QDoubleValidator(float(self.min_depth_xline_tb.text()), np.inf, 2))
        depth_xline_layout = BoxLayout([min_depth_xline_lbl, self.min_depth_xline_tb,
                                        max_depth_xline_lbl, self.max_depth_xline_tb], 'h')
        depth_xline_layout.addStretch()
        self.depth_xline_gb = GroupBox('Depth (m, crosslines)', depth_xline_layout, True, False, 'depth_xline_gb')
        self.depth_xline_gb.setToolTip('Hide crossline data by depth (m, positive down).\n\n'
                                       'Acceptable min/max fall within [0 inf].')

        # add custom swath angle limits (-port, +stbd on [-inf, inf]; not [0, inf] like swath coverage plotter)
        min_angle_xline_lbl = Label('Min:', width=50, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        max_angle_xline_lbl = Label('Max:', width=50, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.min_angle_xline_tb = LineEdit('-75', 40, 20, 'min_angle_xline_tb',
                                           'Set the minimum beam angle for accuracy calculations (port, <= max angle)')
        self.max_angle_xline_tb = LineEdit('75', 40, 20, 'max_angle_xline_tb',
                                           'Set the maximum beam angle for accuracy calulations (stbd, >= min angle)')
        self.min_angle_xline_tb.setValidator(QDoubleValidator(-1*np.inf, float(self.max_angle_xline_tb.text()), 2))
        self.max_angle_xline_tb.setValidator(QDoubleValidator(float(self.min_angle_xline_tb.text()), np.inf, 2))
        # angle_layout = BoxLayout([min_angle_lbl, self.min_angle_tb, max_angle_lbl, self.max_angle_tb], 'h')
        angle_xline_layout = BoxLayout([min_angle_xline_lbl, self.min_angle_xline_tb,
                                        max_angle_xline_lbl, self.max_angle_xline_tb],
                                       'h', add_stretch=True)

        self.angle_xline_gb = GroupBox('Angle (deg)', angle_xline_layout, True, False, 'angle_xline_gb')
        self.angle_xline_gb.setToolTip('Hide soundings based on nominal swath angles calculated from depths and '
                                 'acrosstrack distances; these swath angles may differ slightly from RX beam '
                                 'angles (w.r.t. RX array) due to installation, attitude, and refraction.\n\n'
                                 'Angles are treated as negative to port and positive to starboard.')

        # add custom reported backscatter limits
        min_bs_xline_lbl = Label('Min:', width=50, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.min_bs_xline_tb = LineEdit('-50', 40, 20, 'min_bs_xline_tb',
                                        'Set the minimum reported backscatter (e.g., -50 dB); '
                                        'while backscatter values in dB are inherently negative, the filter range may '
                                        'include positive values to accommodate anomalous reported backscatter data')
        max_bs_xline_lbl = Label('Max:', width=50, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.max_bs_xline_tb = LineEdit('0', 40, 20, 'max_bs_xline_tb',
                                  'Set the maximum reported backscatter of the data (e.g., 0 dB); '
                                  'while backscatter values in dB are inherently negative, the filter range may '
                                  'include positive values to accommodate anomalous reported backscatter data')
        self.min_bs_xline_tb.setValidator(QDoubleValidator(-1 * np.inf, float(self.max_bs_xline_tb.text()), 2))
        self.max_bs_xline_tb.setValidator(QDoubleValidator(float(self.min_bs_xline_tb.text()), np.inf, 2))
        bs_xline_layout = BoxLayout([min_bs_xline_lbl, self.min_bs_xline_tb,
                               max_bs_xline_lbl, self.max_bs_xline_tb], 'h', add_stretch=True)
        self.bs_xline_gb = GroupBox('Backscatter (dB)', bs_xline_layout, True, False, 'bs_xline_gb')
        self.bs_xline_gb.setToolTip('Hide data by reported backscatter amplitude (dB).\n\n'
                                    'Acceptable min/max fall within [-inf inf] to accommodate anomalous data >0.')

        # add depth mode filter for crossline
        # self.depth_mode_list = ['Very Shallow', 'Shallow', 'Medium', 'Deep', 'Deeper',
        #                         'Very Deep', 'Extra Deep', 'Extreme Deep']
        self.depth_mode_list = ['None']
        self.depth_mode_cbox = ComboBox(self.depth_mode_list, 80, 20, 'depth_mode_cbox',
                                        'Filter crosslines by depth mode(s) found in file(s)')
        depth_mode_layout = BoxLayout([Label('Available mode(s):', 50, 20,'depth_mode_lbl',
                                             (Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)), self.depth_mode_cbox], 'h')
        self.depth_mode_gb = GroupBox('Depth mode', depth_mode_layout, True, False, 'depth_mode_gb')


        # add minimum sounding count for binning
        min_bin_count_lbl = Label('Min:', width=50, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.min_bin_count_tb = LineEdit('10', 40, 20, 'min_bin_count_tb',
                                         'Set the min. sounding count per angle bin for running accuracy calculations)')
        self.min_bin_count_tb.setValidator(QDoubleValidator(0, np.inf, 2))

        # add separate absolute depth difference filtering
        max_dz_lbl = Label('Max (m):', width=50, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.max_dz_tb = LineEdit('10', 40, 20, 'max_dz_tb',
                                  'Set the maximum absolute depth difference in meters (>0)')
        self.max_dz_tb.setValidator(QDoubleValidator(0.001, np.inf, 3))
        
        dz_abs_layout = BoxLayout([max_dz_lbl, self.max_dz_tb], 'h', add_stretch=True)
        self.dz_abs_gb = GroupBox('Depth difference (m)', dz_abs_layout, True, False, 'dz_abs_gb')
        self.dz_abs_gb.setToolTip('Hide crossline data by absolute depth difference from reference surface (>0 m)')

        # add separate percentage of water depth filtering
        max_dz_wd_lbl = Label('Max (%WD):', width=50, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.max_dz_wd_tb = LineEdit('5', 40, 20, 'max_dz_wd_tb',
                                     'Set the maximum depth difference as percentage of water depth (1-99%)')
        self.max_dz_wd_tb.setValidator(QDoubleValidator(1, 99, 0))
        
        dz_pct_layout = BoxLayout([max_dz_wd_lbl, self.max_dz_wd_tb], 'h', add_stretch=True)
        self.dz_pct_gb = GroupBox('Depth difference (%WD)', dz_pct_layout, True, False, 'dz_pct_gb')
        self.dz_pct_gb.setToolTip('Hide crossline data by depth difference as percentage of water depth (1-99%)')

        bin_count_layout = BoxLayout([min_bin_count_lbl, self.min_bin_count_tb], 'h', add_stretch=True)
        self.bin_count_gb = GroupBox('Bin count (soundings/bin)', bin_count_layout, True, True, 'bin_count_gb')

        # add plotted point max count and decimation factor control in checkable groupbox
        max_count_lbl = Label('Max. plotted points (0-inf):', width=140, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.max_count_tb = LineEdit(str(self.n_points_max_default), 50, 20, 'max_count_tb',
                                     'Set the maximum number of plotted points for each data set')
        self.max_count_tb.setValidator(QDoubleValidator(0, np.inf, 2))
        max_count_layout = BoxLayout([max_count_lbl, self.max_count_tb], 'h')
        dec_fac_lbl = Label('Decimation factor (1-inf):', width=140, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.dec_fac_tb = LineEdit(str(self.dec_fac_default), 50, 20, 'dec_fac_tb', 'Set the custom decimation factor')
        self.dec_fac_tb.setValidator(QDoubleValidator(1, np.inf, 2))
        dec_fac_layout = BoxLayout([dec_fac_lbl, self.dec_fac_tb], 'h')
        pt_count_layout = BoxLayout([max_count_layout, dec_fac_layout], 'v')
        self.pt_count_gb = GroupBox('Limit plotted point count (plot faster)', pt_count_layout, True, False, 'pt_ct_gb')
        
        # add bin-based decimation controls
        max_points_per_bin_lbl = Label('Max points per beam angle bin:', width=180, alignment=(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter))
        self.max_points_per_bin_tb = LineEdit('350', 50, 20, 'max_points_per_bin_tb',
                                              'Set the maximum number of points per beam angle bin for plotting')
        self.max_points_per_bin_tb.setValidator(QIntValidator(1, 100000))
        bin_decimation_layout = BoxLayout([max_points_per_bin_lbl, self.max_points_per_bin_tb], 'h')
        self.bin_decimation_gb = GroupBox('Bin-based decimation (uniform coverage)', bin_decimation_layout, True, True, 'bin_decimation_gb')
        self.bin_decimation_gb.setToolTip('Decimate points by beam angle bins to ensure uniform coverage across all beam angles. '
                                          'This prevents overplotting in dense beam angle regions while preserving sparse regions. '
                                          'Uses the same beam angle bins as the accuracy calculation.')
        self.pt_count_gb.setToolTip('To maintain reasonable plot and refresh speeds, the display will be limited '
                                    'by default to a total of ' + str(self.n_points_max_default) + ' soundings.  '
                                    'The limit is applied to new and archive datasets separately.  If needed, the user '
                                    'may specify a custom maximum point count.\n\n'
                                    'Reduction of each dataset is accomplished by simple decimation as a final step '
                                    'after all user-defined filtering (depth, angle, backscatter, etc.).  Non-integer '
                                    'decimation factors are handled using nearest-neighbor interpolation; soundings '
                                    'are not altered, just downsampled to display the maximum count allowed by the '
                                    'user parameters.'
                                    '\n\nAlternatively, the user may also specify a custom decimation factor.  '
                                    'Each dataset will be downsampled according to the more aggressive of the two '
                                    'inputs (max. count or dec. fac.) to achieve the greatest reduction in total '
                                    'displayed sounding count.  Unchecking these options will revert to the default.  '
                                    'In any case, large sounding counts may significantly slow the plotting process.')

        # set up layout and groupbox for tabs
        tab2_xline_filter_layout = BoxLayout([self.angle_xline_gb, self.depth_xline_gb, self.dz_abs_gb, self.dz_pct_gb,
                                              self.bs_xline_gb, self.depth_mode_gb, self.bin_count_gb], 'v')
        tab2_xline_filter_gb = GroupBox('Crosslines', tab2_xline_filter_layout,
                                        False, True, 'tab2_xline_filter_gb')
        tab2_xline_filter_gb.setEnabled(True)

        # add filter management buttons
        self.save_filters_btn = PushButton('Store Filters', 95, 25, 'save_filters_btn',
                                           'Save current Reference Surface and Crossline filter settings to a parameter file')
        self.load_last_filters_btn = PushButton('Load Filters', 95, 25, 'load_last_filters_btn',
                                                'Load the last saved set of filter settings')
        self.default_filters_btn = PushButton('Restore Filters', 95, 25, 'default_filters_btn',
                                              'Reset all filters to their default values')
        
        filter_buttons_layout = BoxLayout([self.save_filters_btn, self.load_last_filters_btn, self.default_filters_btn], 'h')
        filter_buttons_gb = GroupBox('Filter Management', filter_buttons_layout, False, False, 'filter_buttons_gb')

        # set up tabs
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setStyleSheet("background-color: none")

        # add plot settings management buttons
        self.save_plot_settings_btn = PushButton('Store Settings', 95, 25, 'save_plot_settings_btn',
                                                 'Save current Plot tab settings to a parameter file')
        self.load_last_plot_settings_btn = PushButton('Load Settings', 95, 25, 'load_last_plot_settings_btn',
                                                      'Load the last saved set of plot settings')
        self.default_plot_settings_btn = PushButton('Restore Defaults', 95, 25, 'default_plot_settings_btn',
                                                    'Reset all plot settings to their default values')
        
        plot_settings_buttons_layout = BoxLayout([self.save_plot_settings_btn, self.load_last_plot_settings_btn, self.default_plot_settings_btn], 'h')
        plot_settings_buttons_gb = GroupBox('Plot Settings Management', plot_settings_buttons_layout, False, False, 'plot_settings_buttons_gb')

        # set up tab 1: plot options
        self.tab1 = QtWidgets.QWidget()
        self.tab1.layout = BoxLayout([self.custom_info_gb, self.data_ref_gb, pt_param_gb,
                                      self.plot_lim_gb, toggle_gb, self.flatten_mean_gb, plot_settings_buttons_gb], 'v')
        self.tab1.layout.addStretch()
        self.tab1.setLayout(self.tab1.layout)

        # set up tab 2: filtering options
        self.tab2 = QtWidgets.QWidget()
        label_in_prog = Label('FILTERS IN PROGRESS\nNOT APPLIED TO DATA')
        label_in_prog.setStyleSheet("color: red")
        # self.tab2.layout = BoxLayout([label_in_prog, self.angle_gb, self.depth_ref_gb, self.depth_xline_gb,
        #                               self.bs_gb, self.slope_gb, self.density_gb], 'v')
        self.tab2.layout = BoxLayout([tab2_ref_filter_gb, tab2_xline_filter_gb, self.pt_count_gb, self.bin_decimation_gb, filter_buttons_gb], 'v')

        self.tab2.layout.addStretch()
        self.tab2.setLayout(self.tab2.layout)

        # add tabs to tab layout
        self.tabs.addTab(self.tab1, 'Plot')
        self.tabs.addTab(self.tab2, 'Filter')

        self.tabw = 340  # set fixed tab width
        self.tabs.setFixedWidth(self.tabw)

        self.right_layout = BoxLayout([self.tabs], 'v')
        self.right_layout.addStretch()

    def _detect_ref_filter_changes(self, old_filters):
        """Detect if reference surface filters have changed"""
        try:
            # Check if any reference surface filter values or states have changed
            current_filters = {
                'depth_enabled': self.depth_ref_gb.isChecked(),
                'depth_min': self.min_depth_ref_tb.text(),
                'depth_max': self.max_depth_ref_tb.text(),
                'slope_enabled': self.slope_gb.isChecked(),
                'slope_max': self.max_slope_tb.text(),
                'slope_window': self.slope_win_cbox.currentText(),
                'density_enabled': self.density_gb.isChecked(),
                'density_min': self.min_dens_tb.text(),
                'uncertainty_enabled': self.uncertainty_gb.isChecked(),
                'uncertainty_max': self.max_u_tb.text()
            }
            
            print(f"DEBUG: Old ref filters: {old_filters}")
            print(f"DEBUG: Current ref filters: {current_filters}")
            
            if (old_filters['depth_enabled'] != current_filters['depth_enabled'] or
                old_filters['depth_min'] != current_filters['depth_min'] or
                old_filters['depth_max'] != current_filters['depth_max'] or
                old_filters['slope_enabled'] != current_filters['slope_enabled'] or
                old_filters['slope_max'] != current_filters['slope_max'] or
                old_filters['slope_window'] != current_filters['slope_window'] or
                old_filters['density_enabled'] != current_filters['density_enabled'] or
                old_filters['density_min'] != current_filters['density_min'] or
                old_filters['uncertainty_enabled'] != current_filters['uncertainty_enabled'] or
                old_filters['uncertainty_max'] != current_filters['uncertainty_max']):
                print("DEBUG: Reference filters changed!")
                return True
            print("DEBUG: No reference filter changes detected")
            return False
        except Exception as e:
            print(f"Error detecting reference filter changes: {str(e)}")
            return True  # Assume changes occurred if detection fails

    def _detect_xline_filter_changes(self, old_filters):
        """Detect if crossline filters have changed"""
        try:
            # Check if any crossline filter values or states have changed
            if (old_filters['depth_enabled'] != self.depth_xline_gb.isChecked() or
                old_filters['depth_min'] != self.min_depth_xline_tb.text() or
                old_filters['depth_max'] != self.max_depth_xline_tb.text() or
                old_filters['angle_enabled'] != self.angle_xline_gb.isChecked() or
                old_filters['angle_min'] != self.min_angle_xline_tb.text() or
                old_filters['angle_max'] != self.max_angle_xline_tb.text() or
                old_filters['bs_enabled'] != self.bs_xline_gb.isChecked() or
                old_filters['bs_min'] != self.min_bs_xline_tb.text() or
                old_filters['bs_max'] != self.max_bs_xline_tb.text() or
                old_filters['dz_abs_enabled'] != (hasattr(self, 'dz_abs_gb') and self.dz_abs_gb.isChecked()) or
                old_filters['dz_max'] != self.max_dz_tb.text() or
                old_filters['dz_pct_enabled'] != (hasattr(self, 'dz_pct_gb') and self.dz_pct_gb.isChecked()) or
                old_filters['dz_wd_max'] != self.max_dz_wd_tb.text() or
                old_filters['depth_mode_enabled'] != self.depth_mode_gb.isChecked() or
                old_filters['depth_mode'] != self.depth_mode_cbox.currentText() or
                old_filters['bin_count_enabled'] != self.bin_count_gb.isChecked() or
                old_filters['bin_count_min'] != self.min_bin_count_tb.text()):
                return True
            return False
        except Exception as e:
            print(f"Error detecting crossline filter changes: {str(e)}")
            return True  # Assume changes occurred if detection fails
        
    def set_main_layout(self):
        # set the main layout with file controls on left and swath figure on right
        main_layout = QtWidgets.QHBoxLayout()
        main_layout.addLayout(self.left_layout)
        # main_layout.addLayout(self.swath_layout)
        main_layout.addLayout(self.center_layout)
        main_layout.addLayout(self.right_layout)
        
        self.mainWidget.setLayout(main_layout)

    def _setup_density_tooltip(self):
        """Setup mouse motion event handler for density plot tooltip"""
        try:
            # Create a text annotation for the tooltip
            self.density_tooltip = self.density_final_ax.text(0.02, 0.98, '', 
                                                              transform=self.density_final_ax.transAxes,
                                                              bbox=dict(boxstyle='round,pad=0.3', 
                                                                       facecolor='white', 
                                                                       alpha=0.8,
                                                                       edgecolor='black'),
                                                              fontsize=8,
                                                              verticalalignment='top',
                                                              horizontalalignment='left')
            self.density_tooltip.set_visible(False)
            
            # Connect mouse motion event
            self.density_final_canvas.mpl_connect('motion_notify_event', self._on_density_mouse_move)
            
        except Exception as e:
            print(f"Error setting up density tooltip: {str(e)}")

    def _on_density_mouse_move(self, event):
        """Handle mouse motion on density plot to show tooltip"""
        try:
            if event.inaxes != self.density_final_ax:
                self.density_tooltip.set_visible(False)
                self.density_final_canvas.draw_idle()
                return
            
            # Check if density data is available
            if 'c_grid' not in self.ref or 'c_mask' not in self.ref:
                return
            
            # Get mouse position in data coordinates
            x, y = event.xdata, event.ydata
            
            if x is None or y is None:
                self.density_tooltip.set_visible(False)
                self.density_final_canvas.draw_idle()
                return
            
            # Convert data coordinates to grid indices
            extent = self.ref['z_ref_extent']
            grid_shape = self.ref['c_grid'].shape
            
            # Calculate grid indices
            # Note: imshow flips the y-axis, so we need to account for this
            x_idx = int((x - extent[0]) / (extent[1] - extent[0]) * grid_shape[1])
            y_idx = int((extent[3] - y) / (extent[3] - extent[2]) * grid_shape[0])
            
            # Check bounds
            if 0 <= x_idx < grid_shape[1] and 0 <= y_idx < grid_shape[0]:
                # Get density value
                density_value = self.ref['c_grid'][y_idx, x_idx]
                mask_value = self.ref['c_mask'][y_idx, x_idx]
                
                # Only show tooltip if the point is not masked (has a valid density value)
                if not np.isnan(density_value) and not np.isnan(mask_value):
                    # Format the tooltip text
                    tooltip_text = f'Density: {density_value:.1f} soundings/cell\nEasting: {x:.0f} m\nNorthing: {y:.0f} m'
                    self.density_tooltip.set_text(tooltip_text)
                    self.density_tooltip.set_visible(True)
                else:
                    self.density_tooltip.set_visible(False)
            else:
                self.density_tooltip.set_visible(False)
            
            # Redraw the canvas
            self.density_final_canvas.draw_idle()
            
        except Exception as e:
            print(f"Error in density mouse move handler: {str(e)}")
            # Hide tooltip on error
            if hasattr(self, 'density_tooltip'):
                self.density_tooltip.set_visible(False)
                self.density_final_canvas.draw_idle()

    def _setup_depth_tooltip(self):
        """Setup mouse motion event handler for depth plot tooltip"""
        try:
            # Create a text annotation for the tooltip
            self.depth_tooltip = self.depth_ax.text(0.02, 0.98, '', 
                                                   transform=self.depth_ax.transAxes,
                                                   bbox=dict(boxstyle='round,pad=0.3', 
                                                            facecolor='white', 
                                                            alpha=0.8,
                                                            edgecolor='black'),
                                                   fontsize=8,
                                                   verticalalignment='top',
                                                   horizontalalignment='left')
            self.depth_tooltip.set_visible(False)
            
            # Connect mouse motion event
            self.depth_canvas.mpl_connect('motion_notify_event', self._on_depth_mouse_move)
            
        except Exception as e:
            print(f"Error setting up depth tooltip: {str(e)}")

    def _on_depth_mouse_move(self, event):
        """Handle mouse motion on depth plot to show tooltip"""
        try:
            if event.inaxes != self.depth_ax:
                self.depth_tooltip.set_visible(False)
                self.depth_canvas.draw_idle()
                return
            
            # Check if depth data is available
            if 'z_grid' not in self.ref or 'z_mask' not in self.ref:
                return
            
            # Get mouse position in data coordinates
            x, y = event.xdata, event.ydata
            
            if x is None or y is None:
                self.depth_tooltip.set_visible(False)
                self.depth_canvas.draw_idle()
                return
            
            # Convert data coordinates to grid indices
            extent = self.ref['z_ref_extent']
            grid_shape = self.ref['z_grid'].shape
            
            # Calculate grid indices
            # Note: imshow flips the y-axis, so we need to account for this
            x_idx = int((x - extent[0]) / (extent[1] - extent[0]) * grid_shape[1])
            y_idx = int((extent[3] - y) / (extent[3] - extent[2]) * grid_shape[0])
            
            # Check bounds
            if 0 <= x_idx < grid_shape[1] and 0 <= y_idx < grid_shape[0]:
                # Get depth value
                depth_value = self.ref['z_grid'][y_idx, x_idx]
                mask_value = self.ref['z_mask'][y_idx, x_idx]
                
                # Only show tooltip if the point is not masked (has a valid depth value)
                if not np.isnan(depth_value) and not np.isnan(mask_value):
                    # Format the tooltip text
                    tooltip_text = f'Depth: {depth_value:.1f} m\nEasting: {x:.0f} m\nNorthing: {y:.0f} m'
                    self.depth_tooltip.set_text(tooltip_text)
                    self.depth_tooltip.set_visible(True)
                else:
                    self.depth_tooltip.set_visible(False)
            else:
                self.depth_tooltip.set_visible(False)
            
            # Redraw the canvas
            self.depth_canvas.draw_idle()
            
        except Exception as e:
            print(f"Error in depth mouse move handler: {str(e)}")
            # Hide tooltip on error
            if hasattr(self, 'depth_tooltip'):
                self.depth_tooltip.set_visible(False)
                self.depth_canvas.draw_idle()

    def _setup_uncertainty_tooltip(self):
        """Setup mouse motion event handler for uncertainty plot tooltip"""
        try:
            # Create a text annotation for the tooltip
            self.uncertainty_tooltip = self.uncertainty_final_ax.text(0.02, 0.98, '', 
                                                                     transform=self.uncertainty_final_ax.transAxes,
                                                                     bbox=dict(boxstyle='round,pad=0.3', 
                                                                              facecolor='white', 
                                                                              alpha=0.8,
                                                                              edgecolor='black'),
                                                                     fontsize=8,
                                                                     verticalalignment='top',
                                                                     horizontalalignment='left')
            self.uncertainty_tooltip.set_visible(False)
            
            # Connect mouse motion event
            self.uncertainty_final_canvas.mpl_connect('motion_notify_event', self._on_uncertainty_mouse_move)
            
        except Exception as e:
            print(f"Error setting up uncertainty tooltip: {str(e)}")

    def _on_uncertainty_mouse_move(self, event):
        """Handle mouse motion on uncertainty plot to show tooltip"""
        try:
            if event.inaxes != self.uncertainty_final_ax:
                self.uncertainty_tooltip.set_visible(False)
                self.uncertainty_final_canvas.draw_idle()
                return
            
            # Check if uncertainty data is available
            if 'u_grid' not in self.ref or 'u_mask' not in self.ref:
                return
            
            # Get mouse position in data coordinates
            x, y = event.xdata, event.ydata
            
            if x is None or y is None:
                self.uncertainty_tooltip.set_visible(False)
                self.uncertainty_final_canvas.draw_idle()
                return
            
            # Convert data coordinates to grid indices
            extent = self.ref['z_ref_extent']
            grid_shape = self.ref['u_grid'].shape
            
            # Calculate grid indices
            # Note: imshow flips the y-axis, so we need to account for this
            x_idx = int((x - extent[0]) / (extent[1] - extent[0]) * grid_shape[1])
            y_idx = int((extent[3] - y) / (extent[3] - extent[2]) * grid_shape[0])
            
            # Check bounds
            if 0 <= x_idx < grid_shape[1] and 0 <= y_idx < grid_shape[0]:
                # Get uncertainty value
                uncertainty_value = self.ref['u_grid'][y_idx, x_idx]
                mask_value = self.ref['u_mask'][y_idx, x_idx]
                
                # Only show tooltip if the point is not masked (has a valid uncertainty value)
                if not np.isnan(uncertainty_value) and not np.isnan(mask_value):
                    # Format the tooltip text
                    tooltip_text = f'Uncertainty: {uncertainty_value:.2f} m\nEasting: {x:.0f} m\nNorthing: {y:.0f} m'
                    self.uncertainty_tooltip.set_text(tooltip_text)
                    self.uncertainty_tooltip.set_visible(True)
                else:
                    self.uncertainty_tooltip.set_visible(False)
            else:
                self.uncertainty_tooltip.set_visible(False)
            
            # Redraw the canvas
            self.uncertainty_final_canvas.draw_idle()
            
        except Exception as e:
            print(f"Error in uncertainty mouse move handler: {str(e)}")
            # Hide tooltip on error
            if hasattr(self, 'uncertainty_tooltip'):
                self.uncertainty_tooltip.set_visible(False)
                self.uncertainty_final_canvas.draw_idle()

    def _setup_slope_tooltip(self):
        """Setup mouse motion event handler for slope plot tooltip"""
        try:
            # Create a text annotation for the tooltip
            self.slope_tooltip = self.slope_final_ax.text(0.02, 0.98, '', 
                                                         transform=self.slope_final_ax.transAxes,
                                                         bbox=dict(boxstyle='round,pad=0.3', 
                                                                  facecolor='white', 
                                                                  alpha=0.8,
                                                                  edgecolor='black'),
                                                         fontsize=8,
                                                         verticalalignment='top',
                                                         horizontalalignment='left')
            self.slope_tooltip.set_visible(False)
            
            # Connect mouse motion event
            self.slope_final_canvas.mpl_connect('motion_notify_event', self._on_slope_mouse_move)
            
        except Exception as e:
            print(f"Error setting up slope tooltip: {str(e)}")

    def _on_slope_mouse_move(self, event):
        """Handle mouse motion on slope plot to show tooltip"""
        try:
            if event.inaxes != self.slope_final_ax:
                self.slope_tooltip.set_visible(False)
                self.slope_final_canvas.draw_idle()
                return
            
            # Check if slope data is available
            if 's_grid' not in self.ref or 's_mask' not in self.ref:
                return
            
            # Get mouse position in data coordinates
            x, y = event.xdata, event.ydata
            
            if x is None or y is None:
                self.slope_tooltip.set_visible(False)
                self.slope_final_canvas.draw_idle()
                return
            
            # Convert data coordinates to grid indices
            extent = self.ref['z_ref_extent']
            grid_shape = self.ref['s_grid'].shape
            
            # Calculate grid indices
            # Note: imshow flips the y-axis, so we need to account for this
            x_idx = int((x - extent[0]) / (extent[1] - extent[0]) * grid_shape[1])
            y_idx = int((extent[3] - y) / (extent[3] - extent[2]) * grid_shape[0])
            
            # Check bounds
            if 0 <= x_idx < grid_shape[1] and 0 <= y_idx < grid_shape[0]:
                # Get slope value
                slope_value = self.ref['s_grid'][y_idx, x_idx]
                mask_value = self.ref['s_mask'][y_idx, x_idx]
                
                # Only show tooltip if the point is not masked (has a valid slope value)
                if not np.isnan(slope_value) and not np.isnan(mask_value):
                    # Format the tooltip text
                    tooltip_text = f'Slope: {slope_value:.2f}°\nEasting: {x:.0f} m\nNorthing: {y:.0f} m'
                    self.slope_tooltip.set_text(tooltip_text)
                    self.slope_tooltip.set_visible(True)
                else:
                    self.slope_tooltip.set_visible(False)
            else:
                self.slope_tooltip.set_visible(False)
            
            # Redraw the canvas
            self.slope_final_canvas.draw_idle()
            
        except Exception as e:
            print(f"Error in slope mouse move handler: {str(e)}")
            # Hide tooltip on error
            if hasattr(self, 'slope_tooltip'):
                self.slope_tooltip.set_visible(False)
                self.slope_final_canvas.draw_idle()

    def handle_save_all_plots(self):
        """Handle save all plots with validation for ship name and description"""
        from PyQt6.QtWidgets import QMessageBox
        
        # Check if ship name and description are filled out
        ship_name = self.ship_tb.text().strip()
        description = self.cruise_tb.text().strip()
        
        missing_fields = []
        if not ship_name:
            missing_fields.append("Ship Name")
        if not description:
            missing_fields.append("Description")
        
        # If any fields are missing, show warning dialog
        if missing_fields:
            missing_text = " and ".join(missing_fields)
            msg = f"The following information is missing: {missing_text}\n\nDo you really want to continue with exporting the plots?"
            reply = QMessageBox.question(self, 'Missing Information', msg, 
                                        QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                                        QMessageBox.StandardButton.No)
            
            if reply == QMessageBox.StandardButton.No:
                return  # User chose not to continue
        
        # If we get here, either all fields are filled or user chose to continue
        # Call the original save_all_plots function
        from multibeam_tools.libs.swath_accuracy_lib import save_all_plots
        save_all_plots(self)


class NewPopup(QtWidgets.QWidget): # new class for additional plots
    def __init__(self):
        QtWidgets.QWidget.__init__(self)
        
        
if __name__ == '__main__':
    print("Starting Swath Accuracy Plotter...")
    print("Creating QApplication...")
    app = QtWidgets.QApplication(sys.argv)
    print("Creating MainWindow...")
    main = MainWindow()
    main.resize(1700, 1050)  # Set initial size to 1700x1050 pixels
    main.setFixedSize(1700, 1050)  # Prevent window resizing
    print("Showing MainWindow...")
    main.show()
    print("Entering event loop...")
    sys.exit(app.exec())
