"""Functions for swath accuracy plotting in NOAA / MAC echosounder assessment tools"""
import numpy as np
import utm
import os
import json
from copy import deepcopy

from PyQt6 import QtWidgets, QtGui
from PyQt6.QtGui import QDoubleValidator
from PyQt6.QtCore import Qt, QSize, QTimer

from libs.file_fun import *
from libs.swath_fun import *
from libs.readEM import convertXYZ, sort_active_pos_system

import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
from matplotlib import colors
from matplotlib import colorbar
from matplotlib.colors import LightSource
# from matplotlib import patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from scipy.interpolate import griddata
from time import process_time
import pyproj
import re
from scipy.spatial import cKDTree as KDTree
from scipy.ndimage import uniform_filter
from scipy.interpolate import interp1d
from datetime import timedelta, datetime
import matplotlib.dates as mdates


def setup(self):
	# initialize other necessities
	self.print_updates = False
	# self.print_updates = False
	self.xline = {}
	self.ref = {}
	self.xline_track = {}
	self.tide = {}
	self.ref_utm_str = 'N/A'
	# self.det = {}  # detection dict (new data)
	# self.det_archive = {}  # detection dict (archive data)
	# self.spec = {}  # dict of theoretical coverage specs
	# self.filenames = ['']  # initial file list
	# self.input_dir = ''  # initial input dir
	
	# Load last used directory from config file
	self.output_dir = load_last_directory()
	if not self.output_dir or not os.path.exists(self.output_dir):
		self.output_dir = os.getcwd()  # fallback to current working directory
	
	# self.clim_last_user = {'depth': [0, 1000], 'backscatter': [-50, -20]}
	# self.last_cmode = 'depth'
	# self.refresh_ref_plot = False
	# self.refresh_acc_plot = False
	self.clim_z = [0, 1]
	self.clim_c = [0, 1]
	self.clim_s = [0, 1]
	self.clim_u = [0, 1]

	self.cbar_ax1 = None  # initial colorbar for ref surf depth subplot
	self.cbar_ax2 = None  # initial colorbar for ref surf slope subplot
	self.cbar_ax3 = None  # initial colorbar for ref surf density subplot
	self.cbar_ax4 = None  # initial colorbar for ref surf uncertainty subplot
	self.cbar_ax5 = None  # initial colorbar for ref surf filtered depth subplot
	self.cbar_ax6 = None  # initial colorbar for ref surf final large plot
	self.tide_cbar_ax = None  # initial colorbar for tide plot
	self.legendbase = None  # initial legend
	self.cbar_font_size = 8  # colorbar/legend label size
	self.cbar_title_font_size = 8  # colorbar/legend title size
	self.cbar_loc = 1  # set upper right as default colorbar/legend location
	self.n_points_max_default = 50000  # default maximum number of points to plot in order to keep reasonable speed
	self.n_points_max = 50000
	self.n_points_plotted = 0
	# # self.n_points_plotted_arc = 0
	self.dec_fac_default = 1  # default decimation factor for point count
	self.dec_fac = 1
	# self.rtp_angle_buffer_default = 0  # default runtime angle buffer
	# self.rtp_angle_buffer = 0  # +/- deg from runtime parameter swath angle limit to filter RX angles
	# self.x_max = 0.0
	# self.z_max = 0.0
	self.model_list = ['EM 2040', 'EM 2042', 'EM 302', 'EM 304', 'EM 710', 'EM 712', 'EM 122', 'EM 124']
	# self.cmode_list = ['Depth', 'Backscatter', 'Ping Mode', 'Pulse Form', 'Swath Mode', 'Solid Color']
	# self.top_data_list = []
	# self.clim_list = ['All data', 'Filtered data', 'Fixed limits']
	self.data_ref_list = ['Waterline']  # , 'Origin', 'TX Array', 'Raw Data']
	self.tide_unit_dict = {'Meter': float(1), 'Foot': float(0.3048), 'Smoot': float(1/1.702)}  # conversion fac to m
	self.unit_mode = '%WD'  # default plot as % Water Depth; option to toggle alternative meters
	self.tide_applied = False
	# plt.margins(x=0.1, y=0.1)
	self.model_updated = False
	self.ship_name_updated = False
	self.cruise_name_updated = False
	self.sn_updated = False
	self.N_plotted = 0
	
	# Initialize flashing timer for Calc Accuracy button
	# self.flash_timer = QTimer()
	# self.flash_timer.timeout.connect(self.flash_calc_accuracy_button)
	# self.is_flashing = False


def init_all_axes(self):
	init_swath_ax(self)
	init_surf_ax(self)
	init_tide_ax(self)
	add_grid_lines(self)
	update_axes(self)


def init_swath_ax(self):  # set initial swath parameters
	self.pt_size = np.square(float(self.pt_size_cbox.currentText()))  # swath plot point size
	self.pt_size_cov = np.square(float(self.pt_size_cov_cbox.currentText()))  # coverage plot point size
	self.pt_alpha = np.divide(float(self.pt_alpha_acc_tb.text()), 100)  # accuracy plot opacity
	self.pt_alpha_cov = np.divide(float(self.pt_alpha_cov_tb.text()), 100)  # coverage plot opacity

	self.ax1 = self.swath_figure.add_subplot(211)
	self.ax2 = self.swath_figure.add_subplot(212)

	self.x_max_default = 75
	self.x_spacing_default = 15
	self.y_max_std_default = 0.5  # max y range of depth st. dev. plot (top subplot)
	self.y_max_bias_default = 1  # max +/- y range of depth bias (raw, mean, +/- 1 sigma, bottom subplot)
	self.axis_margin_default = 5  # % axis margin

	self.x_max_custom = self.x_max_default  # store future custom entries
	self.x_spacing_custom = self.x_spacing_default
	self.y_max_bias_custom = self.y_max_bias_default
	self.y_max_std_custom = self.y_max_std_default
	self.axis_margin_custom = self.axis_margin_default

	self.max_beam_angle_tb.setText(str(self.x_max_default))
	self.angle_spacing_tb.setText(str(self.x_spacing_default))
	self.max_bias_tb.setText(str(self.y_max_bias_default))
	self.max_std_tb.setText(str(self.y_max_std_default))
	self.axis_margin_tb.setText(str(self.axis_margin_default))

	self.cruise_name = ''
	self.swath_ax_margin = 1.1  # scale axes to multiple of max data in each direction
	self.fsize_title = 12
	self.fsize_label = 10
	self.lwidth = 1  # line width
	self.color = QtGui.QColor(0, 0, 0)  # set default solid color to black for new data
	self.archive_color = QtGui.QColor('darkGray')


def init_surf_ax(self):  # set initial ref surf parameters
	self.surf_ax1 = self.surf_figure.add_subplot(221)
	self.surf_ax2 = self.surf_figure.add_subplot(222, sharex=self.surf_ax1, sharey=self.surf_ax1)
	self.surf_ax3 = self.surf_figure.add_subplot(223, sharex=self.surf_ax1, sharey=self.surf_ax1)
	self.surf_ax4 = self.surf_figure.add_subplot(224, sharex=self.surf_ax1, sharey=self.surf_ax1)
	self.surf_ax5 = self.surf_final_figure.add_subplot(111)
	self.surf_ax5.set_aspect(1)
	
	# Add axis for density final tab
	self.density_final_ax = self.density_final_figure.add_subplot(111)
	self.density_final_ax.set_aspect(1)
	
	# Add axis for slope final tab
	self.slope_final_ax = self.slope_final_figure.add_subplot(111)
	self.slope_final_ax.set_aspect(1)
	
	# Add axis for depth tab
	self.depth_ax = self.depth_figure.add_subplot(111)
	self.depth_ax.set_aspect(1)
	
	# Add axis for uncertainty final tab
	self.uncertainty_final_ax = self.uncertainty_final_figure.add_subplot(111)
	self.uncertainty_final_ax.set_aspect(1)

	# setup dict of colorbar parameters; u and z_final are separate colorbars (u added later) plotted on same axis
	self.cbar_dict = {'z': {'cax': self.cbar_ax1, 'ax': self.surf_ax1, 'clim': self.clim_z, 'label': 'Depth (m)'},
			  'c': {'cax': self.cbar_ax2, 'ax': self.surf_ax2, 'clim': self.clim_c, 'label': 'Soundings/Cell'},
			  's': {'cax': self.cbar_ax3, 'ax': self.surf_ax3, 'clim': self.clim_s, 'label': 'Slope (deg)'},
			  'u': {'cax': self.cbar_ax4, 'ax': self.surf_ax4, 'clim': self.clim_u, 'label': 'Uncertainty (m)'},
			  'z_filt': {'cax': self.cbar_ax5, 'ax': self.surf_ax4, 'clim': self.clim_z, 'label': 'Depth (m)'},
			  'z_final': {'cax': self.cbar_ax6, 'ax': self.surf_ax5, 'clim': self.clim_z, 'label': 'Depth (m)'},
			  'z_depth': {'cax': None, 'ax': self.depth_ax, 'clim': self.clim_z, 'label': 'Depth (m)'},
			  'c_final': {'cax': None, 'ax': self.density_final_ax, 'clim': self.clim_c, 'label': 'Soundings/Cell'},
			  's_final': {'cax': None, 'ax': self.slope_final_ax, 'clim': self.clim_s, 'label': 'Slope (deg)'},
			  'u_final': {'cax': None, 'ax': self.uncertainty_final_ax, 'clim': self.clim_u, 'label': 'Uncertainty (m)'}}


def init_tide_ax(self):  # set initial tide plot parameters
	self.tide_ax = self.tide_figure.add_subplot(111)
	self.tide_cbar_dict = {'tide_amp': {'cax': self.tide_cbar_ax, 'ax': self.tide_ax, 'label': 'Amplitude (m)'}}





def update_buttons(self, recalc_acc=False):
	# enable or disable file selection and calc_accuracy buttons depending on loaded files
	print('updating buttons...')
	get_current_file_list(self)
	fnames_ref = [f for f in self.filenames if '.xyz' in f]
	# fnames_xline = get_new_file_list(self, ['.all', '.kmall'], [])  # list new .all files not in det dict
	fnames_xline = get_new_file_list(self, ['.all', '.kmall', 'ASCII.txt'], [])  # list new .all files not in det dict
	fnames_xline_all = [f for f in self.filenames if any(ext in f for ext in ['.all', '.kmall', 'ASCII.txt'])]  # all crossline files

	fnames_tide = [f for f in self.filenames if '.tid' in f]

	self.add_ref_surf_btn.setEnabled(len(fnames_ref) == 0)  # enable ref surf selection only if none loaded
	self.add_dens_surf_btn.setEnabled(('z' in self.ref.keys() and 'c' not in self.ref.keys()))  # enable if ref z avail
	# Only enable Add Tide if at least one crossline and no tide file loaded
	self.add_tide_btn.setEnabled(len(fnames_xline_all) > 0 and len(fnames_tide) == 0)
	
	# Style Add Tide button when crossline files are loaded
	# Check if accuracy has been calculated (indicates tide has been processed)
	accuracy_calculated = (hasattr(self, 'xline') and 
	                      len(self.xline) > 0 and 
	                      any([k in self.xline.keys() for k in ['dz_ref_wd', 'dz_ref']]))
	
	if len(fnames_xline_all) > 0 and len(fnames_tide) == 0 and not accuracy_calculated:
		self.add_tide_btn.setStyleSheet("background-color: lightyellow; color: black; font-weight: bold")
	else:
		self.add_tide_btn.setStyleSheet("color: black; font-weight: normal")

	# enable calc_accuracy button only if one ref surf and at least one crossline are loaded
	if (len(fnames_ref) == 1 and len(fnames_xline) > 0):
		self.calc_accuracy_btn.setEnabled(True)

		if recalc_acc:
			self.calc_accuracy_btn.setStyleSheet("background-color: yellow; color: black; font-weight: bold")
			# Start flashing the button text
			# if not self.is_flashing:
			# 	self.flash_timer.start(500)  # Flash every 500ms
			# 	self.is_flashing = True

	else:
		self.calc_accuracy_btn.setEnabled(False)
		self.calc_accuracy_btn.setStyleSheet("color: black; font-weight: normal")
		# Stop flashing when button is disabled
		# if self.is_flashing:
		# 	self.stop_flashing_calc_accuracy_button()


def add_ref_file(self, ftype_filter, input_dir='HOME', include_subdir=False):
	# add single reference surface file with extensions in ftype_filter
	fname = add_files(self, ftype_filter, input_dir, include_subdir, multiselect=False)
	
	update_file_list(self, fname)
	# try to get UTM zone from filename; zone can be, e.g,, 'UTM-11S', '14N',  w/ or w/o UTM preceding and -, _, or ' '
	# get decimal and hemisphere, strip zero padding and remove spaces for comparison to UTM combobox list
	if fname and len(fname) > 0:
		fname_str = fname[0]
		fname_str = fname_str[fname_str.rfind('/') + 1:].rstrip()

		try:
			utm_idx = -1
			utm_str = re.search(r"[_-]*\s*[0-9]{1,2}[NS]", fname_str).group()
			utm_str = re.search(r'\d+[NS]', utm_str).group().strip('0').replace(' ', '')
			utm_idx = self.ref_proj_cbox.findText(utm_str)
			print('found utm_str, utm_idx =', utm_str, utm_idx)

		except:
			update_log(self, 'Please select the reference surface UTM zone')
			self.ref_proj_cbox.setCurrentIndex(0)

		if utm_idx > -1:
			update_log(self, 'Found UTM zone from filename: ' + utm_str)
			self.ref_proj_cbox.setCurrentIndex(utm_idx)
			self.ref_utm_str = utm_str
			parse_ref_depth(self)

			if any([f for f in self.filenames if '.xyd' in f]) and 'c' not in self.ref:
				print('processing density file already loaded but not added in self.ref yet (')
				process_dens_file(self)

			make_ref_surf(self)
			refresh_plot(self, refresh_list=['ref'], set_active_tab=1, sender='add_ref_file')

	update_buttons(self, recalc_acc=True)


def update_ref_utm_zone(self):
	# update ref surf UTM zone after user selection, transform crossline data into current zone, and recalc accuracy
	self.ref_utm_str = self.ref_proj_cbox.currentText()  # update with current UTM zone
	self.ref['utm_zone'] = self.ref_proj_cbox.currentText()

	if self.xline:  # crossline data already processed, need to recalc accuracy using updated ref surf
		update_log(self, 'Reference surface UTM zone has been updated; recalculating accuracy with loaded crosslines')
		convert_track_utm(self)  # update track
		calc_accuracy(self, recalc_utm_only=True)  # skip parsing, convert_crossline_utm, recalc stats, refresh plots

	
	else:
		update_log(self, 'Reference surface has been updated; no crossline data available')
		refresh_plot(self, refresh_list=['ref'], sender='udate_ref_utm_zone (self.xline=FALSE)')


def add_dens_file(self, ftype_filter, input_dir='HOME', include_subdir=False):
	# add single density surface file with extensions in ftype_filter
	fname = add_files(self, ftype_filter, input_dir, include_subdir, multiselect=False)
	
	update_file_list(self, fname)
	process_dens_file(self)


def process_dens_file(self):
	# process the loaded density file
	# update_buttons(self)
	try:
		parse_ref_dens(self)
		make_ref_surf(self)
		update_buttons(self)  # turn off density button if loaded
		refresh_plot(self, refresh_list=['ref'], set_active_tab=1, sender='process_dens_file')
	except Exception as e:
		update_log(self, f'ERROR: Failed to process density file: {str(e)}', font_color="red")
		print(f'Process density file error: {e}')


def add_tide_file(self, ftype_filter, input_dir='HOME', include_subdir=False):
	fname = add_files(self, ftype_filter, input_dir, include_subdir, multiselect=False)
	
	update_file_list(self, fname)
	update_buttons(self, recalc_acc=True)  # turn off tide button if loaded
	process_tide(self)


def process_tide(self, unit_set_by_user=False):
	parse_tide(self, unit_set_by_user=unit_set_by_user)
	plot_tide(self, set_active_tab=True)

	refresh_plot(self, refresh_list=['tide'], set_active_tab=3, sender='add_tide_file')


def add_acc_files(self, ftype_filter, input_dir='HOME', include_subdir=False):
	# add accuracy crossline files with extensions in ftype_filter from input_dir and subdir if desired
	fnames = add_files(self, ftype_filter, input_dir, include_subdir)
	
	update_file_list(self, fnames)
	update_buttons(self, recalc_acc=True)


def remove_acc_files(self):  # remove selected files only
	# recalc_acc = False
	get_current_file_list(self)
	selected_files = self.file_list.selectedItems()
	# fnames_ref = [f for f in self.filenames if '.xyz' in f]
	# fnames_dens = [f for f in self.filenames if '.xyd' in f]
	# fnames_tide = [f for f in self.filenames if '.tid' in f]
	# fnames_xline = [f for f in self.filenames if f.rsplit('.')[-1] in ['.all', '.kmall']]

	# print('in remove_acc_files, fnames_xline is', fnames_xline)

	# if len(fnames_xline) + len(fnames_ref) == 0:  # all .all and .xyz files have been removed, reset det dicts
	# 	self.xline = {}
	# 	self.ref = {}
	# self.bin_beamwise()  # call bin_beamwise with empty xline results to clear plots
	#            self.xline_archive = {}

	if not selected_files:  # files exist but nothing is selected
		update_log(self, 'No files selected for removal.')
		return

	# elif len(fnames_xline) ==

	else:  # remove only the files that have been selected
		for f in selected_files:
			fname = f.text().split('/')[-1]
			print('working on fname', fname)
			self.file_list.takeItem(self.file_list.row(f))
			update_log(self, 'Removed ' + fname)

			try:  # try to remove detections associated with this file
				if fname.rsplit('.')[-1] in ['all', 'kmall']:
					# get indices of soundings in det dict with matching filenames
					i = [j for j in range(len(self.xline['fname'])) if self.xline['fname'][j] == fname]

					for k in self.xline.keys():  # loop through all keys and remove values at these indices
						print(k)
						self.xline[k] = np.delete(self.xline[k], i).tolist()

					# remove trackline associated with this file
					for k in self.xline_track.keys():
						print('in remove_acc_files, removing xline_track key = ', k)
						# if self.xline_track[k]['fname'] == fname:
						if k == fname:
							print('k == fname, trying to pop')
							# self.xline_track[k].pop()
							self.xline_track.pop(k, None)
							print('success')
						else:
							print('skipping popping this k =', fname)

					# self.xline_track[fname].pop()

				elif '.xyz' in fname:  # remove the single reference surface and reset ref dict

					self.ref = {}
					self.add_ref_surf_btn.setEnabled(True)  # enable button to add replacement reference surface
					self.ref_proj_cbox.setCurrentIndex(0)

				elif '.xyd' in fname:  # remove the single reference density and reset associated ref dict keys
					# self.ref = {}
					# remove density data from ref dict, if it exists
					print('in remove_acc_files, self.ref.keys = ', self.ref.keys())
					for k in ['fname_dens', 'c', 'c_grid']:
						print('in remove_acc_files, popping density key =', k)
						self.ref.pop(k, None)

					self.add_dens_surf_btn.setEnabled(True)  # enable button to add replacement density surface

				elif '.tid' in fname:  # remove the singe tide and reset tide dict
					self.tide = {}
					self.tide_ax.clear()
					self.add_tide_btn.setEnabled(True)
					self.tide_applied = False

			except:  # will fail if detection dict has not been created yet (e.g., if calc_coverage has not been run)
				#                    update_log(self, 'Failed to remove soundings stored from ' + fname)
				pass

		# call bin_beamwise() to update results if any crosslines files are removed (not required to recalc ref surf)
		f_exts = [f.text().split('/')[-1].rsplit('.', 1)[-1] for f in selected_files]
		print('got f_exts =', f_exts)

		# recalculate beamwise results if crossline soundings were removed
		if any([ext in ['all', 'kmall'] for ext in f_exts]) and self.xline:
			print('found all or kmall extension, calling bin_beamwise from remove_files')
			bin_beamwise(self)

		# recalculate accuracy if crossline data exists and the tide files change
		if any([ext in ['tid'] for ext in f_exts]) and self.xline:
			print('found tid extension, calling calc_accuracy from remove_files')
			calc_accuracy(self)

		# recalculate reference surface if depth or density
		if any([ext in ['xyz', 'xyd'] for ext in f_exts]):
			print('at end of remove_acc_files, found xyz or xyd, caling make_ref_surf')
			make_ref_surf(self)

	# reset xline dict if all files have been removed
	get_current_file_list(self)
	print('in remove_acc_files, self.filenames after file removal is', self.filenames)
	fnames_xline = [f for f in self.filenames if f.rsplit('.')[-1] in ['all', 'kmall']]
	print('in remove_acc_files, fnames_xline after file removal is', fnames_xline)
	if len(fnames_xline) == 0:  # all .all and .xyz files have been removed, reset det dicts
		self.xline = {}
		self.xline_track = {}

	update_buttons(self)
	refresh_plot(self, refresh_list=['acc', 'ref', 'tide'], sender='remove_acc_files')  # refresh with updated (reduced or cleared) detection data


def clear_files(self):
	self.file_list.clear()  # clear the file list display
	self.filenames = []  # clear the list of (paths + files) passed to calc_coverage
	self.xline = {}  # clear current non-archive detections
	self.ref = {}  # clear ref surf data
	self.tide = {}
	bin_beamwise(self)  # call bin_beamwise with empty self.xline to reset all other binned results
	remove_acc_files(self)  # remove files and refresh plot
	update_log(self, 'Cleared all files')
	self.current_file_lbl.setText('Current File [0/0]:')
	self.calc_pb.setValue(0)
	self.add_ref_surf_btn.setEnabled(True)
	self.ref_proj_cbox.setCurrentIndex(0)
	
	# Update button styling after clearing files
	update_buttons(self)
	
	# Stop flashing when files are cleared
	# if hasattr(self, 'is_flashing') and self.is_flashing:
	# 	self.stop_flashing_calc_accuracy_button()


def calc_accuracy(self, recalc_utm_only=False, recalc_bins_only=False, recalc_dz_only=False, recalc_ref_only=False):
	# calculate accuracy of soundings from at least one crossline over exactly one reference surface
	# calc_accuracy is called after all filter updates; skip calc attempt if no crossline files are loaded
	
	from time import process_time
	start_time = process_time()
	
	# Handle reference surface recalculation if needed (can be done without crossline files)
	if recalc_ref_only:
		update_log(self, '=== REFERENCE SURFACE RECALCULATION ===')
		update_log(self, f'Current ref filters - Depth: {self.min_depth_ref_tb.text()}-{self.max_depth_ref_tb.text()}, Slope: {self.max_slope_tb.text()}, Density: {self.min_dens_tb.text()}, Uncertainty: {self.max_u_tb.text()}')
		update_log(self, f'Current ref filter states - Depth: {self.depth_ref_gb.isChecked()}, Slope: {self.slope_gb.isChecked()}, Density: {self.density_gb.isChecked()}, Uncertainty: {self.uncertainty_gb.isChecked()}')
		
		ref_start = process_time()
		make_ref_surf(self)  # rebuild reference surface from raw data
		ref_time = process_time() - ref_start
		update_log(self, f'Reference surface generation completed in {ref_time:.2f} seconds')
		
		plot_start = process_time()
		plot_ref_surf(self)  # update reference surface plots
		plot_time = process_time() - plot_start
		update_log(self, f'Reference surface plotting completed in {plot_time:.2f} seconds')
		
		refresh_plot(self, refresh_list=['ref'], set_active_tab=2, sender='calc_accuracy')
		total_time = process_time() - start_time
		update_log(self, f'=== REFERENCE SURFACE RECALCULATION COMPLETED in {total_time:.2f} seconds ===')
		return
	
	# For other calculations, require crossline files
	if not get_new_file_list(self, ['.all', '.kmall', 'ASCII.txt'], []):
		return
	else:
		update_log(self, '=== STARTING ACCURACY CALCULATIONS ===')

	force_sys_info_update = not any([recalc_utm_only, recalc_bins_only, recalc_dz_only, recalc_ref_only]) # force sys info update unless skipping

	print('in calc_accuracy with kwargs ', recalc_utm_only, recalc_bins_only, recalc_dz_only, recalc_ref_only)

	if not recalc_bins_only:  # skip if just recalculating bins;
		if not recalc_dz_only:  # skip if just recalculating dz from ref after filtering
			if not recalc_utm_only:  # skip if simply updating UTM zone for existing data; parse crosslines and calc z_final
				if 'c_grid' not in self.ref:
					update_log(self, 'Parsing reference density data...')
					parse_ref_dens(self)  # parse density data if not available

				update_log(self, 'Parsing crossline files...')
				num_new_xlines = parse_crosslines(self)  # parse the crossline(s)
				update_log(self, f'Parsed {num_new_xlines} new crossline files')

				if num_new_xlines > 0 or not self.tide_applied:  # (re)calc z_final if new files or tide has not been applied
					update_log(self, 'Calculating final depths with tide and reference adjustments...')
					z_final_start = process_time()
					calc_z_final(self)  # adjust sounding depths to desired reference, adjust for tide, and flip sign as necessary
					z_final_time = process_time() - z_final_start
					update_log(self, f'Final depth calculations completed in {z_final_time:.2f} seconds')

			update_log(self, 'Converting crossline coordinates to reference UTM zone...')
			utm_start = process_time()
			convert_crossline_utm(self)  # convert crossline X,Y to UTM zone of reference surface
			utm_time = process_time() - utm_start
			update_log(self, f'UTM conversion completed in {utm_time:.2f} seconds')

		update_log(self, 'Interpolating reference surface onto crossline positions...')
		interp_start = process_time()
		calc_dz_from_ref_interp(self)  # interpolate ref surf onto sounding positions, take difference
		interp_time = process_time() - interp_start
		update_log(self, f'Reference interpolation completed in {interp_time:.2f} seconds')

	# after parsing, converting, and calculating dz as necessary, bin results (apply xline filters in binning step)
	update_log(self, 'Binning soundings by beam angle...')
	bin_start = process_time()
	bin_beamwise(self)  # bin the results by beam angle
	bin_time = process_time() - bin_start
	update_log(self, f'Beam angle binning completed in {bin_time:.2f} seconds')
	
	update_log(self, 'Updating system information...')
	sys_start = process_time()
	update_system_info(self, self.xline, force_update=force_sys_info_update, fname_str_replace='_trimmed')
	sys_time = process_time() - sys_start
	update_log(self, f'System info update completed in {sys_time:.2f} seconds')

	if not recalc_dz_only:  # update mode list if this is the first calc; otherwise, leave mode list unchanged
		update_depth_mode_list(self)

	update_log(self, 'Generating plots...')
	plot_start = process_time()
	refresh_plot(self, refresh_list=['acc', 'ref', 'tide'], set_active_tab=0, sender='calc_accuracy')
	plot_time = process_time() - plot_start
	update_log(self, f'Plot generation completed in {plot_time:.2f} seconds')
	
	total_time = process_time() - start_time
	update_log(self, f'=== ACCURACY CALCULATIONS COMPLETED in {total_time:.2f} seconds ===')
	
	# Update button styling after accuracy calculations are complete
	update_buttons(self)
	
	# Stop flashing the Calc Accuracy button after calculation is complete
	# if self.is_flashing:
	# 	self.stop_flashing_calc_accuracy_button()

def update_depth_mode_list(self):  # determine set of available depth modes and update combobox for filtering if desired
	# Check if xline data exists and has ping_mode key
	if not self.xline or 'ping_mode' not in self.xline:
		# No crossline data loaded yet, set default list
		self.depth_mode_list = ['None']
		self.depth_mode_cbox.clear()
		self.depth_mode_cbox.addItems(self.depth_mode_list)
		return
	
	# Extract unique depth modes from ping_mode data
	try:
		self.depth_mode_list = [mode.split('(')[0].strip() for mode in set([pm for pm in self.xline['ping_mode']])]
		self.depth_mode_cbox.clear()
		self.depth_mode_cbox.addItems(self.depth_mode_list)
	except (KeyError, AttributeError, TypeError) as e:
		print(f"Error updating depth mode list: {e}")
		# Fallback to default
		self.depth_mode_list = ['None']
		self.depth_mode_cbox.clear()
		self.depth_mode_cbox.addItems(self.depth_mode_list)

def filter_xline(self, print_updates=True):
	# Check if xline data exists
	if not self.xline or 'z' not in self.xline:
		if print_updates:
			print('No crossline data available for filtering')
		return
	
	# set up indices for optional masking on angle, depth, bs; all idx true until fail optional filter settings
	# all soundings masked for nans (e.g., occasional nans in EX0908 data)
	# print('in filter_xline with fields:', self.xline.keys())

	# set up indices for each filter parameter
	idx_shape = np.shape(np.asarray(self.xline['z']))
	angle_idx = np.ones(idx_shape)  # inx of beam angle
	depth_idx = np.ones(idx_shape)  # idx of depth final (meters) after tide and z reference adjustments
	bs_idx = np.ones(idx_shape)  # idz of reported backscatter (dB)
	dz_ref_idx = np.ones(idx_shape)  # idx of dz from ref (meters)
	dz_ref_wd_idx = np.ones(idx_shape)  # idx of dz from ref (% WD)
	depth_mode_idx = np.ones(idx_shape)  # idx of depth mode matching combo box
	real_idx = np.logical_not(np.logical_or(np.isnan(self.xline['z_final']),
											np.isnan(self.xline['z_final'])))  # idx true for NON-NAN soundings

	if print_updates:
		print('number of xline soundings:', len(self.xline['z_final']))
		print('number of nans found in xline z_final=', np.sum(np.logical_not(real_idx)))

	if self.angle_xline_gb.isChecked():  # get idx satisfying current swath angle filter based on depth/acrosstrack angle
		lims = [float(self.min_angle_xline_tb.text()), float(self.max_angle_xline_tb.text())]
		# swath coverage method uses abs(angle) to accomplish symmetrical angle filtering (e.g., excluding nadir)
		# angle_idx = np.logical_and(np.abs(np.asarray(self.xline['beam_angle'])) >= lims[0],
		# 						   np.abs(np.asarray(self.xline['beam_angle'])) <= lims[1])

		# swath accuracy method uses angle including sign to filter actual beam angles (e.g., port and stbd outer swath)
		angle_idx = np.logical_and(np.asarray(self.xline['beam_angle']) >= lims[0],
								   np.asarray(self.xline['beam_angle']) <= lims[1])

	if self.depth_xline_gb.isChecked():  # get idx satisfying current depth filter
		lims = [-1*float(self.max_depth_xline_tb.text()), -1*float(self.min_depth_xline_tb.text())]
		depth_idx = np.logical_and(np.asarray(self.xline['z_final']) >= lims[0],
								   np.asarray(self.xline['z_final']) <= lims[1])

	# Separate absolute depth difference filter
	if hasattr(self, 'dz_abs_gb') and self.dz_abs_gb.isChecked():
		dz_abs_idx = np.abs(self.xline['dz_ref']) <= float(self.max_dz_tb.text())
	else:
		dz_abs_idx = None
	
	# Separate percentage of water depth filter
	if hasattr(self, 'dz_pct_gb') and self.dz_pct_gb.isChecked():
		dz_pct_idx = np.abs(self.xline['dz_ref_wd']) <= float(self.max_dz_wd_tb.text())
	else:
		dz_pct_idx = None

	if self.bs_xline_gb.isChecked():  # get idx satisfying current backscatter filter
		lims = [float(self.min_bs_xline_tb.text()), float(self.max_bs_xline_tb.text())]  # parsed BS is converted to dB
		bs_idx = np.logical_and(np.asarray(self.xline['bs']) >= lims[0],
								np.asarray(self.xline['bs']) <= lims[1])

	if self.depth_mode_gb.isChecked():  # get idx satisfying current depth mode filter
		depth_mode_filter = self.depth_mode_cbox.currentText().strip().lower()
		update_log(self, 'Crosslines will be filtered by depth mode: ' + depth_mode_filter)
		# ping mode may be, e.g., 'Deeper (Manual)', need to compare just mode and not whether Manual
		depth_mode_idx = np.asarray([pm.split('(')[0].strip().lower() == depth_mode_filter for
									 pm in self.xline['ping_mode']])

	# # apply filter masks to x, z, angle, and bs fields
	# self.xline['filter_idx'] = np.logical_and.reduce((angle_idx, depth_idx, bs_idx, dz_ref_idx, dz_ref_wd_idx))
	
	# Ensure all filter indices are boolean arrays with the same shape
	try:
		# Convert all indices to boolean arrays and ensure they have the same shape
		filter_arrays = []
		for idx_array in [angle_idx, depth_idx, bs_idx, dz_abs_idx, dz_pct_idx, depth_mode_idx]:
			if idx_array is not None:
				# Convert to numpy array and ensure boolean type
				bool_array = np.asarray(idx_array, dtype=bool)
				filter_arrays.append(bool_array)
			else:
				# If any index is None, create a True array of the correct shape
				bool_array = np.ones(idx_shape, dtype=bool)
				filter_arrays.append(bool_array)
		
		# Apply all filters
		self.xline['filter_idx'] = np.logical_and.reduce(filter_arrays)
		
	except Exception as e:
		print(f"Error in filter_xline: {e}")
		# Fallback: create a simple True array if filtering fails
		self.xline['filter_idx'] = np.ones(idx_shape, dtype=bool)

	# print('number of soundings passing filter: ', np.sum(self.xline['filter_idx']))

def parse_tide(self, unit_set_by_user=False):  # force_refresh=False):
	# add tide file if available -
	fnames_tide = get_new_file_list(self, ['.tid', '.txt'], [])  # list .tid and .txt files
	# print('fnames_tide is', fnames_tide)

	if len(fnames_tide) != 1:  # warn user to add exactly one tide file
		update_log(self, 'Add one tide text (.tid or .txt) file corresponding to the accuracy crossline location and '
						 'covering the entire crossline data duration.\n\n'
						 'Format should be ''YYYY/MM/DD[/]HH:MM[:SS] tide_amplitude'' with time in the same zone as '
						 'the crossline data (e.g., UTC/GMT)\n\n'
						 'Amplitude is positive up in meters or feet; the vertical datum for the tide file must '
						 'match that of the reference surface for meaningful comparison (e.g., same tide applied during'
						 'processing of the refernce surface should be applied here for the crosslines).')
		pass

	else:
		tic = process_time()
		update_log(self, 'Processing tide file')
		fname_tide = fnames_tide[0]

		self.tide = {}
		self.tide['fname'] = fname_tide.rsplit('/', 1)[1]
		# print('storing tide fname =', self.tide['fname'])
		self.tide['time_obj'], self.tide['amplitude'] = [], []

		# try to get units from filename (assumed meters, look for feet, etc.)
		tide_unit = ''
		try:  # look for any units specified in the filename
			if unit_set_by_user:  # use unit selected by user (override
				update_log(self, 'Tide file amplitude unit set by user to ' + self.tide_unit_cbox.currentText())
				tide_unit = self.tide_unit_cbox.currentText()

			else:  # look for unit in filename, update cbox for user to review/confirm units
				if any([self.tide['fname'].lower().find(u) > -1 for u in ['ft', 'feet', 'foot']]):
					tide_unit = 'Foot'

				elif any([self.tide['fname'].lower().find(u) > -1 for u in ['meter', '_m.']]):
					tide_unit = 'Meter'  # found meters, try to avoid simple m character

				elif self.tide['fname'].lower().find('smoot') > -1 and \
						self.tide['fname'].lower().find('smooth') == -1:
					tide_unit = 'Smoot'  # last case, found smoot but not smooth

				# try to update the combobox if a unit was found in the file and not manually set by user
				if tide_unit:
					update_log(self, 'Found tide unit ' + tide_unit.upper() + ' in tide filename; please verify or update in tide '
																			  'unit selection menu and update units as necessary')
					if tide_unit != 'Meter':
						update_log(self, 'All data will be converted to meters for analysis')

				else:
					tide_unit = 'Meter'
					update_log(self, 'Tide unit not detected in filename; assumed tide amplitude in METERS')

				print('Looking for matching cbox index...')
				tide_unit_idx = self.tide_unit_cbox.findText(tide_unit)
				print('got index = ', tide_unit_idx)
				if tide_unit_idx >= 0:
					print('trying to set tide unit cbox to detected unit')
					self.tide_unit_cbox.setCurrentIndex(tide_unit_idx)
					print('great success')

		except:
			print('Tide unit search failed in filename: ', self.tide['fname'])

		with open(fname_tide, 'r') as fid_tide:  # read each line of the tide file, strip newline
			tide_list = [line.strip().rstrip() for line in fid_tide]

		fid_tide.close()
		temp_time = []
		temp_tide = []

		time_fmt = ''  # start without assuming tide time format
		print('starting tide time format loop')

		for l in tide_list:
			print('working on l =', l)
			try:  # try parsing and converting each line before adding to temp time and tide lists
				# print('in first try statement')
				# print('l.rsplit = ', l.replace('\t', ' ').rsplit(' ', 1))

				# time_str, amp = l.replace('\t', ' ').rsplit(' ', 1)  # separate tide amp at end from uncertain time str
				# time_str, amp = l.replace('\t', ' ').rsplit('.', 1)

				# TPXO format
				# 37.3730 -123.1850	06.26.2021	01:00:00 -0.056 1299.549

				# some NOAA exports have spaces mixed into the amplitude decimals, which need to be removed to get the
				# proper amplitude

				# print('len(l.split( ) = ', len(l.split(' ')), 'with fields', l.split(' '))
				print('len(l.split( ) = ', len(l.split()), 'with fields', l.split())

				n_fields = len(l.split())
				tpxo_format = False
				is_yyyy_mm_dd_format = False  # Initialize the flag

				if n_fields == 6:  # probably TPXO download format if lots of fields (more than date time amp)
					# TPXO format
					# 37.3730 -123.1850	06.26.2021	01:00:00 -0.056 1299.549
					print('checking TPXO download format')
					try:
						lat, lon, date, time, amp, depth = l.split()
						time_str = ' '.join([date, time]).strip()
						amp = amp.strip()
						print('got results: ', lat, lon, date, time, amp, depth)
						tpxo_format = True
					except:
						print('failed TPXO format')

				elif n_fields == 3:  # probably standard .tid format, but may have many different date time separators (/, :, ., etc.)
					# Check if it's YYYY/MM/DD HH:MM:SS amplitude format
					fields = l.split()
					is_yyyy_mm_dd_format = len(fields) == 3 and '/' in fields[0] and ':' in fields[1]
					if is_yyyy_mm_dd_format:
						# Format: 2025/07/04 17:00:00 0.176
						time_str = ' '.join([fields[0], fields[1]]).strip()
						amp = fields[2].strip()
						print('detected YYYY/MM/DD HH:MM:SS format')
					else:
						# Original logic for decimal split
						part1, part2 = l.rsplit('.', 1)  # split line at the last decimal
						# print('got part1, part2 = ', part1, part2)
						time_str, amp_int = part1.replace('\t', ' ').rsplit(' ', 1)  # split time str from amplitude whole num
						print('time str after first split is', time_str)
						# time_str2 = ' '.join(time_str.rsplit(' ')[-2:])  # exclude lat lon prior to time in TPXO download format
						# print('time str is now ', time_str2)
						# print('got time_str, amp_int =', time_str, amp_int)
						amp_dec = part2.replace(' ', '')  # remove any spaces in decimal part
						# print('got amp_dec =', amp_dec)
						amp = '.'.join([amp_int, amp_dec])  # rejoin whole number and decimal part of amplitude
						# print('got time_str, amp_int, and amp_dec = ', time_str, amp)

				else:
					update_log(self, 'Tide format not recognized. Please load .tid format or TPXO download tide file.')

				# For YYYY/MM/DD HH:MM:SS format, parse date and time separately
				if is_yyyy_mm_dd_format:
					# Parse YYYY/MM/DD HH:MM:SS format by splitting date and time
					date_part, time_part = time_str.split(' ')
					year, month, day = date_part.split('/')
					hour, minute, second = time_part.split(':')
					# Create datetime object directly
					from datetime import datetime
					dt = datetime(int(year), int(month), int(day), int(hour), int(minute), int(second))
					time_str_reduced = dt.strftime('%Y%m%d%H:%M:%S')  # Convert to standard format without space
					print(f'Parsed YYYY/MM/DD HH:MM:SS format: {time_str} -> {time_str_reduced}')
				else:
					# time_str_reduced = re.sub('[^0-9^.]', '', time_str).strip()  # remove alpha (day name) and :, leave . and ms
					time_str_reduced = re.sub('[^0-9^.^:^ ]', '', time_str).strip()  # remove alpha (day name), leave . and : and space and ms

				print('got time_str --> reduced =', time_str, ' --> ', time_str_reduced)

				# if time_fmt:  # if successful time format is known, try that
				try:  # try parsing time string w/ last successful format (first line will fail until format is found)
					dt = datetime.strptime(time_str_reduced, time_fmt)
					# print('got datetime =', dt)

				except:  # figure out time format and use this format for future lines until it fails
					# list of formats in increasing complexity (may parse incorrectly if starting with most complex)
					# fmt_list = ['%Y%m%d%H%M', '%Y%m%d%H%M%S', '%Y%m%d%H%M%S.%f']  # reduced time_str should match one...
					# fmt_list = ['%Y%m%d%H%M', '%Y%m%d%H%M%S', '%Y%m%d%H%M%S.%f', '%m.%d.%Y%H%M%S']  # reduced time_str should match one...
					fmt_list = ['%Y%m%d%H:%M', '%Y%m%d%H:%M:%S', '%Y%m%d%H:%M:%S.%f', '%m.%d.%Y%H:%M:%S', '%Y/%m/%d %H:%M:%S', '%Y%m%d%H:%M:%S']  # reduced time_str should match one...

					print('time format = ', time_fmt, ' did not work; trying other formats')
					for fmt in fmt_list:
						print('trying tide format', fmt)
						try:
							print('trying to parse time with ', fmt)

							# TPXO download format includes 24:00 that should be converted to 00:00 (and add a day)
							# if fmt == '%m.%d.%Y%H:%M:%S' and time_str_reduced.find('24:') > -1:
							if tpxo_format and time_str_reduced.find('24:') > -1:
								print('found a 24: in time_str_reduced =', time_str_reduced)
								dt = datetime.strptime(time_str_reduced.replace('24:', '00:'), fmt)
								dt = dt + timedelta(days=1)

							else:
								dt = datetime.strptime(time_str_reduced, fmt)

							print('successfully parsed with format: ', fmt)
							time_fmt = fmt  # store format if successfully parsed this time string
							print('set time_fmt =', time_fmt)

							update_log(self, 'Found tide time format: ' + time_fmt)
							print('breaking!')

							break

						except:
							print('Tide time format ', fmt, ' did not work for time_str_reduced =', time_str,
								  '--> trying next format in list')

				if dt:  # try to parse the amplitude only if time was successfully parse
					try:
						amp = float(amp)
						temp_time.append(dt)
						temp_tide.append(amp)

					except:
						print('failed to convert amp to float: ', amp)

				else:
					print('Time was not parsed; skipping amplitude')

			except:
				print('failed to parse tide file line =', l, '(possible header)')


		tide_amp_fac = self.tide_unit_dict[self.tide_unit_cbox.currentText()]  # apply selected/updated amplitude unit
		print('for tide units = ', self.tide_unit_cbox.currentText(), ' got amp fac =', tide_amp_fac)

		print('temp_tide =', temp_tide)

		tide_m = np.multiply(temp_tide, tide_amp_fac)

		self.tide['time_obj'] = deepcopy(temp_time)
		self.tide['amplitude'] = deepcopy(tide_m)
		self.tide_applied = False
		# print('final self.tide time and amp are ', self.tide['time_obj'], self.tide['amplitude'])

		# if force_refresh:
		# 	plot_tide(self, set_active_tab=True)
		# 	refresh_plot(self, refresh_list=['tide'], set_active_tab=3, sender='add_tide_file')


def parse_ref_depth(self):
	# parse the loaded reference surface .xyz file; assumes all units are meters in UTM projection
	# .xyz file is assumed comma or space delimited with minimum fields: east, north, depth (+Z up)
	# optional fourth field for uncertainty (e.g., Qimera CUBE surface export), set to 0 if not included

	self.ref = {}
	fnames_xyz = get_new_file_list(self, ['.xyz'], [])  # list .xyz files
	print('fnames_xyz is', fnames_xyz)

	if len(fnames_xyz) != 1:  # warn user to add exactly one ref grid
		update_log(self, 'Please add one reference grid .xyz file in a UTM projection')
		pass

	else:
		fname_ref = fnames_xyz[0]
		self.ref['fname'] = fname_ref.rsplit('/', 1)[1]
		print(fname_ref)
		fid_ref = open(fname_ref, 'r')
		e_ref, n_ref, z_ref, u_ref = [], [], [], []

		for line in fid_ref:
			# strip and split space- or comma-delimited line; append '0' to list in case uncertainty field is not avail
			temp = line.replace(',', ' ').strip().rstrip().split() + ['0']
			e_ref.append(temp[0])  # easting
			n_ref.append(temp[1])  # northing
			z_ref.append(temp[2])  # up
			u_ref.append(temp[3])  # uncertainty (fourth value in line if included; 0 if not)

	print('*** finished parsing .xyz, got first ten uncertainty values:', u_ref[0:10])

	# update log about uncertainty
	update_log(self, 'Uncertainty ' + ('not ' if len(set(u_ref)) == 1 and u_ref[0] == '0' else '') + 'parsed from .xyz')

	# convert to arrays with Z positive up; vertical datum for ref grid and crosslines is assumed same for now
	self.ref['e'] = np.array(e_ref, dtype=np.float32)
	self.ref['n'] = np.array(n_ref, dtype=np.float32)
	self.ref['z'] = -1 * np.abs(np.array(z_ref, dtype=np.float32))  # ensure grid is +Z UP (neg. depths)
	self.ref['utm_zone'] = self.ref_proj_cbox.currentText()
	self.ref['u'] = np.array(u_ref, dtype=np.float32)


	update_log(self, 'Imported ref grid: ' + fname_ref.split('/')[-1] + ' with ' +
			   str(len(self.ref['z'])) + ' nodes')
	update_log(self, 'Ref grid is assigned UTM zone ' + self.ref['utm_zone'])

	ref_de, ref_dn = check_cell_size(self, self.ref['e'], self.ref['n'])

	if ref_de == ref_dn:
		self.ref_cell_size = ref_de
		update_log(self, 'Imported ref grid has uniform cell size: ' + str(self.ref_cell_size) + ' m')

	else:
		self.ref_cell_size = np.max([ref_de, ref_dn])
		update_log(self, 'WARNING: Unequal ref grid cell sizes (E: ' + str(ref_de) + ' m , N: ' + str(ref_dn) + ' m)',
				   font_color="red")

	print('leaving parse_ref_depth with self.ref.keys =', self.ref.keys())


def check_cell_size(self, easting, northing):
	de = np.mean(np.diff(np.sort(np.unique(easting))))
	dn = np.mean(np.diff(np.sort(np.unique(northing))))

	return de, dn


def validate_density_surface_compatibility(self, e_dens, n_dens, c_dens):
	"""
	Validate that the density surface is compatible with the loaded depth grid.
	
	Args:
		e_dens: List of easting coordinates from density file
		n_dens: List of northing coordinates from density file  
		c_dens: List of count values from density file
		
	Returns:
		tuple: (is_compatible, error_message)
	"""
	try:
		# Check if reference surface exists
		if 'z' not in self.ref or 'e' not in self.ref or 'n' not in self.ref:
			return False, "Reference surface must be loaded before density surface"
		
		# Check if we have density data
		if len(e_dens) == 0 or len(n_dens) == 0 or len(c_dens) == 0:
			return False, "No valid data found in density file"
		
		# Check cell size compatibility
		dens_de, dens_dn = check_cell_size(self, e_dens, n_dens)
		
		# Check if cell sizes match (within 1% tolerance)
		cell_size_tolerance = 0.01  # 1% tolerance
		if abs(dens_de - self.ref_cell_size) / self.ref_cell_size > cell_size_tolerance or \
		   abs(dens_dn - self.ref_cell_size) / self.ref_cell_size > cell_size_tolerance:
			return False, f"Density grid cell size ({dens_de:.1f}m x {dens_dn:.1f}m) does not match reference grid cell size ({self.ref_cell_size:.1f}m). Please use a density surface with the same cell size."
		
		# Check spatial extent compatibility
		ref_e_min, ref_e_max = min(self.ref['e']), max(self.ref['e'])
		ref_n_min, ref_n_max = min(self.ref['n']), max(self.ref['n'])
		
		dens_e_min, dens_e_max = min(e_dens), max(e_dens)
		dens_n_min, dens_n_max = min(n_dens), max(n_dens)
		
		# Calculate overlap
		e_overlap = min(ref_e_max, dens_e_max) - max(ref_e_min, dens_e_min)
		n_overlap = min(ref_n_max, dens_n_max) - max(ref_n_min, dens_n_min)
		
		# Check if there's significant overlap (at least 50% of the smaller grid)
		ref_e_extent = ref_e_max - ref_e_min
		ref_n_extent = ref_n_max - ref_n_min
		dens_e_extent = dens_e_max - dens_e_min
		dens_n_extent = dens_n_max - dens_n_min
		
		min_e_extent = min(ref_e_extent, dens_e_extent)
		min_n_extent = min(ref_n_extent, dens_n_extent)
		
		if e_overlap <= 0 or n_overlap <= 0:
			return False, f"Density surface does not overlap with reference surface. Reference extent: E({ref_e_min:.0f}-{ref_e_max:.0f}), N({ref_n_min:.0f}-{ref_n_max:.0f}). Density extent: E({dens_e_min:.0f}-{dens_e_max:.0f}), N({dens_n_min:.0f}-{dens_n_max:.0f})"
		
		# Check if overlap is significant (at least 50% of the smaller grid)
		e_overlap_ratio = e_overlap / min_e_extent
		n_overlap_ratio = n_overlap / min_n_extent
		
		if e_overlap_ratio < 0.5 or n_overlap_ratio < 0.5:
			return False, f"Insufficient overlap between density and reference surfaces. Overlap ratios: E={e_overlap_ratio:.1%}, N={n_overlap_ratio:.1%}. Minimum required: 50%"
		
		return True, "Density surface is compatible with reference surface"
		
	except Exception as e:
		return False, f"Error validating density surface compatibility: {str(e)}"


def remove_incompatible_density_surface(self, fname_dens):
	"""
	Remove an incompatible density surface from the sources list and clean up any associated data.
	
	Args:
		fname_dens: Path to the incompatible density surface file
	"""
	try:
		# Find the file in the file list and remove it
		for i in range(self.file_list.count()):
			item = self.file_list.item(i)
			if item and fname_dens in item.text():
				self.file_list.takeItem(i)
				update_log(self, f'Removed incompatible density surface: {os.path.basename(fname_dens)}')
				break
		
		# Remove from filenames list
		if fname_dens in self.filenames:
			self.filenames.remove(fname_dens)
		
		# Clean up any density data from ref dict
		for k in ['fname_dens', 'c', 'c_grid', 'c_ref_extent']:
			if k in self.ref:
				self.ref.pop(k, None)
		
		# Re-enable the density surface button
		self.add_dens_surf_btn.setEnabled(True)
		
		# Refresh the plot to remove density-related displays
		refresh_plot(self, refresh_list=['ref'], sender='remove_incompatible_density_surface')
		
	except Exception as e:
		update_log(self, f'WARNING: Failed to remove incompatible density surface: {str(e)}', font_color="orange")
		print(f'Error removing incompatible density surface: {e}')


def parse_ref_dens(self):
	# add density surface if available - this is useful for Qimera .xyz files that do not include density
	fnames_xyd = get_new_file_list(self, ['.xyd'], [])  # list .xyz files
	print('fnames_xyd is', fnames_xyd)

	if len(fnames_xyd) == 0:  # no density files
		update_log(self, 'No density surface files found')
		return
	elif len(fnames_xyd) > 1:  # multiple density files - use the most recent one
		update_log(self, f'Multiple density files found ({len(fnames_xyd)}). Using the most recent one.')
		# Sort by modification time and use the most recent
		fnames_xyd.sort(key=lambda x: os.path.getmtime(x), reverse=True)
		fname_dens = fnames_xyd[0]
		update_log(self, f'Using density file: {os.path.basename(fname_dens)}')
	else:  # exactly one density file
		fname_dens = fnames_xyd[0]
	
	try:
		tic = process_time()
		update_log(self, 'Processing reference surface density')
		print('parsing/matching density data')
		
		# Check if reference surface exists
		if 'z' not in self.ref or 'e' not in self.ref or 'n' not in self.ref:
			update_log(self, 'ERROR: Reference surface must be loaded before density surface', font_color="red")
			return
		
		# Check if file exists and is readable
		if not os.path.exists(fname_dens):
			update_log(self, f'ERROR: Density file not found: {fname_dens}', font_color="red")
			return
			
		self.ref['fname_dens'] = fname_dens
		e_dens, n_dens, c_dens = [], [], []
		
		# Parse XYD file with error handling
		line_count = 0
		try:
			with open(fname_dens, 'r') as fid_dens:
				for line in fid_dens:
					line_count += 1
					try:
						# Skip empty lines
						if not line.strip():
							continue
							
						temp = line.replace(',', ' ').strip().rstrip().split()
						if len(temp) >= 3:  # Ensure we have at least 3 columns
							e_dens.append(float(temp[0]))  # easting
							n_dens.append(float(temp[1]))  # northing
							c_dens.append(float(temp[2]))  # count
						else:
							print(f'Warning: Skipping line {line_count} - insufficient columns')
					except (ValueError, IndexError) as e:
						print(f'Warning: Skipping line {line_count} - parsing error: {e}')
						continue
		except Exception as e:
			update_log(self, f'ERROR: Failed to read density file: {str(e)}', font_color="red")
			return
		
		if len(e_dens) == 0:
			update_log(self, 'ERROR: No valid data found in XYD file', font_color="red")
			return

		# Convert to numpy arrays with error handling
		try:
			e_dens = np.array(e_dens, dtype=np.float32)
			n_dens = np.array(n_dens, dtype=np.float32)
			c_dens = np.array(c_dens, dtype=np.float32)
		except Exception as e:
			update_log(self, f'ERROR: Failed to convert density data to arrays: {str(e)}', font_color="red")
			return

		# Validate density surface compatibility
		is_compatible, error_message = validate_density_surface_compatibility(self, e_dens, n_dens, c_dens)
		if not is_compatible:
			update_log(self, f'ERROR: {error_message}. Please remove the file from the sources list.', font_color="red")
			# Remove the incompatible density surface from sources
			remove_incompatible_density_surface(self, fname_dens)
			return

		# Log successful validation
		update_log(self, 'Density surface validation passed')

		# loop through all nodes in reference surface, attach density value if matching E and N
		try:
			self.ref['c'] = np.empty_like(self.ref['z'])
			self.ref['c'][:] = np.nan
			len_ref = len(self.ref['z'])
			parse_prog_old = -1
			matched_count = 0

			update_log(self, f'Starting density matching for {len_ref} reference nodes...')
			
			for i in range(len(self.ref['z'])):
				parse_prog = round(10 * i / len_ref)
				if parse_prog > parse_prog_old:
					print("%s%%" % (parse_prog * 10) + ('\n' if parse_prog == 10 else ''), end=" ", flush=True)
					parse_prog_old = parse_prog

				try:
					idx_e = np.where(e_dens == self.ref['e'][i])
					idx_n = np.where(n_dens == self.ref['n'][i])
					idx_match = np.intersect1d(idx_e, idx_n)
					if len(idx_match) == 1:
						self.ref['c'][i] = c_dens[idx_match][0]
						matched_count += 1
				except Exception as e:
					print(f'Warning: Error matching node {i}: {e}')
					continue

			update_log(self, f'Imported density grid: {fname_dens.split("/")[-1]} with {len(self.ref["c"])} nodes, {matched_count} matched')

			toc = process_time()
			refresh_time = toc - tic
			print('assigning density nodes took', refresh_time, ' s')
			
		except Exception as e:
			update_log(self, f'ERROR: Failed to match density data: {str(e)}', font_color="red")
			# Clean up any partial data
			if 'c' in self.ref:
				del self.ref['c']
			if 'fname_dens' in self.ref:
				del self.ref['fname_dens']
			return
		
	except Exception as e:
		update_log(self, f'ERROR: Failed to process density file: {str(e)}', font_color="red")
		print(f'Density processing error: {e}')
		# Clean up any partial data
		if 'c' in self.ref:
			del self.ref['c']
		if 'fname_dens' in self.ref:
			del self.ref['fname_dens']


def update_ref_slope(self):
	# update slope and plot after changing slope calc params
	calc_ref_slope(self)
	refresh_plot(self, refresh_list=['ref'], sender='update_ref_slope')


def calc_ref_slope(self):
	# calculate representative maximum slope for each node in reference grid
	# 0. make z grid with nearest neighbor interpolation (so empty cells and edges do not cause NaN gradients; this is
	# for gradient/slope calculation only, as these cells will be masked later to match shape of depth grid)
	print('in calc_ref_slope')
	z_grid_temp, z_grid_temp_extent = make_grid(self, self.ref['e'], self.ref['n'], self.ref['z'], self.ref_cell_size,
												mask_original_shape=False, method='nearest')

	# 1. take moving average of z grid with user-selected window size
	z_grid_smooth = uniform_filter(z_grid_temp, size=int(self.slope_win_cbox.currentText()[0]))

	# 2. calculate north-south and east-west gradients; convert to slopes with known cell spacing; assume the maximum
	# absolute value of slope in either direction is representative of the maximum slope in any direction for each node
	grad = np.gradient(z_grid_smooth)
	slope_y = np.rad2deg(np.arctan2(grad[0], self.ref_cell_size))
	slope_x = np.rad2deg(np.arctan2(grad[1], self.ref_cell_size))
	slope_max = np.maximum(np.abs(slope_x), np.abs(slope_y))

	# 3. apply mask from original z_grid (NaN wherever no depth node exists
	slope_max[np.isnan(self.ref['z_grid'])] = np.nan
	self.ref['s_grid'] = slope_max  # store slope grid with same shape as depth and density grids
	self.ref['s'] = slope_max[~np.isnan(slope_max)][:]  # store slopes like easting, northing, depth, dens for filtering
	print('leaving calc_ref_slope')


def make_ref_surf(self):
	# make grids of reference surface depth, density, slope
	# note: this is separate from applying masks from user limits, e.g., adjustable slope filters
	from time import process_time
	start_time = process_time()
	
	print('calling make_grid from make_ref_surf')
	update_log(self, '=== GENERATING REFERENCE SURFACE GRIDS ===')

	try:
		# Check if reference surface data exists
		if 'z' not in self.ref or 'e' not in self.ref or 'n' not in self.ref:
			update_log(self, 'ERROR: Reference surface data not available', font_color="red")
			return
			
		# Check if cell size is defined
		if not hasattr(self, 'ref_cell_size') or self.ref_cell_size is None:
			update_log(self, 'ERROR: Reference cell size not defined', font_color="red")
			return

		update_log(self, f'Cell size: {self.ref_cell_size} meters')
		update_log(self, f'Reference surface contains {len(self.ref["z"])} nodes')

		# for dim in ['n', 'e', 'z', 'c']:
		for dim in ['n', 'e', 'z', 'c', 'u']:

			grid_str = dim + '_grid'
			extent_str = dim + '_ref_extent'

			if dim in self.ref and not grid_str in self.ref:  # generate grid, extent for parameter if not done so already
				update_log(self, f'Generating grid for {dim} dimension...')
				print('in make_ref_surf, calling make_grid for dim =', dim, 'with self.ref_cell_size =', self.ref_cell_size)
				
				grid_start = process_time()
				try:
					grid, extent = make_grid(self, self.ref['e'], self.ref['n'], self.ref[dim], self.ref_cell_size)
					
					# Check if make_grid succeeded
					if grid is not None and extent is not None:
						self.ref[grid_str] = deepcopy(grid)
						self.ref[extent_str] = deepcopy(extent)
						grid_time = process_time() - grid_start
						update_log(self, f'Grid generation for {dim} completed in {grid_time:.2f} seconds')
					else:
						print(f'ERROR: Failed to generate grid for dimension {dim}')
						update_log(self, f'ERROR: Failed to generate grid for {dim}', font_color="red")
						return
				except Exception as e:
					print(f'ERROR: Exception in make_grid for dimension {dim}: {e}')
					update_log(self, f'ERROR: Exception in make_grid for {dim}: {str(e)}', font_color="red")
					return
			else:
				if grid_str in self.ref:
					update_log(self, f'Grid for {dim} already exists, skipping generation')
				else:
					update_log(self, f'No data available for {dim} dimension, skipping grid generation')
				print('grid not generated for dim =', dim, 'because either grid already exists or data not available')

		if 'z_grid' in self.ref and not 's_grid' in self.ref:  # make slope grid masked with original shape of depth grid
			update_log(self, 'Generating reference slope grid...')
			slope_start = process_time()
			try:
				calc_ref_slope(self)  # calculate slope using z_grid
				slope_time = process_time() - slope_start
				update_log(self, f'Slope grid generation completed in {slope_time:.2f} seconds')
			except Exception as e:
				print(f'ERROR: Exception in calc_ref_slope: {e}')
				update_log(self, f'ERROR: Exception in calc_ref_slope: {str(e)}', font_color="red")
				return

		total_time = process_time() - start_time
		update_log(self, f'=== REFERENCE SURFACE GRID GENERATION COMPLETED in {total_time:.2f} seconds ===')
		print('leaving make_ref_surf with self.ref.keys = ', self.ref.keys())
		
	except Exception as e:
		print(f'ERROR: Exception in make_ref_surf: {e}')
		update_log(self, f'ERROR: Exception in make_ref_surf: {str(e)}', font_color="red")


def make_grid(self, ax1, ax2, val, spacing, mask_original_shape=True, fill_value=np.nan, method='linear'):
	# generate a grid covering the ax1 and ax2 extents with NaN wherever there is no val
	# ax1, ax2, and val are lists of, e.g., easting, northing, and depth values with same len
	# spacing is the assumed fixed grid size desired in the output grid
	# by default, return Nans where there is no data in the gridded result; set mask_original_shape to False to
	# keep the interpolated result within the convex hull and NaNs outside; set method to 'nearest' to fill grid with
	# nearest value

	try:
		# Check input data validity
		if len(ax1) == 0 or len(ax2) == 0 or len(val) == 0:
			print('ERROR: Empty input arrays to make_grid')
			return None, None
			
		if len(ax1) != len(ax2) or len(ax1) != len(val):
			print('ERROR: Input arrays to make_grid have different lengths')
			return None, None

		# generate a z_grid covering the E and N extents with NaN wherever there is no imported ref depth
		ax1_min, ax1_max, ax2_min, ax2_max = min(ax1), max(ax1), min(ax2), max(ax2)
		
		try:
			print('in make_grid, ax1_min, max =', ax1_min, ax1_max, ' and ax2_min, max =', ax2_min, ax2_max,
				  'and spacing =', spacing)
			num_ax1 = round((ax1_max - ax1_min) / spacing) + 1  # number of easting nodes in grid if all are spaced equally
			num_ax2 = round((ax2_max - ax2_min) / spacing) + 1  # number of northing nodes in grid if all are spaced equally
			print('in make_grid, got rounded num_ax1, num_ax2 =', num_ax1, num_ax2)
			
			# Check for reasonable grid size to prevent memory issues
			if num_ax1 * num_ax2 > 1000000:  # More than 1 million grid points
				print(f'WARNING: Large grid size ({num_ax1} x {num_ax2} = {num_ax1 * num_ax2} points). This may cause memory issues.')
				update_log(self, f'WARNING: Large grid size may cause memory issues', font_color="orange")
				
		except Exception as e:
			print('in make_grid, failed to get num_ax1 and/or num_ax2 from ax min/max and spacing info:', e)
			return None, None

		ax1_nodes, ax2_nodes = np.linspace(ax1_min, ax1_max, num_ax1), np.linspace(ax2_min, ax2_max, num_ax2)
		ax1_grid, ax2_grid = np.meshgrid(ax1_nodes, ax2_nodes, indexing='xy')  # generate arrays for all nodes for griddata

		# grid data (griddata interps within convex hull; remove results outside 'original' shape in next step if desired)
		try:
			val_grid = griddata((ax1, ax2), val, (ax1_grid, ax2_grid), method=method, fill_value=fill_value)
		except Exception as e:
			print(f'ERROR: griddata failed: {e}')
			return None, None

		if mask_original_shape:
			try:
				# mask griddata results that include convex hull areas outside the imported ref surface data
				# (adapted from SO answer by pv. for fast masking of griddata when convex hull is insufficient)
				tree = KDTree(np.c_[ax1, ax2])
				dist, _ = tree.query(np.c_[ax1_grid.ravel(), ax2_grid.ravel()], k=1)
				dist = dist.reshape(ax1_grid.shape)
				# val_grid[dist > 0.1] = np.nan  # old fixed method lead to strange masking results for imperfectly spaced grids
				print('masking val_grid for original grid shape, setting nans for dist > 1% of grid spacing')
				val_grid[dist > spacing/100] = np.nan
			except Exception as e:
				print(f'WARNING: Masking failed, continuing without masking: {e}')

		val_grid = np.flipud(val_grid)  # flipud required to get same orientation as ref surf in proc software
		extent = ax1_min, ax1_max, ax2_min, ax2_max

		print('returning from make_grid with size of val_grid = ', np.size(val_grid), 'and extents =', extent)
		print('first ten values of val_grid are', val_grid[0:10])

		return val_grid, extent
		
	except Exception as e:
		print(f'ERROR: make_grid failed: {e}')
		return None, None


def calc_ref_mask(self):
	# calculate the final reference surface based on depth, density, and slope filters
	z_min = (-1*float(self.max_depth_ref_tb.text()) if self.depth_ref_gb.isChecked() else -1*np.inf)  # +Z down in GUI
	z_max = (-1*float(self.min_depth_ref_tb.text()) if self.depth_ref_gb.isChecked() else np.inf)  # +Z down in GUI
	c_min, c_max = (float(self.min_dens_tb.text()) if self.density_gb.isChecked() else 0), np.inf  # sounding count
	s_min, s_max = 0, (float(self.max_slope_tb.text()) if self.slope_gb.isChecked() else np.inf)  # slope
	u_min, u_max = 0, (float(self.max_u_tb.text()) if self.uncertainty_gb.isChecked() else np.inf)  # uncertainty

	# print('MASKING WITH min/max of z, c, s, and u=', z_min, z_max, c_min, c_max, s_min, s_max, u_min, u_max)

	try:
		print('IN CALC_REF_MASK, self.ref.keys =', self.ref.keys())
		if all([k in self.ref.keys() for k in ['z_grid', 's_grid']]):  # for k in self.ref.keys()]):
			print('in calc_ref_mask, found z_grid and s_grid in self.ref.keys')
			# initialize masks for each grid, true unless filtered
			for mask in ['z_mask', 's_mask', 'c_mask', 'u_mask', 'final_mask']:
				self.ref[mask] = np.nan*np.ones_like(self.ref['z_grid'])

			print('in calc_ref_mask, applying depth and slope criteria to mask')
			self.ref['z_mask'][np.logical_and(self.ref['z_grid'] >= z_min, self.ref['z_grid'] <= z_max)] = 1
			self.ref['s_mask'][np.logical_and(self.ref['s_grid'] >= s_min, self.ref['s_grid'] <= s_max)] = 1

			# uncertainty is set to 0 if not parsed from .xyz file
			print('working on u_mask...')
			self.ref['u_mask'][np.logical_and(self.ref['u_grid'] >= u_min, self.ref['u_grid'] <= u_max)] = 1

			if 'c_grid' in self.ref:  # sounding count may not be loaded
				self.ref['c_mask'][np.logical_and(self.ref['c_grid'] >= c_min, self.ref['c_grid'] <= c_max)] = 1

			else:
				self.ref['c_mask'] = np.ones_like(self.ref['z_grid'])

			self.ref['final_mask'] = self.ref['z_mask']*self.ref['s_mask']*self.ref['c_mask']*self.ref['u_mask']

			for mask in ['z_mask', 's_mask', 'c_mask', 'u_mask', 'final_mask']:
				print('num in ', mask, '=', np.sum(~np.isnan(self.ref[mask])))

		else:
			print('z_grid and/or s_grid not found in self.ref')

	except:
		print('error in calc_ref_mask, possibly because self.ref is nonexistent')


def add_shaded_relief(self, ax, data_grid, extent, alpha=0.3):
	"""
	Add shaded relief visualization to a reference surface plot using multidirectional hillshade technique.
	
	Args:
		ax: matplotlib axis to add shaded relief to
		data_grid: 2D numpy array of elevation/depth data
		extent: extent tuple (left, right, bottom, top)
		alpha: transparency of the shaded relief overlay (default 0.3)
	"""
	if data_grid is None or np.all(np.isnan(data_grid)):
		return
		
	# Prepare data for hillshade calculation
	masked_data = np.ma.array(data_grid, mask=np.isnan(data_grid))
	data_for_hillshade = np.nan_to_num(masked_data.filled(np.nanmin(data_grid)))
	
	# Use multidirectional hillshade: average of 4 azimuths
	azimuths = [45, 135, 225, 315]
	altitude = 45
	hillshades = []
	
	for az in azimuths:
		ls = LightSource(azdeg=az, altdeg=altitude)
		hillshade = ls.hillshade(data_for_hillshade, vert_exag=0.6)
		hillshades.append(hillshade)
	
	# Average the hillshades
	hillshade_array = np.mean(hillshades, axis=0)
	
	# Add the hillshade as a gray overlay
	ax.imshow(hillshade_array, extent=extent, cmap='gray', origin='upper', 
			  alpha=alpha, zorder=1)


def plot_ref_surf(self):
	# plot reference depth, density, slope, and final masked grids
	print('in plot_ref_surf, calling calc_ref_mask')
	calc_ref_mask(self)  # update masks before plotting
	print('back in plot_ref_surf after calling calc_ref_mask')

	if 'final_mask' not in self.ref.keys():
		update_log(self, 'Final mask not computed; ref surf will not be plotted')
		return

	ones_mask = np.ones_like(self.ref['final_mask'])  # alternative mask to show all data in each grid

	# update subplots with reference surface and masks; all subplots use same extent from depth grid
	if 'z_grid' in self.ref:  # plot depth and final depths if available
		print('plotting depth)')
		self.clim_z = [np.nanmin(self.ref['z_grid']), np.nanmax(self.ref['z_grid'])]
		self.cbar_dict['z']['clim'] = self.clim_z
		self.cbar_dict['z_filt']['clim'] = self.clim_z
		self.cbar_dict['z_final']['clim'] = self.clim_z
		self.cbar_dict['z_depth']['clim'] = self.clim_z

		# subplot depth as parsed or as filtered
		plot_mask = (self.ref['z_mask'] if self.update_ref_plots_chk.isChecked() else ones_mask)
		self.surf_ax1.imshow(self.ref['z_grid']*plot_mask, interpolation='none', cmap='rainbow',
							 vmin=self.clim_z[0], vmax=self.clim_z[1], extent=self.ref['z_ref_extent'])
		
		# Add shaded relief if enabled
		if hasattr(self, 'show_shaded_relief_chk') and self.show_shaded_relief_chk.isChecked():
			add_shaded_relief(self, self.surf_ax1, self.ref['z_grid']*plot_mask, self.ref['z_ref_extent'])

						# large plot of final masked depth surface
		print('*******plotting final masked depth surface on final tab')
		print('the final masked depth grid will be =', self.ref['z_grid']*self.ref['final_mask'])
		self.surf_ax5.imshow(self.ref['z_grid']*self.ref['final_mask'], interpolation='none', cmap='rainbow',
					vmin=self.clim_z[0], vmax=self.clim_z[1], extent=self.ref['z_ref_extent'])
		
		# Add shaded relief if enabled
		if hasattr(self, 'show_shaded_relief_chk') and self.show_shaded_relief_chk.isChecked():
			add_shaded_relief(self, self.surf_ax5, self.ref['z_grid']*self.ref['final_mask'], self.ref['z_ref_extent'])
		
		# Add crossline coverage to final surface tab if enabled
		print(f"DEBUG: Checking crossline coverage condition: 'e' in self.xline = {'e' in self.xline}, 'n' in self.xline = {'n' in self.xline}, show_xline_cov_chk.isChecked() = {self.show_xline_cov_chk.isChecked()}")
		if 'e' in self.xline and 'n' in self.xline and self.show_xline_cov_chk.isChecked():
			print(f"DEBUG: Crossline coverage condition met! Proceeding to plot crossline coverage.")
			# Apply filter to crossline data
			filter_idx = np.where(self.xline['filter_idx'])[0]
			print(f"DEBUG: filter_idx length = {len(filter_idx)}")
			real_e = [self.xline['e'][i] for i in filter_idx]
			real_n = [self.xline['n'][i] for i in filter_idx]
			print(f"DEBUG: real_e length = {len(real_e)}, real_n length = {len(real_n)}")
			
			# Decimate data for plotting
			# Pass beam angles as first element for bin-based decimation
			beam_angles = [self.xline['beam_angle'][i] for i in filter_idx]
			dec_data = decimate_data(self, data_list=[beam_angles, real_e, real_n])
			real_e_dec = dec_data[1]  # Second element is easting
			real_n_dec = dec_data[2]  # Third element is northing
			# Defensive: Only print debug info if decimated data are list/array/np.ndarray and not int
			if (isinstance(real_e_dec, (list, tuple, np.ndarray)) and not isinstance(real_e_dec, int) and isinstance(real_n_dec, (list, tuple, np.ndarray)) and not isinstance(real_n_dec, int)):
				len_e_dec = len(real_e_dec)
				len_n_dec = len(real_n_dec)
				print(f"DEBUG: After decimation: real_e_dec length = {len_e_dec}, real_n_dec length = {len_n_dec}")
				print(f"DEBUG: Sample real_e_dec values: {real_e_dec[:5] if len_e_dec > 0 else 'No data'}")
				print(f"DEBUG: Sample real_n_dec values: {real_n_dec[:5] if len_n_dec > 0 else 'No data'}")
				print(f"DEBUG: real_e_dec range: {min(real_e_dec) if len_e_dec > 0 else 'No data'} to {max(real_e_dec) if len_e_dec > 0 else 'No data'}")
				print(f"DEBUG: real_n_dec range: {min(real_n_dec) if len_n_dec > 0 else 'No data'} to {max(real_n_dec) if len_n_dec > 0 else 'No data'}")
			else:
				print(f"DEBUG: real_e_dec is not a sequence, type: {type(real_e_dec)}, value: {real_e_dec}")
				print(f"DEBUG: real_n_dec is not a sequence, type: {type(real_n_dec)}, value: {real_n_dec}")
			
			# Plot crossline soundings on final surface tab
			print(f"DEBUG: About to plot crossline coverage on final surface tab with {len(real_e_dec)} points")
			self.surf_ax5.scatter(real_e_dec, real_n_dec,
								 s=self.pt_size_cov, c='lightgray', marker='o', alpha=self.pt_alpha_cov, linewidths=0)
			print(f"DEBUG: Finished plotting crossline coverage on final surface tab")
			# Plot tracklines on final surface tab
			for f in self.xline_track.keys():
				self.surf_ax5.scatter(self.xline_track[f]['e'], self.xline_track[f]['n'],
									 s=2, c='black', marker='o', linewidths=2)
			# Force redraw of the canvas to ensure points are visible
			print(f"DEBUG: Checking if surf_final_canvas exists: {hasattr(self, 'surf_final_canvas')}")
			if hasattr(self, 'surf_final_canvas'):
				print(f"DEBUG: Forcing redraw of surf_final_canvas")
				self.surf_final_canvas.draw()
			else:
				print(f"DEBUG: surf_final_canvas does not exist")
		
		# Plot depth surface on depth tab
		self.depth_ax.clear()
		self.depth_ax.imshow(self.ref['z_grid'] * self.ref['z_mask'], interpolation='none', cmap='rainbow',
					vmin=self.clim_z[0], vmax=self.clim_z[1], extent=self.ref['z_ref_extent'])
		
		# Add shaded relief to depth tab if enabled
		if hasattr(self, 'show_shaded_relief_chk') and self.show_shaded_relief_chk.isChecked():
			add_shaded_relief(self, self.depth_ax, self.ref['z_grid'] * self.ref['z_mask'], self.ref['z_ref_extent'])
		

		
		# Add labels to depth tab
		self.depth_ax.set_xlabel('Easting (m, UTM ' + self.ref_utm_str + ')', fontsize=8)
		self.depth_ax.set_ylabel('Northing (m, UTM ' + self.ref_utm_str + ')', fontsize=8)
		self.depth_ax.use_sticky_edges = False
		self.depth_ax.margins(float(self.axis_margin/100))
		self.depth_ax.autoscale(True)
		self.depth_ax.tick_params(axis='x', labelsize=6)
		self.depth_ax.tick_params(axis='y', labelsize=6)
		self.depth_ax.set_title('Reference Surface (Depth)', fontsize=10)
		
		# Apply gridlines to depth tab
		if self.grid_lines_toggle_chk.isChecked():
			self.depth_ax.grid()
			self.depth_ax.minorticks_on()
			self.depth_ax.grid(which='minor', linestyle='-', linewidth='0.5', color='black')
			self.depth_ax.grid(which='major', linestyle='-', linewidth='1.0', color='black')
		else:
			self.depth_ax.grid(False)
			self.depth_ax.minorticks_off()

		# Add depth tooltip functionality
		self._setup_depth_tooltip()

		print('finished plotting final masked ref surf')
		print(f"DEBUG: self.xline keys = {list(self.xline.keys()) if hasattr(self, 'xline') and self.xline else 'No xline data'}")

	# plot density if available
	if 'c_grid' in self.ref and 'z_ref_extent' in self.ref:
		self.clim_c = [np.nanmin(self.ref['c_grid']), np.nanmax(self.ref['c_grid'])]
		self.cbar_dict['c']['clim'] = self.clim_c
		self.cbar_dict['c_final']['clim'] = self.clim_c
		# plot density grid as parsed or as filtered
		plot_mask = (self.ref['c_mask'] if self.update_ref_plots_chk.isChecked() else ones_mask)
		self.surf_ax2.imshow(self.ref['c_grid']*plot_mask, interpolation='none', cmap='rainbow',
							 vmin=self.clim_c[0], vmax=self.clim_c[1], extent=self.ref['z_ref_extent'])
		
		# Plot density final masked surface on density final tab
		self.density_final_ax.clear()
		self.density_final_ax.imshow(self.ref['c_grid'] * self.ref['c_mask'], interpolation='none', cmap='rainbow',
									vmin=self.clim_c[0], vmax=self.clim_c[1], extent=self.ref['z_ref_extent'])
		
		
		
		# Add labels to density final tab
		self.density_final_ax.set_xlabel('Easting (m, UTM ' + self.ref_utm_str + ')', fontsize=8)
		self.density_final_ax.set_ylabel('Northing (m, UTM ' + self.ref_utm_str + ')', fontsize=8)
		self.density_final_ax.use_sticky_edges = False
		self.density_final_ax.margins(float(self.axis_margin/100))
		self.density_final_ax.autoscale(True)
		self.density_final_ax.tick_params(axis='x', labelsize=6)
		self.density_final_ax.tick_params(axis='y', labelsize=6)
		self.density_final_ax.set_title('Reference Surface (Density Final)', fontsize=10)
		
		# Apply gridlines to density final tab
		if self.grid_lines_toggle_chk.isChecked():
			self.density_final_ax.grid()
			self.density_final_ax.minorticks_on()
			self.density_final_ax.grid(which='minor', linestyle='-', linewidth='0.5', color='black')
			self.density_final_ax.grid(which='major', linestyle='-', linewidth='1.0', color='black')
		else:
			self.density_final_ax.grid(False)
			self.density_final_ax.minorticks_off()

		# Add density tooltip functionality
		self._setup_density_tooltip()

	else:
		# self.surf_ax2.clear()
		update_log(self, 'No sounding density data available for plotting/filtering (load .xyd text file of density '
						 'corresponding to reference depth .xyz file)')

	# plot slope if available
	if 's_grid' in self.ref and 'z_ref_extent' in self.ref:
		# Use actual data range instead of hardcoded [0, 5]
		slope_data = self.ref['s_grid']
		if self.update_ref_plots_chk.isChecked():
			slope_data = slope_data * self.ref['s_mask']
		# Calculate actual min/max from valid (non-NaN) data
		valid_slope = slope_data[~np.isnan(slope_data)]
		if len(valid_slope) > 0:
			self.clim_s = [np.nanmin(valid_slope), np.nanmax(valid_slope)]
		else:
			self.clim_s = [0, 1]  # fallback if no valid data
		self.cbar_dict['s']['clim'] = self.clim_s
		self.cbar_dict['s_final']['clim'] = self.clim_s
		# plot max slope as calculated or as filtered
		plot_mask = (self.ref['s_mask'] if self.update_ref_plots_chk.isChecked() else ones_mask)
		self.surf_ax3.imshow(self.ref['s_grid']*plot_mask, interpolation='none', cmap='rainbow',
							 vmin=self.clim_s[0], vmax=self.clim_s[1], extent=self.ref['z_ref_extent'])
		
		# Plot slope final masked surface on slope final tab
		self.slope_final_ax.clear()
		self.slope_final_ax.imshow(self.ref['s_grid'] * self.ref['s_mask'], interpolation='none', cmap='rainbow',
								  vmin=self.clim_s[0], vmax=self.clim_s[1], extent=self.ref['z_ref_extent'])
		

		
		# Add labels to slope final tab
		self.slope_final_ax.set_xlabel('Easting (m, UTM ' + self.ref_utm_str + ')', fontsize=8)
		self.slope_final_ax.set_ylabel('Northing (m, UTM ' + self.ref_utm_str + ')', fontsize=8)
		self.slope_final_ax.use_sticky_edges = False
		self.slope_final_ax.margins(float(self.axis_margin/100))
		self.slope_final_ax.autoscale(True)
		self.slope_final_ax.tick_params(axis='x', labelsize=6)
		self.slope_final_ax.tick_params(axis='y', labelsize=6)
		self.slope_final_ax.set_title('Reference Surface (Slope Final)', fontsize=10)
		
		# Apply gridlines to slope final tab
		if self.grid_lines_toggle_chk.isChecked():
			self.slope_final_ax.grid()
			self.slope_final_ax.minorticks_on()
			self.slope_final_ax.grid(which='minor', linestyle='-', linewidth='0.5', color='black')
			self.slope_final_ax.grid(which='major', linestyle='-', linewidth='1.0', color='black')
		else:
			self.slope_final_ax.grid(False)
			self.slope_final_ax.minorticks_off()

		# Add slope tooltip functionality
		self._setup_slope_tooltip()

		print('survived plotting slope')

	# plot uncertainty if available and selected by user (replaces small "final" masked surface subplot)
	if 'u_grid' in self.ref and 'z_ref_extent' in self.ref and self.show_u_plot_chk.isChecked():
		print('plotting uncertainty')
		# print('self.ref[u_grid] =', self.ref['u_grid'])
		# print('self.ref[z_ref_extent] =', self.ref['z_ref_extent'])
		print('in u_grid, setting clim_u')
		self.clim_u = [np.nanmin(self.ref['u_grid']), np.nanmax(self.ref['u_grid'])]
		print('in plot_ref_surf, plotting uncertainty with self.clim_u = ', self.clim_u)
		self.cbar_dict['u']['clim'] = self.clim_u
		self.cbar_dict['u_final']['clim'] = self.clim_u
		print('assigned self.cbar_dict =', self.cbar_dict)
		print('calling plot_mask')
		plot_mask = (self.ref['u_mask'] if self.update_ref_plots_chk.isChecked() else ones_mask)
		print('survived plot_mask')
		self.surf_ax4.imshow(self.ref['u_grid'] * plot_mask, interpolation='none', cmap='rainbow',
							 vmin=self.clim_u[0], vmax=self.clim_u[1], extent=self.ref['z_ref_extent'])
		
		# Plot uncertainty final masked surface on uncertainty final tab
		self.uncertainty_final_ax.clear()
		self.uncertainty_final_ax.imshow(self.ref['u_grid'] * self.ref['u_mask'], interpolation='none', cmap='rainbow',
										vmin=self.clim_u[0], vmax=self.clim_u[1], extent=self.ref['z_ref_extent'])
		

		
		# Add labels to uncertainty final tab
		self.uncertainty_final_ax.set_xlabel('Easting (m, UTM ' + self.ref_utm_str + ')', fontsize=8)
		self.uncertainty_final_ax.set_ylabel('Northing (m, UTM ' + self.ref_utm_str + ')', fontsize=8)
		self.uncertainty_final_ax.use_sticky_edges = False
		self.uncertainty_final_ax.margins(float(self.axis_margin/100))
		self.uncertainty_final_ax.autoscale(True)
		self.uncertainty_final_ax.tick_params(axis='x', labelsize=6)
		self.uncertainty_final_ax.tick_params(axis='y', labelsize=6)
		self.uncertainty_final_ax.set_title('Reference Surface (Uncertainty Final)', fontsize=10)
		
		# Apply gridlines to uncertainty final tab
		if self.grid_lines_toggle_chk.isChecked():
			self.uncertainty_final_ax.grid()
			self.uncertainty_final_ax.minorticks_on()
			self.uncertainty_final_ax.grid(which='minor', linestyle='-', linewidth='0.5', color='black')
			self.uncertainty_final_ax.grid(which='major', linestyle='-', linewidth='1.0', color='black')
		else:
			self.uncertainty_final_ax.grid(False)
			self.uncertainty_final_ax.minorticks_off()

		# Add uncertainty tooltip functionality
		self._setup_uncertainty_tooltip()
		
	elif 'z_grid' in self.ref:  # otherwise, plot masked final depths
		self.surf_ax4.imshow(self.ref['z_grid'] * self.ref['final_mask'], interpolation='none', cmap='rainbow',
							 vmin=self.clim_z[0], vmax=self.clim_z[1], extent=self.ref['z_ref_extent'])
		
		# Add shaded relief if enabled
		if hasattr(self, 'show_shaded_relief_chk') and self.show_shaded_relief_chk.isChecked():
			add_shaded_relief(self, self.surf_ax4, self.ref['z_grid'] * self.ref['final_mask'], self.ref['z_ref_extent'])

	# add labels to all subplots (update uncertainty title later, if plotted)
	for ax, t in {self.surf_ax1: 'Reference Surface (Depth)', self.surf_ax2: 'Reference Surface (Density)',
				  self.surf_ax3: 'Reference Surface (Slope)', self.surf_ax4: 'Reference Surface (Final)',
				  self.surf_ax5: 'Reference Surface (Final)'}.items():
		ax.set_xlabel('Easting (m, UTM ' + self.ref_utm_str + ')', fontsize=8)
		ax.set_ylabel('Northing (m, UTM ' + self.ref_utm_str + ')', fontsize=8)
		ax.use_sticky_edges = False
		# ax.margins(float(self.axis_margin_tb.text())/100)
		ax.margins(float(self.axis_margin/100))
		ax.autoscale(True)
		# ticks = ax.xaxis.get_major_ticks()

		ax.tick_params(axis='x', labelsize=6)
		ax.tick_params(axis='y', labelsize=6)
		# for tick_ax in [ax.xaxis, ax.yaxis]:
		# 	print('working on tick_ax =', tick_ax)
		# 	ticks = tick_ax.get_major_ticks()
		# 	print('got ticks =', ticks)
		# 	for tick in ticks:
		# 		print('setting fontsize for tick =', tick)
		# 		# tick.label.set_fontsize(6) ##### the label attribute no longer exists in newer Matplotlib

		ax.set_title(t, fontsize=10)

		if ax == self.surf_ax4 and self.show_u_plot_chk.isChecked():  # update subplot4 uncertainty title
			ax.set_title('Reference Surface (Uncertainty)', fontsize=10)

	# sort out colorbars
	for subplot, params in self.cbar_dict.items():  # set colorbars for each ref surf subplot
		if params['cax'] != None:  # remove all colorbars
			params['cax'].remove()

		if subplot == ['u', 'z_filt'][int(self.show_u_plot_chk.isChecked())]:  # skip u or z_filt plot in subplot4
			params['cax'] = None  # save placeholder cax = None for this ax so removal step doesn't fail on next refresh
			continue

		# Skip creating colorbars for tabs that don't need them
		# Create colorbars for both surface filter plots and final tabs
		if subplot not in ['z', 'c', 's', 'u', 'z_filt', 'z_final', 'z_depth', 'c_final', 's_final', 'u_final']:
			continue

		clim = params['clim']
		cbaxes = inset_axes(params['ax'], width="2%", height="30%", loc=self.cbar_loc)
		
		# Use 5 tick marks for reference surface plots (z, c, s, u), 11 for others
		if subplot in ['z', 'c', 's', 'u', 'z_filt', 'z_final']:
			tval = np.linspace(clim[0], clim[1], 5)
		else:
			tval = np.linspace(clim[0], clim[1], 11)
			
		cbar = colorbar.ColorbarBase(cbaxes, cmap='rainbow', orientation='vertical',
									norm=colors.Normalize(clim[0], clim[1]),
									ticks=tval, ticklocation='left')

		cbar.ax.tick_params(labelsize=self.cbar_font_size)  # set font size for entries
		cbar.set_label(label=params['label'], size=self.cbar_title_font_size)
		# cbar.ax.set_facecolor('white')
		# cbar.set_facecolor('white')
		# cbar.patch.set_color('white')
		# cbar.patch.set_facecolor((1.0,1.0,1.0,1.0))
		# cbar.patch.set_facecolor('white')
		# cbar.patch.set_fill(True)
		# cbar.filled(True)
		# cbar.set_color('white')

		if subplot in ['c', 'c_final']:
			tlab = ['%d' % tick for tick in tval]  # integer sounding count tick labels

		elif subplot in ['s', 's_final']:
			tlab = ['%0.1f' % float(tick) for tick in tval]  # slope tick labels

		else:
			tlab = ['%0.1f' % float((1 if subplot in ['u'] else -1)*tick) for tick in tval]

		cbar.set_ticklabels(tlab)
		params['cax'] = cbar

	# plot crossline soundings and trackline if checked
	if 'e' in self.xline and 'n' in self.xline and self.show_xline_cov_chk.isChecked():
		# nan_idx = np.isnan(self.xline['e'])
		# print('got nan_idx with sum of non-nans=', np.sum(~nan_idx))
		# real_e = np.asarray(self.xline['e'])[~nan_idx].tolist()
		# real_n = np.asarray(self.xline['n'])[~nan_idx].tolist()
		# print('got real_e with len=', len(real_e))
		# print('got real_n with len=', len(real_n))

		# original method with all coverage
		# dec_data = decimate_data(self, data_list=[deepcopy(self.xline['e']), deepcopy(self.xline['n'])])
		# real_e_dec = dec_data[0]
		# real_n_dec = dec_data[1]

		# new method with coverage based on filters
		filter_idx = np.where(self.xline['filter_idx'])[0]
		real_e = (np.asarray(self.xline['e'])[filter_idx]).tolist()
		real_n = (np.asarray(self.xline['n'])[filter_idx]).tolist()

		dec_data = decimate_data(self, data_list=[real_e, real_n])
		real_e_dec = dec_data[0]
		real_n_dec = dec_data[1]

		self.surf_ax5.scatter(real_e_dec, real_n_dec,
							  s=self.pt_size_cov, c='lightgray', marker='o', alpha=self.pt_alpha_cov, linewidths=0)
		
		for f in self.xline_track.keys():  # plot soundings on large final surface plot
			self.surf_ax5.scatter(self.xline_track[f]['e'], self.xline_track[f]['n'],
								  s=2, c='black', marker='o', linewidths=2)

def plot_tide(self, set_active_tab=False):
	# plot imported tide data
	if not all([k in self.tide.keys() for k in ['time_obj', 'amplitude']]):
		print('Tide time and amplitude not available for plotting')
		return

	update_log(self, 'Plotting tide')
	plot_tide_optimized(self, set_active_tab)


def plot_tide_optimized(self, set_active_tab=False):
	"""Optimized tide plotting - shows crossline time ranges as colored regions over tide curve"""
	from time import process_time
	from datetime import timedelta
	start_time = process_time()
	
	self.tide_ax.clear()  # Clear the axes before plotting
	
	# Find crossline time ranges first to determine plot limits
	plot_start = None
	plot_end = None
	
	if all([k in self.xline.keys() for k in ['datetime', 'fname']]):
		# Group by filename to get time ranges for each crossline file
		file_times = {}
		for i, fname in enumerate(self.xline['fname']):
			if fname not in file_times:
				file_times[fname] = []
			file_times[fname].append(self.xline['datetime'][i])
		
		# Find the overall time range of crossline data
		all_crossline_times = []
		for times in file_times.values():
			all_crossline_times.extend(times)
		
		if all_crossline_times:
			crossline_start = min(all_crossline_times)
			crossline_end = max(all_crossline_times)
			
			# Get tide range from UI parameter
			try:
				tide_range_hours = float(self.tide_range_tb.text())
			except (ValueError, AttributeError):
				tide_range_hours = 12.0  # Default fallback
			
			# Set plot limits to ±tide_range_hours around crossline time range
			plot_start = crossline_start - timedelta(hours=tide_range_hours)
			plot_end = crossline_end + timedelta(hours=tide_range_hours)
			
			update_log(self, f'Crossline time range: {crossline_start} to {crossline_end}')
			update_log(self, f'Plot time range: {plot_start} to {plot_end} (±{tide_range_hours} hours)')
	
	# Use the complete original tide data for plotting
	tide_times = self.tide['time_obj']
	tide_amplitudes = self.tide['amplitude']
	
	# Plot the complete tide curve
	self.tide_ax.plot(tide_times, tide_amplitudes,
					  color='black', marker='o', markersize=self.pt_size/10, alpha=self.pt_alpha,
					  label='Tide Curve')
	
	# Plot colored regions for each crossline file
	if all([k in self.xline.keys() for k in ['datetime', 'fname']]):
		colors = ['red', 'orange', 'green', 'blue', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
		color_idx = 0
		
		for fname, times in file_times.items():
			if times:
				start_time_file = min(times)
				end_time_file = max(times)
				
				# Create shaded region
				color = colors[color_idx % len(colors)]
				self.tide_ax.axvspan(start_time_file, end_time_file, 
									alpha=0.3, color=color, 
									label=f'{fname} ({len(times)} pings)')
				
				color_idx += 1
		
		# Add legend if we have multiple files
		if len(file_times) > 1:
			self.tide_ax.legend(loc='upper left', fontsize=8, ncol=1)
	
	# Set title and labels
	if plot_start and plot_end:
		self.tide_ax.set_title(f'Tide Curve with Crossline Time Ranges (±{tide_range_hours}h)')
	else:
		self.tide_ax.set_title('Tide Curve with Crossline Time Ranges')
	self.tide_ax.set_xlabel('Time')
	self.tide_ax.set_ylabel('Tide Amplitude (m)')
	
	# Format x-axis to show dates properly
	self.tide_ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
	self.tide_ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))  # More frequent ticks for focused view
	self.tide_figure.autofmt_xdate()
	
	# Set x-axis limits if we have a specific range
	if plot_start and plot_end:
		self.tide_ax.set_xlim(plot_start, plot_end)
	
	# Force autoscale for y-axis only
	self.tide_ax.relim()
	self.tide_ax.autoscale_view(scalex=False, scaley=True)  # Don't autoscale x-axis if we set limits
	
	# Always draw the canvas
	self.tide_canvas.draw()
	
	# Force a repaint
	self.tide_canvas.repaint()
	
	# Force the canvas to update
	self.tide_canvas.update()
	
	self.tide_canvas.show()
	
	plot_time = process_time() - start_time
	update_log(self, f'Tide plotting completed in {plot_time:.2f} seconds')
	
	if set_active_tab:
		self.plot_tabs.setCurrentIndex(3)  # Set to tide tab


def plot_tide2(self, set_active_tab=False):
    # plot imported tide data - ALL DATA since sorting takes longer than its worth
    if not all([k in self.tide.keys() for k in ['time_obj', 'amplitude']]):
        print('Tide time and amplitude not available for plotting')
        return

    update_log(self, 'Plotting tide')
    
    self.tide_ax.clear()  # Clear the axes before plotting
    
    self.tide_ax.plot(self.tide['time_obj'], self.tide['amplitude'],
                      color='black', marker='o', markersize=self.pt_size/10, alpha=self.pt_alpha)
    
    # Add the missing functionality from original plot_tide() - plot applied tide data
    if all([k in self.xline.keys() for k in ['tide_applied', 'datetime']]):
        # get unique ping times by finding where applied tide diff != 0, rather than resorting
        testing = set(self.xline['datetime'])
        ping_idx = [self.xline['datetime'].index(t) for t in set(self.xline['datetime'])]  # get unique ping times
        ping_time_set = [self.xline['datetime'][i] for i in ping_idx]
        tide_ping_set = [self.xline['tide_applied'][i] for i in ping_idx]
        sort_idx = np.argsort(ping_time_set)
        self.tide_ax.plot(np.asarray(ping_time_set)[sort_idx], np.asarray(tide_ping_set)[sort_idx], 'ro',
                          markersize=self.pt_size / 10)
    
    # Set title and labels
    self.tide_ax.set_title('Tide Applied to Accuracy Crosslines')
    self.tide_ax.set_xlabel('Time')
    self.tide_ax.set_ylabel('Tide Amplitude (m)')
    
    # Format x-axis to show dates properly
    self.tide_ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
    self.tide_ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
    self.tide_figure.autofmt_xdate()
    
    # Force autoscale
    self.tide_ax.relim()
    self.tide_ax.autoscale_view()
    
    # Always draw the canvas
    self.tide_canvas.draw()
    
    # Force a repaint
    self.tide_canvas.repaint()
    
    # Force the canvas to update
    self.tide_canvas.update()
    
    self.tide_canvas.show()
    
    if set_active_tab:
        self.plot_tabs.setCurrentIndex(3)  # Set to tide tab





def parse_crosslines(self):
	# parse crosslines
	from time import process_time
	start_time = process_time()
	
	update_log(self, '=== PARSING CROSSLINE FILES ===')
	try:
		fnames_xline = list(set(self.xline['fname']))  # make list of unique filenames already in detection dict
		update_log(self, f'Found {len(fnames_xline)} previously parsed files')

	except:
		fnames_xline = []  # self.xline has not been created yet; initialize this and self.xline detection dict
		self.xline = {}
		self.xline_track = {}
		update_log(self, 'Initializing new crossline data structures')

	# fnames_new_all = self.get_new_file_list('.all', fnames_xline)  # list new .all files not included in det dict
	# fnames_new = get_new_file_list(self, ['.all', '.kmall'], fnames_xline)  # list all files not in xline dict
	fnames_new = get_new_file_list(self, ['.all', '.kmall', 'ASCII.txt'], fnames_xline)  # list all files not in xline dict

	num_new_files = len(fnames_new)
	# update_log(self, 'Found ' + str(len(fnames_new)) + ' new crossline .all files')

	if num_new_files == 0:
		update_log(self, 'No new crossline files to parse')
		return 0

	else:
		update_log(self, f'Found {len(fnames_new)} new crossline files to parse')
		update_log(self, f'Calculating accuracy from {num_new_files} new file(s)')
		QtWidgets.QApplication.processEvents()  # try processing and redrawing the GUI to make progress bar update
		data_new = {}
		track_new ={}

		# update progress bar and log
		self.calc_pb.setValue(0)  # reset progress bar to 0 and max to number of files
		self.calc_pb.setMaximum(len(fnames_new))

		for f in range(len(fnames_new)):
			# for f in range(len(fnames_new_all)):         # read previously unparsed files
			fname_str = fnames_new[f].rsplit('/')[-1]
			self.current_file_lbl.setText(
				'Parsing new file [' + str(f + 1) + '/' + str(num_new_files) + ']:' + fname_str)
			QtWidgets.QApplication.processEvents()
			ftype = fname_str.rsplit('.', 1)[-1]
			
			file_start = process_time()
			update_log(self, f'Parsing file {f+1}/{num_new_files}: {fname_str} ({ftype.upper()})')

			if ftype == 'all':
				# parse IPSTART73, RRA78, POS80, RTP82, XYZ88
				# data_new[f] = multibeam_tools.libs.readEM.parseEMfile(fnames_new[f],
				#                                                       parse_list=[73, 78, 80, 82, 88],
				#                                                       print_updates=False)

				# WORKING METHOD
				# data_new[f] = readALLswath(self, fnames_new[f], print_updates=False, parse_outermost_only=False)

				# TESTING METHOD
				print('sending .all file to readALLswath')
				data = readALLswath(self, fnames_new[f], print_updates=False, parse_outermost_only=False)
				# print('got data back from readAllswath with type =', type(data))
				# print('now sending dictionary = {0:data} to convertXYZ')
				converted_data = convertXYZ({0: data}, print_updates=False)  # convertXYZ for dict of parsed .all data
				# converted_data = convertXYZ({0: data}, print_updates=True)  # convertXYZ for dict of parsed .all data
				# print('got back converted_data with type =', type(converted_data))
				# print('now trying to store converted data in data_new[f]')
				data_new[f] = converted_data[0]
				# print('stored converted data, data_new is now', data_new)

				# store xline track data separately from detection data for plotting
				# WARNING: THIS IS OVERWRITING EXISTING TRACKS (e.g., f=0, 1, 2, etc.) WHEN NEW FILES ARE ADDED
				# self.xline_track[f] = {k: data_new[f][k] for k in ['POS', 'IP']}  # store POS and IP for track
				# self.xline_track[f]['fname'] = fname_str
				track_new[f] = {k: data_new[f][k] for k in ['POS', 'IP']}  # store POS and IP for track
				track_new[f]['fname'] = fname_str
				# self.xline_track[fname_str] = {k: data_new[f][k] for k in ['POS', 'IP']}  # store POS and IP for track
				# print('stored track=', self.xline_track[fname_str])

			elif ftype == 'kmall':
				data_new[f] = readKMALLswath(self, fnames_new[f], read_mode='plot')

				# store RTP with pingInfo lat/lon
				print('storing pingInfo lat/lon as ship track for this kmall file')
				# self.xline_track[f] = {k: data_new[f][k] for k in ['RTP', 'IP']}
				# self.xline_track[f]['fname'] = fname_str
				track_new[f] = {k: data_new[f][k] for k in ['HDR', 'RTP', 'IP']}
				track_new[f]['fname'] = fname_str
				# print('data_new[IP]=', data_new[f]['IP'])
				# print('IP text =', data_new[f]['IP']['install_txt'])

			elif ftype == 'txt':
				print('********** GOING TO TRY PARSING ASCII SOUNDINGS')
				zone = self.ref_proj_cbox.currentText()
				data_new[f] = readASCIIswath(self, fnames_new[f], utm_zone=zone)  # read soundings from ASCII
				track_new[f] = track_new[f] = {k: data_new[f][k] for k in ['HDR', 'RTP', 'IP']}  # read track from ASCII (vessel X, vessel Y)
				track_new[f]['fname'] = fname_str
				track_new[f]['e'] = data_new[f]['ping_e']
				track_new[f]['n'] = data_new[f]['ping_n']
				track_new[f]['datetime'] = data_new[f]['datetime']

			else:
				update_log(self, 'Warning: Skipping unrecognized file type for ' + fname_str)

			file_time = process_time() - file_start
			update_log(self, f'Parsed file {fname_str} in {file_time:.2f} seconds')
			update_prog(self, f + 1)

		parse_time = process_time() - start_time
		update_log(self, f'File parsing completed in {parse_time:.2f} seconds')
		print('finished parsing, data_new has keys=', data_new.keys(), ' and data_new[0].keys = ', data_new[0].keys())

		# convert XYZ to lat, lon using active pos sensor; maintain depth as reported in file; interpret/verify modes
		update_log(self, 'Interpreting modes and converting coordinates...')
		interpret_start = process_time()
		# self.data_new = convertXYZ(data_new, print_updates=True)  # for .all files only?
		self.data_new = interpretMode(self, data_new, print_updates=self.print_updates)  # True)
		interpret_time = process_time() - interpret_start
		update_log(self, f'Mode interpretation completed in {interpret_time:.2f} seconds')
		print('survived interpretMode, self.data_new has keys')
		
		# print('testing calling sortDetectionsAccuracy before verifying modes')
		update_log(self, 'Sorting detections for accuracy analysis...')
		sort_start = process_time()
		det_new = sortDetectionsAccuracy(self, self.data_new, print_updates=True)  # sort new accuracy soundings
		sort_time = process_time() - sort_start
		update_log(self, f'Detection sorting completed in {sort_time:.2f} seconds')
		print('survived sortDetectionsAccuracy, det_new has keys', det_new.keys())

		# files_OK, EM_params = verifyMode(self.data_new)  # check install and runtime params
		# if not files_OK:  # warn user if inconsistencies detected (perhaps add logic later for sorting into user-selectable lists for archiving and plotting)
		# 	update_log(self, 'WARNING! CROSSLINES HAVE INCONSISTENT MODEL, S/N, or RUNTIME PARAMETERS')


		# sort ship track
		track_new = sort_xline_track(self, track_new)  # sort crossline track after adding any new files  ---> UPDATE TO KEEP EARLIER TRACKS
		print('got track_new with keys = ', track_new.keys())

		if len(self.xline) == 0 and len(self.xline_track) == 0:  # if detection dict is empty, store all new detections
			print('len of self.xline and self.xline_track == 0, so setting equal to det_new and track_new')

			self.xline = det_new
			self.xline_track = track_new

		# print('det_new =', det_new)
			# print('track_new = ', track_new)

		else:  # otherwise, append new detections to existing detection dict
			for key, value in det_new.items():  # loop through the new data and append to existing self.xline
				print('extending self.xline with key =', key)
				# print('value has type', type(value), 'and = ', value)
				self.xline[key].extend(value)

			# for key, value in track_new[0].items():  # loop through the new track and append to existing self.xline_track
			for key, value in track_new.items():  # loop through the new track fnames append to existing self.xline_track

				print('extending self.xline_track with key =', key)
				# print('value has type', type(value), ' and = ', value)
				self.xline_track[key] = value
				# print('extending self.xline_track with key =', key)

		# sort_xline_track(self)  # sort crossline track after adding any new files  ---> UPDATE TO KEEP EARLIER TRACKS
		# track_new = sort_xline_track(self, track_new)  # sort crossline track after adding any new files  ---> UPDATE TO KEEP EARLIER TRACKS

		# verify consistent system info and runtime params
		sys_info = verifyModelAndModes(self.xline)
		if len(sys_info['datetime']) == 1:
			update_log(self, 'Consistent model, serial number, and modes in all files')

		else:
			log_str = ''
			for i in range(len(sys_info['datetime'])):
				# temp_str = '\n' + ','.join([k + ': ' + str(sys_info[k][i]) for k in sys_info.keys()])
				temp_str = '\n' + ', '.join([str(sys_info[k][i]) for k in \
											 ['fname', 'model', 'sn', 'ping_mode', 'swath_mode', 'pulse_form']])
				log_str += temp_str

			update_log(self, 'Found inconsistent parameters:\n' + \
							 'Filename: Model, S/N, Ping Mode, Swath Mode, Pulse Form\n' + log_str)


		total_time = process_time() - start_time
		update_log(self, f'=== CROSSLINE PARSING COMPLETED in {total_time:.2f} seconds ===')
		update_log(self, f'Total soundings loaded: {len(self.xline.get("z", []))}')
		update_log(self, f'Finished parsing {num_new_files} new file(s)')
		self.current_file_lbl.setText('Current File [' + str(f + 1) + '/' + str(num_new_files) +
									  ']: Finished parsing crosslines')

	self.calc_accuracy_btn.setStyleSheet("color: black; font-weight: normal")  # reset the button text color to default

	return num_new_files


def sortDetectionsAccuracy(self, data, print_updates=False):
	# sort through .all and .kmall data dict and store valid soundings, BS, and modes
	# note: .all data must be converted from along/across/depth data to lat/lon with convertXYZ before sorting
	det_key_list = ['fname',  'model', 'datetime', 'date', 'time', 'sn',
					'lat', 'lon', 'x', 'y', 'z', 'z_re_wl', 'n', 'e', 'utm_zone', 'bs', 'rx_angle',
					'ping_mode', 'pulse_form', 'swath_mode', 'frequency',
					'max_port_deg', 'max_stbd_deg', 'max_port_m', 'max_stbd_m',
					'tx_x_m', 'tx_y_m', 'tx_z_m', 'aps_x_m', 'aps_y_m', 'aps_z_m', 'wl_z_m',
					'ping_e', 'ping_n', 'ping_utm_zone']  # mode_bin

	det = {k: [] for k in det_key_list}

	# examine detection info across swath, find outermost valid soundings for each ping
	for f in range(len(data)):  # loop through all data
		if print_updates:
			print('Sorting detections for accuracy, f =', f, ' and data[f] keys =', data[f].keys())
		# set up keys for dict fields of interest from parsers for each file type (.all or .kmall)
		ftype = data[f]['fname'].rsplit('.', 1)[1]
		key_idx = int(ftype == 'kmall')  # keys in data dicts depend on parser used, get index to select keys below
		det_int_threshold = [127, 0][key_idx]  # threshold for valid sounding (.all  <128 and .kmall == 0)
		det_int_key = ['RX_DET_INFO', 'detectionType'][key_idx]  # key for detect info depends on ftype
		depth_key = ['RX_DEPTH', 'z_reRefPoint_m'][key_idx]  # key for depth
		across_key = ['RX_ACROSS', 'y_reRefPoint_m'][key_idx]  # key for acrosstrack distance
		along_key = ['RX_ALONG', 'z_reRefPoint_m'][key_idx]  # key for alongtrack distance
		bs_key = ['RX_BS', 'reflectivity1_dB'][key_idx]  # key for backscatter in dB
		bs_scale = [0.1, 1][key_idx]  # backscatter scale in X dB; multiply parsed value by this factor for dB
		angle_key = ['RX_ANGLE', 'beamAngleReRx_deg'][key_idx]  # key for RX angle re RX array
		lat_key = ['SOUNDING_LAT', 'lat'][key_idx]
		lon_key = ['SOUNDING_LON', 'lon'][key_idx]
		e_key = ['SOUNDING_E', 'e'][key_idx]
		n_key = ['SOUNDING_N', 'n'][key_idx]
		utm_key = ['SOUNDING_UTM_ZONE', 'utm_zone'][key_idx]

		for p in range(len(data[f]['XYZ'])):  # loop through each ping
			# print('working on ping number ', p)
			det_int = data[f]['XYZ'][p][det_int_key]  # get detection integers for this ping
			det_idx = [i for i, v in enumerate(det_int) if v <= det_int_threshold]  # indices of all valid detections

			# extend swath data from appropriate keys/values in data dicts
			# future general sorter: accuracy, keep all valid det_int; coverage, reduce for outermost valid det_int
			det['fname'].extend([data[f]['fname'].rsplit('/')[-1]] * len(det_idx))  # store fname for each det
			det['x'].extend([data[f]['XYZ'][p][along_key][i] for i in det_idx])  # as parsed
			det['y'].extend([data[f]['XYZ'][p][across_key][i] for i in det_idx])  # as parsed
			det['z'].extend([data[f]['XYZ'][p][depth_key][i] for i in det_idx])  # as parsed

			det['lat'].extend([data[f]['XYZ'][p][lat_key][i] for i in det_idx])
			det['lon'].extend([data[f]['XYZ'][p][lon_key][i] for i in det_idx])

			det['n'].extend([data[f]['XYZ'][p][n_key][i] for i in det_idx])
			det['e'].extend([data[f]['XYZ'][p][e_key][i] for i in det_idx])
			# det['bs'].extend([data[f]['XYZ'][p][bs_key][i] for i in det_idx])
			det['bs'].extend([data[f]['XYZ'][p][bs_key][i] * bs_scale for i in det_idx])
			det['ping_mode'].extend([data[f]['XYZ'][p]['PING_MODE']] * len(det_idx))
			det['pulse_form'].extend([data[f]['XYZ'][p]['PULSE_FORM']] * len(det_idx))
			# det['z_re_wl'].extend([data[f]['XYZ'][p]['SOUNDING_Z'][i] for i in det_idx])  # corrected to waterline
			# det['ping_utm_zone'].extend([data[f]['XYZ'][p]['PING_UTM_ZONE']] * len(det_idx))
			# det['ping_e'].extend([data[f]['XYZ'][p]['PING_E']] * len(det_idx))
			# det['ping_n'].extend([data[f]['XYZ'][p]['PING_N']] * len(det_idx))

			if ftype == 'all':  # .all store date and time from ms from midnight
				det['model'].extend([data[f]['XYZ'][p]['MODEL']] * len(det_idx))
				det['sn'].extend([data[f]['XYZ'][p]['SYS_SN']] * len(det_idx))
				dt = datetime.strptime(str(data[f]['XYZ'][p]['DATE']), '%Y%m%d') + \
					 timedelta(milliseconds=data[f]['XYZ'][p]['TIME'])
				det['datetime'].extend([dt] * len(det_idx))
				det['date'].extend([dt.strftime('%Y-%m-%d')] * len(det_idx))
				det['time'].extend([dt.strftime('%H:%M:%S.%f')] * len(det_idx))
				det['utm_zone'].extend([data[f]['XYZ'][p][utm_key]] * len(det_idx))  # convertXYZ --> one utmzone / ping
				# det['rx_angle'].extend([data[f]['RRA_78'][p][angle_key][i] for i in det_idx])
				det['swath_mode'].extend([data[f]['XYZ'][p]['SWATH_MODE']] * len(det_idx))
				det['max_port_deg'].extend([data[f]['XYZ'][p]['MAX_PORT_DEG']] * len(det_idx))
				det['max_stbd_deg'].extend([data[f]['XYZ'][p]['MAX_STBD_DEG']] * len(det_idx))
				det['max_port_m'].extend([data[f]['XYZ'][p]['MAX_PORT_M']] * len(det_idx))
				det['max_stbd_m'].extend([data[f]['XYZ'][p]['MAX_STBD_M']] * len(det_idx))
				# print('in ping', p, 'with data[f][IP_START] =', data[f]['IP_start'])
				det['tx_x_m'].extend([data[f]['XYZ'][p]['TX_X_M']] * len(det_idx))
				det['tx_y_m'].extend([data[f]['XYZ'][p]['TX_Y_M']] * len(det_idx))
				det['tx_z_m'].extend([data[f]['XYZ'][p]['TX_Z_M']] * len(det_idx))
				det['aps_x_m'].extend([data[f]['XYZ'][p]['APS_X_M']] * len(det_idx))
				det['aps_y_m'].extend([data[f]['XYZ'][p]['APS_Y_M']] * len(det_idx))
				det['aps_z_m'].extend([data[f]['XYZ'][p]['APS_Z_M']] * len(det_idx))
				det['wl_z_m'].extend([data[f]['XYZ'][0]['WL_Z_M']] * len(det_idx))

			elif ftype == 'kmall':  # .kmall store date and time from datetime object
				det['model'].extend([data[f]['HDR'][p]['echoSounderID']] * len(det_idx))
				det['datetime'].extend([data[f]['HDR'][p]['dgdatetime']] * len(det_idx))
				det['date'].extend([data[f]['HDR'][p]['dgdatetime'].strftime('%Y-%m-%d')] * len(det_idx))
				det['time'].extend([data[f]['HDR'][p]['dgdatetime'].strftime('%H:%M:%S.%f')] * len(det_idx))
				det['utm_zone'].extend([data[f]['XYZ'][p][utm_key][i] for i in det_idx])  # readKMALLswath 1 utm/sounding
				det['aps_x_m'].extend([0] * len(det_idx))  # not needed for KMALL; append 0 as placeholder
				det['aps_y_m'].extend([0] * len(det_idx))  # not needed for KMALL; append 0 as placeholder
				det['aps_z_m'].extend([0] * len(det_idx))  # not needed for KMALL; append 0 as placeholder

				# get first installation parameter datagram in file for s/n and offsets, assume no changes within file
				ip_text = data[f]['IP']['install_txt'][0]
				# get TX array offset text: EM304 = 'TRAI_TX1' and 'TRAI_RX1', EM2040P = 'TRAI_HD1', not '_TX1' / '_RX1'
				ip_tx1 = ip_text.split('TRAI_')[1].split(',')[0].strip()  # all heads/arrays split by comma; use 1st hd
				det['tx_x_m'].extend([float(ip_tx1.split('X=')[1].split(';')[0].strip())] * len(det_idx))  # TX array X
				det['tx_y_m'].extend([float(ip_tx1.split('Y=')[1].split(';')[0].strip())] * len(det_idx))  # TX array Y
				det['tx_z_m'].extend([float(ip_tx1.split('Z=')[1].split(';')[0].strip())] * len(det_idx))  # TX array Z
				det['wl_z_m'].extend([float(ip_text.split('SWLZ=')[-1].split(',')[0].strip())] * len(det_idx))  # WL Z

				# get serial number from installation parameter: 'SN=12345'
				sn = ip_text.split('SN=')[1].split(',')[0].strip()
				det['sn'].extend([sn] * len(det_idx))

				# get index of latest runtime parameter timestamp prior to ping of interest; default to 0 for cases
				# where earliest pings in file might be timestamped earlier than first runtime parameter datagram
				IOP_headers = data[f]['IOP']['header']  # get list of IOP header dicts in new kmall module output
				IOP_datetimes = [IOP_headers[d]['dgdatetime'] for d in range(len(IOP_headers))]
				# print('got IOP datetimes =', IOP_datetimes)

				MRZ_headers = data[f]['HDR']
				MRZ_datetimes = [MRZ_headers[d]['dgdatetime'] for d in range(len(MRZ_headers))]

				# find index of last IOP datagram before current ping, default to first if
				IOP_idx = max([i for i, t in enumerate(IOP_datetimes) if t <= MRZ_datetimes[p]], default=0)

				if IOP_datetimes[IOP_idx] > MRZ_datetimes[p]:
					print('*****ping', p, 'occurred before first runtime datagram; using first RTP dg in file')

				# get runtime text from applicable IOP datagram, split and strip at keywords and append values
				rt = data[f]['IOP']['runtime_txt'][IOP_idx]  # get runtime text for splitting NEW KMALL FORMAT
				# print('IOP_idx = ', IOP_idx)
				print('rt = ', rt)

				# dict of keys for detection dict and substring to split runtime text at entry of interest
				rt_dict = {'max_port_deg': 'Max angle Port:', 'max_stbd_deg': 'Max angle Starboard:',
						   'max_port_m': 'Max coverage Port:', 'max_stbd_m': 'Max coverage Starboard:'}

				# iterate through rt_dict and append coverage limits from split/stripped runtime text
				for k, v in rt_dict.items():
					try:
						det[k].extend([float(rt.split(v)[-1].split('\n')[0].strip())] * len(det_idx))

					except:
						det[k].extend(['NA'] * len(det_idx))

				# parse swath mode text
				try:
					dual_swath_mode = rt.split('Dual swath:')[-1].split('\n')[0].strip()

					# print('kmall dual_swath_mode =', dual_swath_mode)

					depth_mode = rt.split('Depth setting:')[-1].split('\n')[0].strip().lower()
					print('found depth_mode =', depth_mode)
					print('depth_mode in [very deep, extra deep, extreme deep] = ',
						  depth_mode in ['very deep', 'extra deep', 'extreme deep'])

					# if dual_swath_mode == 'Off':
					if dual_swath_mode == 'Off' or depth_mode in ['very deep', 'extra deep', 'extreme deep']:
						swath_mode = 'Single Swath'

					else:
						swath_mode = 'Dual Swath (' + dual_swath_mode + ')'

				except:
					swath_mode = 'NA'

				det['swath_mode'].extend([swath_mode] * len(det_idx))

				# parse frequency from runtime parameter text, if available
				try:
					# print('trying to split runtime text')
					frequency_rt = rt.split('Frequency:')[-1].split('\n')[0].strip().replace('kHz', ' kHz')
					# print('frequency string from runtime text =', frequency_rt)

				except:  # use default frequency stored from interpretMode
					pass

				# store parsed freq if not empty, otherwise store default
				frequency = frequency_rt if frequency_rt else data[f]['XYZ'][p]['FREQUENCY']
				det['frequency'].extend([frequency] * len(det_idx))

				if print_updates:
					# print('found IOP_idx=', IOP_idx, 'with IOP_datetime=', data[f]['IOP']['dgdatetime'][IOP_idx])
					print('found IOP_idx=', IOP_idx, 'with IOP_datetime=', IOP_datetimes[IOP_idx])
					print('max_port_deg=', det['max_port_deg'][-1])
					print('max_stbd_deg=', det['max_stbd_deg'][-1])
					print('max_port_m=', det['max_port_m'][-1])
					print('max_stbd_m=', det['max_stbd_m'][-1])
					print('swath_mode=', det['swath_mode'][-1])

			elif ftype == 'txt':  #ASCII text
				print('\n\n**** NEED TO SORT ASCII TEXT SOUNDINGS ****\n\n')

			else:
				print('UNSUPPORTED FTYPE --> NOT SORTING DETECTION!')

	if print_updates:
		print('\nDone sorting detections...')

	# print('leaving sortDetectionsCoverage with det[frequency] =', det['frequency'])

	return det


def calc_z_final(self):
	# adjust sounding depths to desired reference and flip sign as necessary for comparison to ref surf (positive up)
	_, _, dz_ping = adjust_depth_ref(self.xline, depth_ref=self.ref_cbox.currentText().lower())
	_, _, dz_ping_sonar = adjust_depth_ref(self.xline, depth_ref='TX array'.lower())
	# print('got dz_ping =', dz_ping)
	# print('got dz_ping_sonar =', dz_ping_sonar)
	# print('dz_ping has len', len(dz_ping))
	# print('first 20 of xline[z]=', self.xline['z'][0:20])
	# print('first 20 of dz_ping =', dz_ping[0:20])
	# z_final = [z + dz for z, dz in zip(self.xline['z'], dz_ping)]  # add dz

	# adjust reported depths (+DOWN) for tide at ping time (+UP); e.g., if the waterline-adjusted sounding was reported
	# as +10 m (down) when the tide was +1 m (up), the sounding is +9 m (down) from the tide datum
	# print('the first few dz_ping values are', dz_ping[0:10])
	# interpolate tide onto ping time
	# print('in calc_z_final, self.xlines.keys =', self.xline.keys())
	# print('the first few dates are ', self.xline['date'][0:10])
	# print('the first few times are ', self.xline['time'][0:10])
	# ping_times = [datetime.datetime.strptime(d + ' ' + t, '%Y-%m-%d %H:%M:%S.%f') for d, t in
	# 			  zip(self.xline['date'], self.xline['time'])]
	ping_times = self.xline['datetime']
	# print('the first few ping times are now', ping_times[0:10])

	self.tide_applied = False
	# ping_times = self.xlines
	if all([k in self.tide.keys() for k in ['time_obj', 'amplitude']]):
		print('working on tide interpolation onto ping times')
		# interp step needs non-datetime-object axis; use UNIX, assume UTC; sub-second precision not necessary!
		epoch = datetime.utcfromtimestamp(0)
		print('got epoch = ', epoch)
		tide_times_s = [(dt - epoch).total_seconds() for dt in self.tide['time_obj']]  # tide time in s from 1970
		
		# Fix timezone issue: make ping_times timezone-naive if they are timezone-aware
		ping_times_naive = []
		for dt in ping_times:
			if dt.tzinfo is not None:
				# Convert timezone-aware datetime to naive UTC
				ping_times_naive.append(dt.replace(tzinfo=None))
			else:
				ping_times_naive.append(dt)
		
		ping_times_s = [(dt - epoch).total_seconds() for dt in ping_times_naive]

		if ping_times_s[0] < tide_times_s[0] or ping_times_s[-1] > tide_times_s[-1]:
			update_log(self, 'WARNING: ping times found outside the tide record; zero tide will be applied',
					   font_color="red")
			tide_ping = np.zeros_like(self.xline['z'])

		else:
			update_log(self, 'Interpolating tide onto final crossline sounding times')
			tide_ping = np.interp(ping_times_s, tide_times_s, self.tide['amplitude'])
			self.tide_applied = True

		# print('got first few tide_time_s =', tide_times_s[0:10])
		# print('got first few ping_time_s =', ping_times_s[0:10])
		self.xline['tide_applied'] = deepcopy(tide_ping.tolist())
		# self.xline['time_obj'] = deepcopy(ping_times)
		# print('just stored self.xline[tide_applied] --> first few =', self.xline['tide_applied'][0:10])

	else:
		update_log(self, 'WARNING: No tide applied during final Z calculation', font_color="red")
		tide_ping = np.zeros_like(self.xline['z'])

	# print('got tide_ping first few values:', tide_ping[0:10])

	# add dz to bring soundings to waterline and then subtract tide to bring soundings to tide datum; results are +DOWN
	print('first couple Z values prior to adjustment: ', self.xline['z'][0:10])
	z_final = [z + dz - tide for z, dz, tide in zip(self.xline['z'], dz_ping, tide_ping)]  # z from user ref, w/ tide
	z_sonar = [z + dz for z, dz in zip(self.xline['z'], dz_ping_sonar)]  # z obs. from sonar ref for angle calc

	# add waterline adjustment (m +DOWN); a pos. value shifts WL down/lower in the mapping frame and results in
	# shallower final xline depths for processing
	wl_adjustment = [float(self.waterline_tb.text())]*len(z_final)  # get current waterline adjustment (m +down)
	z_final_wl = [z - dwl for z, dwl in zip(z_final, wl_adjustment)]  # z +DOWN after WL adjustment (also z+ DOWN)
	print('first ten wl_adjustment, z_final, and z_final_wl = ', wl_adjustment[:10], z_final[:10], z_final_wl[:10])

	# report final z re waterline and z re sonar
	self.xline['z_final'] = (-1 * np.asarray(z_final_wl)).tolist()  # flip sign to neg down and store 'final' soundings
	self.xline['z_sonar'] = (-1 * np.asarray(z_sonar)).tolist()  # store z from TX array for nominal beam angle calc
	z_sonar_pos_down = (-1 * np.asarray(self.xline['z_sonar'])).tolist()
	self.xline['beam_angle'] = np.rad2deg(np.arctan2(self.xline['y'], z_sonar_pos_down)).tolist()

	print('first couple Z_final values after adjustment: ', self.xline['z_final'][0:10])
	print('first couple Z_sonar values after adjustment: ', self.xline['z_sonar'][0:10])


def convert_crossline_utm(self):
	# if necessary, convert crossline X, Y to UTM zone of reference surface
	update_log(self, 'Checking UTM zones of ref grid and crossline(s)')
	ref_utm = self.ref['utm_zone']
	# format xline UTM zone for comparison with ref_utm and use with pyproj; replace zone letter with S if southern
	# hemisphere (UTM zone letter C-M) or N if northern hemisphere (else)
	xline_utm = [utm_str.replace(" ", "") for utm_str in self.xline['utm_zone']]
	xline_utm = [utm_str[:-1] + 'S' if utm_str[-1] <= 'M' else utm_str[:-1] + 'N' for utm_str in xline_utm]
	self.xline['utm_zone'] = xline_utm  # replace with new format
	print('detected ref surf UTM =', ref_utm, ' and set of xline utm zones =', set(xline_utm))
	xline_utm_list = [u for u in set(xline_utm) if u != ref_utm]  # get list of xline utm zones != ref surf utm zone
	print('non-matching utm zones:', xline_utm_list)

	if len(xline_utm_list) > 0:  # transform soundings from non-matching xline utm zone(s) into ref utm zone
		update_log(self, 'Found crossline soundings in UTM zone (' + ', '.join(xline_utm_list) + \
				   ') other than selected ref. surface UTM zone (' + ref_utm +'); transforming soundings to ' + ref_utm)
		utm_zone = int(''.join([c for c in self.ref['utm_zone'] if c.isnumeric()]))  # need integer for pyproj.Proj
		south_hem = True if self.ref['utm_zone'][-1].lower() == 's' else False  # need hemisphere for pyproj.Proj
		print('in convert_crossline_utm, using zone = ', utm_zone, ' and south_hem =', south_hem, ' for pyproj.Proj')

		# define projection of reference surface and numpy array for easier indexing
		# print('in convert_crossline_utm, sending utm zone to pyproj = ', utm_zone)
		# p2 = pyproj.Proj(proj='utm', zone=ref_utm, ellps='WGS84')
		p2 = pyproj.Proj(proj='utm', zone=utm_zone, ellps='WGS84', south=south_hem)  # update UTM zone with hemisphere

		xline_e = np.asarray(self.xline['e'])
		xline_n = np.asarray(self.xline['n'])
		N_soundings = len(self.xline['utm_zone'])
		print('N_soundings is originally', N_soundings)

		for u in xline_utm_list:  # for each non-matching xline utm zone, convert those soundings to ref utm
			print('working on non-matching utm zone', u)
			p1_zone = int(''.join([c for c in u if c.isnumeric()]))  # need integer for pyproj.Proj
			p1_south = True if u[-1].lower() == 's' else False  # need hemisphere for pyproj.Proj
			print('in convert_crossline_utm, working on non-matching zone =', p1_zone, ' and south =', p1_south)
			# p1 = pyproj.Proj(proj='utm', zone=u, ellps='WGS84')  # define proj of xline soundings
			p1 = pyproj.Proj(proj='utm', zone=p1_zone, ellps='WGS84', south=p1_south)  # define proj of xline with hem

			print('first ten xline_utm are:', xline_utm[0:10])

			idx = [s for s in range(N_soundings) if xline_utm[s] == u]  # get indices of soundings with this zone
			# print('length of idx is', str(len(idx)))
			print('first ten xline_utm for idx matches are:', [xline_utm[i] for i in idx[0:10]])

			(xline_e_new, xline_n_new) = pyproj.transform(p1, p2, xline_e[idx], xline_n[idx])  # transform
			xline_e[idx] = xline_e_new
			xline_n[idx] = xline_n_new

			update_log(self, 'Transformed ' + str(len(idx)) + ' soundings (out of '
					   + str(N_soundings) + ') from ' + u + ' to ' + ref_utm)

			print('fixed the eastings to:', xline_e_new)

		# reassign the final coordinates
		self.xline['e'] = xline_e.tolist()
		self.xline['n'] = xline_n.tolist()
		self.xline['utm_zone'] = [ref_utm] * N_soundings

		print('new xline_easting is', self.xline['e'][0:30])
		print('new xline_northing is', self.xline['n'][0:30])
		print('new utm_zone is', self.xline['utm_zone'][0:30])


def convert_track_utm(self):
	# if necessary, convert crossline track X,Y to UTM zone of reference surface
	update_log(self, 'Checking UTM zones of ref grid and crossline track(s)')
	ref_utm = self.ref['utm_zone']
	print('in convert_track_utm, ref_utm=', ref_utm)
	# get set of
	track_utm_set = [u for u in set([self.xline_track[f]['utm_zone'] for f in self.xline_track.keys()]) if u != ref_utm]
	print('track_utm_set =', track_utm_set)

	if len(track_utm_set) == 0:  # return if no conversions are needed
		return

	# update user and continue with conversions as necessary
	update_log(self, 'Found tracklines in UTM zone(s) (' + ', '.join([u for u in track_utm_set]) + ') ' + \
			   'other than selected ref. surface UTM zone (' + ref_utm + '); transforming track to ' + ref_utm)

	# # format xline UTM zone for comparison with ref_utm and use with pyproj; replace zone letter with S if southern
	# # hemisphere (UTM zone letter C-M) or N if northern hemisphere (else)
	# xline_utm = [utm_str.replace(" ", "") for utm_str in self.xline['utm_zone']]
	# xline_utm = [utm_str[:-1] + 'S' if utm_str[-1] <= 'M' else utm_str[:-1] + 'N' for utm_str in xline_utm]

	p2_zone = int(''.join([c for c in self.ref['utm_zone'] if c.isnumeric()]))  # need integer for pyproj.Proj
	p2_south = True if self.ref['utm_zone'][-1].lower() == 's' else False  # need hemisphere for pyproj.Proj

	print('in convert_track_utm, sending the output utm zone to pyproj: ', p2_zone, ' and south hem = ', p2_south)

	# p2 = pyproj.Proj(proj='utm', zone=ref_utm, ellps='WGS84')  # define ref surf projection for desired transform output
	p2 = pyproj.Proj(proj='utm', zone=p2_zone, ellps='WGS84', south=p2_south)  # define ref surf proj for output

	for f in self.xline_track.keys():  # check track utm zone for each fname (key) and transform to ref UTM zone if nec.
		print('in file f=', f, 'the xline_track utm zone is', self.xline_track[f]['utm_zone'])
		track_utm = self.xline_track[f]['utm_zone']

		if track_utm != ref_utm:  # one utm zone assigned to each track dict (key = fname)
			track_e = np.asarray(self.xline_track[f]['e'])
			track_n = np.asarray(self.xline_track[f]['n'])
			print('in convert_track_utm loop, sending the track_utm = ', track_utm)
			p1_zone = int(''.join([c for c in track_utm if c.isnumeric()]))  # need integer for pyproj.Proj
			p1_south = True if self.xline_track[f]['utm_zone'][-1].lower() == 's' else False
			# p1 = pyproj.Proj(proj='utm', zone=track_utm, ellps='WGS84')  # define proj of current track line
			p1 = pyproj.Proj(proj='utm', zone=p1_zone, ellps='WGS84', south=p1_south)  # define proj of current track line
			(track_e_new, track_n_new) = pyproj.transform(p1, p2, track_e, track_n)  # transform all track points
			update_log(self, 'Transformed ' + str(len(track_e_new)) + ' track points from ' + \
					   track_utm + ' to ' + ref_utm)

			# reassign the final coordinates
			self.xline_track[f]['e'] = track_e_new.tolist()
			self.xline_track[f]['n'] = track_n_new.tolist()
			self.xline_track[f]['utm_zone'] = ref_utm  # store single UTM zone for whole file

			print('new track easting is', self.xline_track[f]['e'][0:30])
			print('new track northing is', self.xline_track[f]['n'][0:30])
			print('new track utm_zone is', self.xline_track[f]['utm_zone'])


def calc_dz_from_ref_interp(self):
	# calculate the difference of each sounding from the reference grid (interpolated onto sounding X, Y position)
	update_log(self, 'Calculating ref grid depths at crossline sounding positions')
	print('N ref_surf nodes e =', len(self.ref['e']), 'with first ten =', self.ref['e'][0:10])
	print('N ref_surf nodes n =', len(self.ref['n']), 'with first ten =', self.ref['n'][0:10])
	print('N ref_surf nodes z =', len(self.ref['z']), 'with first ten =', self.ref['z'][0:10])
	print('N xline soundings e =', len(self.xline['e']), 'with first ten =', self.xline['e'][0:10])
	print('N xline soundings n =', len(self.xline['n']), 'with first ten =', self.xline['n'][0:10])
	print('N xline soundings z =', len(self.xline['z']), 'with first ten =', self.xline['z'][0:10])
	print('N xline soundings z_final =', len(self.xline['z_final']), 'with first ten =', self.xline['z_final'][0:10])

	print('\n\n ******** MASKING REF GRID PRIOR TO DZ CALC *************')

	print('')

	# get all nodes in masked final reference grid in shape, set nans to inf for interpolating xline soundings
	print('starting calc_dz_from_ref_interp with self.ref[e_grid]=', self.ref['e_grid'])
	print('starting calc_dz_from_ref_interp with self.ref[n_grid]=', self.ref['n_grid'])
	print('starting calc_dz_from_ref_interp with self.ref_cell_size =', self.ref_cell_size)

	e_ref = self.ref['e_grid'].flatten()  # available reference grid (includes nans)
	e_ref_range = np.arange(np.nanmin(e_ref), np.nanmax(e_ref) + self.ref_cell_size, self.ref_cell_size)
	n_ref = self.ref['n_grid'].flatten()
	n_ref_range = np.arange(np.nanmin(n_ref), np.nanmax(n_ref) + self.ref_cell_size, self.ref_cell_size)
	z_ref = np.multiply(self.ref['z_grid'], self.ref['final_mask']).flatten()  # mask final ref z_grid
	print('e_final, n_final, z_final have shape', np.shape(e_ref), np.shape(n_ref), np.shape(z_ref))

	# create fully populated easting and northing grids without NANs (ref Z NANs are OK) in order to use griddata with
	# nearest interp method; this is ultimately used to identify / mask soundings that occur over NaN ref grid cells
	e_ref_full = np.tile(e_ref_range, (np.size(n_ref_range), 1))  # full grid eastings
	n_ref_full = np.flipud(np.transpose(np.tile(n_ref_range, (np.size(e_ref_range), 1))))  # full grid northings



	if self.xline:
		print('number of INFs in e_final, n_final, z_final, e_xline, n_xline =',
			  [np.sum(np.isinf(thing)) for thing in [e_ref, n_ref, z_ref, self.xline['e'], self.xline['n']]])
		print('number of NANs in e_final, n_final, z_final, e_xline, n_xline =',
			  [np.sum(np.isnan(thing)) for thing in [e_ref, n_ref, z_ref, self.xline['e'], self.xline['n']]])
	else:
		print('number of INFs in e_final, n_final, z_final =',
			  [np.sum(np.isinf(thing)) for thing in [e_ref, n_ref, z_ref]])
		print('number of NANs in e_final, n_final, z_final =',
			  [np.sum(np.isnan(thing)) for thing in [e_ref, n_ref, z_ref]])
		print('***Returning early from calc_dz_from_ref_interp***')
		return

	# linearly interpolate masked reference grid onto xline sounding positions, get mask with NaNs wherever the closest
	# reference grid node is a NaN, apply mask to interpolated xline ref depths such that all xline soundings 'off' the
	# masked ref grid will be nan and excluded from further analysis
	print('e_ref =', e_ref)
	print('n_ref =', n_ref)
	print('z_ref =', z_ref)

	z_ref_interp = griddata((self.ref['e'], self.ref['n']), self.ref['z'], (self.xline['e'], self.xline['n']), method='linear')

	print('made z_ref_interp')

	# linearly interpolate masked reference grid onto xline sounding positions, returns NANs for soundings over masked cells
	#### this behavior has changed since running on SEABREAM - exact code works and griddata accepts nans in the input,
	# but no longer works on POLLACK
	# z_ref_interp_mask = griddata((e_ref, n_ref), z_ref, (self.xline['e'], self.xline['n']), method='nearest')
	z_ref_interp_mask = griddata((e_ref_full.flatten(), n_ref_full.flatten()), z_ref,
								 (self.xline['e'], self.xline['n']), method='nearest')

	# goal is to identify xline soundings over empty/nan ref grid cells
	# z_ref_interp_mask = griddata((self.ref['e'], self.ref['n']), self.ref['z'], (self.xline['e'], self.xline['n']), method='nearest')
	# z_ref_interp_mask = griddata((e_ref, n_ref), z_ref, (self.xline['e'], self.xline['n']), method='linear')

	print('made z_ref_interp_mask')

	z_ref_interp_mask[~np.isnan(z_ref_interp_mask)] = 1.0  # set all non-nan values to 1 for masking
	self.xline['z_ref_interp'] = z_ref_interp*z_ref_interp_mask
	self.xline['z_ref_interp_mask'] = z_ref_interp_mask

	print('NUM NANS and INFS in z_ref_interp_mask =',
		  np.sum(np.isnan(self.xline['z_ref_interp'])), np.sum(np.isinf(self.xline['z_ref_interp'])))

	print('NUM XLINE ON FINAL REF SURF: ', np.sum(~np.isnan(self.xline['z_ref_interp'])), '/', np.size(self.xline['z']))
	# print('z_ref_interp looks like', self.xline['z_ref_interp'][0:30])
	# print('xline z_final after flip looks like', self.xline['z_final'][0:30])
	self.xline['num_on_ref'] = np.sum(~np.isnan(self.xline['z_ref_interp']))  # count non-Nan interp values
	update_log(self, 'Found ' + str(self.xline['num_on_ref']) +
			   ' crossline soundings on reference grid (after filtering, if applicable)')

	if self.xline['num_on_ref'] == 0:  # warn user if zero soundings found
		update_log(self, 'WARNING: Verify reference surface UTM zone and update if necessary (0 crossline soundings '
						 'found on reference grid in selected UTM zone ' + self.ref_utm_str + ')', font_color="red")

	# calculate dz for xline soundings with non-NaN interpolated reference grid depths
	# note that xline['z'] is positive down as returned from parser; flip sign for differencing from ref surf
	update_log(self, 'Calculating crossline sounding differences from filtered reference grid')
	self.xline['dz_ref'] = np.subtract(self.xline['z_final'], self.xline['z_ref_interp'])
	# print('xline dz_ref looks like', self.xline['dz_ref'][0:100])

	# store dz as percentage of water depth, with positive dz_ref_wd meaning shoal-bias crossline soundings to retain
	# intuitive plotting appearance, with shallower soundings above deeper soundings
	# e.g., if xline z = -98 and z_ref_interp = -100, then dz_ref = +2; dz_ref_wd should be positive; division of
	# positive bias (up) by reference depth (always negative) yields negative, so must be flipped in sign for plot
	dz_ref_wd = np.array(-1*100*np.divide(np.asarray(self.xline['dz_ref']), np.asarray(self.xline['z_ref_interp'])))
	self.xline['dz_ref_wd'] = dz_ref_wd.tolist()
	print('xline dz_ref_wd looks like', self.xline['dz_ref_wd'][0:100])
	self.ref['z_mean'] = np.nanmean(self.xline['z_ref_interp'])  # mean of ref grid interp values used


def bin_beamwise(self, refresh_plot=False):
	# bin by angle, calc mean and std of sounding differences in that angular bin
	from time import process_time
	start_time = process_time()
	
	print('starting bin_beamwise')
	update_log(self, '=== BEAM ANGLE BINNING ===')
	
	self.beam_bin_size = 1  # beam angle bin size (deg)
	self.beam_bin_lim = 75  # max angle (deg)
	update_log(self, f'Bin size: {self.beam_bin_size}°, Max angle: ±{self.beam_bin_lim}°')

	self.beam_bin_dz_mean = []  # declare dz mean, std, and sample count
	self.beam_bin_dz_std = []
	self.beam_bin_dz_N = []
	self.beam_bin_dz_wd_mean = []
	self.beam_bin_dz_wd_std = []
	self.beam_bin_dz_wd_zero = []
	self.dz_wd_bin_zero_mean = []  # bin mean value for soundings in this bin for plotting adjustment to zero mean
	self.beam_range = range(-1 * self.beam_bin_lim, self.beam_bin_lim, self.beam_bin_size)

	# if crossline data AND reference surface are available, convert soundings with meaningful reference surface
	# nodes to array for binning; otherwise, continue to refresh plot with empty results
	if self.xline == {}:  # skip if crossline data dict is empty (bin_beamwise was called only to reset stats)
		print('self.xline == {}; bin_beamwise called to reset stats')
		update_log(self, 'No crossline data available for binning')

	elif 'z_final' in self.xline and 'z' in self.ref:
		update_log(self, 'Applying crossline filters...')
		filter_start = process_time()
		filter_xline(self)
		filter_time = process_time() - filter_start
		update_log(self, f'Crossline filtering completed in {filter_time:.2f} seconds')
		
		# Count filtered soundings
		total_soundings = len(self.xline['z_final'])
		filtered_soundings = np.sum(self.xline['filter_idx'])
		update_log(self, f'Total soundings: {total_soundings}, Filtered soundings: {filtered_soundings} ({filtered_soundings/total_soundings*100:.1f}%)')
		
		update_log(self, 'Binning soundings by angle...')
		# calculate simplified beam angle from acrosstrack distance and depth
		# depth is used here as negative down re WL, consistent w/ %WD results
		# Kongsberg angle convention is right-hand-rule about +X axis (fwd), so port angles are + and stbd are -
		# however, for plotting purposes with conventional X axes, use neg beam angles to port and pos to stbd, per
		# plotting conventions used elsewhere (e.g., Qimera)
		# z_final_pos_down = (-1 * np.asarray(self.xline['z_final'])).tolist()
		# z_sonar_pos_down = (-1 * np.asarray(self.xline['z_sonar'])).tolist()

		# calculate nominal beam angle referenced from a horizontal plane at the arrays (sonar vertical height)
		# FUTURE: adjust y to sonar ref rather than raw ref (sonar in .all, mapping origin in .kmall) for more accurate
		# beam angle calcs; this method will suffice for .kmall as long as the athwartship offset from sonar to origin
		# is small compared to y and z
		# self.xline['beam_angle'] = np.rad2deg(np.arctan2(self.xline['y'], z_final_pos_down)).tolist()
		# self.xline['beam_angle'] = np.rad2deg(np.arctan2(self.xline['y'], z_sonar_pos_down)).tolist()

		# print('size of beam_angle is now', len(self.xline['beam_angle']))
		# print('first 30 beam angles are', self.xline['beam_angle'][0:30])
		# self.xline['beam_angle'] = (-1*np.rad2deg(np.arctan2(self.xline['y'], self.xline['z_re_wl']))).tolist()

		# print('found beam_angle in self.xline and z in self.ref')
		beam_array = np.asarray(self.xline['beam_angle'])
		dz_array = np.asarray(self.xline['dz_ref'])
		dz_wd_array = np.asarray(self.xline['dz_ref_wd'])
		dz_wd_zero_mean_array = np.zeros_like(dz_wd_array)  # adjustments for plotting zero mean

		filter_idx = self.xline['filter_idx']
		# print('filter_idx has size', np.size(filter_idx))

		# Check if filter_idx is empty (all crosslines removed)
		if len(filter_idx) == 0 or np.sum(filter_idx) == 0:
			print('filter_idx is empty or all False; no data to bin')
			# Initialize empty results
			for b in self.beam_range:
				self.beam_bin_dz_N.append(0)
				self.beam_bin_dz_mean.append(np.nan)
				self.beam_bin_dz_std.append(np.nan)
				self.beam_bin_dz_wd_mean.append(np.nan)
				self.beam_bin_dz_wd_std.append(np.nan)
				self.beam_bin_dz_wd_zero.append(np.nan)
			self.dz_wd_bin_zero_mean = []
			self.bin_count_mask = np.ones_like(beam_array, dtype=bool)  # All soundings valid when no data
			return

		min_bin_count = 0  # minimum number of soundings per bin; zero if not selected by user

		if self.bin_count_gb.isChecked():
			try:
				min_bin_count = float(self.min_bin_count_tb.text())

			except:
				min_bin_count = 0

		print('prior to binning, set min_bin_count = ', min_bin_count)
		print('DEBUG: bin_count_gb.isChecked() =', self.bin_count_gb.isChecked())
		print('DEBUG: min_bin_count_tb.text() =', self.min_bin_count_tb.text())

		for b in self.beam_range:  # loop through all beam bins, calc mean and std for dz results within each bin
			beam_idx = (beam_array >= b) & (beam_array < b + self.beam_bin_size)  # find indices of angles in this bin

			# print('for beam b, beam_idx has size', np.size(beam_idx))
			idx = np.logical_and(filter_idx, beam_idx)  # filtered soundings in this bin

			# print('Found', str(np.sum(beam_idx)), str(np.sum(idx)),
			# 	  '(unfiltered / filtered) soundings between', str(b), 'and', str(b + self.beam_bin_size), 'deg')
			# print('Found', str(np.sum(idx)), 'soundings between', str(b), 'and', str(b + self.beam_bin_size), 'deg')

			self.beam_bin_dz_N.append(np.sum(idx))

			# if np.sum(idx) > 0:  # calc only if at least one sounding on ref surf within this bin
			if np.sum(idx) > min_bin_count:  # calc only if N soundings on ref surf in this bin

				# THIS CAN BE UPDATED WITH A MINIMUM BIN COUNT SETTING
				# no filter_idx applied
				# self.beam_bin_dz_mean.append(np.nanmean(dz_array[idx]))
				# self.beam_bin_dz_std.append(np.nanstd(dz_array[idx]))
				# self.beam_bin_dz_wd_mean.append(np.nanmean(dz_wd_array[idx]))  # simple mean of WD percentages
				# self.beam_bin_dz_wd_std.append(np.nanstd(dz_wd_array[idx]))

				# apply filter idx
				self.beam_bin_dz_mean.append(np.nanmean(dz_array[idx]))
				self.beam_bin_dz_std.append(np.nanstd(dz_array[idx]))
				self.beam_bin_dz_wd_mean.append(np.nanmean(dz_wd_array[idx]))  # simple mean of WD percentages
				self.beam_bin_dz_wd_std.append(np.nanstd(dz_wd_array[idx]))
				self.beam_bin_dz_wd_zero.append(0)  # store zero for bins with N>0 soundings for plotting zero mean
				dz_wd_zero_mean_array[idx] = self.beam_bin_dz_wd_mean[-1]  # store current bin mean for these soundings

				# if b >= float(self.)

			else:  # store NaN if no dz results are available for this bin
				self.beam_bin_dz_mean.append(np.nan)
				self.beam_bin_dz_std.append(np.nan)
				self.beam_bin_dz_wd_mean.append(np.nan)  # this is the simple mean of WD percentages
				self.beam_bin_dz_wd_std.append(np.nan)
				self.beam_bin_dz_wd_zero.append(np.nan)  # store nan for bins with N=0 soundings for plotting zero mean

		self.dz_wd_bin_zero_mean = dz_wd_zero_mean_array.tolist()  # store adjustments to zero mean for plotting
		
		# Create a mask for individual soundings based on bin count filter
		# This will be used to filter out soundings that belong to bins with insufficient counts
		self.bin_count_mask = np.ones_like(beam_array, dtype=bool)  # Initialize all soundings as valid
		
		if min_bin_count > 0:
			for b in self.beam_range:
				beam_idx = (beam_array >= b) & (beam_array < b + self.beam_bin_size)
				idx = np.logical_and(filter_idx, beam_idx)
				
				# If this bin has insufficient soundings, mark all soundings in this bin as invalid
				if np.sum(idx) <= min_bin_count:
					self.bin_count_mask = np.logical_and(self.bin_count_mask, ~beam_idx)

	else:
		print('Error in bin_beamwise')


def plot_soundings_per_bin(self):
    """
    Plot a bar graph showing the number of soundings per bin.
    This helps visualize the bin count distribution and the effect of the bin count filter.
    """
    print("DEBUG: plot_soundings_per_bin called")
    
    if not hasattr(self, 'xline') or not self.xline:
        print("DEBUG: No xline data available")
        update_log(self, 'No crossline data loaded. Please load crossline files and calculate accuracy to show soundings per bin.')
        return
    
    if not hasattr(self, 'beam_bin_dz_N') or not self.beam_bin_dz_N:
        print("DEBUG: No binned data available")
        update_log(self, 'No binned data available. Please calculate accuracy to show soundings per bin.')
        return
    
    # Clear the figure
    self.soundings_figure.clear()
    
    # Create the main axis
    ax = self.soundings_figure.add_subplot(111)
    
    # Get the bin centers and counts
    beam_bin_centers = np.asarray([b + self.beam_bin_size / 2 for b in self.beam_range])
    bin_counts = np.asarray(self.beam_bin_dz_N)
    
    # Create the bar plot
    bars = ax.bar(beam_bin_centers, bin_counts, width=self.beam_bin_size * 0.8, 
                  color='skyblue', edgecolor='navy', alpha=0.7)
    
    # Add a horizontal line for the minimum bin count if the filter is enabled
    if hasattr(self, 'bin_count_gb') and self.bin_count_gb.isChecked():
        try:
            min_bin_count = float(self.min_bin_count_tb.text())
            ax.axhline(y=min_bin_count, color='red', linestyle='--', linewidth=2, 
                      label=f'Min bin count ({min_bin_count})')
            ax.legend()
        except (ValueError, AttributeError):
            pass
    
    # Set labels and title
    ax.set_xlabel('Beam Angle (deg)')
    ax.set_ylabel('Number of Soundings')
    ax.set_title('Soundings per Bin')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Set x-axis limits to match the beam range
    ax.set_xlim(self.beam_range[0], self.beam_range[-1] + self.beam_bin_size)
    
    # Add statistics text
    total_soundings = np.sum(bin_counts)
    valid_bins = np.sum(bin_counts > 0)
    max_soundings = np.max(bin_counts)
    mean_soundings = np.mean(bin_counts[bin_counts > 0]) if valid_bins > 0 else 0
    
    stats_text = f'Total soundings: {total_soundings:,}\n'
    stats_text += f'Valid bins: {valid_bins}\n'
    stats_text += f'Max per bin: {max_soundings}\n'
    stats_text += f'Mean per bin: {mean_soundings:.1f}'
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Refresh the canvas
    self.soundings_canvas.draw()
    print("DEBUG: plot_soundings_per_bin completed")


def plot_bathymetric_orders(self):
    """
    Plot bathymetric order accuracy limits on the accuracy plot.
    These are the IHO S-44 standards for hydrographic surveys.
    """
    print("DEBUG: plot_bathymetric_orders called")
    
    if not hasattr(self, 'xline'):
        print("DEBUG: No xline data available")
        update_log(self, 'No crossline data loaded. Please load crossline files and calculate accuracy to show bathymetric orders.')
        return
    
    if not self.xline:
        print("DEBUG: xline is empty")
        update_log(self, 'No crossline data loaded. Please load crossline files and calculate accuracy to show bathymetric orders.')
        return
        
    # Use filtered data for plotting
    if 'beam_angle' not in self.xline:
        print("DEBUG: beam_angle not in xline")
        update_log(self, 'Beam angle data not found. Please load crossline files and calculate accuracy to show bathymetric orders.')
        return
        
    if 'z_ref' in self.xline:
        ref_depths = np.abs(self.xline['z_ref'])
        print(f"DEBUG: Using z_ref, found {len(ref_depths)} depths")
    elif 'z_final' in self.xline:
        ref_depths = np.abs(self.xline['z_final'])
        print(f"DEBUG: Using z_final, found {len(ref_depths)} depths")
    else:
        print("DEBUG: Neither z_ref nor z_final found in xline")
        update_log(self, 'Depth data not found. Please load crossline files and calculate accuracy to show bathymetric orders.')
        return
        
    beam_angles = np.asarray(self.xline['beam_angle'])
    print(f"DEBUG: Found {len(beam_angles)} beam angles")
    
    # Use filter_idx to get valid points
    if 'filter_idx' in self.xline:
        filt_idx = self.xline['filter_idx']
        print(f"DEBUG: Using filter_idx, original data: {len(beam_angles)} points")
        beam_angles = beam_angles[filt_idx]
        ref_depths = ref_depths[filt_idx]
        print(f"DEBUG: After filtering: {len(beam_angles)} points")
    else:
        print("DEBUG: No filter_idx found, using all data")
        
    if len(ref_depths) == 0 or len(beam_angles) == 0:
        print("DEBUG: No valid data after filtering")
        update_log(self, 'No valid data after filtering. Please check your filter settings.')
        return

    orders = {
        'special_order': {'a': 0.25, 'b': 0.0075, 'color': 'red', 'linestyle': '--', 'label': 'Special Order'},
        'order_1a': {'a': 0.5, 'b': 0.013, 'color': 'orange', 'linestyle': '--', 'label': 'Order 1a'},
        'order_1b': {'a': 0.5, 'b': 0.013, 'color': 'yellow', 'linestyle': '--', 'label': 'Order 1b'},
        'order_2': {'a': 1.0, 'b': 0.023, 'color': 'green', 'linestyle': '--', 'label': 'Order 2'},
        'order_3': {'a': 2.0, 'b': 0.05, 'color': 'blue', 'linestyle': '--', 'label': 'Order 3'},
    }
    checkbox_map = {
        'special_order': 'show_special_order_chk',
        'order_1a': 'show_order_1a_chk',
        'order_1b': 'show_order_1b_chk',
        'order_2': 'show_order_2_chk',
        'order_3': 'show_order_3_chk'
    }
    # Plot on ax2 (bias vs angle)
    ax2 = getattr(self, 'ax2', None)
    if ax2 is None:
        print("DEBUG: ax2 not found")
        return
    print("DEBUG: Found ax2")
    
    # Check if any orders are enabled
    any_orders_enabled = False
    print("DEBUG: Checking for enabled orders:")
    for order_name, checkbox_attr in checkbox_map.items():
        if hasattr(self, checkbox_attr):
            checkbox = getattr(self, checkbox_attr)
            is_checked = checkbox.isChecked()
            print(f"DEBUG: {checkbox_attr}: {is_checked}")
            if is_checked:
                any_orders_enabled = True
        else:
            print(f"DEBUG: {checkbox_attr} not found")
    
    if not any_orders_enabled:
        print("DEBUG: No orders enabled")
        return
    print("DEBUG: Found enabled orders")
    
    # Calculate representative depth for the order lines
    # Use the median depth to avoid outliers affecting the lines
    median_depth = np.median(ref_depths)
    
    # Get the beam angle range for plotting
    beam_angle_range = np.linspace(np.min(beam_angles), np.max(beam_angles), 100)
    
    for order_name, checkbox_attr in checkbox_map.items():
        if hasattr(self, checkbox_attr) and getattr(self, checkbox_attr).isChecked():
            order_params = orders[order_name]
            
            # Calculate limits for the representative depth
            limit_meters = order_params['a'] + order_params['b'] * median_depth
            limit_percent = (limit_meters / median_depth) * 100
            
            # Plot positive and negative limits as horizontal lines
            ax2.axhline(y=limit_percent, color=order_params['color'], linestyle=order_params['linestyle'], 
                       label=order_params['label'] + ' (+)', linewidth=2, alpha=0.8)
            ax2.axhline(y=-limit_percent, color=order_params['color'], linestyle=order_params['linestyle'], 
                       label=order_params['label'] + ' (-)', linewidth=2, alpha=0.8)
    
    # Add legend if any lines were plotted
    handles, labels = ax2.get_legend_handles_labels()
    if handles:
        # Remove duplicate labels
        unique_labels = []
        unique_handles = []
        for handle, label in zip(handles, labels):
            if label not in unique_labels:
                unique_labels.append(label)
                unique_handles.append(handle)
        ax2.legend(unique_handles, unique_labels, loc='upper right', fontsize=8)
    
    # Add order compliance statistics to the plot
    _add_order_compliance_stats(self, ax2, ref_depths, beam_angles)


def _add_order_compliance_stats(self, ax, ref_depths, beam_angles):
    """
    Add compliance statistics text box to the plot showing how many points
    meet each bathymetric order standard.
    """
    if not hasattr(self, 'xline') or 'dz_ref_wd' not in self.xline:
        return
    
    # Get the accuracy data
    dz_wd = np.asarray(self.xline['dz_ref_wd'])
    if 'filter_idx' in self.xline:
        filt_idx = self.xline['filter_idx']
        dz_wd = dz_wd[filt_idx]
        ref_depths = ref_depths[filt_idx]
    
    if len(dz_wd) == 0:
        return
    
    # Calculate compliance for each order
    orders = {
        'Special Order': {'a': 0.25, 'b': 0.0075},
        'Order 1a': {'a': 0.5, 'b': 0.013},
        'Order 1b': {'a': 0.5, 'b': 0.013},
        'Order 2': {'a': 1.0, 'b': 0.023},
        'Order 3': {'a': 2.0, 'b': 0.05},
    }
    
    # Calculate limits for each depth
    limits_meters = {}
    for order_name, params in orders.items():
        limits_meters[order_name] = params['a'] + params['b'] * ref_depths
        limits_percent = (limits_meters[order_name] / ref_depths) * 100
        
        # Count compliant points
        compliant = np.abs(dz_wd) <= limits_percent
        compliance_rate = np.sum(compliant) / len(dz_wd) * 100
        
        # Add text to plot
        text = f"{order_name}: {compliance_rate:.1f}%"
        ax.text(0.02, 0.98 - 0.05 * list(orders.keys()).index(order_name), 
                text, transform=ax.transAxes, fontsize=8, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


def plot_accuracy(self, set_active_tab=False):  # plot the accuracy results
	# set point size; slider is on [1-11] for small # of discrete steps
	# print('in plot_accuracy with xline keys=', self.xline.keys())

	if not all([k in self.xline.keys() for k in ['beam_angle', 'dz_ref_wd']]):
		update_log(self, 'Beam angle and depth difference not found; crossline results will not be plotted')
		# Still call plot_bathymetric_orders to show appropriate message if orders are enabled
		plot_bathymetric_orders(self)
		return

	beam_bin_centers = np.asarray([b + self.beam_bin_size / 2 for b in self.beam_range])  # bin centers for plot
	beam_bin_dz_wd_std = np.asarray(self.beam_bin_dz_wd_std)  # std as a percent of wd across all bins
	self.N_plotted = 0

	# print('before calling decimate_data, the lens of beam_angle and dz_ref_wd = ', len(self.xline['beam_angle']),
	# 	  len(self.xline['dz_ref_wd']))

	nan_idx = np.isnan(self.xline['dz_ref_wd'])  # index of nans for all sounding dz values
	# print('got nan_idx with sum of non-nans=', np.sum(~nan_idx))
	# real_beam_angle = np.asarray(self.xline['beam_angle'])[~nan_idx].tolist()
	# real_dz_ref_wd = np.asarray(self.xline['dz_ref_wd'])[~nan_idx].tolist()
	# real_dz_ref = np.asarray(self.xline['dz_ref'])[~nan_idx].tolist()

	# apply crossline filter to data
	print('applying crossline filters')
	real_filt_idx = np.logical_and(~nan_idx, self.xline['filter_idx'])  # indices of real, filtered soundings
	
	# Apply bin count filter if available
	if hasattr(self, 'bin_count_mask'):
		print('DEBUG: bin_count_mask exists, shape:', self.bin_count_mask.shape)
		print('DEBUG: bin_count_mask sum (valid soundings):', np.sum(self.bin_count_mask))
		real_filt_idx = np.logical_and(real_filt_idx, self.bin_count_mask)
		print('Applied bin count filter, remaining soundings:', np.sum(real_filt_idx))
	else:
		print('DEBUG: bin_count_mask does not exist')
	real_beam_angle = np.asarray(self.xline['beam_angle'])[real_filt_idx].tolist()  # real, filtered beam angles
	real_dz_ref_wd = np.asarray(self.xline['dz_ref_wd'])[real_filt_idx].tolist()  # real, filtered dz results as %WD
	real_dz_ref = np.asarray(self.xline['dz_ref'])[real_filt_idx].tolist()  # real, filtered dz results in meters

	# calculate mean bias of filtered results across whole swath
	# offsets to apply to each real, filtered sounding to plot with zero mean
	print('size of real_filt_idx is', np.size(real_filt_idx))
	print('len of self.dz_wd_bin_zero_mean is', len(self.dz_wd_bin_zero_mean))
	real_dz_wd_bin_zero_mean = np.asarray(self.dz_wd_bin_zero_mean)[real_filt_idx].tolist()  # real, filtered dz results as %WD when forcing mean to zero
	# print('got real_dz_wd_bin_zero_mean =', real_dz_wd_bin_zero_mean)
	print('len of real_dz_wd_bin_zero_mean is', len(real_dz_wd_bin_zero_mean))
	### THIS WORKS FOR CALCULATING MEAN ACROSS ALL BINS
	real_dz_wd_bin_mean_mean = np.mean(real_dz_wd_bin_zero_mean)  # mean of filtered bin biases
	print('got initial value for real_dz_wd_bin_mean_mean = ', real_dz_wd_bin_mean_mean)

	# calculate mean bias of filtered results within user defined swath portion
	# print('getting real_dz_wd_bin_mean as array of')
	# real_filt_user_mean_angle_idx = np.logical_and(real_filt_idx, user_mean_angle_idx)
	# real_dz_wd_bin_user_mean = np.asarray(self.dz_wd_bin_zero_mean)[user_mean_angle_idx].tolist()  # real, filtered dz results as %WD with mean from user-defined swath portion

	# if self.flatten_mean_gb.isChecked():  # calculate mean within user defined swath limits, if checked
	# idx of angles that fall within user defined swath portion to use for mean bias calcs, if desired / checked
	if self.flatten_mean_gb.isChecked():
		print('*** flatten_mean_gb IS CHECKED****')

		if not self.set_zero_mean_chk.isChecked():
			print('***************** set_zero_mean_chk is NOT CHECKED')
			print('\n\n***getting indices of soundings within user defined beam angles')
			user_mean_angle_idx = np.logical_and(np.greater_equal(self.xline['beam_angle'], -1*float(self.mean_bias_angle_lim_tb.text())),
												 np.less_equal(self.xline['beam_angle'], float(self.mean_bias_angle_lim_tb.text())))
			print('got user_mean_angle_idx!')
			real_filt_user_mean_angle_idx = np.logical_and(real_filt_idx, user_mean_angle_idx)
			print('got real_filt_user_mean_angle_idx!')
			real_dz_wd_bin_user_mean = np.asarray(self.dz_wd_bin_zero_mean)[real_filt_user_mean_angle_idx].tolist()  # real, filtered dz results as %WD with mean from user-defined swath portion
			print('got real_dz_wd_bin_user_mean!')
			real_dz_wd_bin_mean_mean = np.mean(real_dz_wd_bin_user_mean)  # store mean within user selected swath limits
			print('survived user defined swath portion mean calc...')

	print('*******got real_dz_wd_bin_mean_mean = ', real_dz_wd_bin_mean_mean)

	# print('got real_beam_angle with len=', len(real_beam_angle), ' and = ', real_beam_angle)
	# print('got real_dz_ref_wd with len=', len(real_dz_ref_wd), ' and = ', real_dz_ref_wd)
	# print('got real_dz_ref with len=', len(real_dz_ref), ' and = ', real_dz_ref)
	# print('got real_dz_wd_bin_mean with len=', len(real_dz_wd_bin_mean), ' and = ', real_dz_wd_bin_mean)
	# print('got real_dz_wd_bin_mean_mean = ', real_dz_wd_bin_mean_mean)

	dec_data = decimate_data(self, data_list=[real_beam_angle, real_dz_ref_wd, real_dz_wd_bin_zero_mean])

	# print('back from decimate_data, dz_dec has len=', len(dec_data))
	real_beam_angle_dec = dec_data[0]
	real_dz_ref_wd_dec = dec_data[1]
	real_dz_wd_bin_zero_mean_dec = dec_data[2]

	real_dz_wd_bin_zero_mean_plot = [x - y for x, y in zip(real_dz_ref_wd_dec, real_dz_wd_bin_zero_mean_dec)]  # offsets to real filtered soundings for flat curve forced to zero mean
	real_dz_wd_bin_flat_mean_plot = [x + real_dz_wd_bin_mean_mean for x in real_dz_wd_bin_zero_mean_plot]  # offsets to real filtered soundings for flat curve preserving mean bias

	# print('beam_angle_dec and dz_ref_wd_dec have lens=', len(real_beam_angle_dec), len(real_dz_ref_wd_dec))
	self.N_plotted = len(real_beam_angle_dec)

	# plot standard deviation as %WD versus beam angle
	self.ax1.plot(beam_bin_centers, beam_bin_dz_wd_std, '-', linewidth=self.lwidth, color='b')  # bin mean + st. dev.

	# plot mean and std trend lines over decimated raw soundings (as %WD)
	if self.flatten_mean_gb.isChecked():  # flatten the mean bias curve

		print('\n\n************* FLATTEN MEAN GROUPBBOX IS CHECKED***************\n\n')

		if self.set_zero_mean_chk.isChecked():  # plot adjustment to zero mean (remove all mean biases)

			print('\n\n################# ZERO MEAN IS CHECKED ######################\n\n')
			# plot the raw differences, mean, and +/- 1 sigma as %wd versus beam angle, forcing the mean bias to zero
			self.ax2.scatter(real_beam_angle_dec, real_dz_wd_bin_zero_mean_plot,
							 marker='o', color='0.75', s=self.pt_size, alpha=self.pt_alpha)

			self.ax2.plot(beam_bin_centers, np.asarray(self.beam_bin_dz_wd_zero), '-',
						  linewidth=self.lwidth, color='r')  # beamwise bin mean diff

			self.ax2.plot(beam_bin_centers, self.beam_bin_dz_wd_std, '-',
						  linewidth=self.lwidth, color='b')  # beamwise bin mean + st. dev.
			self.ax2.plot(beam_bin_centers, np.multiply(-1, self.beam_bin_dz_wd_std), '-',
						  linewidth=self.lwidth, color='b')  # beamwise bin mean - st. dev.

			print('survived zero mean steps...')

		else:  # flatten curve but preserve the total mean across all bins (e.g., an across-the-board waterline error)
			# plot the raw differences, mean, and +/- 1 sigma as %wd versus beam angle, keeping the mean bias

			print('**** ZERO mean is NOT checked... keeping mean from user-defined portion ')
			self.ax2.scatter(real_beam_angle_dec, real_dz_wd_bin_flat_mean_plot,
							 marker='o', color='0.75', s=self.pt_size, alpha=self.pt_alpha)

			self.ax2.plot(beam_bin_centers, np.add(np.asarray(self.beam_bin_dz_wd_zero), real_dz_wd_bin_mean_mean), '-',
						  linewidth=self.lwidth, color='r')  # beamwise bin mean diff

			self.ax2.plot(beam_bin_centers, np.add(real_dz_wd_bin_mean_mean, self.beam_bin_dz_wd_std), '-',
						  linewidth=self.lwidth, color='b')  # beamwise bin mean + st. dev.
			self.ax2.plot(beam_bin_centers, np.subtract(real_dz_wd_bin_mean_mean, self.beam_bin_dz_wd_std), '-',
						  linewidth=self.lwidth, color='b')  # beamwise bin mean - st. dev.

			print('survived flattening without zero mean...')

	else:  # plot the raw differences with no flattening or adjustment of mean biases

		print('\n\n----------> FLATTEN MEAN IS NOT CHECKED\n\n')

		# plot the raw differences, mean, and +/- 1 sigma as %wd versus beam angle
		self.ax2.scatter(real_beam_angle_dec, real_dz_ref_wd_dec,
						 marker='o', color='0.75', s=self.pt_size, alpha=self.pt_alpha)

		# raw differences from reference grid, small gray points
		self.ax2.plot(beam_bin_centers, self.beam_bin_dz_wd_mean, '-',
					  linewidth=self.lwidth, color='r')  # beamwise bin mean diff
		self.ax2.plot(beam_bin_centers, np.add(self.beam_bin_dz_wd_mean, self.beam_bin_dz_wd_std), '-',
					  linewidth=self.lwidth, color='b')  # beamwise bin mean + st. dev.
		self.ax2.plot(beam_bin_centers, np.subtract(self.beam_bin_dz_wd_mean, self.beam_bin_dz_wd_std), '-',
					  linewidth=self.lwidth, color='b')  # beamwise bin mean - st. dev.

		print('survived normal plotting!')

			# Plot bathymetric order lines if enabled
		plot_bathymetric_orders(self)
		
		# Plot soundings per bin
		plot_soundings_per_bin(self)

	# update system info with xline detections
	# update_system_info(self, self.det, force_update=True, fname_str_replace='_trimmed')

	if set_active_tab:
		self.plot_tabs.setCurrentIndex(0)  # show accuracy results tab

def decimate_data_by_bins(self, data_list=[]):
	# decimate data by beam angle bins to achieve point count limit per bin
	# optional input: list of data to decimate (each item is a list with same length)
	# otherwise, return indices to apply for decimation of
	
	if not hasattr(self, 'xline') or not self.xline or 'beam_angle' not in self.xline:
		print('No beam angle data available for bin-based decimation, returning original data')
		return data_list if data_list else list(range(len(self.xline['dz_ref_wd'])))
	
	# Get the maximum points per bin from UI
	max_points_per_bin = 500  # default
	if hasattr(self, 'max_points_per_bin_tb') and self.max_points_per_bin_tb.text():
		try:
			max_points_per_bin = int(self.max_points_per_bin_tb.text())
		except ValueError:
			max_points_per_bin = 500
	
	print(f'Applying bin-based decimation with max {max_points_per_bin} points per bin')
	
	# Use the same beam bin parameters as the accuracy calculation
	beam_bin_size = 1  # beam angle bin size (deg)
	beam_bin_lim = 75  # max angle (deg)
	beam_range = range(-1 * beam_bin_lim, beam_bin_lim, beam_bin_size)
	
	if data_list:
		# When data_list is provided, we're working with already filtered data
		# The first element of data_list should be the beam angles for the filtered data
		self.n_points = len(data_list[0])
		
		# Use the beam angles from the data_list (first element)
		beam_angles = np.asarray(data_list[0])
		print(f'Using beam angles from data_list, length: {len(beam_angles)}')
	else:
		# When no data_list, work with full dataset
		self.n_points = len(self.xline['dz_ref_wd'])
		beam_angles = np.asarray(self.xline['beam_angle'])
		if len(beam_angles) != self.n_points:
			print('Beam angle array length does not match data length, returning original data')
			return list(range(len(self.xline['dz_ref_wd'])))
	
	# Initialize output indices (all points initially included)
	selected_indices = []
	
	# Process each beam angle bin
	for b in beam_range:
		# Find indices of angles in this bin
		beam_idx = (beam_angles >= b) & (beam_angles < b + beam_bin_size)
		bin_indices = np.where(beam_idx)[0]
		
		if len(bin_indices) > max_points_per_bin:
			# Randomly sample down to max_points_per_bin
			selected_bin_indices = np.random.choice(bin_indices, max_points_per_bin, replace=False)
			selected_indices.extend(selected_bin_indices)
			print(f'Bin {b}° to {b+beam_bin_size}°: {len(bin_indices)} -> {max_points_per_bin} points')
		else:
			# Keep all points in this bin
			selected_indices.extend(bin_indices)
			if len(bin_indices) > 0:
				print(f'Bin {b}° to {b+beam_bin_size}°: {len(bin_indices)} points (no decimation)')
	
	selected_indices = sorted(selected_indices)  # Sort for consistent ordering
	
	if data_list:
		# Return decimated data lists
		data_list_out = []
		for data in data_list:
			data_list_out.append([data[i] for i in selected_indices])
		return data_list_out
	else:
		# Return indices for decimation
		return selected_indices


def decimate_data(self, data_list=[]):
	# decimate data to achieve point count limit or apply user input (if selected)
	# optional input: list of data to decimate (each item is a list with same length)
	# otherwise, return indices to apply for decimation of
	# self.n_points = len(self.xline['dz_ref_wd'])
	
	# Check if bin-based decimation is enabled
	if hasattr(self, 'bin_decimation_gb') and self.bin_decimation_gb.isChecked():
		print('Using bin-based decimation')
		return decimate_data_by_bins(self, data_list)

	if data_list:
		self.n_points = len(data_list[0])
		data_list_out = data_list  # data_list_out will be decimated later if necessary
		# print('in decimate_data, got len(data_list[0]) = ', self.n_points)

	else:
		# print('in decimate_data, no data_list provided')
		self.n_points = len(self.xline['dz_ref_wd'])
		idx_out = [int(i) for i in range(self.n_points)]  # idx_out will be reduced later if necessary
		# print('got n_points = ', self.n_points)

	# print(1)
	self.n_points_max = self.n_points_max_default

	if self.pt_count_gb.isChecked() and self.max_count_tb.text():  # override default only if explicitly set by user
		print('setting n_points_max = max_count_tb.text')
		self.n_points_max = float(self.max_count_tb.text())

	# print(2)
	# default dec fac to meet n_points_max, regardless of whether user has checked box for plot point limits
	if self.n_points_max == 0:
		update_log(self, 'WARNING: Max plotting sounding count set equal to zero', font_color='red')
		self.dec_fac_default = np.inf
		# print(3)
	else:
		# print('setting dec_fac_default!')
		self.dec_fac_default = float(self.n_points / self.n_points_max)

		# print('self.dec_fac_default =', self.dec_fac_default)

	# print(4)
	if self.dec_fac_default > 1 and not self.pt_count_gb.isChecked():  # warn user if large count may slow down plot
		update_log(self, 'Large filtered sounding count (' + str(self.n_points) + ') may slow down plotting')

	# print(5)

	# get user dec fac as product of whether check box is checked (default 1)
	self.dec_fac_user = max(self.pt_count_gb.isChecked() * float(self.dec_fac_tb.text()), 1)
	self.dec_fac = max(self.dec_fac_default, self.dec_fac_user)

	if self.dec_fac_default > self.dec_fac_user:  # warn user if default max limit was reached
		update_log(self, 'Decimating crossline data (for plotting only) by factor of ' + "%.1f" % self.dec_fac +
				   ' to keep plotted point count under ' + "%.0f" % self.n_points_max +
				   '; accuracy results include all soundings, as filtered according to user input')

	elif self.pt_count_gb.isChecked() and self.dec_fac_user > self.dec_fac_default and self.dec_fac_user > 1:
		# otherwise, warn user if their manual dec fac was applied because it's more aggressive than max count
		update_log(self, 'Decimating crossline data (for plotting only) by factor of ' + "%.1f" % self.dec_fac +
				   ' per user input; accuracy results are include all soundings, as filtered according to user input')

	# print('before decimation, c_all=', c_all)

	if self.dec_fac > 1:
		print('dec fac > 1, trying to determine idx_dec')
		# print('dec_fac > 1 --> attempting interp1d')
		# n_points = len(self.xline['dz_ref_wd'])
		idx_all = np.arange(self.n_points)  # integer indices of all filtered data
		idx_dec = np.arange(0, self.n_points - 1, self.dec_fac)  # desired decimated indices, may be non-integer

		# print(6)
		# interpolate indices of colors, not color values directly
		f_dec = interp1d(idx_all, idx_all, kind='nearest')  # nearest neighbor interpolation function of all indices
		idx_out = [int(i) for i in f_dec(idx_dec)]  # list of decimated integer indices

		# print(7)

		# print(8)
		if data_list:
			# print('yes, data_list exists')
			for i, d in enumerate(data_list):
				# print('enumerating i=', i, 'applying idx_out to data_list')
				# print('d = ', d)
				# print('idx_out =', idx_out)
				data_list_out[i] = [d[j] for j in idx_out]
				# print('survived')

		# else:  # if no
		# 	data_list_out = data_list

	# print(9)

	if data_list:  # return list of decimated data if data were provided
		print('returning data_list_out from decimate_data')
		return data_list_out

	else:  # return indices for decimation if data were not provided
		print('returning idx_new from decimate_data')
		return idx_out

	# self.n_points = len(y_all)
	# print('self n_points = ', self.n_points)


def sort_xline_track(self, new_track):
	# pull ship track from dict of parsed crossline track information (different data from .all and .kmall files) and
	# convert to current UTM zone
	# .all files: use active position sensor, full position time series
	# .kmall files: use ping position (assumes active position sensor), ping time only

	track_out = {}  # simplified trackline dict with lat, lon, easting, northing, utm_zone

	print('\n\n\n******in sort_xline_track with utm_zone =', self.ref['utm_zone'])

	# refProj = pyproj.Proj(proj='utm', zone=self.ref['utm_zone'], ellps='WGS84')  # define output projection
	utm_zone = int(''.join([c for c in self.ref['utm_zone'] if c.isnumeric()]))  # need integer for pyproj.Proj
	south_hem = True if self.ref['utm_zone'][-1].lower() == 's' else False  # need hemisphere for pyproj.Proj

	print('in sort_xline_track, sending zone = ', utm_zone, ' and south_hem =', south_hem, ' to pyproj.Proj')
	refProj = pyproj.Proj(proj='utm', zone=utm_zone, ellps='WGS84', south=south_hem)  # update UTM zone with hemisphere

	for f in range(len(new_track)):
		lat, lon = [], []
		fname = new_track[f]['fname']
		# track_out[fname] = {}

		if '.all' in new_track[f]['fname']:
			temp = {0: dict(new_track[f])}  # reformat dict with key=0 for sort_active_pos_system
			dt_pos, lat, lon, sys_num = sort_active_pos_system(temp, print_updates=True)  # use only active pos

		elif '.kmall' in new_track[f]['fname']:
			# print('new_track[f].keys() = ', new_track[f].keys())
			pingInfo = new_track[f]['RTP']  # temp method using pingInfo lat/lon stored in RTP
			headerInfo = new_track[f]['HDR']
			# print('pingInfo has len =', len(pingInfo))
			# print('headerInfo has len =', len(headerInfo))

			lat = [pingInfo[p]['latitude_deg'] for p in range(len(pingInfo))]
			lon = [pingInfo[p]['longitude_deg'] for p in range(len(pingInfo))]
			dt_pos = [headerInfo[p]['dgdatetime'] for p in range(len(headerInfo))]

			# print('pingInfo[0].keys() = ', pingInfo[0].keys())
			# print('headerInfo[0].keys() = ', headerInfo[0].keys())

			# dt_pos = [pingInfo[p]['dgdatetime'] for p in range(len(pingInfo))]

		elif 'ASCII.txt' in new_track[f]['fname']:
			lat = np.unique(new_track[f]['ping_lat']).tolist()
			lon = np.unique(new_track[f]['ping_lon']).tolist()
			dt_pos = np.unique(new_track[f]['datetime']).tolist()

		else:
			print('in sort_xline_track, not sorting track for f =', f, '-->', new_track[f]['fname'])
			continue

		print('first couple lat, lon are', lat[0:10], lon[0:10])
		temp_out = {}
		temp_out['lat'] = lat
		temp_out['lon'] = lon
		temp_out['datetime'] = dt_pos
		temp_out['e'], temp_out['n'] = refProj(lon, lat)
		temp_out['utm_zone'] = self.ref['utm_zone']

		track_out[fname] = temp_out

		print('for fname =', fname, 'the first 10 track e, n =',
			  track_out[fname]['e'][0:10], track_out[fname]['n'][0:10])
		print('for fname =', fname, ' sorted trackline temp_out[utm_zone] =', temp_out['utm_zone'])

	return track_out


def refresh_plot(self, refresh_list=['ref', 'acc', 'tide'], sender=None, set_active_tab=None):
	# update swath plot with new data and options
	from time import process_time
	start_time = process_time()
	
	print('refresh_plot called from sender=', sender, ', refresh_list=', refresh_list, ', active_tab=', set_active_tab)
	update_log(self, f'=== REFRESHING PLOTS === (sender: {sender}, plots: {refresh_list})')
	
	print('calling clear_plot from refresh_plot')
	clear_start = process_time()
	clear_plot(self, refresh_list)
	clear_time = process_time() - clear_start
	update_log(self, f'Plot clearing completed in {clear_time:.2f} seconds')
	print('cleared plot in refresh_plot')
	
	self.pt_size = np.square(float(self.pt_size_cbox.currentText()))  # swath plot point size
	self.pt_size_cov = np.square(float(self.pt_size_cov_cbox.currentText()))  # coverage plot point size
	self.pt_alpha = np.divide(float(self.pt_alpha_acc_tb.text()), 100)  # accuracy plot opacity
	self.pt_alpha_cov = np.divide(float(self.pt_alpha_cov_tb.text()), 100)  # coverage plot opacity
	update_log(self, f'Plot settings - Point size: {self.pt_size}, Coverage size: {self.pt_size_cov}, Alpha (acc): {self.pt_alpha:.2f}, Alpha (cov): {self.pt_alpha_cov:.2f}')
	print('got pt_size, pt_size_cov, pt_alpha (acc), and pt_alpha_cov = ', self.pt_size, self.pt_size_cov, self.pt_alpha, self.pt_alpha_cov)

	try:
		print('in refresh_plot, calling update_axes')
		axes_start = process_time()
		update_axes(self)
		axes_time = process_time() - axes_start
		update_log(self, f'Axes update completed in {axes_time:.2f} seconds')
		
		print('in refresh_plot, calling add_grid_lines')
		grid_start = process_time()
		add_grid_lines(self)  # add grid lines
		grid_time = process_time() - grid_start
		update_log(self, f'Grid lines added in {grid_time:.2f} seconds')
		print('in refresh_plot, survived add_grid_lines')

		if 'ref' in refresh_list:
			update_log(self, 'Generating reference surface plots...')
			ref_start = process_time()
			print('refreshing ref plot')
			update_plot_limits(self)
			print('1')
			plot_ref_surf(self)
			print('2')
			self.surf_canvas.draw()
			print('3')
			self.surf_final_canvas.draw()
			print('4')
			# Also refresh the density final tab canvas if it exists
			if hasattr(self, 'density_final_canvas'):
				self.density_final_canvas.draw()
			print('5')
			# Also refresh the slope final tab canvas if it exists
			if hasattr(self, 'slope_final_canvas'):
				self.slope_final_canvas.draw()
			print('6')
			# Also refresh the uncertainty final tab canvas if it exists
			if hasattr(self, 'uncertainty_final_canvas'):
				self.uncertainty_final_canvas.draw()
			print('7')
			# Also refresh the depth tab canvas if it exists
			if hasattr(self, 'depth_canvas'):
				self.depth_canvas.draw()
			print('8')
			plt.show()
			ref_time = process_time() - ref_start
			update_log(self, f'Reference surface plots completed in {ref_time:.2f} seconds')
			print('finished refreshing ref in refresh_list')

		if 'acc' in refresh_list:
			update_log(self, 'Generating accuracy plots...')
			acc_start = process_time()
			print('refreshing acc plot')
			print('a')
			update_plot_limits(self)
			print('b')
			plot_accuracy(self)
			print('c)')
			self.swath_canvas.draw()
			print('d')
			plt.show()
			acc_time = process_time() - acc_start
			update_log(self, f'Accuracy plots completed in {acc_time:.2f} seconds')
			print('finished refreshing acc in refresh_list')

		if 'tide' in refresh_list:
			update_log(self, 'Generating tide plots...')
			tide_start = process_time()
			# plot_tide(self)
			print('refreshing tide plot')
			print('a1')
			plot_tide(self)
			print('b2')
			self.tide_canvas.draw()
			print('c3')
			plt.show()
			tide_time = process_time() - tide_start
			update_log(self, f'Tide plots completed in {tide_time:.2f} seconds')
			print('finished refreshing tide in refresh_list')
				
		total_plot_time = process_time() - start_time
		update_log(self, f'=== PLOT REFRESH COMPLETED in {total_plot_time:.2f} seconds ===')
		print('survived calling plot steps from refresh_plot')

	except Exception as e:
		import traceback
		error_msg = f'Error in refreshing plot: {str(e)}\n{traceback.format_exc()}'
		update_log(self, error_msg, font_color="red")
		print(error_msg)


	# Determine which tab to show based on whether accuracy data has been calculated
	if set_active_tab != None:  # if explicitly set, use that value
		self.plot_tabs.setCurrentIndex(set_active_tab)
		update_log(self, f'Set active tab to index {set_active_tab}')
	else:
		# Check if accuracy data has been calculated
		accuracy_calculated = (hasattr(self, 'xline') and 
							  self.xline and 
							  'dz_ref_wd' in self.xline and 
							  len(self.xline['dz_ref_wd']) > 0)
		
		if accuracy_calculated:
			# Show Accuracy tab (index 0) when accuracy data is available
			self.plot_tabs.setCurrentIndex(0)
			update_log(self, 'Set active tab to Accuracy (data available)')
		else:
			# Show Surface Filters tab (index 1) when no accuracy data is available
			self.plot_tabs.setCurrentIndex(1)
			update_log(self, 'Set active tab to Surface Filters (no accuracy data)')

	update_buttons(self)


def update_axes(self):
	# udpate axes for swath and tide plots; ref surf axes are handled in plot_ref_surf
	# update top subplot axes (std. dev. as %WD)
	# update_system_info(self)
	update_system_info(self, self.xline, force_update=False, fname_str_replace='_trimmed')
	update_plot_limits(self)

	# set x axis limits and ticks for both swath subplots
	plt.setp((self.ax1, self.ax2),
			 xticks=np.arange(-1 * self.x_max, self.x_max + self.x_spacing, self.x_spacing),
			 xlim=(-1 * self.x_max, self.x_max))

	# set y axis limits for both swath subplots
	self.ax1.set_ylim(0, self.y_max_std)  # set y axis for std (0 to max, top plot)
	self.ax2.set_ylim(-1 * self.y_max_bias, self.y_max_bias)  # set y axis for total bias+std (bottom plot)

	title_str = 'Swath Accuracy vs. Beam Angle'
	title_str_surf = 'Reference Surface'
	title_str_tide = 'Tide Applied to Accuracy Crosslines'

	# update plot title with default or custom combination of system info fields
	if self.custom_info_gb.isChecked():  # include custom system info that is checked on
		sys_info_list = [['', self.model_name][self.show_model_chk.isChecked()],
						 ['', self.ship_name][self.show_ship_chk.isChecked()],
						 ['', self.cruise_name][self.show_cruise_chk.isChecked()]]
		print('got sys_info_list = ', sys_info_list)
		sys_info_str = ' - '.join([str for str in sys_info_list if str != ''])

	else:  # otherwise, default to all system info in the title
		sys_info_str = ' - '.join([self.model_name, self.ship_name, self.cruise_name])

	# get set of modes in these crosslines
	if self.xline:  # get set of modes in these crosslines and add to title string
		try:
			# default set of modes based on all available in crossline files
			# modes = [' / '.join([self.xline['ping_mode'][i],
			# 					 self.xline['swath_mode'][i],
			# 					 self.xline['pulse_form'][i]]) for i in range(len(self.xline['ping_mode']))]

			# set of modes based on filters (e.g., if filtered for mode(s)); np.where returns tuple, use [0] as list
			print('trying to get mode set: np.where(self.xline[filter_idx]) =', np.where(self.xline['filter_idx'])[0])
			modes = [' / '.join([self.xline['ping_mode'][i],
								 self.xline['swath_mode'][i],
								 self.xline['pulse_form'][i]]) for i in np.where(self.xline['filter_idx'])[0]]

			modes_str = ' + '.join(list(set(modes)))

		except:
			modes_str = 'Modes N/A'
	else:
		modes_str = 'Modes N/A'

	if self.ref:
		fname_ref = self.ref['fname']
	else:
		fname_ref = 'Reference file N/A'

	if self.tide:
		fname_tide = self.tide['fname']
	else:
		fname_tide = 'Tide file N/A'

	self.title_str = '\n'.join([title_str, sys_info_str, modes_str])
	self.title_str_ref = '\n'.join([title_str_surf, sys_info_str, fname_ref])
	self.title_str_tide = '\n'.join([title_str_tide, sys_info_str, fname_tide])

	# set axis labels
	self.ax1.set(xlabel='Beam Angle (deg, pos. stbd.)', ylabel='Depth Bias Std. Dev (% Water Depth)')
	self.ax2.set(xlabel='Beam Angle (deg, pos. stbd.)', ylabel='Depth Bias (% Water Depth)')
	self.tide_ax.set(xlabel='Time', ylabel='Tide Ampitude (m from tide file datum, positive up)')

	# set super titles
	self.swath_figure.suptitle(self.title_str)
	self.surf_figure.suptitle(self.title_str_ref)
	# Use ax.set_title() for final surface to avoid large gaps
	self.surf_ax5.set_title(self.title_str_ref)
	self.tide_figure.suptitle(self.title_str_tide)

	# add processing text boxes
	add_ref_proc_text(self)
	add_xline_proc_text(self)


	try:
		plt.show()  # need show() after update; failed until matplotlib.use('qt5agg') added at start

	except:
		print('in update_axes, failed plt.show()')

	print('survived update_axes')
	print('*******')


def update_plot_limits(self):
	# update plot limits if custom limits are selected
	if self.plot_lim_gb.isChecked():  # use custom plot limits if checked, store custom values in text boxes
		# self.max_gb.setEnabled(True)
		self.x_max_custom = float(self.max_beam_angle_tb.text())
		self.x_spacing_custom = float(self.angle_spacing_tb.text())
		self.y_max_std_custom = float(self.max_std_tb.text())
		self.y_max_bias_custom = float(self.max_bias_tb.text())
		self.axis_margin_custom = float(self.axis_margin_tb.text())

		# assign to current plot limits
		self.x_max = self.x_max_custom
		self.x_spacing = self.x_spacing_custom
		self.y_max_std = self.y_max_std_custom
		self.y_max_bias = self.y_max_bias_custom
		self.axis_margin = self.axis_margin_custom

	else:  # revert to default limits from the data if unchecked, but keep the custom numbers in text boxes
		# self.max_gb.setEnabled(False)
		self.x_max = self.x_max_default
		self.x_spacing = self.x_spacing_default
		self.y_max_std = self.y_max_std_default
		self.y_max_bias = self.y_max_bias_default
		self.axis_margin = self.axis_margin_default

		print('axis margin is now ', self.axis_margin)

		# set text boxes to latest custom values for easier toggling between custom/default limits
		self.max_beam_angle_tb.setText(str(float(self.x_max_custom)))
		self.angle_spacing_tb.setText(str(float(self.x_spacing_custom)))
		self.max_bias_tb.setText(str(float(self.y_max_bias_custom)))
		self.max_std_tb.setText(str(float(self.y_max_std_custom)))
		self.axis_margin_tb.setText(str(float(self.axis_margin_custom)))


def add_grid_lines(self):
	for ax in [self.ax1, self.ax2, self.surf_ax1, self.surf_ax2, self.surf_ax3, self.surf_ax4,
		   self.surf_ax5, self.depth_ax, self.density_final_ax, self.slope_final_ax, self.uncertainty_final_ax, self.tide_ax]:
		if self.grid_lines_toggle_chk.isChecked():
			ax.grid()
			ax.minorticks_on()
			ax.grid(which='minor', linestyle='-', linewidth='0.5', color='black')
			ax.grid(which='major', linestyle='-', linewidth='1.0', color='black')

		else:
			ax.grid(False)
			ax.minorticks_off()


def add_xline_proc_text(self):
	# add text for depth ref and filters applied
	proc_str = 'Crossline processing parameters'
	proc_str += '\nSounding reference: ' + self.ref_cbox.currentText()
	proc_str += '\nReference surface: ' + (self.ref['fname'].split('/')[-1] if self.ref else 'N/A')
	proc_str += '\nUTM Zone: ' + self.ref_proj_cbox.currentText()
	proc_str += '\nNum. crosslines: ' + (str(len(list(set(self.xline['fname'])))) if self.xline else '0') + ' files'
	proc_str += '\nNum. soundings: ' + (str(len(self.xline['z'])) + ' (' + str(self.xline['num_on_ref']) +
										' on filtered ref. surf.)' if self.xline else '0')
	proc_str += '\nNum. plotted: ' + (str(self.N_plotted) + ' (filtered soundings)' if self.xline else '0')
	proc_str += '\nTide applied: ' + (self.tide['fname'] if self.tide and self.tide_applied else 'None')
	# make dict of text to include based on user input
	depth_fil_xline = ['None', self.min_depth_xline_tb.text() + ' to ' + self.max_depth_xline_tb.text() + ' m']
	dz_abs_fil = ['None', self.max_dz_tb.text() + ' m']
	dz_pct_fil = ['None', self.max_dz_wd_tb.text() + ' %WD']
	angle_fil = ['None', self.min_angle_xline_tb.text() + ' to ' + self.max_angle_xline_tb.text() + '\u00b0']
	bs_fil = ['None', ('+' if float(self.min_bs_xline_tb.text()) > 0 else '') + self.min_bs_xline_tb.text() + ' to ' +
			  ('+' if float(self.max_bs_xline_tb.text()) > 0 else '') + self.max_bs_xline_tb.text() + ' dB']

	fil_dict = {'Angle filter: ': angle_fil[self.angle_xline_gb.isChecked()],
				'Depth filter (crossline): ': depth_fil_xline[self.depth_xline_gb.isChecked()],
				'Backscatter filter: ': bs_fil[self.bs_xline_gb.isChecked()],
				'Max. diff. (m):': dz_abs_fil[hasattr(self, 'dz_abs_gb') and self.dz_abs_gb.isChecked()],
				'Max. diff. (%WD):': dz_pct_fil[hasattr(self, 'dz_pct_gb') and self.dz_pct_gb.isChecked()]}

	for fil in fil_dict.keys():
		proc_str += '\n' + fil + fil_dict[fil]

	if self.show_acc_proc_text_chk.isChecked():
		self.ax1.text(0.02, 0.98, proc_str,
					  ha='left', va='top', fontsize=6, transform=self.ax1.transAxes,
					  bbox=dict(facecolor='white', edgecolor=None, linewidth=0, alpha=0.8))


def add_ref_proc_text(self):
	print('made it to add_ref_proc_text')
	# add text for depth ref and filters applied
	proc_str = 'Reference surface processing parameters'
	# proc_str = 'Vertical reference: ' + self.ref_cbox.currentText()
	proc_str += '\nReference surface: ' + (self.ref['fname'].split('/')[-1] if self.ref else 'N/A')
	proc_str += '\nUTM Zone: ' + self.ref_utm_str

	# make dict of text to include based on user input
	# depth_fil_xline = ['None', self.min_depth_xline_tb.text() + ' to ' + self.max_depth_xline_tb.text() + ' m']
	depth_fil_ref = ['None', self.min_depth_ref_tb.text() + ' to ' + self.max_depth_ref_tb.text() + ' m']
	slope_win = ['None', self.slope_win_cbox.currentText() + ' cells']
	slope_fil = ['None', '0 to ' + self.max_slope_tb.text() + ' deg']

	fil_dict = {'Depth filter (reference): ': depth_fil_ref[self.depth_ref_gb.isChecked()],
				'Slope smoothing window: ': slope_win[self.slope_gb.isChecked()],
				'Slope filter (reference): ': slope_fil[self.slope_gb.isChecked()]}

	print('made it to end of fil_dict')

	for fil in fil_dict.keys():
		proc_str += '\n' + fil + fil_dict[fil]

	if self.show_ref_proc_text_chk.isChecked():
		self.surf_ax1.text(0.02, 0.98, proc_str,ha='left', va='top', fontsize=6, transform=self.surf_ax1.transAxes,
						   bbox=dict(facecolor='white', edgecolor=None, linewidth=0, alpha=0.8))


def save_plot(self):
	# option 1: try to save a .PNG of the swath plot
	plot_path = QtWidgets.QFileDialog.getSaveFileName(self, 'Save plot as...', os.getenv('HOME'),
													  ".PNG file (*.PNG);; .JPG file (*.JPG);; .TIF file (*.TIF)")
	fname_out = plot_path[0]

	if fname_out:  # Only save if a filename was selected
		try:
			self.swath_figure.savefig(fname_out,
						dpi=600, facecolor='w', edgecolor='k',
						transparent=False, bbox_inches='tight', pad_inches=0.1)
			update_log(self, 'Saved figure ' + fname_out.rsplit('/')[-1])
		except Exception as e:
			update_log(self, f'Error saving plot: {str(e)}')
			print(f"Error saving plot: {str(e)}")


def save_all_plots(self):
	# Save all plot tabs as separate files
	# Load last used plot settings
	config = load_session_config()
	last_basename = config.get("last_plot_basename", "swath_accuracy_plots")
	last_directory = config.get("last_plot_directory", os.getcwd())
	last_extension = config.get("last_plot_extension", ".PNG")
	
	# First, prompt for parent directory
	parent_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select parent directory for plots...', last_directory)
	
	if not parent_dir:  # User cancelled
		return
	
	# Then, prompt for basename with a custom dialog
	from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QDialogButtonBox
	
	# Create custom dialog for basename input
	dialog = QDialog(self)
	dialog.setWindowTitle('Enter basename for plots')
	dialog.setModal(True)
	dialog.setMinimumWidth(500)  # Make dialog wider
	
	# Create layout
	layout = QVBoxLayout()
	
	# Add label
	label = QLabel('Basename (without extension):')
	layout.addWidget(label)
	
	# Add input field
	input_field = QLineEdit()
	input_field.setText(last_basename)
	input_field.setMinimumWidth(400)  # Make input field wider
	layout.addWidget(input_field)
	
	# Add button box
	button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
	button_box.accepted.connect(dialog.accept)
	button_box.rejected.connect(dialog.reject)
	layout.addWidget(button_box)
	
	dialog.setLayout(layout)
	
	# Show dialog and get result
	result = dialog.exec()
	ok = result == QDialog.DialogCode.Accepted
	basename = input_field.text() if ok else ""
	
	if not ok or not basename.strip():  # User cancelled or entered empty string
		return
	
	basename = basename.strip()
	
	# Set default extension
	file_ext = last_extension
	
	# Debug: Print what we're using
	print(f"DEBUG: Parent directory: {parent_dir}")
	print(f"DEBUG: Basename: {basename}")
	print(f"DEBUG: Extension: {file_ext}")

	# --- NEW: Create a directory with the base name ---
	output_subdir = os.path.join(parent_dir, basename)
	if not os.path.exists(output_subdir):
		os.makedirs(output_subdir)

	# Define plot names and their corresponding figures
	plot_configs = [
		('accuracy', self.swath_figure),
		('surface_filters', self.surf_figure),
		('final_surface', self.surf_final_figure),
		('depth', self.depth_figure),
		('uncertainty', self.uncertainty_final_figure),
		('density', self.density_final_figure),
		('slope', self.slope_final_figure),
		('soundings', self.soundings_figure),
		('tide', self.tide_figure)
	]
	
	# Check for existing files and ask for overwrite permission
	existing_files = []
	for plot_name, figure in plot_configs:
		if figure and hasattr(figure, 'savefig'):
			filename = f"{basename}_{plot_name}{file_ext}"
			filepath = os.path.join(output_subdir, filename)
			if os.path.exists(filepath):
				existing_files.append(filename)
	
	# Ask user for permission to overwrite if files exist
	if existing_files:
		from PyQt6.QtWidgets import QMessageBox
		file_list = "\n".join([f"  • {f}" for f in existing_files])
		msg = f"The following files already exist:\n{file_list}\n\nDo you want to overwrite them?"
		reply = QMessageBox.question(self, 'Overwrite Files?', msg, 
									QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
									QMessageBox.StandardButton.No)
		if reply == QMessageBox.StandardButton.No:
			update_log(self, 'Save operation cancelled by user')
			return

	# Save the new settings for next session
	config["last_plot_basename"] = basename
	config["last_plot_directory"] = parent_dir  # Save the parent directory where the subdirectory was created
	config["last_plot_extension"] = file_ext
	
	# Debug: Print what we're saving
	print(f"DEBUG: Saving config - basename: {basename}, directory: {parent_dir}, extension: {file_ext}")
	
	save_session_config(config)
	
	saved_count = 0
	for plot_name, figure in plot_configs:
		if figure and hasattr(figure, 'savefig'):
			try:
				# Create filename by appending plot type to base name
				filename = f"{basename}_{plot_name}{file_ext}"
				filepath = os.path.join(output_subdir, filename)
				
				# Save the figure
				figure.savefig(filepath,
					  dpi=600, facecolor='w', edgecolor='k',
					  transparent=False, bbox_inches='tight', pad_inches=0.1)
				saved_count += 1
				update_log(self, f'Saved {plot_name} plot: {filename}')
			except Exception as e:
				update_log(self, f'Error saving {plot_name} plot: {str(e)}')
				print(f"Error saving {plot_name} plot: {str(e)}")
	
	# Save info text file
	try:
		info_filename = f"{basename}_info.txt"
		info_filepath = os.path.join(output_subdir, info_filename)
		
		with open(info_filepath, 'w') as f:
			f.write("SWATH ACCURACY PLOTTER - EXPORT INFORMATION\n")
			f.write("=" * 50 + "\n\n")
			
			# System Information
			f.write("SYSTEM INFORMATION:\n")
			f.write("-" * 20 + "\n")
			f.write(f"Ship Name: {getattr(self, 'ship_name', 'N/A')}\n")
			f.write(f"Cruise Name: {getattr(self, 'cruise_name', 'N/A')}\n")
			f.write(f"Model Name: {getattr(self, 'model_name', 'N/A')}\n")
			
			# Get ping mode, swath mode, and pulse form from crossline data
			if hasattr(self, 'xline') and 'ping_mode' in self.xline and len(self.xline['ping_mode']) > 0:
				ping_modes = list(set(self.xline['ping_mode']))
				f.write(f"Ping Mode: {', '.join(ping_modes)}\n")
			else:
				f.write("Ping Mode: N/A\n")
				
			if hasattr(self, 'xline') and 'swath_mode' in self.xline and len(self.xline['swath_mode']) > 0:
				swath_modes = list(set(self.xline['swath_mode']))
				f.write(f"Swath Mode: {', '.join(swath_modes)}\n")
			else:
				f.write("Swath Mode: N/A\n")
				
			if hasattr(self, 'xline') and 'pulse_form' in self.xline and len(self.xline['pulse_form']) > 0:
				pulse_forms = list(set(self.xline['pulse_form']))
				f.write(f"Pulse Form: {', '.join(pulse_forms)}\n")
			else:
				f.write("Pulse Form: N/A\n")
			
			f.write("\n")
			
			# File Information
			f.write("FILE INFORMATION:\n")
			f.write("-" * 20 + "\n")
			
			# Reference surface
			if hasattr(self, 'ref') and 'fname' in self.ref:
				ref_fname = os.path.basename(self.ref['fname']) if '/' in self.ref['fname'] else self.ref['fname']
				f.write(f"Reference Surface: {ref_fname}\n")
			else:
				f.write("Reference Surface: N/A\n")
			
			# Density surface
			if hasattr(self, 'ref') and 'fname_dens' in self.ref:
				dens_fname = os.path.basename(self.ref['fname_dens']) if '/' in self.ref['fname_dens'] else self.ref['fname_dens']
				f.write(f"Density Surface: {dens_fname}\n")
			else:
				f.write("Density Surface: N/A\n")
			
			# Crosslines
			if hasattr(self, 'xline') and 'fname' in self.xline and len(self.xline['fname']) > 0:
				crossline_fnames = list(set(self.xline['fname']))
				f.write(f"Crosslines ({len(crossline_fnames)} files):\n")
				for fname in crossline_fnames:
					f.write(f"  - {os.path.basename(fname) if '/' in fname else fname}\n")
			else:
				f.write("Crosslines: N/A\n")
			
			f.write("\n")
			
			# Filter Settings
			f.write("FILTER SETTINGS:\n")
			f.write("-" * 20 + "\n")
			
			# Reference Surface Filters
			f.write("Reference Surface Filters:\n")
			# Depth filter
			if hasattr(self, 'depth_ref_gb') and self.depth_ref_gb.isChecked():
				min_depth = getattr(self, 'min_depth_ref_tb', None)
				max_depth = getattr(self, 'max_depth_ref_tb', None)
				if min_depth and max_depth:
					f.write(f"  - Depth: {min_depth.text()} to {max_depth.text()} m\n")
			else:
				f.write("  - Depth: Not applied\n")
			# Slope filter
			if hasattr(self, 'slope_gb') and self.slope_gb.isChecked():
				max_slope = getattr(self, 'max_slope_tb', None)
				if max_slope:
					f.write(f"  - Slope: 0 to {max_slope.text()} degrees\n")
			else:
				f.write("  - Slope: Not applied\n")
			# Density filter
			if hasattr(self, 'density_gb') and self.density_gb.isChecked():
				min_dens = getattr(self, 'min_dens_tb', None)
				if min_dens:
					f.write(f"  - Density: {min_dens.text()}+ soundings/cell\n")
			else:
				f.write("  - Density: Not applied\n")
			# Uncertainty filter
			if hasattr(self, 'uncertainty_gb') and self.uncertainty_gb.isChecked():
				max_u = getattr(self, 'max_u_tb', None)
				if max_u:
					f.write(f"  - Uncertainty: 0 to {max_u.text()} m\n")
			else:
				f.write("  - Uncertainty: Not applied\n")

			f.write("\n")

			# Crossline Filters
			f.write("Crossline Filters:\n")
			# Angle filter
			if hasattr(self, 'angle_xline_gb') and self.angle_xline_gb.isChecked():
				min_angle = getattr(self, 'min_angle_xline_tb', None)
				max_angle = getattr(self, 'max_angle_xline_tb', None)
				if min_angle and max_angle:
					f.write(f"  - Beam Angle: {min_angle.text()} to {max_angle.text()} degrees\n")
			else:
				f.write("  - Beam Angle: Not applied\n")
			# Depth filter (crossline)
			if hasattr(self, 'depth_xline_gb') and self.depth_xline_gb.isChecked():
				min_depth_x = getattr(self, 'min_depth_xline_tb', None)
				max_depth_x = getattr(self, 'max_depth_xline_tb', None)
				if min_depth_x and max_depth_x:
					f.write(f"  - Depth: {min_depth_x.text()} to {max_depth_x.text()} m\n")
			else:
				f.write("  - Depth: Not applied\n")
			# Depth difference filter
			if hasattr(self, 'dz_gb') and self.dz_gb.isChecked():
				max_dz = getattr(self, 'max_dz_tb', None)
				max_dz_wd = getattr(self, 'max_dz_wd_tb', None)
				if max_dz and max_dz_wd:
					f.write(f"  - Depth Difference: 0 to {max_dz.text()} m AND {max_dz_wd.text()}% of water depth (both criteria applied)\n")
			else:
				f.write("  - Depth Difference: Not applied\n")
			# Depth mode filter
			if hasattr(self, 'depth_mode_gb') and self.depth_mode_gb.isChecked():
				depth_mode = getattr(self, 'depth_mode_cbox', None)
				if depth_mode:
					f.write(f"  - Depth Mode: {depth_mode.currentText()}\n")
			else:
				f.write("  - Depth Mode: Not applied\n")

			f.write("\n")
			f.write(f"Export Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
	except Exception as e:
		update_log(self, f'Error saving info file: {str(e)}')
		print(f"Error saving info file: {str(e)}")


def clear_plot(self, refresh_list=['ref', 'acc', 'tide']):
	# clear plots in refresh_list
	print('in clear_plot with refresh_list=', refresh_list)
	if 'acc' in refresh_list:
		for ax in [self.ax1, self.ax2]:
			ax.clear()

		self.swath_canvas.draw()

	if 'ref' in refresh_list:
		for ax in [self.surf_ax1, self.surf_ax2, self.surf_ax3, self.surf_ax4, self.surf_ax5, self.depth_ax, self.density_final_ax, self.slope_final_ax, self.uncertainty_final_ax]:
			ax.clear()

		self.surf_canvas.draw()

	if 'tide' in refresh_list:
		self.tide_ax.clear()
		self.tide_canvas.draw()


def export_all_to_geotiff(self):
	"""Export all reference surface data as GeoTIFF files"""
	try:
		# Check if reference surface data is available
		if not hasattr(self, 'ref') or not self.ref:
			update_log(self, 'ERROR: No reference surface data available for export', font_color="red")
			return
		
		# Check if grids are available
		required_grids = ['z_grid', 'e_grid', 'n_grid']
		available_grids = []
		for grid in required_grids:
			if grid in self.ref and self.ref[grid] is not None:
				available_grids.append(grid)
		
		if not available_grids:
			update_log(self, 'ERROR: No reference surface grids available for export', font_color="red")
			return
		
		# Load last used directory
		config = load_session_config()
		last_directory = config.get("last_export_directory", os.getcwd())
		
		# First, prompt for parent directory
		parent_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select parent directory for GeoTIFF export...', last_directory)
		
		if not parent_dir:  # User cancelled
			return
		
		# Then, prompt for export name with a custom dialog
		from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QDialogButtonBox
		
		# Create custom dialog for export name input
		dialog = QDialog(self)
		dialog.setWindowTitle('Enter export name')
		dialog.setModal(True)
		dialog.setMinimumWidth(500)
		
		# Create layout
		layout = QVBoxLayout()
		
		# Add label
		label = QLabel('Export name (will create directory with this name):')
		layout.addWidget(label)
		
		# Add input field
		input_field = QLineEdit()
		input_field.setText('reference_surface_export')
		input_field.setMinimumWidth(400)
		layout.addWidget(input_field)
		
		# Add button box
		button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
		button_box.accepted.connect(dialog.accept)
		button_box.rejected.connect(dialog.reject)
		layout.addWidget(button_box)
		
		dialog.setLayout(layout)
		
		# Show dialog and get result
		result = dialog.exec()
		ok = result == QDialog.DialogCode.Accepted
		export_name = input_field.text() if ok else ""
		
		if not ok or not export_name.strip():  # User cancelled or entered empty string
			return
		
		export_name = export_name.strip()
		
		# Create output directory
		output_dir = os.path.join(parent_dir, export_name)
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)
		
		# Save last directory
		save_last_directory(parent_dir)
		update_last_directory("last_export_directory", parent_dir)
		
		update_log(self, f'Exporting reference surface data to: {output_dir}')
		
		# Import required libraries for GeoTIFF export
		try:
			import rasterio
			from rasterio.transform import from_bounds
			from rasterio.crs import CRS
		except ImportError:
			update_log(self, 'ERROR: rasterio library not available. Please install it with: pip install rasterio', font_color="red")
			return
		
		# Get the extent and CRS information
		if 'z_ref_extent' in self.ref:
			extent = self.ref['z_ref_extent']
			# extent format: [xmin, xmax, ymin, ymax]
			bounds = [extent[0], extent[2], extent[1], extent[3]]  # [xmin, ymin, xmax, ymax]
		else:
			update_log(self, 'ERROR: Reference surface extent not available', font_color="red")
			return
		
		# Get UTM zone for CRS
		utm_zone = self.ref.get('utm_zone', '1N')  # Default to 1N if not available
		if 'N' in utm_zone:
			epsg_code = 32600 + int(utm_zone.replace('N', ''))
		elif 'S' in utm_zone:
			epsg_code = 32700 + int(utm_zone.replace('S', ''))
		else:
			epsg_code = 32601  # Default to UTM 1N
		
		crs = CRS.from_epsg(epsg_code)
		transform = from_bounds(*bounds, self.ref['z_grid'].shape[1], self.ref['z_grid'].shape[0])
		
		# Define the surfaces to export
		surfaces_to_export = []
		
		# Final Surface (masked depth)
		if 'z_final_grid' in self.ref and self.ref['z_final_grid'] is not None:
			surfaces_to_export.append(('final_surface', self.ref['z_final_grid'], 'Final Surface (m)'))
		
		# Depth
		if 'z_grid' in self.ref and self.ref['z_grid'] is not None:
			surfaces_to_export.append(('depth', self.ref['z_grid'], 'Depth (m)'))
		
		# Uncertainty
		if 'u_grid' in self.ref and self.ref['u_grid'] is not None:
			surfaces_to_export.append(('uncertainty', self.ref['u_grid'], 'Uncertainty (m)'))
		
		# Density
		if 'c_grid' in self.ref and self.ref['c_grid'] is not None:
			surfaces_to_export.append(('density', self.ref['c_grid'], 'Density (count)'))
		
		# Slope
		if 's_grid' in self.ref and self.ref['s_grid'] is not None:
			surfaces_to_export.append(('slope', self.ref['s_grid'], 'Slope (degrees)'))
		
		if not surfaces_to_export:
			update_log(self, 'ERROR: No valid surfaces available for export', font_color="red")
			return
		
		# Export each surface
		exported_files = []
		for surface_name, data_grid, description in surfaces_to_export:
			try:
				# Prepare data for export (replace NaN with nodata value)
				export_data = data_grid.copy()
				nodata_value = -9999
				export_data[np.isnan(export_data)] = nodata_value
				
				# Create filename
				filename = f"{export_name}_{surface_name}.tif"
				filepath = os.path.join(output_dir, filename)
				
				# Write GeoTIFF
				with rasterio.open(
					filepath,
					'w',
					driver='GTiff',
					height=export_data.shape[0],
					width=export_data.shape[1],
					count=1,
					dtype=export_data.dtype,
					crs=crs,
					transform=transform,
					nodata=nodata_value
				) as dst:
					dst.write(export_data, 1)
					dst.update_tags(TIFFTAG_DATETIME=description)
				
				exported_files.append(filename)
				update_log(self, f'Exported: {filename}')
				
			except Exception as e:
				update_log(self, f'ERROR: Failed to export {surface_name}: {str(e)}', font_color="red")
				print(f"Error exporting {surface_name}: {str(e)}")
		
		if exported_files:
			update_log(self, f'Successfully exported {len(exported_files)} GeoTIFF files to: {output_dir}')
			
			# Show success message
			from PyQt6.QtWidgets import QMessageBox
			file_list = "\n".join([f"  • {f}" for f in exported_files])
			msg = f"Successfully exported {len(exported_files)} GeoTIFF files:\n{file_list}\n\nLocation: {output_dir}"
			QMessageBox.information(self, 'Export Complete', msg)
		else:
			update_log(self, 'ERROR: No files were successfully exported', font_color="red")
		
	except Exception as e:
		update_log(self, f'ERROR: Export failed: {str(e)}', font_color="red")
		print(f"Error in export_all_to_geotiff: {str(e)}")

def save_last_directory(directory):
	"""Save the last used directory to configuration file"""
	config_file = get_config_file_path()
	config = {}
	
	# Load existing config if it exists
	if os.path.exists(config_file):
		try:
			with open(config_file, 'r') as f:
				config = json.load(f)
		except (json.JSONDecodeError, IOError):
			config = {}
	
	# Update with new directory
	config['last_output_directory'] = directory
	
	# Save config
	try:
		with open(config_file, 'w') as f:
			json.dump(config, f, indent=2)
	except IOError:
		print(f"Warning: Could not save configuration to {config_file}")

def load_last_directory():
	"""Load the last used directory from configuration file"""
	config_file = get_config_file_path()
	
	if not os.path.exists(config_file):
		return None
	
	try:
		with open(config_file, 'r') as f:
			config = json.load(f)
		return config.get('last_output_directory', None)
	except (json.JSONDecodeError, IOError):
		return None

def load_session_config():
	"""Load session configuration including last used directories"""
	config_file = os.path.join(os.path.expanduser("~"), ".swath_accuracy_session.json")
	
	default_config = {
		"last_ref_surface_dir": os.getcwd(),
		"last_density_surface_dir": os.getcwd(), 
		"last_crossline_dir": os.getcwd(),
		"last_output_dir": os.getcwd(),
		"last_plot_basename": "swath_accuracy_plots",
		"last_plot_directory": os.getcwd(),
		"last_plot_extension": ".PNG",
		"last_session_dir": os.getcwd(),
		"last_xyz_dir": os.getcwd(),
		"last_xyd_dir": os.getcwd(),
		"last_tide_dir": os.getcwd()
	}
	
	try:
		if os.path.exists(config_file):
			with open(config_file, 'r') as f:
				config = json.load(f)
				# Update with any missing keys from default
				for key, value in default_config.items():
					if key not in config:
						config[key] = value
				print(f"DEBUG: Loaded session config: {config}")
				return config
		else:
			print(f"DEBUG: No session config file found, using defaults")
			return default_config
	except Exception as e:
		print(f"Warning: Could not load session config: {e}")
		return default_config

def save_session_config(config):
	"""Save session configuration including last used directories"""
	config_file = os.path.join(os.path.expanduser("~"), ".swath_accuracy_session.json")
	
	try:
		with open(config_file, 'w') as f:
			json.dump(config, f, indent=2)
	except Exception as e:
		print(f"Warning: Could not save session config: {e}")

def update_last_directory(config_key, directory):
	"""Update the last used directory for a specific file type"""
	print(f"DEBUG: update_last_directory called with config_key={config_key}, directory={directory}")
	if directory and os.path.exists(directory):
		config = load_session_config()
		config[config_key] = directory
		save_session_config(config)
		print(f"DEBUG: Saved {config_key} = {directory}")
	else:
		print(f"DEBUG: Directory not saved - directory={directory}, exists={os.path.exists(directory) if directory else 'None'}")


def save_current_filters(self):
	"""Save current filter settings to a parameter file"""
	try:
		# Create filter settings dictionary
		filter_settings = {
			# Reference surface filters
			'ref_depth_min': self.min_depth_ref_tb.text(),
			'ref_depth_max': self.max_depth_ref_tb.text(),
			'ref_slope_max': self.max_slope_tb.text(),
			'ref_density_min': self.min_dens_tb.text(),
			'ref_uncertainty_max': self.max_u_tb.text(),
			'ref_slope_window': self.slope_win_cbox.currentText(),
			
			# Crossline filters
			'xline_depth_min': self.min_depth_xline_tb.text(),
			'xline_depth_max': self.max_depth_xline_tb.text(),
			'xline_angle_min': self.min_angle_xline_tb.text(),
			'xline_angle_max': self.max_angle_xline_tb.text(),
			'xline_bs_min': self.min_bs_xline_tb.text(),
			'xline_bs_max': self.max_bs_xline_tb.text(),
			'xline_dz_max': self.max_dz_tb.text(),
			'xline_dz_wd_max': self.max_dz_wd_tb.text(),
			'xline_bin_count_min': self.min_bin_count_tb.text(),
			'xline_depth_mode': self.depth_mode_cbox.currentText(),
			
			# Filter enable/disable states
			'ref_depth_enabled': self.depth_ref_gb.isChecked(),
			'ref_slope_enabled': self.slope_gb.isChecked(),
			'ref_density_enabled': self.density_gb.isChecked(),
			'ref_uncertainty_enabled': self.uncertainty_gb.isChecked(),
			'xline_depth_enabled': self.depth_xline_gb.isChecked(),
			'xline_angle_enabled': self.angle_xline_gb.isChecked(),
			'xline_bs_enabled': self.bs_xline_gb.isChecked(),
			'xline_dz_abs_enabled': hasattr(self, 'dz_abs_gb') and self.dz_abs_gb.isChecked(),
			'xline_dz_pct_enabled': hasattr(self, 'dz_pct_gb') and self.dz_pct_gb.isChecked(),
			'xline_depth_mode_enabled': self.depth_mode_gb.isChecked(),
			'xline_bin_count_enabled': self.bin_count_gb.isChecked(),
			
			# Point count settings
			'pt_count_enabled': self.pt_count_gb.isChecked(),
			'pt_count_max': self.max_count_tb.text(),
			'pt_count_dec_fac': self.dec_fac_tb.text(),
			
			# Bin decimation settings
			'bin_decimation_enabled': hasattr(self, 'bin_decimation_gb') and self.bin_decimation_gb.isChecked(),
			'bin_decimation_max_points': hasattr(self, 'max_points_per_bin_tb') and self.max_points_per_bin_tb.text() or '1000',
			
			# Timestamp
			'saved_timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
		}
		
		# Get config file path
		config_path = get_config_file_path()
		filter_file = os.path.join(os.path.dirname(config_path), 'last_filters.json')
		
		# Save to JSON file
		with open(filter_file, 'w') as f:
			json.dump(filter_settings, f, indent=2)
		
		update_log(self, f'Filter settings saved to: {os.path.basename(filter_file)}')
		
	except Exception as e:
		update_log(self, f'Error saving filter settings: {str(e)}')
		print(f"Error saving filter settings: {str(e)}")


def load_last_filters(self):
	"""Load the last saved filter settings"""
	try:
		# Get config file path
		config_path = get_config_file_path()
		filter_file = os.path.join(os.path.dirname(config_path), 'last_filters.json')
		
		if not os.path.exists(filter_file):
			update_log(self, 'No saved filter settings found')
			return
		
		# Load from JSON file
		with open(filter_file, 'r') as f:
			filter_settings = json.load(f)
		
		# Store current filter states to detect changes
		old_ref_filters = {
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
		
		old_xline_filters = {
			'depth_enabled': self.depth_xline_gb.isChecked(),
			'depth_min': self.min_depth_xline_tb.text(),
			'depth_max': self.max_depth_xline_tb.text(),
			'angle_enabled': self.angle_xline_gb.isChecked(),
			'angle_min': self.min_angle_xline_tb.text(),
			'angle_max': self.max_angle_xline_tb.text(),
			'bs_enabled': self.bs_xline_gb.isChecked(),
			'bs_min': self.min_bs_xline_tb.text(),
			'bs_max': self.max_bs_xline_tb.text(),
			'dz_abs_enabled': hasattr(self, 'dz_abs_gb') and self.dz_abs_gb.isChecked(),
			'dz_max': self.max_dz_tb.text(),
			'dz_pct_enabled': hasattr(self, 'dz_pct_gb') and self.dz_pct_gb.isChecked(),
			'dz_wd_max': self.max_dz_wd_tb.text(),
			'depth_mode_enabled': self.depth_mode_gb.isChecked(),
			'depth_mode': self.depth_mode_cbox.currentText(),
			'bin_count_enabled': self.bin_count_gb.isChecked(),
			'bin_count_min': self.min_bin_count_tb.text()
		}
		
		# Apply reference surface filters
		self.min_depth_ref_tb.setText(filter_settings.get('ref_depth_min', '0'))
		self.max_depth_ref_tb.setText(filter_settings.get('ref_depth_max', '10000'))
		self.max_slope_tb.setText(filter_settings.get('ref_slope_max', '5'))
		self.min_dens_tb.setText(filter_settings.get('ref_density_min', '10'))
		self.max_u_tb.setText(filter_settings.get('ref_uncertainty_max', '10'))
		
		# Set slope window
		slope_window = filter_settings.get('ref_slope_window', '3x3')
		slope_idx = self.slope_win_cbox.findText(slope_window)
		if slope_idx >= 0:
			self.slope_win_cbox.setCurrentIndex(slope_idx)
		
		# Apply crossline filters
		self.min_depth_xline_tb.setText(filter_settings.get('xline_depth_min', '0'))
		self.max_depth_xline_tb.setText(filter_settings.get('xline_depth_max', '10000'))
		self.min_angle_xline_tb.setText(filter_settings.get('xline_angle_min', '-75'))
		self.max_angle_xline_tb.setText(filter_settings.get('xline_angle_max', '75'))
		self.min_bs_xline_tb.setText(filter_settings.get('xline_bs_min', '-50'))
		self.max_bs_xline_tb.setText(filter_settings.get('xline_bs_max', '0'))
		self.max_dz_tb.setText(filter_settings.get('xline_dz_max', '10'))
		self.max_dz_wd_tb.setText(filter_settings.get('xline_dz_wd_max', '5'))
		self.min_bin_count_tb.setText(filter_settings.get('xline_bin_count_min', '10'))
		
		# Set depth mode
		depth_mode = filter_settings.get('xline_depth_mode', 'None')
		depth_mode_idx = self.depth_mode_cbox.findText(depth_mode)
		if depth_mode_idx >= 0:
			self.depth_mode_cbox.setCurrentIndex(depth_mode_idx)
		
		# Apply point count settings
		self.max_count_tb.setText(filter_settings.get('pt_count_max', '50000'))
		self.dec_fac_tb.setText(filter_settings.get('pt_count_dec_fac', '1'))
		
		# Apply bin decimation settings
		if hasattr(self, 'max_points_per_bin_tb'):
			self.max_points_per_bin_tb.setText(filter_settings.get('bin_decimation_max_points', '350'))
		
		# Set enable/disable states
		self.depth_ref_gb.setChecked(filter_settings.get('ref_depth_enabled', True))
		self.slope_gb.setChecked(filter_settings.get('ref_slope_enabled', True))
		self.density_gb.setChecked(filter_settings.get('ref_density_enabled', True))
		self.uncertainty_gb.setChecked(filter_settings.get('ref_uncertainty_enabled', True))
		self.depth_xline_gb.setChecked(filter_settings.get('xline_depth_enabled', True))
		self.angle_xline_gb.setChecked(filter_settings.get('xline_angle_enabled', True))
		self.bs_xline_gb.setChecked(filter_settings.get('xline_bs_enabled', True))
		if hasattr(self, 'dz_abs_gb'):
			self.dz_abs_gb.setChecked(filter_settings.get('xline_dz_abs_enabled', True))
		if hasattr(self, 'dz_pct_gb'):
			self.dz_pct_gb.setChecked(filter_settings.get('xline_dz_pct_enabled', True))
		self.depth_mode_gb.setChecked(filter_settings.get('xline_depth_mode_enabled', True))
		self.bin_count_gb.setChecked(filter_settings.get('xline_bin_count_enabled', True))
		self.pt_count_gb.setChecked(filter_settings.get('pt_count_enabled', True))
		
		# Set bin decimation enable/disable state
		if hasattr(self, 'bin_decimation_gb'):
			self.bin_decimation_gb.setChecked(filter_settings.get('bin_decimation_enabled', False))
		
		# Get timestamp
		timestamp = filter_settings.get('saved_timestamp', 'Unknown')
		update_log(self, f'Filter settings loaded (saved: {timestamp})')
		
		# Force UI to process events so widget values are updated
		from PyQt6 import QtWidgets
		QtWidgets.QApplication.processEvents()
		
		# Determine which filters changed and update appropriate data
		ref_filters_changed = self._detect_ref_filter_changes(old_ref_filters)
		xline_filters_changed = self._detect_xline_filter_changes(old_xline_filters)
		
		# Update appropriate data based on what changed
		if ref_filters_changed:
			update_log(self, 'Reference surface filters changed - updating reference surface data')
			calc_accuracy(self, recalc_ref_only=True)  # Update reference surface only
		elif xline_filters_changed:
			update_log(self, 'Crossline filters changed - updating crossline data')
			calc_accuracy(self, recalc_bins_only=True)  # Only binning needs update
		else:
			update_log(self, 'No filter changes detected - no recalculation needed')
		
	except Exception as e:
		update_log(self, f'Error loading filter settings: {str(e)}')
		print(f"Error loading filter settings: {str(e)}")


def load_default_filters(self):
	"""Reset all filters to their default values"""
	try:
		# Store current filter states to detect changes
		old_ref_filters = {
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
		
		old_xline_filters = {
			'depth_enabled': self.depth_xline_gb.isChecked(),
			'depth_min': self.min_depth_xline_tb.text(),
			'depth_max': self.max_depth_xline_tb.text(),
			'angle_enabled': self.angle_xline_gb.isChecked(),
			'angle_min': self.min_angle_xline_tb.text(),
			'angle_max': self.max_angle_xline_tb.text(),
			'bs_enabled': self.bs_xline_gb.isChecked(),
			'bs_min': self.min_bs_xline_tb.text(),
			'bs_max': self.max_bs_xline_tb.text(),
			'dz_abs_enabled': hasattr(self, 'dz_abs_gb') and self.dz_abs_gb.isChecked(),
			'dz_max': self.max_dz_tb.text(),
			'dz_pct_enabled': hasattr(self, 'dz_pct_gb') and self.dz_pct_gb.isChecked(),
			'dz_wd_max': self.max_dz_wd_tb.text(),
			'depth_mode_enabled': self.depth_mode_gb.isChecked(),
			'depth_mode': self.depth_mode_cbox.currentText(),
			'bin_count_enabled': self.bin_count_gb.isChecked(),
			'bin_count_min': self.min_bin_count_tb.text()
		}
		
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
		
		# Bin decimation defaults (use 350 as default, matching the decimation function)
		if hasattr(self, 'max_points_per_bin_tb'):
			self.max_points_per_bin_tb.setText('350')
		
		# Disable all filters by default (except bin count and decimation settings)
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
		self.bin_count_gb.setChecked(True)  # Keep bin count enabled by default
		
		# Use bin-based decimation by default (point count limit disabled)
		self.pt_count_gb.setChecked(False)
		
		# Enable bin decimation by default
		if hasattr(self, 'bin_decimation_gb'):
			self.bin_decimation_gb.setChecked(True)
		
		update_log(self, 'All filters reset to default values')
		
		# Force UI to process events so widget values are updated
		from PyQt6 import QtWidgets
		QtWidgets.QApplication.processEvents()
		
		# Determine which filters changed and update appropriate data
		ref_filters_changed = self._detect_ref_filter_changes(old_ref_filters)
		xline_filters_changed = self._detect_xline_filter_changes(old_xline_filters)
		
		# Update appropriate data based on what changed
		if ref_filters_changed:
			update_log(self, 'Reference surface filters changed - updating reference surface data')
			calc_accuracy(self, recalc_ref_only=True)  # Update reference surface only
		elif xline_filters_changed:
			update_log(self, 'Crossline filters changed - updating crossline data')
			calc_accuracy(self, recalc_bins_only=True)  # Only binning needs update
		else:
			update_log(self, 'No filter changes detected - no recalculation needed')
		
	except Exception as e:
		update_log(self, f'Error loading default filters: {str(e)}')
		print(f"Error loading default filters: {str(e)}")


def get_session_default_directory(self):
	"""Get the best default directory for session files - either last session directory or last plot directory"""
	config = load_session_config()
	last_session_dir = config.get("last_session_dir", "")
	last_plot_dir = config.get("last_plot_directory", "")
	
	# Prefer the plot directory if it exists and is more recent, otherwise use session directory
	if last_plot_dir and os.path.exists(last_plot_dir):
		return last_plot_dir
	elif last_session_dir and os.path.exists(last_session_dir):
		return last_session_dir
	else:
		return os.getcwd()

def save_session(self):
	"""Save complete session including files, plot settings, and filter settings"""
	try:
		from PyQt6.QtWidgets import QFileDialog
		import json
		from datetime import datetime
		
		# Get save file path from user
		default_dir = get_session_default_directory(self)
		initial_path = os.path.join(default_dir, "analysis_session.json")
		file_path, _ = QFileDialog.getSaveFileName(
			self, 
			"Save Session", 
			initial_path, 
			"Session Files (*.json);;All Files (*)"
		)
		
		# Update last used directory if file was selected
		if file_path:
			directory = os.path.dirname(file_path)
			update_last_directory("last_session_dir", directory)
		
		if not file_path:
			return  # User cancelled
		
		# Get current file list
		from file_fun import get_current_file_list
		get_current_file_list(self)
		
		# Collect all session data
		session_data = {
			'saved_timestamp': datetime.now().isoformat(),
			'version': '1.0',
			
			# File paths
			'files': {
				'crossline_files': getattr(self, 'filenames', []),
				'reference_surface_file': getattr(self, 'ref_surf_file', ''),
				'density_surface_file': getattr(self, 'dens_surf_file', ''),
				'tide_file': getattr(self, 'tide_file', ''),
				'output_directory': self.output_dir
			},
			
			# Plot settings (from Plot tab)
			'plot_settings': {
				# Custom info
				'custom_info_enabled': self.custom_info_gb.isChecked(),
				'model': self.model_cbox.currentText(),
				'ship_name': self.ship_tb.text(),
				'cruise_name': self.cruise_tb.text(),
				'show_model': self.show_model_chk.isChecked(),
				'show_ship': self.show_ship_chk.isChecked(),
				'show_cruise': self.show_cruise_chk.isChecked(),
				
				# Data reference
				'reference_data': self.ref_cbox.currentText(),
				'tide_units': self.tide_unit_cbox.currentText(),
				'waterline_adjustment': self.waterline_tb.text(),
				
				# Point style
				'point_size': self.pt_size_cbox.currentText(),
				'point_size_coverage': self.pt_size_cov_cbox.currentText(),
				'point_opacity_accuracy': self.pt_alpha_acc_tb.text(),
				'point_opacity_coverage': self.pt_alpha_cov_tb.text(),
				
				# Plot limits
				'custom_plot_limits': self.plot_lim_gb.isChecked(),
				'max_beam_angle': self.max_beam_angle_tb.text(),
				'angle_spacing': self.angle_spacing_tb.text(),
				'max_bias': self.max_bias_tb.text(),
				'max_std': self.max_std_tb.text(),
				'axis_margin': self.axis_margin_tb.text(),
				
				# Toggle options
				'show_acc_proc_text': self.show_acc_proc_text_chk.isChecked(),
				'show_ref_proc_text': self.show_ref_proc_text_chk.isChecked(),
				'show_grid_lines': self.grid_lines_toggle_chk.isChecked(),
				'show_IHO_lines': self.IHO_lines_toggle_chk.isChecked(),
				'update_ref_plots': self.update_ref_plots_chk.isChecked(),
				'show_xline_coverage': self.show_xline_cov_chk.isChecked(),
				'show_uncertainty_plot': self.show_u_plot_chk.isChecked(),
				'show_shaded_relief': self.show_shaded_relief_chk.isChecked(),
				'show_special_order': self.show_special_order_chk.isChecked(),
				'show_order_1a': self.show_order_1a_chk.isChecked(),
				'show_order_1b': self.show_order_1b_chk.isChecked(),
				'show_order_2': self.show_order_2_chk.isChecked(),
				'show_order_3': self.show_order_3_chk.isChecked(),
				
				# Flatten swath
				'flatten_mean_enabled': self.flatten_mean_gb.isChecked(),
				'force_zero_mean': self.set_zero_mean_chk.isChecked(),
				'mean_bias_angle_limit': self.mean_bias_angle_lim_tb.text()
			},
			
			# Filter settings (from Filter tab)
			'filter_settings': {
				# Reference surface filters
				'ref_depth_enabled': self.depth_ref_gb.isChecked(),
				'ref_depth_min': self.min_depth_ref_tb.text(),
				'ref_depth_max': self.max_depth_ref_tb.text(),
				'ref_slope_enabled': self.slope_gb.isChecked(),
				'ref_slope_max': self.max_slope_tb.text(),
				'ref_slope_window': self.slope_win_cbox.currentText(),
				'ref_density_enabled': self.density_gb.isChecked(),
				'ref_density_min': self.min_dens_tb.text(),
				'ref_uncertainty_enabled': self.uncertainty_gb.isChecked(),
				'ref_uncertainty_max': self.max_u_tb.text(),
				
				# Crossline filters
				'xline_depth_enabled': self.depth_xline_gb.isChecked(),
				'xline_depth_min': self.min_depth_xline_tb.text(),
				'xline_depth_max': self.max_depth_xline_tb.text(),
				'xline_angle_enabled': self.angle_xline_gb.isChecked(),
				'xline_angle_min': self.min_angle_xline_tb.text(),
				'xline_angle_max': self.max_angle_xline_tb.text(),
				'xline_bs_enabled': self.bs_xline_gb.isChecked(),
				'xline_bs_min': self.min_bs_xline_tb.text(),
				'xline_bs_max': self.max_bs_xline_tb.text(),
				'xline_dz_abs_enabled': hasattr(self, 'dz_abs_gb') and self.dz_abs_gb.isChecked(),
				'xline_dz_max': self.max_dz_tb.text(),
				'xline_dz_pct_enabled': hasattr(self, 'dz_pct_gb') and self.dz_pct_gb.isChecked(),
				'xline_dz_wd_max': self.max_dz_wd_tb.text(),
				'xline_depth_mode_enabled': self.depth_mode_gb.isChecked(),
				'xline_depth_mode': self.depth_mode_cbox.currentText(),
				'xline_bin_count_enabled': self.bin_count_gb.isChecked(),
				'xline_bin_count_min': self.min_bin_count_tb.text(),
				
				# Point count settings
				'pt_count_enabled': self.pt_count_gb.isChecked(),
				'pt_count_max': self.max_count_tb.text(),
				'pt_count_dec_fac': self.dec_fac_tb.text(),
				
				# Bin decimation settings
				'bin_decimation_enabled': hasattr(self, 'bin_decimation_gb') and self.bin_decimation_gb.isChecked(),
				'bin_decimation_max_points': hasattr(self, 'max_points_per_bin_tb') and self.max_points_per_bin_tb.text() or '350'
			}
		}
		
		# Save to file
		with open(file_path, 'w') as f:
			json.dump(session_data, f, indent=2)
		
		update_log(self, f'Session saved to: {file_path}')
		
	except Exception as e:
		update_log(self, f'Error saving session: {str(e)}')
		print(f"Error saving session: {str(e)}")


def load_session(self):
	"""Load complete session including files, plot settings, and filter settings"""
	try:
		from PyQt6.QtWidgets import QFileDialog
		import json
		import os
		
		# Get load file path from user
		default_dir = get_session_default_directory(self)
		file_path, _ = QFileDialog.getOpenFileName(
			self, 
			"Load Session", 
			default_dir, 
			"Session Files (*.json);;All Files (*)"
		)
		
		# Update last used directory if file was selected
		if file_path:
			directory = os.path.dirname(file_path)
			update_last_directory("last_session_dir", directory)
		
		if not file_path:
			return  # User cancelled
		
		# Load session data
		with open(file_path, 'r') as f:
			session_data = json.load(f)
		
		# Clear current files
		clear_files(self)
		
		# Load files
		files = session_data.get('files', {})
		
		# Load crossline files
		crossline_files = files.get('crossline_files', [])
		for file_path in crossline_files:
			if os.path.exists(file_path):
				self.file_list.add_file(file_path)
			else:
				update_log(self, f'Warning: Crossline file not found: {file_path}')
		
		# Load reference surface file
		ref_file = files.get('reference_surface_file', '')
		if ref_file and os.path.exists(ref_file):
			self.ref_surf_file = ref_file
			parse_ref_depth(self)
		elif ref_file:
			update_log(self, f'Warning: Reference surface file not found: {ref_file}')
		
		# Load density surface file
		dens_file = files.get('density_surface_file', '')
		if dens_file and os.path.exists(dens_file):
			self.dens_surf_file = dens_file
			parse_ref_dens(self)
		elif dens_file:
			update_log(self, f'Warning: Density surface file not found: {dens_file}')
		
		# Load tide file
		tide_file = files.get('tide_file', '')
		if tide_file and os.path.exists(tide_file):
			self.tide_file = tide_file
			process_tide(self)
		elif tide_file:
			update_log(self, f'Warning: Tide file not found: {tide_file}')
		
		# Set output directory
		output_dir = files.get('output_directory', '')
		if output_dir and os.path.exists(output_dir):
			self.output_dir = output_dir
			self.current_outdir_lbl.setText('Current output directory:\n' + self.output_dir)
		
		# Load plot settings
		plot_settings = session_data.get('plot_settings', {})
		
		# Custom info
		if plot_settings.get('custom_info_enabled', False):
			self.custom_info_gb.setChecked(True)
		self.model_cbox.setCurrentText(plot_settings.get('model', ''))
		self.ship_tb.setText(plot_settings.get('ship_name', ''))
		self.cruise_tb.setText(plot_settings.get('cruise_name', ''))
		self.show_model_chk.setChecked(plot_settings.get('show_model', True))
		self.show_ship_chk.setChecked(plot_settings.get('show_ship', True))
		self.show_cruise_chk.setChecked(plot_settings.get('show_cruise', True))
		
		# Data reference
		self.ref_cbox.setCurrentText(plot_settings.get('reference_data', ''))
		self.tide_unit_cbox.setCurrentText(plot_settings.get('tide_units', ''))
		self.waterline_tb.setText(plot_settings.get('waterline_adjustment', '0.00'))
		
		# Point style
		self.pt_size_cbox.setCurrentText(plot_settings.get('point_size', '1'))
		self.pt_size_cov_cbox.setCurrentText(plot_settings.get('point_size_coverage', '5'))
		# Handle both old and new opacity settings for backward compatibility
		if 'point_opacity_accuracy' in plot_settings:
			self.pt_alpha_acc_tb.setText(plot_settings.get('point_opacity_accuracy', '100'))
			self.pt_alpha_cov_tb.setText(plot_settings.get('point_opacity_coverage', '6'))
		else:
			# Old format: single opacity value applies to both
			old_opacity = plot_settings.get('point_opacity', '100')
			self.pt_alpha_acc_tb.setText(old_opacity)
			self.pt_alpha_cov_tb.setText(old_opacity)
		
		# Plot limits
		self.plot_lim_gb.setChecked(plot_settings.get('custom_plot_limits', False))
		self.max_beam_angle_tb.setText(plot_settings.get('max_beam_angle', ''))
		self.angle_spacing_tb.setText(plot_settings.get('angle_spacing', ''))
		self.max_bias_tb.setText(plot_settings.get('max_bias', ''))
		self.max_std_tb.setText(plot_settings.get('max_std', ''))
		self.axis_margin_tb.setText(plot_settings.get('axis_margin', ''))
		
		# Toggle options
		self.show_acc_proc_text_chk.setChecked(plot_settings.get('show_acc_proc_text', False))
		self.show_ref_proc_text_chk.setChecked(plot_settings.get('show_ref_proc_text', False))
		self.grid_lines_toggle_chk.setChecked(plot_settings.get('show_grid_lines', True))
		self.IHO_lines_toggle_chk.setChecked(plot_settings.get('show_IHO_lines', True))
		self.update_ref_plots_chk.setChecked(plot_settings.get('update_ref_plots', True))
		self.show_xline_cov_chk.setChecked(plot_settings.get('show_xline_coverage', True))
		self.show_u_plot_chk.setChecked(plot_settings.get('show_uncertainty_plot', True))
		self.show_shaded_relief_chk.setChecked(plot_settings.get('show_shaded_relief', True))
		self.show_special_order_chk.setChecked(plot_settings.get('show_special_order', False))
		self.show_order_1a_chk.setChecked(plot_settings.get('show_order_1a', False))
		self.show_order_1b_chk.setChecked(plot_settings.get('show_order_1b', False))
		self.show_order_2_chk.setChecked(plot_settings.get('show_order_2', False))
		self.show_order_3_chk.setChecked(plot_settings.get('show_order_3', False))
		
		# Flatten swath
		self.flatten_mean_gb.setChecked(plot_settings.get('flatten_mean_enabled', False))
		self.set_zero_mean_chk.setChecked(plot_settings.get('force_zero_mean', False))
		self.mean_bias_angle_lim_tb.setText(plot_settings.get('mean_bias_angle_limit', '40'))
		
		# Load filter settings
		filter_settings = session_data.get('filter_settings', {})
		
		# Reference surface filters
		self.depth_ref_gb.setChecked(filter_settings.get('ref_depth_enabled', False))
		self.min_depth_ref_tb.setText(filter_settings.get('ref_depth_min', '0'))
		self.max_depth_ref_tb.setText(filter_settings.get('ref_depth_max', '10000'))
		self.slope_gb.setChecked(filter_settings.get('ref_slope_enabled', False))
		self.max_slope_tb.setText(filter_settings.get('ref_slope_max', '5'))
		slope_window = filter_settings.get('ref_slope_window', '3x3')
		slope_idx = self.slope_win_cbox.findText(slope_window)
		if slope_idx >= 0:
			self.slope_win_cbox.setCurrentIndex(slope_idx)
		self.density_gb.setChecked(filter_settings.get('ref_density_enabled', False))
		self.min_dens_tb.setText(filter_settings.get('ref_density_min', '10'))
		self.uncertainty_gb.setChecked(filter_settings.get('ref_uncertainty_enabled', False))
		self.max_u_tb.setText(filter_settings.get('ref_uncertainty_max', '10'))
		
		# Crossline filters
		self.depth_xline_gb.setChecked(filter_settings.get('xline_depth_enabled', False))
		self.min_depth_xline_tb.setText(filter_settings.get('xline_depth_min', '0'))
		self.max_depth_xline_tb.setText(filter_settings.get('xline_depth_max', '10000'))
		self.angle_xline_gb.setChecked(filter_settings.get('xline_angle_enabled', False))
		self.min_angle_xline_tb.setText(filter_settings.get('xline_angle_min', '-75'))
		self.max_angle_xline_tb.setText(filter_settings.get('xline_angle_max', '75'))
		self.bs_xline_gb.setChecked(filter_settings.get('xline_bs_enabled', False))
		self.min_bs_xline_tb.setText(filter_settings.get('xline_bs_min', '-50'))
		self.max_bs_xline_tb.setText(filter_settings.get('xline_bs_max', '0'))
		if hasattr(self, 'dz_abs_gb'):
			self.dz_abs_gb.setChecked(filter_settings.get('xline_dz_abs_enabled', False))
		self.max_dz_tb.setText(filter_settings.get('xline_dz_max', '10'))
		if hasattr(self, 'dz_pct_gb'):
			self.dz_pct_gb.setChecked(filter_settings.get('xline_dz_pct_enabled', False))
		self.max_dz_wd_tb.setText(filter_settings.get('xline_dz_wd_max', '5'))
		self.depth_mode_gb.setChecked(filter_settings.get('xline_depth_mode_enabled', False))
		depth_mode = filter_settings.get('xline_depth_mode', 'None')
		depth_mode_idx = self.depth_mode_cbox.findText(depth_mode)
		if depth_mode_idx >= 0:
			self.depth_mode_cbox.setCurrentIndex(depth_mode_idx)
		self.bin_count_gb.setChecked(filter_settings.get('xline_bin_count_enabled', True))
		self.min_bin_count_tb.setText(filter_settings.get('xline_bin_count_min', '10'))
		
		# Point count settings
		self.pt_count_gb.setChecked(filter_settings.get('pt_count_enabled', False))
		self.max_count_tb.setText(filter_settings.get('pt_count_max', '50000'))
		self.dec_fac_tb.setText(filter_settings.get('pt_count_dec_fac', '1'))
		
		# Bin decimation settings
		if hasattr(self, 'bin_decimation_gb'):
			self.bin_decimation_gb.setChecked(filter_settings.get('bin_decimation_enabled', True))
		if hasattr(self, 'max_points_per_bin_tb'):
			self.max_points_per_bin_tb.setText(filter_settings.get('bin_decimation_max_points', '350'))
		
		# Force UI to process events
		from PyQt6 import QtWidgets
		QtWidgets.QApplication.processEvents()
		
		# Update buttons and recalculate if needed
		update_buttons(self)
		
		# Recalculate accuracy if files are loaded
		if self.file_list.count() > 0:
			calc_accuracy(self)
		
		# Refresh plots
		refresh_plot(self)
		
		# Get timestamp
		timestamp = session_data.get('saved_timestamp', 'Unknown')
		update_log(self, f'Session loaded (saved: {timestamp})')
		
	except Exception as e:
		update_log(self, f'Error loading session: {str(e)}')
		print(f"Error loading session: {str(e)}")


def save_current_plot_settings(self):
	"""Save current plot settings to a parameter file"""
	try:
		# Create plot settings dictionary
		plot_settings = {
			# Custom info settings
			'custom_info_enabled': self.custom_info_gb.isChecked(),
			'model': self.model_cbox.currentText(),
			'ship_name': self.ship_tb.text(),
			'cruise_name': self.cruise_tb.text(),
			'show_model': self.show_model_chk.isChecked(),
			'show_ship': self.show_ship_chk.isChecked(),
			'show_cruise': self.show_cruise_chk.isChecked(),
			
			# Data reference settings
			'reference_data': self.ref_cbox.currentText(),
			'tide_units': self.tide_unit_cbox.currentText(),
			'waterline_adjustment': self.waterline_tb.text(),
			
			# Point style settings
			'point_size': self.pt_size_cbox.currentText(),
			'point_size_coverage': self.pt_size_cov_cbox.currentText(),
			'point_opacity_accuracy': self.pt_alpha_acc_tb.text(),
			'point_opacity_coverage': self.pt_alpha_cov_tb.text(),
			
			# Plot limits settings
			'custom_plot_limits': self.plot_lim_gb.isChecked(),
			'max_beam_angle': self.max_beam_angle_tb.text(),
			'angle_spacing': self.angle_spacing_tb.text(),
			'max_bias': self.max_bias_tb.text(),
			'max_std': self.max_std_tb.text(),
			'axis_margin': self.axis_margin_tb.text(),
			'tide_range': self.tide_range_tb.text(),
			
			# Toggle options
			'show_acc_proc_text': self.show_acc_proc_text_chk.isChecked(),
			'show_ref_proc_text': self.show_ref_proc_text_chk.isChecked(),
			'show_grid_lines': self.grid_lines_toggle_chk.isChecked(),
			'show_IHO_lines': self.IHO_lines_toggle_chk.isChecked(),
			'update_ref_plots': self.update_ref_plots_chk.isChecked(),
			'show_xline_coverage': self.show_xline_cov_chk.isChecked(),
			'show_uncertainty_plot': self.show_u_plot_chk.isChecked(),
			'show_shaded_relief': self.show_shaded_relief_chk.isChecked(),
			'show_special_order': self.show_special_order_chk.isChecked(),
			'show_order_1a': self.show_order_1a_chk.isChecked(),
			'show_order_1b': self.show_order_1b_chk.isChecked(),
			'show_order_2': self.show_order_2_chk.isChecked(),
			'show_order_3': self.show_order_3_chk.isChecked(),
			
			# Flatten swath settings
			'flatten_mean_enabled': self.flatten_mean_gb.isChecked(),
			'force_zero_mean': self.set_zero_mean_chk.isChecked(),
			'mean_bias_angle_limit': self.mean_bias_angle_lim_tb.text(),
			
			# Timestamp
			'saved_timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
		}
		
		# Get config file path
		config_path = get_config_file_path()
		plot_file = os.path.join(os.path.dirname(config_path), 'last_plot_settings.json')
		
		# Save to JSON file
		with open(plot_file, 'w') as f:
			json.dump(plot_settings, f, indent=2)
		
		update_log(self, f'Plot settings saved to: {os.path.basename(plot_file)}')
		
	except Exception as e:
		update_log(self, f'Error saving plot settings: {str(e)}')
		print(f"Error saving plot settings: {str(e)}")


def load_last_plot_settings(self):
	"""Load the last saved plot settings"""
	try:
		# Get config file path
		config_path = get_config_file_path()
		plot_file = os.path.join(os.path.dirname(config_path), 'last_plot_settings.json')
		
		if not os.path.exists(plot_file):
			update_log(self, 'No saved plot settings found')
			return
		
		# Load from JSON file
		with open(plot_file, 'r') as f:
			plot_settings = json.load(f)
		
		# Apply custom info settings
		self.custom_info_gb.setChecked(plot_settings.get('custom_info_enabled', False))
		self.model_cbox.setCurrentText(plot_settings.get('model', ''))
		self.ship_tb.setText(plot_settings.get('ship_name', 'R/V Unsinkable II'))
		self.cruise_tb.setText(plot_settings.get('cruise_name', 'A 3-hour tour'))
		self.show_model_chk.setChecked(plot_settings.get('show_model', True))
		self.show_ship_chk.setChecked(plot_settings.get('show_ship', True))
		self.show_cruise_chk.setChecked(plot_settings.get('show_cruise', True))
		
		# Apply data reference settings
		ref_data = plot_settings.get('reference_data', 'Waterline')
		ref_idx = self.ref_cbox.findText(ref_data)
		if ref_idx >= 0:
			self.ref_cbox.setCurrentIndex(ref_idx)
		
		tide_units = plot_settings.get('tide_units', 'Meter')
		tide_idx = self.tide_unit_cbox.findText(tide_units)
		if tide_idx >= 0:
			self.tide_unit_cbox.setCurrentIndex(tide_idx)
		
		self.waterline_tb.setText(plot_settings.get('waterline_adjustment', '0.00'))
		
		# Apply point style settings
		pt_size = plot_settings.get('point_size', '1')
		pt_size_idx = self.pt_size_cbox.findText(pt_size)
		if pt_size_idx >= 0:
			self.pt_size_cbox.setCurrentIndex(pt_size_idx)
		
		pt_size_cov = plot_settings.get('point_size_coverage', '5')
		pt_size_cov_idx = self.pt_size_cov_cbox.findText(pt_size_cov)
		if pt_size_cov_idx >= 0:
			self.pt_size_cov_cbox.setCurrentIndex(pt_size_cov_idx)
		
		# Handle both old and new opacity settings for backward compatibility
		if 'point_opacity_accuracy' in plot_settings:
			pt_alpha_acc = plot_settings.get('point_opacity_accuracy', '100')
			self.pt_alpha_acc_tb.setText(pt_alpha_acc)
			
			pt_alpha_cov = plot_settings.get('point_opacity_coverage', '6')
			self.pt_alpha_cov_tb.setText(pt_alpha_cov)
		else:
			# Old format: single opacity value applies to both
			pt_alpha = plot_settings.get('point_opacity', '100')
			self.pt_alpha_acc_tb.setText(pt_alpha)
			self.pt_alpha_cov_tb.setText(pt_alpha)
		
		# Apply plot limits settings
		self.plot_lim_gb.setChecked(plot_settings.get('custom_plot_limits', False))
		self.max_beam_angle_tb.setText(plot_settings.get('max_beam_angle', ''))
		self.angle_spacing_tb.setText(plot_settings.get('angle_spacing', ''))
		self.max_bias_tb.setText(plot_settings.get('max_bias', ''))
		self.max_std_tb.setText(plot_settings.get('max_std', ''))
		self.axis_margin_tb.setText(plot_settings.get('axis_margin', ''))
		self.tide_range_tb.setText(plot_settings.get('tide_range', '12'))
		
		# Apply toggle options
		self.show_acc_proc_text_chk.setChecked(plot_settings.get('show_acc_proc_text', False))
		self.show_ref_proc_text_chk.setChecked(plot_settings.get('show_ref_proc_text', False))
		self.grid_lines_toggle_chk.setChecked(plot_settings.get('show_grid_lines', True))
		self.IHO_lines_toggle_chk.setChecked(plot_settings.get('show_IHO_lines', True))
		self.update_ref_plots_chk.setChecked(plot_settings.get('update_ref_plots', True))
		self.show_xline_cov_chk.setChecked(plot_settings.get('show_xline_coverage', True))
		self.show_u_plot_chk.setChecked(plot_settings.get('show_uncertainty_plot', True))
		self.show_shaded_relief_chk.setChecked(plot_settings.get('show_shaded_relief', True))
		self.show_special_order_chk.setChecked(plot_settings.get('show_special_order', False))
		self.show_order_1a_chk.setChecked(plot_settings.get('show_order_1a', False))
		self.show_order_1b_chk.setChecked(plot_settings.get('show_order_1b', False))
		self.show_order_2_chk.setChecked(plot_settings.get('show_order_2', False))
		self.show_order_3_chk.setChecked(plot_settings.get('show_order_3', False))
		
		# Apply flatten swath settings
		self.flatten_mean_gb.setChecked(plot_settings.get('flatten_mean_enabled', False))
		self.set_zero_mean_chk.setChecked(plot_settings.get('force_zero_mean', False))
		self.mean_bias_angle_lim_tb.setText(plot_settings.get('mean_bias_angle_limit', '40'))
		
		# Get timestamp
		timestamp = plot_settings.get('saved_timestamp', 'Unknown')
		update_log(self, f'Plot settings loaded (saved: {timestamp})')
		
		# Force UI to process events so widget values are updated
		from PyQt6 import QtWidgets
		QtWidgets.QApplication.processEvents()
		
		# Refresh plots to apply the loaded settings
		refresh_plot(self, refresh_list=['acc', 'ref'], sender='load_plot_settings')
		
	except Exception as e:
		update_log(self, f'Error loading plot settings: {str(e)}')
		print(f"Error loading plot settings: {str(e)}")


def load_default_plot_settings(self):
	"""Reset all plot settings to their default values"""
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
		
		# Flatten swath defaults
		self.flatten_mean_gb.setChecked(False)
		self.set_zero_mean_chk.setChecked(False)
		self.mean_bias_angle_lim_tb.setText('40')
		
		update_log(self, 'All plot settings reset to default values')
		
		# Force UI to process events so widget values are updated
		from PyQt6 import QtWidgets
		QtWidgets.QApplication.processEvents()
		
		# Refresh plots to apply the default settings
		refresh_plot(self, refresh_list=['acc', 'ref'], sender='load_default_plot_settings')
		
	except Exception as e:
		update_log(self, f'Error loading default plot settings: {str(e)}')
		print(f"Error loading default plot settings: {str(e)}")


def get_config_file_path():
	"""Get the path to the configuration file"""
	return os.path.join(os.path.expanduser("~"), ".swath_accuracy_config.json")


def export_all_to_geotiff(self):
	"""Export all reference surface data as GeoTIFF files"""
	try:
		# Check if reference surface data is available
		if not hasattr(self, 'ref') or not self.ref:
			update_log(self, 'ERROR: No reference surface data available for export', font_color="red")
			return
		
		# Load last used directory
		config = load_session_config()
		last_directory = config.get("last_export_directory", os.getcwd())
		
		# First, prompt for parent directory
		parent_dir = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select parent directory for GeoTIFF export...', last_directory)
		
		if not parent_dir:  # User cancelled
			return
		
		# Then, prompt for export name with a custom dialog
		from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QDialogButtonBox
		
		# Create custom dialog for export name input
		dialog = QDialog(self)
		dialog.setWindowTitle('Enter export name')
		dialog.setModal(True)
		dialog.setMinimumWidth(500)
		
		# Create layout
		layout = QVBoxLayout()
		
		# Add label
		label = QLabel('Export name (will create directory with this name):')
		layout.addWidget(label)
		
		# Add input field
		input_field = QLineEdit()
		input_field.setText('reference_surface_export')
		input_field.setMinimumWidth(400)
		layout.addWidget(input_field)
		
		# Add button box
		button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
		button_box.accepted.connect(dialog.accept)
		button_box.rejected.connect(dialog.reject)
		layout.addWidget(button_box)
		
		dialog.setLayout(layout)
		
		# Show dialog and get result
		result = dialog.exec()
		ok = result == QDialog.DialogCode.Accepted
		export_name = input_field.text() if ok else ""
		
		if not ok or not export_name.strip():  # User cancelled or entered empty string
			return
		
		export_name = export_name.strip()
		
		# Create output directory
		output_dir = os.path.join(parent_dir, export_name)
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)
		
		# Save last directory
		save_last_directory(parent_dir)
		update_last_directory("last_export_directory", parent_dir)
		
		update_log(self, f'Exporting reference surface data to: {output_dir}')
		
		# Import required libraries for GeoTIFF export
		try:
			import rasterio
			from rasterio.transform import from_bounds
			from rasterio.crs import CRS
		except ImportError:
			update_log(self, 'ERROR: rasterio library not available. Please install it with: pip install rasterio', font_color="red")
			return
		
		# Get the extent and CRS information
		if 'z_ref_extent' in self.ref:
			extent = self.ref['z_ref_extent']
			# extent format: [xmin, xmax, ymin, ymax]
			bounds = [extent[0], extent[2], extent[1], extent[3]]  # [xmin, ymin, xmax, ymax]
		else:
			update_log(self, 'ERROR: Reference surface extent not available', font_color="red")
			return
		
		# Get UTM zone for CRS
		utm_zone = self.ref.get('utm_zone', '1N')  # Default to 1N if not available
		if 'N' in utm_zone:
			epsg_code = 32600 + int(utm_zone.replace('N', ''))
		elif 'S' in utm_zone:
			epsg_code = 32700 + int(utm_zone.replace('S', ''))
		else:
			epsg_code = 32601  # Default to UTM 1N
		
		crs = CRS.from_epsg(epsg_code)
		transform = from_bounds(*bounds, self.ref['z_grid'].shape[1], self.ref['z_grid'].shape[0])
		
		# Define the surfaces to export
		surfaces_to_export = []
		
		# Final Surface (masked depth)
		if 'z_final_grid' in self.ref and self.ref['z_final_grid'] is not None:
			surfaces_to_export.append(('final_surface', self.ref['z_final_grid'], 'Final Surface (m)'))
		
		# Depth
		if 'z_grid' in self.ref and self.ref['z_grid'] is not None:
			surfaces_to_export.append(('depth', self.ref['z_grid'], 'Depth (m)'))
		
		# Uncertainty
		if 'u_grid' in self.ref and self.ref['u_grid'] is not None:
			surfaces_to_export.append(('uncertainty', self.ref['u_grid'], 'Uncertainty (m)'))
		
		# Density
		if 'c_grid' in self.ref and self.ref['c_grid'] is not None:
			surfaces_to_export.append(('density', self.ref['c_grid'], 'Density (count)'))
		
		# Slope
		if 's_grid' in self.ref and self.ref['s_grid'] is not None:
			surfaces_to_export.append(('slope', self.ref['s_grid'], 'Slope (degrees)'))
		
		if not surfaces_to_export:
			update_log(self, 'ERROR: No valid surfaces available for export', font_color="red")
			return
		
		# Export each surface
		exported_files = []
		for surface_name, data_grid, description in surfaces_to_export:
			try:
				# Prepare data for export (replace NaN with nodata value)
				export_data = data_grid.copy()
				nodata_value = -9999
				export_data[np.isnan(export_data)] = nodata_value
				
				# Create filename
				filename = f"{export_name}_{surface_name}.tif"
				filepath = os.path.join(output_dir, filename)
				
				# Write GeoTIFF
				with rasterio.open(
					filepath,
					'w',
					driver='GTiff',
					height=export_data.shape[0],
					width=export_data.shape[1],
					count=1,
					dtype=export_data.dtype,
					crs=crs,
					transform=transform,
					nodata=nodata_value
				) as dst:
					dst.write(export_data, 1)
					dst.update_tags(TIFFTAG_DATETIME=description)
				
				exported_files.append(filename)
				update_log(self, f'Exported: {filename}')
				
			except Exception as e:
				update_log(self, f'ERROR: Failed to export {surface_name}: {str(e)}', font_color="red")
				print(f"Error exporting {surface_name}: {str(e)}")
		
		if exported_files:
			update_log(self, f'Successfully exported {len(exported_files)} GeoTIFF files to: {output_dir}')
			
			# Show success message
			from PyQt6.QtWidgets import QMessageBox
			file_list = "\n".join([f"  • {f}" for f in exported_files])
			msg = f"Successfully exported {len(exported_files)} GeoTIFF files:\n{file_list}\n\nLocation: {output_dir}"
			QMessageBox.information(self, 'Export Complete', msg)
		else:
			update_log(self, 'ERROR: No files were successfully exported', font_color="red")
		
	except Exception as e:
		update_log(self, f'ERROR: Export failed: {str(e)}', font_color="red")
		print(f"Error in export_all_to_geotiff: {str(e)}")