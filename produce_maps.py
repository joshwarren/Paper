## Routine to produce maps in a single grid for the MNRAS paper.

from checkcomp import checkcomp
cc = checkcomp()
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import matplotlib.pyplot as plt # used for plotting
import cPickle as pickle
import numpy as np 
from plot_velfield_nointerp import plot_velfield_nointerp
from prefig import Prefig
from astropy.io import fits 
from sauron_colormap import sauron#2 as sauron


from Bin import myArray
class Ds(object):
	def __init__(self):
		self.x=np.array([0,0,0,1,1,1,2,2,40])
		self.y=np.array([0,1,2,0,1,2,0,1,40])
		self.bin_num = np.array([0,0,1,0,1,1,2,2,3])
		self.xBar = np.array([0.5,1.5,2,40])
		self.yBar = np.array([0.5,1.5,1,40])
		self.SNRatio = np.array([0,1,1,2])
		if 'home' in cc.device:
			self.unbinned_flux = np.zeros((40,40))
		else:
			self.unbinned_flux = np.zeros((150,150))
		self.number_of_bins = 4
		self.components = {'stellar':comp(), '[OIII]5007d':comp()}
		self.flux = np.array([0,1,1,2])

	def absorption_line(self, p, uncert=False, remove_badpix=False):
		if uncert:
			return np.array([0,1,1,2]), np.array([0,1,1,2])
		else:
			return np.array([0,1,1,2])

class comp(object):
	def __init__(self):
		self.plot = {'vel':myArray(np.array([0,1,1,2]), np.array([0,1,1,2])), 
			'sigma': myArray(np.array([0,1,1,2]), np.array([0,1,1,2]))}
		self.flux = np.array([0,1,1,2])


def plot(galaxy, instrument='vimos', debug=True):
	print galaxy, instrument

	# opt = 'kin'
	overplot={'CO':'c', 'radio':'g'}
	if instrument=='vimos':
		from plot_results import add_, set_lims
		from errors2 import get_dataCubeDirectory
		ext = 0
		rows = 5

		plots = [[
			'D.flux',
			"D2.components['[OIII]5007d'].flux",
			"D.components['stellar'].plot['vel']",
			"D.components['stellar'].plot['vel'].uncert",
			"D.components['stellar'].plot['sigma']",
			"D.components['stellar'].plot['sigma'].uncert"
			],[
			"D2.absorption_line('G4300')",
			"D2.absorption_line('G4300',uncert=True)[1]",
			"D2.absorption_line('Fe4383')",
			"D2.absorption_line('Fe4383',uncert=True)[1]",
			"D2.absorption_line('Ca4455')",
			"D2.absorption_line('Ca4455',uncert=True)[1]"
			],[
			"D2.absorption_line('Fe4531')",
			"D2.absorption_line('Fe4531',uncert=True)[1]",
			"D2.absorption_line('H_beta')",
			"D2.absorption_line('H_beta',uncert=True)[1]",
			"D2.absorption_line('Fe5015')",
			"D2.absorption_line('Fe5015',uncert=True)[1]"
			],[
			"D2.absorption_line('Mg_b')",
			"D2.absorption_line('Mg_b',uncert=True)[1]",
			'',
			'',
			'',
			''
			]]

		str_plots = [[
			'Flux',
			r'[OIII]$\lambda\lambda$4959,5007 Flux',
			"Stellar mean velocity",
			"Stellar mean velocity\nuncertainty",
			"Stellar velocity dispersion",
			"Stellar velocity dispersion\nuncertainty"
			],[
			"G4300",
			"G4300 uncertainty",
			"Fe4383",
			"Fe4383 uncertainty",
			"Ca4455",
			"Ca4455 uncertainty"
			],[
			"Fe4531",
			"Fe4531 uncertainty",
			r'H$\,\beta$',
			r'H$\,\beta$ uncertainty',
			"Fe5015",
			"Fe5015 uncertainty"
			],[
			r'Mg$\,$b',
			r'Mg$\,$b uncertainty',
			'',
			'',
			'',
			''
			]]
		units = [[
			'',
			'',
			r'km s$^{-1}$',
			r'km s$^{-1}$',
			r'km s$^{-1}$',
			r'km s$^{-1}$'
			],[
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$'
			],[
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$'
			],[
			r'$\AA$',
			r'$\AA$',
			'',
			'',
			'',
			''
			]]

	elif instrument == 'muse':
		from plot_results_muse import add_, set_lims
		from errors2_muse import get_dataCubeDirectory
		ext = 1
		rows = 6
		
		plots = [[
			'D.flux',
			"D2.components['[OIII]5007d'].flux",
			"D.components['stellar'].plot['vel']",
			"D.components['stellar'].plot['vel'].uncert",
			"D.components['stellar'].plot['sigma']",
			"D.components['stellar'].plot['sigma'].uncert"
			],[
			"D2.absorption_line('H_beta')",
			"D2.absorption_line('H_beta',uncert=True)[1]",
			"D2.absorption_line('Fe5015')",
			"D2.absorption_line('Fe5015',uncert=True)[1]",
			"D2.absorption_line('Mg_b')",
			"D2.absorption_line('Mg_b',uncert=True)[1]"
			],[
			"D2.absorption_line('Fe5270')",
			"D2.absorption_line('Fe5270',uncert=True)[1]",
			"D2.absorption_line('Fe5335')",
			"D2.absorption_line('Fe5335',uncert=True)[1]",
			"D2.absorption_line('Fe5406')",
			"D2.absorption_line('Fe5406',uncert=True)[1]"
			],[
			"D2.absorption_line('Fe5709')",
			"D2.absorption_line('Fe5709',uncert=True)[1]",
			"D2.absorption_line('Fe5782')",
			"D2.absorption_line('Fe5782',uncert=True)[1]",
			"D2.absorption_line('NaD')",
			"D2.absorption_line('NaD',uncert=True)[1]"
			],[
			"D2.absorption_line('TiO1',remove_badpix=True)",
			"D2.absorption_line('TiO1',uncert=True,remove_badpix=True)[1]",
			"D2.absorption_line('TiO2',remove_badpix=True)",
			"D2.absorption_line('TiO2',uncert=True,remove_badpix=True)[1]",
			'',
			''
			]]

		str_plots = [[
			'Flux',
			r'[OIII]$\lambda\lambda$4959,5007 Flux',
			"Stellar mean velocity",
			"Stellar mean velocity\nuncertainty",
			"Stellar velocity dispersion",
			"Stellar velocity dispersion\nuncertainty"
			],[
			r'H$\,\beta$',
			r'H$\,\beta$ uncertainty',
			"Fe5015",
			"Fe5015 uncertainty",
			r'Mg$\,$b',
			r'Mg$\,$b uncertainty'
			],[
			'Fe5270',
			'Fe5270 uncertainty',
			'Fe5335',
			'Fe5335 uncertainty',
			'Fe5406',
			'Fe5406 uncertainty'
			],[
			'Fe5709',
			'Fe5709 uncertainty',
			'Fe5782',
			'Fe5782 uncertainty',
			'NaD',
			'NaD uncertainty'
			],[
			'TiO1',
			'TiO1 uncertainty',
			'TiO2',
			'TiO2 uncertainty',
			'',
			''
			]]
		units = [[
			'',
			'',
			r'km s$^{-1}$',
			r'km s$^{-1}$',
			r'km s$^{-1}$',
			r'km s$^{-1}$'
			],[
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$'
			],[
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$'
			],[
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$',
			r'$\AA$'
			],[
			'mag',
			'mag',
			'mag',
			'mag',
			'',
			''
			]]
	if galaxy in ['ngc0612', 'pks0718-34']:
		rows -= 1 # No Mgb
	if galaxy in ['ic1459', 'ngc0612', 'ngc1316', 'ngc3100']:
		rows += 1
		plots.insert(1, [
			"D2.components['[OIII]5007d'].plot['vel']",
			"D2.components['[OIII]5007d'].plot['vel'].uncert",
			"D2.components['[OIII]5007d'].plot['sigma']",
			"D2.components['[OIII]5007d'].plot['sigma'].uncert",
			'',
			''
			])
		str_plots.insert(1, [
			r"[OIII]$\lambda\lambda$4959,5007"+"\nmean velocity",
			r"[OIII]$\lambda\lambda$4959,5007"+"\nmean velocity uncertainty",
			r"[OIII]$\lambda\lambda$4959,5007"+"\nvelocity dispersion",
			r"[OIII]$\lambda\lambda$4959,5007"+"\nvelocity dispersion\nuncertainty",
			'',
			''
			])
		units.insert(1,[
			r'km s$^{-1}$',
			r'km s$^{-1}$',
			r'km s$^{-1}$',
			r'km s$^{-1}$',
			'',
			''
			])

	Prefig(size=np.array((6, rows))*6)
	fig, axs = plt.subplots(rows, 6)
	
	out_dir = '%s/Documents/paper/' % (cc.home_dir)

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[ext].header
	f.close()


	if debug: 
		D = Ds()
		D2 = Ds()

	vin_dir = '%s/Data/%s/analysis' % (cc.base_dir, instrument)
	data_file =  "%s/galaxies.txt" % (vin_dir)
	file_headings = np.loadtxt(data_file, dtype=str)[0]
	col = np.where(file_headings=='SN_kin')[0][0]
	col2 = np.where(file_headings=='SN_pop')[0][0]
	SN_target_kin_gals, SN_target_pop_gals = np.loadtxt(data_file, 
		unpack=True, skiprows=1, usecols=(col,col2))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	SN_target_kin=SN_target_kin_gals[i_gal]
	SN_target_pop=SN_target_pop_gals[i_gal]

	# attr, vmin, vmax = np.loadtxt('%s/lims.txt' % (vin_dir), dtype=str, 
	# 	usecols=(0,1,2), skiprows=1, unpack=True)
	# vmin, vmax = vmin.astype(float), vmax.astype(float)

	vin_dir2 = str(vin_dir + '/%s/pop' % (galaxy)) 
	vin_dir += '/%s/kin' % (galaxy) 

	if not debug:
		pickle_file = '%s/pickled' % (vin_dir)
		pickleFile = open("%s/dataObj.pkl" % (pickle_file), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()
		pickle_file2 = '%s/pickled' % (vin_dir2)
		pickleFile2 = open("%s/dataObj.pkl" % (pickle_file2), 'rb')
		D2 = pickle.load(pickleFile2)
		pickleFile2.close()		

	for i, row in enumerate(plots):
		for j, p in enumerate(row):
			print i, j, p
			if galaxy in ['ngc0612', 'pks0718-34'] and 'Mg_b' in p:
				break # break out of for-loop
			elif 'flux' in p:
				if p == 'D.flux':
					axs[i,j] = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
						D.xBar, D.yBar, eval(p), header, cmap='gist_yarg', 
						flux_unbinned=D.unbinned_flux, galaxy=str_plots[i][j],
						ax=axs[i,j])
				elif p == "D2.components['[OIII]5007d'].flux":
					axs[i,j] = plot_velfield_nointerp(D2.x, D2.y, D2.bin_num, 
						D2.xBar, D2.yBar, eval(p), header, cmap='gist_yarg', 
						flux_unbinned=D2.unbinned_flux, galaxy=str_plots[i][j],
						ax=axs[i,j])

				if overplot and p != '' and 'uncert' not in p:
					for o, color in overplot.iteritems():
						scale = 'log' if o == 'radio' else 'lin'
						add_(o, color, axs[i,j], galaxy, nolegend=True, scale=scale)
				if i > 0:
					if plots[i-1][j] != '':
						axs[i-1,j].ax_dis.set_xticklabels([])
						axs[i-1,j].ax_dis.set_xlabel('')
				if j > 0:
					if plots[i][j-1] != '':
						axs[i,j].ax_dis.set_yticklabels([])
						axs[i,j].ax_dis.set_ylabel('')

			elif p != '':
				vmin, vmax = set_lims(eval(p), 
					symmetric = 'vel' in p and 'uncert' not in p,
					positive = not ('vel' in p and 'uncert' not in p))
				if 'D2' in p:
					axs[i,j] = plot_velfield_nointerp(D2.x, D2.y, D2.bin_num, 
						D2.xBar, D2.yBar, eval(p), header, vmin=vmin, vmax=vmax, 
						flux_unbinned=D2.unbinned_flux, galaxy=str_plots[i][j],
						signal_noise=D2.SNRatio, signal_noise_target=SN_target_pop, 
						ax=axs[i,j], lim_labels=True, 
						lim_labels_units=units[i][j])
				else:
					axs[i,j] = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
						D.xBar, D.yBar, eval(p), header, vmin=vmin, vmax=vmax, 
						flux_unbinned=D.unbinned_flux, galaxy=str_plots[i][j],
						signal_noise=D.SNRatio, signal_noise_target=SN_target_kin, 
						ax=axs[i,j], lim_labels=True,
						lim_labels_units=units[i][j])

				if overplot and p != '' and 'uncert' not in p:
					for o, color in overplot.iteritems():
						scale = 'log' if o == 'radio' else 'lin'
						add_(o, color, axs[i,j], galaxy, nolegend=True, scale=scale)
				if i > 0:
					if plots[i-1][j] != '':
						axs[i-1,j].ax_dis.set_xticklabels([])
						axs[i-1,j].ax_dis.set_xlabel('')
				if j > 0:
					if plots[i][j-1] != '':
						axs[i,j].ax_dis.set_yticklabels([])
						axs[i,j].ax_dis.set_ylabel('')
			
				
			else:
				axs[i,j].remove()


	age = np.zeros(D2.number_of_bins)
	met = np.zeros(D2.number_of_bins)
	alp = np.zeros(D2.number_of_bins)
	unc_age = np.zeros(D2.number_of_bins)
	unc_met = np.zeros(D2.number_of_bins)
	unc_alp = np.zeros(D2.number_of_bins)

	if not debug:
		for j in xrange(D2.number_of_bins):
			ag, me, al = np.loadtxt('%s/pop/distribution/%i.dat' % (
				vin_dir2, j), unpack=True)

			for plot, unc_plot, pop in zip([age,met,alp],
				[unc_age,unc_met,unc_alp], [ag,me,al]):

				hist = np.histogram(pop, bins=40)
				x = (hist[1][0:-1]+hist[1][1:])/2
				hist = hist[0]
				plot[j] = x[np.argmax(hist)]

				gt_fwhm = hist >= np.max(hist)/2
				unc_plot[j] = np.max(x[gt_fwhm]) - np.min(x[gt_fwhm])

	str_plots = ['Age', 'Age uncertainty', 'Metalicity', 'Metalicity uncertainty', 
		'Alpha enhancement', 'Alpha enhancement\nuncertainty']
	units = ['Gyr', 'Gyr', None, None, None, None]
	vmin = [0, 0, -2.25, 0, -0.3, 0]
	vmax = [15, 2, 0.67, 0.4, 0.5, 0.25]
	for j, p in enumerate(['age', 'unc_age', 'met', 'unc_met', 'alp', 'unc_alp']):
		axs[-1, j] = plot_velfield_nointerp(D2.x, D2.y, D2.bin_num, 
			D2.xBar, D2.yBar, eval(p), header,  
			vmin=vmin[j], vmax=vmax[j], 
			# cmap='inferno', 
			flux_unbinned=D2.unbinned_flux, galaxy=str_plots[j],
			signal_noise=D2.SNRatio, signal_noise_target=SN_target_pop, 
			ax=axs[-1,j], lim_labels=True, lim_labels_units=units[j])
		if j > 0:
			axs[-1,j].ax_dis.set_yticklabels([])
			axs[-1,j].ax_dis.set_ylabel('')

		if plots[-1][j] != '' or galaxy in ['ngc0612', 'pks0718-34']:
			axs[-2, j].ax_dis.set_xticklabels([])
			axs[-2, j].ax_dis.set_xlabel('')
		if overplot and 'unc' not in p:
			for o, color in overplot.iteritems():
				scale = 'log' if o == 'radio' else 'lin'
				add_(o, color, axs[-1,j], galaxy, nolegend=True, scale=scale)

	for a in axs.flatten():
		if hasattr(a, 'ax_dis'):
			a.ax_dis.tick_params(top=True, bottom=True, left=True, 
				right=True, direction='in', which='major', length=20,
				width=3, labelsize='large', zorder=10)
			a.ax_dis.tick_params(top=True, bottom=True, left=True, 
				right=True, direction='in', which='minor', length=10,
				width=3, zorder=10)
			a.ax_dis.xaxis.label.set_size(22)
			a.ax_dis.yaxis.label.set_size(22)

	# Create gap between pairs of columns
	for j in range(1, 6, 2):
		for a in axs[:, j].flatten():
			ax_loc = a.get_position()
			ax_loc.x0 -= 0.007
			ax_loc.x1 -= 0.007

			a.set_position(ax_loc)
			if hasattr(a, 'ax_dis'):
				a.ax_dis.set_position(ax_loc)
	# for j in range(0, 6, 2):
	# 	for a in axs[:, j+1].flatten():
	# 		ax_loc = a.get_position()
	# 		ax_loc.x0 += 0.01
	# 		ax_loc.x1 += 0.01

	# 		a.set_position(ax_loc)
	# 		if hasattr(a, 'ax_dis'):
	# 			a.ax_dis.set_position(ax_loc)


	# fig.text(0.24, 0.9, r'Flux', va='top', ha='center', size='xx-large')
	# fig.text(0.51, 0.9, r'Velocty', va='top', ha='center', size='xx-large')
	# fig.text(0.8, 0.9, r'Velocty Dispersion', va='top', ha='center', 
	# 	size='xx-large')

	
	# Add colorbar
	ax_loc = axs[0,-1].get_position()
	cax = fig.add_axes([ax_loc.x1+0.03, ax_loc.y0, 0.02, ax_loc.height])
	cbar = plt.colorbar(axs[1,1].cs, cax=cax)
	cbar.ax.set_yticklabels([])

	# plt.show()
	fig.savefig('%s/%s.png' % (out_dir, galaxy), bbox_inches='tight',
		dpi=200)



if __name__=='__main__':
	# plot('ngc1316', instrument='muse', debug=True)
	if 'home' in cc.device:
		for galaxy in ['eso443-g024', 'ic1531', 'ngc0612', 'ngc3100', 'ngc3557', 
			'ngc7075', 'pks0718-34']:
			plot(galaxy, instrument='vimos', debug=False)
	elif cc.device == 'uni':
		for galaxy in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
			plot(galaxy, instrument='muse', debug=False)