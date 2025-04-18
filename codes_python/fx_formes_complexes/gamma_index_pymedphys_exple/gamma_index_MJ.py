import pymedphys
import numpy as np
import matplotlib.pyplot as plt


axis_reference = data_exp[0] #1st column is x in mm
dose_reference = data_exp[1] #2nd column is dose in Gy

axis_evaluation = data_simu[0]
dose_evaluation = data_simu[1]
gamma_options = {
    'dose_percent_threshold': 1,
    'distance_mm_threshold': 1,
    'lower_percent_dose_cutoff': 10,
    'interp_fraction': 10,  # Should be 10 or more for more accurate results
    'max_gamma': 2,
    'random_subset': None,
    'local_gamma': True, # False indicates global gamma is calculated
    'ram_available': 5*10**9 # 1/2 GB
}

#for global dose normalization, the maximum reference dose is used
#but in TG218, it said usually the prescribed dose or the maximum dose in a plan (evaluation) is used
gamma = pymedphys.gamma(axis_reference, dose_reference,axis_evaluation, dose_evaluation,**gamma_options)

valid_gamma = gamma[~np.isnan(gamma)]
print('# of reference points with a valid gamma {0}'.format( len(valid_gamma)) )
num_bins = (
    gamma_options['interp_fraction'] * gamma_options['max_gamma'])
bins = np.linspace(0, gamma_options['max_gamma'], num_bins + 1)

if gamma_options['local_gamma']:
    gamma_norm_condition = 'Local gamma'
else:
    gamma_norm_condition = 'Global gamma'

pass_ratio = np.sum(valid_gamma <= 1) / len(valid_gamma)
print(pass_ratio)

width = 10
height = (width /(1.618))
fig.set_size_inches(width, height)


fig4,ax4 = plt.subplots()

ax4.hist(valid_gamma, bins, density=True) #y value is probability density in each bin
#plt.hist(valid_gamma, bins, density=False) #y value is counts in each bin
ax4.set(xlabel = 'gamma index of reference point', ylabel ='probability density', xlim = [0, gamma_options['max_gamma']]  )


max_ref_dose = np.max(dose_reference) #max reference dose
max_eva_dose = np.max(dose_evaluation) #max evaluation dose
lower_dose_cutoff = gamma_options['lower_percent_dose_cutoff']/100*max_ref_dose

fig11,ax11 = plt.subplots()
ax11.set(xlabel = r'$\textbf{Depth [mm]}$', ylabel = r'$\textbf{Dose [u.a]}$', ylim = [0, max(max_ref_dose,max_eva_dose) * 1.1] )
ax11.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


ax11_2 = ax11.twinx()
ax11_2.set(ylabel =r'$\textbf{Gamma index}$',ylim = [-0.25, gamma_options['max_gamma'] * 2.0 ], xlim = [-10,499])
curve_3 = ax11_2.plot(np.linspace(-500,500,20),np.zeros(len(np.linspace(-500,500,20))) + 1, ls = '--', lw = 2, label = r'$\textbf{Passing rate: ' +str(int(np.round(pass_ratio,2)*100)) + '\%}$' )

curve_0 = ax11.plot(axis_reference, dose_reference,'k-',label=r'$\textbf{Reference}$')
curve_1 = ax11.plot(axis_evaluation[1:], dose_evaluation[1:],'bo', mfc='none', markersize=5, label=r'$\textbf{Simulation}$')
str_gamma_p_dose = str(gamma_options['dose_percent_threshold'])
str_gamma_mm_dist = str(gamma_options['distance_mm_threshold'])
curve_2 = ax11_2.plot(axis_reference[1:], gamma[1:], 'r*', label=r'$\textbf{gamma ('+str_gamma_p_dose+'\%/'+str_gamma_p_dose +' mm)}$')
curves = curve_0 + curve_1 + curve_2 + curve_3

labels = [l.get_label() for l in curves]
fig11.subplots_adjust(top=0.9,
                         bottom=0.129,
                         left=0.09,
                         right=0.902,
                         hspace=0.2,
                         wspace=0.2)
ax11.legend(curves, labels)
ax11.grid(True)

fig11.set_size_inches(width, height)
plt.show()
