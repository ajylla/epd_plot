from epd_loader import *
from matplotlib.ticker import LinearLocator, MultipleLocator, AutoMinorLocator
import matplotlib.dates as mdates

def extract_data(df_electrons, plotstart, plotend, searchstart, searchend, bgstart, bgend, channels = 
                 [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33], averaging_mode='none', averaging=2):
    
    # Energy bin boundaries.
    e_low = [0.0312, 0.0330, 0.0348, 0.0380, 0.0406, 0.0432, 0.0459, 0.0497, 0.0533, 0.0580, 0.0627, 0.0673, 0.0731, 0.0788, 0.0856, 0.0934, 0.1011, 0.1109, 0.1197, 0.1305, 0.1423, 0.1541, 0.1679, 0.1835, 0.1995, 0.2181, 0.2371, 0.2578, 0.2817, 0.3061, 0.3339, 0.3661, 0.3989, 0.4348]
    e_high = [0.0348, 0.0369, 0.0380, 0.0406, 0.0432, 0.0459, 0.0497, 0.0533, 0.0580, 0.0627, 0.0673, 0.0731, 0.0788, 0.0856, 0.0934, 0.1011, 0.1109, 0.1197, 0.1305, 0.1423, 0.1541, 0.1679, 0.1835, 0.1995, 0.2181, 0.2371, 0.2578, 0.2817, 0.3061, 0.3339, 0.3661, 0.3989, 0.4348, 0.4714]
    
    # Takes electron flux and uncertainty values from original data.
    df_electron_fluxes = df_electrons['Electron_Flux']
    df_electron_uncertainties = df_electrons['Electron_Uncertainty']
    
    # Does averaging (if selected) on flux and uncertainty data and selects only the plot area range of data.
    if(averaging_mode=='none'):
        df_electron_fluxes = df_electron_fluxes[plotstart:plotend]
        df_electron_uncertainties = df_electron_uncertainties[plotstart:plotend]
    elif(averaging_mode=='rolling_window'):
        df_electron_fluxes = df_electron_fluxes.rolling(window=averaging, min_periods=1).mean()[plotstart:plotend]
        df_electron_uncertainties = df_electron_uncertainties.rolling(window=averaging, min_periods=1).mean()[plotstart:plotend]
    elif(averaging_mode=='rolling_mean'):
        df_electron_fluxes = df_electron_fluxes.resample('{}min'.format(averaging)).mean()[plotstart:plotend]
        df_electron_uncertainties = df_electron_uncertainties.resample('{}min'.format(averaging)).mean()[plotstart:plotend]
    
    # Main information dataframe containing most of the required data.
    df_info = pd.DataFrame({'Plot_start':[], 'Plot_end':[], 'Bg_start':[], 'Bg_end':[], 'Averaging':[], 'Energy_channel':[], 'Primary_energy':[], 'Energy_error_low':[], 'Energy_error_high':[], 'Peak_timestamp':[], 'Flux_peak':[], 'Peak_electron_uncertainty':[], 'Background_flux':[],'Bg_electron_uncertainty':[], 'Bg_subtracted_peak':[]})
    
    # Adds basic metadata to main info df.
    df_info['Plot_start'] = [plotstart]+['']*(len(channels)-1)
    df_info['Plot_end'] = [plotend]+['']*(len(channels)-1)
    df_info['Bg_start'] = [bgstart]+['']*(len(channels)-1)
    df_info['Bg_end'] = [bgend]+['']*(len(channels)-1)
    if(averaging_mode == 'none'):
        df_info['Averaging'] = ['No averaging']+['']*(len(channels)-1)
    elif(averaging_mode == 'rolling_window'):
        df_info['Averaging'] = ['Rolling window', 'Window size = ' + str(averaging)] + ['']*(len(channels)-2)
    elif(averaging_mode == 'rolling_mean'):
        df_info['Averaging'] = ['Rolling mean', 'Resampled to ' + str(averaging) + 'min'] + ['']*(len(channels)-2)
    
    # Energy bin primary energies; geometric mean.
    primary_energies = []
    for i in range(0,len(e_low)):
        primary_energies.append(np.sqrt(e_low[i]*e_high[i]))
    df_info['Primary_energy'] = primary_energies
        
    # Next blocks of code calculate information from data and append them to main info df.
    list_bg_fluxes = []
    list_flux_peaks = []
    list_peak_timestamps = []
    list_bg_subtracted_peaks = []
    list_peak_electron_uncertainties = []
    list_average_bg_uncertainties = []
    
    for channel in channels:
        bg_flux = df_electron_fluxes['Electron_Flux_{}'.format(channel)][bgstart:bgend].mean()
        list_bg_fluxes.append(bg_flux)
        
        flux_peak = df_electron_fluxes['Electron_Flux_{}'.format(channel)][searchstart:searchend].max()
        list_flux_peaks.append(flux_peak)
        
        peak_timestamp = df_electron_fluxes['Electron_Flux_{}'.format(channel)][searchstart:searchend].idxmax()
        list_peak_timestamps.append(peak_timestamp)
        
        peak_electron_uncertainty = df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)][peak_timestamp]
        list_peak_electron_uncertainties.append(peak_electron_uncertainty)
        
        average_bg_uncertainty = np.sqrt((df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)]
                                          [bgstart:bgend]**2).sum(axis=0))/len(df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)][bgstart:bgend])
        list_average_bg_uncertainties.append(average_bg_uncertainty)
        
    
    for i in range(0,len(list_flux_peaks)):
        list_bg_subtracted_peaks.append(list_flux_peaks[i]-list_bg_fluxes[i])
        
    df_info['Energy_channel'] = channels
    df_info['Background_flux'] = list_bg_fluxes
    df_info['Flux_peak'] = list_flux_peaks
    df_info['Peak_timestamp'] = list_peak_timestamps
    df_info['Bg_subtracted_peak'] = list_bg_subtracted_peaks
    df_info['Peak_electron_uncertainty'] = list_peak_electron_uncertainties
    df_info['Bg_electron_uncertainty'] = list_average_bg_uncertainties
    
    # Calculates energy errors for spectrum plot.
    energy_error_low = []
    energy_error_high = []
    for i in range(0,len(primary_energies)):
        energy_error_low.append(primary_energies[i]-e_low[i])
        energy_error_high.append(e_high[i]-primary_energies[i])
    df_info['Energy_error_low'] = energy_error_low
    df_info['Energy_error_high'] = energy_error_high
    
    return df_electron_fluxes, df_info, [searchstart, searchend], [e_low, e_high]

def plot_channels(args, bg_subtraction=False, savefig=False, key=''):
    
    hours = mdates.HourLocator(interval = 1)
    df_electron_fluxes = args[0]
    df_info = args[1]
    search_area = args[2]
    energy_bin = args[3]
    
    # If background subtraction is enabled, subtracts bg_flux from all observations. If flux value is negative, changes it to 0.
    if(bg_subtraction == False):
        pass
    elif(bg_subtraction == True):
        df_electron_fluxes = df_electron_fluxes.sub(df_info['Background_flux'].values, axis=1)
        df_electron_fluxes[df_electron_fluxes<0] = 0
        
    # Plotting part.
    # Initialized the main figure.
    fig = plt.figure()
    plt.xticks([])
    plt.yticks([])
    plt.ylabel("Flux \n [1/s cm$^2$ sr MeV]", labelpad=40)
    plt.xlabel("Time", labelpad=45)
    if(bg_subtraction==False):
        plt.title(str(df_info['Plot_start'][0][:-5]) + ", bg subtraction off")
    elif(bg_subtraction==True):
        plt.title(str(df_info['Plot_start'][0][:-5]) + ", bg subtraction on")
        
    # Loop through selected energy channels and creates a subplot for each.
    n=1
    for channel in df_info['Energy_channel']:
        ax = fig.add_subplot(len(df_info['Energy_channel']),1,n)
        ax = df_electron_fluxes['Electron_Flux_{}'.format(channel)].plot(logy=True, figsize=(20,40), color='red')
        
        plt.text(0.025,0.55, str(energy_bin[0][channel]) + " - " + str(energy_bin[1][channel]) + " MeV", transform=ax.transAxes)
        
        # Search area vertical lines.
        ax.axvline(search_area[0], color='black')
        ax.axvline(search_area[1], color='black')
        
        # Peak vertical line.
        ax.axvline(df_info['Peak_timestamp'][n-1], color='green')
        
        # Background measurement area.
        ax.axvspan(df_info['Bg_start'][0], df_info['Bg_end'][0], color='gray', alpha=0.25)
        
        if(n != len(df_info['Energy_channel'])):
            ax.get_xaxis().set_visible(False)
           
        plt.xlabel("")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
        ax.xaxis.set_minor_locator(hours)
        ax.get_xaxis().set_visible(True)
        
        n=n+1
        
    # Saves figure, if enabled.
    if(savefig):
        plt.savefig('/home/smurf/srl/images/channels-' + str(df_info['Plot_start'][0][:-5]) + '-' + str(key) + '.jpg', bbox_inches='tight')
        
    plt.show()
        
def plot_spectrum(args, bg_subtraction=True, savefig=False, key=''):
    
    df_info = args[1]
    
    # Plots either the background subtracted or raw flux peaks depending on choice.
    if(bg_subtraction):
        ax = df_info.plot.scatter(x='Primary_energy', y='Bg_subtracted_peak', c='red', label='Flux peaks', figsize=(13,10))
        ax.errorbar(x=df_info['Primary_energy'], y=df_info['Bg_subtracted_peak'], yerr=df_info['Peak_electron_uncertainty'], 
                    xerr=[df_info['Energy_error_low'], df_info['Energy_error_high']], fmt='.', ecolor='red', alpha=0.5)
        plt.title(str(df_info['Plot_start'][0][:-5]) + ' flux peaks, bg subtraction on')
    elif(bg_subtraction == False):
        ax = df_info.plot.scatter(x='Primary_energy', y='Flux_peak', c='red', label='Flux peaks', figsize=(13,10))
        ax.errorbar(x=df_info['Primary_energy'], y=df_info['Flux_peak'], yerr=df_info['Peak_electron_uncertainty'], 
                    xerr=[df_info['Energy_error_low'], df_info['Energy_error_high']], fmt='.', ecolor='red', alpha=0.5)
        plt.title(str(df_info['Plot_start'][0][:-5]) + ' flux peaks, bg subtraction off')
    
    # Plots background flux and background errorbars in same scatterplot.
    df_info.plot(kind='scatter', x='Primary_energy', y='Background_flux', c='red', alpha=0.25, ax=ax, label='Background flux')
    ax.errorbar(x=df_info['Primary_energy'], y=df_info['Background_flux'], yerr=df_info['Bg_electron_uncertainty'], xerr=[df_info['Energy_error_low'],df_info['Energy_error_high']],
                fmt='.', ecolor='red', alpha=0.15)
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Energy [MeV]')
    ax.set_ylabel('Flux \n [1/s cm$^2$ sr MeV]')
    plt.grid()
    
    if(savefig):
        plt.savefig('/home/smurf/srl/images/spectrum-' + str(df_info['Plot_start'][0][:-5]) + '-' + str(key) + '.jpg', dpi=300, bbox_inches='tight')
        
    plt.show()
    
def write_to_csv(args, key=''):
    
    df_info = args[1]
    
    df_info.to_csv('/home/smurf/srl/csv/electron_data-' + str(df_info['Plot_start'][0][:-5]) + '-' + str(key) + '.csv', index=False)