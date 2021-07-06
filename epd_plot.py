from epd_loader import *
from matplotlib.ticker import LinearLocator, MultipleLocator, AutoMinorLocator
import matplotlib.dates as mdates
from matplotlib.ticker import FormatStrFormatter
from adjustText import adjust_text

def extract_data(df_electrons, plotstart, plotend, searchstart, searchend, bgstart, bgend, instrument = 'ept', data_type = 'l2', averaging_mode='none', averaging=2):

    if(instrument != 'step'):
        df_electron_fluxes = df_electrons['Electron_Flux']
        df_electron_uncertainties = df_electrons['Electron_Uncertainty']

    # Takes electron flux and uncertainty values from original data.
    if(instrument != 'step'):
        df_electron_fluxes = df_electrons['Electron_Flux'][plotstart:plotend]
        df_electron_uncertainties = df_electrons['Electron_Uncertainty'][plotstart:plotend]

    # Pretty bad workarounds for implementing low-latency data.
    if(instrument == 'ept'):

        if(data_type == 'll'):

            channels = [0,1,2,3,4,5,6,7]

            for i in channels:
                df_electron_fluxes = df_electron_fluxes.rename(columns={'Ele_Flux_{}'.format(i):'Electron_Flux_{}'.format(i)})
                df_electron_uncertainties = df_electron_uncertainties.rename(columns={'Ele_Flux_Sigma_{}'.format(i):'Electron_Uncertainty_{}'.format(i)})

            e_low = [0.0329, 0.0411, 0.0537, 0.0733, 0.1013, 0.1425, 0.1997, 0.2821]
            e_high = [0.0411, 0.0537, 0.0733, 0.1013, 0.1425, 0.1997, 0.2821, 0.3977]

        elif(data_type == 'l2'):

            channels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33]

            e_low = [0.0312, 0.0330, 0.0348, 0.0380, 0.0406, 0.0432, 0.0459, 0.0497, 0.0533, 0.0580, 0.0627, 0.0673, 0.0731, 0.0788, 0.0856, 0.0934, 0.1011, 0.1109, 0.1197, 0.1305, 0.1423, 0.1541, 0.1679, 0.1835, 0.1995, 0.2181, 0.2371, 0.2578, 0.2817, 0.3061, 0.3339, 0.3661, 0.3989, 0.4348]
            e_high = [0.0348, 0.0369, 0.0380, 0.0406, 0.0432, 0.0459, 0.0497, 0.0533, 0.0580, 0.0627, 0.0673, 0.0731, 0.0788, 0.0856, 0.0934, 0.1011, 0.1109, 0.1197, 0.1305, 0.1423, 0.1541, 0.1679, 0.1835, 0.1995, 0.2181, 0.2371, 0.2578, 0.2817, 0.3061, 0.3339, 0.3661, 0.3989, 0.4348, 0.4714]

    elif(instrument == 'het'):

        if(data_type == 'll'):

            channels = [0,1,2,3]

            for i in channels:
                df_electron_fluxes = df_electron_fluxes.rename(columns={'Ele_Flux_{}'.format(i):'Electron_Flux_{}'.format(i)})
                df_electron_uncertainties = df_electron_uncertainties.rename(columns={'Ele_Flux_Sigma_{}'.format(i):'Electron_Uncertainty_{}'.format(i)})

            e_low = [0.4533, 1.0530, 2.4010, 5.9930]
            e_high = [1.0380, 2.4010, 5.9930, 18.8300]

        elif(data_type == 'l2'):

            channels = [0,1,2,3]

            e_low = [0.4533, 1.0530, 2.4010, 5.9930]
            e_high = [1.0380, 2.4010, 5.9930, 18.8300]

    elif(instrument == 'step'):

        if(data_type == 'l2'):

            channels = range(0,48)

            step_data = make_step_electron_flux(df_electrons)
            df_electron_fluxes = step_data[0][plotstart:plotend]
            df_electron_uncertainties = step_data[1][plotstart:plotend]

            e_low = step_data[2]
            e_high = step_data[3]

    if(averaging_mode == 'mean'):
        df_electron_fluxes = df_electron_fluxes.resample('{}min'.format(averaging)).mean()
        df_electron_uncertainties = df_electron_uncertainties.resample('{}min'.format(averaging)).apply(average_flux_error)
    elif(averaging_mode == 'rolling_window'):
        df_electron_fluxes = df_electron_fluxes.rolling(window=averaging, min_periods=1).mean()


    # Main information dataframe containing most of the required data.
    df_info = pd.DataFrame({'Plot_start':[], 'Plot_end':[], 'Bg_start':[], 'Bg_end':[], 'Averaging':[], 'Energy_channel':[], 'Primary_energy':[], 'Energy_error_low':[], 'Energy_error_high':[], 'Peak_timestamp':[], 'Flux_peak':[], 'Peak_significance':[], 'Peak_electron_uncertainty':[], 'Background_flux':[],'Bg_electron_uncertainty':[], 'Bg_subtracted_peak':[]})

    # Adds basic metadata to main info df.
    df_info['Plot_start'] = [plotstart]+['']*(len(channels)-1)
    df_info['Plot_end'] = [plotend]+['']*(len(channels)-1)
    df_info['Bg_start'] = [bgstart]+['']*(len(channels)-1)
    df_info['Bg_end'] = [bgend]+['']*(len(channels)-1)
    if(averaging_mode == 'none'):
        df_info['Averaging'] = ['No averaging']+['']*(len(channels)-1)
    elif(averaging_mode == 'rolling_window'):
        df_info['Averaging'] = ['Rolling window', 'Window size = ' + str(averaging)] + ['']*(len(channels)-2)
    elif(averaging_mode == 'mean'):
        df_info['Averaging'] = ['Mean', 'Resampled to ' + str(averaging) + 'min'] + ['']*(len(channels)-2)

    # Energy bin primary energies; geometric mean.
    primary_energies = []
    for i in range(0,len(e_low)):
        primary_energies.append(np.sqrt(e_low[i]*e_high[i]))

    primary_energies_channels = []
    for energy in channels:
        primary_energies_channels.append(primary_energies[energy])
    df_info['Primary_energy'] = primary_energies_channels

    # Next blocks of code calculate information from data and append them to main info df.
    list_bg_fluxes = []
    list_flux_peaks = []
    list_peak_timestamps = []
    list_bg_subtracted_peaks = []
    list_peak_electron_uncertainties = []
    list_average_bg_uncertainties = []
    list_bg_std = []
    list_peak_significance = []

    for channel in channels:
        bg_flux = df_electron_fluxes['Electron_Flux_{}'.format(channel)][bgstart:bgend].mean()
        list_bg_fluxes.append(bg_flux)

        flux_peak = df_electron_fluxes['Electron_Flux_{}'.format(channel)][searchstart:searchend].max()
        list_flux_peaks.append(flux_peak)

        peak_timestamp = df_electron_fluxes['Electron_Flux_{}'.format(channel)][searchstart:searchend].idxmax()
        list_peak_timestamps.append(peak_timestamp)

        # First finding the index location of the peak timestamp in uncertainty dataframe and the getting value of that index location.
        timestamp_loc = df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)].index.get_loc(peak_timestamp, method='nearest')
        peak_electron_uncertainty = df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)].iloc[timestamp_loc]

        #peak_electron_uncertainty = df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)][peak_timestamp]
        list_peak_electron_uncertainties.append(peak_electron_uncertainty)

        average_bg_uncertainty = np.sqrt((df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)]
                                          [bgstart:bgend]**2).sum(axis=0))/len(df_electron_uncertainties['Electron_Uncertainty_{}'.format(channel)][bgstart:bgend])
        list_average_bg_uncertainties.append(average_bg_uncertainty)

        bg_std = df_electron_fluxes['Electron_Flux_{}'.format(channel)][bgstart:bgend].std()
        list_bg_std.append(bg_std)


    for i in range(0,len(list_flux_peaks)):
        list_bg_subtracted_peaks.append(list_flux_peaks[i]-list_bg_fluxes[i])
        list_peak_significance.append(list_bg_subtracted_peaks[i]/list_bg_std[i])

    df_info['Energy_channel'] = channels
    df_info['Background_flux'] = list_bg_fluxes
    df_info['Flux_peak'] = list_flux_peaks
    df_info['Peak_timestamp'] = list_peak_timestamps
    df_info['Bg_subtracted_peak'] = list_bg_subtracted_peaks
    df_info['Peak_electron_uncertainty'] = list_peak_electron_uncertainties
    df_info['Bg_electron_uncertainty'] = list_average_bg_uncertainties
    df_info['Peak_significance'] = list_peak_significance

    # Calculates energy errors for spectrum plot.
    energy_error_low = []
    energy_error_high = []
    for i in range(0,len(primary_energies)):
        energy_error_low.append(primary_energies[i]-e_low[i])
        energy_error_high.append(e_high[i]-primary_energies[i])

    energy_error_low_channels = []
    energy_error_high_channels = []
    for i in channels:
        energy_error_low_channels.append(energy_error_low[i])
        energy_error_high_channels.append(energy_error_high[i])

    df_info['Energy_error_low'] = energy_error_low_channels
    df_info['Energy_error_high'] = energy_error_high_channels

    return df_electron_fluxes, df_info, [searchstart, searchend], [e_low, e_high], [instrument, data_type]

def make_step_electron_flux(stepdata, mask_conta=False):
    '''
    here we use the calibration factors from Paco (Alcala) to calculate the electron flux out of the (integral - magnet) fluxes (we now use level2 data files to get these)
    we also check if the integral counts are sufficiently higher than the magnet counts so that we can really assume it's electrons (otherwise we mask the output arrays)
    As suggested by Alex Kollhoff & Berger use a 5 sigma threshold:
    C_INT >> C_MAG:
    C_INT - C_MAG > 5*sqrt(C_INT)
    Alex: die count rates und fuer die uebrigen Zeiten gebe ich ein oberes Limit des Elektronenflusses an, das sich nach 5*sqrt(C_INT) /(E_f - E_i) /G_e berechnet.
    '''
    # calculate electron flux from F_INT - F_MAG:
    colnames = ["ch_num", "E_low", "E_hi", "factors"]
    paco = pd.read_csv('step_electrons_calibration.csv', names=colnames, skiprows=1)
    paco.E_low = round(paco.E_low/1000, 5)
    paco.E_hi = round(paco.E_hi/1000, 5)

    F_INT = stepdata['Integral_Flux']
    F_MAG = stepdata['Magnet_Flux']
    step_flux =  (F_INT - F_MAG) * paco.factors.values
    U_INT = stepdata['Integral_Uncertainty']
    U_MAG = stepdata['Magnet_Uncertainty']
    # from Paco:
    # Ele_Uncertainty = k * sqrt(Integral_Uncertainty^2 + Magnet_Uncertainty^2)
    step_unc = np.sqrt(U_INT*2 + U_MAG*2) * paco.factors.values
    param_list = ['Electron_Flux', 'Electron_Uncertainty']

    if mask_conta:
        C_INT = stepdata['Integral_Rate']
        C_MAG = stepdata['Magnet_Rate']
        clean = (C_INT - C_MAG) > 5*np.sqrt(C_INT)
        step_flux = step_flux.mask(clean)
        step_unc = step_unc.mask(clean)
    step_data = pd.concat([step_flux, step_unc], axis=1, keys=param_list)

    df_electron_fluxes = step_data['Electron_Flux']
    df_electron_uncertainties = step_data['Electron_Uncertainty']

    for channel in df_electron_fluxes:
        df_electron_fluxes = df_electron_fluxes.rename(columns={channel:'Electron_Flux_{}'.format(channel)})
    for channel in df_electron_uncertainties:
        df_electron_uncertainties = df_electron_uncertainties.rename(columns={channel:'Electron_Uncertainty_{}'.format(channel)})

    return df_electron_fluxes, df_electron_uncertainties, paco.E_low, paco.E_hi

def average_flux_error(flux_err: pd.DataFrame) -> pd.Series:
    return np.sqrt((flux_err ** 2).sum(axis=0)) / len(flux_err.values)

def plot_channels(args, bg_subtraction=False, savefig=False, key=''):

    hours = mdates.HourLocator(interval = 1)
    df_electron_fluxes = args[0]
    df_info = args[1]
    search_area = args[2]
    energy_bin = args[3]
    instrument = args[4][0]
    data_type = args[4][1]

    # If background subtraction is enabled, subtracts bg_flux from all observations. If flux value is negative, changes it to 0.
    if(bg_subtraction == False):
        pass
    elif(bg_subtraction == True):
        df_electron_fluxes = df_electron_fluxes.sub(df_info['Background_flux'].values, axis=1)
        df_electron_fluxes[df_electron_fluxes<0] = np.NaN

    # Plotting part.
    # Initialized the main figure.
    fig = plt.figure()
    plt.xticks([])
    plt.yticks([])
    plt.ylabel("Flux \n [1/s cm$^2$ sr MeV]", labelpad=40)
    plt.xlabel("Time", labelpad=45)
    if(bg_subtraction==False):
        plt.title(str(instrument) + ', ' + str(data_type) + ', ' + str(df_info['Plot_start'][0][:-5]) + ", bg subtraction off, " + str(df_info['Averaging'][0]) + ', ' + str(df_info['Averaging'][1]))
    elif(bg_subtraction==True):
        plt.title(str(instrument) + ', ' + str(data_type) + ', ' + str(df_info['Plot_start'][0][:-5]) + ", bg subtraction on, " + str(df_info['Averaging'][0]) + ', ' + str(df_info['Averaging'][1]))

    # Loop through selected energy channels and creates a subplot for each.
    n=1
    for channel in df_info['Energy_channel']:
        ax = fig.add_subplot(len(df_info['Energy_channel']),1,n)
        ax = df_electron_fluxes['Electron_Flux_{}'.format(channel)].plot(logy=True, figsize=(20,25), color='red')

        plt.text(0.025,0.7, str(energy_bin[0][channel]) + " - " + str(energy_bin[1][channel]) + " MeV", transform=ax.transAxes, size=13)

        # Search area vertical lines.
        ax.axvline(search_area[0], color='black')
        ax.axvline(search_area[1], color='black')

        # Peak vertical line.
        ax.axvline(df_info['Peak_timestamp'][n-1], color='green')

        # Background measurement area.
        ax.axvspan(df_info['Bg_start'][0], df_info['Bg_end'][0], color='gray', alpha=0.25)

        ax.get_xaxis().set_visible(False)

        if(n == len(df_info['Energy_channel'])):
            ax.get_xaxis().set_visible(True)

        plt.xlabel("")
        #ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
        #ax.xaxis.set_minor_locator(hours)

        n+=1

    # Saves figure, if enabled.
    if(savefig):
        if(bg_subtraction):
            plt.savefig('/home/smurf/srl/images/channels-' + str(df_info['Plot_start'][0][:-5]) + '-' + instrument + '-' + data_type + '-' + str(df_info['Averaging'][1]) + '-bg_subtracted' + str(key) +'.jpg', bbox_inches='tight')
        else:
            plt.savefig('/home/smurf/srl/images/channels-' + str(df_info['Plot_start'][0][:-5]) + '-' + instrument + '-' + data_type + '-' + str(df_info['Averaging'][1]) + str(key) +'.jpg', bbox_inches='tight')

    plt.show()

def plot_spectrum(args, bg_subtraction=True, savefig=False, key=''):

    df_info = args[1]
    instrument = args[4][0]
    data_type = args[4][1]

    # Plots either the background subtracted or raw flux peaks depending on choice.
    if(bg_subtraction):
        ax = df_info.plot.scatter(x='Primary_energy', y='Bg_subtracted_peak', c='red', label='Flux peaks', figsize=(13,10))
        ax.errorbar(x=df_info['Primary_energy'], y=df_info['Bg_subtracted_peak'], yerr=df_info['Peak_electron_uncertainty'],
                    xerr=[df_info['Energy_error_low'], df_info['Energy_error_high']], fmt='.', ecolor='red', alpha=0.5)
        plt.title(str(instrument) + ', ' + str(data_type) + ', ' + str(df_info['Plot_start'][0][:-5]) + ' flux peaks, bg subtraction on, ' + str(df_info['Averaging'][0]) + ', ' + str(df_info['Averaging'][1]), size=18)
    elif(bg_subtraction == False):
        ax = df_info.plot.scatter(x='Primary_energy', y='Flux_peak', c='red', label='Flux peaks', figsize=(13,10))
        ax.errorbar(x=df_info['Primary_energy'], y=df_info['Flux_peak'], yerr=df_info['Peak_electron_uncertainty'],
                    xerr=[df_info['Energy_error_low'], df_info['Energy_error_high']], fmt='.', ecolor='red', alpha=0.5)
        plt.title(str(instrument) + ', ' + str(data_type) + ', ' + str(df_info['Plot_start'][0][:-5]) + ' flux peaks, bg subtraction off, ' + str(df_info['Averaging'][0]) + ', ' + str(df_info['Averaging'][1]), size=18)

    # Plots background flux and background errorbars in same scatterplot.
    df_info.plot(kind='scatter', x='Primary_energy', y='Background_flux', c='red', alpha=0.25, ax=ax, label='Background flux')
    ax.errorbar(x=df_info['Primary_energy'], y=df_info['Background_flux'], yerr=df_info['Bg_electron_uncertainty'], xerr=[df_info['Energy_error_low'],df_info['Energy_error_high']],
                fmt='.', ecolor='red', alpha=0.15)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Energy [MeV]', size=20)
    ax.set_ylabel('Flux \n [1/s cm$^2$ sr MeV]', size=20)
    plt.tick_params(axis='x', which='minor', labelsize=16)
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))
    #plt.tick_params(axis='y', which='minor')
    #ax.yaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
    plt.legend(prop={'size': 18})
    plt.xticks(size=16)
    plt.yticks(size=16)
    plt.grid()

    for label in ax.xaxis.get_ticklabels(which='minor')[1::2]:
        label.set_visible(False)

    if(savefig):
        if(bg_subtraction):
            plt.savefig('/home/smurf/srl/images/spectrum-' + str(df_info['Plot_start'][0][:-5]) + '-' + instrument + '-' + data_type + '-' + str(df_info['Averaging'][1]) + '-bg_subtracted' + str(key) +'.jpg', dpi=300, bbox_inches='tight')
        else:
            plt.savefig('/home/smurf/srl/images/spectrum-' + str(df_info['Plot_start'][0][:-5]) + '-' + instrument + '-' + data_type + '-' + str(df_info['Averaging'][1]) + str(key) +'.jpg', dpi=300, bbox_inches='tight')

    plt.show()

def write_to_csv(args, key=''):

    df_info = args[1]
    instrument = args[4][0]
    data_type = args[4][1]

    df_info.to_csv('/home/smurf/srl/csv/electron_data-' + str(df_info['Plot_start'][0][:-5]) + '-' + instrument + '-' + data_type + '-' + str(df_info['Averaging'][1]) + str(key) + '.csv', index=False)

def acc_flux(args, time=[]):

    df_electron_fluxes = args[0]
    df_info = args[1]

    # If no timeframe specified, use search area.
    if(time==[]):
        time = args[2]

    # Calculates average fluxes for each enery channel from given timeframe and appends to list.
    list_flux_averages = []
    for channel in df_info['Energy_channel']:
        list_flux_averages.append(df_electron_fluxes['Electron_Flux_{}'.format(channel)][time[0]:time[1]].mean())

    df_acc = pd.DataFrame({'Primary_energy':[], 'Acc_flux':[]})
    df_acc['Primary_energy'] = df_info['Primary_energy']
    df_acc['Acc_flux'] = list_flux_averages

    ax = df_acc.plot(kind='scatter', x='Primary_energy', y='Acc_flux', logy=True, logx=True, color='green', figsize=(13,10))
