from epd_loader import *
from matplotlib.ticker import LinearLocator, MultipleLocator, AutoMinorLocator
import matplotlib.dates as mdates

"""
Function to plot electron fluxes.

Function arguments:
df_electrons: Electron data acquired from epd_loader script.
plotstart: Start time of the plot.
plotend: End time of the plot.
searchstart: Start time of the peak search area.
searchend: End time of the peak search area.
bgstart: Start time of the background area.
bgend: End time of the background area.
!!! All start and end times are give in format "YYYY-MM-DD-HHMM" (string) !!!
flux_list: List of flux energy channels to be plotted. Defaults to all channels if missing. Parameter format is list. Example: [0,2,4] would plot channels 0, 2 and 4. 

Returns:
Function returns a list of plot data.
"""
def electron_flux_plot(df_electrons, plotstart, plotend, searchstart, searchend, bgstart, bgend, 
                       flux_list = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33], bg_subtraction=False):
    
    plot_area = [plotstart, plotend]
    search_area = [searchstart, searchend]
    bg_area = [bgstart, bgend]
    hours = mdates.HourLocator(interval = 1)
    
    # Gets the uncertainty values
    df_electron_uncertainty = df_electrons['Electron_Uncertainty']
    df_electron_uncertainty = df_electron_uncertainty[plot_area[0]:plot_area[1]]
    
    df_bg_electron_uncertainty = df_electron_uncertainty[bg_area[0]:bg_area[1]]
    #df_bg_electron_uncertainty = df_bg_electron_uncertainty.resample(averaging,label='left').apply(average_flux_error)
    
    # Gets just the electron flux data from specified time range.
    df_electrons_trim = df_electrons['Electron_Flux']
    df_electrons_fluxes = df_electrons_trim[plot_area[0]:plot_area[1]]
    
    electron_energies = [['0.0312 - 0.0348 MeV'], \
       ['0.0330 - 0.0369 MeV'], \
       ['0.0348 - 0.0380 MeV'], \
       ['0.0380 - 0.0406 MeV'], \
       ['0.0406 - 0.0432 MeV'], \
       ['0.0432 - 0.0459 MeV'], \
       ['0.0459 - 0.0497 MeV'], \
       ['0.0497 - 0.0533 MeV'], \
       ['0.0533 - 0.0580 MeV'], \
       ['0.0580 - 0.0627 MeV'], \
       ['0.0627 - 0.0673 MeV'], \
       ['0.0673 - 0.0731 MeV'], \
       ['0.0731 - 0.0788 MeV'], \
       ['0.0788 - 0.0856 MeV'], \
       ['0.0856 - 0.0934 MeV'], \
       ['0.0934 - 0.1011 MeV'], \
       ['0.1011 - 0.1109 MeV'], \
       ['0.1109 - 0.1197 MeV'], \
       ['0.1197 - 0.1305 MeV'], \
       ['0.1305 - 0.1423 MeV'], \
       ['0.1423 - 0.1541 MeV'], \
       ['0.1541 - 0.1679 MeV'], \
       ['0.1679 - 0.1835 MeV'], \
       ['0.1835 - 0.1995 MeV'], \
       ['0.1995 - 0.2181 MeV'], \
       ['0.2181 - 0.2371 MeV'], \
       ['0.2371 - 0.2578 MeV'], \
       ['0.2578 - 0.2817 MeV'], \
       ['0.2817 - 0.3061 MeV'], \
       ['0.3061 - 0.3339 MeV'], \
       ['0.3339 - 0.3661 MeV'], \
       ['0.3661 - 0.3989 MeV'], \
       ['0.3989 - 0.4348 MeV'], \
       ['0.4348 - 0.4714 MeV']] 
    
    # Average electron uncertainty for bg area.
    df_average_bg_electron_uncertainty = pd.DataFrame({'Energy channel':[],
                                                       'Bg electron uncertainty':[]})
    
    # Electron uncertainty for the peaks.
    df_peak_electron_uncertainty = pd.DataFrame({'Energy channel':[],
                                                 'Electron uncertainty':[]})
    
    # Observations from background area
    df_bg_obs = df_electrons_fluxes[bg_area[0]:bg_area[1]]
    
    # Background flux dataframe
    df_bg_fluxes = pd.DataFrame({'Energy channel':[],
                                 'Background flux':[]})
    
    #sqrt( (unc**2).sum(axis=0) ) / len(unc)
    
    
    for i in range(0,34):
        # Appends background data to background dataframe.
        df_bg_fluxes = df_bg_fluxes.append(pd.DataFrame({'Energy channel':[i],
                                                         'Background flux':[df_bg_obs.iloc[:,i].mean()]}), ignore_index=False)
        
        # Averages the bg electron uncertaintys and appends them to df. Probably not the final implementation.
        df_average_bg_electron_uncertainty = df_average_bg_electron_uncertainty.append(pd.DataFrame({'Energy channel':[i],
                                                                                                     'Bg electron uncertainty':[np.sqrt((df_bg_electron_uncertainty.iloc[:,i]**2).sum(axis=0))/len(df_bg_electron_uncertainty.iloc[:,i])]}), ignore_index=False)
    
    # Peak dataframe
    df_peaks = pd.DataFrame({'Energy channel':[],
                             'Flux peak':[],
                             'Timestamp':[]})
    
    df_peaks_subtracted = pd.DataFrame({'Energy channel':[],
                                        'Subtracted flux peak':[],
                                        'Timestamp':[]})
    
    # Background subtracted fluxes
    df_electrons_fluxes_subtracted = df_electrons_fluxes.sub(df_bg_fluxes.iloc[:,1].values,axis=1)
    
    # Main figure.
    fig = plt.figure(1)
    plt.ylabel("Flux \n [1/s cm$^2$ sr MeV]")
    plt.title(plotstart[:-5])
    
    # Indexing variable for the plot order.
    n = 1
    
    for i in flux_list:
        
        if(bg_subtraction==False):
            ax = fig.add_subplot(len(flux_list),1,n)
            ax = df_electrons_fluxes.iloc[:, i].plot(logy=True, figsize=(20,20), color='red')
        elif(bg_subtraction==True):
            ax = fig.add_subplot(len(flux_list),1,n)
            ax = df_electrons_fluxes_subtracted.iloc[:, i].plot(logy=True, figsize=(20,20), color='red')
    
        # Shaded area.
        ax.axvspan(bg_area[0], bg_area[1], color='gray', alpha=0.25)

        # Two vertical lines, or the "search area".
        ax.axvline(search_area[0], color='black')
        ax.axvline(search_area[1], color='black')

        # Finds peaks timestamp within search area and draws a vertical line.
        peak = df_electrons_fluxes.iloc[:, i].loc[search_area[0]:search_area[1]].idxmax()
        ax.axvline(peak, color='red')
        
        df_peak_electron_uncertainty = df_peak_electron_uncertainty.append(pd.DataFrame({'Energy channel':[i],
                                                                                         'Electron uncertainty':[df_electron_uncertainty.iloc[:, i].loc[peak]]}))
        
        # Appends peak data to peak dataframe.
        df_peaks = df_peaks.append(pd.DataFrame({'Energy channel':[i],
                                                 'Flux peak':[str(df_electrons_fluxes.iloc[:, i].loc[search_area[0]:search_area[1]].max())],
                                                 'Timestamp':[peak]}), ignore_index=False)

        #print("Flux peak: " + str(df_electrons_fluxes.iloc[:, i].loc[search_area[0]:search_area[1]].max()) + " at timestamp: " + \
        #     str(df_electrons_fluxes.iloc[:, i].loc[search_area[0]:search_area[1]].idxmax()))

        # Labels. X-axis label and ticks only visible for the last plot. One y-label for the whole figure, not for each plot separately.
        plt.xlabel("")
        ax.get_xaxis().set_visible(False)

        if(n==len(flux_list)):
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
            ax.xaxis.set_minor_locator(hours)
            plt.xlabel("Time")
            ax.get_xaxis().set_visible(True)

        # Adds the energy channels as text on top of plots.
        plt.text(0.025,0.55, electron_energies[i][0], transform=ax.transAxes)
        
        n=n+1
        
    df_peaks['Flux peak'] = pd.to_numeric(df_peaks['Flux peak'])
                                                         
    df_peaks_subtracted['Energy channel'] = df_peaks['Energy channel']
    df_peaks_subtracted['Subtracted flux peak'] = df_peaks['Flux peak'].sub(df_bg_fluxes['Background flux'])
    df_peaks_subtracted['Timestamp'] = df_peaks['Timestamp']
    
    
    # Spectrum plot error bar things.
    e_low = [0.0312, 0.0330, 0.0348, 0.0380, 0.0406, 0.0432, 0.0459, 0.0497, 0.0533, 0.0580, 0.0627, 0.0673, 0.0731, 0.0788, 0.0856, 0.0934, 0.1011, 0.1109, 0.1197, 0.1305, 0.1423, 0.1541, 0.1679, 0.1835, 0.1995, 0.2181, 0.2371, 0.2578, 0.2817, 0.3061, 0.3339, 0.3661, 0.3989, 0.4348]
    
    e_high = [0.0348, 0.0369, 0.0380, 0.0406, 0.0432, 0.0459, 0.0497, 0.0533, 0.0580, 0.0627, 0.0673, 0.0731, 0.0788, 0.0856, 0.0934, 0.1011, 0.1109, 0.1197, 0.1305, 0.1423, 0.1541, 0.1679, 0.1835, 0.1995, 0.2181, 0.2371, 0.2578, 0.2817, 0.3061, 0.3339, 0.3661, 0.3989, 0.4348, 0.4714]
    
    # Geometric mean
    #primary_energies = [np.sqrt(0.0312*0.0348), np.sqrt(0.033*0.0369), np.sqrt(0.0348*0.038), np.sqrt(0.038*0.0406), np.sqrt(0.0406*0.0432), np.sqrt(0.0432*0.0459), np.sqrt(0.0459*0.0497), np.sqrt(0.0497*0.0533), np.sqrt(0.0533*0.058), np.sqrt(0.058*0.0627), np.sqrt(0.0627*0.0673), np.sqrt(0.0673*0.0731), np.sqrt(0.0731*0.0788), np.sqrt(0.0788*0.0856), np.sqrt(0.0856*0.0934), np.sqrt(0.0934*0.1011), np.sqrt(0.1011*0.1109), np.sqrt(0.1109*0.1197), np.sqrt(0.1197*0.1305), np.sqrt(0.1305*0.1423), np.sqrt(0.1423*0.1541), np.sqrt(0.1541*0.1679), np.sqrt(0.1679*0.1835), np.sqrt(0.1835*0.1995), np.sqrt(0.1995*0.2181), np.sqrt(0.2181*0.2371), np.sqrt(0.2371*0.2578), np.sqrt(0.2578*0.2817), np.sqrt(0.2817*0.3061), np.sqrt(0.3061*0.3339), np.sqrt(0.3339*0.3661), np.sqrt(0.3661*0.3989), np.sqrt(0.3989*0.4348), np.sqrt(0.4348*0.4714)]
    
    primary_energies = []
    
    for i in range(0,len(e_low)):
        primary_energies.append(np.sqrt(e_low[i]*e_high[i]))
    
    df_peaks['Primary energy'] = primary_energies
    df_bg_fluxes['Primary energy'] = primary_energies
    df_peaks_subtracted['Primary energy'] = primary_energies
    
    # X errorbars.
    xerr_low = []
    xerr_high = []
    for i in range(0,len(primary_energies)):
        xerr_low.append(primary_energies[i]-e_low[i])
        xerr_high.append(e_high[i]-primary_energies[i])
    
    df_peaks['Flux peak'] = pd.to_numeric(df_peaks['Flux peak'])
    
    plt.savefig('/home/smurf/srl/images/energy channels-{}.jpg'.format(plotstart[:-5]), dpi=500)
    plt.show()
    
    return [df_electrons, plotstart, plotend, searchstart, searchend, bgstart, bgend, flux_list, df_peaks, df_bg_fluxes, df_peak_electron_uncertainty, df_average_bg_electron_uncertainty, xerr_low, xerr_high, primary_energies, df_peaks_subtracted]


"""
Function for previewing the plotting data. Mainly for checking that the data is correct before writing to a file.

Takes the same arguments as electron_flux_plot function. Recommended way is to store electron_flux_plot functions return as a variable and pass it on to this function.
"""
def preview_data(data_list):
    
    plotstart = data_list[1]
    plotend = data_list[2]
    searchstart = data_list[3]
    searchend = data_list[4]
    bgstart = data_list[5]
    bgend = data_list[6]
    df_peaks = data_list[8]
    
    print("--------------------------------------------------")
    print("Plot start: " + str(plotstart))
    print("Plot end: " + str(plotend))
    print("--------------------------------------------------")
    print("Search start: " + str(searchstart))
    print("Search end: " + str(searchend))
    print("--------------------------------------------------")
    print("Background start: " + str(bgstart))
    print("Background end: " + str(bgend))
    print("--------------------------------------------------")
            
    print(df_peaks)

    
    
    
"""
Function to write plotting data to a .csv file.

Takes the same arguments as electron_flux_plot function + a filepath (no filename).
"""
def write_to_csv(data_list, path):
    
    df_electrons = data_list[0]
    plotstart = data_list[1]
    plotend = data_list[2]
    search_area = [data_list[3],data_list[4]]
    bgstart = data_list[5]
    bgend = data_list[6]
    flux_list = data_list[7]
    df_peaks = data_list[8]
    df_bg_fluxes = data_list[9]
    df_peak_electron_uncertainty = data_list[10]
    df_average_bg_electron_uncertainty = data_list[11]
    xerr_low = data_list[12]
    xerr_high = data_list[13]
    primary_energies = data_list[14]
    
    df_peaks['Flux peak'] = pd.to_numeric(df_peaks['Flux peak'])
    df_bg_fluxes = df_bg_fluxes.apply(pd.to_numeric)
    
    df_to_write = pd.DataFrame({'Energy channel':[],
                                'Primary energy':[],
                                'Energy error low':[],
                                'Energy error high':[],
                                'Peak timestamp':[],
                                'Flux peak':[],
                                'Peak electron uncertainty':[],
                                'Background flux':[],
                                'Background electron uncertainty':[],
                                'Plot start':[],
                                'Plot end':[],
                                'Background start':[],
                                'Background end':[]})
    
    df_to_write['Energy channel'] = df_peaks['Energy channel']
    df_to_write['Primary energy'] = primary_energies
    df_to_write['Energy error low'] = xerr_low
    df_to_write['Energy error high'] = xerr_high
    df_to_write['Peak timestamp'] = df_peaks['Timestamp']
    df_to_write['Flux peak'] = df_peaks['Flux peak']
    df_to_write['Peak electron uncertainty'] = df_peak_electron_uncertainty['Electron uncertainty']
    df_to_write['Background flux'] = df_bg_fluxes['Background flux']
    df_to_write['Background electron uncertainty'] = df_average_bg_electron_uncertainty['Bg electron uncertainty']
    df_to_write['Plot start'] = plotstart
    df_to_write['Plot end'] = plotend
    df_to_write['Background start'] = bgstart
    df_to_write['Background end'] = bgend
    
    
    # Only specify the path, the filename is generated automatically.
    try:
        df_to_write.to_csv(path+"electron_data-"+str(plotstart)[:-5]+".csv", index=False)
        
    except FileNotFoundError:
        print("Not a valid file path. Subfolders have to be manually created for now.")
        
        
def spectrum_plot(data_list, bg_subtraction=False):
    
    df_bg_fluxes = data_list[9]
    df_peak_electron_uncertainty = data_list[10]
    df_average_bg_electron_uncertainty = data_list[11]
    xerr_low = data_list[12]
    xerr_high = data_list[13]
    df_peaks_subtracted = data_list[15]
    df_peaks = data_list[8]
    
    
    if(bg_subtraction):
        ax = df_peaks_subtracted.plot.scatter(x='Primary energy', y='Subtracted flux peak', c='red', label='Peaks')
        df_bg_fluxes.plot(kind='scatter', x='Primary energy', y='Background flux', c='red', alpha=0.25, ax=ax, label='Background')

        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Energy [MeV]')
        ax.set_ylabel('Flux \n [1/s cm$^2$ sr MeV]')

        # Observation errorbars.
        ax.errorbar(x=df_peaks_subtracted['Primary energy'], y=df_peaks_subtracted['Subtracted flux peak'], xerr=[xerr_low, xerr_high], yerr=df_peak_electron_uncertainty['Electron uncertainty'], fmt='.', ecolor='red', alpha=0.5)

        # Bg errorbars.
        ax.errorbar(x=df_bg_fluxes['Primary energy'], y=df_bg_fluxes['Background flux'], yerr=df_average_bg_electron_uncertainty['Bg electron uncertainty'], fmt='.', ecolor='red', alpha=0.15)
        
        plt.title('Flux peaks, ' + data_list[1][:-5] + '\n background subtracted')
        
        plt.savefig('/home/smurf/srl/images/spectrum-{}.jpg'.format(data_list[1][:-5]),dpi=500)
    
    else:
        ax = df_peaks.plot.scatter(x='Primary energy', y='Flux peak', c='red', label='Peaks')
        df_bg_fluxes.plot(kind='scatter', x='Primary energy', y='Background flux', c='red', alpha=0.25, ax=ax, label='Background')

        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Energy [MeV]')
        ax.set_ylabel('Flux \n [/s cm$^2$ sr MeV]')

        # Observation errorbars.
        ax.errorbar(x=df_peaks['Primary energy'], y=df_peaks['Flux peak'], xerr=[xerr_low, xerr_high], yerr=df_peak_electron_uncertainty['Electron uncertainty'], fmt='.', ecolor='red', alpha=0.5)

        # Bg errorbars.
        ax.errorbar(x=df_bg_fluxes['Primary energy'], y=df_bg_fluxes['Background flux'], yerr=df_average_bg_electron_uncertainty['Bg electron uncertainty'], fmt='.', ecolor='red', alpha=0.15)
        
        plt.title('Flux peaks, ' + str(data_list[1])[:-5] + '\n background subtraction off')
    
    plt.grid()
    plt.show()
    
    #print(df_peaks)
    
        
# Lines to make things faster.
# 2020-12-10
df_protons, df_electrons, energies = read_epd_cdf('ept', 'sun', 'l2', 20201210, 20201211, path='/home/smurf/solo/data/', autodownload=True)
data = electron_flux_plot(df_electrons, "2020-12-10-2000", "2020-12-11-0600", "2020-12-10-2300", "2020-12-11-0200", "2020-12-10-2100", "2020-12-10-2200")

# 2020-11-17
#df_protons, df_electrons, energies = read_epd_cdf('ept', 'sun', 'l2', 20201117, 20201118, path='/home/smurf/solo/data/', autodownload=True)
#data = electron_flux_plot(df_electrons, "2020-11-17-1600", "2020-11-17-2200", "2020-11-17-1800", "2020-11-17-2000", "2020-11-17-1630", "2020-11-17-1730")

# 2020-11-18
#df_protons, df_electrons, energies = read_epd_cdf('ept', 'sun', 'l2', 20201118, 20201118, path='/home/smurf/solo/data/', autodownload=True)
#data = electron_flux_plot(df_electrons, "2020-11-18-1100", "2020-11-18-1800", "2020-11-18-1300", "2020-11-18-1600", "2020-11-18-1130", "2020-11-18-1230")

# 2021-04-17
#df_protons, df_electrons, energies = read_epd_cdf('ept', 'sun', 'l2', 20210417, 20210418, path='/home/smurf/solo/data/', autodownload=True)
#data = electron_flux_plot(df_electrons, "2021-04-17-1500", "2021-04-17-2100", "2021-04-17-1700", "2021-04-17-1900", "2021-04-17-1530", "2021-04-17-1630")



# to-do
# - bg electron uncertainty implementation is not very good.