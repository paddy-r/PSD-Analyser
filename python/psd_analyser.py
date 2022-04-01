''' HR 08/02/22 Particle size distribution (PSD) analyser
    To:     (1) Parse Mastersizer-derived spreadsheet and extract PSD data,
            (2) Fit PSD data to log-normal distribution, and
            (3) Apply product difference algorithm to PSD data to give n-component distribution with same moments
    References in methods '''

import pandas as pd
# from numpy import diff
import os
import matplotlib as mpl
# from scipy.stats import lognorm
from scipy.stats import linregress
from scipy.special import erf, erfinv
from scipy.optimize import curve_fit
import numpy as np
# from numpy import exp, sqrt, log, pi, linspace, diff, zeros, diag, cos
from numpy.linalg import eig
# from matplotlib.figure import Figure
import matplotlib as mpl



''' Returns value of lognormal CDF at x
    with lognormal parameters M, S '''
def lognormal_CDF(x, M, S, called_by = None):
    if called_by:
        print('Called by: ', called_by)
    # print('LN_CDF: ', x, M, S)
    CDF = 0.5*(1.0+erf((np.log(x)-M)/(S*np.sqrt(2.0))))
    return CDF

''' Returns value of (normalised) lognormal PDF at x
    with lognormal parameters M, S '''
def lognormal_PDF(x, M, S):
    PDF = (1.0/(S*np.sqrt(2.0*np.pi)*x))*np.exp(-((np.log(x)-M)**2)/(2*(S**2)))
    return PDF

''' To fit to lognormal using normalised PDF '''
def fit_lognormal_PDF(d, P_norm, p0 = None):
    results = curve_fit(lognormal_PDF, d, P_norm, p0)
    M,S = results[0]
    return M,S

''' To fit to lognormal using CDF '''
def fit_lognormal_CDF(d, C, p0 = None):
    results = curve_fit(lognormal_CDF, d, C, p0)
    M,S = results[0]
    return M,S

''' To fit to lognormal using linearised CDF '''
def fit_lognormal_CDF_linear(d, C):
    x = [erfinv(2*el - 1.0) for el in C]
    y = [np.log(el) for el in d]
    # print(x,y)
    results = linregress(x, y)
    M = results[1]
    S = results[0]/np.sqrt(2.0)
    # R = results[2]
    # print('Linear regression R = ', R)
    return M,S

def mu(M, S):
    _mu = np.exp(M + S**2)
    return _mu

def sigma(M, S):
    _sigma = np.sqrt(np.exp(S**2 + 2.0*M)*(np.exp(S**2)-1.0))
    return _sigma

''' Returns x value corresponding to pth quantile (0 < p < 100)
    of lognormal distribution with parameters M, S '''
def lognormal_quantile(M, S, p):
    x_p = np.exp(M + np.sqrt(2*S**2)*erfinv(2*(p/100)-1.0))
    return x_p

''' Calculate nth or (n/m)th moment of distribution
    f = frequencies at discrete values of x
    m = 0 by default, i.e. denominator is unity '''
def moment(x, f, n, m = 0):
    M1 = sum([(a**n) * b for a,b in zip(x,f)])
    if m == 0:
        M2 = 1.0
    else:
        M2 = sum([(a**m) * b for a,b in zip(x,f)])
    result = M1/M2
    return result

''' Get nth or (n/m)th moment of lognormal distribution
    m = 0 by default, i.e. denominator is unity '''
def lognormal_moment(M, S, n, m = 0):
    M1 = np.exp(n*M) + 0.5*n**2*S**2
    if m == 0:
        M2 = 1.0
    else:
        M2 = np.exp(m*M) + 0.5*m**2*S**2
    result = M1/M2
    return result

''' To compute moment from PSD (raw or from PDA) and lognormal parameters '''
def moments_compare(x, P, M, S, n, m):
    ''' Get moments from PDF and lognormal fit '''
    m_PSD = moment(x,P,n,m)
    m_lognormal = lognormal_moment(M,S,n,m)
    moments = (m_PSD, m_lognormal)
    return moments

''' To use product difference algorithm to create new PSD with N arbitrary elements
    having same moments as input PSD
    - McGraw (1997), DOI: https://doi.org/10.1080/02786829708965471
    - Marchisio et al. (2003), DOI: https://doi.org/10.1016/S0009-2509(03)00211-2
    - Gordon (1968), DOI: https://doi.org/10.1063/1.1664624
    - Wheeler and Gordon (1971), Bounds for averages using moment constraints,
        In: Baker and Gammel (eds.), The Padé Approximant in Theoretical Physics,
        New York and London: Elsevier Science. ISBN: 9780080955803 '''
def product_difference_algorithm(PSD, N = None):
    ''' Check if binned; if so, unbin, etc. '''
    x = PSD[0]
    P = PSD[1]
    if len(x) != len(P):
        # Do something here
        pass

    ''' Default to 3 for output PSD '''
    if not N:
        N = 3

    ''' Get moments '''
    m = []
    for k in range(2*N):
        mk = moment(x,P,k)
        m.append(mk)

    ''' Populate b array with zeros '''
    B = np.zeros([2*N + 1,2*N + 1])
    ''' First column; only first element is non-zero '''
    B[0,0] = 1
    ''' Second column '''
    for i in range(2*N):
        B[i,1] = (-1)**i * m[i];
    ''' Third column onwards '''
    for j in range(2,2*N+1):
        for i in range(2*N+1-j):
            B[i,j] = B[0,j-1]*B[i+1,j-2] - B[0,j-2]*B[i+1,j-1]

    ''' Populate alpha; first element is zero '''
    alpha = [0]
    for i in range(1,2*N):
        alpha_i = B[0,i+1]/(B[0,i]*B[0,i-1])
        alpha.append(alpha_i)

    ''' Populate a '''
    a = []
    for i in range(N):
        ai = alpha[2*i+1] + alpha[2*i]
        a.append(ai)

    ''' Populate b '''
    b = []
    for i in range(N-1):
        bi = -np.sqrt(alpha[2*i+2]*alpha[2*i+1])
        b.append(bi)

    ''' Get result:
        1. Eigenvectors are weights, w
        2. Eigenvalues are values, d '''
    J = np.diag(b,-1) + np.diag(a) + np.diag(b,1)
    x_new,v = eig(J)
    # P_new = [m[0]*el[0]**2 for el in v.transpose()]
    P_new = [m[0]*el**2 for el in v[0]]

    ''' Create new PSD '''
    PSD_new = (x_new,P_new)
    return PSD_new

''' Estimate packing fraction of particle species;
    assumes non-interacting, spherical particles with lognormal sizes
    - Farr (2013), DOI: https://doi.org/10.1016/j.powtec.2013.04.009 '''
def packing_fraction(M, S):
    pf = 1.0 - 0.57*np.exp(-S) + 0.2135*np.exp(-0.57*S/0.2135) + 0.0019*(np.cos(2.0*np.pi*(1 - np.exp(-0.75*S**(0.7) - 0.025*S**4))) - 1)
    return pf

def load_csv(file, *args, **kwargs):
    loaded_data = pd.read_csv(file, *args, **kwargs)
    return loaded_data

def load_excel(file, *args, **kwargs):
    loaded_data = pd.read_excel(file, *args, **kwargs)
    return loaded_data

''' Plot CDF, optionally with fit curve '''
def plot_CDF(x, CDF, ax, fit_data = None, xlabel = None, ylabel = None):
    ax.plot(x, CDF, 'kx')
    if fit_data:
        M = fit_data[0]
        S = fit_data[1]
        # spread = max(x) - min(x)
        ''' Setting to zero as otherwise borks plotting,
            i.e. log(-number) if fit_data specified '''
        spread = 0
        x_fit = np.linspace(min(x)-0.05*spread, max(x)+0.05*spread)
        y_fit = lognormal_CDF(x_fit, M, S, called_by = 'plot_CDF')
        ax.plot(x_fit, y_fit, 'k--')
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

''' Plot PDF, optionally with fit curve '''
def plot_PDF(x, PDF, ax, fit_data = None, xlabel = None, ylabel = None):
    ax.plot(x, PDF, 'kx')
    if fit_data:
        M = fit_data[0]
        S = fit_data[1]
        # spread = max(x) - min(x)
        ''' Setting to zero as otherwise borks plotting,
            i.e. log(-number) if fit_data specified '''
        spread = 0
        x_fit = np.linspace(min(x)-0.05*spread, max(x)+0.05*spread)
        y_fit = lognormal_PDF(x_fit, M, S)
        ax.plot(x_fit, y_fit, 'k--')
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

def get_valids(x_raw, C_raw, P_raw):
    valids = [i for i,el in enumerate(P_raw) if (el > 0.0) and C_raw[i] < 1.0]
    # print('nonzeros: \n', valids)
    x = [x_raw[i] for i in valids]
    C = [C_raw[i] for i in valids]
    P = [P_raw[i] for i in valids]
    return x,C,P


class PSDAnalyser():

    ''' Constructor '''
    def __init__(self, ax = None, *args, **kwargs):
        self.suffix_dict = {'.csv': load_csv,
                            '.xls': load_excel,
                            '.xlsx': load_excel}
        self._loaded = False
        self._file = None
        self.SAVE_TEXT_DEFAULT = "_computed"
        self.IMAGE_DUMP_DEFAULT = "plot"
        self.IMAGE_DUMP_EXT_DEFAULT = ".jpg"
        
        ''' Default fit and plot modes (True) are CDF '''
        # self.fit_mode_dict = {True: "Fit by CDF", False: "Fit by PDF"}
        self.FIT_MODE_DEFAULT = True
        self.FIT_MODE_PREFIT_DEFAULT = True
        self.PLOT_MODE_DEFAULT = True

        self.COLUMN_DICT_DEFAULT = {}

        ''' Set up figure and plot axes '''
        if ax:
            self.ax = ax
            self.fig = self.ax.get_figure()
        if not ax:
            self.fig = mpl.figure.Figure()
            self.ax = self.fig.add_subplot()

    ''' Raw data to DataFrame '''
    def load_spreadsheet(self, file):
        file_ext = os.path.splitext(file)[-1]
        if file_ext in self.suffix_dict:
            try:
                self.loaded_data = self.suffix_dict[file_ext](file)
                print('Loaded data from file: ', file)
                self._loaded = True
                self._file = file
            except Exception as e:
                print('Could not load data from file: ', file)
                print('Exception: ', e)
                return None
        else:
            print('File type not supported; could not load data')
            return None

        # return self.loaded_data

    # def get_all_size_data(self, loaded_data, *args, **kwargs):
    #     ''' Check data loaded from file '''
    #     if not (self._file and self._loaded):
    #         print('No file loaded; cannot get size data')
    #         return None

        ld = self.loaded_data

        start_text = 'Result Between'
        end_text = 'Operator Notes'
        self.start_column_index, self.start_column = [(i,el) for i,el in enumerate(ld) if start_text in el][0]
        self.end_column_index, self.end_column = [(i,el) for i,el in enumerate(ld) if end_text in el][0]
        self.start_rows = [0]
        for i,el in enumerate(ld[self.start_column]):
            if type(el) == str:
                self.start_rows.append(i)

        ''' Create dict containing bins and sizing data rows corresponding to them '''
        self.bin_groups = {0:[]}
        for row in range(len(ld.index)):
            if ld.loc[row][0] == "Sample Name":
                self.bin_groups[row + 1] = []

        ''' Populate bin groups '''
        bin_list = list(self.bin_groups)
        bin_list.sort()
        bin_list.append(len(ld.index) + 1)
        for i,group in enumerate(bin_list[:-1]):
            for row in range(group + 1, bin_list[i+1] - 1):
                if type(ld.loc[row][0]) == float:
                    ''' Exclude if "nan", i.e. empty float '''
                    if np.isnan(ld.loc[row][0]):
                        continue
                    # row_name = str(ld.loc[row][0])
                # elif type(ld.loc[row][0]) == str:
                #     row_name = ld.loc[row][0]
                # self.bin_groups[group].append((row, row_name))
                self.bin_groups[group].append(row)

        ''' Create datasets for each data row '''
        self.datasets = {}
        for group in self.bin_groups:
            for row in self.bin_groups[group]:
                dataset = self.get_dataset(row)
                ''' Add to dataset dict if dataset returned correctly '''
                if dataset:
                    self.datasets[row] = dataset

    ''' Return dataset by row; bin group optional '''
    def get_dataset(self, row, group = None, fit_mode = None, prefit = None, column_dict = None):
        ''' Check if row has corresponding PSD dataset, if none specified;
            else check correct group has been given, i.e. row is in group '''
        # print('Getting data for row ', row)
        if not group:
            try:
                bin_group = [k for k,v in self.bin_groups.items() if row in v][0]
                # print('Found bin group (row', bin_group, ')')
            except:
                # print('Row does not have PSD dataset; returning None')
                return None
        else:
            if not row in self.bin_groups[group]:
                print('Row not present in specified bin group; returning None')
                return None

        ''' Check fit type '''
        if (fit_mode is None) or (fit_mode not in (True,False)):
            fit_mode = self.FIT_MODE_DEFAULT
        # fit_by = self.fit_mode_dict[fit_mode]
        fit_by_CDF = fit_mode

        ''' Create dictionary for output '''
        dataset = {}
        if not column_dict:
            column_dict = self.COLUMN_DICT_DEFAULT
        ld = self.loaded_data

        ''' Get whole rows '''
        row_bins = list(ld.loc[bin_group])
        row_data = list(ld.loc[row])
        # print(row_bins)
        # print('row data:\n', list(row_data))

        ''' Get bins, PDF, CDF, etc. '''
        bins = row_bins[self.start_column_index: self.end_column_index]
        # print('bins:\n', bins)
        bin_widths = np.diff(bins)
        bin_centres = [np.sqrt(bins[i] * bins[i+1]) for i in range(0, len(bins)-1)]
        # print('bin centres:\n', bin_centres)
        PDF = row_data[self.start_column_index + 1: self.end_column_index]
        PDF = [el/100 for el in PDF]
        ''' Must normalise, as PDF is strictly density function
            N.B. Mastersizer gives percentages for raw PDF so sum(PDF) is ~1
            but this isn't true in general for all sources of data,
            so sum(PDF) is retained here '''
        sum_PDF = sum(PDF)
        PDF_normalised = [a/(b*sum_PDF) for a,b in zip(PDF, bin_widths)]
        # print('PDF:\n', PDF)
        CDF = np.cumsum(PDF)
        # print('CDF:\n', CDF)

        ''' Get all constant-like data '''
        sample_name = ld['Sample Name'][row]
        # print('Sample name: ', sample_name)
        D43 = ld['D [4, 3] - Volume weighted mean'][row]
        # print('D43: ', D43)
        D32 = ld['D [3, 2] - Surface weighted mean '][row]
        # print('D32: ', D32)
        d10 = ld['d (0.1)'][row]
        # print('d10: ', d10)
        d50 = ld['d (0.5)'][row]
        # print('d50: ', d50)
        d90 = ld['d (0.9)'][row]
        # print('d90: ', d90)

        ''' Filter by non-zero PDF entries to avoid numerical problems with fitting '''
        # valids = [i for i,el in enumerate(PDF_normalised) if (el > 0.0) and CDF[i] < 1.0]
        # # print('nonzeros: \n', valids)
        # x = [bin_centres[i] for i in valids]
        # C = [CDF[i] for i in valids]
        # P = [PDF_normalised[i] for i in valids]

        x,C,P = get_valids(bin_centres, CDF, PDF_normalised)

        ''' Get all fitted data '''
        if not prefit:
            prefit = self.FIT_MODE_PREFIT_DEFAULT
        if prefit:
            # print('x,c:\n',x,C)
            p0 = fit_lognormal_CDF_linear(x, C)
        else:
            p0 = None

        if (fit_mode is None) or (fit_mode not in (True,False)):
            # print('Fit mode not given; defaulting')
            fit_mode = self.FIT_MODE_DEFAULT
        if fit_mode:
            # print('Fitting by CDF')
            # print('x,c:\n',x,C)
            M,S = fit_lognormal_CDF(x, C, p0 = p0)
        else:
            # print('Fitting by PDF')
            M,S = fit_lognormal_PDF(x, P, p0 = p0)
        # print('M,S:\n', M, S)
        D43_fit = lognormal_moment(M,S,4,3)
        D32_fit = lognormal_moment(M,S,3,2)
        d10_fit = lognormal_quantile(M,S,10)
        d50_fit = lognormal_quantile(M,S,50)
        d90_fit = lognormal_quantile(M,S,90)
        phi = packing_fraction(M,S)

        ''' Populate dataset '''
        data_list = ['sample_name',
                     'D43',
                     'D32',
                     'd10',
                     'd50',
                     'd90',
                     'bins',
                     'bin_widths',
                     'bin_centres',
                     'PDF',
                     'PDF_normalised',
                     'CDF',
                     'M',
                     'S',
                     'D43_fit',
                     'D32_fit',
                     'd10_fit',
                     'd50_fit',
                     'd90_fit',
                     'phi',
                     'fit_by_CDF']
        # dataset['sample_name'] = sample_name
        # dataset['D43'] = D43
        # dataset['D32'] = D32
        # dataset['d10'] = d10
        # dataset['d50'] = d50
        # dataset['d90'] = d90
        # dataset['bins'] = bins
        # dataset['bin_widths'] = bin_widths
        # dataset['bin_centres'] = bin_centres
        # dataset['PDF'] = PDF
        # dataset['PDF_normalised'] = PDF_normalised
        # dataset['CDF'] = CDF
        # dataset['M'] = M
        # dataset['S'] = S
        # dataset['D43_fit'] = D43_fit
        # dataset['D32_fit'] = D32_fit
        # dataset['d10_fit'] = d10_fit
        # dataset['d50_fit'] = d50_fit
        # dataset['d90_fit'] = d90_fit
        # dataset['phi'] = phi
        # dataset['fit_by_CDF'] = fit_by_CDF

        ''' Alternate method for populating "dataset" dict: by "eval",
            instead of by manual addition, as above '''
        for item in data_list:
            dataset[item] = eval(item)

        return dataset

    ''' To toggle dataset between CDF and PDF fit modes; and prefit if specified '''
    def toggle_fit_mode(self, row, prefit = True):
        if not hasattr(self, 'datasets'):
            print('Datasets not initialised; aborting "toggle_fit_mode"')
            return

        ''' Get and check all data from dataset '''
        ds = self.datasets[row]
        x_all = ds['bin_centres']
        C_all = ds['CDF']
        P_all = ds['PDF_normalised']
        x,C,P = get_valids(x_all, C_all, P_all)

        ''' Prefit if specified '''
        if prefit:
            p0 = fit_lognormal_CDF_linear(x, C)
        else:
            p0 = None

        ''' Get new fit results and calculate new derived values'''
        if ds['fit_by_CDF'] == True:
            ds['fit_by_CDF'] = False
            M,S = fit_lognormal_PDF(x, P, p0 = p0)
            # print('M,S after toggle: ', M,S)
        else:
            ds['fit_by_CDF'] = True
            M,S = fit_lognormal_CDF(x, C, p0 = p0)
            # print('M,S after toggle: ', M,S)
        
        D43_fit = lognormal_moment(M,S,4,3)
        D32_fit = lognormal_moment(M,S,3,2)
        d10_fit = lognormal_quantile(M,S,10)
        d50_fit = lognormal_quantile(M,S,50)
        d90_fit = lognormal_quantile(M,S,90)
        phi = packing_fraction(M,S)

        ''' Update dataset with new values from fit '''
        ds['M'] = M
        ds['S'] = S
        ds['D43_fit'] = D43_fit
        ds['D32_fit'] = D32_fit
        ds['d10_fit'] = d10_fit
        ds['d50_fit'] = d50_fit
        ds['d90_fit'] = d90_fit
        ds['phi'] = phi

    def export_to_spreadsheet(self, file = None):
        print('Exporting all original and processed PSD data; in original format unless specified')
        if not file:
            stub,ext = os.path.splitext(self._file)
            file = stub + self.SAVE_TEXT_DEFAULT + ext
        print('-> File full path: ', file)
        ''' Get original and computed data, to go in separate sheets
            - Original data are exactly as in loaded file
            - Computed data are all PSDs, lognormal params, etc., computed from original data '''
        original_data = self.loaded_data
        computed_data = pd.DataFrame()
        try:
            with pd.ExcelWriter(file) as writer:  
                original_data.to_excel(writer, sheet_name = 'Original data')
                computed_data.to_excel(writer, sheet_name = 'Computed data')
            return True
        except Exception as e:
            print('Could not save file; exception follows...')
            print(e)
            return False

    ''' Dump current figure to image file '''
    def dump_figure(self, file = None):
        ''' Create filename if not specified based on defaults (JPG) '''
        if not file:
            cwd = os.get_cwd()
            file_stub = self.IMAGE_DUMP_DEFAULT
            if self._CDF:
                plot_type = "_CDF"
            else:
                plot_type = "_PDF"
            file_ext = self.IMAGE_DUMP_EXT_DEFAULT
            file = file_stub + plot_type + file_ext
            file = os.path.join(cwd, file)

        try:
            print('Saving image file; full path: ', file)
            self.fig.savefig(file)
            print('Image file saved')
        except Exception as e:
            print('Could not save file; exception follows')
            print(e)



if __name__ == "__main__":

    ''' Basic calculations testing; PSDs from old MATLAB scripts '''
    d_h22 = ["7.585776    8.709636    10    11.481536    13.182567    15.135612    17.378008    19.952623    22.908677    26.30268    30.199517    34.673685    39.810717    45.708819    52.480746    60.255959    69.183097    79.432823",
              "0.007635    0.135005    0.310442    0.475434    0.558515    0.617157    0.910635    1.939477    4.28946    8.158915    12.980964    17.173874    18.739872    16.548617    11.479466    5.072049    0.602481",
              "0.007635    0.14264    0.453082    0.928516    1.487031    2.104188    3.014823    4.9543    9.24376    17.402675    30.383639    47.557513    66.297385    82.846002    94.325468    99.397517    99.999998"]

    d_h16 = ["30.199517    34.673685    39.810717    45.708819    52.480746    60.255959    69.183097    79.432823    91.201084    104.712855    120.226443    138.038426    158.489319",
              "0.012362    0.42421    2.158819    6.105918    12.005703    17.790159    20.449133    18.41713    12.907406    6.753966    2.47866    0.46569",
              "0.012362    0.436572    2.595391    8.701309    20.707012    38.497171    58.946304    77.363434    90.27084    97.024806    99.503466    99.969156"]

    raw_data = d_h16
    # raw_data = d_h22

    bins = [float(el) for el in raw_data[0].split()]
    PDF = [float(el)/100 for el in raw_data[1].split()]
    CDF = [float(el)/100 for el in raw_data[2].split()]

    bin_widths = np.diff(bins)
    # x = [(bins[i] + bins[i+1])/2 for i in range(0, len(bins)-1)]
    x = [np.sqrt(bins[i] * bins[i+1]) for i in range(0, len(bins)-1)]

    sum_PDF = sum(PDF)
    PDF_normalised = [a/(b*sum_PDF) for a,b in zip(PDF,bin_widths)]

    ''' Get estimate from linearised version '''
    M,S = fit_lognormal_CDF_linear(x, CDF)
    ''' ...then use estimates for nonlinear fit '''
    print('Estimates of M,S from linearised fit: ', M,S)
    # M,S = fit_lognormal_CDF(x, CDF, p0 = [M,S])
    M,S = fit_lognormal_PDF(x, PDF_normalised, p0 = [M,S])
    print('Final results for M,S from nonlinear fit: ', M,S)

    ''' Plotting; simplest possible MPL plot axes! '''
    ax = mpl.pyplot.axes()

    plot_CDF(x, CDF, ax, (M,S),
             xlabel = "Particle diameter ($\mu$m)",
             ylabel = 'CDF')
    plot_PDF(x, PDF_normalised, ax, (M,S),
             xlabel = "Particle diameter ($\mu$m)",
             ylabel = 'PDF (normalised)')

    ''' File load testing '''
    ps = PSDAnalyser()
    file = 'C:\\Users\\hughr\\MATLAB Drive\\PSD-Analyser\\MS2000 Export.xlsx'
    ps.load_spreadsheet(file)
