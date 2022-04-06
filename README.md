## **PSD Analyser: A set of Python/MATLAB tools for particle size distribution (PSD) analysis and visualisation**

### Hugh Patrick Rice, 2022

[![DOI](https://zenodo.org/badge/328645883.svg)](https://zenodo.org/badge/latestdoi/328645883)

## Quick-start guide

**Run the standalone executable (currently Windows only)**
- Go to the *Releases* on the right of the page
- Download the executable file (ending *.exe*) and any of the spreadsheet example files from the main page
- Once downloaded, double-click the executable file and it will run

**Run the code directly**
- Clone the repository or download files individually
- (Python) Import into your Python environment
  - Ensure all the necessary libraries (see imports in *PSD-Analyser* an *psd_analyser*) are installed in your environment
  - Instructions for creating a miminal environment in which *PSDAnalyser* can be run is provided in *env_setup.txt*, using the *Anaconda* distribution of Python (details [here](https://www.anaconda.com/)) or the standard Python environment management library *venv*
  - There are several examples given in the *psd_analyser* script to demonstrate functionality, specifically: loading and saving data, plotting PSDs, fitting PSDs to log-normal distributions; and reducing them to *N* components using the product difference algorithm (PDA); comment out the code as appropriate
- (MATLAB) Import each function and class individually (**to be completed**)
  - Download and call the functions individually (more functionality under development)

## Detailed guide

**How *PSD Analyser* works (Python)**

The app (*PSD-Analyser*) creates a user interface for loading, interacting with and saving PSD data. However, the hard work is done by the *PSDAnalyser* class within the *psd_analyser* library, which can be used without the app and contains some functionality not available in the app (see below).

**Summary of *PSD Analyser* functionality**

1. Parsing of *Mastersizer* files in spreadsheet (CSV, Excel) formats
2. Log-normal modelling of PSDs, where fitting can be performed to eithe the cumulative distribution function (CDF) or probability density function (PDF) of the PSD data
3. Application to PSDs of the product difference algorithm (PDA) to PSDs, which computes the discrete distribution of *N* elements with the same statistical moments (mean, etc.)
4. (With standalone Python app) Visualisation of results in an interactive viewer, allowing output of figures in various formats (currently Python only)
5. Saving PSD data, with derived log-normal models, into a single spreadsheet for further analysis

**Running the *PSD Analyser* standalone app (Python)**

- See the quick start guide above to get started
- Once *PSD Analyser* is running, load one of the spreadsheets provided as examples; these are all outputs from the *Mastersizer* laser-diffraction-based particle sizing device manufactured by *Malvern Panalytical* (formerly *Malvern Instruments*)
- A summary of the user interface is provided below
  - *Load file* Select a spreadsheet from which to load PSD data; several example files are provided
  - *Save data* Dump all loaded and modelled data to a single spreadsheet; by default the file format is the same as the original load-file, but the user can specify otherwise; the index of each dataset corresponds to the row in which it appeared in the original data file
  - *Export figure* Dump current plot to an image file for later use; user can specify image format
  - *Forward and back arrows* If the file loaded contains multiple PSD datasets, these arrow buttons scroll through them; the dataset number is shown in the information panel (bottom left)
  - *Fit to PDF/CDF buttons* These buttons toggle between the two methods for fitting of a log-normal distribution to the dataset displayed; see code for more information, and note that each can give different results for the log-normal fit parameters
  - *Plot PDF/CDF buttons* These buttons toggle between plotting the loaded and log-normal fitted data in CDF and PDF form, depending on the user's preference
  - *Plot log/linear buttons* Toggle between logarithmic and linear scale on the *x*-axis (*i.e.* particle size)

**Using *PSD Analyser* in code form (Python)**

- Running the code directly provides additional flexibility not available in the app, in particular:
  1. The product difference algorithm (PDA) can be used to model a given PSD as another with an arbitrary number of elements and the same moments, *N*; this is intended as a tool in computationally expensive applications that rely on size-fraction-dependent calculations
  2. Log-normal fitting can be executed with or without pre-fitting (pre-fitting is used by default in the app) *via* linear regression using a linearised version of the equation for the log-normal CDF; pre-fitting has a computational cost but makes it more likely that a subsequent non-linear fit will be successful

**Creating your own standalone app using the *PSD Analyser* source code (Python)**

The app was created using the Python library *pyinstaller*, and you can do the same. A rough outline of the necessary code for doing so is below.

1. In the command line in your Python environment, install *pyinstaller* (see [here](https://pypi.org/project/pyinstaller/)):

```
pip install pyinstaller
```

2. Then navigate to the directory containing your Python scripts (*i.e.* *PSD-Analyser.py* and *psd_analyser.py*) and create the executables, where: the ```noconsole``` option creates an executable without a console window (this option can be removed if you wish to use the console window for debugging); and the ```onefile``` option creates a single executable file (if this option is removed, the resulting executable may be smaller but a number of directories will be created alongside the executable)

```
pyinstaller --onefile --noconsole PSD-Analyser.py
```

**Important equations used by *PSD Analyser***

1. The log-normal cumulative distribution function (CDF). See [here](https://en.wikipedia.org/wiki/Log-normal_distribution#Cumulative_distribution_function) and [here](https://mathworld.wolfram.com/LogNormalDistribution.html).
2. The log-normal probability density function (PDF). See [here](https://en.wikipedia.org/wiki/Log-normal_distribution#Probability_density_function) and [here](https://mathworld.wolfram.com/LogNormalDistribution.html).
3. The product difference algorithm (PDA). See references below.

## Contact and issue reporting

Please either raise an issue here at Github or contact me directly.

*Contact:* Hugh Rice, h.p.rice@leeds.ac.uk

## How to cite this repository

- Copy or click the Zenodo link above, which has a corresponding DOI attached, and construct your own citation that contains it
- Depending on your style, your citation should look something like this: Rice HP (2022), *PSD Analyser: A set of Python/MATLAB tools for particle size distribution (PSD) analysis and visualisation*, Github code repository, DOI: 10.5281/zenodo.6408352
- If you're unsure, please contact me

## References and notes

*Note on history, development and previous implementation*
- The concepts behind *PSD Analyser* were developed during research that led to the production of references (1) and (2). Please refer to those for more details and examples.
- The standalone app was constructed using the Python GUI library *Tkinter* and relies heavily on standard Python numerical and mathematical libraries *scipy*, *numpy* and *pandas*, as well as *matplotlib* for visualisation.

(1) Rice HP, Fairweather M, Peakall J, Hunter TN, Mahmoud B and Biggs SR (2015), *Constraints on the functional form of the critical deposition velocity in solid–liquid pipe flow at low solid volume fractions*, Chemical Engineering Science **126** 759-770, DOI: https://doi.org/10.1016/j.ces.2014.12.039

(2) Rice HP, Peakall J, Fairweather M and Hunter TN (2020), *Extending estimation of the critical deposition velocity in solid–liquid pipe flow to ideal and non-ideal particles at low and intermediate solid volume fractions*, Chemical Engineering Science **211** 115308 (9 p.), DOI: https://doi.org/10.1016/j.ces.2019.115308

(3) McGraw R (1997), *Description of Aerosol Dynamics by the Quadrature Method of Moments*, Aerosol Science and Technology **27** (2) 255-265, DOI: https://doi.org/10.1080/02786829708965471

(4) Marchisio DL, Vigil RD and Fox RO (2003), *Implementation of the quadrature method of moments in CFD codes for aggregation–breakage problems*, Chemical Engineering Science **58** (15) 3337-3351, DOI: https://doi.org/10.1016/S0009-2509(03)00211-2

(5) Gordon RG (1968), *Error Bounds in Equilibrium Statistical Mechanics*, Journal of Mathematical Physics **9** 655-663, DOI: https://doi.org/10.1063/1.1664624

(6) Wheeler and Gordon (1971), *Bounds for averages using moment constraints*, In: Baker and Gammel (eds.), *The Padé Approximant in Theoretical Physics*, New York and London: Elsevier Science, ISBN: 9780080955803

(7) Farr SF (2013), *Random close packing fractions of lognormal distributions of hard spheres*, Powder Technology **245** 28-34, DOI: https://doi.org/10.1016/j.powtec.2013.04.009
