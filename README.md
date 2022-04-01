## **PSD Analyser: A set of Python/MATLAB tools for particle size distribution (PSD) analysis and visualisation**

### Hugh Patrick Rice, January 2021 onwards

Zenodo link here

## Quick start user guide

**Run the standalone executable (currently Windows only)**
- Go to the *Releases* on the right of the right of the page
- Download the executable file (ending *.exe*) and any of the spreadsheet example files from the main page

**Run the Python or MATLAB code directly**
- Clone the repository or download files individually
- Import into your Python environment or import into MATLAB
- All dependencies are listed in the *Dependencies* document

## Detailed guide with examples

*Summary of PSD Analyser functionality*
1. Parsing of *Mastersizer* files in spreadsheet (CSV, Excel) formats
2. Log-normal modelling of PSDs
3. Application to PSDs of the product difference algorithm (PDA) to PSDs, which computes the discrete distribution of *N* elements with the same statistical moments (mean, etc.)
4. Visualisation of results in interactive viewer, allowing output of figures in various formats (currently Python only)

*Note on history and previous implementation*
- The concepts behind *PSD Analyser* were developed during research that led to the production of references (1) and (2). Please refer to those for more details and examples.

## Contact and issue reporting

Please either raise an issue here at Github or contact me directly.

*Contact:* Hugh Rice, h.p.rice@leeds.ac.uk

## How to cite this repository

- Copy or click the Zenodo link above, which has a corresponding DOI attached, and construct your own citation that contains it
- Depending on your style, your citation should look something like this: Rice HP (2022), *PSD Analyser: A set of Python/MATLAB tools for particle size distribution (PSD) analysis and visualisation*, Github code repository, DOI: 
- If you're unsure, please contact me

## References

(1) Rice HP, Fairweather M, Peakall J, Hunter TN, Mahmoud B and Biggs SR (2015), *Constraints on the functional form of the critical deposition velocity in solid–liquid pipe flow at low solid volume fractions*, Chemical Engineering Science **126** 759-770, DOI: https://doi.org/10.1016/j.ces.2014.12.039

(2) Rice HP, Peakall J, Fairweather M and Hunter TN (2020), *Extending estimation of the critical deposition velocity in solid–liquid pipe flow to ideal and non-ideal particles at low and intermediate solid volume fractions*, Chemical Engineering Science **211** 115308 (9 p.), DOI: https://doi.org/10.1016/j.ces.2019.115308

(3) McGraw R (1997), *Description of Aerosol Dynamics by the Quadrature Method of Moments*, Aerosol Science and Technology **27** (2) 255-265, DOI: https://doi.org/10.1080/02786829708965471

(4) Marchisio DL, Vigil RD and Fox RO (2003), *Implementation of the quadrature method of moments in CFD codes for aggregation–breakage problems*, Chemical Engineering Science **58** (15) 3337-3351, DOI: https://doi.org/10.1016/S0009-2509(03)00211-2

(5) Gordon RG (1968), *Error Bounds in Equilibrium Statistical Mechanics*, Journal of Mathematical Physics **9** 655-663, DOI: https://doi.org/10.1063/1.1664624

(6) Wheeler and Gordon (1971), *Bounds for averages using moment constraints*, In: Baker and Gammel (eds.), *The Padé Approximant in Theoretical Physics*, New York and London: Elsevier Science, ISBN: 9780080955803

(7) Farr SF (2013), *Random close packing fractions of lognormal distributions of hard spheres*, Powder Technology **245** 28-34, DOI: https://doi.org/10.1016/j.powtec.2013.04.009
