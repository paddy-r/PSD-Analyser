# -*- coding: utf-8 -*-
''' HR 29/03/22 '''
import tkinter as tk
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import os

# ''' Workaround to avoid "module not found error" when importing psd_analyser '''
# import sys
# sys.path.append(os.path.join(os.path.dirname(__file__)))
import psdanalyser.psd_analyser as psda



class TkApp(tk.Tk):

    ''' Constructor '''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.args = args
        self.kwargs = kwargs

        ''' Fit and plot preferences
                - Default to CDF, PDF otherwise
                - Default to linear x scale, log scale otherwise '''
        self._plot_mode = True
        # self._fit_mode = True
        self._log_mode = False

        ''' All GUI set-up '''
        self.title("PSD-Analyser")

        self.frame_edge_col = 'black'
        self.frame_edge_thickness = 1

        self.grid_rowconfigure(0, weight = 1)
        self.grid_columnconfigure(0, weight = 1)

        self.main_frame = tk.Frame(self,
                           highlightbackground = self.frame_edge_col,
                           highlightthickness = self.frame_edge_thickness)
        self.main_frame.grid(column = 0, row = 0, sticky = "nsew")

        self.main_frame.grid_rowconfigure(0, weight = 3)
        self.main_frame.grid_rowconfigure(1, weight = 1)
        self.main_frame.grid_columnconfigure(0, weight = 1)
        self.main_frame.grid_columnconfigure(1, weight = 3)

        self.topleft = tk.Frame(self.main_frame,
                           highlightbackground = self.frame_edge_col,
                           highlightthickness = self.frame_edge_thickness)
        self.topright = tk.Frame(self.main_frame,
                            highlightbackground = self.frame_edge_col,
                            highlightthickness = self.frame_edge_thickness)
        self.bottomleft = tk.Frame(self.main_frame,
                              highlightbackground = self.frame_edge_col,
                              highlightthickness = self.frame_edge_thickness)
        self.bottomright = tk.Frame(self.main_frame,
                               highlightbackground = self.frame_edge_col,
                               highlightthickness = self.frame_edge_thickness)

        self.topleft.grid(row = 0, column = 0, sticky = 'nesw')
        self.topright.grid(row = 0, column = 1, sticky = 'nesw')
        self.bottomleft.grid(row = 1, column = 0, sticky = 'nesw')
        self.bottomright.grid(row = 1, column = 1, sticky = 'nesw')

        ''' File ops frame '''
        self.button_load_file = tk.Button(self.topleft,
                                          text = 'Load file',
                                          command = self.load_file)
        self.button_save_data = tk.Button(self.topleft,
                                          text = 'Save data',
                                          command = self.save_data)
        self.button_export_figure = tk.Button(self.topleft,
                                              text = 'Export figure',
                                              command = self.export_figure)

        self.button_load_file.grid(row = 0, sticky = 'nsew')
        self.button_save_data.grid(row = 1, sticky = 'nsew')
        self.button_export_figure.grid(row = 2, sticky = 'nsew')

        self.topleft.grid_columnconfigure(0, weight = 1)
        self.topleft.grid_rowconfigure(0, weight = 1)
        self.topleft.grid_rowconfigure(1, weight = 1)
        self.topleft.grid_rowconfigure(2, weight = 1)

        ''' Plot control frame '''
        self.button_previous = tk.Button(self.bottomright,
                                         text = '<--',
                                         command = self.previous_dataset)
        self.button_next = tk.Button(self.bottomright,
                                     text = '-->',
                                     command = self.next_dataset)
        self.button_fit_PDF = tk.Button(self.bottomright,
                                        text = 'Fit to PDF',
                                        command = self.fit_mode_PDF)
        self.button_fit_CDF = tk.Button(self.bottomright,
                                        text = 'Fit to CDF',
                                        command = self.fit_mode_CDF)
        self.button_plot_PDF = tk.Button(self.bottomright,
                                         text = 'Plot PDF',
                                         command = self.plot_mode_PDF)
        self.button_plot_CDF = tk.Button(self.bottomright,
                                         text = 'Plot CDF',
                                         command = self.plot_mode_CDF)
        self.button_plot_log = tk.Button(self.bottomright,
                                         text = 'Plot log scale',
                                         command = self.plot_mode_log)
        self.button_plot_linear = tk.Button(self.bottomright,
                                         text = 'Plot linear scale',
                                         command = self.plot_mode_linear)

        self.button_previous.grid(row = 0, column = 0, sticky = 'nsew')
        self.button_next.grid(row = 0, column = 1, sticky = 'nsew')
        self.button_fit_PDF.grid(row = 1, column = 0, sticky = 'nsew')
        self.button_fit_CDF.grid(row = 1, column = 1, sticky = 'nsew')
        self.button_plot_PDF.grid(row = 2, column = 0, sticky = 'nsew')
        self.button_plot_CDF.grid(row = 2, column = 1, sticky = 'nsew')
        self.button_plot_log.grid(row = 3, column = 0, sticky = 'nsew')
        self.button_plot_linear.grid(row = 3, column = 1, sticky = 'nsew')

        self.bottomright.grid_rowconfigure(0, weight = 1)
        self.bottomright.grid_rowconfigure(1, weight = 1)
        self.bottomright.grid_rowconfigure(2, weight = 1)
        self.bottomright.grid_rowconfigure(3, weight = 1)
        self.bottomright.grid_columnconfigure(0, weight = 1)
        self.bottomright.grid_columnconfigure(1, weight = 1)

        ''' Dataset info frame '''
        self.info_label = tk.Label(self.bottomleft,
                                   text = '',
                                   background = 'white')

        self.info_label.grid(sticky = 'nsew')
        self.bottomleft.grid_columnconfigure(0, weight = 1)
        self.bottomleft.grid_rowconfigure(0, weight = 1)

        ''' Plot frame '''
        self.fig = Figure()
        self.ax = self.fig.add_subplot()
        self.plot_canvas = FigureCanvasTkAgg(self.fig,
                                             master = self.topright)

        # plot_canvas.get_tk_widget().grid(sticky = "nsew")
        self.plot_canvas.get_tk_widget().pack(fill = 'both', expand = True)
        # self.topright.grid_columnconfigure(0, weight = 1)
        # self.topright.grid_rowconfigure(0, weight = 1)

        ''' Set up PSDAnalyser for data management '''
        self.manager = psda.PSDAnalyser(ax = self.ax)

        ''' Get plot label defaults from manager '''
        self.XLABEL_DEFAULT = self.manager.XLABEL_DEFAULT
        self.YLABEL_DEFAULT_CDF = self.manager.YLABEL_DEFAULT_CDF
        self.YLABEL_DEFAULT_PDF = self.manager.YLABEL_DEFAULT_PDF

        self._current = None

        ''' Initialise info frame '''
        self.update_info_frame()

    ''' All bindings '''
    def load_file(self):
        print('Loading file...')
        import_file = tk.filedialog.askopenfilename(defaultextension = ".xlsx",
                                                    filetypes = (("Excel files (.xlsx)", "*.xlsx"),
                                                                 ("Excel files (.xls)", "*.xls"),
                                                                 ("Comma-separated files (.csv)", "*.csv")))
        if import_file:
            print('Load file full path: ', import_file)
        else:
            print('No file loaded')
            return

        ''' Load file through manager '''
        self.manager.load_spreadsheet(import_file)

        ''' Grab first dataset '''
        self._current = min(self.manager.datasets)

        self.update_plot()
        self.update_info_frame()

    def save_data(self):
        print('Saving data...')

        ''' Get original file format of data '''
        original_file = self.manager._file
        if original_file:
            stub,ext = os.path.splitext(original_file)
        else:
            print('No data loaded')
            return

        ''' Default filename (full path) to add to dialog '''
        initialfile = stub + self.manager.SAVE_TEXT_DEFAULT

        save_file = tk.filedialog.asksaveasfilename(defaultextension = ext,
                                                    initialfile = initialfile,
                                                    filetypes = (("Excel files (.xlsx)", "*.xlsx"),
                                                                 ("Excel files (.xls)", "*.xls"),
                                                                 ("Comma-separated files (.csv)", "*.csv")))

        ''' Can pass None here, as defaults to loaded file name + "_computed" in that case '''
        self.manager.export_to_spreadsheet(file = save_file)

    def export_figure(self):
        print('Outputting figure...')
        initialfile = self.manager.IMAGE_DUMP_DEFAULT
        figure_file = tk.filedialog.asksaveasfilename(defaultextension = ".jpg",
                                                      initialfile = initialfile,
                                                      filetypes = (("JPEG (.jpg)", "*.jpg"),
                                                                   ("PNG (.png)", "*.png")))
        if figure_file:
            print('Figure export file full path: ', figure_file)
        else:
            print('No figure export file selected; returning')
            return

        self.manager.dump_figure(figure_file)

    def previous_dataset(self):
        if not self.manager._loaded:
            print('Cannot display previous dataset: no file loaded')
            return

        print('Moving to previous dataset...')
        ''' Check if first dataset; if so, loop to last '''
        if self._current == min(self.manager.datasets):
            self._current = max(self.manager.datasets)
        else:
            self._current -= 1
            while self._current not in self.manager.datasets:
                ''' Traverse datasets until valid one found '''
                self._current -= 1

        ''' Update other views '''
        self.update_plot()
        self.update_info_frame()

    def next_dataset(self):
        if not self.manager._loaded:
            print('Cannot display next dataset: no file loaded')
            return

        print('Moving to next dataset...')
        ''' Check if last dataset; if so, loop to first '''
        if self._current == max(self.manager.datasets):
            self._current = min(self.manager.datasets)
        else:
            self._current += 1
            while self._current not in self.manager.datasets:
                ''' Traverse datasets until valid one found '''
                self._current += 1

        ''' Update other views '''
        self.update_plot()
        self.update_info_frame()

    def fit_mode_CDF(self):
        if not self.manager._loaded:
            print('Cannot change fit mode: no file loaded')
            return

        print('Fitting to CDF and replotting...')
        # if self._fit_mode:
        if self.manager.datasets[self._current]['fit_by_CDF'] == True:
            print('Already in CDF fit mode; aborting...')
            return
        # self._fit_mode = True
        self.manager.toggle_fit_mode(self._current)

        ''' Update other views '''
        self.update_plot()
        self.update_info_frame()

    def fit_mode_PDF(self):
        if not self.manager._loaded:
            print('Cannot change fit mode: no file loaded')
            return

        print('Fitting to PDF and replotting...')
        # if self._fit_mode:
        if self.manager.datasets[self._current]['fit_by_CDF'] == False:
            print('Already in PDF fit mode; aborting...')
            return
        # self._fit_mode = True
        self.manager.toggle_fit_mode(self._current)

        ''' Update other views '''
        self.update_plot()
        self.update_info_frame()

    def plot_mode_CDF(self):
        if not self.manager._loaded:
            print('Cannot change plot mode: no file loaded')
            return

        print('Plotting CDF...')
        if self._plot_mode:
            print('Already in CDF plot mode; aborting...')
            return
        self._plot_mode = True

        ''' Update other views '''
        self.update_plot()
        self.update_info_frame()

    def plot_mode_PDF(self):
        if not self.manager._loaded:
            print('Cannot change plot mode: no file loaded')
            return

        print('Plotting PDF...')
        if not self._plot_mode:
            print('Already in PDF plot mode; aborting...')
            return
        self._plot_mode = False

        ''' Update other views '''
        self.update_plot()
        self.update_info_frame()

    def plot_mode_log(self):
        if not self.manager._loaded:
            print('Cannot change plot mode: no file loaded')
            return

        print('Plotting in log mode...')
        if self._log_mode:
            print('Already in log plot mode; aborting...')
            return
        self._log_mode = True

        ''' Update other views '''
        self.update_plot()
        self.update_info_frame()

    def plot_mode_linear(self):
        if not self.manager._loaded:
            print('Cannot change plot mode: no file loaded')
            return

        print('Plotting in log mode...')
        if not self._log_mode:
            print('Already in linear plot mode; aborting...')
            return
        self._log_mode = False

        ''' Update other views '''
        self.update_plot()
        self.update_info_frame()

    def update_plot(self):
        if not self.manager._loaded:
            return

        print('Updating plot...')
        if self.manager.datasets[self._current]['fit_by_CDF']:
            text1 = "CDF"
        else:
            text1 = "PDF"
        print('Fit mode: ', text1)

        if self._plot_mode:
            text2 = "CDF"
        else:
            text2 = "PDF"
        print('Plot mode: ', text2)

        ''' Clear plot axes '''
        try:
            print('Clearing plot axes...')
            self.ax.clear()
            print('Done')
        except Exception as e:
            print('Could not clear axes; exception follows')
            print(e)

        ''' Get common data to plot '''
        ds = self.manager.datasets[self._current]
        x_all = ds['bin_centres']
        C_all = ds['CDF']
        P_all = ds['PDF_normalised']
        x,C,P = psda.get_valids(x_all, C_all, P_all)
        # print('x_all:\n', x_all)
        # print('x:\n', x)
        fit_data = (ds['M'], ds['S'])

        ''' Get mode-specific data to plot '''
        if self._plot_mode:
            psda.plot_CDF(x, C, self.ax, fit_data = fit_data,
                          log_mode = self._log_mode,
                          xlabel = self.XLABEL_DEFAULT,
                          ylabel = self.YLABEL_DEFAULT_CDF)
        else:
            psda.plot_PDF(x, P, self.ax, fit_data = fit_data,
                          log_mode = self._log_mode,
                          xlabel = self.XLABEL_DEFAULT,
                          ylabel = self.YLABEL_DEFAULT_PDF)

        self.fig.canvas.draw()
        self.update_info_frame()

    def update_info_frame(self):
        if not self.manager._loaded:
            info_text = "No data loaded"
            self.info_label.configure(text = info_text)
            return

        print('Updating info frame...')
        if self.manager.datasets[self._current]['fit_by_CDF']:
            text1 = "CDF"
        else:
            text1 = "PDF"
        print('Fit mode: ', text1)

        if self._plot_mode:
            text2 = "CDF"
        else:
            text2 = "PDF"
        print('Plot mode: ', text2)

        ''' Get dataset index and display in info frame '''
        info_text = ("Dataset", str(self._current))
        self.info_label.configure(text = info_text)

        # ''' Dataset index text '''
        # text3 = self._current
        # ''' Dataset name text '''
        # text4 = ds['sample_name']
        # ''' Quantiles (measured) '''
        # text5 = (ds['d10'], ds['d50'], ds['d90'])
        # ''' Quantiles (lognormal fit) '''
        # text6 = (ds['d10_fit'], ds['d50_fit'], ds['d90_fit'])
        # ''' Means (measured) '''
        # text7 = (ds['D43'], ds['D32'])
        # ''' Means (lognormal fit) '''
        # text8 = (ds['D43_fit'], ds['D32_fit'])
        # ''' Lognormal parameters '''
        # text9 = (ds['M'], ds['S'])
        # ''' Dataset packing fraction '''
        # text10 = ds['phi']
        # info_text = str(("Plot mode: ", text1,
        #                 '\nFit mode: ', text2,
        #                 '\nDataset index: ', text3,
        #                 '\nDataset name: ', text4,
        #                 '\nQuantiles (measured; d10, d50, d90): ', text5,
        #                 '\nQuantiles (lognormal fit; d10, d50, d90): ', text6,
        #                 '\nMeans (measured; D43, D32): ', text7,
        #                 '\nMeans (lognormal; D43, D32); ', text8,
        #                 '\nLognormal params: ', text9,
        #                 '\nEstimated packing fraction: ', text10))
        # self.info_label.configure(text = info_text)



if __name__ == "__main__":
    app = TkApp()
    app.mainloop()