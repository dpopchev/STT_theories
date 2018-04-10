#!/usr/bin/env python


class plot_result:

    def __init__(self):

        self.path_res = ""
        self.path_daniela = ""

        self.data_res = []
        self.headline_res = []
        self.label_res = []

        self.data_daniela = []
        self.headline_daniela = []
        self.label_daniela = []

        self.units = self._units_coef_clac()

        self._set_default_paths(self)

        self.daniela_mapping = {
            "rho_c": 0,
            "AR": 1,
            "M": 2,
            "J": 3,
            "phiScal_c": 4,
            "p_c": 5
        }

        return

    @staticmethod
    def _units_coef_clac():
        """
        Calculates the unit coefficients

        Parameters
        ----------

        Returns
        -------
            : dictionary
                "density" in (double) g cm^-3
                "pressure" in (double) dyns cm^-3
                "r" in (double) km
                "j" in (double) m^2 kg
        """

        # mas of sun in kg
        const_msun = 1.9891e30
        # gravitational const in m^3kg^-1s^-2
        const_g = 6.67384e-11
        # speed of light in ms^-1
        const_c = 299792458

        units = {}

        # units of density in g cm^-3
        units["density"] = 1e-3 * const_c**6 / (const_g**3 * const_msun**2)

        # units of pressure in dyns cm^-3
        units["pressure"] = const_c**8 / (const_g**3 * const_msun**2) * 10

        # units of rad coordinate in km
        units["rad"] = 1e-3 * const_g * const_msun / const_c**2

        # units of moment of inertia
        units["j"] = 1e7 * const_g**2 * const_msun**3 / const_c**4

        return units

    def _set_default_paths(self):
        """
        sets default paths in self to the result directories

        Parameters
        ----------

        Returns
        -------
        """

        self.path_res = "/home/dimitar/projects/STT_theories/results"

        self.path_daniela = "/home/dimitar/Documents/Teaching_Materials/" \
            + "University/Uvod.Fizika.Cherni.Dupki/Doktorant/" \
            + "Daniela_Static_NS_Massive_SlowRot"

        print(
            "\n mine path_res: {} \n path_daniela: {} \n".format(
                self.path_res,
                self.path_daniela
            )
        )

        return

    def _clear_data(self):
        """
        clear the data in self:
            data_*, which contains the plot data
            headline_*, which is the naming of each column of each data
            label_*, which is name of the file the data come from

        Parameters
        ----------

        Returns
        -------
        """

        self.data_res.clear()
        self.headline_res.clear()
        self.label_res.clear()

        self.data_daniela.clear()
        self.headline_daniela.clear()
        self.label_daniela.clear()

        return

    @staticmethod
    def _check_if_StrNotBlank(string):
        """
        check if a sting is blank/empty

        Parameters
        ----------

        Returns
        -------
            True if string is not blank/empty
            False if string is blank/empty
        """

        return bool(string and string.strip())

    def _get_latest_file(self, path, model):
        """
        return the name of the latest modified file in path
        whos name is based on model

        Parameters
        ----------
            path: string
                the path to the directory of interest
            model: string
                how the file name looks like, it serves as wildcard
                of the type model + "*"

        Returns
        -------
            :string
                the name of latest modified file
        """

        import glob
        import os

        if (self._check_if_StrNotBlank(path) and
            self._check_if_StrNotBlank(model)):

            filename = max(
                glob.glob(path + "/" + model + "*"),
                key=os.path.getctime
            )

            filename = filename.split("/")[-1]

            print(
                "\n latest file in \n\t {} is: \n {} \n".format(
                    path, filename
                )
            )
        else:

            filename = ""

            print("\n path and name not set, some filename var!!! \n")

        return filename

    def _load_file(self, path, filename, data, headline, label):
        """
        load the data from <filename> in <path> by appending them to
        <data>, append the headline (the name of each column) in
        <headline> and the name of the file in <label>

        Parameters
        ----------
            path: string
                the path to the directory containing the file with data
            filename: string
                the name of the file which will be loaded
            data: list
                the loaded data will be appended here
            headline: list
                the name of each column will be appended here as a list
            label: list
                the name of the file will be apended here

        Returns
        ------
        """

        if not( self._check_if_StrNotBlank(_path) and
            self._check_if_StrNotBlank(filename) ):

            print(
                "\n from where to load a file!!! \n\t {} / {} \n".format(
                    path, filename
                )
            )

        with open(path + "/" + filename, "r") as f:
            all_data = f.readlines()

        headline.append(
            [
                i.strip() for i in
                all_data.pop(0).strip().split(" ")
                if
                "#" not in i and
                len(i.strip())
            ]
        )

        label.append(filename.split("/")[-1][-41:])

        data.append(
            [
                [] for i in all_data[0].strip().split(" ") if len(i.strip())
            ]
        )

        for cnt, line in enumerate(all_data):
            for apnd, num in zip(
                self.data_res[-1], [
                    float(i) for i in line.strip().split(" ") if len(i.strip())
                ]
            ):
                apnd.append(num)

        return

    @staticmethod
    def _single_grid_placement(fig, gs):
        """
            create single axes in <fig> using GridSpec <gs>
            and fig.add_subplot() method

            Parameters
            ---------
                fig: Figure
                    the figure whos geometry will be specified

                gs: class
                    predefined class of how many rows and columns the
                    figure should have

            Returns
            ------
                all_ax: list
                    list containing all axes defined by the <gs> in <fig>
        """

        all_ax = []
        all_ax.append(fig.add_subplot(gs[:]))

        return all_ax

    @staticmethod
    def _get_figure(nrows, ncols, grid_placement):
        """
            create Figure and assign specific subplot geometry defined by
            <nrows> amount of rows and <ncols> amount of columns
            defined using the function <grid_placement>

            Parameters
            ----------
                nrows: int
                    amount of rows the figure will have

                ncols: int
                    amount of columns the figure will have

                grid_placement: class
                    function to define the placement of the subplots in the
                    figure using <nrows> and <ncols> and return
                    list of axes

            Returns
            -------
                fig: Figure
                    the figure itself

                : list
                    the list of axes which are associated with the subplot
                    placement defined by <grid_placement> function using
                    GridSpec class
        """

        from matplotlib import pyplot as plt
        from matplotlib import style
        from matplotlib.gridspec import GridSpec

        style.use("seaborn-poster")
        gs = GridSpec(nrows=nrows, ncols=ncols)

        fig = plt.figure()
        fig.set_tight_layout(True)

        return fig, grid_placement(fig, gs)

    @staticmethod
    def _get_index(max_index):
        """
            get the user input choice of an index and return it
            by making sure it is an positive integer not greater than <max_index>

            Parameters
            ----------
                max_index: int
                    the maximum index available

            Return
            ------
                index: int
                    the index which the user chose
        """

        while "waiting for the number":
            try:
                index = int(input("\n\t\t index =... "))
            except ValueError:
                print("\n Not a number is you typed \n")
                continue

            if index < 0 or index >= max_index:
                print("\n {} not in allowed boundaries \n".format(index))
            else:
                break

        return index

    def _get_xy_current_file(self, headline):
        """
            list the enumerated <headlines> and let the user choose
            which will be plotted as Y and as X on the fig

            Parameters
            ----------
                headline: list
                    list of available headlines to choose from


            Return
            ------
                index_x: int
                    the index of the headline to be used on the abscissa

                index_y: int
                    the index of the headline to be used on the ordinate
        """

        print("\n choose columns to plot Y vs X \n")

        for cnt, value in enumerate(headline):
            print("\t {}: {}".format(cnt, value))

        print("\n\t Y... ")
        index_y = self._get_index(len(headline))

        print("\n\t X... ")
        index_x = self._get_index(len(headline))

        return index_x, index_y

    @staticmethod
    def _set_parms(ax, label_x, label_y):
        """
            set the parameters of the provided axes <ax>; it will modify to
            specific size the fonts and tick of the ax; also sets the
            format of the numbers on each axis and the labels using <label_x> and
            <label_y>

            Parameter
            ---------
                <ax>: class
                    the specific ax whos parameters will be set

                <label_x>: string
                    the abscissa label

                <label_y>: string
                    the ordinate label

            Return
            ------

        """

        from matplotlib.ticker import FormatStrFormatter

        fontsize = 12
        ticksize = 10

        ax.clear()

        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax.set_xlabel(label_x, fontsize=fontsize)
        ax.xaxis.set_tick_params(labelsize=ticksize)

        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax.set_ylabel(label_y, fontsize=fontsize)
        ax.yaxis.set_tick_params(labelsize=ticksize)

        return

    def load_latest_ResultF(self, filename_type="STT_phiScal_"):
        """
            load the latest modified file from self.path

            Parameters
            ----------
                filename_type: string
                    this string will be used as wildcard filename_type + "*"
                    to search in self.path_res; if no value is supplied, the
                    default "STT_phiScal_" will be used

            Return
            ------
        """

        fname_res = self._get_latest_file(
            self, self.path_res, filename_type
        )

        self._load_file(
            self,
            self.path_res,
            fname_res,
            self.data_res,
            self.headline_res,
            self.label_res
        )

        return

    def plot_last_ResultFile(self):
        """
            creates figure with single subplot on which it will plot the
            last appended data in <self.data_res> by letting choose which column
            to be placed on X and Y axis

            Parameters
            ---------

            Return
            ------
        """

        from matplotlib import pyplot as plt
        from IPython import get_ipython

        # if we start it in jupyter qtconsole
        # this will force open it in new window not inside
        get_ipython().run_line_magic('matplotlib', "qt5")

        index_x, index_y = self._get_xy_current_file(
            self, self.headline_res[-1]
        )

        fig, all_ax = self._get_figure(self, 1, 1, self._single_grid_placement)

        ax = all_ax[-1]

        self._set_parms(
            self, ax,
            self.headline_res[-1][index_x],
            self.headline_res[-1][index_y]
        )

        ax.plot(
            self.data_res[-1][index_x],
            self.data_res[-1][index_y],
            linewidth=1.5,
            label=self.label_res[-1]
        )

        ax.legend(loc="best", fontsize=8)
        plt.show()

        return

    def _get_filename_from_dir(self, path):
        """
            lists all the files in <path> and lets the user choose
            which file to be loaded

            Parameters
            ----------
            path: string
                the full path to the directory whos files will be listed

            Return
            ------
                : string
                    the filename of our choice
        """

        import os

        list_dir = os.listdir(path)

        for ind, val in enumerate(list_dir):
            print("{}: {}".format(ind, val))

        return list_dir[self._get_index(len(list_dir))]

    @staticmethod
    def _yes_or_no_answer():
        """
            get the answer of yes-no question

            Parameters
            ----------

            Return
            ------
        """

        while "Answer is invalid":
            ans = str(input("\n\t y/n answer... ")).lower().strip()

            if ans[:1] == "y":
                return True
            elif ans[:1] == "n":
                return False
            else:
                print("\n y or n answer only !\n")
                continue

    def load_many_mine_results(self):
        """
            clear the <self.data_res>, loads as many as I like data results from <self.path_res>

            Parameters
            ----------

            Return
            ------
        """



    def plot_compare_mine(self):
        """

        """

        from matplotlib import pyplot as plt
        from IPython import get_ipython

        get_ipython().run_line_magic('matplotlib', "qt5")

        self._clear_data()

        print("\n File to load...")
        while self.yes_or_no_answer():

            fname_res = self.get_filename_from_dir(self, self.path_res)
            print("\n now will load: \n\t {} \n".format(fname_res))

            self._load_file(
                self,
                self.path_res,
                fname_res,
                self.data_res,
                self.headline_res,
                self.label_res
            )

            print("\n Append another file? \n")


        index_x, index_y = self._get_xy_current_file(
            self, self.headline_res[-1]
        )

        fig, all_ax = self._get_figure(self, 1, 1, self._single_grid_placement)

        ax = all_ax[-1]

        self._set_parms(
            self, ax,
            self.headline_res[-1][index_x],
            self.headline_res[-1][index_y]
        )

        linestyle_color = self.get_prod_iter_ls_c()

        for d, l, lsc in zip(self.data_res, self.label_res, linestyle_color):

            ax.plot(
                d[index_x],
                d[index_y],
                linewidth=1.5,
                label=l,
                linestyle=lsc[0],
                color=lsc[1]
            )

        ax.legend(loc="best", fontsize=8)
        plt.show()

    @staticmethod
    def load_file_daniela(self, _path, _filename, _data, _label):

        if not(self._check_if_StrNotBlank(_path) and
        self._check_if_StrNotBlank(_filename)):

            print(
                "from where to load a daniela file!!! {} / {}".format(
                    _path, _filename
                )
            )

        with open(_path + "/" + _filename, "r") as f:
            all_data = f.readlines()

        # get rid of inline text
        for i, line in enumerate(all_data):
            if "rho" in line:
                del all_data[i]

        # i am not interested in all data, so only the indexes here
        # will be saved
        # since in some files the central pressure is in different column
        # we check it and take into account
        if _filename in [
            "models_APR_beta-4.5_mphi5e-3.dat",
            "models_APR_beta-6.0_mphi0.dat",
            "models_APR_beta-6.0_mphi1e-3.dat"
        ]:
            indexes = [ 0, 1, 2, 3, 4, 10 ]
        else:
            indexes = [ 0, 1, 2, 3, 4, 11 ]

        _units = [
            self.units["density"], self.units["rad"],
            1, self.units["j"], 1, 1
        ]

        _data.append( [ [] for i in indexes ] )

        for cnt, line in enumerate(all_data):
            for i, u, d in zip(indexes, _units, _data[-1]):
                try:
                    d.append(float(line.split(" ")[i]) / u)
                except ValueError:
                    print(
                        "\n ValueError: line {}: {} \n".format(
                            cnt, line
                        )
                    )
                    break

        _data[-1][-2] = \
            [ (-1)*i for i in _data[-1][-2] ]

        _label.append(_filename)

        return

    @staticmethod
    def get_max_y_coord(x, y):

        max_y = max(y, key=abs)
        max_x = x[y.index(max_y)]

        return max_x, max_y

    @staticmethod
    def get_center_coordinates_axes(ax):

        center_x = (max(ax.get_xlim(), key=abs) + min(ax.get_xlim(), key=abs))/2
        center_y = (max(ax.get_ylim(), key=abs) + min(ax.get_ylim(), key=abs))/2

        return center_x, center_y

    @staticmethod
    def get_axes_lim(ax):

        return max(ax.get_xlim(), key=abs), max(ax.get_ylim(), key=abs), \
        min(ax.get_xlim(), key=abs), min(ax.get_ylim(), key=abs)

    def plot_compare_OneSet_mine_daniela(self):

        from matplotlib import pyplot as plt
        from IPython import get_ipython

        mine_color = "y"
        daniela_color = "m"

        get_ipython().run_line_magic('matplotlib', "qt5")

        self._clear_data()

        fname_res = self.get_filename_from_dir(self, self.path_res)
        print("\n now to load... {} \n".format(fname_res))

        self._load_file(
            self,
            self.path_res,
            fname_res,
            self.data_res,
            self.headline_res,
            self.label_res
        )

        fname_res_daniela = self.get_filename_from_dir(self, self.path_daniela)
        print("\n now to load... {} \n".format(fname_res))

        self.load_file_daniela(
            self,
            self.path_daniela,
            fname_res_daniela,
            self.data_daniela,
            self.label_daniela
        )

        index_x, index_y = self._get_xy_current_file(
            self, self.headline_res[-1]
        )

        fig, all_ax = self._get_figure(self, 1, 1, self._single_grid_placement)

        ax = all_ax[-1]

        self._set_parms(
            self, ax,
            self.headline_res[-1][index_x],
            self.headline_res[-1][index_y]
        )

        ax.plot(
            self.data_res[-1][index_x],
            self.data_res[-1][index_y],
            linewidth=1.5,
            label=self.label_res[-1],
            color=mine_color
        )

        ax.plot(
            self.data_daniela[-1][
                self.daniela_mapping[self.headline_res[-1][index_x]]
            ],
            self.data_daniela[-1][
                self.daniela_mapping[self.headline_res[-1][index_y]]
            ],
            marker = "o",
            markersize = 5,
            alpha = 0.4,
            linewidth=0,
            color=daniela_color,
            label=fname_res_daniela
        )

        ax.legend(loc="best", fontsize=8)

        plt.show()

    def plot_compare_OneSet_mine_daniela_with_annotations(self):

        from matplotlib import pyplot as plt
        from IPython import get_ipython

        mine_color = "y"
        daniela_color = "m"

        get_ipython().run_line_magic('matplotlib', "qt5")

        self._clear_data()

        fname_res = self.get_filename_from_dir(self, self.path_res)
        print("\n now to load... {} \n".format(fname_res))

        self._load_file(
            self,
            self.path_res,
            fname_res,
            self.data_res,
            self.headline_res,
            self.label_res
        )

        fname_res_daniela = self.get_filename_from_dir(self, self.path_daniela)
        print("\n now to load... {} \n".format(fname_res))

        self.load_file_daniela(
            self,
            self.path_daniela,
            fname_res_daniela,
            self.data_daniela,
            self.label_daniela
        )

        index_x, index_y = self._get_xy_current_file(
            self, self.headline_res[-1]
        )

        fig, all_ax = self._get_figure(self, 1, 1, self._single_grid_placement)

        ax = all_ax[-1]

        self._set_parms(
            self, ax,
            self.headline_res[-1][index_x],
            self.headline_res[-1][index_y]
        )

        ax.plot(
            self.data_res[-1][index_x],
            self.data_res[-1][index_y],
            linewidth=1.5,
            label=self.label_res[-1],
            color=mine_color
        )

        max_val_mine_x, max_val_mine_y = self.get_max_y_coord(
            self.data_res[-1][index_x],
            self.data_res[-1][index_y]
        )

        ax.plot(
            self.data_daniela[-1][
                self.daniela_mapping[self.headline_res[-1][index_x]]
            ],
            self.data_daniela[-1][
                self.daniela_mapping[self.headline_res[-1][index_y]]
            ],
            marker = "o",
            markersize = 5,
            alpha = 0.4,
            linewidth=0,
            color=daniela_color,
            label=fname_res_daniela
        )

        max_val_danieal_x, max_val_danieal_y = self.get_max_y_coord(
            self.data_daniela[-1][
                self.daniela_mapping[self.headline_res[-1][index_x]]
            ],
            self.data_daniela[-1][
                self.daniela_mapping[self.headline_res[-1][index_y]]
            ]
        )

        ax.legend(loc="best", fontsize=8)

        xmax, ymax, xmin, ymin = self.get_axes_lim(ax)
        xcenter, ycenter = self.get_center_coordinates_axes(ax)

        ax.annotate(
            "max_daniela --> ({:.3e}, {:.3e})".format(
                max_val_danieal_x, max_val_danieal_y
            ),
            xy=(max_val_danieal_x,max_val_danieal_y),
            xytext=(
                (xcenter + xmin)/2,
                ycenter
            ),
            arrowprops=dict(
                arrowstyle="-|>",connectionstyle="arc", lw=1, color=daniela_color
            ),
            size=10, ha="center"
        )

        ax.annotate(
            "max_mine --> ({:.3e}, {:.3e})".format(
                max_val_mine_x, max_val_mine_y
            ),
            xy=(max_val_mine_x,max_val_mine_y),
            xytext=(
                (xcenter + xmax)/2,
                ycenter
            ),
            arrowprops=dict(
                arrowstyle="-|>",connectionstyle="arc", lw=1, color=mine_color
            ),
            size=10, ha="center"
        )

        plt.show()

    @staticmethod
    def get_prod_iter_ls_c():

        from itertools import product
        from itertools import cycle
        from random import shuffle

        styles_line = [ "-", "--", "-." ]
        styles_color = [ "b", "g", "r", "c", "m", "y", "k" ]

        shuffle(styles_line)
        shuffle(styles_color)

        return cycle(product( styles_line, styles_color ))

    def plot_compare_mine_one_daniela(self):

        from matplotlib import pyplot as plt
        from IPython import get_ipython

        get_ipython().run_line_magic('matplotlib', "qt5")

        self._clear_data()

        print("\n File to load...")
        while self.yes_or_no_answer():

            fname_res = self.get_filename_from_dir(self, self.path_res)
            print("\n now will load: \n\t {} \n".format(fname_res))

            self._load_file(
                self,
                self.path_res,
                fname_res,
                self.data_res,
                self.headline_res,
                self.label_res
            )

            print("\n Append another file? \n")

        print("\n Daniela file to load \n")
        fname_res_daniela = self.get_filename_from_dir(self, self.path_daniela)
        print("\n now to load... {} \n".format(fname_res))

        self.load_file_daniela(
            self,
            self.path_daniela,
            fname_res_daniela,
            self.data_daniela,
            self.label_daniela
        )

        index_x, index_y = self._get_xy_current_file(
            self, self.headline_res[-1]
        )

        fig, all_ax = self._get_figure(self, 1, 1, self._single_grid_placement)

        ax = all_ax[-1]

        self._set_parms(
            self, ax,
            self.headline_res[-1][index_x],
            self.headline_res[-1][index_y]
        )

        linestyle_color = self.get_prod_iter_ls_c()

        for d, l, lsc in zip(self.data_res, self.label_res, linestyle_color):

            ax.plot(
                d[index_x],
                d[index_y],
                linewidth=1.5,
                label=l,
                linestyle=lsc[0],
                color=lsc[1]
            )

        ls_c = next(linestyle_color)
        ax.plot(
            self.data_daniela[-1][
                self.daniela_mapping[self.headline_res[-1][index_x]]
            ],
            self.data_daniela[-1][
                self.daniela_mapping[self.headline_res[-1][index_y]]
            ],
            marker = "o",
            markersize = 5,
            alpha = 0.4,
            linewidth=0,
            linestyle=ls_c[0],
            color=ls_c[1],
            label=fname_res_daniela
        )

        ax.legend(loc="best", fontsize=8)
        plt.show()

if __name__ == "__main__":

    import sys

    print("\n Hello from {} \n", sys.argv[0])
