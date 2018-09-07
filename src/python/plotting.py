#!/usr/bin/env python

import os
import glob
import shutil
import itertools
import pathlib
import random

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from IPython import get_ipython
from matplotlib import style
from matplotlib.gridspec import GridSpec
from numpy.polynomial.polynomial import polyfit

get_ipython().run_line_magic("matplotlib", "qt")

class plot_result:

    def __init__(self):

        self.my_ResPath = None
        self.my_EOSname = None
        self.my_fname_starts = None
        self_my_file = None
        self.my_headline = None
        self.my_data = None
        self.my_label = None

        #~ self.amount_of_points = 10

        self.kalin_path = None
        self.kalin_file = None
        self.kalin_headline = None
        self.kalin_data = None
        self.kalin_label = None
        self.kalin_mapping = {
            "rho_c": 0,
            "AR": 1,
            "M": 2,
            "J": 3,
            "phiScal_c": 4,
            "p_c": 5
        }

        self.units = self._units_coef_clac()

        self.specific_ms = None
        self.specific_ls = None
        self.specific_c = None

        #~ self.data_daniela = []
        #~ self.headline_daniela = []
        #~ self.label_daniela = []

        #~ self.daniela_mapping = {
            #~ "rho_c": 0,
            #~ "AR": 1,
            #~ "M": 2,
            #~ "J": 3,
            #~ "phiScal_c": 4,
            #~ "p_c": 5
        #~ }

        #~ self.data_kalin = []
        #~ self.headline_kalin = []
        #~ self.label_kalin = []

        #~ self._set_default_paths()

        return

    def set_my_ResPath(
        self,
        my_ResPath = "~/projects/STT_theories/results"
    ):
        """
        path to the results of the shootings
        """

        self.my_ResPath = os.path.expanduser(my_ResPath)

        return

    def set_severalEOSs_ms_ls_c(self, severalEOSs):
        """
        for consistent marking on all graphs for the current instance

        EXAMPLE INPUT
        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        self.specific_ms = self._get_specific_ms( list(
            set( [ _["name"] for _ in severalEOSs ] )
        ) )

        self.specific_ls = self._get_specific_ls( list(
            set( [ _["m"] for _ in severalEOSs ] )
        ) )

        self.specific_c = lambda_linestyle = self._get_specific_c( list(
            set( [ _["lambda"] for _ in severalEOSs ] )
        ) )

        return


    def get_my_ResPath(self):

        return self.my_ResPath

    def set_my_EOSname(self, EOS_name):
        """
        set the name of EOS to be plotted, as it is a directory in my_ResPath
        """
        self.my_EOSname = EOS_name

        return

    def get_my_EOSname(self):

        return self.my_EOSname

    def get_my_latests_res(self, fname = "STT_phiScal_J_"):
        """
        return the latest result file in my_ResPath
        """

        try:
            return max(
                glob.glob(os.path.join(self.my_ResPath, fname + "*")),
                key = os.path.getctime
            )

        except ValueError:
            return None

    def move_my_latest_res(self):
        """
        move the latest result file from my_ResPath to the EOS_model name dir
        """

        full_latest_res = self.get_my_latests_res()
        if full_latest_res:
            latest_result = os.path.basename(full_latest_res)
        else:
            print("\n No latest result at \n\t {} \n".format(
                    self.my_ResPath
                )
            )
            return

        src = os.path.join( self.my_ResPath, latest_result )
        dst = os.path.join( self.my_ResPath, self.my_EOSname, latest_result )

        print(
            "\n moving \n\t from {} \n\t to {} \n".format(
                src,
                dst
            )
        )

        shutil.move(src, dst)

        return

    def get_my_latests_res_file(self, fname = "STT_phiScal_J_"):
        """
        return the latest result file in my_ResPath for current EOS
        """

        try:
            return max(
                glob.glob(os.path.join(
                        self.my_ResPath, fname + "*"
                ) ),
                key = os.path.getctime
            )

        except ValueError:

            return None

    def plot_latest_resEOSname_severalEOSs(self, severalEOSs):
        """
        get the latest result from the root folder and plot it alongside
        severalEOSs, which again is list of dicsts in the expected form
        """

        label, headline, data = self.get_resEOSname_data(
            self.get_my_latests_res_file()
        )

        fig, all_axes = self._get_figure(2,2, self._4by4_grid_placement)

        ax_M_AR = all_axes[0]
        ax_J_M = all_axes[1]
        ax_phiScal_c_p_c = all_axes[2]
        ax_rho_c_p_c = all_axes[3]

        self._set_parms( ax_M_AR, "AR", "M" )
        self._set_parms( ax_J_M, "M", "J" )
        self._set_parms( ax_phiScal_c_p_c, "$p_c$", "$\\varphi_c$" )
        self._set_parms( ax_rho_c_p_c, "$p_c$", "$\\rho_c$" )

        ax_M_AR.plot(
            data[3],
            data[2],
            linewidth=0,
            linestyle="",
            marker = "o",
            markersize = 3
        )

        ax_J_M.plot(
            data[2],
            data[5],
            linewidth=0,
            linestyle="",
            marker = "o",
            markersize = 3
        )

        ax_phiScal_c_p_c.plot(
            data[0],
            data[1],
            linewidth=0,
            linestyle="",
            marker = "o",
            markersize = 3
        )

        ax_rho_c_p_c.plot(
            data[0],
            data[4],
            linewidth=0,
            linestyle="",
            marker = "o",
            markersize = 3
        )

        val_EOSname, val_beta, val_m, val_lambda = self._get_parameter_values(label)

        plt.suptitle(
            "EOS = {}; beta = {:.1f}; m = {:.1e}; lambda = {:.1e}".format(
                val_EOSname, val_beta, val_m, val_lambda
            ),
            fontsize=10, y=0.998
        )

        all_label, all_headline, all_data = self.get_severalEOS_data(severalEOSs)

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            ls, lc, ms, mc = self._get_ls_lc_ms_mc()

            ax_M_AR.plot(
                data[3],
                data[2],
                label = "{}"
                    "\n\t $\\beta$ = {:.1f}"
                    "\n\t m = {:.1e}"
                    "\n\t $\\lambda$ = {:.1e}".format(
                    eos["name"], eos["beta"], eos["m"], eos["lambda"]
                ),
                color = lc,
                linestyle = ls,
                marker = ms,
                markerfacecolor = mc,
                markeredgecolor = mc,
                markersize = 2.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[3]),
                alpha = 0.5
            )

            ax_J_M.plot(
                data[2],
                data[5],
                label = "{}"
                    "\n\t $\\beta$ = {:.1f}"
                    "\n\t m = {:.1e}"
                    "\n\t $\\lambda$ = {:.1e}".format(
                    eos["name"], eos["beta"], eos["m"], eos["lambda"]
                ),
                color = lc,
                linestyle = ls,
                marker = ms,
                markerfacecolor = mc,
                markeredgecolor = mc,
                markersize = 2.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[3]),
                alpha = 0.5
            )

            ax_phiScal_c_p_c.plot(
                data[0],
                data[1],
                label = "{}"
                    "\n\t $\\beta$ = {:.1f}"
                    "\n\t m = {:.1e}"
                    "\n\t $\\lambda$ = {:.1e}".format(
                    eos["name"], eos["beta"], eos["m"], eos["lambda"]
                ),
                color = lc,
                linestyle = ls,
                marker = ms,
                markerfacecolor = mc,
                markeredgecolor = mc,
                markersize = 2.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[3]),
                alpha = 0.5
            )

            ax_rho_c_p_c.plot(
                data[0],
                data[4],
                label = "{}"
                    "\n\t $\\beta$ = {:.1f}"
                    "\n\t m = {:.1e}"
                    "\n\t $\\lambda$ = {:.1e}".format(
                    eos["name"], eos["beta"], eos["m"], eos["lambda"]
                ),
                color = lc,
                linestyle = ls,
                marker = ms,
                markerfacecolor = mc,
                markeredgecolor = mc,
                markersize = 2.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[3]),
                alpha = 0.5
            )


        ax_rho_c_p_c.legend(
            loc="best",
            fontsize=8,
            handlelength=3.2,
            numpoints=1,
            fancybox=True,
            markerscale = 1.5
        )

        plt.show()

        return

    def get_resEOSname_data(self, fpath):
        """
        get the following data form fpath, which is the full path to file
            label - the name of the file as string
            headline - the name for each column, as a list
            data - the data itself as a list, each sublist is different column
        """

        with open(fpath, "r") as f:
            all_data = f.readlines()

        label = os.path.basename(fpath)

        headline = [
            _.strip() for _ in all_data.pop(0).strip().split(" ")
            if
            "#" not in _ and
            len(_.strip())
        ]

        data = [
            [] for _ in all_data[0].strip().split(" ") if len(_.strip())
        ]

        for line in all_data:

            for d, n in zip(
                data, [ float(_) for _ in line.strip().split(" ") if len(_.strip()) ]
            ):
                d.append(n)

            if abs(data[1][-1]) > 1e-5 and data[1][-1] > 0:
                data[1][-1] *= (-1)

        return label, headline, data

    def get_severalEOS_data(
        self,
        severalEOSs,
        fname = "STT_phiScal_J"
    ):
        """
        for provided list of dictionaries called severalEOSs get the data
        the dict has following structure
        {
            "name": EOSname_string,
            "beta": Value_Beta,
            "m": Value_M,
            "lambda": Value_lambda
        }
        """

        all_label = []
        all_headline = []
        all_data = []

        for eos in severalEOSs:

            EOSname = "_".join( [
                fname,
                eos["name"],
                "beta{:.3e}".format(eos["beta"]),
                "m{:.3e}".format(eos["m"]),
                "lambda{:.3e}".format(eos["lambda"])
            ] )

            EOSpath = os.path.join( self.my_ResPath, eos["name"], EOSname)

            _label, _headline, _data = self.get_resEOSname_data(EOSpath)

            all_label.append(_label)
            all_headline.append(_headline)
            all_data.append(_data)

        return all_label, all_headline, all_data

    def get_uniEOSname_data_uniI(self, fpath):
        """
        for provided filepaths to tilde I return the entries as nested list
        labels and headline
        """

        with open(fpath, "r") as f:
            all_data = f.readlines()

        label = os.path.basename(fpath)

        headline = [
            _.strip() for _ in all_data.pop(0).strip().split(" ")
            if
            "#" not in _ and
            len(_.strip())
        ]

        data = [
            [] for _ in all_data[0].strip().split(" ") if len(_.strip())
        ]

        for line in all_data:

            for d, n in zip(
                data, [ float(_) for _ in line.strip().split(" ") if len(_.strip()) ]
            ):
                d.append(n)

        return label, headline, data

    def get_severalEOS_uniTildeI_data(
        self,
        severalEOSs,
        fname = "STT_phiScal_J"
    ):
        """
        for provided list of dictionaries called severalEOSs get the data for
        universal I, which are in the "Fitting" directory of each EOS

        one of the is *_tildeI for tilde I
        other is *_barI for bar I

        for each entyr in severalEOS return two nested lists barI and tildeI

        {
            "name": EOSname_string,
            "beta": Value_Beta,
            "m": Value_M,
            "lambda": Value_lambda
        }
        """

        all_label = []
        all_headline = []
        all_data = []

        for eos in severalEOSs:

            EOSname_tildeI = "_".join( [
                fname,
                eos["name"],
                "beta{:.3e}".format(eos["beta"]),
                "m{:.3e}".format(eos["m"]),
                "lambda{:.3e}".format(eos["lambda"]),
                "tildeI"
            ] )

            EOSpath_tildeI = os.path.join(
                self.my_ResPath, eos["name"], "Fitting", EOSname_tildeI
            )

            _label, _headline, _data = self.get_uniEOSname_data_uniI(EOSpath_tildeI)

            all_label.append(_label)
            all_headline.append(_headline)
            all_data.append(_data)

        return all_label, all_headline, all_data

    def get_severalEOS_uniBarI_data(
        self,
        severalEOSs,
        fname = "STT_phiScal_J"
    ):
        """
        for provided list of dictionaries called severalEOSs get the data for
        universal I, which are in the "Fitting" directory of each EOS

        one of the is *_tildeI for tilde I
        other is *_barI for bar I

        for each entyr in severalEOS return two nested lists barI and tildeI

        {
            "name": EOSname_string,
            "beta": Value_Beta,
            "m": Value_M,
            "lambda": Value_lambda
        }
        """

        all_label = []
        all_headline = []
        all_data = []

        for eos in severalEOSs:

            EOSname_barI = "_".join( [
                fname,
                eos["name"],
                "beta{:.3e}".format(eos["beta"]),
                "m{:.3e}".format(eos["m"]),
                "lambda{:.3e}".format(eos["lambda"]),
                "barI"
            ] )

            EOSpath_tildeI = os.path.join(
                self.my_ResPath, eos["name"], "Fitting", EOSname_barI
            )

            _label, _headline, _data = self.get_uniEOSname_data_uniI(EOSpath_tildeI)

            all_label.append(_label)
            all_headline.append(_headline)
            all_data.append(_data)

        return all_label, all_headline, all_data

    def plot_severalEOSs_MvsR(self, severalEOSs):
        """
        plot several EOSs by listing them in <severalEOSs> with dictionaries
        see get_severalEOS_data for the format

        EXAMPLE INPUT
        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        all_label, all_headline, all_data = self.get_severalEOS_data(severalEOSs)

        fig, all_axes = self._get_figure(1,1,self._1by1_grid_placement)

        ax = all_axes[0]

        self._set_parms(ax, "R [km]", "$M/M_{\odot}$")

        EOS_markers = self._get_specific_ms( list(
            set( [ _["name"] for _ in severalEOSs ] )
        ) )

        m_colors = self._get_specific_c( list(
            set( [ _["m"] for _ in severalEOSs ] )
        ) )

        lambda_linestyle = self._get_specific_ls( list(
            set( [ _["lambda"] for _ in severalEOSs ] )
        ) )

        for label, data, eos in zip( all_label, all_data, severalEOSs ):
            data[3] = [ _*self.units["R"] for _ in data[3] ]

            #~ old way it set random at every iiteration for every entry
            #~ ls, lc, ms, mc = self._get_ls_lc_ms_mc()

            #~ label entry for each line, but now it will specific
            #~ label_entry = "{}"
                #~ "\n\t $\\beta$ = {:.1f}"
                #~ "\n\t m = {:.1e}"
                #~ "\n\t $\\lambda$ = {:.1e}".format(
                #~ eos["name"], eos["beta"], eos["m"], eos["lambda"]
            #~ )

            ax.plot(
                data[3],
                data[2],
                label = None,
                color = m_colors.get(eos["m"], None),
                linestyle = lambda_linestyle.get(eos["lambda"], None),
                marker = EOS_markers.get(eos["name"], None),
                markerfacecolor = m_colors.get(eos["m"], None),
                markeredgecolor = m_colors.get(eos["m"], None),
                markersize = 5.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[2])
            )

        legend_eos = [
            Line2D(
                [0], [0],
                color="b",
                marker = __,
                linewidth = 0,
                label = _
            )
            for _, __ in EOS_markers.items()
        ]

        legend_lambda = [
            Line2D(
                [0], [0],
                color="b",
                marker = None,
                linewidth = 1.5,
                linestyle = __,
                label = "$\\lambda =$ {:.2e}".format(_)
            )
            for _, __ in lambda_linestyle.items()
        ]

        legend_m = [
            Line2D(
                [0], [0],
                color=__,
                marker = None,
                linewidth = 1.5,
                linestyle = "-",
                label = "$m =$ {:.2e}".format(_)
            )
            for _, __ in m_colors.items()
        ]

        ax.legend(
            handles = [*legend_eos, *legend_lambda, *legend_m],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol=3,
            frameon = False,
        )

        plt.show()

        return

    def plot_severalEOSs_phiScal_cVSp_c(self, severalEOSs):
        """
        plot several EOSs by listing them in <severalEOSs> with dictionaries
        see get_severalEOS_data for the format

        EXAMPLE INPUT
        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        all_label, all_headline, all_data = self.get_severalEOS_data(severalEOSs)

        fig, all_axes = self._get_figure(1,1,self._1by1_grid_placement)

        ax = all_axes[0]

        self._set_parms(ax, "$p_c$", "$\\varphi_c$")

        EOS_markers = self._get_specific_ms( list(
            set( [ _["name"] for _ in severalEOSs ] )
        ) )

        m_colors = self._get_specific_c( list(
            set( [ _["m"] for _ in severalEOSs ] )
        ) )

        lambda_linestyle = self._get_specific_ls( list(
            set( [ _["lambda"] for _ in severalEOSs ] )
        ) )

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ ls, lc, ms, mc = self._get_ls_lc_ms_mc()

            ax.plot(
                data[0],
                data[1],
                label = None,
                color = m_colors.get(eos["m"], None),
                linestyle = lambda_linestyle.get(eos["lambda"], None),
                marker = EOS_markers.get(eos["name"], None),
                markerfacecolor = m_colors.get(eos["m"], None),
                markeredgecolor = m_colors.get(eos["m"], None),
                markersize = 5.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[1])
            )

        legend_eos = [
            Line2D(
                [0], [0],
                color="b",
                marker = __,
                linewidth = 0,
                label = _
            )
            for _, __ in EOS_markers.items()
        ]

        legend_lambda = [
            Line2D(
                [0], [0],
                color="b",
                marker = None,
                linewidth = 1.5,
                linestyle = __,
                label = "$\\lambda =$ {:.2e}".format(_)
            )
            for _, __ in lambda_linestyle.items()
        ]

        legend_m = [
            Line2D(
                [0], [0],
                color=__,
                marker = None,
                linewidth = 1.5,
                linestyle = "-",
                label = "$m =$ {:.2e}".format(_)
            )
            for _, __ in m_colors.items()
        ]

        ax.legend(
            handles = [*legend_eos, *legend_lambda, *legend_m],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol=3,
            frameon = False,
        )

        plt.show()

        return

    def plot_severalEOSs_phiScal_cvsrho_c(self, severalEOSs):
        """
        plot several EOSs by listing them in <severalEOSs> with dictionaries
        see get_severalEOS_data for the format

        EXAMPLE INPUT
        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        all_label, all_headline, all_data = self.get_severalEOS_data(severalEOSs)

        fig, all_axes = self._get_figure(1,1,self._1by1_grid_placement)

        ax = all_axes[0]

        self._set_parms(ax, "$\\rho_c [g/cm^3]$", "$\\varphi_c$")

        EOS_markers = self._get_specific_ms( list(
            set( [ _["name"] for _ in severalEOSs ] )
        ) )

        m_colors = self._get_specific_c( list(
            set( [ _["m"] for _ in severalEOSs ] )
        ) )

        lambda_linestyle = self._get_specific_ls( list(
            set( [ _["lambda"] for _ in severalEOSs ] )
        ) )

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            data[-2] = [ _*self.units["density"] for _ in data[-2] ]

            #~ ls, lc, ms, mc = self._get_ls_lc_ms_mc()
            #~ label_entry = "{}"
                #~ "\n\t $\\beta$ = {:.1f}"
                #~ "\n\t m = {:.1e}"
                #~ "\n\t $\\lambda$ = {:.1e}".format(
                #~ eos["name"], eos["beta"], eos["m"], eos["lambda"]
            #~ )

            ax.plot(
                data[-2],
                data[1],
                label = None,
                color = m_colors.get(eos["m"], None),
                linestyle = lambda_linestyle.get(eos["lambda"], None),
                marker = EOS_markers.get(eos["name"], None),
                markerfacecolor = m_colors.get(eos["m"], None),
                markeredgecolor = m_colors.get(eos["m"], None),
                markersize = 5.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[1])
            )

        legend_eos = [
            Line2D(
                [0], [0],
                color="b",
                marker = __,
                linewidth = 0,
                label = _
            )
            for _, __ in EOS_markers.items()
        ]

        legend_lambda = [
            Line2D(
                [0], [0],
                color="b",
                marker = None,
                linewidth = 1.5,
                linestyle = __,
                label = "$\\lambda =$ {:.2e}".format(_)
            )
            for _, __ in lambda_linestyle.items()
        ]

        legend_m = [
            Line2D(
                [0], [0],
                color=__,
                marker = None,
                linewidth = 1.5,
                linestyle = "-",
                label = "$m =$ {:.2e}".format(_)
            )
            for _, __ in m_colors.items()
        ]

        ax.legend(
            handles = [*legend_eos, *legend_lambda, *legend_m],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol=3,
            frameon = False,
        )

        plt.show()

        return

    def plot_severalEOSs_Mvsrho_c(self, severalEOSs):
        """
        plot several EOSs by listing them in <severalEOSs> with dictionaries
        see get_severalEOS_data for the format

        EXAMPLE INPUT
        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        all_label, all_headline, all_data = self.get_severalEOS_data(severalEOSs)

        fig, all_axes = self._get_figure(1,1,self._1by1_grid_placement)

        ax = all_axes[0]

        self._set_parms(ax, "$\\rho_c [g/cm^3]$", "$M/M_{\odot}$")

        EOS_markers = self._get_specific_ms( list(
            set( [ _["name"] for _ in severalEOSs ] )
        ) )

        m_colors = self._get_specific_c( list(
            set( [ _["m"] for _ in severalEOSs ] )
        ) )

        lambda_linestyle = self._get_specific_ls( list(
            set( [ _["lambda"] for _ in severalEOSs ] )
        ) )

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            data[-2] = [ _*self.units["density"] for _ in data[-2] ]

            #~ old way it set random at every iiteration for every entry
            #~ ls, lc, ms, mc = self._get_ls_lc_ms_mc()
            #~ label entry for each line, but now it will specific
            #~ label_entry = "{}"
                #~ "\n\t $\\beta$ = {:.1f}"
                #~ "\n\t m = {:.1e}"
                #~ "\n\t $\\lambda$ = {:.1e}".format(
                #~ eos["name"], eos["beta"], eos["m"], eos["lambda"]
            #~ )

            ax.plot(
                data[-2],
                data[2],
                label = None,
                color = m_colors.get(eos["m"], None),
                linestyle = lambda_linestyle.get(eos["lambda"], None),
                marker = EOS_markers.get(eos["name"], None),
                markerfacecolor = m_colors.get(eos["m"], None),
                markeredgecolor = m_colors.get(eos["m"], None),
                markersize = 5.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[2])
            )

        legend_eos = [
            Line2D(
                [0], [0],
                color="b",
                marker = __,
                linewidth = 0,
                label = _
            )
            for _, __ in EOS_markers.items()
        ]

        legend_lambda = [
            Line2D(
                [0], [0],
                color="b",
                marker = None,
                linewidth = 1.5,
                linestyle = __,
                label = "$\\lambda =$ {:.2e}".format(_)
            )
            for _, __ in lambda_linestyle.items()
        ]

        legend_m = [
            Line2D(
                [0], [0],
                color=__,
                marker = None,
                linewidth = 1.5,
                linestyle = "-",
                label = "$m =$ {:.2e}".format(_)
            )
            for _, __ in m_colors.items()
        ]

        ax.legend(
            handles = [*legend_eos, *legend_lambda, *legend_m],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol=3,
            frameon = False,
        )
        plt.show()

        return

    def plot_severalEOSs_JvsM(self, severalEOSs):
        """
        plot several EOSs by listing them in <severalEOSs> with dictionaries
        see get_severalEOS_data for the format

        EXAMPLE INPUT
        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        all_label, all_headline, all_data = self.get_severalEOS_data(severalEOSs)

        fig, all_axes = self._get_figure(1,1,self._1by1_grid_placement)

        ax = all_axes[0]

        self._set_parms(ax, "$M/M_{\odot}$", "$J 10^{45} [g cm^3]$")

        EOS_markers = self._get_specific_ms( list(
            set( [ _["name"] for _ in severalEOSs ] )
        ) )

        m_colors = self._get_specific_c( list(
            set( [ _["m"] for _ in severalEOSs ] )
        ) )

        lambda_linestyle = self._get_specific_ls( list(
            set( [ _["lambda"] for _ in severalEOSs ] )
        ) )

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            data[-1] = [ _*self.units["J"]*1e-45 for _ in data[-1] ]

            #~ old way it set random at every iiteration for every entry
            #~ ls, lc, ms, mc = self._get_ls_lc_ms_mc()
            #~ label entry for each line, but now it will specific
            #~ label_entry = "{}"
                #~ "\n\t $\\beta$ = {:.1f}"
                #~ "\n\t m = {:.1e}"
                #~ "\n\t $\\lambda$ = {:.1e}".format(
                #~ eos["name"], eos["beta"], eos["m"], eos["lambda"]
            #~ )

            ax.plot(
                data[2],
                data[-1],
                label = None,
                color = m_colors.get(eos["m"], None),
                linestyle = lambda_linestyle.get(eos["lambda"], None),
                marker = EOS_markers.get(eos["name"], None),
                markerfacecolor = m_colors.get(eos["m"], None),
                markeredgecolor = m_colors.get(eos["m"], None),
                markersize = 5.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[2])
            )
        legend_eos = [
            Line2D(
                [0], [0],
                color="b",
                marker = __,
                linewidth = 0,
                label = _
            )
            for _, __ in EOS_markers.items()
        ]

        legend_lambda = [
            Line2D(
                [0], [0],
                color="b",
                marker = None,
                linewidth = 1.5,
                linestyle = __,
                label = "$\\lambda =$ {:.2e}".format(_)
            )
            for _, __ in lambda_linestyle.items()
        ]

        legend_m = [
            Line2D(
                [0], [0],
                color=__,
                marker = None,
                linewidth = 1.5,
                linestyle = "-",
                label = "$m =$ {:.2e}".format(_)
            )
            for _, __ in m_colors.items()
        ]

        ax.legend(
            handles = [*legend_eos, *legend_lambda, *legend_m],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol=3,
            frameon = False,
        )

        plt.show()

        return

    def convert_to_fitting(self, severalEOSs, fname = "STT_phiScal_J"):
        """
        for the provided list of dics of EOSs go over their results and create, by
        appending [name of result]_tildeI and [name of result]_barI, the
        neaceassery ceofficients for the fitting

        IT WILL OVERWRITE EXISTING !!!!

        EXAMPLE
        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        for eos in severalEOSs:

            result = os.path.join(
                self.my_ResPath,
                eos["name"],
                "_".join( [
                    fname,
                    eos["name"],
                    "beta{:.3e}".format(eos["beta"]),
                    "m{:.3e}".format(eos["m"]),
                    "lambda{:.3e}".format(eos["lambda"])
                ] )
            )

            with open(result, "r") as f:
                src_data_lines = f.readlines()

            #~ first line is just headline
            src_data_lines.pop(0)

            current_convert = "_tildeI"

            target_path = os.path.join(
                self.my_ResPath, eos["name"], "Fitting"
            )

            pathlib.Path( target_path ).mkdir(parents=True, exist_ok=True)

            target = os.path.join(
                target_path,
                os.path.basename(result) + current_convert
            )

            print(
                "\n will convert \n\t from {} \n\t to {} \n\t as {}".format(
                    result, target, current_convert
                )
            )

            with open( target, "w" ) as f:

                f.write("# M/R I/(MR**2) \n")

                for line in src_data_lines:

                    if not line.strip():
                        continue

                    tmp = [
                        float(_) for _ in line.strip().split(" ") if len(_.strip())
                    ]

                    f.write(
                        "{:.6e} {:.6e} \n".format(
                            tmp[2]/tmp[3], tmp[5]/(tmp[2]*tmp[3]**2)
                        )
                    )

            current_convert = "_barI"

            target = os.path.join(
                target_path,
                os.path.basename(result) + current_convert
            )


            print(
                "\n will convert \n\t from {} \n\t to {} \n\t as {}".format(
                    result, target, current_convert
                )
            )

            with open( target, "w" ) as f:

                f.write("# (M/R)**-1 I/M**3 \n")

                for line in src_data_lines:

                    if not line.strip():
                        continue


                    tmp = [
                        float(_) for _ in line.strip().split(" ") if len(_.strip())
                    ]

                    f.write(
                        "{:.6e} {:.6e} \n".format(
                            tmp[2]/tmp[3], tmp[5]/(tmp[2]**3)
                        )
                    )

        return

    def plot_severalEOSs_uniTildeI(self, severalEOSs ):
        """
        plot severalEOS unifersal Tilde I relationships
        <severalEOSs> with dictionaries see get_severalEOS_data for the format

        EXAMPLE INPUT

        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        all_label, all_headline, all_data = self.get_severalEOS_uniTildeI_data(severalEOSs)

        fig, all_axes = self._get_figure(1,1,self._1by1_grid_placement)

        ax = all_axes[0]

        self._set_parms(ax, "M/R", "$I/(MR^2)$")

        EOS_markers = self._get_specific_ms( list(
            set( [ _["name"] for _ in severalEOSs ] )
        ) )

        m_colors = self._get_specific_c( list(
            set( [ _["m"] for _ in severalEOSs ] )
        ) )

        lambda_linestyle = self._get_specific_ls( list(
            set( [ _["lambda"] for _ in severalEOSs ] )
        ) )

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ old way it set random at every iiteration for every entry
            #~ ls, lc, ms, mc = self._get_ls_lc_ms_mc()

            #~ label entry for each line, but now it will specific
            #~ label_entry = "{}"
                #~ "\n\t $\\beta$ = {:.1f}"
                #~ "\n\t m = {:.1e}"
                #~ "\n\t $\\lambda$ = {:.1e}".format(
                #~ eos["name"], eos["beta"], eos["m"], eos["lambda"]
            #~ )

            ax.plot(
                data[0],
                data[1],
                label = None,
                color = m_colors.get(eos["m"], None),
                linestyle = lambda_linestyle.get(eos["lambda"], None),
                marker = EOS_markers.get(eos["name"], None),
                markerfacecolor = m_colors.get(eos["m"], None),
                markeredgecolor = m_colors.get(eos["m"], None),
                markersize = 5.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[1])
            )

        legend_eos = [
            Line2D(
                [0], [0],
                color="b",
                marker = __,
                linewidth = 0,
                label = _
            )
            for _, __ in EOS_markers.items()
        ]

        legend_lambda = [
            Line2D(
                [0], [0],
                color="b",
                marker = None,
                linewidth = 1.5,
                linestyle = __,
                label = "$\\lambda =$ {:.2e}".format(_)
            )
            for _, __ in lambda_linestyle.items()
        ]

        legend_m = [
            Line2D(
                [0], [0],
                color=__,
                marker = None,
                linewidth = 1.5,
                linestyle = "-",
                label = "$m =$ {:.2e}".format(_)
            )
            for _, __ in m_colors.items()
        ]

        ax.legend(
            handles = [*legend_eos, *legend_lambda, *legend_m],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol=3,
            frameon = False,
        )

        plt.show()

        return

    def plot_severalEOSs_uniTildeI_polyFit(self, severalEOSs ):
        """
        plot severalEOS unifersal Tilde I relationships
        <severalEOSs> with dictionaries see get_severalEOS_data for the format

        EXAMPLE INPUT

        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        all_label, all_headline, all_data = self.get_severalEOS_uniTildeI_data(severalEOSs)

        #~ fig, all_axes = self._get_figure(1,1,self._1by1_grid_placement)
        fig, all_axes  = self._get_figure(
            2,  1,  self._3by1_shareX_grid_placement, height_ratios = [2,1]
        )

        ax_up = all_axes[0]
        ax_down = all_axes[1]

        #~ self._set_parms(ax_up, "M/R", "$I/(MR^3)$")
        self._set_parms(ax_up, "", r"$I/(MR^2)$")

        EOS_markers = self._get_specific_ms( list(
            set( [ _["name"] for _ in severalEOSs ] )
        ) )

        m_colors = self._get_specific_c( list(
            set( [ _["m"] for _ in severalEOSs ] )
        ) )

        lambda_linestyle = self._get_specific_ls( list(
            set( [ _["lambda"] for _ in severalEOSs ] )
        ) )

        max_x = 0
        min_x = 1e9

        all_colors = []

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ old way it set random at every iiteration for every entry
            #~ all_colors.append(self._get_ls_lc_ms_mc())

            #~ label entry for each line, but now it will specific
            #~ label_entry = "{}"
                #~ "\n\t $\\beta$ = {:.1f}"
                #~ "\n\t m = {:.1e}"
                #~ "\n\t $\\lambda$ = {:.1e}".format(
                #~ eos["name"], eos["beta"], eos["m"], eos["lambda"]
            #~ ),

            if max(data[0]) > max_x:
                max_x = max(data[0])

            if min(data[0]) < min_x:
                min_x = min(data[0])

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                color = m_colors.get(eos["m"], None),
                linestyle = lambda_linestyle.get(eos["lambda"], None),
                marker = EOS_markers.get(eos["name"], None),
                markerfacecolor = m_colors.get(eos["m"], None),
                markeredgecolor = m_colors.get(eos["m"], None),
                markersize = 5.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[1])
            )

        coef, rest  = polyfit(
            x = [ __ for _ in all_data for __ in _[0] ],
            y = [ __ for _ in all_data for __ in _[1] ],
            deg = [ 0, 1, 4 ],
            w = np.sqrt(np.array([ __ for _ in all_data for __ in _[1] ])),
            full = True,
        )

        p = lambda x: coef[0] + coef[1]*x + coef[4]*x**4

        p_x = np.linspace(min_x, max_x, 100)
        p_y = [ p(_) for _ in p_x ]

        lines_eos = [
            Line2D(
                [0], [0],
                color="b",
                marker = __,
                linewidth = 0,
                label = _
            )
            for _, __ in EOS_markers.items()
        ]

        lines_lambda = [
            Line2D(
                [0], [0],
                color="b",
                marker = None,
                linewidth = 1.5,
                linestyle = __,
                label = "$\\lambda =$ {:.2e}".format(_)
            )
            for _, __ in lambda_linestyle.items()
        ]

        lines_m = [
            Line2D(
                [0], [0],
                color=__,
                marker = None,
                linewidth = 1.5,
                linestyle = "-",
                label = "$m =$ {:.2e}".format(_)
            )
            for _, __ in m_colors.items()
        ]

        lines_polyfit = [
            Line2D(
                [0], [0],
                color = "#ff81c0",
                marker = None,
                linewidth = 2,
                linestyle = "-",
                label = "poly fit, $\chi_r^2$ = {:.3e}"
                    "\n {:.3f} + {:.3f}x + {:.3f}$x^4$".format(
                    rest[0][0]/(len([ __ for _ in all_data for __ in _[0] ]) - 3),
                    coef[0],
                    coef[1],
                    coef[4]
                )
            )
        ]

        ax_up.plot(
            p_x,
            p_y,
            color = "#ff81c0",
            marker = None,
            linewidth = 2,
            linestyle = "-",
            label = None
        )

        ax_up.get_shared_x_axes().join(ax_up, ax_down)
        ax_up.set_xticklabels([])

        ax_down.set_yscale("log")
        self._set_parms(ax_down, "M/R", r"$\left| 1 - \tilde I/\tilde I_{fit} \right|$  ")

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            _data = [
                abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
            ]

            ax_down.plot(
                data[0],
                _data,
                label = None,
                color = m_colors.get(eos["m"], None),
                linestyle = lambda_linestyle.get(eos["lambda"], None),
                marker = EOS_markers.get(eos["name"], None),
                markerfacecolor = m_colors.get(eos["m"], None),
                markeredgecolor = m_colors.get(eos["m"], None),
                markersize = 5.5,
                linewidth = 1.5,
                markevery = self._get_markevry(_data)
            )

        ax_up.legend(
            handles = [*lines_eos, *lines_lambda, *lines_m, *lines_polyfit],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol = 3,
            frameon = False,
            mode = None
        )

        plt.show()

        return

    def plot_severalEOSs_uniBarI(self, severalEOSs ):
        """
        plot severalEOS unifersal Tilde I relationships
        <severalEOSs> with dictionaries see get_severalEOS_data for the format

        EXAMPLE INPUT

        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        all_label, all_headline, all_data = self.get_severalEOS_uniBarI_data(severalEOSs)

        fig, all_axes = self._get_figure(1,1,self._1by1_grid_placement)

        ax = all_axes[0]

        self._set_parms(ax, "M/R", "$I/(M^3)$")

        EOS_markers = self._get_specific_ms( list(
            set( [ _["name"] for _ in severalEOSs ] )
        ) )

        m_colors = self._get_specific_c( list(
            set( [ _["m"] for _ in severalEOSs ] )
        ) )

        lambda_linestyle = self._get_specific_ls( list(
            set( [ _["lambda"] for _ in severalEOSs ] )
        ) )

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ old way it set random at every iiteration for every entry
            #~ ls, lc, ms, mc = self._get_ls_lc_ms_mc()

            #~ label entry for each line, but now it will specific
            #~ label_entry = "{}"
                #~ "\n\t $\\beta$ = {:.1f}"
                #~ "\n\t m = {:.1e}"
                #~ "\n\t $\\lambda$ = {:.1e}".format(
                #~ eos["name"], eos["beta"], eos["m"], eos["lambda"]
            #~ )

            ax.plot(
                data[0],
                data[1],
                label = None,
                color = m_colors.get(eos["m"], None),
                linestyle = lambda_linestyle.get(eos["lambda"], None),
                marker = EOS_markers.get(eos["name"], None),
                markerfacecolor = m_colors.get(eos["m"], None),
                markeredgecolor = m_colors.get(eos["m"], None),
                markersize = 5.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[1])
            )

        legend_eos = [
            Line2D(
                [0], [0],
                color="b",
                marker = __,
                linewidth = 0,
                label = _
            )
            for _, __ in EOS_markers.items()
        ]

        legend_lambda = [
            Line2D(
                [0], [0],
                color="b",
                marker = None,
                linewidth = 1.5,
                linestyle = __,
                label = "$\\lambda =$ {:.2e}".format(_)
            )
            for _, __ in lambda_linestyle.items()
        ]

        legend_m = [
            Line2D(
                [0], [0],
                color=__,
                marker = None,
                linewidth = 1.5,
                linestyle = "-",
                label = "$m =$ {:.2e}".format(_)
            )
            for _, __ in m_colors.items()
        ]

        ax.legend(
            handles = [*legend_eos, *legend_lambda, *legend_m],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol=3,
            frameon = False,
        )

        plt.show()

        return

    def plot_severalEOSs_uniBarI_polyFit(self, severalEOSs ):
        """
        plot severalEOS unifersal Tilde I relationships
        <severalEOSs> with dictionaries see get_severalEOS_data for the format

        EXAMPLE INPUT

        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        all_label, all_headline, all_data = self.get_severalEOS_uniBarI_data(severalEOSs)

        #~ fig, all_axes = self._get_figure(1,1,self._1by1_grid_placement)
        fig, all_axes  = self._get_figure(
            2,  1,  self._3by1_shareX_grid_placement, height_ratios = [2,1]
        )

        ax_up = all_axes[0]
        ax_down = all_axes[1]

        self._set_parms(ax_up, "", "$I/(M^3)$")

        EOS_markers = self._get_specific_ms( list(
            set( [ _["name"] for _ in severalEOSs ] )
        ) )

        m_colors = self._get_specific_c( list(
            set( [ _["m"] for _ in severalEOSs ] )
        ) )

        lambda_linestyle = self._get_specific_ls( list(
            set( [ _["lambda"] for _ in severalEOSs ] )
        ) )

        max_x = 0
        min_x = 1e9

        all_colors = []

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ old way it set random at every iiteration for every entry
            #~ all_colors.append(self._get_ls_lc_ms_mc())
            #~ old way it set random at every iiteration for every entry
            #~ label_entry =  "{}"
                #~ "\n\t $\\beta$ = {:.1f}"
                #~ "\n\t m = {:.1e}"
                #~ "\n\t $\\lambda$ = {:.1e}".format(
                #~ eos["name"], eos["beta"], eos["m"], eos["lambda"]
            #~ )

            if max(data[0]) > max_x:
                max_x = max(data[0])

            if min(data[0]) < min_x:
                min_x = min(data[0])

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                color = m_colors.get(eos["m"], None),
                linestyle = lambda_linestyle.get(eos["lambda"], None),
                marker = EOS_markers.get(eos["name"], None),
                markerfacecolor = m_colors.get(eos["m"], None),
                markeredgecolor = m_colors.get(eos["m"], None),
                markersize = 5.5,
                linewidth = 1.5,
                markevery = self._get_markevry(data[1])
            )

        coef, rest = polyfit(
            x = [ __**-1 for _ in all_data for __ in _[0] ],
            y = [ __ for _ in all_data for __ in _[1] ],
            deg = [ 1, 2, 3, 4 ],
            w = np.sqrt(np.array([ __ for _ in all_data for __ in _[1] ])),
            full = True
        )

        p = lambda x: coef[1]*x**-1 + coef[2]*x**-2 + coef[3]*x**-3 + coef[4]*x**-4

        p_x = np.linspace(min_x, max_x, 100)
        p_y = [ p(_) for _ in p_x ]

        lines_eos = [
            Line2D(
                [0], [0],
                color="b",
                marker = __,
                linewidth = 0,
                label = _
            )
            for _, __ in EOS_markers.items()
        ]

        lines_lambda = [
            Line2D(
                [0], [0],
                color="b",
                marker = None,
                linewidth = 1.5,
                linestyle = __,
                label = "$\\lambda =$ {:.2e}".format(_)
            )
            for _, __ in lambda_linestyle.items()
        ]

        lines_m = [
            Line2D(
                [0], [0],
                color=__,
                marker = None,
                linewidth = 1.5,
                linestyle = "-",
                label = "$m =$ {:.2e}".format(_)
            )
            for _, __ in m_colors.items()
        ]

        lines_polyfit = [
            Line2D(
                [0], [0],
                color = "#ff81c0",
                marker = None,
                linewidth = 2,
                linestyle = "-",
                label = "poly fit, $\chi_r^2$ = {:.3e}"
                    "\n {:.3f}$(1/x)^1$ + {:.3f}$(1/x)^2$ + {:.3f}$(1/x)^3$ + {:.3f}$(1/x)^4$".format(
                        rest[0][0]/(len([ __ for _ in all_data for __ in _[0] ]) - 4),
                        coef[1],
                        coef[2],
                        coef[3],
                        coef[4]
                )
            )
        ]

        ax_up.plot(
            p_x,
            p_y,
            color = "#ff81c0",
            marker = None,
            linewidth = 2,
            linestyle = "-",
            label = None
        )

        ax_up.get_shared_x_axes().join(ax_up, ax_down)
        ax_up.set_xticklabels([])

        ax_down.set_yscale("log")
        self._set_parms(ax_down, "M/R", r"$\left| 1 - \bar I/\bar I_{fit} \right|$")

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            _data = [
                abs(1 - _/p(__)) for _, __ in zip(data[1], data[0])
            ]

            ax_down.plot(
                data[0],
                _data,
                label = None,
                color = m_colors.get(eos["m"], None),
                linestyle = lambda_linestyle.get(eos["lambda"], None),
                marker = EOS_markers.get(eos["name"], None),
                markerfacecolor = m_colors.get(eos["m"], None),
                markeredgecolor = m_colors.get(eos["m"], None),
                markersize = 5.5,
                linewidth = 1.5,
                markevery = self._get_markevry(_data)
            )

        #~ ax_up.set_xlim(0.1)

        ax_up.legend(
            handles = [*lines_eos, *lines_lambda, *lines_m, *lines_polyfit],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol = 3,
            frameon = False,
            mode = None
        )

        plt.show()

        return

    def foo(self,severalEOSs):
        self.plot_severalEOSs_MvsR(severalEOSs)
        self.plot_severalEOSs_phiScal_cVSp_c(severalEOSs)
        return

    def foo2(self,severalEOSs):
        self.plot_latest_resEOSname_severalEOSs(severalEOSs)
        return

    def foo_uniI(self,severalEOSs):
        self.plot_severalEOSs_uniTildeI_polyFit(severalEOSs)
        self.plot_severalEOSs_uniBarI_polyFit(severalEOSs)
        return

    @staticmethod
    def _get_figure(nrows, ncols, grid_placement, height_ratios=None):
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

        #~ style.use("seaborn-poster")
        gs = GridSpec(nrows=nrows, ncols=ncols, height_ratios = height_ratios)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        fig = plt.figure()
        fig.set_tight_layout(True)

        return fig, grid_placement(fig, gs)

    @staticmethod
    def _4by4_grid_placement(fig, gs):
        """
        create grid which looks like

        +---------------------------------+
        |                  |              |
        |       Y0 vs X0   |   Y1 vs X1   |
        |       ax[0]      |   ax[1]      |
        |------------------+--------------|
        |                  |              |
        |  Y2 vs X2        |   Y3 vs X3   |
        |  ax[2]           |   ax[3]      |
        +---------------------------------+
        """
        return [
            fig.add_subplot( gs[0,0]),
            fig.add_subplot( gs[0,1]),
            fig.add_subplot( gs[1,0]),
            fig.add_subplot( gs[1,1]),
        ]

    @staticmethod
    def _1by1_grid_placement(fig, gs):
        """
        create grid which looks like

        +----------------------+
        |                      |
        |       Y0 vs X0       |
        |                      |
        +----------------------+
        """
        return [
            fig.add_subplot(gs[:]),
        ]

    @staticmethod
    def _3by1_shareX_grid_placement(fig, gs):
        """
        create grid which looks like ( proportion 1 : 2 )

        +----------------------+
        |                      |
        |       Y0 vs X0       |
        |                      |
        |                      |
        |                      |
        |                      |
        +----------------------+ shared Ox
        |                      |
        |       Y1 vs X1       |
        |                      |
        +----------------------+
        """

        up = fig.add_subplot( gs[0:1,0] )
        down = fig.add_subplot( gs[1,0] )

        return [ up, down  ]

    @staticmethod
    def _set_parms(
        ax, label_x, label_y, fontsize=10, x_ticksize = 8, y_ticksize = 8,
        x_format_str = "%.1e", y_format_str = "%.1e"
    ):
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

        if label_x:
            ax.set_xlabel(label_x, fontsize=fontsize)

        if x_ticksize:
            ax.xaxis.set_tick_params(labelsize=x_ticksize)

        if x_format_str:
            ax.xaxis.set_major_formatter(FormatStrFormatter(x_format_str))

        if label_y:
            ax.set_ylabel(label_y, fontsize=fontsize)

        if y_ticksize:
            ax.yaxis.set_tick_params(labelsize=y_ticksize)

        if y_format_str:
            ax.yaxis.set_major_formatter(FormatStrFormatter(y_format_str))

        return

    @staticmethod
    def _get_parameter_values(label):

        val_eosName = ""
        val_beta = 0
        val_m = 0
        val_lambda = 0

        for _ in label.split("_"):

            if "beta" in _:
                val_beta = float(_.replace("beta", ""))
            elif "lambda" in _:
                val_lambda = float(_.replace("lambda", ""))
            elif "m" in _:
                val_m = float(_.replace("m", ""))
            elif "STT" not in _ and "phiScal" not in _ and "J" not in _:
                val_eosName = _

        return val_eosName, val_beta, val_m, val_lambda

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
        units["R"] = 1e-3 * const_g * const_msun / const_c**2

        # units of moment of inertia
        units["J"] = 1e7 * const_g**2 * const_msun**3 / const_c**4

        return units

    @staticmethod
    def _get_ls_lc_ms_mc():
        """
        get tuple for
            ls ---> line style
            lc ---> line color
            ms ---> marker style
            mc----> marker color
        """

        markerstyles = [
            "s", ">", "<", "^", "v", "o", "X", "P", "d", "D", "H", "*", "p",
        ]

        linestyles = [
            ":", "-.", "--", "-"
        ]

        #~ colors = [
            #~ "b", "g", "r", "c", "m", "y", "k"
        #~ ]
        colors = [
            "b", "g", "r", "c", "y", "k"
        ]

        return random.sample(
            set(
                itertools.product(linestyles,colors,markerstyles, colors)
            ), 1
        )[0]

    #~ @staticmethod
    def _get_specific_ms(self, map_me):
        """
        for provided list return list of dictionaries, each having the
        entry of the provided list as key and marker style as value

        it is here because to be easier to find
        """

        if self.specific_ms:
            return self.specific_ms

        all_makrers_styles = [
            "s", "8", ">", "<", "^", "v", "o",
            "X", "P", "d", "D", "H", "h", "*", "p",
        ]

        if len(map_me) > len(all_makrers_styles):
            print(
                "\n not enough markers, amount of markers are {}\n".format(
                    len(all_makrers_styles)
                )
            )

            return

        #~ shuffle the list just for fun
        for _ in range(5):
            random.shuffle(all_makrers_styles)

        return {
            _: __ for _, __ in zip(map_me, all_makrers_styles)
        }

    #~ @staticmethod
    def _get_specific_ls(self, map_me):
        """
        for provided list return list of dictionaries, each having the
        entry of the provided list as key and line style as value

        it is here because to be easier to find
        """

        if self.specific_ls:
            return self.specific_ls

        all_line_styles = [
            ":", "-.", "--"
        ]

        if len(map_me) > len(all_line_styles):
            print(
                "\n not enough markers, amount of lines are {}\n".format(
                    len(all_line_styles)
                )
            )

            return

        #~ shuffle the list just for fun
        for _ in range(5):
            random.shuffle(all_line_styles)

        return {
            _: __ for _, __ in zip(map_me, all_line_styles)
        }

    #~ @staticmethod
    def _get_specific_c(self, map_me):
        """
        for provided list return list of dictionaries, each having the
        entry of the provided list as key and colour as value

        it is here because to be easier to find
        """

        if self.specific_c:
            return self.specific_c

        all_colors = [
            "b", "g", "r", "c", "m", "y", "k"
        ]

        if len(map_me) > len(all_colors):
            print(
                "\n not enough markers, amount of lines are {}\n".format(
                    len(all_line_styles)
                )
            )

            return

        #~ shuffle the list just for fun
        for _ in range(5):
            random.shuffle(all_colors)

        return {
            _: __ for _, __ in zip(map_me, all_colors)
        }

    #~ @staticmethod
    #~ def _get_val(ldict, key):
        #~ """
        #~ for the provided ldict list of dictionaries return the value of the
        #~ key, or None if not found dictionary with such key
        #~ """

        #~ for _ in map(lambda _: _.get(key, None), ldict):
            #~ if _:
                #~ return _
            #~ else:
                #~ continue

    @staticmethod
    def _get_markevry(data, amount_points = 20):
        """
        for provided data list get the difference between the max and min to
        evluate the stem for getting amount_points and return list of
        indexes of points closes to the step

        solution used from
        https://stackoverflow.com/questions/9873626/choose-m-evenly-spaced-elements-from-a-sequence-of-length-n
        """

        #~ f = lambda m, n: [i*n//m + n//(2*m) for i in range(m)]

        #~ step = abs(max(data) - min(data))/amount_points

        #~ idx = [ 0 ]
        #~ last = data[0]
        #~ for _, __ in enumerate(data):
            #~ if abs(__ - last) > step:
                #~ idx.append(_)
                #~ last = __

        #~ idx.append(len(data) - 1)
        #~ return f(amount_points, len(data))

        #~ abs(max(data) - min(data))/amount_points

        #~ def _get_pairwise(_iter):
            #~ "s -> (s0,s1), (s1,s2), (s2, s3), ..."
            #~ a, b = itertools.tee(_iter)
            #~ next(b, None)
            #~ return zip(a, b)

        #~ data_diff = np.diff(data)

        #~ idxs = [ 0 ]
        #~ for _, __ in enumerate(data_diff):
            #~ if __ * data_diff[idxs[-1]] > 0:
                #~ continue
            #~ else:
                #~ idxs.append(_)

        #~ idxs_p = [
            #~ int(abs(_[0] - _[1])/len(data)*amount_points)
            #~ for _ in _get_pairwise(idxs)
        #~ ]



        #~ idx_max = data.index(max(data, key=abs))

        #~ evenly_l = np.linspace(
            #~ data[0], data[idx_max], num = amount_points/2, endpoint = True
        #~ )

        #~ idx_l = [
            #~ data.index(
                #~ min(data[:idx_max], key=lambda x: abs(x-_))
            #~ ) for _ in evenly_l
        #~ ]

        #~ evenly_r = np.linspace(
            #~ data[idx_max], data[-1], num = amount_points/2, endpoint = True
        #~ )

        #~ idx_r = [
            #~ data.index(
                #~ min(data[idx_max:], key=lambda x: abs(x-_))
            #~ ) for _ in evenly_r
        #~ ]

        #~ return [ *idx_l, *idx_r ]

        evenly_spaced = np.linspace(
            min(data, key=abs), max(data, key=abs),
            num = amount_points,
            endpoint = True
        )

        return [
            data.index(min(data, key=lambda x: abs(x-_)))
            for _ in evenly_spaced
        ]


#################################################################################
##### PAST VERSION TO BE DELETED ################################################
#################################################################################
    #@staticmethod
    #def _check_if_StrNotBlank(string):
        #"""
        #check if a sting is blank/empty

        #Parameters
        #----------

        #Returns
        #-------
        #: boolean
            #True if string is not blank/empty
            #False if string is blank/empty
        #"""

        #return bool(string and string.strip())

    #@staticmethod
    #def _get_max_xy_coord(x, y):
        #"""
        #return the maximum coordinates of axes; this method is receiving
        #the axes class method get_?lim(), who is an list with the min/max
        #values

        #Parameters
        #----------
        #x: matplotlib.axes.get_xlim method
            #list of the min and max value of the abscissa of an axes

        #y: matplotlib.axes.get_ylim method
            #list of the min and max value of the ordinate of an axes

        #Returns
        #-------
        #max_x, max_y: tuple
            #the max values of x and y of an axes
        #"""

        #max_y = max(y, key=abs)
        #max_x = max(x, key=abs)

        #return max_x, max_y

    #@staticmethod
    #def _get_min_xy_coord(x, y):
        #"""
        #return the maximum coordinates of axes; this method is receiving
        #the axes class method get_?lim(), who is an list with the min/max
        #values

        #Parameters
        #----------
        #x: matplotlib.axes.get_xlim method
            #list of the min and max value of the abscissa of an axes

        #y: matplotlib.axes.get_ylim method
            #list of the min and max value of the ordinate of an axes

        #Returns
        #-------
        #max_x, max_y: tuple
            #the max values of x and y of an axes
        #"""

        #max_y = min(y, key=abs)
        #max_x = min(x, key=abs)

        #return max_x, max_y

    #@staticmethod
    #def _get_center_coordinates_axes(ax):
        #"""
        #return the center coordinates of an axes

        #Parameters
        #----------
        #ax: matplotlib.axes.Axes class

        #Returns:
        #center_x, center_y: double
            #the center coordinates of an axes <ax>
        #"""
        #max_x, max_y = get_max_xy_coord(ax.get_xlim(), ax.get_ylim())
        #min_x, min_y = get_min_xy_coord(ax.get_xlim(), ax.get_ylim())

        #center_x = (max_x + min_x) / 2
        #center_y = (max_y + min_y) / 2

        #return center_x, center_y

    #@staticmethod
    #def _get_iter_ls_c():
        #"""
        #get iterator who provides all possible combinations of predefined
        #list of line styles and colour to be used for plotting

        #the lists are shuffled before creating the generator

        #the generator cycles, thus it will never end

        #Parameters
        #----------

        #Returns
        #-------
        #: cycle object
            #generator who will repeat the product(styles_line, styles_colours)
            #indefinitely
        #"""
        #from itertools import product
        #from itertools import cycle
        #from random import shuffle

        #styles_line = ["-", "--", "-."]
        #styles_colours = ["b", "g", "r", "c", "m", "y", "k"]

        #shuffle(styles_line)
        #shuffle(styles_colours)

        #return cycle(product(styles_line, styles_colours))

    #@staticmethod
    #def _set_single_grid_placement(fig, gs):
        #"""
        #create single axes in <fig> using GridSpec <gs>
        #and fig.add_subplot() method

        #Parameters
        #---------
        #fig: Figure
            #the figure whos geometry will be specified

        #gs: class
            #predefined class of how many rows and columns the
            #figure should have

        #Returns
        #------
        #all_ax: list
            #list containing all axes defined by the <gs> in <fig>
        #"""

        #all_ax = []
        #all_ax.append(fig.add_subplot(gs[:]))

        #return all_ax

    #@staticmethod
    #def _get_index(max_index):
        #"""
        #get the user input choice of an index and return it
        #by making sure it is an positive integer not greater than <max_index>

        #Parameters
        #----------
        #max_index: int
            #the maximum index available

        #Return
        #------
        #index: int
            #the index which the user chose
        #"""

        #while "waiting for the number":
            #try:
                #index = int(input("\n\t\t index =... "))
            #except ValueError:
                #print("\n Not a number is you typed \n")
                #continue

            #if index < 0 or index >= max_index:
                #print("\n {} not in allowed boundaries \n".format(index))
            #else:
                #return index

    #@staticmethod
    #def _get_yes_or_no_answer():
        #"""
        #get the answer of yes-no question

        #Parameters
        #----------

        #Return
        #------
        #"""

        #while "Answer is invalid":
            #ans = str(input("\n\t y/n answer... ")).lower().strip()

            #print("\n")

            #if ans[:1] == "y":
                #return True
            #elif ans[:1] == "n":
                #return False
            #else:
                #print("\n y or n answer only !\n")
                #continue

    #def _get_latest_file(self, path, model):
        #"""
        #return the name of the latest modified file in path
        #by creating wildcard of the type model + "*"

        #Parameters
        #----------
        #path: string
            #the path to the directory of interest
        #model: string
            #how the file name looks like, it serves as wildcard
            #of the type model + "*"

        #Returns
        #-------
        #: string
            #the name of latest modified file
        #"""

        #import glob
        #import os

        #if (self._check_if_StrNotBlank(path) and
            #self._check_if_StrNotBlank(model)):

            #filename = max(
                #glob.glob(path + "/" + model + "*"),
                #key=os.path.getctime
            #)

            #filename = filename.split("/")[-1]

            #print(
                #"\n latest file in \n\t {} is: {} \n".format(
                    #path, filename
                #)
            #)
        #else:

            #filename = ""

            #print("\n path and name not set, some filename var!!! \n")

        #return filename

    #def _get_xy_current_file(self, headline):
        #"""
        #list the enumerated <headlines> and let the user choose
        #which will be plotted as Y and as X on the fig

        #Parameters
        #----------
        #headline: list
            #list of available headlines to choose from

        #Return
        #------
        #index_x: int
            #the index of the headline to be used on the abscissa

        #index_y: int
            #the index of the headline to be used on the ordinate
        #"""

        #print("\n choose columns to plot Y vs X \n")

        #for cnt, value in enumerate(headline):
            #print("\t {}: {}".format(cnt, value))

        #print("\n\t Y... ")
        #index_y = self._get_index(len(headline))

        #print("\n\t X... ")
        #index_x = self._get_index(len(headline))

        #return index_x, index_y

    #def _get_filename_from_dir(self, path):
        #"""
        #lists all the files in <path> and lets the user choose
        #which file to be loaded

        #Parameters
        #----------
        #path: string
            #the full path to the directory whos files will be listed

        #Return
        #------
            #: string
                #the filename of our choice
        #"""

        #import os

        #list_dir = os.listdir(path)

        #for ind, val in enumerate(list_dir):
            #print("{}. {}".format(ind, val))

        #return list_dir[self._get_index(len(list_dir))]

    #def _set_default_paths(self):
        #"""
        #sets default paths in self to the result directories

        #Parameters
        #----------

        #Returns
        #-------
        #"""

        #self.path_res = "/home/dimitar/projects/STT_theories/results"

        #print(
            #"\n mine path_res: {} \n".format(
                #self.path_res
            #)
        #)

        #self.path_daniela = "/home/dimitar/Documents/Teaching_Materials/" \
            #+ "University/Uvod.Fizika.Cherni.Dupki/Doktorant/" \
            #+ "Daniela_Static_NS_Massive_SlowRot"

        #print(
            #"\n daniela path_res: {} \n".format(
                #self.path_daniela
            #)
        #)

        #self.path_kalin = "/home/dimitar/Documents/Teaching_Materials/" \
            #+ "University/Uvod.Fizika.Cherni.Dupki/Doktorant/" \
            #+ "Kalin_Static_NS_SlowRot_Massive_Lambda/beta-6"

        #print(
            #"\n kalin path_res: {} \n".format(
                #self.path_kalin
            #)
        #)

        #return

    #def _load_file(self, path, filename, data, headline, label):
        #"""
        #load the data from <filename> in <path> by appending them to
        #<data>, append the headline (the name of each column) in
        #<headline> and the name of the file in <label>

        #Parameters
        #----------
        #path: string
            #the path to the directory containing the file with data
        #filename: string
            #the name of the file which will be loaded
        #data: list
            #the loaded data will be appended here
        #headline: list
            #the name of each column will be appended here as a list
        #label: list
            #the name of the file will be apended here

        #Returns
        #------
        #"""

        #if not(self._check_if_StrNotBlank(path) and
               #self._check_if_StrNotBlank(filename)):

            #print(
                #"\n from where to load a file!!! \n\t {} / {} \n".format(
                    #path, filename
                #)
            #)

        #with open(path + "/" + filename, "r") as f:
            #all_data = f.readlines()

        #headline.append(
            #[
                #i.strip() for i in
                #all_data.pop(0).strip().split(" ")
                #if
                #"#" not in i and
                #len(i.strip())
            #]
        #)

        #label.append(filename.split("/")[-1][-41:])

        #data.append(
            #[
                #[] for i in all_data[0].strip().split(" ") if len(i.strip())
            #]
        #)

        #for cnt, line in enumerate(all_data):

            #for apnd, num in zip(
                #data[-1], [
                    #float(i) for i in line.strip().split(" ") if len(i.strip())
                #]
            #):
                #apnd.append(num)

        #return

    #def _load_daniela_file(self, path, filename, data, headline, label):
        #"""
        #load danieala result file with name <filename> by appending the data
        #after taking care of coefficients to <data> and saving the <headlines>(
        #what is the name of each column) and label - the filename itself

        #Parameters
        #----------
        #path: string
            #the full path of the directory where the file is
        #filename: string
            #the name of the file itself
        #data: list
            #where to append the data of interest
        #headline: list
            #what each column is
        #label: list
            #the name the given data, expect it will be the filename

        #Returns
        #-------
        #"""

        #if not(self._check_if_StrNotBlank(path) and
               #self._check_if_StrNotBlank(filename)):

            #print(
                #"from where to load a daniela file!!! {} / {}".format(
                    #path, filename
                #)
            #)

        #with open(path + "/" + filename, "r") as f:
            #all_data = f.readlines()

        ##~ get rid of inline text
        #all_data = [
            #line for line in all_data
            #if "rho" not in line
        #]

        ##~ i am not interested in all data, so only the indexes here
        ##~ will be saved
        ##~ since in some files the central pressure is in different column
        ##~ we check it and take into account
        #if filename in [
            #"models_APR_beta-4.5_mphi5e-3.dat",
            #"models_APR_beta-6.0_mphi0.dat",
            #"models_APR_beta-6.0_mphi1e-3.dat"
        #]:
            #indexes = [0, 1, 2, 3, 4, 10]
        #else:
            #indexes = [0, 1, 2, 3, 4, 11]

        #units = [
            #self.units["density"], self.units["rad"],
            #1, self.units["j"], 1, 1
        #]

        #data.append([[] for i in indexes])

        #for cnt, line in enumerate(all_data):
            #for i, u, d in zip(indexes, units, data[-1]):
                #try:
                    #d.append(float(line.split(" ")[i]) / u)
                #except ValueError:
                    #print(
                        #"\n ValueError: line {}: {} \n".format(
                            #cnt, line
                        #)
                    #)
                    #break

        #data[-1][-2] = \
            #[(-1)*i for i in data[-1][-2]]

        #label.append(filename)

        #return

    #def _load_kalin_file(self, path, filename, data, headline, label):
        #"""
        #load kalin result file with name <filename> by appending the data
        #after taking care of coefficients to <data> and saving the <headlines>(
        #what is the name of each column) and label - the filename itself

        #Parameters
        #----------
        #path: string
            #the full path of the directory where the file is
        #filename: string
            #the name of the file itself
        #data: list
            #where to append the data of interest
        #headline: list
            #what each column is
        #label: list
            #the name the given data, expect it will be the filename

        #Returns
        #-------
        #"""

        #if not(self._check_if_StrNotBlank(path) and
               #self._check_if_StrNotBlank(filename) ):

            #print(
                #"from where to load a kalin file!!! {} / {}".format(
                    #path, filename
                #)
            #)

        #with open(path + "/" + filename, "r") as f:
            #all_data = f.readlines()

        ##~ get rid of inline text
        #all_data = [
            #line for line in all_data
            #if "lambda" not in line and "f_rot" not in line
        #]

        ##~ i am not interested in all data, so only the indexes here
        ##~ will be saved
        ##~ since in some files the central pressure is in different column
        ##~ we check it and take into account
        #indexes = [3, 5, 7, 10, 6, 2]

        #units = [
            #self.units["density"], self.units["rad"],
            #1, 1, 1, 1
        #]

        #data.append([[] for i in indexes])

        #for cnt, line in enumerate(all_data):
            #for i, u, d in zip(indexes, units, data[-1]):
                #try:
                    #d.append(float(line.split(" ")[i]) / u)
                #except ValueError:
                    #print(
                        #"\n ValueError: line {}: {} \n".format(
                            #cnt, line
                        #)
                    #)
                    #break

        #label.append(filename)

        #return

    #def load_ResultFile_latest(self, filename_type="STT_phiScal_"):
        #"""
        #load the latest modified file from <self.path_res>

        #Parameters
        #----------
        #filename_type: string
            #this string will be used as wildcard filename_type + "*"
            #to search in self.path_res; if no value is supplied, the
            #default "STT_phiScal_" will be used

        #Return
        #------
        #"""

        #fname_res = self._get_latest_file(
            #self.path_res, filename_type
        #)

        #print("\n now will load:\n\t {} \n".format(fname_res))

        #self._load_file(
            #self.path_res,
            #fname_res,
            #self.data_res,
            #self.headline_res,
            #self.label_res
        #)

        #return

    #def load_ResultFile_argument(self, filename_type=""):
        #"""
        #load <filename_type> and append the results in self.data_res and etc
        #it will search the file at <self.path>

        #Parameters
        #----------
        #filename_type: string
            #the file name at <self.path> which will be loaded

        #Return
        #------
        #"""

        #if self._check_if_StrNotBlank(filename_type):

            #self._load_file(
                #self.path_res,
                #filename_type,
                #self.data_res,
                #self.headline_res,
                #self.label_res
            #)

        #else:
            #print("\n you forgot to give filename \n")

        #return

    #def load_ResultFile(self):
        #"""
        #load <filename_type> and append the results in self.data_res and etc
        #it will search the file at <self.path>

        #Parameters
        #----------
        #filename_type: string
            #the file name at <self.path> which will be loaded

        #Return
        #------
        #"""

        #filename = self._get_filename_from_dir(self.path_res)

        #print("\n now will load:\n\t {} \n".format(filename))

        #self._load_file(
            #self.path_res,
            #filename,
            #self.data_res,
            #self.headline_res,
            #self.label_res
        #)

        #return

    #def load_DanielaFile(self):
        #"""
        #load daniela data from <self.path_daniela>

        #Parameters
        #----------

        #Returns
        #-------
        #"""

        #filename = self._get_filename_from_dir(self.path_daniela)

        #print("\n now will load:\n\t {} \n".format(filename))

        #self._load_daniela_file(
            #self.path_daniela,
            #filename,
            #self.data_daniela,
            #self.headline_daniela,
            #self.label_daniela
        #)

    #def load_KalinFile(self):
        #"""
        #load kalin data from <self.path_kalin>

        #Parameters
        #----------

        #Returns
        #-------
        #"""

        #filename = self._get_filename_from_dir(self.path_kalin)

        #print("\n now will load:\n\t {} \n".format(filename))

        #self._load_kalin_file(
            #self.path_kalin,
            #filename,
            #self.data_kalin,
            #self.headline_kalin,
            #self.label_kalin
        #)

    #def clear_data(self):
        #"""
        #clear the data in self:
            #data_*, which contains the plot data
            #headline_*, which is the naming of each column of each data
            #label_*, which is name of the file the data come from

        #Parameters
        #----------

        #Returns
        #-------
        #"""

        #self.data_res.clear()
        #self.headline_res.clear()
        #self.label_res.clear()

        #self.data_daniela.clear()
        #self.headline_daniela.clear()
        #self.label_daniela.clear()

        #self.data_kalin.clear()
        #self.headline_kalin.clear()
        #self.label_kalin.clear()

        #return

    #def plot_data_res_last(self):
        #"""
        #creates figure with single subplot to plot the last list in
        #<self.data_res> with asking which column will be used for
        #X and Y axis

        #Parameters
        #----------

        #Return
        #------
        #"""

        #from matplotlib import pyplot as plt
        #from IPython import get_ipython

        ##~ if we start it in jupyter qtconsole
        ##~ this will force open it in new window not inside
        #get_ipython().run_line_magic("matplotlib", "qt5")

        ##~ lets check if we have anything to plot
        #if not self.data_res:
            #self.load_ResultFile_latest()

        #index_x, index_y = self._get_xy_current_file(
            #self.headline_res[-1]
        #)

        #fig, all_ax = self._get_figure(1, 1, self._set_single_grid_placement)

        #ax = all_ax[-1]

        #self._set_parms(
            #ax,
            #self.headline_res[-1][index_x],
            #self.headline_res[-1][index_y]
        #)

        #ax.plot(
            #self.data_res[-1][index_x],
            #self.data_res[-1][index_y],
            #linewidth=1.5,
            #label=self.label_res[-1]
        #)

        #ax.legend(loc="best", fontsize=8)
        #plt.show()

        #return

    #def plot_data_res_all(self):
        #"""
        #creates figure with single subplot to plot all data in
        #<self.data_res> with asking which column will be used for all
        #X and Y axis (it presumes that all share the same headline)

        #Parameters
        #----------

        #Return
        #------
        #"""

        #from matplotlib import pyplot as plt
        #from IPython import get_ipython

        ##~ if we start it in jupyter qtconsole
        ##~ this will force open it in new window not inside
        #get_ipython().run_line_magic('matplotlib', "qt5")

        ##~ lets check if we have anything to plot
        #if not self.data_res:
            #self.load_ResultFile()

        #index_x, index_y = self._get_xy_current_file(
            #self.headline_res[-1]
        #)

        #fig, all_ax = self._get_figure(1, 1, self._set_single_grid_placement)

        #ax = all_ax[-1]

        #self._set_parms(
            #ax,
            #self.headline_res[-1][index_x],
            #self.headline_res[-1][index_y]
        #)

        #for data, label in zip(self.data_res, self.label_res):
            #ax.plot(
                #data[index_x],
                #data[index_y],
                #linewidth=1.5,
                #label=label
            #)

        #ax.legend(loc="best", fontsize=8)
        #plt.show()

        #return

    #def plot_compare_last_res_daniela(self):

        #from matplotlib import pyplot as plt
        #from IPython import get_ipython

        #get_ipython().run_line_magic('matplotlib', "qt5")

        ##~ lets check if we have anything to plot
        #if not self.data_res:
            #self.load_ResultFile()

        #if not self.data_daniela:
            #self.load_DanielaFile()

        #index_x, index_y = self._get_xy_current_file(self.headline_res[-1])

        #fig, all_ax = self._get_figure(1, 1, self._set_single_grid_placement)

        #ax = all_ax[-1]

        #self._set_parms(
            #ax,
            #self.headline_res[-1][index_x],
            #self.headline_res[-1][index_y]
        #)

        #ax.plot(
            #self.data_res[-1][index_x],
            #self.data_res[-1][index_y],
            #linewidth=1.5,
            #label=self.label_res[-1],
        #)

        #ax.plot(
            #self.data_daniela[-1][
                #self.daniela_mapping[self.headline_res[-1][index_x]]
            #],
            #self.data_daniela[-1][
                #self.daniela_mapping[self.headline_res[-1][index_y]]
            #],
            #marker="o",
            #markersize=5,
            #alpha=0.4,
            #linewidth=0,
            #label=self.label_daniela[-1]
        #)

        #ax.legend(loc="best", fontsize=8)
        #plt.show()

        #return

    #def plot_compare_last_res_daniela_kalin(self):

        #from matplotlib import pyplot as plt
        #from IPython import get_ipython

        #get_ipython().run_line_magic('matplotlib', "qt5")

        ##~ lets check if we have anything to plot
        #if not self.data_res:
            #self.load_ResultFile()

        #if not self.data_daniela:
            #self.load_DanielaFile()

        #if not self.data_kalin:
            #self.load_KalinFile()

        #index_x, index_y = self._get_xy_current_file( self.headline_res[-1] )

        #fig, all_ax = self._get_figure(1, 1, self._set_single_grid_placement)

        #ax = all_ax[-1]

        #self._set_parms(
            #ax,
            #self.headline_res[-1][index_x],
            #self.headline_res[-1][index_y]
        #)

        #ax.plot(
            #self.data_res[-1][index_x],
            #self.data_res[-1][index_y],
            #linewidth=1.5,
            #label=self.label_res[-1],
        #)

        #ax.plot(
            #self.data_daniela[-1][
                #self.daniela_mapping[self.headline_res[-1][index_x]]
            #],
            #self.data_daniela[-1][
                #self.daniela_mapping[self.headline_res[-1][index_y]]
            #],
            #marker="o",
            #markersize=5,
            #alpha=0.4,
            #linewidth=0,
            #label=self.label_daniela[-1]
        #)

        #ax.plot(
            #self.data_kalin[-1][
                #self.kalin_mapping[self.headline_res[-1][index_x]]
            #],
            #self.data_kalin[-1][
                #self.kalin_mapping[self.headline_res[-1][index_y]]
            #],
            #marker="X",
            #markersize=5,
            #alpha=0.4,
            #linewidth=0,
            #label=self.label_kalin[-1]
        #)

        #ax.legend(loc="best", fontsize=8)
        #plt.show()

        #return

    #def plot_compare_last_res_kalin(self):

        #from matplotlib import pyplot as plt
        #from IPython import get_ipython

        #get_ipython().run_line_magic('matplotlib', "qt5")

        ##~ lets check if we have anything to plot
        #if not self.data_res:
            #self.load_ResultFile()

        #if not self.data_kalin:
            #self.load_KalinFile()

        #index_x, index_y = self._get_xy_current_file(self.headline_res[-1])

        #fig, all_ax = self._get_figure(1, 1, self._set_single_grid_placement)

        #ax = all_ax[-1]

        #self._set_parms(
            #ax,
            #self.headline_res[-1][index_x],
            #self.headline_res[-1][index_y]
        #)

        #ax.plot(
            #self.data_res[-1][index_x],
            #self.data_res[-1][index_y],
            #linewidth=1.5,
            #label=self.label_res[-1],
        #)

        #ax.plot(
            #self.data_kalin[-1][
                #self.kalin_mapping[self.headline_res[-1][index_x]]
            #],
            #self.data_kalin[-1][
                #self.kalin_mapping[self.headline_res[-1][index_y]]
            #],
            #marker="X",
            #markersize=5,
            #alpha=0.5,
            #linewidth=0,
            #label=self.label_kalin[-1]
        #)

        #ax.legend(loc="best", fontsize=8)
        #plt.show()

        #return

    #def remove_from_data_res(self):
        #"""
        #user says who from <self.data_res> based on their labels saved in
        #<self.label_res>

        #Parameters
        #----------

        #Returns
        #-------
        #"""

        #while "want to remove something":
            #print("\n Want to remove some data from self.data_res? \n")

            #if not self._get_yes_or_no_answer():
                #return

            #for cnt, val in enumerate(self.label_res):

                #print("{}. {} ".format(cnt, val))

            #print("\n Choose index... ")
            #i = self._get_index(len(self.label_res))

            #del self.data_res[i]
            #del self.label_res[i]
            #del self.headline_res[i]


if __name__ == "__main__":

    print("\n asd \n")
#~ cd projects/STT_theories/src/python/
#~ import imp
#~ import plotting as myplt_M
#~ severalEOSs = [
    #~ {'name': 'SLy4', 'beta': 0, 'm': 0, 'lambda': 0},
    #~ {'name': 'SLy4', 'beta': -6, 'm': 0, 'lambda': 1}
#~ ]
#~ myplt = myplt_M.plot_result()
#~ myplt.set_my_ResPath()
#~ myplt.plot_severalEOSs_MvsR(severalEOSs)
#~ myplt.plot_severalEOSs_phiScal_cVSp_c(severalEOSs)

#~ STEP LIKE M
#~ eos = "APR4"
#~ severalEOSs = [
    #~ {'name': eos, 'beta': 0, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 0},
#~ ]
#~ myplt.foo(severalEOSs)
#~ severalEOSs = [
    #~ {'name': eos, 'beta': 0, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 5e-3, 'lambda': 0},
#~ ]
#~ myplt.foo(severalEOSs)
#~ severalEOSs = [
    #~ {'name': eos, 'beta': 0, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 5e-3, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 1e-2, 'lambda': 0},
#~ ]
#~ myplt.foo(severalEOSs)
#~ severalEOSs = [
    #~ {'name': eos, 'beta': 0, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 5e-3, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 1e-2, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 5e-2, 'lambda': 0},
#~ ]
#~ myplt.foo(severalEOSs)

#~ STEP LIKE LAMBDA
#~ eos = "APR4"
#~ severalEOSs = [
    #~ {'name': eos, 'beta': 0, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 0},
#~ ]
#~ myplt.foo(severalEOSs)
#~ severalEOSs = [
    #~ {'name': eos, 'beta': 0, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 1e-1},
#~ ]
#~ myplt.foo(severalEOSs)
#~ severalEOSs = [
    #~ {'name': eos, 'beta': 0, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 1e-1},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 1e0},
#~ ]
#~ myplt.foo(severalEOSs)
#~ severalEOSs = [
    #~ {'name': eos, 'beta': 0, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 1e-1},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 1e0},
    #~ {'name': eos, 'beta': -6, 'm': 0, 'lambda': 1e1},
#~ ]
#~ myplt.foo(severalEOSs)
