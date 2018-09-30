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

#~ get_ipython().run_line_magic("matplotlib", "qt")
#~ import matplotlib
#~ matplotlib.use('Qt5Agg')

class plot_result:

    def __init__(self):

        self.my_ResPath = None
        self.my_EOSname = None
        self.my_fname_starts = None
        self_my_file = None
        self.my_headline = None
        self.my_data = None
        self.my_label = None

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

        return

    def _luminosity_color(self, color, amount=0.5):
        """
        Lightens the given color by multiplying (1-luminosity) by the given amount.
        Input can be matplotlib color string, hex string, or RGB tuple.

        Examples:
        >> lighten_color('g', 0.3)
        >> lighten_color('#F034A3', 0.6)
        >> lighten_color((.3,.55,.1), 0.5)
        """
        import matplotlib.colors as mc
        import colorsys
        try:
            c = mc.cnames[color]
        except:
            c = color

        c = colorsys.rgb_to_hls(*mc.to_rgb(c))

        return colorsys.hls_to_rgb(abs(c[0]), abs(1 - amount * (1 - c[1])), abs(c[2]))

    def set_my_ResPath(
        self,
        my_ResPath = "~/projects/STT_theories/results/"
    ):
        """
        path to the results of the shootings
        """

        self.my_ResPath = os.path.expanduser(my_ResPath)

        return

    def set_severalEOSs_ms_ls_c(
        self, severalEOSs, m_lambda_style = { "ls": "lambda", "c": "m" }
    ):
        """
        for consistent marking on all graphs for the current instance

        provide which will be the marker for m and lambda
        EOS is always the marker style with dictionary

        {"ls": "lambda", "c": "m"} will map linestyle to lambda and color to m

        EXAMPLE INPUT
        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        self.specific_ms = None
        self.specific_ms = self._get_specific_ms( list(
            set( [ _["name"] for _ in severalEOSs ] )
        ), "name" )

        self.specific_ls = None
        self.specific_ls = self._get_specific_ls( list(
            set( [ _[m_lambda_style["ls"]] for _ in severalEOSs ] )
        ), m_lambda_style["ls"] )

        self.specific_c = None
        self.specific_c = self._get_specific_c( list(
            set( [ _[m_lambda_style["c"]] for _ in severalEOSs ] )
        ), m_lambda_style["c"] )

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
        fname = "STT_phiScal_J",
        append = ""
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
                "tildeI",
                append
            ] ) if append else "_".join( [
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
        fname = "STT_phiScal_J",
        append = ""
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
                "barI",
                append
            ] ) if append else "_".join( [
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

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        for label, data, eos in zip( all_label, all_data, severalEOSs ):
            data[3] = [ _*self.units["R"] for _ in data[3] ]

            ax.plot(
                data[3],
                data[2],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(data[3], data[2]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        ax.legend(
            handles = [*lines_markers, *lines_colors, *lines_linestyles],
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

    def plot_severalEOSs_MvsR_stable(self, severalEOSs):
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

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ we are interested only in stable modules - part of the graphs
            #~ where the mass is increasing
            #~ after the maximum mass we can cut them off
            #~ the min mass we are interested is the one which is  at least 0.5 M sun
            min_mass = 0.5
            little_offset = 2

            max_m_i = data[2].index(max(data[2]))

            max_m_i = max_m_i + little_offset \
            if max_m_i + little_offset < len(data[2]) - 1 else max_m_i

            min_m_i = list(
                map(lambda _: _ >= min_mass,  data[2][:max_m_i])
            ).index(True)

            min_m_i = min_m_i - little_offset \
            if min_m_i - little_offset > 0 else min_m_i

            stable_R = [ _*self.units["R"] for _ in data[3][min_m_i:max_m_i] ]

            ax.plot(
                stable_R,
                data[2][min_m_i:max_m_i],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(stable_R, data[2][min_m_i:max_m_i]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        ax.legend(
            handles = [*lines_markers, *lines_colors, *lines_linestyles],
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

    def plot_severalEOSs_MvsR_GR(self, severalEOSs):
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

        all_label_GR, all_headline_GR, all_data_GR = self.get_severalEOS_data( [
                { "name": _, "beta": 0, "m": 0, "lambda": 0 }
                for _ in set( [ _["name"] for _ in severalEOSs ] )
            ]
        )

        fig, all_axes = self._get_figure(1,1,self._1by1_grid_placement)

        ax = all_axes[0]

        self._set_parms(ax, "R [km]", "$M/M_{\odot}$")

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ we are interested only in stable modules - part of the graphs
            #~ where the mass is increasing
            #~ after the maximum mass we can cut them off
            #~ the min mass we are interested is the one which is  at least 0.5 M sun
            min_mass = 0.5
            little_offset = 2

            max_m_i = data[2].index(max(data[2]))

            max_m_i = max_m_i + little_offset \
            if max_m_i + little_offset < len(data[2]) - 1 else max_m_i

            min_m_i = list(
                map(lambda _: _ >= min_mass,  data[2][:max_m_i])
            ).index(True)

            min_m_i = min_m_i - little_offset \
            if min_m_i - little_offset > 0 else min_m_i

            stable_R = [ _*self.units["R"] for _ in data[3][min_m_i:max_m_i] ]

            ax.plot(
                stable_R,
                data[2][min_m_i:max_m_i],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(stable_R, data[2][min_m_i:max_m_i]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

            _tmp = self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )

            ax.plot(
                [ _*self.units["R"] for _ in data[3][max_m_i:] ],
                data[2][max_m_i:],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(
                    [ _*self.units["R"] for _ in data[3][max_m_i:] ],
                    data[2][max_m_i:],
                    amount_points=5
                ),
                marker = _tmp["marker"],
                markerfacecolor = _tmp["color"],
                markeredgecolor = _tmp["color"],
                color = _tmp["color"],
                linestyle = (0, (5,10)),
            )

        #~ GR_color_markers = "#ef4026"
        GR_color_markers = "#c0022f"
        #~ GR_color_markers = "#a9f971"
        #~ GR_color_fit = "#ed0dd9"
        GR_color_fit = GR_color_markers

        for label, data, eos in zip(
            all_label_GR,
            all_data_GR,
            [ _ for _ in set( [ _["name"] for _ in severalEOSs ] ) ]
        ):

            #~ we are interested only in stable modules - part of the graphs
            #~ where the mass is increasing
            #~ after the maximum mass we can cut them off
            #~ the min mass we are interested is the one which is  at least 0.5 M sun
            min_mass = 0.5
            little_offset = 2

            #~ max_m_i = data[2].index(max(data[2]))
            max_m_i = data[2].index(data[2][-1])

            max_m_i = max_m_i + little_offset \
            if max_m_i + little_offset < len(data[2]) - 1 else max_m_i

            min_m_i = list(
                map(lambda _: _ >= min_mass,  data[2][:max_m_i])
            ).index(True)

            min_m_i = min_m_i - little_offset \
            if min_m_i - little_offset > 0 else min_m_i

            stable_R = [ _*self.units["R"] for _ in data[3][min_m_i:max_m_i] ]

            ax.plot(
                stable_R,
                data[2][min_m_i:max_m_i],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(
                    stable_R,
                    data[2][min_m_i:max_m_i]
                ),
                marker = markers.get(eos, None),
                color = GR_color_markers,
                markerfacecolor = GR_color_markers,
                markeredgecolor = GR_color_markers,
            )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles, severalEOSs
        )

        #~ ax.set_xlim(9,15)
        #~ ax.set_ylim(0.5,3.2)

        ax.legend(
            handles = [
                *lines_markers, *lines_colors, *lines_linestyles,
                Line2D(
                    [0], [0], color = GR_color_fit, marker = None, linestyle
                    = "-", linewidth = 1.5, label = "GR"
                )
            ],
            loc="best",
            fontsize=8,
            handlelength=3,
            numpoints=1,
            fancybox=True,
            markerscale = 1.25,
            ncol = 4,
            frameon = False,
            mode = None
        )

        plt.show()

        return

    def plot_severalEOSs_MvsR_GR_Ifunc(self, severalEOSs):
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

        all_label_GR, all_headline_GR, all_data_GR = self.get_severalEOS_data( [
                { "name": _, "beta": 0, "m": 0, "lambda": 0 }
                for _ in set( [ _["name"] for _ in severalEOSs ] )
            ]
        )

        fig, all_axes = self._get_figure(1,1,self._1by1_grid_placement)

        ax = all_axes[0]

        self._set_parms(ax, "R [km]", "$M/M_{\odot}$")

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ we are interested only in stable modules - part of the graphs
            #~ where the mass is increasing
            #~ after the maximum mass we can cut them off
            #~ the min mass we are interested is the one which is  at least 0.5 M sun
            min_mass = 0.5
            little_offset = 2

            max_m_i = data[2].index(max(data[2]))

            max_m_i = max_m_i + little_offset \
            if max_m_i + little_offset < len(data[2]) - 1 else max_m_i

            min_m_i = list(
                map(lambda _: _ >= min_mass,  data[2][:max_m_i])
            ).index(True)

            min_m_i = min_m_i - little_offset \
            if min_m_i - little_offset > 0 else min_m_i

            stable_R = [ _*self.units["R"] for _ in data[3][min_m_i:max_m_i] ]

            ax.plot(
                stable_R,
                data[2][min_m_i:max_m_i],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(stable_R, data[2][min_m_i:max_m_i]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

            _tmp = self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )

            ax.plot(
                [ _*self.units["R"] for _ in data[3][max_m_i:] ],
                data[2][max_m_i:],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(
                    [ _*self.units["R"] for _ in data[3][max_m_i:] ],
                    data[2][max_m_i:],
                    amount_points = 5
                ),
                marker = _tmp["marker"],
                markerfacecolor = _tmp["color"],
                markeredgecolor = _tmp["color"],
                color = _tmp["color"],
                linestyle = (0, (5,10)),
            )

        #~ GR_color_markers = "#ef4026"
        GR_color_markers = "#c0022f"
        #~ GR_color_markers = "#a9f971"
        #~ GR_color_fit = "#ed0dd9"
        GR_color_fit = GR_color_markers

        for label, data, eos in zip(
            all_label_GR,
            all_data_GR,
            [ _ for _ in set( [ _["name"] for _ in severalEOSs ] ) ]
        ):

            #~ we are interested only in stable modules - part of the graphs
            #~ where the mass is increasing
            #~ after the maximum mass we can cut them off
            #~ the min mass we are interested is the one which is  at least 0.5 M sun
            min_mass = 0.5
            little_offset = 2

            #~ max_m_i = data[2].index(max(data[2]))
            max_m_i = data[2].index(data[2][-1])

            max_m_i = max_m_i + little_offset \
            if max_m_i + little_offset < len(data[2]) - 1 else max_m_i

            min_m_i = list(
                map(lambda _: _ >= min_mass,  data[2][:max_m_i])
            ).index(True)

            min_m_i = min_m_i - little_offset \
            if min_m_i - little_offset > 0 else min_m_i

            stable_R = [ _*self.units["R"] for _ in data[3][min_m_i:max_m_i] ]

            ax.plot(
                stable_R,
                data[2][min_m_i:max_m_i],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(
                    stable_R, data[2][min_m_i:max_m_i]
                ),
                marker = markers.get(eos, None),
                color = GR_color_markers,
                markerfacecolor = GR_color_markers,
                markeredgecolor = GR_color_markers,
            )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        #~ reduce the uncessasery markers, colors and labels
        #~ as they may be not be part of severalEOSs
        lines_markers = [
            _ for _ in lines_markers if _.get_label() in [
                _["name"] for _ in severalEOSs
            ]
        ]

        lines_colors = [
            _ for _ in lines_colors if _.get_color() in [
                colors[_[colors["label"]]] for _ in severalEOSs
            ]
        ]

        lines_linestyles = [
            _ for _ in lines_linestyles if str(_.get_linestyle()) in [
                str(linestyles[_[linestyles["label"]]]) for _ in severalEOSs
            ]
        ]

        ax.set_xlim(9,15)
        ax.set_ylim(0.5,3.2)

        ax.legend(
            handles = [
                *lines_markers, *lines_colors, *lines_linestyles,
                Line2D(
                    [0], [0], color = GR_color_fit, marker = None, linestyle
                    = "-", linewidth = 1.5, label = "GR fit"
                )
            ],
            loc="best",
            fontsize=8,
            handlelength=2.5,
            numpoints=1,
            fancybox=True,
            markerscale = 1.25,
            ncol = 4,
            frameon = False,
            mode = None
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

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ ls, lc, ms, mc = self._get_ls_lc_ms_mc()

            ax.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        ax.legend(
            handles = [*lines_markers, *lines_colors, *lines_linestyles],
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

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            data[-2] = [ _*self.units["density"] for _ in data[-2] ]

            ax.plot(
                data[-2],
                data[1],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(data[-2],data[1]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        ax.legend(
            handles = [*lines_markers, *lines_colors, *lines_linestyles],
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

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            data[-2] = [ _*self.units["density"] for _ in data[-2] ]

            ax.plot(
                data[-2],
                data[2],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(data[-2], data[2]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        ax.legend(
            handles = [*lines_markers, *lines_colors, *lines_linestyles],
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

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ the last multiplication is added as it is scaled out in the plot
            data[-1] = [ _*self.units["J"]*1e-45 for _ in data[-1] ]

            ax.plot(
                data[2],
                data[-1],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(data[2],data[-1]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        ax.legend(
            handles = [*lines_markers, *lines_colors, *lines_linestyles],
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

    def convert_to_fitting_stable(self, severalEOSs, fname = "STT_phiScal_J"):
        """
        for the provided list of dics of EOSs go over their results and create, by
        appending [name of result]_tildeI and [name of result]_barI, the
        neaceassery ceofficients for the fitting

        IT WILL OVERWRITE EXISTING !!!!

        Since I want to convert only the stable models I will get only those
        entries in the data which have increasing Mass,
        Strictly speaking - those entries up to the line with maximum mass

        EXAMPLE
        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        def _get_max_min_M_i(src_data_lines):
            """
            load all netries of column 2 - the mass
            find the line with maximum mass and use only those
            since we want only stable models
            """

            min_mass = 0.5

            stable_m = []
            for line in src_data_lines:
                stable_m.append(
                    [
                        float(_) for _ in line.strip().split(" ") if len(_.strip())
                    ][2]
                )

            little_offset = 2

            max_m_i = stable_m.index(max(stable_m))

            max_m_i = max_m_i + little_offset \
            if max_m_i + little_offset < len(stable_m) - 1 else max_m_i

            min_m_i = list(
                map(lambda _: _ >= min_mass, stable_m[:max_m_i])
            ).index(True)

            min_m_i = min_m_i - little_offset \
            if min_m_i - little_offset > 0 else min_m_i

            return max_m_i, min_m_i

        #~ min_M*Mass of Sun
        min_M = 0.5
        print(
            "\n Since the minimum measured mass is 1 times the mass of Sun"
            " will cut out masses smaller than {} \n".format(min_M)
        )

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

            convert_tilde = "_tildeI_stable"
            convert_bar = "_barI_stable"
            convert_dir = "Fitting"

            target_path = os.path.join(
                self.my_ResPath, eos["name"], convert_dir
            )
            pathlib.Path( target_path ).mkdir(parents=True, exist_ok=True)

            with open(result, "r") as src, \
            open(
                os.path.join( target_path, os.path.basename(result) + convert_tilde ),
                "w"
            ) as dst_tilde, \
            open(
                os.path.join( target_path, os.path.basename(result) + convert_bar ),
                "w"
            ) as dst_bar:

                src_data_lines = src.readlines()

                #~ first line is just headline
                src_data_lines.pop(0)

                max_M_i, min_M_i = _get_max_min_M_i(src_data_lines)

                print(
                    "\n Now converting from {}"
                    "\n\t EOS {}"
                    "\n\t model {}"
                    "\n\t to tilde {}"
                    "\n\t to bar {}".format(
                        self.my_ResPath,
                        eos["name"],
                        result,
                        dst_tilde.name,
                        dst_bar.name
                    )
                )

                dst_tilde.write("# M/R I/(MR**2) \n")
                dst_bar.write("# (M/R)**-1 I/M**3 \n")

                for line in src_data_lines[min_M_i:max_M_i]:

                    if not line.strip():
                        continue

                    tmp = [
                        float(_) for _ in line.strip().split(" ") if len(_.strip())
                    ]

                    dst_tilde.write(
                        "{:.6e} {:.6e} \n".format(
                            tmp[2]/tmp[3], tmp[5]/(tmp[2]*tmp[3]**2)
                        )
                    )
                    dst_bar.write(
                        "{:.6e} {:.6e} \n".format(
                            tmp[2]/tmp[3], tmp[5]/(tmp[2]**3)
                        )
                    )
        print(
            "\n Since the minimum measured mass is 1 times the mass of Sun"
            " will cut out masses smaller than {} \n".format(min_M)
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

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            ax.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        ax.legend(
            handles = [*lines_markers, *lines_colors, *lines_linestyles],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1.25,
            ncol = 3,
            frameon = False,
            mode = None
        )

        ax_up.set_xlim(0.09)

        plt.show()

        return

    def plot_severalEOSs_uniTildeI_polyFit(
        self, severalEOSs, append_stable = "stable"
    ):
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

        _get_max = lambda data, _: max(data) if max(data) > _ else _
        _get_min = lambda data, _: min(data) if min(data) < _ else _

        all_label, all_headline, all_data = self.get_severalEOS_uniTildeI_data(
            severalEOSs, append = append_stable
        )

        fig, all_axes  = self._get_figure(
            2,  1,  self._3by1_shareX_grid_placement, height_ratios = [2,1]
        )

        ax_up = all_axes[0]
        ax_down = all_axes[1]

        self._set_parms(ax_up, "", r"$I/(MR^2)$")

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        max_x = max_y = 0
        min_x = min_y = 1e9
        min_compactness = 0.09

        #~ color_fit = "#fcc006"
        color_fit = "#fa4224"
        #~ color_avg = "#9f2305"
        #~ color_avg_worst = "#dfc5fe"

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ we set mimimal compactenss threshold and cut out all entries
            #~ who are below it only if we are interested in stable solutions
            if append_stable:
                _min_x = list(
                    map(lambda _: _ >= min_compactness, data[0])
                ).index(True)

                data[0] = [ _ for _ in data[0][_min_x:] ]
                data[1] = [ _ for _ in data[1][_min_x:] ]

            max_x = _get_max(data[0], max_x)
            min_x = _get_min(data[0], min_x)

            max_y = _get_max(data[1], max_y)
            min_y = _get_min(data[1], min_y)

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

        coef, rest  = polyfit(
            x = [ __ for _ in all_data for __ in _[0] ],
            y = [ __ for _ in all_data for __ in _[1] ],
            deg = [ 0, 1, 4 ],
            w = np.sqrt(np.array([ __ for _ in all_data for __ in _[1] ])),
            full = True,
        )
        chi_red = rest[0][0]/(len([ __ for _ in all_data for __ in _[0] ]) - 3)

        p = lambda x: coef[0] + coef[1]*x + coef[4]*x**4

        p_x = np.linspace(min_x, max_x, 100)
        p_y = [ p(_) for _ in p_x ]

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        #~ polyfit is on its own
        lines_polyfit = [
            Line2D(
                [0], [0],
                color = color_fit,
                marker = None,
                linewidth = 2,
                linestyle = "-",
                label = "poly fit, $\chi_r^2$ = {:.3e}"
                    "\n {:.3f} + {:.3f}x + {:.3f}$x^4$".format(
                    chi_red,
                    coef[0],
                    coef[1],
                    coef[4]
                )
            )
        ]

        ax_up.get_shared_x_axes().join(ax_up, ax_down)
        ax_up.set_xticklabels([])

        ax_down.set_yscale("log")
        self._set_parms(ax_down, "M/R", r"$\left| 1 - \tilde I/\tilde I_{fit} \right|$  ")

        max_y_down = 0
        min_y_down = 1e9

        #~ average over all EOSs of all the residuals
        delta_all = 0
        n_all = 0

        #~ average over all EOSs of largest residual
        delta_all_max = 0
        n_all_max = 0

        #~ the largest residual across all EOSs
        delta_max = 0

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            _data = [
                abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
            ]

            delta_all += sum(_data)
            n_all += len(_data)

            delta_all_max += max(_data)
            n_all_max += 1

            delta_max = _get_max(_data, delta_max)

            max_y_down = _get_max(_data, max_y_down)
            min_y_down = _get_min(_data, min_y_down)

            ax_down.plot(
                data[0],
                _data,
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], _data),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

        avg_L_1 = delta_all/n_all
        avg_L_inf = delta_all_max/n_all_max
        L_inf_worst = delta_max

        ax_up.plot(
            p_x,
            p_y,
            color = color_fit,
            marker = None,
            linewidth = 2,
            linestyle = "-",
            label = None
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*(1 + avg_L_inf),
            np.array(p_y)*(1 - avg_L_inf),
            facecolor=color_avg,
            alpha=0.75
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*( 1 + L_inf_worst ),
            np.array(p_y)*( 1 - L_inf_worst ),
            facecolor=color_avg_worst,
            alpha=0.5
        )

        ax_up.legend(
            handles = [
                *lines_markers, *lines_colors, *lines_linestyles, *lines_polyfit
            ],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1.25,
            ncol = 3,
            frameon = False,
            mode = None
        )

        ax_up.set_xlim(min_x, max_x)
        ax_up.set_ylim(min_y, max_y)

        ax_down.set_ylim(1e-3, 1e0)

        print(
            "\n All fit information"
            "\n\t $\chi_r^2$ = {:.3e}"
            "\n\t $a_0$ = {:.3e}"
            "\n\t $a_1$ = {:.3e}"
            "\n\t $a_4$ = {:.3e}"
            "\n\t $< L_1 >$ = {:.3e}"
            "\n\t $< L_\inf >$ = {:.3e}"
            "\n\t $ L_\inf $ = {:.3e}\n".format(
                chi_red,
                coef[0],
                coef[1],
                coef[4],
                avg_L_1,
                avg_L_inf,
                L_inf_worst
            )
        )

        if n_all_max != len(severalEOSs):
            print(
                "\n SOMETHING FISSHY, n_all = {}, len EOSs = {} \n".format(
                n_all_max, len(severalEOSs))
            )

        plt.show()

        return

    def plot_severalEOSs_uniTildeI_polyFitAll_GR(
        self, severalEOSs, append_stable = "stable"
    ):
        """
        plot severalEOS unifersal Tilde I relationships
        <severalEOSs> with dictionaries see get_severalEOS_data for the format

        EXAMPLE INPUT

        append_stable = "stable" is set as default, so only stable solutions
        will be calculated

        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        #~ expand the data and do the polyfit
        def _get_polyfit_res(xp, yp):
            xp = [ __ for _ in xp for __ in _ ]
            yp = [ __ for _ in yp for __ in _ ]

            coef, rest = polyfit(
                x = xp, y = yp,
                deg = [ 0, 1, 4 ],
                w = np.sqrt(yp),
                full = True
            )

            #~ calcualte the chi reduce
            chi_red = rest[0][0]/(len(xp) - 3)
            p = lambda x: coef[0] + coef[1]*x + coef[4]*x**4

            return coef, chi_red, p

        _get_max = lambda data, _: max(data) if max(data) > _ else _
        _get_min = lambda data, _: min(data) if min(data) < _ else _

        all_label, all_headline, all_data = self.get_severalEOS_uniTildeI_data(
            severalEOSs, append = append_stable
        )

        all_label_GR, all_headline_GR, all_data_GR = self.get_severalEOS_uniTildeI_data( [
                { "name": _, "beta": 0, "m": 0, "lambda": 0 }
                for _ in set( [ _["name"] for _ in severalEOSs ] )
            ],
            append = append_stable
        )

        fig, all_axes  = self._get_figure(
            2,  1,  self._3by1_shareX_grid_placement, height_ratios = [2,1]
        )

        ax_up = all_axes[0]
        ax_down = all_axes[1]

        ax_up.get_shared_x_axes().join(ax_up, ax_down)
        ax_up.set_xticklabels([])

        ax_down.set_yscale("log")
        self._set_parms(ax_down, "M/R", r"$\left| 1 - \tilde I/\tilde I_{fit} \right|$  ")

        self._set_parms(ax_up, "", r"$I/(MR^2)$")

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        max_x = max_y = 0
        min_x = min_y = 1e9
        min_compactness = 0.09

        #~ GR_color_markers = "#ef4026"
        GR_color_markers = "#c0022f"
        #~ GR_color_markers = "#a9f971"
        #~ GR_color_fit = "#ed0dd9"
        GR_color_fit = GR_color_markers

        plot_alpha = 0.75

        #~ lets plot severEOSs on the up plot and eventually cut out data
        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ we set mimimal compactenss threshold and cut out all entries
            #~ who are below it only if we are interested in stable solutions
            if append_stable:
                _min_x = list(
                    map(lambda _: _ >= min_compactness, data[0])
                ).index(True)

                data[0] = [ _ for _ in data[0][_min_x:] ]
                data[1] = [ _ for _ in data[1][_min_x:] ]

            max_x = _get_max(data[0], max_x)
            min_x = _get_min(data[0], min_x)

            max_y = _get_max(data[1], max_y)
            min_y = _get_min(data[1], min_y)

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                ),
                alpha = plot_alpha
            )

        #~ will use this foreach to fill the polyfit for GR
        polyfit_res = []
        xp = []
        yp = []

        #~ plot all GR data and gather all x and y for evaluating the polyfit
        for label, data, eos in zip(
            all_label_GR,
            all_data_GR,
            [ _ for _ in set( [ _["name"] for _ in severalEOSs ] ) ]
        ):

            #~ we set mimimal compactenss threshold and cut out all entries
            #~ who are below it only if we are interested in stable solutions
            if append_stable:
                _min_x = list(
                    map(lambda _: _ >= min_compactness, data[0])
                ).index(True)

                data[0] = [ _ for _ in data[0][_min_x:] ]
                data[1] = [ _ for _ in data[1][_min_x:] ]

            max_x = _get_max(data[0], max_x)
            min_x = _get_min(data[0], min_x)

            max_y = _get_max(data[1], max_y)
            min_y = _get_min(data[1], min_y)

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                marker = markers.get(eos, None),
                color = GR_color_markers,
                markerfacecolor = GR_color_markers,
                markeredgecolor = GR_color_markers,
                alpha = plot_alpha
            )

            xp.append(data[0])
            yp.append(data[1])

        coef, chi_red, p = _get_polyfit_res(xp, yp)

        max_y_down = 0
        min_y_down = 1e9

        #~ average over all EOSs of all the residuals
        delta_all = 0
        n_all = 0

        #~ average over all EOSs of largest residual
        delta_all_max = 0
        n_all_max = 0

        #~ the largest residual across all EOSs
        delta_max = 0

        #~ for the generated polyfit function calcualte the
        #~ relative error and plot it donw for GR
        for label, data, eos in zip(
            all_label_GR,
            all_data_GR,
            [ _ for _ in set( [ _["name"] for _ in severalEOSs ] ) ]
        ):

            _data = [
                abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
            ]

            delta_all += sum(_data)
            n_all += len(_data)

            delta_all_max += max(_data)
            n_all_max += 1

            delta_max = _get_max(_data, delta_max)

            max_y_down = _get_max(_data, max_y_down)
            min_y_down = _get_min(_data, min_y_down)

            ax_down.plot(
                data[0],
                _data,
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], _data),
                marker = markers.get(eos, None),
                color = GR_color_markers,
                markerfacecolor = GR_color_markers,
                markeredgecolor = GR_color_markers
            )

        avg_L_1 = delta_all/n_all
        avg_L_inf = delta_all_max/n_all_max
        L_inf_worst = delta_max

        print(
            "\n GR fit"
            "\n\t $\chi_r^2$ = {:.3e}"
            "\n\t $a_0$ = {:.3e}"
            "\n\t $a_1$ = {:.3e}"
            "\n\t $a_4$ = {:.3e}"
            "\n\t $< L_1 >$ = {:.3e}"
            "\n\t $< L_\inf >$ = {:.3e}"
            "\n\t $ L_\inf $ = {:.3e}\n".format(
                chi_red,
                coef[0],
                coef[1],
                coef[4],
                avg_L_1,
                avg_L_inf,
                L_inf_worst
            )
        )

        lines_polyfit = [
            Line2D(
                [0], [0],
                color = GR_color_fit,
                marker = None,
                linestyle = "-",
                linewidth = 1.5,
                label = "GR fit"
            )
        ]

        p_x = np.linspace(min_x, max_x, 100)
        p_y = [ p(_) for _ in p_x ]

        #~ generate 100 points between min and max of x and plot it values
        ax_up.plot(
            p_x,
            p_y,
            label = None,
            linewidth = 2,
            linestyle = "-",
            markersize = 0,
            markevery = 0,
            marker = None,
            color = GR_color_fit,
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*(1 + avg_L_inf),
            np.array(p_y)*(1 - avg_L_inf),
            facecolor=GR_color_markers,
            alpha= plot_alpha - 0.25
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*( 1 + L_inf_worst ),
            np.array(p_y)*( 1 - L_inf_worst ),
            facecolor=GR_color_markers,
            alpha= plot_alpha - 0.5
        )

        #~ now do the same for each color if there are more than 1
        for k, v in colors.items():

            #~ the colors will have label key containing the name of parameter
            #~ which they represetn
            if k == "label":
                continue

            xp = []
            yp = []

            for data, eos in zip(all_data, severalEOSs):
                #~ if the current eos has parameter value equal to the current one
                #~ lets append its data
                if eos[colors["label"]] == k:
                    xp.append( data[0] )
                    yp.append( data[1] )

            #~ expand all the data into flat list to calculate the polyfit
            coef, chi_red, p = _get_polyfit_res(xp, yp)

            max_y_down = 0
            min_y_down = 1e9

            #~ average over all EOSs of all the residuals
            delta_all = 0
            n_all = 0

            #~ average over all EOSs of largest residual
            delta_all_max = 0
            n_all_max = 0

            #~ the largest residual across all EOSs
            delta_max = 0

            for data, eos in zip(all_data, severalEOSs):
                #~ if the current eos has parameter value equal to the current one
                #~ lets append its data
                if eos[colors["label"]] == k:

                    _data = [
                        abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
                    ]

                    delta_all += sum(_data)
                    n_all += len(_data)

                    delta_all_max += max(_data)
                    n_all_max += 1

                    delta_max = _get_max(_data, delta_max)

                    max_y_down = _get_max(_data, max_y_down)
                    min_y_down = _get_min(_data, min_y_down)

                    ax_down.plot(
                        data[0],
                        _data,
                        label = None,
                        linewidth = 0,
                        markersize = 5.5,
                        markevery = self._get_markevry(data[0], _data),
                        marker = markers.get(eos["name"], None),
                        color = v,
                        markerfacecolor = v,
                        markeredgecolor = v
                    )

            avg_L_1 = delta_all/n_all
            avg_L_inf = delta_all_max/n_all_max
            L_inf_worst = delta_max

            print(
                "\n {} {:.3e} fit"
                "\n\t $\chi_r^2$ = {:.3e}"
                "\n\t $a_0$ = {:.3e}"
                "\n\t $a_1$ = {:.3e}"
                "\n\t $a_4$ = {:.3e}"
                "\n\t $< L_1 >$ = {:.3e}"
                "\n\t $< L_\inf >$ = {:.3e}"
                "\n\t $ L_\inf $ = {:.3e}\n".format(
                    "$\\lambda =$ " if colors["label"] == "lambda" else "m =",
                    k,
                    chi_red,
                    coef[0],
                    coef[1],
                    coef[4],
                    avg_L_1,
                    avg_L_inf,
                    L_inf_worst
                )
            )

            lines_polyfit.append(
                Line2D(
                    [0], [0],
                    color = v,
                    marker = None,
                    linestyle = "-",
                    linewidth = 1.5,
                    label = "{} {:.3e} fit".format(
                        "$\\lambda =$ " if colors["label"] == "lambda" else "m =",
                        k
                    )
                )
            )

            p_x = np.linspace(min_x, max_x, 100)
            p_y = [ p(_) for _ in p_x ]

            ax_up.plot(
                p_x,
                p_y,
                label = None,
                linewidth = 2,
                linestyle = "-",
                markersize = 0,
                markevery = 0,
                marker = None,
                color = v,
            )

            ax_up.fill_between(
                p_x,
                np.array(p_y)*(1 + avg_L_inf),
                np.array(p_y)*(1 - avg_L_inf),
                facecolor=v,
                alpha= plot_alpha - 0.25
            )

            ax_up.fill_between(
                p_x,
                np.array(p_y)*( 1 + L_inf_worst ),
                np.array(p_y)*( 1 - L_inf_worst ),
                facecolor=v,
                alpha= plot_alpha - 0.5
            )

        #####################################################################
        #~ THE CODE BELOW CAN WORK FOR LINSTEYLES DO NOT DELETE
        #####################################################################
        #~ now do the same for each linestyle if there are more than 1
        #for k, v in linestyles.items():

            ##~ the linestyle will have label key containing the name of parameter
            ##~ which they represetn
            #if k == "label":
                #continue

            #xp = []
            #yp = []

            #for data, eos in zip(all_data, severalEOSs):
                ##~ if the current eos has parameter value equal to the current one
                ##~ lets append its data
                #if eos[linestyles["label"]] == k:
                    #xp.append( data[0] )
                    #yp.append( data[1] )

            ##~ expand all the data into flat list to calculate the polyfit
            #coef, chi_red, p = _get_polyfit_res(xp, yp)

            #max_y_down = 0
            #min_y_down = 1e9

            ##~ average over all EOSs of all the residuals
            #delta_all = 0
            #n_all = 0

            ##~ average over all EOSs of largest residual
            #delta_all_max = 0
            #n_all_max = 0

            ##~ the largest residual across all EOSs
            #delta_max = 0

            #for data, eos in zip(all_data, severalEOSs):
                ##~ if the current eos has parameter value equal to the current one
                ##~ lets append its data
                #if eos[linestyles["label"]] == k:

                    #_data = [
                        #abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
                    #]

                    #delta_all += sum(_data)
                    #n_all += len(_data)

                    #delta_all_max += max(_data)
                    #n_all_max += 1

                    #delta_max = _get_max(_data, delta_max)

                    #max_y_down = _get_max(_data, max_y_down)
                    #min_y_down = _get_min(_data, min_y_down)

                    #ax_down.plot(
                        #data[0],
                        #_data,
                        #label = None,
                        #linewidth = 0,
                        #markersize = 5.5,
                        #markevery = self._get_markevry(data[0], _data),
                        #marker = markers.get(eos["name"], None),
                        #color = "k",
                        #markerfacecolor = "k",
                        #markeredgecolor = "k"
                    #)

            #avg_L_1 = delta_all/n_all
            #avg_L_inf = delta_all_max/n_all_max
            #L_inf_worst = delta_max

            #print(
                #"\n {} {:.3e} fit"
                #"\n\t $\chi_r^2$ = {:.3e}"
                #"\n\t $a_0$ = {:.3e}"
                #"\n\t $a_1$ = {:.3e}"
                #"\n\t $a_4$ = {:.3e}"
                #"\n\t $< L_1 >$ = {:.3e}"
                #"\n\t $< L_\inf >$ = {:.3e}"
                #"\n\t $ L_\inf $ = {:.3e}\n".format(
                    #"$\\lambda =$ " if linestyles["label"] == "lambda" else "m =",
                    #k,
                    #chi_red,
                    #coef[0],
                    #coef[1],
                    #coef[4],
                    #avg_L_1,
                    #avg_L_inf,
                    #L_inf_worst
                #)
            #)

            #lines_polyfit.append(
                #Line2D(
                    #[0], [0],
                    #color = "k",
                    #marker = None,
                    #linestyle = v,
                    #linewidth = 1.5,
                    #label = "{} {:.3e} fit".format(
                        #"$\\lambda =$ " if linestyles["label"] == "lambda" else "m =",
                        #k
                    #)
                #)
            #)

            #p_x = np.linspace(min_x, max_x, 100)
            #p_y = [ p(_) for _ in p_x ]

            #ax_up.plot(
                #p_x,
                #p_y,
                #label = None,
                #linewidth = 2,
                #linestyle = v,
                #markersize = 0,
                #markevery = 0,
                #marker = None,
                #color = "k",
            #)

            #ax_up.fill_between(
                #p_x,
                #np.array(p_y)*(1 + avg_L_inf),
                #np.array(p_y)*(1 - avg_L_inf),
                #facecolor="k",
                #alpha= plot_alpha - 0.25
            #)

            #ax_up.fill_between(
                #p_x,
                #np.array(p_y)*( 1 + L_inf_worst ),
                #np.array(p_y)*( 1 - L_inf_worst ),
                #facecolor="k",
                #alpha= plot_alpha - 0.5
            #)

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        #~ I think the code below is copied by mistake
        #~ commeting out for clearence
        #~ ax_up.legend(
            #~ handles = [
                #~ *lines_markers, *lines_colors, *lines_linestyles, *lines_polyfit
            #~ ],
            #~ loc="best",
            #~ fontsize=8,
            #~ handlelength=2,
            #~ numpoints=1,
            #~ fancybox=True,
            #~ markerscale = 1.25,
            #~ ncol = 3,
            #~ frameon = False,
            #~ mode = None
        #~ )

        ax_up.legend(
            handles = [
                *lines_markers, *lines_colors, *lines_polyfit
            ],
            loc="best",
            fontsize=8,
            handlelength=2.5,
            numpoints=1,
            fancybox=True,
            markerscale = 1.25,
            ncol = 4,
            frameon = False,
            mode = None
        )

        ax_up.set_xlim(min_x, max_x)
        ax_up.set_ylim(min_y, max_y)

        ax_down.set_ylim(1e-3, 1e0)

        plt.show()

        return

    def plot_severalEOSs_uniTildeI_polyFitAll_GR_ParmsProduct(
        self, severalEOSs, append_stable = "stable"
    ):
        """
        plot severalEOS unifersal Tilde I relationships
        <severalEOSs> with dictionaries see get_severalEOS_data for the format

        EXAMPLE INPUT

        append_stable = "stable" is set as default, so only stable solutions
        will be calculated

        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        #~ expand the data and do the polyfit
        def _get_polyfit_res(xp, yp):
            xp = [ __ for _ in xp for __ in _ ]
            yp = [ __ for _ in yp for __ in _ ]

            coef, rest = polyfit(
                x = xp, y = yp,
                deg = [ 0, 1, 4 ],
                w = np.sqrt(yp),
                full = True
            )

            #~ calcualte the chi reduce
            chi_red = rest[0][0]/(len(xp) - 3)
            p = lambda x: coef[0] + coef[1]*x + coef[4]*x**4

            return coef, chi_red, p

        _get_max = lambda data, _: max(data) if max(data) > _ else _
        _get_min = lambda data, _: min(data) if min(data) < _ else _

        all_label, all_headline, all_data = self.get_severalEOS_uniTildeI_data(
            severalEOSs, append = append_stable
        )

        all_label_GR, all_headline_GR, all_data_GR = self.get_severalEOS_uniTildeI_data( [
                { "name": _, "beta": 0, "m": 0, "lambda": 0 }
                for _ in set( [ _["name"] for _ in severalEOSs ] )
            ],
            append = append_stable
        )

        fig, all_axes  = self._get_figure(
            2,  1,  self._3by1_shareX_grid_placement, height_ratios = [2,1]
        )

        #~ turn off the tight layout since it does not permit sticking the two axes
        #~ but can do it manually
        fig.set_tight_layout(False)
        fig.set_rasterized(True)
        fig.subplots_adjust(wspace=0, hspace=0)

        ax_up = all_axes[0]
        ax_down = all_axes[1]

        ax_up.get_shared_x_axes().join(ax_up, ax_down)
        ax_up.set_xticklabels([])

        ax_down.set_yscale("log")
        self._set_parms(
            ax_down,
            "M/R",
            "$\\left| 1 - \\tilde I/\\tilde I_{{fit}} \\right|$"
        )

        self._set_parms(ax_up, "", r"$I/(MR^2)$")

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        max_x = max_y = 0
        min_x = min_y = 1e9
        min_compactness = 0.09

        #~ GR_color_markers = "#ef4026"
        GR_color_markers = "#469990"
        #~ GR_color_markers = "#a9f971"
        #~ GR_color_fit = "#ed0dd9"
        GR_color_fit = self._luminosity_color(GR_color_markers, 1.1)

        plot_alpha = 0.6

        #~ lets plot severEOSs on the up plot and eventually cut out data
        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ we set mimimal compactenss threshold and cut out all entries
            #~ who are below it only if we are interested in stable solutions
            if append_stable:
                _min_x = list(
                    map(lambda _: _ >= min_compactness, data[0])
                ).index(True)

                data[0] = [ _ for _ in data[0][_min_x:] ]
                data[1] = [ _ for _ in data[1][_min_x:] ]

            max_x = _get_max(data[0], max_x)
            min_x = _get_min(data[0], min_x)

            max_y = _get_max(data[1], max_y)
            min_y = _get_min(data[1], min_y)

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                ),
                alpha = plot_alpha
            )

        #~ will use this foreach to fill the polyfit for GR
        polyfit_res = []
        xp = []
        yp = []

        #~ plot all GR data and gather all x and y for evaluating the polyfit
        for label, data, eos in zip(
            all_label_GR,
            all_data_GR,
            [ _ for _ in set( [ _["name"] for _ in severalEOSs ] ) ]
        ):

            #~ we set mimimal compactenss threshold and cut out all entries
            #~ who are below it only if we are interested in stable solutions
            if append_stable:
                _min_x = list(
                    map(lambda _: _ >= min_compactness, data[0])
                ).index(True)

                data[0] = [ _ for _ in data[0][_min_x:] ]
                data[1] = [ _ for _ in data[1][_min_x:] ]

            max_x = _get_max(data[0], max_x)
            min_x = _get_min(data[0], min_x)

            max_y = _get_max(data[1], max_y)
            min_y = _get_min(data[1], min_y)

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                marker = markers.get(eos, None),
                color = GR_color_markers,
                markerfacecolor = GR_color_markers,
                markeredgecolor = GR_color_markers,
                alpha = plot_alpha
            )

            xp.append(data[0])
            yp.append(data[1])

        coef, chi_red, p = _get_polyfit_res(xp, yp)

        max_y_down = 0
        min_y_down = 1e9

        #~ average over all EOSs of all the residuals
        delta_all = 0
        n_all = 0

        #~ average over all EOSs of largest residual
        delta_all_max = 0
        n_all_max = 0

        #~ the largest residual across all EOSs
        delta_max = 0

        #~ for the generated polyfit function calcualte the
        #~ relative error and plot it donw for GR
        for label, data, eos in zip(
            all_label_GR,
            all_data_GR,
            [ _ for _ in set( [ _["name"] for _ in severalEOSs ] ) ]
        ):

            _data = [
                abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
            ]

            delta_all += sum(_data)
            n_all += len(_data)

            delta_all_max += max(_data)
            n_all_max += 1

            delta_max = _get_max(_data, delta_max)

            max_y_down = _get_max(_data, max_y_down)
            min_y_down = _get_min(_data, min_y_down)

            ax_down.plot(
                data[0],
                _data,
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], _data),
                marker = markers.get(eos, None),
                color = GR_color_markers,
                markerfacecolor = GR_color_markers,
                markeredgecolor = GR_color_markers,
                alpha = plot_alpha
            )

        avg_L_1 = delta_all/n_all
        avg_L_inf = delta_all_max/n_all_max
        L_inf_worst = delta_max

        print(
            "\n GR fit"
            "\n\t $\chi_r^2$ = {:.3e}"
            "\n\t $a_0$ = {:.3e}"
            "\n\t $a_1$ = {:.3e}"
            "\n\t $a_4$ = {:.3e}"
            "\n\t $< L_1 >$ = {:.3e}"
            "\n\t $< L_\inf >$ = {:.3e}"
            "\n\t $ L_\inf $ = {:.3e}\n".format(
                chi_red,
                coef[0],
                coef[1],
                coef[4],
                avg_L_1,
                avg_L_inf,
                L_inf_worst
            )
        )

        lines_polyfit = [
            Line2D(
                [0], [0],
                color = GR_color_fit,
                marker = None,
                linestyle = "-",
                linewidth = 1.5,
                label = "GR fit"
            )
        ]

        p_x = np.linspace(min_x, max_x, 100)
        p_y = [ p(_) for _ in p_x ]

        #~ generate 100 points between min and max of x and plot it values
        ax_up.plot(
            p_x,
            p_y,
            label = None,
            linewidth = 2.5,
            linestyle = "-",
            markersize = 0,
            markevery = 0,
            marker = None,
            color = GR_color_fit,
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*(1 + avg_L_inf),
            np.array(p_y)*(1 - avg_L_inf),
            facecolor=GR_color_markers,
            alpha= plot_alpha - 0.2
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*( 1 + L_inf_worst ),
            np.array(p_y)*( 1 - L_inf_worst ),
            facecolor=GR_color_markers,
            alpha= plot_alpha - 0.4
        )

        #~ now do the same for each color if there are more than 1
        for c, l in itertools.product(colors.items(), linestyles.items()):

            #~ the colors will have label key containing the name of parameter
            #~ which they represetn
            if c[0] == "label" or l[0] == "label":
                continue

            xp = []
            yp = []

            for data, eos in zip(all_data, severalEOSs):
                #~ if the current eos has parameter value equal to the current one
                #~ lets append its data
                if eos[colors["label"]] == c[0] and eos[linestyles["label"]] == l[0]:
                    xp.append( data[0] )
                    yp.append( data[1] )

            #~ colors and linestyles have all possible combinations of EOSs
            #~ but we may not need all of them but only porsion determineed
            #~ by severalEOSs, so a quick fix, if no points added we just continue
            if not xp:
                continue

            #~ expand all the data into flat list to calculate the polyfit
            coef, chi_red, p = _get_polyfit_res(xp, yp)

            max_y_down = 0
            min_y_down = 1e9

            #~ average over all EOSs of all the residuals
            delta_all = 0
            n_all = 0

            #~ average over all EOSs of largest residual
            delta_all_max = 0
            n_all_max = 0

            #~ the largest residual across all EOSs
            delta_max = 0

            for data, eos in zip(all_data, severalEOSs):
                #~ if the current eos has parameter value equal to the current one
                #~ lets append its data
                if eos[colors["label"]] == c[0] and eos[linestyles["label"]] == l[0]:

                    _data = [
                        abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
                    ]

                    delta_all += sum(_data)
                    n_all += len(_data)

                    delta_all_max += max(_data)
                    n_all_max += 1

                    delta_max = _get_max(_data, delta_max)

                    max_y_down = _get_max(_data, max_y_down)
                    min_y_down = _get_min(_data, min_y_down)

                    ax_down.plot(
                        data[0],
                        _data,
                        label = None,
                        linewidth = 0,
                        markersize = 5.5,
                        markevery = self._get_markevry(data[0], _data),
                        marker = markers.get(eos["name"], None),
                        color = c[1],
                        markerfacecolor = c[1],
                        markeredgecolor = c[1],
                        alpha = plot_alpha - 0.1
                    )

            avg_L_1 = delta_all/n_all
            avg_L_inf = delta_all_max/n_all_max
            L_inf_worst = delta_max

            print(
                "\n lambda = {:.3e}, m = {:.3e} fit"
                "\n\t $\chi_r^2$ = {:.3e}"
                "\n\t $a_0$ = {:.3e}"
                "\n\t $a_1$ = {:.3e}"
                "\n\t $a_4$ = {:.3e}"
                "\n\t $< L_1 >$ = {:.3e}"
                "\n\t $< L_\inf >$ = {:.3e}"
                "\n\t $ L_\inf $ = {:.3e}\n".format(
                    l[0] if linestyles["label"] == "lambda" else c[0],
                    c[0] if colors["label"] == "m" else l[0],
                    chi_red,
                    coef[0],
                    coef[1],
                    coef[4],
                    avg_L_1,
                    avg_L_inf,
                    L_inf_worst
                )
            )

            #~ lines_polyfit.append(
                #~ Line2D(
                    #~ [0], [0],
                    #~ color = c[1],
                    #~ marker = None,
                    #~ linestyle = l[1],
                    #~ linewidth = 1.5,
                    #~ label = "$\\lambda$ = {:.3e},\n m = {:.3e} fit".format(
                        #~ l[0] if linestyles["label"] == "lambda" else c[0],
                        #~ c[0] if colors["label"] == "m" else l[0]
                    #~ )
                #~ )
            #~ )

            p_x = np.linspace(min_x, max_x, 100)
            p_y = [ p(_) for _ in p_x ]

            ax_up.plot(
                p_x,
                p_y,
                label = None,
                linewidth = 2.5,
                linestyle = l[1],
                markersize = 0,
                markevery = 0,
                marker = None,
                color = self._luminosity_color(c[1], 1.1)
            )

            #~ ax_up.fill_between(
                #~ p_x,
                #~ np.array(p_y)*(1 + avg_L_inf),
                #~ np.array(p_y)*(1 - avg_L_inf),
                #~ facecolor=c[1],
                #~ alpha= plot_alpha - 0.25
            #~ )

            #~ ax_up.fill_between(
                #~ p_x,
                #~ np.array(p_y)*( 1 + L_inf_worst ),
                #~ np.array(p_y)*( 1 - L_inf_worst ),
                #~ facecolor=c[1],
                #~ alpha= plot_alpha - 0.5
            #~ )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles, severalEOSs
        )

        ax_up.add_artist( ax_up.legend(
            handles = [
                *lines_markers
            ],
            loc="upper left",
            fontsize=8,
            handlelength=3,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol = 3,
            frameon = False,
            mode = None
        ) )

        ax_up.add_artist( ax_up.legend(
            handles = [
                *lines_colors, *lines_polyfit, *lines_linestyles
            ],
            loc="lower right",
            fontsize=8,
            handlelength=3,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol = 2,
            frameon = False,
            mode = None
        ) )

        ax_up.set_xlim(min_x, max_x)
        ax_up.set_ylim(min_y, max_y)

        ax_down.set_ylim(1e-3, 1.5e0)

        plt.savefig(
            'uniTilde.eps', format="eps",
            bbox_inches='tight',
            dpi=1200,
            pad_inches=0
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

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            ax.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        ax.legend(
            handles = [*lines_markers, *lines_colors, *lines_linestyles],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1.25,
            ncol = 3,
            frameon = False,
            mode = None
        )

        ax_up.set_xlim(0.09)

        plt.show()

        return

    def plot_severalEOSs_uniBarI_polyFit(
        self, severalEOSs, append_stable = "stable"
    ):
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

        _get_max = lambda data, _: max(data) if max(data) > _ else _
        _get_min = lambda data, _: min(data) if min(data) < _ else _

        all_label, all_headline, all_data = self.get_severalEOS_uniBarI_data(
            severalEOSs, append = append_stable
        )

        fig, all_axes  = self._get_figure(
            2,  1,  self._3by1_shareX_grid_placement, height_ratios = [2,1]
        )

        ax_up = all_axes[0]
        ax_down = all_axes[1]

        self._set_parms(ax_up, "", "$I/(M^3)$")

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        max_x = max_y = 0
        min_x = min_y = 1e9
        min_compactness = 0.09

        color_fit = "#fcc006"
        color_avg = "#9f2305"
        color_avg_worst = "#dfc5fe"

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ we set mimimal compactenss threshold and cut out all entries
            #~ who are below it only if we are interested only in stable solutions
            if append_stable:
                _min_x = list(
                    map(lambda _: _ >= min_compactness, data[0])
                ).index(True)

                data[0] = [ _ for _ in data[0][_min_x:] ]
                data[1] = [ _ for _ in data[1][_min_x:] ]

            max_x = _get_max(data[0], max_x)
            min_x = _get_min(data[0], min_x)

            max_y = _get_max(data[1], max_y)
            min_y = _get_min(data[1], min_y)

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 1.5,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )

        coef, rest = polyfit(
            x = [ __**-1 for _ in all_data for __ in _[0] ],
            y = [ __ for _ in all_data for __ in _[1] ],
            deg = [ 1, 2, 3, 4 ],
            #~ w = np.sqrt(np.array([ __ for _ in all_data for __ in _[1] ])),
            full = True
        )
        chi_red = rest[0][0]/(len([ __ for _ in all_data for __ in _[0] ]) - 4)

        p = lambda x: \
        coef[1]*(1/x)**1 + coef[2]*(1/x)**2 + coef[3]*(1/x)**3 + coef[4]*(1/x)**4

        p_x = np.linspace(min_x, max_x, 100)
        p_y = [ p(_) for _ in p_x ]

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        lines_polyfit = [
            Line2D(
                [0], [0],
                color = color_fit,
                marker = None,
                linewidth = 2,
                linestyle = "-",
                label = "poly fit, $\chi_r^2$ = {:.3e}"
                    "\n {:.3f}$(1/x)^1$ + {:.3f}$(1/x)^2$ + {:.3f}$(1/x)^3$ + {:.3f}$(1/x)^4$".format(
                        chi_red,
                        coef[1],
                        coef[2],
                        coef[3],
                        coef[4]
                )
            )
        ]

        ax_up.get_shared_x_axes().join(ax_up, ax_down)
        ax_up.set_xticklabels([])

        ax_down.set_yscale("log")
        self._set_parms(ax_down, "M/R", r"$\left| 1 - \bar I/\bar I_{fit} \right|$")

        max_y_down = 0
        min_y_down = 1e9

        #~ average over all EOSs of all the residuals
        delta_all = 0
        n_all = 0

        #~ average over all EOSs of largest residual
        delta_all_max = 0
        n_all_max = 0

        #~ the largest residual across all EOSs
        delta_max = 0

        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            _data = [
                abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
            ]

            delta_all += sum(_data)
            n_all += len(_data)

            delta_all_max += max(_data)
            n_all_max += 1

            delta_max = _get_max(_data, delta_max)

            max_y_down = _get_max(_data, max_y_down)
            min_y_down = _get_min(_data, min_y_down)

            ax_down.plot(
                data[0],
                _data,
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], _data),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                )
            )
        avg_L_1 = delta_all/n_all
        avg_L_inf = delta_all_max/n_all_max
        L_inf_worst = delta_max

        ax_up.plot(
            p_x,
            p_y,
            color = color_fit,
            marker = None,
            linewidth = 2,
            linestyle = "-",
            label = None
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*( 1 + L_inf_worst ),
            np.array(p_y)*( 1 - L_inf_worst ),
            facecolor=color_avg_worst,
            alpha=0.5
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*(1 + avg_L_inf),
            np.array(p_y)*(1 - avg_L_inf),
            facecolor=color_avg,
            alpha=0.5
        )

        ax_up.legend(
            handles = [
                *lines_markers, *lines_colors, *lines_linestyles, *lines_polyfit
            ],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1.25,
            ncol = 3,
            frameon = False,
            mode = None
        )

        ax_up.set_xlim(min_x, max_x)
        ax_up.set_ylim(min_y, max_y)

        ax_down.set_ylim(1e-3, 1e0)

        print(
            "\n All fit information"
            "\n\t $\chi_r^2$ = {:.3e}"
            "\n\t $a_1$ = {:.3e}"
            "\n\t $a_2$ = {:.3e}"
            "\n\t $a_3$ = {:.3e}"
            "\n\t $a_4$ = {:.3e}"
            "\n\t $< L_1 >$ = {:.3e}"
            "\n\t $< L_\inf >$ = {:.3e}"
            "\n\t $ L_\inf $ = {:.3e}\n".format(
                chi_red,
                coef[1],
                coef[2],
                coef[3],
                coef[4],
                avg_L_1,
                avg_L_inf,
                L_inf_worst
            )
        )

        if n_all_max != len(severalEOSs):
            print(
                "\n SOMETHING FISSHY, n_all = {}, len EOSs = {} \n".format(
                n_all_max, len(severalEOSs))
            )

        plt.show()

        return

    def plot_severalEOSs_uniBarI_polyFitAll_GR(
        self, severalEOSs, append_stable = "stable"
    ):
        """
        plot severalEOS unifersal Tilde I relationships
        <severalEOSs> with dictionaries see get_severalEOS_data for the format

        EXAMPLE INPUT

        append_stable = "stable" is set as default, so only stable solutions
        will be calculated

        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        #~ expand the data and do the polyfit
        def _get_polyfit_res(xp, yp):
            xp = [ __**-1 for _ in xp for __ in _ ]
            yp = [ __ for _ in yp for __ in _ ]

            coef, rest = polyfit(
                x = xp, y = yp,
                deg = [ 1, 2, 3, 4 ],
                #~ w = np.sqrt(yp),
                full = True
            )

            #~ calcualte the chi reduce
            chi_red = rest[0][0]/(len(xp) - 4)
            p = lambda x: coef[1]*x**-1 + coef[2]*x**-2 + coef[3]*x**-3 + coef[4]*x**-4

            return coef, chi_red, p

        _get_max = lambda data, _: max(data) if max(data) > _ else _
        _get_min = lambda data, _: min(data) if min(data) < _ else _

        all_label, all_headline, all_data = self.get_severalEOS_uniBarI_data(
            severalEOSs, append = append_stable
        )

        all_label_GR, all_headline_GR, all_data_GR = self.get_severalEOS_uniBarI_data( [
                { "name": _, "beta": 0, "m": 0, "lambda": 0 }
                for _ in set( [ _["name"] for _ in severalEOSs ] )
            ],
            append = append_stable
        )

        fig, all_axes  = self._get_figure(
            2,  1,  self._3by1_shareX_grid_placement, height_ratios = [2,1]
        )

        ax_up = all_axes[0]
        ax_down = all_axes[1]

        ax_up.get_shared_x_axes().join(ax_up, ax_down)
        ax_up.set_xticklabels([])

        ax_down.set_yscale("log")
        self._set_parms(ax_down, "M/R", r"$\left| 1 - \bar I/\bar I_{fit} \right|$  ")

        self._set_parms(ax_up, "", r"$I/(M^3)$")

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        max_x = max_y = 0
        min_x = min_y = 1e9
        min_compactness = 0.09

        #~ GR_color_markers = "#ef4026"
        GR_color_markers = "#c0022f"
        #~ GR_color_markers = "#a9f971"
        #~ GR_color_fit = "#ed0dd9"
        GR_color_fit = GR_color_markers

        plot_alpha = 0.75

        #~ lets plot severEOSs on the up plot and eventually cut out data
        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ we set mimimal compactenss threshold and cut out all entries
            #~ who are below it only if we are interested in stable solutions
            if append_stable:
                _min_x = list(
                    map(lambda _: _ >= min_compactness, data[0])
                ).index(True)

                data[0] = [ _ for _ in data[0][_min_x:] ]
                data[1] = [ _ for _ in data[1][_min_x:] ]

            max_x = _get_max(data[0], max_x)
            min_x = _get_min(data[0], min_x)

            max_y = _get_max(data[1], max_y)
            min_y = _get_min(data[1], min_y)

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                ),
                alpha = plot_alpha
            )

        #~ will use this foreach to fill the polyfit for GR
        polyfit_res = []
        xp = []
        yp = []

        #~ plot all GR data and gather all x and y for evaluating the polyfit
        for label, data, eos in zip(
            all_label_GR,
            all_data_GR,
            [ _ for _ in set( [ _["name"] for _ in severalEOSs ] ) ]
        ):

            #~ we set mimimal compactenss threshold and cut out all entries
            #~ who are below it only if we are interested in stable solutions
            if append_stable:
                _min_x = list(
                    map(lambda _: _ >= min_compactness, data[0])
                ).index(True)

                data[0] = [ _ for _ in data[0][_min_x:] ]
                data[1] = [ _ for _ in data[1][_min_x:] ]

            max_x = _get_max(data[0], max_x)
            min_x = _get_min(data[0], min_x)

            max_y = _get_max(data[1], max_y)
            min_y = _get_min(data[1], min_y)

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                marker = markers.get(eos, None),
                color = GR_color_markers,
                markerfacecolor = GR_color_markers,
                markeredgecolor = GR_color_markers,
                alpha = plot_alpha
            )

            xp.append(data[0])
            yp.append(data[1])

        coef, chi_red, p = _get_polyfit_res(xp, yp)

        max_y_down = 0
        min_y_down = 1e9

        #~ average over all EOSs of all the residuals
        delta_all = 0
        n_all = 0

        #~ average over all EOSs of largest residual
        delta_all_max = 0
        n_all_max = 0

        #~ the largest residual across all EOSs
        delta_max = 0

        #~ for the generated polyfit function calcualte the
        #~ relative error and plot it donw for GR
        for label, data, eos in zip(
            all_label_GR,
            all_data_GR,
            [ _ for _ in set( [ _["name"] for _ in severalEOSs ] ) ]
        ):

            _data = [
                abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
            ]

            delta_all += sum(_data)
            n_all += len(_data)

            delta_all_max += max(_data)
            n_all_max += 1

            delta_max = _get_max(_data, delta_max)

            max_y_down = _get_max(_data, max_y_down)
            min_y_down = _get_min(_data, min_y_down)

            ax_down.plot(
                data[0],
                _data,
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], _data),
                marker = markers.get(eos, None),
                color = GR_color_markers,
                markerfacecolor = GR_color_markers,
                markeredgecolor = GR_color_markers
            )

        avg_L_1 = delta_all/n_all
        avg_L_inf = delta_all_max/n_all_max
        L_inf_worst = delta_max

        print(
            "\n GR fit"
            "\n\t $\chi_r^2$ = {:.3e}"
            "\n\t $a_1$ = {:.3e}"
            "\n\t $a_2$ = {:.3e}"
            "\n\t $a_3$ = {:.3e}"
            "\n\t $a_4$ = {:.3e}"
            "\n\t $< L_1 >$ = {:.3e}"
            "\n\t $< L_\inf >$ = {:.3e}"
            "\n\t $ L_\inf $ = {:.3e}\n".format(
                chi_red,
                coef[1],
                coef[2],
                coef[3],
                coef[4],
                avg_L_1,
                avg_L_inf,
                L_inf_worst
            )
        )

        #~ add the polyfit lien as entry in polyfit lines to be displayed in legend
        lines_polyfit = [
            Line2D(
                [0], [0],
                color = GR_color_fit,
                marker = None,
                linestyle = "-",
                linewidth = 1.5,
                label = "GR fit"
            )
        ]

        p_x = np.linspace(min_x, max_x, 100)
        p_y = [ p(_) for _ in p_x ]

        #~ generate 100 points between min and max of x and plot it values
        ax_up.plot(
            p_x,
            p_y,
            label = None,
            linewidth = 2,
            linestyle = "-",
            markersize = 0,
            markevery = 0,
            marker = None,
            color = GR_color_fit,
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*(1 + avg_L_inf),
            np.array(p_y)*(1 - avg_L_inf),
            facecolor=GR_color_markers,
            alpha= plot_alpha - 0.25
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*( 1 + L_inf_worst ),
            np.array(p_y)*( 1 - L_inf_worst ),
            facecolor=GR_color_markers,
            alpha= plot_alpha - 0.5
        )

        #~ now do the same for each color if there are more than 1
        for k, v in colors.items():

            #~ the colors will have label key containing the name of parameter
            #~ which they represetn
            if k == "label":
                continue

            xp = []
            yp = []

            for data, eos in zip(all_data, severalEOSs):
                #~ if the current eos has parameter value equal to the current one
                #~ lets append its data
                if eos[colors["label"]] == k:
                    xp.append( data[0] )
                    yp.append( data[1] )

            #~ expand all the data into flat list to calculate the polyfit
            coef, chi_red, p = _get_polyfit_res(xp, yp)

            max_y_down = 0
            min_y_down = 1e9

            #~ average over all EOSs of all the residuals
            delta_all = 0
            n_all = 0

            #~ average over all EOSs of largest residual
            delta_all_max = 0
            n_all_max = 0

            #~ the largest residual across all EOSs
            delta_max = 0

            for data, eos in zip(all_data, severalEOSs):
                #~ if the current eos has parameter value equal to the current one
                #~ lets append its data
                if eos[colors["label"]] == k:

                    _data = [
                        abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
                    ]

                    delta_all += sum(_data)
                    n_all += len(_data)

                    delta_all_max += max(_data)
                    n_all_max += 1

                    delta_max = _get_max(_data, delta_max)

                    max_y_down = _get_max(_data, max_y_down)
                    min_y_down = _get_min(_data, min_y_down)

                    ax_down.plot(
                        data[0],
                        _data,
                        label = None,
                        linewidth = 0,
                        markersize = 5.5,
                        markevery = self._get_markevry(data[0], _data),
                        marker = markers.get(eos["name"], None),
                        color = v,
                        markerfacecolor = v,
                        markeredgecolor = v
                    )

            avg_L_1 = delta_all/n_all
            avg_L_inf = delta_all_max/n_all_max
            L_inf_worst = delta_max

            print(
                "\n {} {:.3e} fit"
                "\n\t $\chi_r^2$ = {:.3e}"
                "\n\t $a_1$ = {:.3e}"
                "\n\t $a_2$ = {:.3e}"
                "\n\t $a_3$ = {:.3e}"
                "\n\t $a_4$ = {:.3e}"
                "\n\t $< L_1 >$ = {:.3e}"
                "\n\t $< L_\inf >$ = {:.3e}"
                "\n\t $ L_\inf $ = {:.3e}\n".format(
                    "$\\lambda =$ " if colors["label"] == "lambda" else "m =",
                    k,
                    chi_red,
                    coef[1],
                    coef[2],
                    coef[3],
                    coef[4],
                    avg_L_1,
                    avg_L_inf,
                    L_inf_worst
                )
            )
            lines_polyfit.append(
                Line2D(
                    [0], [0],
                    color = v,
                    marker = None,
                    linestyle = "-",
                    linewidth = 1.5,
                    label = "{} {:.3e} fit".format(
                        "$\\lambda =$ " if colors["label"] == "lambda" else "m =",
                        k
                    )
                )
            )

            p_x = np.linspace(min_x, max_x, 100)
            p_y = [ p(_) for _ in p_x ]

            ax_up.plot(
                p_x,
                p_y,
                label = None,
                linewidth = 2,
                linestyle = "-",
                markersize = 0,
                markevery = 0,
                marker = None,
                color = v,
            )

            ax_up.fill_between(
                p_x,
                np.array(p_y)*(1 + avg_L_inf),
                np.array(p_y)*(1 - avg_L_inf),
                facecolor=v,
                alpha= plot_alpha - 0.25
            )

            ax_up.fill_between(
                p_x,
                np.array(p_y)*( 1 + L_inf_worst ),
                np.array(p_y)*( 1 - L_inf_worst ),
                facecolor=v,
                alpha= plot_alpha - 0.5
            )

        #~ THERE IS ADDITIONAL CODE FOR LINSTYLES WHICH IS NOT TRANSFERED HERE

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles
        )

        ax_up.legend(
            handles = [
                *lines_markers, *lines_colors, *lines_polyfit
            ],
            loc="best",
            fontsize=8,
            handlelength=2,
            numpoints=1,
            fancybox=True,
            markerscale = 1.25,
            ncol = 4,
            frameon = False,
            mode = None
        )

        ax_up.set_xlim(min_x, max_x)
        ax_up.set_ylim(min_y, max_y)

        ax_down.set_ylim(1e-3, 1e0)

        plt.show()

        return

    def plot_severalEOSs_uniBarI_polyFitAll_GR_ParmsProduct(
        self, severalEOSs, append_stable = "stable"
    ):
        """
        plot severalEOS unifersal Tilde I relationships
        <severalEOSs> with dictionaries see get_severalEOS_data for the format

        EXAMPLE INPUT

        append_stable = "stable" is set as default, so only stable solutions
        will be calculated

        severalEOSs = [
            { "name": "SLy4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "APR4", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "FPS", "beta": 0, "m": 0, "lambda": 0 },
            { "name": "WFF2", "beta": 0, "m": 0, "lambda": 0 }
        ]
        """

        #~ expand the data and do the polyfit
        def _get_polyfit_res(xp, yp):
            xp = [ __**-1 for _ in xp for __ in _ ]
            yp = [ __ for _ in yp for __ in _ ]

            coef, rest = polyfit(
                x = xp, y = yp,
                deg = [ 1, 2, 3, 4 ],
                #~ w = np.sqrt(yp),
                full = True
            )

            #~ calcualte the chi reduce
            chi_red = rest[0][0]/(len(xp) - 4)
            p = lambda x: coef[1]*x**-1 + coef[2]*x**-2 + coef[3]*x**-3 + coef[4]*x**-4

            return coef, chi_red, p

        _get_max = lambda data, _: max(data) if max(data) > _ else _
        _get_min = lambda data, _: min(data) if min(data) < _ else _

        all_label, all_headline, all_data = self.get_severalEOS_uniBarI_data(
            severalEOSs, append = append_stable
        )

        all_label_GR, all_headline_GR, all_data_GR = self.get_severalEOS_uniBarI_data( [
                { "name": _, "beta": 0, "m": 0, "lambda": 0 }
                for _ in set( [ _["name"] for _ in severalEOSs ] )
            ],
            append = append_stable
        )

        fig, all_axes  = self._get_figure(
            2,  1,  self._3by1_shareX_grid_placement, height_ratios = [2,1]
        )

        #~ turn off the tight layout since it does not permit sticking the two axes
        #~ but can do it manually
        fig.set_tight_layout(False)
        fig.set_rasterized(True)
        fig.subplots_adjust(wspace=0, hspace=0)

        ax_up = all_axes[0]
        ax_down = all_axes[1]

        ax_up.get_shared_x_axes().join(ax_up, ax_down)
        ax_up.set_xticklabels([])

        ax_down.set_yscale("log")
        self._set_parms(ax_down, "M/R", r"$\left| 1 - \bar I/\bar I_{fit} \right|$  ")

        self._set_parms(ax_up, "", r"$I/(M^3)$")

        markers, colors, linestyles = self._get_MSs_Cs_LSs(severalEOSs)

        max_x = max_y = 0
        min_x = min_y = 1e9
        min_compactness = 0.09

        #~ GR_color_markers = "#ef4026"
        GR_color_markers = "#469990"
        #~ GR_color_markers = "#a9f971"
        #~ GR_color_fit = "#ed0dd9"
        GR_color_fit = self._luminosity_color(GR_color_markers, 1.1)

        plot_alpha = 0.6

        #~ lets plot severEOSs on the up plot and eventually cut out data
        for label, data, eos in zip( all_label, all_data, severalEOSs ):

            #~ we set mimimal compactenss threshold and cut out all entries
            #~ who are below it only if we are interested in stable solutions
            if append_stable:
                _min_x = list(
                    map(lambda _: _ >= min_compactness, data[0])
                ).index(True)

                data[0] = [ _ for _ in data[0][_min_x:] ]
                data[1] = [ _ for _ in data[1][_min_x:] ]

            max_x = _get_max(data[0], max_x)
            min_x = _get_min(data[0], min_x)

            max_y = _get_max(data[1], max_y)
            min_y = _get_min(data[1], min_y)

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                **self._get_plot_keywords(
                    markers, colors, linestyles,
                    {
                        "name": eos["name"],
                        "m": eos["m"],
                        "lambda": eos["lambda"]
                    }
                ),
                alpha = plot_alpha
            )

        #~ will use this foreach to fill the polyfit for GR
        polyfit_res = []
        xp = []
        yp = []

        #~ plot all GR data and gather all x and y for evaluating the polyfit
        for label, data, eos in zip(
            all_label_GR,
            all_data_GR,
            [ _ for _ in set( [ _["name"] for _ in severalEOSs ] ) ]
        ):

            #~ we set mimimal compactenss threshold and cut out all entries
            #~ who are below it only if we are interested in stable solutions
            if append_stable:
                _min_x = list(
                    map(lambda _: _ >= min_compactness, data[0])
                ).index(True)

                data[0] = [ _ for _ in data[0][_min_x:] ]
                data[1] = [ _ for _ in data[1][_min_x:] ]

            max_x = _get_max(data[0], max_x)
            min_x = _get_min(data[0], min_x)

            max_y = _get_max(data[1], max_y)
            min_y = _get_min(data[1], min_y)

            ax_up.plot(
                data[0],
                data[1],
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], data[1]),
                marker = markers.get(eos, None),
                color = GR_color_markers,
                markerfacecolor = GR_color_markers,
                markeredgecolor = GR_color_markers,
                alpha = plot_alpha
            )

            xp.append(data[0])
            yp.append(data[1])

        coef, chi_red, p = _get_polyfit_res(xp, yp)

        max_y_down = 0
        min_y_down = 1e9

        #~ average over all EOSs of all the residuals
        delta_all = 0
        n_all = 0

        #~ average over all EOSs of largest residual
        delta_all_max = 0
        n_all_max = 0

        #~ the largest residual across all EOSs
        delta_max = 0

        #~ for the generated polyfit function calcualte the
        #~ relative error and plot it donw for GR
        for label, data, eos in zip(
            all_label_GR,
            all_data_GR,
            [ _ for _ in set( [ _["name"] for _ in severalEOSs ] ) ]
        ):

            _data = [
                abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
            ]

            delta_all += sum(_data)
            n_all += len(_data)

            delta_all_max += max(_data)
            n_all_max += 1

            delta_max = _get_max(_data, delta_max)

            max_y_down = _get_max(_data, max_y_down)
            min_y_down = _get_min(_data, min_y_down)

            ax_down.plot(
                data[0],
                _data,
                label = None,
                linewidth = 0,
                markersize = 5.5,
                markevery = self._get_markevry(data[0], _data),
                marker = markers.get(eos, None),
                color = GR_color_markers,
                markerfacecolor = GR_color_markers,
                markeredgecolor = GR_color_markers,
                alpha = plot_alpha
            )

        avg_L_1 = delta_all/n_all
        avg_L_inf = delta_all_max/n_all_max
        L_inf_worst = delta_max

        print(
            "\n GR fit"
            "\n\t $\chi_r^2$ = {:.3e}"
            "\n\t $a_1$ = {:.3e}"
            "\n\t $a_2$ = {:.3e}"
            "\n\t $a_3$ = {:.3e}"
            "\n\t $a_4$ = {:.3e}"
            "\n\t $< L_1 >$ = {:.3e}"
            "\n\t $< L_\inf >$ = {:.3e}"
            "\n\t $ L_\inf $ = {:.3e}\n".format(
                chi_red,
                coef[1],
                coef[2],
                coef[3],
                coef[4],
                avg_L_1,
                avg_L_inf,
                L_inf_worst
            )
        )

        lines_polyfit = [
            Line2D(
                [0], [0],
                color = GR_color_fit,
                marker = None,
                linestyle = "-",
                linewidth = 1.5,
                label = "GR fit"
            )
        ]

        p_x = np.linspace(min_x, max_x, 100)
        p_y = [ p(_) for _ in p_x ]

        #~ generate 100 points between min and max of x and plot it values
        ax_up.plot(
            p_x,
            p_y,
            label = None,
            linewidth = 2.5,
            linestyle = "-",
            markersize = 0,
            markevery = 0,
            marker = None,
            color = GR_color_fit,
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*(1 + avg_L_inf),
            np.array(p_y)*(1 - avg_L_inf),
            facecolor=GR_color_markers,
            alpha= plot_alpha - 0.2
        )

        ax_up.fill_between(
            p_x,
            np.array(p_y)*( 1 + L_inf_worst ),
            np.array(p_y)*( 1 - L_inf_worst ),
            facecolor=GR_color_markers,
            alpha= plot_alpha - 0.4
        )

        #~ now do the same for each color if there are more than 1
        for c, l in itertools.product(colors.items(), linestyles.items()):

            #~ the colors will have label key containing the name of parameter
            #~ which they represetn
            if c[0] == "label" or l[0] == "label":
                continue

            xp = []
            yp = []

            for data, eos in zip(all_data, severalEOSs):
                #~ if the current eos has parameter value equal to the current one
                #~ lets append its data
                if eos[colors["label"]] == c[0] and eos[linestyles["label"]] == l[0]:
                    xp.append( data[0] )
                    yp.append( data[1] )

            #~ colors and linestyles have all possible combinations of EOSs
            #~ but we may not need all of them but only porsion determineed
            #~ by severalEOSs, so a quick fix, if no points added we just continue
            if not xp:
                continue

            #~ expand all the data into flat list to calculate the polyfit
            coef, chi_red, p = _get_polyfit_res(xp, yp)

            max_y_down = 0
            min_y_down = 1e9

            #~ average over all EOSs of all the residuals
            delta_all = 0
            n_all = 0

            #~ average over all EOSs of largest residual
            delta_all_max = 0
            n_all_max = 0

            #~ the largest residual across all EOSs
            delta_max = 0

            for data, eos in zip(all_data, severalEOSs):
                #~ if the current eos has parameter value equal to the current one
                #~ lets append its data
                if eos[colors["label"]] == c[0] and eos[linestyles["label"]] == l[0]:

                    _data = [
                        abs(1 - _/p(__)) for _,__ in zip(data[1], data[0])
                    ]

                    delta_all += sum(_data)
                    n_all += len(_data)

                    delta_all_max += max(_data)
                    n_all_max += 1

                    delta_max = _get_max(_data, delta_max)

                    max_y_down = _get_max(_data, max_y_down)
                    min_y_down = _get_min(_data, min_y_down)

                    ax_down.plot(
                        data[0],
                        _data,
                        label = None,
                        linewidth = 0,
                        markersize = 5.5,
                        markevery = self._get_markevry(data[0], _data),
                        marker = markers.get(eos["name"], None),
                        color = c[1],
                        markerfacecolor = c[1],
                        markeredgecolor = c[1],
                        alpha = plot_alpha - 0.1
                    )

            avg_L_1 = delta_all/n_all
            avg_L_inf = delta_all_max/n_all_max
            L_inf_worst = delta_max

            print(
                "\n lambda = {:.3e}, m = {:.3e} fit"
                "\n\t $\chi_r^2$ = {:.3e}"
                "\n\t $a_1$ = {:.3e}"
                "\n\t $a_2$ = {:.3e}"
                "\n\t $a_3$ = {:.3e}"
                "\n\t $a_4$ = {:.3e}"
                "\n\t $< L_1 >$ = {:.3e}"
                "\n\t $< L_\inf >$ = {:.3e}"
                "\n\t $ L_\inf $ = {:.3e}\n".format(
                    l[0] if linestyles["label"] == "lambda" else c[0],
                    c[0] if colors["label"] == "m" else l[0],
                    chi_red,
                    coef[1],
                    coef[2],
                    coef[3],
                    coef[4],
                    avg_L_1,
                    avg_L_inf,
                    L_inf_worst
                )
            )

            #~ lines_polyfit.append(
                #~ Line2D(
                    #~ [0], [0],
                    #~ color = c[1],
                    #~ marker = None,
                    #~ linestyle = l[1],
                    #~ linewidth = 1.5,
                    #~ label = "$\\lambda$ = {:.3e},\n m = {:.3e} fit".format(
                        #~ l[0] if linestyles["label"] == "lambda" else c[0],
                        #~ c[0] if colors["label"] == "m" else l[0]
                    #~ )
                #~ )
            #~ )

            p_x = np.linspace(min_x, max_x, 100)
            p_y = [ p(_) for _ in p_x ]

            ax_up.plot(
                p_x,
                p_y,
                label = None,
                linewidth = 2.5,
                linestyle = l[1],
                markersize = 0,
                markevery = 0,
                marker = None,
                color = self._luminosity_color(c[1], 1.1)
            )

            #~ ax_up.fill_between(
                #~ p_x,
                #~ np.array(p_y)*(1 + avg_L_inf),
                #~ np.array(p_y)*(1 - avg_L_inf),
                #~ facecolor=c[1],
                #~ alpha= plot_alpha - 0.25
            #~ )

            #~ ax_up.fill_between(
                #~ p_x,
                #~ np.array(p_y)*( 1 + L_inf_worst ),
                #~ np.array(p_y)*( 1 - L_inf_worst ),
                #~ facecolor=c[1],
                #~ alpha= plot_alpha - 0.5
            #~ )

        lines_markers, lines_colors, lines_linestyles = self._get_lines_MSs_Cs_LSs(
            markers, colors, linestyles, severalEOSs
        )

        ax_up.add_artist( ax_up.legend(
            handles = [
                *lines_markers
            ],
            loc="upper right",
            fontsize=8,
            handlelength=3,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol = 3,
            frameon = False,
            mode = None
        ) )

        ax_up.add_artist( ax_up.legend(
            handles = [
                *lines_colors, *lines_polyfit, *lines_linestyles
            ],
            loc="center right",
            fontsize=8,
            handlelength=3,
            numpoints=1,
            fancybox=True,
            markerscale = 1,
            ncol = 2,
            frameon = False,
            mode = None
        ) )

        ax_up.set_xlim(min_x, max_x)
        ax_up.set_ylim(min_y, max_y)

        ax_down.set_ylim(1e-3, 1.5e0)

        plt.savefig(
            'uniBar.eps', format="eps",
            bbox_inches='tight',
            dpi=1200,
            pad_inches=0
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
        fig = plt.figure(figsize=(8,6))
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
        x_format_str = "", y_format_str = ""
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

        ax.tick_params(direction="in")

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

        colors = [
            "b", "g", "r", "c", "y", "k"
        ]

        return random.sample(
            set(
                itertools.product(linestyles,colors,markerstyles, colors)
            ), 1
        )[0]

    @staticmethod
    def _get_specific_ms(map_me, label):
        """
        for provided list return list of dictionaries, each having the
        entry of the provided list as key and marker style as value

        it is here because to be easier to find
        """

        all_makrers_styles = [
            "s", "8", ">", "<", "^", "v", "o",
            "X", "P", "d", "D", "H", "h", "*", "p",
            "$\\bowtie$", "$\\clubsuit$", "$\\diamondsuit$", "$\\heartsuit$",
            "$\\spadesuit$",
            "$\\bigotimes$", "$\\bigoplus$",
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

        res = {
            _: __ for _, __ in zip(map_me, all_makrers_styles)
        }

        res.update({"label": label})

        return res

    def _get_specific_ls(self, map_me, label):
        """
        for provided list return list of dictionaries, each having the
        entry of the provided list as key and line style as value

        it is here because to be easier to find
        """

        if self.specific_ls:
            return self.specific_ls

        #~ all_line_styles = [
            #~ ":", "-.", "--",
            #~ (0, (5, 1, 1, 1, 1, 1)),
            #~ (0, (5, 1, 1, 1, 1, 1, 1, 1)),
            #~ (0, (5, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
            #~ (0, (8, 1, 1, 1, 3, 1, 1, 1))
        #~ ]
        all_line_styles = [
            "-.", "--",
            #~ (0, (5, 1, 1, 1, 1, 1)),
            (0, (5, 1, 1, 1, 1, 1, 1, 1)),
            #~ (0, (5, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
            (0, (8, 1, 1, 1, 3, 1, 1, 1))
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

        res = {
            _: __ for _, __ in zip(map_me, all_line_styles)
        }

        res.update({"label": label})

        return res

    def _get_specific_c(self, map_me, label):
        """
        for provided list return list of dictionaries, each having the
        entry of the provided list as key and colour as value

        it is here because to be easier to find
        """

        if self.specific_c:
            return self.specific_c

        #~ all_colors = [
            #~ "b", "g", "r", "c", "m", "y"
        #~ ]
        #~ all_colors = [
            #~ "#e6194B", "#3cb44b", "#4363d8"
        #~ ]

        all_colors = [
            "#800000", "#4363d8", "#f58231"
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

        res = {
            _: __ for _, __ in zip(map_me, all_colors)
        }

        res.update({"label": label})

        return res

    @staticmethod
    def _get_markevry(data_x, data_y, amount_points = 20):
        """
        for provided data list get the difference between the max and min to
        evluate the stem for getting amount_points and return list of
        indexes of points closes to the step

        solution used from
        https://stackoverflow.com/questions/9873626/choose-m-evenly-spaced-elements-from-a-sequence-of-length-n
        """

        def _get_EvenlySpacedIdxs(data, amount_points):

            evenly_spaced, step = np.linspace(
                min(data, key=abs), max(data, key=abs),
                num = amount_points,
                endpoint = True,
                retstep = True
            )

            step = abs(step)
            evenly_spaced_idx = [ 0 ]
            cumitv_step = 0
            for _, __ in enumerate(np.diff(data)):
                cumitv_step += abs(__)
                if cumitv_step > step:
                    evenly_spaced_idx.append(_)
                    cumitv_step = 0

            if abs(data[-1] - data[evenly_spaced_idx[-1]]) > step:
                evenly_spaced_idx.append(len(data)-1)

            return evenly_spaced_idx

        idx_x = _get_EvenlySpacedIdxs(data_x, amount_points)
        idx_y = _get_EvenlySpacedIdxs(data_y, amount_points)

        #~ return list( set(idx_x).union(idx_y) )
        all_i = list( set(idx_x).union(idx_y) )

        return np.linspace(
            min(all_i),
            max(all_i),
            amount_points,
            endpoint=True,
            dtype=int
        ).tolist()

    @staticmethod
    def _get_plot_keywords( markers, colors, linestyles, current ):
        """
        for provided dictionaries for markers, colors and linestyles

        for provided list <current> which contains the current eos name,
        m value and lambda value

        search where the eos name, m and lambda value are and subscribe its value

        the label is returned last
        """

        marker = None
        color = None
        linestyle = None

        for _ in current.keys():

            if _ == markers["label"]:
                marker = markers.get(current[_])

            elif _ == colors["label"]:
                color = colors.get(current[_])

            elif _ == linestyles["label"]:
                linestyle = linestyles.get(current[_])

        return {
            "color": color,
            "linestyle": linestyle,
            "marker": marker,
            "markerfacecolor": color,
            "markeredgecolor": color
        }

    def _get_MSs_Cs_LSs(self,_severalEOSs):

        tmp_ms = None
        if self.specific_ms:
            tmp_ms = self.specific_ms
        else:

            print("\n did you forgot to set_severalEOSs_ms_ls_c ?? \n")

            tmp_ms = self._get_specific_ms(
                list( set( [ _["name"] for _ in _severalEOSs ] ) ), "name"
            )

        tmp_c = None
        if self.specific_c:
            tmp_c = self.specific_c
        else:

            print("\n did you forgot to set_severalEOSs_ms_ls_c ?? \n")

            tmp_c = self._get_specific_c(
                list( set( [ _["m"] for _ in _severalEOSs ] ) ), "m"
            )

        tmp_ls = None
        if self.specific_ls:
            tmp_ls = self.specific_ls
        else:

            print("\n did you forgot to set_severalEOSs_ms_ls_c ?? \n")

            tmp_ls = self._get_specific_ls(
                list( set( [ _["lambda"] for _ in _severalEOSs ] ) ), "lambda"
            )

        return tmp_ms, tmp_c, tmp_ls

    @staticmethod
    def _get_lines_MSs_Cs_LSs(markers, colors, linestyles, severalEOSs = None ):

        def _convert_sci(num):

            snum = "{:.2e}".format(num)
            mantisa, power = snum.split("e")

            mantisa = "{}".format(mantisa)
            power = "{}".format(power)

            return " ${{{:.0f}}} \\times 10^{{{:.0f}}}$".format(
                float(mantisa), float(power)
            ) if float(mantisa) \
            else " ${{{:.0f}}}$".format(
                float(mantisa)
            )

        chs_lambda_m = lambda _: "$\\lambda =$ " if _ == "lambda" else "m = "

        if not severalEOSs:

            lines_markers = [
                Line2D(
                    [0], [0],
                    color="k",
                    marker = __,
                    linewidth = 0,
                    label = _
                )
                for _, __ in markers.items() if _ != "label"
            ]

            lines_colors = [
                Line2D(
                    [0], [0],
                    color=__,
                    marker = None,
                    linewidth = 1.5,
                    linestyle = "-",
                    label = chs_lambda_m(colors["label"]) + _convert_sci(_)
                )
                for _, __ in colors.items() if _ != "label"
            ]

            lines_linestyles = [
                Line2D(
                    [0], [0],
                    color="k",
                    marker = None,
                    linewidth = 1.5,
                    linestyle = __,
                    label = chs_lambda_m(linestyles["label"]) + _convert_sci(_)
                )
                for _, __ in linestyles.items() if _ != "label"
            ]
        else:
            lines_markers = [
                Line2D(
                    [0], [0],
                    color="k",
                    marker = markers[_],
                    linewidth = 0,
                    label = _
                )
                for _ in list( set( [ _["name"] for _ in severalEOSs ] ) )
            ]

            lines_colors = [
                Line2D(
                    [0], [0],
                    color=colors[_],
                    marker = None,
                    linewidth = 1.5,
                    linestyle = "-",
                    label = chs_lambda_m(colors["label"]) + _convert_sci(_)
                )
                for _ in list( set( [ _[colors["label"]] for _ in severalEOSs ] ) )
            ]

            lines_linestyles = [
                Line2D(
                    [0], [0],
                    color="k",
                    marker = None,
                    linewidth = 1.5,
                    linestyle = linestyles[_],

                    label = chs_lambda_m(linestyles["label"]) + _convert_sci(_)
                )
                for _ in list( set( [ _[linestyles["label"]] for _ in severalEOSs ] ) )
            ]

        return lines_markers, lines_colors, lines_linestyles


if __name__ == "__main__":

    print("\n asd \n")
