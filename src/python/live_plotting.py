#!/usr/bin/env python

DEF_RES_PATH = "/home/dimitar/projects/STT_theories/results/APR4"

PREVIUS_FILE_BETA = 0
PREVIUS_FILE_M = 0
PREVIUS_FILE_LAMBDA = 0

#~ DEF_RES_PATH_PREV = "/home/dimitar/projects/STT_theories/results/" \
    #~ + "Results_Statia_1/" \
    #~ + "beta{:.0f}/".format(PREVIUS_FILE_BETA) \
    #~ + "lambda{:.0e}".format(PREVIUS_FILE_LAMBDA)

#~ DEF_RES_FNAME_PREV = "STT_phiScal_J_AkmalPR_" \
    #~ + "beta{:.3e}_".format(PREVIUS_FILE_BETA) \
    #~ + "m{:.3e}_".format(PREVIUS_FILE_M) \
    #~ + "lambda{:.3e}".format(PREVIUS_FILE_LAMBDA)

DEF_RES_PATH_PREV = "/home/dimitar/projects/STT_theories/results/APR4" \
    + "" \
    + "" \
    + ""

DEF_RES_FNAME_PREV = "STT_phiScal_J_APR4_" \
    + "beta{:.3e}_".format(PREVIUS_FILE_BETA) \
    + "m{:.3e}_".format(PREVIUS_FILE_M) \
    + "lambda{:.3e}".format(PREVIUS_FILE_LAMBDA)


_kalin_beta = "beta-6"

DEF_RES_KALIN = "/home/dimitar/Documents/Teaching_Materials/University/" \
  + "Uvod.Fizika.Cherni.Dupki/Doktorant/Moi_Statii_Prezentacii/Statia_3/Kalin_data/" \
  + _kalin_beta + "/"

KALIN_BETA = _kalin_beta
KALIN_LAMBDA = "lambda{}_".format(PREVIUS_FILE_LAMBDA)
KALIN_M = "m{}".format(PREVIUS_FILE_M)

def _check_if_StrNotBlank(string):
        """
        check if a sting is blank/empty

        Parameters
        ----------

        Returns
        -------
        : boolean
            True if string is not blank/empty
            False if string is blank/empty
        """

        return bool(string and string.strip())

def _get_latest_file(
    path = DEF_RES_PATH,
    model = "STT_phiScal_J_"
):
    """
    return the name of the latest modified file in path
    by creating wildcard of the type model + "*"

    the default values are where the results are plus presume that we will
    plot the result file only

    Parameters
    ----------
    path: string
        the path to the directory of interest
    model: string
        how the file name looks like, it serves as wildcard
        of the type model + "*"

    Returns
    -------
    : string
        the name of latest modified file
    """

    import glob
    import os

    if (_check_if_StrNotBlank(path) and
        _check_if_StrNotBlank(model)):

        filename = max(
            glob.glob(path + "/" + model + "*"),
            key=os.path.getctime
        )

        filename = filename.split("/")[-1]

        #~ print(
            #~ "\n latest file in \n\t {} is: {} \n".format(
                #~ path, filename
            #~ )
        #~ )
    else:

        filename = ""

        print("\n path and name not set, some filename var!!! \n")

    return filename

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

def _load_kalin_file(
    data=[], headline=[], label=[],
    path=DEF_RES_KALIN,
    filename=KALIN_BETA+KALIN_LAMBDA+KALIN_M+".txt"
):
    """
    load kalin result file with name <filename> by appending the data
    after taking care of coefficients to <data> and saving the <headlines>(
    what is the name of each column) and label - the filename itself

    Parameters
    ----------
    path: string
        the full path of the directory where the file is
    filename: string
        the name of the file itself
    data: list
        where to append the data of interest
    headline: list
        what each column is
    label: list
        the name the given data, expect it will be the filename

    Returns
    -------
    """

    import os.path

    all_data = []

    if os.path.isfile(path+filename):
        with open(path + filename, "r") as f:
            all_data = f.readlines()
    else:
        print(
            "from where to load a kalin file!!! {}{} pass ".format(
                path, filename
            )
        )

        pass

    #~ get rid of inline text
    all_data = [
        line.strip() for line in all_data
        if "lambda" not in line and "f_rot" not in line and line.strip()
    ]

    #~ i am not interested in all data, so only the indexes here
    #~ will be saved
    #~ since in some files the central pressure is in different column
    #~ we check it and take into account
    indexes = [3, 5, 7, 10, 6, 2]

    _units = _units_coef_clac()

    units = [
        _units["density"], _units["rad"],
        1, 1, 1, 1
    ]

    data.append([[] for i in indexes])

    for cnt, line in enumerate(all_data):
        for i, u, d in zip(indexes, units, data[-1]):
            try:
                d.append(float(line.split(" ")[i]) / u)
            except ValueError:
                print(
                    "\n ValueError: line {}: {} \n".format(
                        cnt, line
                    )
                )
                break

    for i, val in enumerate(data[-1][4]):
        data[-1][4][i] *= (-1) if val > 0 else 1


    label.append(filename)

    return { "data": data[-1], "label": label[-1] }

def _load_GR(
    data = [],headline = [], label = [], path = DEF_RES_PATH, filename = ""
):
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
    -------
    : dictionary
        dictionary with entries
            { "data" : data, "headline": headline, "label": label }
    """

    data.clear()
    headline.clear()
    label.clear()

    if not(_check_if_StrNotBlank(path) and
           _check_if_StrNotBlank(filename)):

        print(
            "\n from where to load a file!!! \n\t {} /{} \n".format(
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
            data[-1], [
                float(i) for i in line.strip().split(" ") if len(i.strip())
            ]
        ):
            apnd.append(num)

    return { "data": data[-1], "headline": headline[-1], "label": label[-1] }

def _load_previus_data(
    data = [],headline = [], label = [],
    path = DEF_RES_PATH_PREV,
    filename = DEF_RES_FNAME_PREV
):
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
    -------
    : dictionary
        dictionary with entries
            { "data" : data, "headline": headline, "label": label }
    """

    data.clear()
    headline.clear()
    label.clear()

    if not(_check_if_StrNotBlank(path) and
           _check_if_StrNotBlank(filename)):

        print(
            "\n from where to load a file!!! \n\t {} /{} \n".format(
                path, filename
            )
        )

    try:
        with open(path + "/" + filename, "r") as f:
            all_data = f.readlines()

    except FileNotFoundError:
        print(path + "/" + filename," MISSING !!! ")
        return {  }

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
            data[-1], [
                float(i) for i in line.strip().split(" ") if len(i.strip())
            ]
        ):
            apnd.append(num)

    return { "data": data[-1], "headline": headline[-1], "label": label[-1] }

def _load_file(
    filename = "", data = [], headline = [], label = [], path = DEF_RES_PATH
):
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
    -------
    : dictionary
        dictionary with entries
            { "data" : data, "headline": headline, "label": label }
    """

    data.clear()
    headline.clear()
    label.clear()

    if not(_check_if_StrNotBlank(path)):

        print(
            "\n from where to load a file!!! \n\t {} /{} \n".format(
                path, filename
            )
        )

    filename = _get_latest_file()

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

    try:
        data.append(
            [
                [] for i in all_data[0].strip().split(" ") if len(i.strip())
            ]
        )
    except IndexError:
        pass

    for cnt, line in enumerate(all_data):

        for apnd, num in zip(
            data[-1], [
                float(i) for i in line.strip().split(" ") if len(i.strip())
            ]
        ):
            apnd.append(num)

    try:
        return { "data": data[-1], "headline": headline[-1], "label": label[-1] }

    except IndexError:
        return { "data": [], "headline": [], "label": [] }

def _set_parms(ax, label_x, label_y):
    """
    clear the plot <ax> and set its Y and X axis respectively to
    <label_y> and <label_x>
    also set the font and thick size to predefined values

    Parameters
    ----------
    ax: matplotlib.pyplot.axes class
        the associated axes

    label_x: string
        label for X axis

    label_y: string
        label for Y axis

    Returns
    -------
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

def _create_plot_parms_dict():
    """
    creates dictionary for the parameters of interests for the plot

    Parameters
    ----------

    Returns
    -------
    : dictionary
        empty dictionary for the parameters of interest for plots
    """

    keys = [
        "linestyle", "linewidth", "label", "marker", "markersize", "alpha"
    ]

    return { key: None for key in keys }

def _create_data_mapping():
    """
    create mapping dictionary

    Parameters
    ----------

    Returns
    -------
    : dictionary
        dictionary with the mapping for the result data
    """

    return {
        "p_c": 0,
        "phiScal_c": 1,
        "M": 2,
        "AR": 3,
        "rho_c": 4,
        "J": 5
    }

def _create_Kalin_data_mapping():
    """
    create mapping dictionary

    Parameters
    ----------

    Returns
    -------
    : dictionary
        dictionary with the mapping for the result data
    """

    return {
        "rho_c": 0,
        "AR": 1,
        "M": 2,
        "J": 3,
        "phiScal_c": 4,
        "p_c": 5
    }

def _plot_data(ax, x, y, parms):
    """
    function to plot the data in the provided <ax> using the lists
    in <x> and <y> for X and Y axis

    Parameters
    ----------
    ax: matplotlib.pyplot.axes class
        where to be plotted

    x: list
        the data for X axis

    y: list
        the data for the Y axis

    parms: dictionary
        parameters for the plot given in the form of dictionary
            "linestyle" is the line style
            "linewidth" is the line width
            "label" is the label for the graph
            "marker" is the marker style
            "markersize" is the marker size
            "alpha" is the alpha of the line

    Returns
    -------
    """

    ax.plot(
        x,
        y,
        linewidth=parms["linewidth"],
        linestyle=parms["linestyle"],
        marker=parms["marker"],
        markersize=parms["markersize"],
        alpha=parms["alpha"],
        label=parms["label"]
    )

    return

def _set_parms_PhiScalMARJ(ax):
    """
    function which specificity sets the list of axes <ax> to
        p_c vs phiScal_c, AR vs M, M vs J

    Parameters
    ----------
        ax: list
            list of matplotlib.pyplot.axes

    Returns
    -------
    """

    ax_phiScal_c, ax_M_AR, ax_J_M = ax

    _set_parms(ax_phiScal_c, "p_c", "phiScal_c")
    _set_parms(ax_M_AR, "AR", "M")
    _set_parms(ax_J_M, "M", "J")

    return

def _plot_PhiScalMARJ(ax, data, data_GR, data_Kalin, data_previus):
    """
    function which specificity sets the <data> to each axis in the <ax>
        p_c vs phiScal_c, AR vs M, M vs J

    Parameters
    ----------
        ax: list
            list of matplotlib.pyplot.axes

        data: dictionary with data
            nested list of the updated data

        data_GR: dictionary with data
            nested list containing data for GR

    Returns
    -------
    """

    ax_phiScal_c, ax_M_AR, ax_J_M = ax

    mapping = _create_data_mapping()
    mapping_kalin = _create_Kalin_data_mapping()

    parms = _create_plot_parms_dict()

    if data_GR["data"]:

        parms["linestyle"] = "-"
        parms["linewidth"] = 1.5
        parms["label"] = "GR"
        parms["marker"] = ""
        parms["markersize"] = 0
        parms["alpha"] = 1

        _plot_data(
            ax_phiScal_c,
            data_GR["data"][mapping["p_c"]],
            data_GR["data"][mapping["phiScal_c"]],
            parms
        )

        _plot_data(
            ax_M_AR,
            data_GR["data"][mapping["AR"]],
            data_GR["data"][mapping["M"]],
            parms
        )

        _plot_data(
            ax_J_M,
            data_GR["data"][mapping["M"]],
            data_GR["data"][mapping["J"]],
            parms
        )

    if data_previus:

        parms["linestyle"] = "-"
        parms["linewidth"] = 1.5
        parms["label"] = DEF_RES_FNAME_PREV
        parms["marker"] = ""
        parms["markersize"] = 0
        parms["alpha"] = 1

        _plot_data(
            ax_phiScal_c,
            data_previus["data"][mapping["p_c"]],
            data_previus["data"][mapping["phiScal_c"]],
            parms
        )

        _plot_data(
            ax_M_AR,
            data_previus["data"][mapping["AR"]],
            data_previus["data"][mapping["M"]],
            parms
        )

        _plot_data(
            ax_J_M,
            data_previus["data"][mapping["M"]],
            data_previus["data"][mapping["J"]],
            parms
        )

    if data_Kalin:

        parms["linestyle"] = "-"
        parms["linewidth"] = 1.5
        parms["label"] = data_Kalin["label"]
        parms["marker"] = ""
        parms["markersize"] = 0
        parms["alpha"] = 1

        _plot_data(
            ax_phiScal_c,
            data_Kalin["data"][mapping_kalin["p_c"]],
            data_Kalin["data"][mapping_kalin["phiScal_c"]],
            parms
        )

        _plot_data(
            ax_M_AR,
            data_Kalin["data"][mapping_kalin["AR"]],
            data_Kalin["data"][mapping_kalin["M"]],
            parms
        )

        _plot_data(
            ax_J_M,
            data_Kalin["data"][mapping_kalin["M"]],
            data_Kalin["data"][mapping_kalin["J"]],
            parms
        )

    if data["data"]:

        parms["linestyle"] = "--"
        parms["linewidth"] = 1
        parms["label"] = data["label"]
        parms["marker"] = "o"
        parms["markersize"] = 5
        parms["alpha"] = 0.5

        _plot_data(
            ax_phiScal_c,
            data["data"][mapping["p_c"]],
            data["data"][mapping["phiScal_c"]],
            parms
        )

        _plot_data(
            ax_M_AR,
            data["data"][mapping["AR"]],
            data["data"][mapping["M"]],
            parms
        )

        _plot_data(
            ax_J_M,
            data["data"][mapping["M"]],
            data["data"][mapping["J"]],
            parms
        )
    else:
        pass

    ax_phiScal_c.legend(loc="best",fontsize=8)
    ax_M_AR.legend(loc="best",fontsize=8)
    ax_J_M.legend(loc="best",fontsize=8)

    return

def _set_PhiScalMARJ_grid_placement(fig, gs):
    """
    using the add_subplot method create the following fig layout

    +-------------------------------+
    |               |               |
    |               |               |
    |    M vs AR    |   J vs M      |
    |   all_ax[1]   |   all_ax[2]   |
    |               |               |
    |-------------------------------|
    |                               |
    |                               |
    |         phiScal_c vs p_c      |
    |           all_ax[0]           |
    |                               |
    +-------------------------------+

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
    all_ax.append(fig.add_subplot(gs[1:,:]))
    all_ax.append(fig.add_subplot(gs[0,0]))
    all_ax.append(fig.add_subplot(gs[0,1]))

    return all_ax

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

def update(frame, ax, set_parms, plot_me, data_GR, data_Kalin, data_previus):
    """
    the function which will clear the all axes ax_*
    will place the GR reference
    and will place the new data from frame, which is loaded by _load_file

    Parameters
    ----------
    frame: dictionary
        the result of _load_file

    ax: list
        list of matplotlib.pyplot.axes which will be used

    set_parms: function
        function to clear plot and set parameters of each axes

    plot_me: function
        function to plot <frame> to each <ax>

    data_GR: list
        nested list of the GR solutions for reference

    Returns
    -------
    """

    if not frame["data"]:
        pass

    set_parms(ax)

    plot_me(ax, frame, data_GR, data_Kalin, data_previus)

    return

def update_data():

    yield _load_file()

if __name__ == "__main__":

    from matplotlib import pyplot as plt
    from matplotlib.animation import FuncAnimation

    data_GR = _load_GR(
        filename = "STT_phiScal_J_APR4_beta0.000e+00_m0.000e+00_lambda0.000e+00"
    )

    data_Kalin = _load_kalin_file()
    #~ do not want kalin data now
    data_Kalin.clear()

    data_previus = _load_previus_data()
    #~ data_previus.clear()

    fig, ax = _get_figure(
        nrows=2,
        ncols=2,
        grid_placement=_set_PhiScalMARJ_grid_placement
    )

    animation_live = FuncAnimation(
        fig = fig,
        func = update,
        fargs = (
          ax, _set_parms_PhiScalMARJ, _plot_PhiScalMARJ,
          data_GR, data_Kalin, data_previus
        ),
        frames = update_data,
        interval = 500,
        save_count = 0,
        repeat=False
    )

    plt.show()
