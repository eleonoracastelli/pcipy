from matplotlib import pyplot as plt
import matplotlib as mpl


def plotconfig(lbsize=20, lgsize=18, autolayout=True, figsize=[8, 6],
               ticklabelsize=16, style='classic', fsize=None,
               tdir=None, major=None, minor=None, lwidth=None, lhandle=None,
               fontfamily=None):

    plt.style.use('default')
    if style == 'classic':
        
        if fontfamily is None:
            fontfamily = 'STIXGeneral'
        if fsize is None:
            fsize=15

        ticks_font = mpl.font_manager.FontProperties(family='serif',
                                                    style='normal',
                                                    weight='normal',
                                                    stretch='normal',
                                                    size=lbsize)

        mpl.rcParams['xtick.labelsize'] = ticklabelsize
        mpl.rcParams['ytick.labelsize'] = ticklabelsize
        mpl.rcParams['font.size'] = fsize
        mpl.rcParams['figure.autolayout'] = autolayout
        # mpl.rcParams['figure.figsize'] = 7.2, 4.45
        mpl.rcParams['figure.figsize'] = figsize[0], figsize[1]
        mpl.rcParams['axes.titlesize'] = lbsize
        mpl.rcParams['axes.labelsize'] = lbsize
        mpl.rcParams['lines.linewidth'] = 2
        mpl.rcParams['lines.markersize'] = 6
        mpl.rcParams['legend.fontsize'] = lgsize
        mpl.rcParams['mathtext.fontset'] = 'stix'
        mpl.rcParams['font.family'] = 'STIXGeneral'

    elif style == 'dataviz':

        mpl.rcParams['figure.autolayout'] = autolayout
        mpl.rcParams['figure.figsize'] = figsize[0], figsize[1]

        if fsize is None:
            fsize=15
        if tdir is None:
            tdir='in'
        if major is None:
            major=5.0
        if minor is None:
            minor=3.0
        if lwidth is None:
            lwidth=0.8
        if lhandle is None:
            lhandle=2.0
        
        mpl.rcParams['xtick.labelsize'] = ticklabelsize
        mpl.rcParams['ytick.labelsize'] = ticklabelsize
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.size'] = fsize
        plt.rcParams['legend.fontsize'] = lgsize
        plt.rcParams['xtick.direction'] = tdir
        plt.rcParams['ytick.direction'] = tdir
        plt.rcParams['xtick.major.size'] = major
        plt.rcParams['xtick.minor.size'] = minor
        plt.rcParams['ytick.major.size'] = major
        plt.rcParams['ytick.minor.size'] = minor
        plt.rcParams['axes.linewidth'] = lwidth
        plt.rcParams['legend.handlelength'] = lhandle
        

    elif style == 'publication':

        if fsize is None:
            fsize=18
        if tdir is None:
            tdir='in'
        if major is None:
            major=10
        if minor is None:
            minor=7
        if lwidth is None:
            lwidth=2
        if lhandle is None:
            lhandle=2.0

        if fontfamily is None:
            fontfamily = 'Avenir'
        elif fontfamily == 'STIXGeneral':
            mpl.rcParams['mathtext.fontset'] = 'stix'

        mpl.rcParams['font.family'] = fontfamily
        mpl.rcParams['figure.figsize'] = figsize[0], figsize[1]
        mpl.rcParams['figure.autolayout'] = autolayout
        mpl.rcParams['xtick.labelsize'] = ticklabelsize
        mpl.rcParams['ytick.labelsize'] = ticklabelsize
        plt.rcParams['font.size'] = fsize
        plt.rcParams['axes.linewidth'] = lwidth

        mpl.rcParams['axes.titlesize'] = lbsize
        mpl.rcParams['axes.labelsize'] = lbsize

        plt.rcParams['xtick.major.size'] = major
        plt.rcParams['xtick.minor.size'] = minor

        plt.rcParams['ytick.major.size'] = major
        plt.rcParams['ytick.minor.size'] = minor

        plt.rcParams['xtick.direction'] = tdir
        plt.rcParams['ytick.direction'] = tdir

        plt.rcParams['xtick.major.width'] = lwidth
        plt.rcParams['ytick.major.width'] = lwidth

        plt.rcParams['xtick.major.top'] = True
        plt.rcParams['xtick.minor.top'] = True
        plt.rcParams['ytick.major.right'] = True
        plt.rcParams['ytick.minor.right'] = True

def comparison_plot(x_list, matrix_list,
                    labels=['Initial', 'Estimated', 'True', 'Eigenvalues'],
                    colors=['gray', 'blue', 'black', 'green'],
                    linewidths=[2, 2, 2, 2],
                    linestyles=['solid', 'solid', 'solid', 'solid'],
                    xscale='log', 
                    yscale='log',
                    xlabel="Frequency [Hz]",
                    ylabels=["PSD 1", "PSD 2", "PSD 3"],
                    titles=[None, None, None],
                    tight=False,
                    bottom=None, top=None,
                    left=None, right=None,
                    figsize=[8, 6],
                    lbsize=17,
                    lgsize=13,
                    sharex=True,
                    sharey=True,
                    frameon=True,
                    bbox_to_anchor=(0, 1.03, 1, 0.3),
                    rasterized=False,
                    ncol=None):

    n_channels = matrix_list[0].shape[1]

    plotconfig(lbsize=lbsize, lgsize=lgsize, autolayout=True, figsize=figsize)
    # Frequency plot
    fig1, ax1 = plt.subplots(nrows=n_channels, sharex=sharex, sharey=sharey)
    for i in range(n_channels):
        [ax1[i].plot(x_list[j], matrix_list[j][:, i],
                     colors[j], label=labels[j],
                     linewidth=linewidths[j],
                     linestyle=linestyles[j],
                     rasterized=rasterized)
         for j in range(len(matrix_list))]
        ax1[i].set_xscale(xscale)
        ax1[i].set_yscale(yscale)
        ax1[i].set_ylabel(ylabels[i], fontsize=lbsize)
        ax1[i].set_xlim(left=left, right=right)
        ax1[i].set_ylim(bottom=bottom, top=top)
        ax1[i].set_title(titles[i])
        ax1[i].minorticks_on()
        if i == 0:
            if ncol is None:
                # Number of columns should not be larger than 3
                if len(matrix_list) <= 3:
                    ncol = len(matrix_list)
                else:
                    ncol = 3
            # Locating legend (x, y, width, height) from bottom left corner
            ax1[i].legend(bbox_to_anchor=bbox_to_anchor,
                          loc="lower left",
                          mode='expand',
                          shadow=False,
                          ncol=ncol,
                          fontsize=lgsize, frameon=frameon)
    ax1[n_channels-1].set_xlabel(xlabel, fontsize=lbsize)
    if tight:
        plt.tight_layout()
    # plt.show()
    return fig1, ax1
