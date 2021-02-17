import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn


def plot_nxy(name, tfile,
             title = '',
             plim = 1e-2,
             lines = [],
             xscale = 'log',
             ylim = None,
             xlim = None):
    """ Plot a list of trajectories.

    Args:
      name (str): Name of the pdf file.
      tfile (str): Filename of a treekin `nxy` output file.
      title (str, optional): Name of the title for the plot.
      plim (float, optional): Minimal occupancy to plot a trajectory. Defaults to 0.01
      lines ([int,..], optional): Selected list of lines to be plotted.
      xscale (str, optional): *lin* or *log*. Default: *log*.
      xlim ((float,float), optional): matplotlib xlim.
      ylim ((float,float), optional): matplotlib ylim.

    Returns:
      [str]: Name of the output file.
    """

    lines = set(lines)
    title = title if title else name

    nxy = []
    with open(tfile) as tkn:
        for line in tkn:
            if re.match('#', line):
                continue
            nxy.append(list(map(float, line.strip().split())))

    # Do the plotting
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['right'].set_visible(False)
    ax.set_xscale(xscale)

    if ylim: ax.set_ylim(ylim)
    if xlim: ax.set_xlim(xlim)

    for e, traject in enumerate(zip(*nxy)):
        if e == 0:
            time = traject
            continue
        if lines and e not in lines:
            continue
        if plim and max(traject) < plim:
            continue
        p, = ax.plot(time, traject, '-', lw = 1.5)
        p.set_label("ID {:d}".format(e))

    fig.set_size_inches(7, 3)
    fig.text(0.5, 0.95, title, ha='center', va='center')
    ax.set_ylabel('occupancy', fontsize = 11)
    ax.set_xlabel('time [s]', ha = 'center', va = 'center', fontsize = 11)

    # Add ticks for 1 minute, 1 hour, 1 day, 1 year
    #ax.axvline(x = 60, linewidth = 1, color = 'black', linestyle = '--') # 1 minute
    #ax.axvline(x = 3600, linewidth = 1, color = 'black', linestyle = '--') # 1 hour
    #ax.axvline(x = 86400, linewidth = 1, color = 'black', linestyle = '--') # 1 day
    #ax.axvline(x = 31536000, linewidth = 1, color = 'black', linestyle = '--') # 1 year
    plt.legend()

    plt.savefig(name, bbox_inches='tight')
    plt.close()
    return name

def plot_nxy_linlog(name, tfile,
             xlim,
             title = '',
             plim = 1e-2,
             lines = [],
             ylim = None,
             figwidth = 7,
             figheight = 3,
             figdivide = 3.5):
    """ Plot a list of trajectories.

    Args:
      name (str): Name of the pdf file.
      tfile (str): Filename of a treekin `nxy` output file.
      xlim (float, float, float): Time poins for start, 
        lin-to-log switch and end.
      title (str, optional): Name of the title for the plot.
      plim (float, optional): Minimal occupancy to plot a trajectory. Defaults to 0.01
      lines ([int,..], optional): Selected list of lines to be plotted.
      ylim ((float,float), optional): matplotlib xlim.
      figdivide(float, optional): Set the size of the log-part of the plot.

    Returns:
      [str]: Name of the output file.
    """

    lines = set(lines)
    title = title if title else name

    nxy = []
    with open(tfile) as tkn:
        for line in tkn:
            if re.match('#', line):
                continue
            tr = list(map(float, line.strip().split()))
            maxtime = tr[0]
            nxy.append(tr)

    # Get the relevant arguments from args
    lintime = xlim[1]
    logtime = max(xlim[2] - lintime, lintime * 10)
    
    # Do the plotting
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['right'].set_visible(False)
    if ylim:
        ax.set_ylim(ylim)
        #ax.set_ylim([-0.05, 1.05])
    ax.set_xscale('linear')

    # Make the second part of the plot logarithmic
    offset = 0.00001
    ax.set_xlim((xlim[0], xlim[1] + offset))
    divider = make_axes_locatable(ax)
    axLog = divider.append_axes("right", size = figdivide, pad = 0, sharey = ax)
    axLog.set_xscale('log')
    axLog.set_xlim((xlim[1] + offset, logtime))
    if ylim:
        ax.set_ylim(ylim)
        #ax.set_ylim([-0.05, 1.05])
    axLog.yaxis.set_visible(False)
    axLog.spines['left'].set_visible(False)

    for e, traject in enumerate(zip(*nxy)):
        if e == 0:
            time = traject
            continue
        if lines and e not in lines:
            continue
        if plim and max(traject) < plim:
            continue
        p, = ax.plot(time, traject, '-', lw = 0.5)
        L, = axLog.plot(time, traject, '-', lw = 0.5)
        L.set_label("ID {:d}".format(e))

    fig.set_size_inches(figwidth, figheight)
    fig.text(0.5, 0.95, title, ha='center', va='center')
    ax.set_ylabel('occupancy', fontsize = 11)
    ax.set_xlabel('time [s]', ha = 'center', va = 'center', fontsize = 11)
    ax.xaxis.set_label_coords(.9, -0.2)

    # Add ticks for 1 minute, 1 hour, 1 day, 1 year
    axLog.axvline(x = lintime, linewidth = 3, color = 'black', linestyle = '-') 
    axLog.axvline(x = 60, linewidth = 1, color = 'black', linestyle = '--') # 1 minute
    axLog.axvline(x = 3600, linewidth = 1, color = 'black', linestyle = '--') # 1 hour
    axLog.axvline(x = 86400, linewidth = 1, color = 'black', linestyle = '--') # 1 day
    axLog.axvline(x = 31536000, linewidth = 1, color = 'black', linestyle = '--') # 1 year
    plt.legend()

    plt.savefig(name, bbox_inches='tight')
    plt.close()
    return name


def plot_xmgrace(trajectories, filename):
    head = """
@with line
@line on
@line loctype world
@line g0
@line linewidth 2
@line linestyle 1
@line color 7
@line arrow 0
@line arrow type 0
@line arrow length 1.000000
@line arrow layout 1.000000, 1.000000
@line def
"""
    with open(filename, 'w') as gfh:
        gfh.write(head)
        for nid in sorted(trajectories):
            course = trajectories[nid]
            t, o = list(zip(*course))
            for i in range(len(t)):
                gfh.write("{:f} {:f}\n".format(t[i], o[i]))
            gfh.write("&\n")
    return


def plot_simulation(trajectories, basename, formats, 
                    lin_time, log_time, title = ''):
    """
    """
    seaborn.set_style("darkgrid")
    assert log_time >= lin_time * 10

    # Do the plotting
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['right'].set_visible(False)
    ax.set_ylim([-0.05, 1.05])
    ax.set_xscale('linear')

    # Make the second part of the plot logarithmic
    offset = 0.00001
    ax.set_xlim((0, lin_time + offset))
    divider = make_axes_locatable(ax)
    axLog = divider.append_axes("right", size=2.5, pad=0, sharey=ax)
    axLog.set_xscale('log')
    axLog.set_xlim((lin_time + offset, log_time))
    axLog.set_ylim([-0.05, 1.05])
    axLog.yaxis.set_visible(False)
    axLog.spines['left'].set_visible(False)

    for ni in sorted(trajectories):
        course = trajectories[ni]
        t, o = list(zip(*course))
        # Determine which lines are part of the legend:
        # like this, it is only those that are populated
        # at the end of transcription and if they reach
        # an occupancy of 10% or higher
        if t[-1] > lin_time:
            p, = ax.plot(t, o, '-', lw=1.5)
            L, = axLog.plot(t, o, '-', lw=1.5)
            if max(o) >= 0.1:
                L.set_label(f"ID {ni}")
        else:
            p, = ax.plot(t, o, '-', lw=0.5)
            L, = axLog.plot(t, o, '-', lw=0.5)

    fig.set_size_inches(7, 3)
    fig.text(0.5, 0.95, title, ha='center', va='center')

    # Add ticks for 1 minute, 1 hour, 1 day, 1 year
    axLog.axvline(x = lin_time, linewidth = 3, color = 'black', linestyle = '-') 
    axLog.axvline(x = 60, linewidth = 1, color = 'black', linestyle = '--') # 1 minute
    axLog.axvline(x = 3600, linewidth = 1, color = 'black', linestyle = '--') # 1 hour
    axLog.axvline(x = 86400, linewidth = 1, color = 'black', linestyle = '--') # 1 day
    axLog.axvline(x = 31536000, linewidth = 1, color = 'black', linestyle = '--') # 1 year
    plt.legend()

    ax.set_ylabel('occupancy', fontsize = 11)
    ax.set_xlabel('time [s]', ha = 'center', va = 'center', fontsize = 11)
    ax.xaxis.set_label_coords(.9, -0.15)

    for ending in formats:
        pfile = basename + '.' + ending
        plt.savefig(pfile, bbox_inches = 'tight')
    return

def main():
    """ DrPlotter.
    """
    import sys
    llen = 0
    lin_time = 0
    all_courses = dict()
    for e, line in enumerate(sys.stdin):
        if e == 0:
            continue
        [ni, time, occu, ss, en] = line.split()
        time = float(time)
        occu = float(occu)
        if ni not in all_courses:
            all_courses[ni] = [(time, occu)]
        else:
            all_courses[ni].append((time, occu))
        if len(ss) > llen:
            lin_time = time
        llen = len(ss)
    log_time = time if time > 10 * lin_time else 10 * lin_time

    plot_xmgrace(all_courses, 'DrPlotter.gr')
    plot_simulation(all_courses, 'DrPlotter', ['pdf'], lin_time, log_time, title = 'DrPlotter')
 
if __name__ == '__main__':
    main()
