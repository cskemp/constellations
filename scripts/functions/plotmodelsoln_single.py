# plotmodelsoln.py
# show star clusters (including edges, if supplied) 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.lines import Line2D
from matplotlib import colors, cm, gridspec
import functions.makecolormap as makecolormap
import functions.starhelpers as sh
import re
import os
from matplotlib.transforms import Bbox

def mag2plotsize(m):
    return 100 * (1.8**(6-m))/(1.8**4)

def plot_edges(thetas, ras, ris, cis, ax, col=[0.5, 0.5, 0.5, 1], count = 1, lwid=1):
    for s1, s2, c, cnt in zip(ris, cis, col, count):
        # deal with RA wrap around
        #print (s1, s2)
        if thetas[s1] > thetas[s2]:
            s1, s2 = s2, s1
        th1 = thetas[s1]
        th2 = thetas[s2]
        th1a, th2a = 0, 0
        if thetas[s1] < -np.pi / 2 and thetas[s2] > (1 / 2) * np.pi:
            th1a = th1 + 2 * np.pi
            th2a = th2 - 2 * np.pi
        if th1a == 0:
            ax.add_line(Line2D([th1, th2], [ras[s1], ras[s2]], linewidth=cnt, color=c, zorder=1))
        else:
            ax.add_line(Line2D([th1a, th2], [ras[s1], ras[s2]], linewidth=cnt, color=c, zorder=1))
            ax.add_line(Line2D([th1, th2a], [ras[s1], ras[s2]], linewidth=cnt, color=c, zorder=1))


def plot_graph(ax, sg, magthresh = 4.0, edgecount= 0, showlabs=0, ethick = 3):

    df = pd.DataFrame.from_dict(dict(sg.nodes.data()), orient='index')
    # filter out faint stars
    df= df[df['m'] <= magthresh]
    # convert magnitude to plot size
    df = df.assign(m = mag2plotsize(df.m))
    # transform so that Orion on left side of plot
    ramap = lambda x:  np.pi - (x + 1.45 * np.pi) if x < 0.55 *np.pi else np.pi -  (x - 0.55*np.pi)
    df.ra = df.ra.apply(ramap)

    edf = pd.DataFrame(list(sg.edges.data()))
    edf = edf.join(pd.DataFrame(edf[2].to_dict()).T)
    if edgecount:
        edf = edf[[0,1,'col', 'ecount']]
    else:
        edf = edf[[0,1,'col']]
        # for thicker edges set to 5
        # use 1 for empty graph
        edf['ecount']=ethick

    plot_edges(df.ra, df.dec, edf[0], edf[1], ax, col=edf['col'], count=edf['ecount'])
    colmatrix = np.array(df['col'].values.tolist())
    ax.scatter(df.ra, df.dec, s=df.m, c=colmatrix, zorder=2)
    if showlabs:
        for index, row in df.iterrows():
            ax.text(row.ra, row.dec, row['name'])

    return

def annotate_nodes(sg, clusters, cols):
    nx.set_node_attributes(sg, cols[0,], 'col')
    for i, cl in enumerate(clusters):
        icol = cols[i,]
        for clmem in cl:
            sg.add_node(clmem, col=icol)
    return

def annotate_edges(sg, edges, clusters, cols, edgecol):
    if edgecol:
        nx.set_edge_attributes(sg, cols[0,], 'col')
    else:  # set edges to white
        nx.set_edge_attributes(sg, [1,1,1,0], 'col')
    if edges:
        for i, cl in enumerate(clusters):
            icol = cols[i,]
            if i>0:
                for u, v in edges[i].edges:
                    sg.add_edge(u, v, col=icol)
    return



decsat = lambda x : x - 0.3 if x > 0.7 else x
def increasesat(rgbcol):
    hsvcol = colors.rgb_to_hsv(rgbcol[0:3])
    rgbnew = colors.hsv_to_rgb( (hsvcol[0], decsat(hsvcol[1]), hsvcol[2] ) )
    return (rgbnew[0], rgbnew[1], rgbnew[2], rgbcol[3])


def preparecolormap_aligned(stargraph, clusters, samecol, ncluster):
    if not clusters:
        cmap = np.array([colors.to_rgba('gray')])
        return(cmap)
    elif samecol:
        cmap = np.array([colors.to_rgba('red') for c in clusters])
        cmap[0] = colors.to_rgba('gray')
        return(cmap)

    fixedc = sh.colorgroups(stargraph)
    cmap = [colors.to_rgba('gray') for c in clusters]
    cassigned= np.array([0 for c in clusters])
    cassigned[0] = 1 # keep gray in first  position
    for (nodei, col) in fixedc.items():
        nodeics = [(len(c), i) for i, c in enumerate(clusters) if nodei in c and cassigned[i]==0]
        if len(nodeics) > 0:
            nodeics.sort()
            cindex = nodeics[-1][1]
            cmap[cindex] = col
            cassigned[cindex] = 1

    start = 0.0
    stop = 1.0
    numberleft = np.sum(cassigned==0)
    cm_subsection = np.linspace(start, stop, numberleft)
    for i, rind in enumerate(np.where(cassigned==0)[0]):
        jrow = cm.jet(cm_subsection[i])
        jrownew = increasesat(jrow)
        cmap[rind] = jrownew

    return np.array(cmap)


def preparecolormap(samecol, ncluster):
    if samecol:  # same color for all groups
        cmap, norm = makecolormap.makecolormap(2)
    else:
        cmap, norm = makecolormap.makecolormap(ncluster)

    if 1:
        # permute colormap so that adjacent constellations in list have different colors
        R = np.random.RandomState(10)
        cidperm = np.insert(R.permutation(ncluster-1)+1, 0, 0)
    else:
        cidperm = np.arange(ncluster)


    ns = [norm(i) for i in cidperm]
    cmap = cmap(ns)

    return cmap

# CLUSTERS: list of star indices
# EDGES: subgraphs showing edges (optional)
# SG: stargraph including all background stars

def plotmodelsoln(clusters, edges, sg, prefix='', show = 1, edgecol = 1, showlabs=0, samecol=0, title='', magthresh=4.5, edgecount= 0, xlim=(-np.pi-0.04, np.pi + 0.04), ylim=(-0.5*np.pi, 0.03 + 0.5*np.pi), figsize=(18.0, 9), alphaval=0.3, ethick = 3):
    # showlabs = 1 to annotate plot with starname , showlabs=0, samecol=0s
    # samecol = 1 to use same color for all groups

    # put small clusters last
    if clusters:
        cepairs = list(zip(clusters, edges))
        cepairs = sorted(cepairs, key = lambda x: len(x[0]))
        cepairs = cepairs[0:1] + cepairs[-1:0:-1]
        clusters, edges= zip(*cepairs)

    ncluster = len(clusters)
    cols = preparecolormap_aligned(sg, clusters, samecol, ncluster)
    #cols = preparecolormap(samecol, ncluster)
    # grey for pairs
    if 0:
        pairs = np.where([True  if len(c)==2 else False for c in clusters])[0]
        cols[pairs] = cols[0]

    # add color information to sg
    annotate_nodes(sg, clusters, cols)
    edgecols = np.copy(cols)

    # alpha value
    edgecols[:,3] = alphaval
    annotate_edges(sg, edges, clusters, edgecols, edgecol)

    # plot stargraph
    plt.rcParams['figure.figsize'] = figsize
    plt.close()
    plt.figure(1)
    ax = plt.subplot(111, projection="rectilinear")

    plot_graph(ax, sg, magthresh=magthresh, edgecount=edgecount, showlabs=showlabs, ethick = ethick)
   
    #ax.set_xlim(-np.pi, np.pi)
    #mindec = -0.5 * np.pi # now we're considering southern stars
    #ax.set_ylim(mindec, 0.5*np.pi)
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])

    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])

    plt.tight_layout()

    if prefix and re.search('videoframes', prefix):
        plt.savefig(prefix+ '.png')
    elif prefix:
        plt.savefig(prefix+ '.pdf', bbox_inches="tight")
        plt.savefig(prefix+ '.svg', format="svg", bbox_inches="tight")

    if show:
        plt.show(block=False)


def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)


def plotfullwithexamples(clusters, edges, sg, flip = 0, prefix='', show = 1, edgecol = 1, showlabs=0, samecol=0,  title='', magthresh=4.5, edgecount= 0, xlim=(-np.pi-0.04, np.pi + 0.04), ylim=(-0.5*np.pi, 0.03 + 0.5*np.pi), figsize=(18.0, 9.0)):

    plt.rcParams['figure.figsize'] = figsize
    plt.close()
    fig = plt.figure(constrained_layout=True)
    gs = gridspec.GridSpec(7, 10, figure=fig)
    gs.update(left=0.0, right=1.0, wspace=0.01, hspace=0.01)
    if flip:
        ax1 = fig.add_subplot(gs[6, 0])
        ax2 = fig.add_subplot(gs[6, 1])
        ax3 = fig.add_subplot(gs[6, 2])
        ax4 = fig.add_subplot(gs[6, 3])
        ax5 = fig.add_subplot(gs[6, 4])
        ax6 = fig.add_subplot(gs[6, 5])
        ax7 = fig.add_subplot(gs[6, 6])
        ax8 = fig.add_subplot(gs[6, 7])
        ax9 = fig.add_subplot(gs[6, 8])
        ax10 = fig.add_subplot(gs[6, 9])
        ax11 = fig.add_subplot(gs[0:6, :])
    else:
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[0, 2])
        ax4 = fig.add_subplot(gs[0, 3])
        ax5 = fig.add_subplot(gs[0, 4])
        ax6 = fig.add_subplot(gs[0, 5])
        ax7 = fig.add_subplot(gs[0, 6])
        ax8 = fig.add_subplot(gs[0, 7])
        ax9 = fig.add_subplot(gs[0, 8])
        ax10 = fig.add_subplot(gs[0, 9])
        ax11 = fig.add_subplot(gs[1:, :])
    plot_graph_sub(ax1, sg, '1. Pleiades', (-2.52, -2.27), (0.41, 0.435), magthresh, edgecount, showlabs)
    plot_graph_sub(ax2, sg, '2. Orion', (-3.02, -2.62), (-0.2, 0.2), magthresh, edgecount, showlabs)
    plot_graph_sub(ax3, sg, '3. Hyades', (-2.65, -2.5), (0.23, 0.38), magthresh, edgecount, showlabs)
    plot_graph_sub(ax4, sg, '4. Big Dipper', (1.1, 2.1), (0.5, 1.5), magthresh, edgecount, showlabs)
    plot_graph_sub(ax5, sg, '5. Southern Cross', (1.47, 1.71), (-1.17, -0.93), magthresh, edgecount, showlabs)
    plot_graph_sub(ax6, sg, '6. Corona Borealis', (0.67, 0.89), (0.39, 0.61), magthresh, edgecount, showlabs)
    plot_graph_sub(ax7, sg, '7. Castor & Pollux', (2.82, 2.92), (0.47, 0.57), magthresh, edgecount, showlabs)
    plot_graph_sub(ax8, sg, '8. Cassiopeia', (-2.0,-1.4), (0.75, 1.35), magthresh, edgecount, showlabs)
    plot_graph_sub(ax9, sg, '9. Delphinus', (-0.59,-0.49), (0.19, 0.29), magthresh, edgecount, showlabs)
    plot_graph_sub(ax10, sg, '10. Head of Aries',(-1.99, -1.89), (0.32, 0.42), magthresh, edgecount, showlabs)
    plot_graph_sub(ax11, sg, '', xlim, ylim, magthresh=magthresh, edgecount=edgecount, showlabs=showlabs)

    if re.search('consensus_edgethick_egs', prefix):
        ax11.text( -2.32,  0.425, '1')
        ax11.text( -3.07,  0.12, '2')
        ax11.text( -2.71,  0.3, '3')
        ax11.text( 2.05,  1.0, '4')
        ax11.text( 1.71, -1.01, '5')
        ax11.text( 0.71, 0.53, '6')
        ax11.text( 2.8, 0.56, '7')
        ax11.text(-1.7, 1.12, '8')
        ax11.text(-0.54, 0.33, '9')
        ax11.text(-1.880, 0.363, '10')
        ax11.text( 1.188, -1.18, '11')
        ax11.text(-0.28, 0.21, '12')
        ax11.text( 0.988, 1.35, '13')
        ax11.text( 0.72, -0.395, '14')
        ax11.text( 0.318, -0.638, '15')
        ax11.text( 2.07, 0.346, '16')
        ax11.text( 1.69, -0.30, '17')
        ax11.text( -0.58, 0.87, '18')
        ax11.text( -0.00, 0.74, '19')
        ax11.text( -1.32, 0.522 , '20')
        ax11.text( -0.10, -0.67, '21')
        ax11.text( 0.064, 0.889, '22')
        ax11.text( 0.0, -0.39, '23')

    if prefix and re.search('videoframes', prefix):
        ax11.set_title(title)
        plt.savefig(prefix+ '.png')
    elif prefix:
        plt.savefig(prefix+ '.pdf', bbox_inches="tight")
        plt.savefig(prefix+ '.svg', format="svg", bbox_inches="tight")

    if show:
        plt.show(block=False)


def plot_graph_sub(ax, sg, title, xlim, ylim, magthresh = 4.0, edgecount= 0, showlabs=0):
    plot_graph(ax, sg, magthresh=magthresh, edgecount=edgecount, showlabs=showlabs)
    if title:
        ax.set_title(title)
    #plt.tight_layout()
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    ax.set_xticks([])
    ax.set_yticks([])


def plotdemo(clusters, edges, sg,  prefix='', show = 1, edgecol = 1, showlabs=0, samecol=0,  title='', magthresh=4.5, edgecount= 0, xlim=(-np.pi-0.04, np.pi + 0.04), ylim=(-0.5*np.pi, 0.03 + 0.5*np.pi), figsize=(18.0, 9), alphaval=0.3, edgehue='gray'):

    ncluster = len(clusters)
    cols = preparecolormap_aligned(sg, clusters, samecol, ncluster)
    # add color information to sg
    annotate_nodes(sg, clusters, cols)
    edgecols = np.copy(cols)
    edgecols[0] = colors.to_rgba(edgehue)
    # alpha value
    edgecols[:,3] = alphaval
    annotate_edges(sg, edges, clusters, edgecols, edgecol)

    plt.rcParams['figure.figsize'] = figsize
    plt.close()
    fig = plt.figure(constrained_layout=True)
    gs = gridspec.GridSpec(6, 10, figure=fig)
    gs.update(left=0.0, right=1.0, wspace=0.01, hspace=0.01)
    ax2 = fig.add_subplot(gs[0, 1:3])
    plot_graph_sub(ax2, sg, '', (0.96, 1.76), (-1.23, -0.83), magthresh, edgecount, showlabs)

    # Save just the portion _inside_ the second axis's boundaries
    extent = full_extent(ax2).transformed(fig.dpi_scale_trans.inverted())

    if prefix:
        plt.savefig(prefix+ '.svg', format="svg", bbox_inches=extent)

    if show:
        plt.show(block=False)


def full_extent(ax, pad=0.0):
    """Get the full extent of an axes, including axes labels, tick labels, and
    titles."""
    # For text objects, we need to draw the figure first, otherwise the extents
    # are undefined.
    ax.figure.canvas.draw()
    items = ax.get_xticklabels() + ax.get_yticklabels()
#    items += [ax, ax.title, ax.xaxis.label, ax.yaxis.label]
    items += [ax, ax.title]
    bbox = Bbox.union([item.get_window_extent() for item in items])

    return bbox.expanded(1.0 + pad, 1.0 + pad)

