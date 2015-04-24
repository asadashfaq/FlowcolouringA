from __future__ import division
import sys
import copy
import numpy as np
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from EUgrid import EU_Nodes_usage
from functions import link_label, binMaker
from aurespf.tools import get_q, AtoKh_old
from link_namer import link_namer
import networkx as nx

"""
Figures for a paper on flowtracing with unconstrained synchronised flow with
gamma = 1, alpha = .7

These plotting functions are modifications of some from plotting.py and also
from Magnus': /AUESG/constrainedDCflow/figutils.py
"""

if len(sys.argv)<2:
    raise Exception('Not enough inputs!')
else:
    figNum = str(sys.argv[1:])

outPath = './figures/Paper figures feb 15/'

# Node indices and names sorted after descending mean load
loadOrder = [18, 4, 7, 22, 24, 20, 8, 5, 2, 6, 1, 15, 0, 10, 14,
             9, 11, 12, 16, 21, 17, 19, 3, 26, 13, 29, 27, 23, 28, 25]

loadNames = np.array(['DE', 'FR', 'GB', 'IT', 'ES', 'SE', 'PL', 'NO', 'NL',
                      'BE', 'FI', 'CZ', 'AT', 'GR', 'RO', 'BG', 'PT', 'CH',
                      'HU', 'DK', 'RS', 'IE', 'BA', 'SK', 'HR', 'LT', 'EE',
                      'SI', 'LV', 'LU'], dtype='|S4')

def make_europe_graph(link_weights, node_weights, t, figfilename='/fig 1/network', savepath=outPath, title=None):
    G = nx.Graph()

    dcolwidth = (2*3.425+0.236)
    colwidth = (3.425)
    blue = '#134b7c'

    all_countries = ['AUT', 'FIN', 'NLD', 'BIH', 'FRA', 'NOR', 'BEL','GBR', \
                     'POL', 'BGR', 'GRC', 'PRT', 'CHE', 'HRV', 'ROU', 'CZE',\
                     'HUN', 'SRB', 'DEU', 'IRL', 'SWE', 'DNK', 'ITA', 'SVN',\
                     'ESP', 'LUX', 'SVK', 'EST', 'LVA', 'LTU']

    ISO3ISO2dict = {'AUT':'AT', 'FIN':'FI', 'NLD':'NL', 'BIH':'BA',\
            'FRA':'FR', 'NOR':'NO', 'BEL':'BE','GBR':'GB', 'POL':'PL',\
            'BGR':'BG', 'GRC':'GR', 'PRT':'PT', 'CHE':'CH', 'HRV':'HR',\
            'ROU':'RO', 'CZE':'CZ', 'HUN':'HU', 'SRB':'RS', 'DEU':'DE',\
            'IRL':'IE', 'SWE':'SE', 'DNK':'DK', 'ITA':'IT', 'SVN':'SI',\
            'ESP':'ES', 'LUX':'LU', 'SVK':'SK', 'EST':'EE', 'LVA':'LV',\
            'LTU':'LT'}

    all_countries_ISO2 = [ISO3ISO2dict[c] for c in all_countries]

    all_links = ['AUT to CHE',
                 'AUT to CZE',
                 'AUT to HUN',
                 'AUT to DEU',
                 'AUT to ITA',
                 'AUT to SVN',
                 'FIN to SWE',
                 'FIN to EST',
                 'NLD to NOR',
                 'NLD to BEL',
                 'NLD to GBR',
                 'NLD to DEU',
                 'BIH to HRV',
                 'BIH to SRB',
                 'FRA to BEL',
                 'FRA to GBR',
                 'FRA to CHE',
                 'FRA to DEU',
                 'FRA to ITA',
                 'FRA to ESP',
                 'NOR to SWE',
                 'NOR to DNK',
                 'GBR to IRL',
                 'POL to CZE',
                 'POL to DEU',
                 'POL to SWE',
                 'POL to SVK',
                 'BGR to GRC',
                 'BGR to ROU',
                 'BGR to SRB',
                 'GRC to ITA',
                 'PRT to ESP',
                 'CHE to DEU',
                 'CHE to ITA',
                 'HRV to HUN',
                 'HRV to SRB',
                 'HRV to SVN',
                 'ROU to HUN',
                 'ROU to SRB',
                 'CZE to DEU',
                 'CZE to SVK',
                 'HUN to SRB',
                 'HUN to SVK',
                 'DEU to SWE',
                 'DEU to DNK',
                 'DEU to LUX',
                 'SWE to DNK',
                 'ITA to SVN',
                 'EST to LVA',
                 'LVA to LTU']

    nodelist = all_countries_ISO2
    all_links_ISO2 = []
    for l in all_links:
        l_ISO2 = copy.deepcopy(l)
        for c in all_countries:
            if c in l_ISO2:
                l_ISO2 = l_ISO2.replace(c, ISO3ISO2dict[c])
        all_links_ISO2.append(l_ISO2)

    linklist = [[link[0:2], link[-2:], all_links_ISO2.index(link)] \
                for link in all_links_ISO2]

    #print linklist
    for n in nodelist:
        G.add_node(n)

    for l in linklist:
        G.add_edge(l[0], l[1], weight = link_weights[l[2]])

    pos={}

    pos['AT']=[0.55,0.45]
    pos['FI']=[.95,1.1]
    pos['NL']=[0.40,0.85]
    pos['BA']=[0.65,0.15]
    pos['FR']=[0.15,0.60]
    pos['NO']=[0.5,1.1]
    pos['BE']=[0.275,0.775]
    pos['GB']=[0.10,1.05]
    pos['PL']=[0.75,0.8]
    pos['BG']=[0.9,0.0]
    pos['GR']=[0.7,0.0]
    pos['PT']=[0.0,0.15]
    pos['CH']=[0.4,0.45]
    pos['HR']=[0.75,0.3]
    pos['RO']=[1.0,0.15]
    pos['CZ']=[0.75,0.60]
    pos['HU']=[1.0,0.45]
    pos['RS']=[0.85,0.15]
    pos['DE']=[0.45,0.7]
    pos['IE']=[0.0,0.95]
    pos['SE']=[0.75,1.0]
    pos['DK']=[0.5,0.875]
    pos['IT']=[0.4,0.2]
    pos['SI']=[0.55,0.3]
    pos['ES']=[0.15,0.35]
    pos['LU']=[0.325,0.575]
    pos['SK']=[0.90,0.55]
    pos['EE']=[1.03,0.985]
    pos['LV']=[0.99,0.85]
    pos['LT']=[0.925,0.72]


    redgreendict = {'red':[(0.0, .3, .3), (0.2, .7, .7), (0.5, 1.0, 1.0) ,(.8, 0.0, 0.0) ,(1.0, 0.0, 0.0)],
                    'green':[(0.0, 0.0, 0.0), (0.2, 0.0, 0.0), (0.5, 1.0, 1.0),(.8, .7, .7), (1.0, .3, .3)],
                    'blue':[(0.0, 0.0, 0.0), (0.2, 0.0, 0.0), (0.5, 1.0, 1.0), (.8, 0.0, 0.0), (1.0, 0.0, 0.0)]}


    cmap = LinearSegmentedColormap('redgreen', redgreendict, 1000)

    fig = plt.figure(dpi=400, figsize=(0.85*dcolwidth,0.85*dcolwidth*0.8))
    ax1 = fig.add_axes([0.05, 0.08, 0.9, 0.08])
    ax2 = fig.add_axes([-0.05, 0.15, 1.1, 0.95])

    #span = 2*np.max(np.abs(node_weights))
    #norm_node_weights = [w/span+0.5 for w in node_weights]
    # now the scale shows the actual phi normalized to unit length
    # the offset and factor of 2 is because the colormap is defined
    # from 0 to 1
    phi_length = np.sqrt(np.sum(node_weights**2))
    norm_node_weights = [w/(2*phi_length) + 0.5 for w in node_weights]

    # node_colors = [cmap(w) for w in norm_node_weights]
    node_colors = [cmap(w) for w in node_weights]
    nx.draw_networkx_nodes(G, pos, node_size=400, nodelist=nodelist,\
            node_color=node_colors)
    nx.draw_networkx_labels(G, pos, font_size=10.8, font_color='k', \
                            font_family='sans-serif')

    maxflow = np.max(np.abs(link_weights))
    no_arrow_limit = 0.05*maxflow
    arrow_length = 0.05
    for l in linklist:
        if -no_arrow_limit<link_weights[l[2]]<no_arrow_limit: continue
        if link_weights[l[2]]>=no_arrow_limit:
            x0 = pos[l[0]][0]
            y0 = pos[l[0]][1]
            x1 = pos[l[1]][0]
            y1 = pos[l[1]][1]
        if link_weights[l[2]]<=-no_arrow_limit:
            x1 = pos[l[0]][0]
            y1 = pos[l[0]][1]
            x0 = pos[l[1]][0]
            y0 = pos[l[1]][1]
        dist = np.sqrt((x0-x1)**2 + (y0-y1)**2)
        plt.arrow(x0, y0, 0.4*(x1-x0), 0.4*(y1-y0), fc=blue, ec=blue,\
                    head_width= 1.2*arrow_length*(abs(link_weights[l[2]]))/maxflow,\
                    head_length = arrow_length)


    nx.draw_networkx_edges(G, pos, edgelist = G.edges(), edge_color=blue)

    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap, orientation='horizontal')
    cb1.set_ticks([0, 0.5, 1])
    cb1.set_ticklabels(['-1', '0', '1'])
    #ax1.set_xlabel(r'$\Phi_n$' + ' [normalized]')
    ax1.set_xlabel(r'$P_n / \left\langleL_n\right\rangle$')
    ax1.xaxis.set_label_position('top')
    ax1.set_xticks('none')
    ax2.axis('off')
    if title!=None:
        fig.suptitle(title)
    fig.savefig(savepath+figfilename+'-'+str(t)+'.png')
    plt.close()
    return

def draw_static_network(N=None, F=None, tit="1", show_link_size=True, typ=0, mode=''):
    N = EU_Nodes_usage('square.npz')
    F = np.load('./results/square-flows.npy')

    colwidth = (3.425)
    rolor= ( (19)/255. , (75)/255. , (124)/255.)
    G=nx.Graph()
    nodelist=[]

    col = ['#d4d4ff',
           '#adadff',
           '#8585ff',
           '#5e5eff',
           '#3737ff',
           '#1010ff',
           '#0000d4',
           '#000099',
           '#000072',
           '#000037']

    for n in N:
        G.add_node(str(n.label))
        nodelist.append(str(n.label))

    K,h,ListF=AtoKh_old(N)

    quants =  np.load('./results/quantiles_square_70128.npy')
    for q in range(50):
        G.add_edge(ListF[q][0], ListF[q][1] , weight=quants[q])

    pos={}
    pos['AT']=[0.55,0.45]
    pos['FI']=[.95,1.1]
    pos['NL']=[0.40,0.85]
    pos['BA']=[0.65,0.15]
    pos['FR']=[0.15,0.60]
    pos['NO']=[0.5,1.1]
    pos['BE']=[0.275,0.775]
    pos['GB']=[0.10,1.05]
    pos['PL']=[0.75,0.8]
    pos['BG']=[0.9,0.0]
    pos['GR']=[0.7,0.0]
    pos['PT']=[0.0,0.15]
    pos['CH']=[0.4,0.45]
    pos['HR']=[0.75,0.3]
    pos['RO']=[1.0,0.15]
    pos['CZ']=[0.75,0.60]
    pos['HU']=[1.0,0.45]
    pos['RS']=[0.85,0.15]
    pos['DE']=[0.45,0.7]
    pos['IE']=[0.0,0.95]
    pos['SE']=[0.75,1.0]
    pos['DK']=[0.5,0.875]
    pos['IT']=[0.4,0.2]
    pos['SI']=[0.55,0.3]
    pos['ES']=[0.15,0.35]
    pos['LU']=[0.325,0.575]
    pos['SK']=[0.90,0.55]
    pos['EE']=[1.0,0.94]
    pos['LV']=[0.95,0.83]
    pos['LT']=[0.87,0.72]

    fig = plt.figure(dpi=400,figsize=(1.7*colwidth,1.7*colwidth*0.75))

    if show_link_size: ax1= fig.add_axes([-0.125,0.135,1.25,0.975]) #For displaying graph
    else: ax1= fig.add_axes([-0.125,-0.075,1.25,1.25])
    nx.draw_networkx_nodes(G,pos,node_size=600,nodelist=nodelist,node_color=rolor,facecolor=(1,1,1))
    if mode == 'old':
        e0=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=1000]
        e1=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1000 and d['weight']<=1500]
        e2=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1500 and d['weight']<=2000]
        e3=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>2000 and d['weight']<=4000]
        e4=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>4000 and d['weight']<=6000]
        e5=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>6000 and d['weight']<=8000]
        e6=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>8000 and d['weight']<=10000]
        e7=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>10000 and d['weight']<=15000]
        e8=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>15000 and d['weight']<=20000]
        e9=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>20000]
        # ax1.text(-0.05,1.05,"(a)",fontsize=12)
        nx.draw_networkx_edges(G,pos,edgelist=e0,width=1.0,edge_color='k',alpha=1.0,style='dotted')
        nx.draw_networkx_edges(G,pos,edgelist=e1,width=0.5,edge_color='k',alpha=0.7,style='dashed')
        nx.draw_networkx_edges(G,pos,edgelist=e2,width=1.0,edge_color='k',alpha=0.8,style='dashed')
        nx.draw_networkx_edges(G,pos,edgelist=e3,width=2.0,edge_color='k',alpha=.9,style='dashed')
        nx.draw_networkx_edges(G,pos,edgelist=e4,width=3.0,edge_color='k',alpha=1.0,style='dashed')
        nx.draw_networkx_edges(G,pos,edgelist=e5,width=3.0,edge_color='k',alpha=0.6)
        nx.draw_networkx_edges(G,pos,edgelist=e6,width=3.5,edge_color='k',alpha=.7)
        nx.draw_networkx_edges(G,pos,edgelist=e7,width=4.0,edge_color='k',alpha=.8)
        nx.draw_networkx_edges(G,pos,edgelist=e8,width=4.5,edge_color='k',alpha=0.9)
        nx.draw_networkx_edges(G,pos,edgelist=e9,width=5.0,edge_color='k',alpha=1.0)
        nx.draw_networkx_labels(G,pos,font_size=13,font_color='w',font_family='sans-serif')
    else:
        e0=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=1000]
        e1=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1000 and d['weight']<=1500]
        e2=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1500 and d['weight']<=2000]
        e3=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>2000 and d['weight']<=4000]
        e4=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>4000 and d['weight']<=6000]
        e5=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>6000 and d['weight']<=8000]
        e6=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>8000 and d['weight']<=10000]
        e7=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>10000 and d['weight']<=15000]
        e8=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>15000 and d['weight']<=20000]
        e9=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>20000]
        # ax1.text(-0.05,1.05,"(a)",fontsize=12)
        nx.draw_networkx_edges(G,pos,edgelist=e0,width=4.0,edge_color=col[0],alpha=1.0)
        nx.draw_networkx_edges(G,pos,edgelist=e1,width=4.0,edge_color=col[1],alpha=1.0)
        nx.draw_networkx_edges(G,pos,edgelist=e2,width=4.0,edge_color=col[2],alpha=1.0)
        nx.draw_networkx_edges(G,pos,edgelist=e3,width=4.0,edge_color=col[3],alpha=1.0)
        nx.draw_networkx_edges(G,pos,edgelist=e4,width=4.0,edge_color=col[4],alpha=1.0)
        nx.draw_networkx_edges(G,pos,edgelist=e5,width=4.0,edge_color=col[5],alpha=1.0)
        nx.draw_networkx_edges(G,pos,edgelist=e6,width=4.0,edge_color=col[6],alpha=1.0)
        nx.draw_networkx_edges(G,pos,edgelist=e7,width=4.0,edge_color=col[7],alpha=1.0)
        nx.draw_networkx_edges(G,pos,edgelist=e8,width=4.0,edge_color=col[8],alpha=1.0)
        nx.draw_networkx_edges(G,pos,edgelist=e9,width=4.0,edge_color=col[9],alpha=1.0)
        nx.draw_networkx_labels(G,pos,font_size=13,font_color='w',font_family='sans-serif')
    ax1.axis('off')

    if show_link_size:
        ax4= fig.add_axes([-0.075,0.075,1.5,.15]) #For displaying graph
        if mode == 'old':
            ax4.plot((0.03*1.05+0.025,0.09*1.05+0.025),(0.75,0.75),linewidth=1.0,color='k',alpha=1.0,linestyle='dotted')
            ax4.plot((0.09*1.05+0.025,0.15*1.05+0.025),(0.75,0.75),linewidth=0.5,color='k',alpha=0.4,linestyle='dashed')
            ax4.plot((0.15*1.05+0.025,0.21*1.05+0.025),(0.75,0.75),linewidth=1.0,color='k',alpha=0.6,linestyle='dashed')
            ax4.plot((0.21*1.05+0.025,0.27*1.05+0.025),(0.75,0.75),linewidth=2.0,color='k',alpha=0.8,linestyle='dashed')
            ax4.plot((0.27*1.05+0.025,0.33*1.05+0.025),(0.75,0.75),linewidth=3.0,color='k',alpha=1.0,linestyle='dashed')
            ax4.plot((0.33*1.05+0.0295,0.39*1.05+0.025),(0.75,0.75),linewidth=3.0,color='k',alpha=0.6)
            ax4.plot((0.39*1.05+0.03,0.45*1.05+0.025),(0.75,0.75),linewidth=3.5,color='k',alpha=0.7)
            ax4.plot((0.45*1.05+0.0305,0.51*1.05+0.025),(0.75,0.75),linewidth=4.0,color='k',alpha=0.8)
            ax4.plot((0.51*1.05+0.0310,0.57*1.05+0.025),(0.75,0.75),linewidth=4.5,color='k',alpha=0.9)
            ax4.plot((0.57*1.05+0.0315,0.63*1.05+0.025),(0.75,0.75),linewidth=5.0,color='k',alpha=1.0)
        else:
            ax4.plot((0.03*1.05+0.025,0.09*1.05+0.025),(0.75,0.75),linewidth=4.0,color=col[0],alpha=1.0)
            ax4.plot((0.09*1.05+0.025,0.15*1.05+0.025),(0.75,0.75),linewidth=4.0,color=col[1],alpha=1.0)
            ax4.plot((0.15*1.05+0.025,0.21*1.05+0.025),(0.75,0.75),linewidth=4.0,color=col[2],alpha=1.0)
            ax4.plot((0.21*1.05+0.025,0.27*1.05+0.025),(0.75,0.75),linewidth=4.0,color=col[3],alpha=1.0)
            ax4.plot((0.27*1.05+0.025,0.33*1.05+0.025),(0.75,0.75),linewidth=4.0,color=col[4],alpha=1.0)
            ax4.plot((0.33*1.05+0.0295,0.39*1.05+0.025),(0.75,0.75),linewidth=4.0,color=col[5],alpha=1.0)
            ax4.plot((0.39*1.05+0.03,0.45*1.05+0.025),(0.75,0.75),linewidth=4.0,color=col[6],alpha=1.0)
            ax4.plot((0.45*1.05+0.0305,0.51*1.05+0.025),(0.75,0.75),linewidth=4.0,color=col[7],alpha=1.0)
            ax4.plot((0.51*1.05+0.0310,0.57*1.05+0.025),(0.75,0.75),linewidth=4.0,color=col[8],alpha=1.0)
            ax4.plot((0.57*1.05+0.0315,0.63*1.05+0.025),(0.75,0.75),linewidth=4.0,color=col[9],alpha=1.0)
        ax4.text(0.06*1.05+0.01,0.5,"$\leq$ 1.0 GW",fontsize=9,rotation=-60)
        ax4.text(0.12*1.05+0.01,0.5,"$\leq$ 1.5 GW",fontsize=9,rotation=-60)
        ax4.text(0.18*1.05+0.01,0.5,"$\leq$ 2.0 GW",fontsize=9,rotation=-60)
        ax4.text(0.24*1.05+0.01,0.5,"$\leq$ 4.0 GW",fontsize=9,rotation=-60)
        ax4.text(0.30*1.05+0.01,0.5,"$\leq$ 6.0 GW",fontsize=9,rotation=-60)
        ax4.text(0.36*1.05+0.01,0.5,"$\leq$ 8.0 GW",fontsize=9,rotation=-60)
        ax4.text(0.42*1.05+0.01,0.5,"$\leq$ 10 GW",fontsize=9,rotation=-60)
        ax4.text(0.48*1.05+0.01,0.5,"$\leq$ 15 GW",fontsize=9,rotation=-60)
        ax4.text(0.54*1.05+0.01,0.5,"$\leq$ 20 GW",fontsize=9,rotation=-60)
        ax4.text(0.60*1.05+0.01,0.5,"$\leq$ 45 GW",fontsize=9,rotation=-60)
        ax4.axis([0.0,1.0,0.0,1.2])
        ax4.axis('off')
    if mode == 'old':
        plt.savefig(outPath+'fig 2/network-old.pdf')
    else:
        plt.savefig(outPath+'fig 2/network.pdf')

def drawnet_capacities(N=None,scheme='square',direction='combined'):
    """
    Make network figures of a node's usage of links for both import, export and
    combined. Adapted from drawnet() in aurespf.plotting
    """
    colwidth = (3.425)
    dcolwidth = (2*3.425+0.236)

    if N==None:
        N=EU_Nodes_usage()
    G=nx.Graph()
    nodelist=[]

    # Add nodes and labels to networkx object for plotting
    for n in N:
        G.add_node(str(n.label))
        nodelist.append(str(n.label))

    # LF is a list of links
    K,h,LF=AtoKh_old(N)

    for l in LF:
        G.add_edge(l[0], l[1], id=l[2])

    # Define position of nodes
    pos={}
    pos['AT']=[0.55,0.45]
    pos['FI']=[.95,1.1]
    pos['NL']=[0.40,0.85]
    pos['BA']=[0.65,0.15]
    pos['FR']=[0.15,0.60]
    pos['NO']=[0.5,1.1]
    pos['BE']=[0.275,0.775]
    pos['GB']=[0.10,1.05]
    pos['PL']=[0.75,0.8]
    pos['BG']=[0.9,0.0]
    pos['GR']=[0.7,0.0]
    pos['PT']=[0.0,0.15]
    pos['CH']=[0.4,0.45]
    pos['HR']=[0.75,0.3]
    pos['RO']=[1.0,0.15]
    pos['CZ']=[0.75,0.60]
    pos['HU']=[1.0,0.45]
    pos['RS']=[0.85,0.15]
    pos['DE']=[0.45,0.7]
    pos['IE']=[0.0,0.95]
    pos['SE']=[0.75,1.0]
    pos['DK']=[0.5,0.875]
    pos['IT']=[0.4,0.2]
    pos['SI']=[0.55,0.3]
    pos['ES']=[0.15,0.35]
    pos['LU']=[0.325,0.575]
    pos['SK']=[0.90,0.55]
    pos['EE']=[1.0,0.985]
    pos['LV']=[0.975,0.87]
    pos['LT']=[0.925,0.77]

    # Define color scale for links
    blueDict = {'red': ((0.0, 1.0, 1.0),(1.0, 0.0, 0.0)),
             'green': ((0.0, 1.0, 1.0),(1.0, 0.0, 0.0)),
             'blue': ((0.0, 1.0, 1.0),(1.0, 1.0, 1.0))}

    blueDict2 = {'red': ((0.0, 1.0, 1.0),(.2,0,0),(.5, 0.5, 0.5),(1.0, 1, 1)),
             'green': ((0.0, 1.0, 1.0),(.2,0,0),(.5, 0, 0),(1.0, 0.0, 0.0)),
             'blue': ((0.0, 1.0, 1.0),(.2,1,1),(.5, 0.5, 0.5),(1.0, 0, 0))}

    cmap = LinearSegmentedColormap('blue',blueDict2,1000)

    # Load usages for given scheme and direction
    quantiles = np.load('./results/quantiles_'+str(scheme)+'_70128.npy')

    # Calculate colors of links
    col = [(cmap(l)) for l in quantiles/max(quantiles)]

    # Create a new figure and plot network below
    fig = plt.figure(dpi=400,figsize=(0.85*dcolwidth,0.85*0.8*dcolwidth))

    # color bar in bottom of figure
    ax1 = fig.add_axes([0.05,0.08,0.9,.08])
    cbl = mpl.colorbar.ColorbarBase(ax1,cmap,orientation='horizontal')

    # Label just above color bar
    ax1.set_xlabel(r'$\mathcal{K}^T_{l}/\max(\mathcal{K}^T_l)$')
    ax1.xaxis.set_label_position('top')
    cbl.set_ticks(np.linspace(0,1,11))
    cbl.set_ticklabels(np.linspace(0,1,11))

    # Set color of nodes and draw
    node_c = ["#000000" for node in N]
    ax2 = fig.add_axes([-0.05,0.15,1.1,0.95])
    nx.draw_networkx_nodes(G,pos,node_size=500,nodelist=nodelist,node_color=node_c,facecolor=(1,1,1))

    # Draw links colored by usage of node n
    edges=[(u,v) for (u,v,d) in G.edges(data=True)]
    edge_id = [d['id'] for (u,v,d) in G.edges(data=True)]

    color_sort = []
    for i in range(len(col)):
        color_sort.append(col[edge_id[i]])
    nx.draw_networkx_edges(G,pos,edgelist=edges,width=3.5,edge_color=color_sort,alpha=1)

    # Draw country names
    nx.draw_networkx_labels(G,pos,font_size=12,font_color='w',font_family='sans-serif')
    ax2.axis('off')

    # Save figure
    plt.savefig(outPath+'fig 2/network-alternate.pdf')

def scatter_plotter(N, F, Fmax, usage, direction, mode):
    """
    Scatter plots of nodes' import/export usages of links saved to ./figures/.
    """
    nodes = [0, 3, 18, 24]
    links = usage.shape[0]
    for l in range(links):
        diag = []
        diagflow = []
        for n in nodes:
            names = ['diagonal', r'$99\%$ quantile', 'avg. usage', 'usage']
            plt.figure()
            ax = plt.subplot(111)
            linkflow = abs(F[l,:])
            qq = get_q(abs(F[l]),.99)
            if mode == 'old':
                usages = usage[l,n,:]/Fmax[l]
                plt.plot([0, Fmax[l]/qq], [0, 1], '-k',lw=1) # diagonal
                ax.set_ylabel(r'$C_{ln}(t)$')
            if mode == 'new':
                usages = usage[l,n,:]/linkflow
                names = names[1:]
                ax.set_ylabel(r'$C_{ln}(t)$')

            # scatter
            plt.scatter(linkflow/qq, usages, c='#000099', edgecolor='none', alpha=.2)

            # 99 quantile
            # plt.plot([1,1],[0,1],':k')

            # Plot bin avg usage
            F_vert = np.reshape(linkflow,(len(linkflow),1))
            exp_vert = np.reshape(usages,(len(usages),1))
            F_matrix = np.hstack([F_vert,exp_vert]) # [flow, usage]
            F_matrix[F_matrix[:,0].argsort()]
            H,bin_edges = binMaker(F_matrix, qq, lapse=70128)
            plt.plot(bin_edges/qq, H[:,1], '-', c='#aa0000', lw=2)

            label = link_label(l,N)
            # ax.set_title('Synchronised'+' '+str(direction)+' flows on link '+label)
            ax.set_xlabel(r'$|F_l(t)|/\mathcal{K}_l^T$')

            # Plot link name and direction in figure
            plt.text(.9, .95, label)
            plt.text(.9, .9, direction)

            # Shrink x-axis to make room for legend
            # box = ax.get_position()
            # ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
            # ax.legend((names),loc='center left', bbox_to_anchor=(1,0.5))#,title='Contributions')

            plt.axis([0, 1, 0, 1])
            if mode == 'old':
                plt.savefig(outPath+'fig 4/'+str(N[n].label)+'/'+str(l)+'-'+str(direction)+'.png', bbox_inches='tight')
            if mode == 'new':
                plt.savefig(outPath+'fig 4/'+str(N[n].label)+'/'+str(l)+'-'+str(direction)+'-'+mode+'.png', bbox_inches='tight')
            plt.close('all') # fixes memory leak.
    return

def drawnet_import(N, scheme='square', direction='import'):
    colwidth = (3.425)
    dcolwidth = (2*3.425+0.236)

    G=nx.Graph()
    nodelist=[]

    # Add nodes and labels to networkx object for plotting
    for n in N:
        G.add_node(str(n.label))
        nodelist.append(str(n.label))

    # LF is a list of links
    K,h,LF=AtoKh_old(N)

    for l in LF:
        G.add_edge(l[0], l[1], id=l[2])

    # Define position of nodes
    pos={}
    pos['AT']=[0.55,0.45]
    pos['FI']=[.95,1.1]
    pos['NL']=[0.40,0.85]
    pos['BA']=[0.65,0.15]
    pos['FR']=[0.15,0.60]
    pos['NO']=[0.5,1.1]
    pos['BE']=[0.275,0.775]
    pos['GB']=[0.10,1.05]
    pos['PL']=[0.75,0.8]
    pos['BG']=[0.9,0.0]
    pos['GR']=[0.7,0.0]
    pos['PT']=[0.0,0.15]
    pos['CH']=[0.4,0.45]
    pos['HR']=[0.75,0.3]
    pos['RO']=[1.0,0.15]
    pos['CZ']=[0.75,0.60]
    pos['HU']=[1.0,0.45]
    pos['RS']=[0.85,0.15]
    pos['DE']=[0.45,0.7]
    pos['IE']=[0.0,0.95]
    pos['SE']=[0.75,1.0]
    pos['DK']=[0.5,0.875]
    pos['IT']=[0.4,0.2]
    pos['SI']=[0.55,0.3]
    pos['ES']=[0.15,0.35]
    pos['LU']=[0.325,0.575]
    pos['SK']=[0.90,0.55]
    pos['EE']=[1.0,0.985]
    pos['LV']=[0.975,0.87]
    pos['LT']=[0.925,0.77]

    # Define color scale for links
    blueDict = {'red': ((0.0, 1.0, 1.0),(1.0, 0.0, 0.0)),
             'green': ((0.0, 1.0, 1.0),(1.0, 0.0, 0.0)),
             'blue': ((0.0, 1.0, 1.0),(1.0, 1.0, 1.0))}

    blueDict2 = {'red': ((0.0, 1.0, 1.0),(.15,0,0),(.4, 0.5, 0.5),(1.0, .7, .7)),
             'green': ((0.0, 1.0, 1.0),(.15,0,0),(.4, 0, 0),(1.0, 0.0, 0.0)),
             'blue': ((0.0, 1.0, 1.0),(.15,.8,.8),(.4, 0.5, 0.5),(1.0, 0, 0))}

    cmap = LinearSegmentedColormap('blue',blueDict2,1000)

    # color scale for nodes
    greenDict = {'red': ((0.0, 1.0, 1.0),(.2, 0.0, 0.0),(1, 0.0, 0.0)),
             'green': ((0.0, 1.0, 1.0),(.2, .7, .7),(1, .3, .3)),
             'blue': ((0.0, 1.0, 1.0),(.2, 0.0, 0.0),(1, 0.0, 0.0))}
    cmapNodes = LinearSegmentedColormap('green',greenDict,1000)

    # Load usages for given scheme and direction
    N_usages = np.load('./results/Node_contrib_'+scheme+'_'+direction+'_70128.npy')
    quantiles = np.load('./results/quantiles_'+str(scheme)+'_70128.npy')

    # Load power mixes
    pmim = np.load('./results/square_pm.npz', mmap_mode='r')['power_mix']
    # pmex = np.load('./results/square_pm.npz', mmap_mode='r')['power_mix_ex']
    pmim = np.sum(pmim, 2)

    # Pick a particular node
    for n in N:
        # Calculate colors of links
        N_usages[n.id] = (N_usages[n.id]/quantiles)*2
        col = [ (cmap(l)) for l in N_usages[n.id]]

        # Calculate colors of nodes
        pmim[n.id,n.id] = 0
        pm = (pmim[n.id]/sum(pmim[n.id]))*2
        node_col = [ (cmapNodes(nn)) for nn in pm]

        # Create a new figure and plot network below
        fig = plt.figure(dpi=400,figsize=(0.85*dcolwidth,0.85*0.8*dcolwidth))

        # color bar in bottom of figure
        ax1 = fig.add_axes([0.05,0.04,0.4,.08])
        cbl = mpl.colorbar.ColorbarBase(ax1,cmap,orientation='horizontal')
        ax2 = fig.add_axes([0.55,0.04,0.40,.08])
        cbl2 = mpl.colorbar.ColorbarBase(ax2,cmapNodes,orientation='horizontal')

        # Label just above color bar
        if scheme == 'linear':
            xlabel = 'Most localised'
        elif scheme == 'square':
            xlabel = 'Synchronised'

        ax1.set_xlabel(r'$\mathcal{K}^T_{ln}/\mathcal{K}^T_l$')
        ax1.xaxis.set_label_position('top')
        cbl.set_ticks(np.linspace(0,1,6))
        cbl.set_ticklabels(['0','0.1','0.2','0.3','0.4','0.5'])

        ax2.set_xlabel(r'$\mathcal{I}_{n\leftarrow m}$')
        ax2.xaxis.set_label_position('top')
        cbl2.set_ticks(np.linspace(0,1,6))
        cbl2.set_ticklabels(['0','0.1','0.2','0.3','0.4','0.5'])

        ax3 = fig.add_axes([-0.05,0.15,1.1,0.95])

        # Set color of nodes, highlight one and draw all
        # node_c = [ "#000000" for node in N]
        node_c = []
        for d in range(30):
            if n.id == d:
                node_c.append(cmap(.99))
            else:
                node_c.append(node_col[d])
        nx.draw_networkx_nodes(G,pos,node_size=500,nodelist=nodelist,node_color=node_c,facecolor=(1,1,1))

        # Draw links colored by usage of node n
        edges=[(u,v) for (u,v,d) in G.edges(data=True)]
        edge_id = [d['id'] for (u,v,d) in G.edges(data=True)]

        color_sort = []
        for i in range(len(col)):
            color_sort.append(col[edge_id[i]])
        nx.draw_networkx_edges(G,pos,edgelist=edges,width=3.5,edge_color=color_sort,alpha=1)

        # Draw country names
        nx.draw_networkx_labels(G,pos,font_size=12,font_color='w',font_family='sans-serif')
        ax3.axis('off')

        # Save figure
        plt.savefig(outPath+'fig 5/'+str(n.id)+'_'+str(direction)+".png")
        plt.close()


def drawnet_export(N, scheme='square', direction='export'):
    colwidth = (3.425)
    dcolwidth = (2*3.425+0.236)

    G=nx.Graph()
    nodelist=[]

    # Add nodes and labels to networkx object for plotting
    for n in N:
        G.add_node(str(n.label))
        nodelist.append(str(n.label))

    # LF is a list of links
    K,h,LF=AtoKh_old(N)

    for l in LF:
        G.add_edge(l[0], l[1], id=l[2])

    # Define position of nodes
    pos={}
    pos['AT']=[0.55,0.45]
    pos['FI']=[.95,1.1]
    pos['NL']=[0.40,0.85]
    pos['BA']=[0.65,0.15]
    pos['FR']=[0.15,0.60]
    pos['NO']=[0.5,1.1]
    pos['BE']=[0.275,0.775]
    pos['GB']=[0.10,1.05]
    pos['PL']=[0.75,0.8]
    pos['BG']=[0.9,0.0]
    pos['GR']=[0.7,0.0]
    pos['PT']=[0.0,0.15]
    pos['CH']=[0.4,0.45]
    pos['HR']=[0.75,0.3]
    pos['RO']=[1.0,0.15]
    pos['CZ']=[0.75,0.60]
    pos['HU']=[1.0,0.45]
    pos['RS']=[0.85,0.15]
    pos['DE']=[0.45,0.7]
    pos['IE']=[0.0,0.95]
    pos['SE']=[0.75,1.0]
    pos['DK']=[0.5,0.875]
    pos['IT']=[0.4,0.2]
    pos['SI']=[0.55,0.3]
    pos['ES']=[0.15,0.35]
    pos['LU']=[0.325,0.575]
    pos['SK']=[0.90,0.55]
    pos['EE']=[1.0,0.985]
    pos['LV']=[0.975,0.87]
    pos['LT']=[0.925,0.77]

    # Define color scale for links
    blueDict = {'red': ((0.0, 1.0, 1.0),(1.0, 0.0, 0.0)),
             'green': ((0.0, 1.0, 1.0),(1.0, 0.0, 0.0)),
             'blue': ((0.0, 1.0, 1.0),(1.0, 1.0, 1.0))}

    blueDict2 = {'red': ((0.0, 1.0, 1.0),(.15,0,0),(.4, 0.5, 0.5),(1.0, .7, .7)),
             'green': ((0.0, 1.0, 1.0),(.15,0,0),(.4, 0, 0),(1.0, 0.0, 0.0)),
             'blue': ((0.0, 1.0, 1.0),(.15,.8,.8),(.4, 0.5, 0.5),(1.0, 0, 0))}

    cmap = LinearSegmentedColormap('blue',blueDict2,1000)

    # color scale for nodes
    greenDict = {'red': ((0.0, 1.0, 1.0),(.2, 0.0, 0.0),(1, 0.0, 0.0)),
             'green': ((0.0, 1.0, 1.0),(.2, .7, .7),(1, .3, .3)),
             'blue': ((0.0, 1.0, 1.0),(.2, 0.0, 0.0),(1, 0.0, 0.0))}
    cmapGreen = LinearSegmentedColormap('green',greenDict,1000)

    # color scale for nodes
    redDict = {'red': ((0.0, 1.0, 1.0),(.15, .7, .7),(1, .4, .4)),
             'green': ((0.0, 1.0, 1.0),(.15, 0.0, 0.0),(1, 0.0, 0.0)),
             'blue': ((0.0, 1.0, 1.0),(.15, 0.0, 0.0),(1, 0.0, 0.0))}
    cmapNodes = LinearSegmentedColormap('red',redDict,1000)

    # Load usages for given scheme and direction
    N_usages = np.load('./results/Node_contrib_'+scheme+'_'+direction+'_70128.npy')
    quantiles = np.load('./results/quantiles_'+str(scheme)+'_70128.npy')

    # Load power mixes
    pmex = np.load('./results/square_pm.npz', mmap_mode='r')['power_mix_ex']
    pmex = np.sum(pmex, 2)

    # Pick a particular node
    for n in N:
        # Calculate colors of links
        N_usages[n.id] = (N_usages[n.id]/quantiles)*2
        col = [ (cmap(l)) for l in N_usages[n.id]]

        # Calculate colors of nodes
        pmex[n.id,n.id] = 0
        pm = (pmex[n.id]/sum(pmex[n.id]))*2
        node_col = [ (cmapNodes(nn)) for nn in pm]

        # Create a new figure and plot network below
        fig = plt.figure(dpi=400,figsize=(0.85*dcolwidth,0.85*0.8*dcolwidth))

        # color bar in bottom of figure
        ax1 = fig.add_axes([0.05,0.04,0.4,.08])
        cbl = mpl.colorbar.ColorbarBase(ax1,cmap,orientation='horizontal')
        ax2 = fig.add_axes([0.55,0.04,0.40,.08])
        cbl2 = mpl.colorbar.ColorbarBase(ax2,cmapNodes,orientation='horizontal')

        # Label just above color bar
        if scheme == 'linear':
            xlabel = 'Most localised'
        elif scheme == 'square':
            xlabel = 'Synchronised'

        ax1.set_xlabel(r'$\mathcal{K}^T_{ln}/\mathcal{K}^T_l$')
        ax1.xaxis.set_label_position('top')
        cbl.set_ticks(np.linspace(0,1,6))
        cbl.set_ticklabels(['0','0.1','0.2','0.3','0.4','0.5'])

        ax2.set_xlabel(r'$\mathcal{E}_{n\rightarrow m}$')
        ax2.xaxis.set_label_position('top')
        cbl2.set_ticks(np.linspace(0,1,6))
        cbl2.set_ticklabels(['0','0.1','0.2','0.3','0.4','0.5'])

        ax3 = fig.add_axes([-0.05,0.15,1.1,0.95])

        # Set color of nodes, highlight one and draw all
        # node_c = [ "#000000" for node in N]
        node_c = []
        for d in range(30):
            if n.id == d:
                node_c.append(cmapGreen(.15))
            else:
                node_c.append(node_col[d])
        nx.draw_networkx_nodes(G,pos,node_size=500,nodelist=nodelist,node_color=node_c,facecolor=(1,1,1))

        # Draw links colored by usage of node n
        edges=[(u,v) for (u,v,d) in G.edges(data=True)]
        edge_id = [d['id'] for (u,v,d) in G.edges(data=True)]

        color_sort = []
        for i in range(len(col)):
            color_sort.append(col[edge_id[i]])
        nx.draw_networkx_edges(G,pos,edgelist=edges,width=3.5,edge_color=color_sort,alpha=1)

        # Draw country names
        nx.draw_networkx_labels(G,pos,font_size=12,font_color='w',font_family='sans-serif')
        ax3.axis('off')

        # Save figure
        plt.savefig(outPath+'fig 5/'+str(n.id)+'_'+str(direction)+".png")
        plt.close()


def usageMesh(direction):
    scheme = 'square'
    N_usages = np.load('./results/Node_contrib_' + scheme + '_' + direction + '_70128.npy')
    quantiles = np.load('./results/quantiles_' + str(scheme) + '_70128.npy')
    usages = N_usages / quantiles
    usages = usages[loadOrder]
    plt.figure(figsize=(10, 11))
    ax = plt.subplot(1, 1, 1)
    plt.pcolormesh(usages.T, cmap='Blues')
    plt.colorbar().set_label(label=r'$\mathcal{K}_{ln}^T / \mathcal{K}_{l}^T$', size=13)
    ax.set_xticks(np.linspace(1, 30, 30))
    ax.set_xticklabels(loadNames, rotation=60, ha="right", va="top", fontsize=10.5)
    plt.grid(True)
    #ax.xaxis.grid(False)
    ax.xaxis.set_tick_params(width=0)
    ax.set_yticks(np.linspace(1, 50, 50))
    ax.set_yticklabels(lnames, ha="right", va="top", fontsize=10.5)
    #ax.yaxis.grid(False)
    ax.yaxis.set_tick_params(width=0)
    plt.savefig(outPath + 'table 3/usages-' + direction + '.pdf', bbox_inches='tight')

if '1' in figNum:
    print 'Making figure 1'
    N = np.load('./results/square.npz', mmap_mode='r')
    F = np.load('./results/square-flows.npy')
    ranges = [range(24), range(70104,70128)]
    for r in ranges:
        for t in r:
            node_weights = N['mismatch'][:,t] + N['balancing'][:,t]
            means = N['mean']
            node_weights =  node_weights / means
            # Scaling to the colorbar so [-1, 1] becomes [0, 1]
            node_weights += 1
            node_weights /= 2
            make_europe_graph(F[:,t], node_weights, t)

if '2' in figNum:
    print 'Making figure 2'
    draw_static_network(mode='old')
    draw_static_network(mode='new')
    drawnet_capacities()

if '3' in figNum:
    print 'Making table 3'
    directions = ['import', 'export', 'combined']
    lnames = np.array(link_namer())
    p = Pool(len(directions))
    p.map(usageMesh, directions)

if '4' in figNum:
    print 'Making figure 4'
    N = EU_Nodes_usage('square.npz')
    F = np.load('./results/square-flows.npy')
    Fmax = np.max(np.abs(F),1)

    export_usage = np.load('./linkcolouring/old_square_copper_link_mix_export_all_alpha=same.npy')
    scatter_plotter(N, F, Fmax, export_usage, 'export', mode='old')
    scatter_plotter(N, F, Fmax, export_usage, 'export', mode='new')
    export_usage = [] # frees roughly 800 MB RAM before loading import usage.

    import_usage = np.load('./linkcolouring/old_square_copper_link_mix_import_all_alpha=same.npy')
    scatter_plotter(N, F, Fmax, import_usage, 'import', mode='old')
    scatter_plotter(N, F, Fmax, import_usage, 'import', mode='new')
    import_usage = []

if '5' in figNum:
    print 'Making figure 5'
    N = EU_Nodes_usage('square.npz')
    drawnet_import(N)
    drawnet_export(N)

if '6' in figNum:
    print 'Making figure 6'
    N = np.load('./results/square.npz')

    sqr_mix_ex = np.load('./linkcolouring/old_square_copper_link_mix_export_all_alpha=same.npy')
    sqr_mix_im = np.load('./linkcolouring/old_square_copper_link_mix_import_all_alpha=same.npy')
    sqr_net = np.mean(np.sum(sqr_mix_ex,0),1) - np.mean(np.sum(sqr_mix_im,0),1)

    order = N['mean'].argsort()[::-1]

    plt.figure(figsize=(10, 4), facecolor='w', edgecolor='k')
    ax = plt.subplot(111)
    plt.plot([0,30], [0,0], '-k')
    plt.bar(range(30), sqr_net[order]/N['mean'][order], edgecolor='none', color='#000099')
    ax.set_xticks(np.linspace(.8,29.8,30))
    ax.set_xticklabels(N['label'][order], rotation=60, ha="right", va="top", fontsize=10.5)
    ax.xaxis.grid(False)
    ax.xaxis.set_tick_params(width=0)
    ax.set_ylabel('Net exports')
    plt.savefig(outPath+'fig 6/net_exports.png', bbox_inches='tight')
