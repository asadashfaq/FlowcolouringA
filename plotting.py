import sys
import numpy as np
from pylab import plt
from multiprocessing import Pool
from aurespf.tools import *
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx
import matplotlib as mpl
from EUgrid import EU_Nodes_usage, EU_Nodes_regions, EU_Nodes_superRegions
from link_namer import node_namer, link_dict
from scipy.stats import pearsonr
from functions import linkProportional, getLengths

"""
Script that makes network figures of a country's usage for import, export and
the combination

Call the script using only one of the following command line arguments:
- network:                  network figures with coloured links, only works for N=30
- total:                    barplots comparing total network usage for different export schemes
- total sensitivity:        same as above, but for N=8 and N=50 networks
- sensitivity compare:      compare different networks: N=30 & N=8 ; N=30 & N=53
- sensitivity compare all:  as just above for N=53, N=30, N=8

The last 4 can be followed by 'length' to include modelling of link lengths.
"""

if len(sys.argv)>1:
    task = str(sys.argv[1:])
else:
    raise Exception('Not enough inputs!')

schemes = ['linear' ,'square'] #'RND'
directions = ['import', 'export', 'combined']

if 'length' in task: length = True # modelling of link lengths
else: length = False

def simpleMerger(data, merge_dict):
    """
    Function to merge a list of data for regions to a list of data for
    countries. This function depends on the merge_dict which does not exist
    before the function regionMerger has been run.
    """
    merge_index = 0
    merged_data = np.zeros(len(merge_dict))
    for d in merge_dict:
        regions = merge_dict[d]
        for region in regions:
            merged_data[merge_index] += data[region]
        merge_index += 1
    return merged_data

def regionMerger(scheme, direction, lapse, md=None, node_ids=None):
    """
    Merge regional usage to compare with country's usage
    """
    N_countries = EU_Nodes_usage()
    N_regions = EU_Nodes_regions()
    
    # The optional parameter is intended for sorting of merged regions. If the
    # parameter is not set the ordering will be the default ordering as at is in
    # nodes objects.
    if node_ids == None:
        node_ids = np.arange(30)

    # Create a dictionary for merging. All entries are integers referring to IDs.
    # The format of the dictionary is: {country : [region1, region2, region3]}
    if md:
        merge_dict = md
    else:
        merge_dict = {}
        for node in node_ids:
            node = int(node)
            for region in N_regions:
                # Merge Northern Ireland with GB
                if ((region.label.tostring() == 'IE_N') and (node == 7)):
                    thisNode = 7 # GB
                    if thisNode in merge_dict:
                        merge_dict[thisNode].append(region.id)
                    else:
                        merge_dict.update({thisNode: [region.id]})
                elif region.label.tostring()[0:2] == N_countries[node].label.tostring():
                    if node in merge_dict:
                        merge_dict[node].append(region.id)
                    else:
                        merge_dict.update({node: [region.id]})

    if length: region_usages = np.load('./sensitivity/regions-Node_contrib_'+scheme+'_'+direction+'_length_'+str(lapse)+'.npy')
    else: region_usages = np.load('./sensitivity/regions-Node_contrib_'+scheme+'_'+direction+'_'+str(lapse)+'.npy')

    # Merge usages and mean load with the dictionary build above.
    merged_usages = np.zeros((len(merge_dict),region_usages.shape[1]))
    merged_mean_load = np.zeros(len(merge_dict))
    merge_index = 0
    for d in merge_dict:
        regions = merge_dict[d]
        for region in regions:
            merged_usages[merge_index] += region_usages[region]
            merged_mean_load[merge_index] += N_regions[region].mean
        merge_index += 1
    
    return merged_usages, merged_mean_load, merge_dict

def superRegionMerger(N, usages, countryNames, region_mean_load=None, md=None):
    """
    Merge regional or country usage to compare with super regions
    """
    aggregate_dict={\
    'UK':['IE','GB'],\
    'FR':['NL','BE','FR'],\
    'HI':['ES','PT'],\
    'IT':['CH','IT'],\
    'SG':['SI','HR','RS','BA','GR','BG','RO'],\
    'EE':['AT','CZ','SK','HU','PL'],\
    'DE':['DE','LU'], \
    'SC':['DK','SE','NO','FI','EE','LV','LV']}
    merge_dict = {}
    for i,SR in enumerate(['DE', 'EE', 'FR', 'HI', 'IT', 'SC', 'SG', 'UK']):
        for country in countryNames:
            if country in aggregate_dict[SR]:
                if i in merge_dict:
                    merge_dict[i].append(countryNames.index(country))
                else:
                    merge_dict.update({i: [countryNames.index(country)]})

    # Merge usages and mean load with the dictionary build above.
    merged_usages = np.zeros((len(merge_dict),usages.shape[1]))
    merged_mean_load = np.zeros(len(merge_dict))
    merge_index = 0
    for d in merge_dict:
        countries = merge_dict[d]
        for country in countries:
            merged_usages[merge_index] += usages[country]
            if region_mean_load != None:
                merged_mean_load[merge_index] += region_mean_load[country]
            else:
                merged_mean_load[merge_index] += N[country].mean
        merge_index += 1

    return merged_usages, merged_mean_load, merge_dict

def drawnet_usage(N=None,scheme='linear',direction='combined'):
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
    N_usages = np.load('./results/Node_contrib_'+scheme+'_'+direction+'_70128.npy')
    quantiles = np.load('./results/quantiles_'+str(scheme)+'_70128.npy')

    # Pick a particular node
    for n in N:
        # Calculate colors of links
        N_usages[n.id] = N_usages[n.id]/quantiles
        col = [ (cmap(l)) for l in N_usages[n.id]]
        
        # Create a new figure and plot network below
        fig = plt.figure(dpi=400,figsize=(0.85*dcolwidth,0.85*0.8*dcolwidth))
    
        # color bar in bottom of figure
        ax1 = fig.add_axes([0.05,0.08,0.9,.08])
        cbl = mpl.colorbar.ColorbarBase(ax1,cmap,orientation='vhorizontal')
        
        # Label just above color bar
        if scheme == 'linear':
            xlabel = 'Most localised'
        elif scheme == 'square':
            xlabel = 'Synchronised'
        else:
            xlabel = 'Market'
        ax1.set_xlabel(xlabel+' '+direction+r" usage $C_n/C^{\,99\%}$")
        ax1.xaxis.set_label_position('top') 
        
        ax2 = fig.add_axes([-0.05,0.15,1.1,0.95])
        
        # Set color of nodes, highlight one and draw all 
        node_c = [ "#000000" for node in N]
        node_c[n.id] = "#B30000"
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
        plt.savefig("./figures/network_figures/"+scheme+"/"+str(n.id)+'_'+str(direction)+".png")

def bars(scheme, verbose=None):
    """
    Figure to compare link proportional and usage proportional for a single
    scheme and put them in ./sensitivity/figures/scheme/
    """
    # Load data and results
    lengths = False
    if len(N) == 8:
        F = abs(np.load('./sensitivity/superRegions-'+scheme+'-flows.npy'))
        quantiles = np.load('./sensitivity/superRegions-quantiles_'+scheme+'_'+str(lapse)+'.npy')
        network = 'superRegions'
        nNodes = 8
        if length:
            lengths = 'superRegions'
    elif len(N) == 30:
        F = abs(np.load('./results/'+scheme+'-flows.npy'))
        quantiles = np.load('./results/quantiles_'+scheme+'_'+str(lapse)+'.npy')
        network = ''
        nNodes = 30
        if length:
            lengths = 'countries'
    elif len(N) == 53:
        F = abs(np.load('./sensitivity/regions-'+scheme+'-flows.npy'))
        quantiles = np.load('./sensitivity/regions-quantiles_'+scheme+'_'+str(lapse)+'.npy')
        network = 'regions'
        nNodes = 53
        if length:
            lengths = 'regions'
    else:
        raise Exception('Wrong network!')

    names = node_namer(N) # array of node labels
    links = range(len(F))
    nodes = np.linspace(0.5,2*nNodes-1.5,nNodes)
    nodes_shift = nodes+.5

    for direction in directions:
        if network == 'regions':
            if length: N_usages = np.load('./sensitivity/regions-Node_contrib_'+scheme+'_'+direction+'_length_'+str(lapse)+'.npy')
            else: N_usages = np.load('./sensitivity/regions-Node_contrib_'+scheme+'_'+direction+'_'+str(lapse)+'.npy')
        elif network == 'superRegions':
            if length: N_usages = np.load('./sensitivity/superRegions-Node_contrib_'+scheme+'_'+direction+'_length_'+str(lapse)+'.npy')
            else: N_usages = np.load('./sensitivity/superRegions-Node_contrib_'+scheme+'_'+direction+'_'+str(lapse)+'.npy')
        else:
            if length: N_usages = np.load('./results/Node_contrib_'+scheme+'_'+direction+'_length_'+str(lapse)+'.npy')
            else: N_usages = np.load('./results/Node_contrib_'+scheme+'_'+direction+'_'+str(lapse)+'.npy')
    
        # Compare node transmission to mean load
        if verbose:
            print('Plotting node comparison - '+scheme+' - '+direction)
        # sort node names for x-axis
        Total_usage = np.sum(N_usages,1)
        node_ids = np.array(range(len(N))).reshape((len(N),1))
        node_mean_load = [n.mean for n in N]
        
        # Calculate node proportional
        EU_load = np.sum(node_mean_load)
        if length:
            linkLengths = getLengths(lengths)
            Total_caps = sum(quantiles*linkLengths)
        else:
            Total_caps = sum(quantiles)
        Node_proportional = node_mean_load/EU_load*Total_caps/node_mean_load

        # Calculate link proportional
        link_proportional = linkProportional(N, link_dic, quantiles, lengths=lengths)
        link_proportional = [link_proportional[i]/node_mean_load[i] for i in range(nNodes)]
        
        # Calculate usage and sort countries by mean load
        normed_usage = Total_usage/node_mean_load
        normed_usage = np.reshape(normed_usage,(len(normed_usage),1))
        node_mean_load = np.reshape(node_mean_load,(len(node_mean_load),1))
        data = np.hstack([normed_usage,node_ids,node_mean_load,link_proportional])
        data_sort = data[data[:,2].argsort()]
        names_sort = [names[int(i)] for i in data_sort[:,1]]
        # flip order so largest is first
        names_sort = names_sort[::-1]
        link_proportional = data_sort[:,3][::-1]
        data_sort = data_sort[:,0][::-1]

        plt.figure(figsize=(10, 4), facecolor='w', edgecolor='k')
        ax = plt.subplot(111)
        green = '#009900'
        blue = '#000099'

        # Plot node proportional
        plt.rc('lines', lw=2)
        plt.rc('lines', dash_capstyle = 'round')
        plt.plot(np.linspace(0,len(N)*2+2,len(N)),Node_proportional,'--k')
        # Plot link proportional
        plt.bar(nodes,link_proportional,width=1,color=green,edgecolor='none')
        # Plot usage proportional
        plt.bar(nodes_shift,data_sort,width=1,color=blue,edgecolor='none')

        # Magic with ticks and labels
        if network == 'regions':
            ax.set_xticks(np.linspace(2,len(N)*2+2,len(N)+1))
            ax.set_xticklabels(names_sort,rotation=60,ha="right",va="top",fontsize=7)
        elif network == 'superRegions':
            ax.set_xticks(np.linspace(1.5,len(N)*2+1.5,len(N)+1))
            ax.set_xticklabels(names_sort,rotation=60,ha="right",va="top",fontsize=12)
        else:
            ax.set_xticks(np.linspace(2,len(N)*2+2,len(N)+1))
            ax.set_xticklabels(names_sort,rotation=60,ha="right",va="top",fontsize=10.5)

        ax.xaxis.grid(False)
        ax.xaxis.set_tick_params(width=0)
        if length:
            ax.set_ylabel(r'Network usage [MW$_T$km/MW$_L$]')
        else:
            ax.set_ylabel(r'Network usage [MW$_T$/MW$_L$]')
        maxes = [max(link_proportional), max(data_sort)]
        plt.axis([0,nNodes*2+.5,0,1.15*max(maxes)])

        # Legend
        artists = [plt.Line2D([0,0],[0,0],ls='dashed',lw=2.0,c='k'), plt.Rectangle((0,0),0,0,ec=green,fc=green), plt.Rectangle((0,0),0,0,ec=blue,fc=blue)]
        LABS = ['node proportional M$^1$','link proportional M$^2$','usage proportional M$^3$']
        leg = plt.legend(artists, LABS,loc='upper left',ncol=len(artists), columnspacing=0.6,borderpad=0.4, borderaxespad=0.0, handletextpad=0.2, handleheight = 1.2)
        leg.get_frame().set_alpha(0)
        leg.get_frame().set_edgecolor('white')
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=9.5)    # the legend text fontsize

        if (network == 'regions' or network == 'superRegions'):
            if length: plt.savefig('./figures/sensitivity/'+scheme+'/'+network+'-network-usage-'+direction+'-length.png', bbox_inches='tight')
            else: plt.savefig('./figures/sensitivity/'+scheme+'/'+network+'-network-usage-'+direction+'.png', bbox_inches='tight')
        else:
            if length:
                plt.savefig('./figures/'+scheme+'/network-usage-'+direction+'-length.png', bbox_inches='tight')
                plt.savefig('./figures/sensitivity/'+scheme+'/network-usage-'+direction+'-length.png', bbox_inches='tight')
            else:
                plt.savefig('./figures/'+scheme+'/network-usage-'+direction+'.png', bbox_inches='tight')
                plt.savefig('./figures/sensitivity/'+scheme+'/network-usage-'+direction+'.png', bbox_inches='tight')
            if verbose:
                print('Saved figures to ./figures/'+scheme+'/network-usage-'+direction+'.png')

def bars2(scheme, verbose=False):
    """
    Creates figures to compare usage proportional for networks with N = 30 and
    N = 53 and places them in ./figures/sensitivity/scheme/
    """
    # Load data and results
    N_regions = EU_Nodes_regions()
    F_regions = abs(np.load('./sensitivity/regions-'+scheme+'-flows.npy'))
    quantiles_regions = np.load('./sensitivity/regions-quantiles_'+scheme+'_'+str(lapse)+'.npy')

    N = EU_Nodes_usage()
    F = abs(np.load('./results/'+scheme+'-flows.npy'))
    quantiles = np.load('./results/quantiles_'+scheme+'_'+str(lapse)+'.npy')

    mergeDict = None # dictionary for merged regions, see regionMerger()
    names = node_namer(N) # array of node labels
    links = range(len(F))
    nNodes = 30
    nodes = np.linspace(0.5,2*nNodes-1.5,nNodes)
    bw = .5 # bar width
    nodes_shift = nodes+.5*bw

    for direction in directions:
        if length: N_usages = np.load('./results/Node_contrib_'+scheme+'_'+direction+'_length_'+str(lapse)+'.npy')
        else: N_usages = np.load('./results/Node_contrib_'+scheme+'_'+direction+'_'+str(lapse)+'.npy')

        # Calculate node-, link- and usage proportional for countries
        # Compare node transmission to mean load
        # sort node names for x-axis
        Total_usage = np.sum(N_usages,1)
        node_ids = np.array(range(len(N))).reshape((len(N),1))
        node_mean_load = [n.mean for n in N]

        # Calculate node proportional
        EU_load = np.sum(node_mean_load)
        if length:
            linkLengths = getLengths('countries')
            Total_caps = sum(quantiles*linkLengths)
        else:
            Total_caps = sum(quantiles)
        Node_proportional = node_mean_load/EU_load*Total_caps/node_mean_load

        # Calculate link proportional
        if length: link_proportional = linkProportional(N, link_dic, quantiles, lengths='countries')
        else: link_proportional = linkProportional(N, link_dic, quantiles)
        link_proportional = [link_proportional[i]/node_mean_load[i] for i in range(nNodes)]
        
        # Calculate usage and sort countries by mean load
        normed_usage = Total_usage/node_mean_load
        normed_usage = np.reshape(normed_usage,(len(normed_usage),1))
        node_mean_load = np.reshape(node_mean_load,(len(node_mean_load),1))

        # Calculate node-, link- and usage proportional for regions. Variable
        # names are the same as above with the addition of '_merged' or '_regions'
        region_usages, node_mean_load_merged, mergeDict = regionMerger(scheme, direction, lapse, md=mergeDict)
        Total_usage_merged = np.sum(region_usages,1)
        EU_load_merged = np.sum(node_mean_load_merged)
        if length:
            linkLengths = getLengths('regions')
            Total_caps_merged = sum(quantiles_regions*linkLengths)
        else:
            Total_caps_merged = sum(quantiles_regions)
        Node_proportional_merged = node_mean_load_merged/EU_load_merged*Total_caps_merged/node_mean_load_merged

        if length: link_proportional_regions = linkProportional(N_regions, link_dic_regions, quantiles_regions, lengths='regions')
        else: link_proportional_regions = linkProportional(N_regions, link_dic_regions, quantiles_regions)
        link_proportional_merged = simpleMerger(link_proportional_regions, mergeDict)
        link_proportional_merged = [link_proportional_merged[i]/node_mean_load_merged[i] for i in range(nNodes)]
        normed_usage_merged = Total_usage_merged/node_mean_load_merged

        normed_usage_merged = np.reshape(normed_usage_merged,(len(normed_usage_merged),1))
        node_mean_load_merged = np.reshape(node_mean_load_merged,(len(node_mean_load_merged),1))
        link_proportional_merged = np.reshape(link_proportional_merged,(len(link_proportional_merged),1))

        # Sort data for plotting
        data = np.hstack([normed_usage, node_ids, node_mean_load, link_proportional,
            normed_usage_merged, link_proportional_merged])
        data_sort = data[data[:,2].argsort()]
        names_sort = [names[int(i)] for i in data_sort[:,1]]
        # flip order so largest is first
        names_sort = names_sort[::-1]
        link_proportional = data_sort[:,3][::-1]
        usage_proportional = data_sort[:,0][::-1]
        link_proportional_merged = data_sort[:,5][::-1]
        usage_proportional_merged = data_sort[:,4][::-1]


        # PLOTTING
        plt.figure(figsize=(10, 4), facecolor='w', edgecolor='k')
        ax = plt.subplot(111)
        green = '#009900'
        blue = '#000099'

        # Plot node proportional
        plt.rc('lines', lw=2)
        plt.rc('lines', dash_capstyle = 'round')
        plt.plot(np.linspace(0,len(N)*2+2,len(N)),Node_proportional,'--k')
        # Plot link proportional
        plt.bar(nodes,link_proportional,width=bw,color=green,edgecolor='none')
        plt.bar(nodes+1.5*bw,link_proportional_merged,width=bw,color=green,edgecolor='none')
        # Plot usage proportional
        plt.bar(nodes_shift+.2*bw,usage_proportional,width=bw,color=blue,edgecolor='none')
        plt.bar(nodes_shift+1.7*bw,usage_proportional_merged,width=bw,color=blue,edgecolor='none')

        # Magic with ticks and labels
        ax.set_xticks(np.linspace(2,len(N)*2+2,len(N)+1))
        ax.set_xticklabels(names_sort,rotation=60,ha="right",va="top",fontsize=10.5)

        ax.xaxis.grid(False)
        ax.xaxis.set_tick_params(width=0)
        if length:
            ax.set_ylabel(r'Network usage [MW$_T$km/MW$_L$]')
        else:
            ax.set_ylabel(r'Network usage [MW$_T$/MW$_L$]')
        maxes = [max(link_proportional), max(usage_proportional), max(link_proportional_merged), max(usage_proportional_merged)]
        plt.axis([0,nNodes*2+.5,0,1.15*max(maxes)])

        # Legend
        artists = [plt.Line2D([0,0],[0,0],ls='dashed',lw=2.0,c='k'), plt.Rectangle((0,0),0,0,ec=green,fc=green), plt.Rectangle((0,0),0,0,ec=blue,fc=blue)]
        LABS = ['node proportional M$^1$','link proportional M$^2$','usage proportional M$^3$']
        leg = plt.legend(artists, LABS,loc='upper left',ncol=len(artists), columnspacing=0.6,borderpad=0.4, borderaxespad=0.0, handletextpad=0.2, handleheight = 1.2)
        leg.get_frame().set_alpha(0)
        leg.get_frame().set_edgecolor('white')
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=9.5)

        # Calculate correlations between countries and merged regions
        link1 = link_proportional[np.where(link_proportional_merged > 0)]
        link2 = link_proportional_merged[np.where(link_proportional_merged > 0)]
        linkCorr =  '%0.2f' % pearsonr(link1, link2)[0]
        usage1 = usage_proportional[np.where(usage_proportional_merged > 0)]
        usage2 = usage_proportional_merged[np.where(usage_proportional_merged > 0)]
        usageCorr =  '%0.2f' % pearsonr(usage1, usage2)[0]
        s = 'linkCorr = '+linkCorr+', usageCorr = '+usageCorr
        plt.text(6,.9*max(maxes),s)

        if length: plt.savefig('./figures/sensitivity/'+scheme+'/compared-30-53-usage-'+direction+'-length.png', bbox_inches='tight')
        else: plt.savefig('./figures/sensitivity/'+scheme+'/compared-30-53-usage-'+direction+'.png', bbox_inches='tight')
        if verbose:
            print('Saved figures to ./figures/sensitivity/'+scheme+'/compared-30-53-usage-'+direction+'.png')

def bars3(scheme, verbose=False):
    """
    Creates figures to compare usage proportional for networks with N = 30 and
    N = 8 and places them in ./figures/sensitivity/scheme/
    """
    """
    Creates figures to compare usage proportional for networks with N = 53,
    N = 50 and N = 8 and places them in ./figures/sensitivity/scheme/
    """
    # Load data and results
    N_superRegions = EU_Nodes_superRegions()
    F_superRregions = abs(np.load('./sensitivity/superRegions-'+scheme+'-flows.npy'))
    quantiles_superRegions = np.load('./sensitivity/superRegions-quantiles_'+scheme+'_'+str(lapse)+'.npy')

    N = EU_Nodes_usage()
    F = abs(np.load('./results/'+scheme+'-flows.npy'))
    quantiles = np.load('./results/quantiles_'+scheme+'_'+str(lapse)+'.npy')

    country_mergeDict = None
    names = node_namer(N_superRegions)
    countryNames = node_namer(N)
    nNodes = 8
    nodes = np.linspace(0.5,2*nNodes-1.5,nNodes)
    bw = .5 # bar width
    nodes_shift = nodes+.5*bw
    
    for direction in directions:
        if length: SR_usages = np.load('./sensitivity/superRegions-Node_contrib_'+scheme+'_'+direction+'_length_'+str(lapse)+'.npy')
        else: SR_usages = np.load('./sensitivity/superRegions-Node_contrib_'+scheme+'_'+direction+'_'+str(lapse)+'.npy')

        # SUPER REGIONS
        # Calculate node-, link- and usage proportional
        Total_SR_usage = np.sum(SR_usages,1)
        node_ids = np.array(range(nNodes)).reshape((nNodes,1))
        SR_mean_load = [n.mean for n in N_superRegions]

        # Calculate node proportional
        Total_load = np.sum(SR_mean_load)
        if length:
            linkLengths = getLengths('superRegions')
            Total_caps = sum(quantiles_superRegions*linkLengths)
        else:
            Total_caps = sum(quantiles_superRegions)
        SR_node_proportional = SR_mean_load/Total_load*Total_caps/SR_mean_load

        # Calculate link proportional
        if length: SR_link_proportional = linkProportional(N_superRegions, link_dic_superRegions, quantiles_superRegions, lengths='superRegions')
        else: SR_link_proportional = linkProportional(N_superRegions, link_dic_superRegions, quantiles_superRegions)
        SR_link_proportional = [SR_link_proportional[i]/SR_mean_load[i] for i in range(nNodes)]

        # Calculate usage and sort countries by mean load
        normed_usage = Total_SR_usage/SR_mean_load
        normed_usage = np.reshape(normed_usage,(len(normed_usage),1))
        SR_mean_load = np.reshape(SR_mean_load,(len(SR_mean_load),1))


        # COUNTRIES
        if length: country_usages = np.load('./results/Node_contrib_'+scheme+'_'+direction+'_length_'+str(lapse)+'.npy')
        else: country_usages = np.load('./results/Node_contrib_'+scheme+'_'+direction+'_'+str(lapse)+'.npy')
        country_usages, country_mean_load, country_mergeDict = superRegionMerger(N, country_usages, countryNames, md=country_mergeDict)
        
        Total_country_usage = np.sum(country_usages,1)
        Total_country_load = np.sum(country_mean_load)
        if length:
            linkLengths = getLengths('countries')
            Total_country_caps = sum(quantiles*linkLengths)
        else:
            Total_country_caps = sum(quantiles)
        country_node_proportional = country_mean_load/Total_country_load*Total_country_caps/country_mean_load
        if length: country_link_proportional = linkProportional(N, link_dic, quantiles, lengths='countries')
        else: country_link_proportional = linkProportional(N, link_dic, quantiles)
        country_link_proportional = simpleMerger(country_link_proportional, country_mergeDict)
        country_link_proportional = [country_link_proportional[i]/country_mean_load[i] for i in range(nNodes)]
        country_normed_usage = Total_country_usage/country_mean_load

        country_normed_usage = np.reshape(country_normed_usage,(len(country_normed_usage),1))
        country_mean_load = np.reshape(country_mean_load,(len(country_mean_load),1))
        country_link_proportional = np.reshape(country_link_proportional,(len(country_link_proportional),1))


        # Sort data for plotting
        data = np.hstack([normed_usage, node_ids, SR_mean_load, SR_link_proportional,
            country_normed_usage, country_link_proportional])
        data_sort = data[data[:,2].argsort()]
        names_sort = [names[int(i)] for i in data_sort[:,1]]
        # flip order so largest is first
        names_sort = names_sort[::-1]
        SR_link_proportional = data_sort[:,3][::-1]
        SR_usage_proportional = data_sort[:,0][::-1]
        country_link_proportional = data_sort[:,5][::-1]
        country_usage_proportional = data_sort[:,4][::-1]


        # PLOTTING
        plt.figure(figsize=(10, 4), facecolor='w', edgecolor='k')
        ax = plt.subplot(111)
        green = '#009900'
        blue = '#000099'

        # Plot node proportional
        plt.rc('lines', lw=2)
        plt.rc('lines', dash_capstyle = 'round')
        plt.plot(np.linspace(0,len(N_superRegions)*2+2,len(N_superRegions)),SR_node_proportional,'--k')
        # Plot link proportional
        plt.bar(nodes,SR_link_proportional,width=bw,color=green,edgecolor='none')
        plt.bar(nodes+1.5*bw,country_link_proportional,width=bw,color=green,edgecolor='none')
        # Plot usage proportional
        plt.bar(nodes_shift+.2*bw,SR_usage_proportional,width=bw,color=blue,edgecolor='none')
        plt.bar(nodes_shift+1.7*bw,country_usage_proportional,width=bw,color=blue,edgecolor='none')

        # Magic with ticks and labels
        ax.set_xticks(np.linspace(1.5,len(N)*2+1.5,len(N)+1))
        ax.set_xticklabels(names_sort,rotation=60,ha="right",va="top",fontsize=10.5)

        ax.xaxis.grid(False)
        ax.xaxis.set_tick_params(width=0)
        if length:
            ax.set_ylabel(r'Network usage [MW$_T$km/MW$_L$]')
        else:
            ax.set_ylabel(r'Network usage [MW$_T$/MW$_L$]')
        maxes = [max(SR_link_proportional), max(SR_usage_proportional), max(country_link_proportional), max(country_usage_proportional)]
        plt.axis([0,nNodes*2+.5,0,1.15*max(maxes)])

        # Legend
        artists = [plt.Line2D([0,0],[0,0],ls='dashed',lw=2.0,c='k'), plt.Rectangle((0,0),0,0,ec=green,fc=green), plt.Rectangle((0,0),0,0,ec=blue,fc=blue)]
        LABS = ['node proportional M$^1$','link proportional M$^2$','usage proportional M$^3$']
        leg = plt.legend(artists, LABS,loc='upper left',ncol=len(artists), columnspacing=0.6,borderpad=0.4, borderaxespad=0.0, handletextpad=0.2, handleheight = 1.2)
        leg.get_frame().set_alpha(0)
        leg.get_frame().set_edgecolor('white')
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=9.5)

        # Calculate correlations between countries and merged regions
        link1 = SR_link_proportional[np.where(country_link_proportional > 0)]
        link2 = country_link_proportional[np.where(country_link_proportional > 0)]
        linkCorr =  '%0.2f' % pearsonr(link1, link2)[0]
        usage1 = SR_usage_proportional[np.where(country_usage_proportional > 0)]
        usage2 = country_usage_proportional[np.where(country_usage_proportional > 0)]
        usageCorr =  '%0.2f' % pearsonr(usage1, usage2)[0]
        s = 'linkCorr = '+linkCorr+', usageCorr = '+usageCorr
        plt.text(2,.9*max(maxes),s)

        if length: plt.savefig('./figures/sensitivity/'+scheme+'/compared-30-8-usage-'+direction+'-length.png', bbox_inches='tight')
        else: plt.savefig('./figures/sensitivity/'+scheme+'/compared-30-8-usage-'+direction+'.png', bbox_inches='tight')
        if verbose:
            print('Saved figures to ./figures/sensitivity/'+scheme+'/compared-30-8-usage-'+direction+'.png')

def bars4(scheme):
    """
    Creates figures to compare usage proportional for networks with N = 53,
    N = 50 and N = 8 and places them in ./figures/sensitivity/scheme/
    """
    # Load data and results
    N_superRegions = EU_Nodes_superRegions()
    F_superRregions = abs(np.load('./sensitivity/superRegions-'+scheme+'-flows.npy'))
    quantiles_superRegions = np.load('./sensitivity/superRegions-quantiles_'+scheme+'_'+str(lapse)+'.npy')

    N_regions = EU_Nodes_regions()
    F_regions = abs(np.load('./sensitivity/regions-'+scheme+'-flows.npy'))
    quantiles_regions = np.load('./sensitivity/regions-quantiles_'+scheme+'_'+str(lapse)+'.npy')

    N = EU_Nodes_usage()
    F = abs(np.load('./results/'+scheme+'-flows.npy'))
    quantiles = np.load('./results/quantiles_'+scheme+'_'+str(lapse)+'.npy')

    country_mergeDict = None
    region_mergeDict = None
    names = node_namer(N_superRegions)
    countryNames = node_namer(N)
    nNodes = 8
    nodes = np.linspace(0.2,2*nNodes-1.8,nNodes)
    bw = .4 # bar width
    nodes2 = nodes+1.5*bw
    nodes3 = nodes2+1.5*bw
    
    for direction in directions:
        if length: SR_usages = np.load('./sensitivity/superRegions-Node_contrib_'+scheme+'_'+direction+'_length_'+str(lapse)+'.npy')
        else: SR_usages = np.load('./sensitivity/superRegions-Node_contrib_'+scheme+'_'+direction+'_'+str(lapse)+'.npy')

        # SUPER REGIONS
        # Calculate node-, link- and usage proportional
        Total_SR_usage = np.sum(SR_usages,1)
        node_ids = np.array(range(nNodes)).reshape((nNodes,1))
        SR_mean_load = [n.mean for n in N_superRegions]

        # Calculate node proportional
        Total_load = np.sum(SR_mean_load)
        if length:
            linkLengths = getLengths('superRegions')
            Total_caps = sum(quantiles_superRegions*linkLengths)
        else:
            Total_caps = sum(quantiles_superRegions)
        SR_node_proportional = SR_mean_load/Total_load*Total_caps/SR_mean_load

        # Calculate link proportional
        if length: SR_link_proportional = linkProportional(N_superRegions, link_dic_superRegions, quantiles_superRegions, lengths='superRegions')
        else: SR_link_proportional = linkProportional(N_superRegions, link_dic_superRegions, quantiles_superRegions)
        SR_link_proportional = [SR_link_proportional[i]/SR_mean_load[i] for i in range(nNodes)]

        # Calculate usage and sort countries by mean load
        normed_usage = Total_SR_usage/SR_mean_load
        normed_usage = np.reshape(normed_usage,(len(normed_usage),1))
        SR_mean_load = np.reshape(SR_mean_load,(len(SR_mean_load),1))


        # COUNTRIES
        if length: country_usages = np.load('./results/Node_contrib_'+scheme+'_'+direction+'_length_'+str(lapse)+'.npy')
        else: country_usages = np.load('./results/Node_contrib_'+scheme+'_'+direction+'_'+str(lapse)+'.npy')
        country_usages, country_mean_load, country_mergeDict = superRegionMerger(N, country_usages, countryNames, md=country_mergeDict)
        
        Total_country_usage = np.sum(country_usages,1)
        Total_country_load = np.sum(country_mean_load)
        if length:
            linkLengths = getLengths('countries')
            Total_country_caps = sum(quantiles*linkLengths)
        else:
            Total_country_caps = sum(quantiles)
        country_node_proportional = country_mean_load/Total_country_load*Total_country_caps/country_mean_load
        if length: country_link_proportional = linkProportional(N, link_dic, quantiles, lengths='countries')
        else: country_link_proportional = linkProportional(N, link_dic, quantiles)
        country_link_proportional = simpleMerger(country_link_proportional, country_mergeDict)
        country_link_proportional = [country_link_proportional[i]/country_mean_load[i] for i in range(nNodes)]
        country_normed_usage = Total_country_usage/country_mean_load

        country_normed_usage = np.reshape(country_normed_usage,(len(country_normed_usage),1))
        country_mean_load = np.reshape(country_mean_load,(len(country_mean_load),1))
        country_link_proportional = np.reshape(country_link_proportional,(len(country_link_proportional),1))
        
        # REGIONS
        region_usages, region_mean_load, region_mergeDict = regionMerger(scheme, direction, lapse, md=region_mergeDict)
        region_usages, region_mean_load, country_mergeDict = superRegionMerger(N_regions, region_usages, countryNames, region_mean_load, md=country_mergeDict)

        Total_region_usage = np.sum(region_usages,1)
        Total_region_load = np.sum(region_mean_load)
        if length:
            linkLengths = getLengths('regions')
            Total_region_caps = sum(quantiles_regions*linkLengths)
        else:
            Total_region_caps = sum(quantiles_regions)
        region_node_proportional = region_mean_load/Total_region_load*Total_region_caps/region_mean_load
        if length: region_link_proportional = linkProportional(N_regions, link_dic_regions, quantiles_regions, lengths='regions')
        else: region_link_proportional = linkProportional(N_regions, link_dic_regions, quantiles_regions)
        region_link_proportional = simpleMerger(region_link_proportional, country_mergeDict)
        region_link_proportional = [region_link_proportional[i]/region_mean_load[i] for i in range(nNodes)]
        region_normed_usage = Total_region_usage/region_mean_load

        region_normed_usage = np.reshape(region_normed_usage,(len(region_normed_usage),1))
        region_mean_load = np.reshape(region_mean_load,(len(region_mean_load),1))
        region_link_proportional = np.reshape(region_link_proportional,(len(region_link_proportional),1))

        # Sort data for plotting
        data = np.hstack([normed_usage, node_ids, SR_mean_load, SR_link_proportional, country_normed_usage, country_link_proportional, region_normed_usage, region_link_proportional])
        data_sort = data[data[:,2].argsort()]
        names_sort = [names[int(i)] for i in data_sort[:,1]]
        # flip order so largest is first
        names_sort = names_sort[::-1]
        SR_link_proportional = data_sort[:,3][::-1]
        SR_usage_proportional = data_sort[:,0][::-1]
        country_link_proportional = data_sort[:,5][::-1]
        country_usage_proportional = data_sort[:,4][::-1]
        region_link_proportional = data_sort[:,7][::-1]
        region_usage_proportional = data_sort[:,6][::-1]
        
        # PLOTTING
        plt.figure(figsize=(10, 4), facecolor='w', edgecolor='k')
        ax = plt.subplot(111)
        green = ['#009900', '#008800', '#007700']
        blue = ['#000099', '#000088', '#000077']

        # Plot node proportional
        plt.rc('lines', lw=2)
        plt.rc('lines', dash_capstyle = 'round')
        #plt.plot(np.linspace(0,nNodes*2+2,nNodes),SR_node_proportional,'--k')
        # Plot link proportional
        plt.bar(nodes,SR_link_proportional,width=bw,color=green[0],edgecolor='none')
        plt.bar(nodes2,country_link_proportional,width=bw,color=green[1],edgecolor='none')
        plt.bar(nodes3,region_link_proportional,width=bw,color=green[2],edgecolor='none')
        # Plot usage proportional
        plt.bar(nodes+.5*bw,SR_usage_proportional,width=bw,color=blue[0],edgecolor='none')
        plt.bar(nodes2+.5*bw,country_usage_proportional,width=bw,color=blue[1],edgecolor='none')
        plt.bar(nodes3+.5*bw,region_usage_proportional,width=bw,color=blue[2],edgecolor='none')

        # Magic with ticks and labels
        ax.set_xticks(nodes+.9)
        ax.set_xticklabels(names_sort,rotation=0,ha="left",va="top",fontsize=10.5)
        ax.xaxis.grid(False)
        ax.xaxis.set_tick_params(width=0)
        if length:
            ax.set_ylabel(r'Network usage [MW$_T$km/MW$_L$]')
        else:
            ax.set_ylabel(r'Network usage [MW$_T$/MW$_L$]')
        maxes = [max(SR_link_proportional), max(SR_usage_proportional), max(country_usage_proportional), max(country_link_proportional), max(region_usage_proportional), max(region_link_proportional)]
        plt.axis([0,nNodes*2+.5,0,1.11*max(maxes)])

        # Legend
        artists = [plt.Line2D([0,0],[0,0],ls='dashed',lw=2.0,c='k'), plt.Rectangle((0,0),0,0,ec=green[0],fc=green[0]), plt.Rectangle((0,0),0,0,ec=blue[0],fc=blue[0])]
        LABS = ['node proportional M$^1$','link proportional M$^2$','usage proportional M$^3$']
        leg = plt.legend(artists, LABS,loc='upper left',ncol=len(artists), columnspacing=0.6,borderpad=0.4, borderaxespad=0.0, handletextpad=0.2, handleheight = 1.2)
        leg.get_frame().set_alpha(0)
        leg.get_frame().set_edgecolor('white')
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=9.5)

        # Calculate correlations between countries and merged regions
        link1 = SR_link_proportional
        link2 = country_link_proportional
        link3 = region_link_proportional
        linkCorr1 =  '%0.2f' % pearsonr(link1, link2)[0]
        linkCorr2 =  '%0.2f' % pearsonr(link1, link3)[0]
        linkCorr3 =  '%0.2f' % pearsonr(link2, link3)[0]
        usage1 = SR_usage_proportional
        usage2 = country_usage_proportional
        usage3 = region_usage_proportional
        usageCorr1 =  '%0.2f' % pearsonr(usage1, usage2)[0]
        usageCorr2 =  '%0.2f' % pearsonr(usage1, usage3)[0]
        usageCorr3 =  '%0.2f' % pearsonr(usage2, usage3)[0]
        s1 = 'link(8,30) = '+linkCorr1+', link(8,53) = '+linkCorr2+', link(30,53) = '+linkCorr3
        s2 = 'usage(8,30) = '+usageCorr1+', usage(8,53) = '+usageCorr2+', usage(30,53) = '+usageCorr3
        plt.text(.5*bw,.9*max(maxes),s1,fontsize=8)
        plt.text(.5*bw,.82*max(maxes),s2,fontsize=8)
        
        if length: plt.savefig('./figures/sensitivity/'+scheme+'/compared-all-network-usage-'+direction+'-length.png', bbox_inches='tight')
        else: plt.savefig('./figures/sensitivity/'+scheme+'/compared-all-network-usage-'+direction+'.png', bbox_inches='tight')
        print('Saved figures to ./figures/sensitivity/'+scheme+'/compared-all-network-usage-'+direction+'.png')


if 'network' in task:
    print('Plotting network figures')
    N = EU_Nodes_usage()
    for scheme in schemes:
        print('Mode: '+scheme)
        for direction in directions:
            print('Direction: '+direction)
            drawnet_usage(N,scheme,direction)

if (('total' in task) and ('sensitivity' not in task)):
    if not length: print('Plotting total network usage')
    else: print('Plotting total network usage with lengths')
    lapse = 70128
    N = EU_Nodes_usage()
    print('Building link dictionary')
    link_dic = link_dict(N) # dictionary of links directly connected to each node
    print('Plotting')
    p = Pool(len(schemes))
    p.map(bars, schemes)

if (('total' in task) and ('sensitivity' in task)):
    if not length: print('Plotting total network usage for different networks')
    else: print('Plotting total network usage for different networks with lengths')
    lapse = 70128
    N = EU_Nodes_superRegions()
    nLinks = np.zeros(12)
    link_dic = link_dict(N,nLinks)
    for scheme in schemes:
        bars(scheme)

    N = EU_Nodes_regions()
    nLinks = np.zeros(96)
    link_dic = link_dict(N,nLinks)
    for scheme in schemes:
        bars(scheme)

if (('sensitivity' in task) and ('compare' in task) and ('all' not in task)):
    if not length: print('Comparing N=30 with N=8 and N=30 with N=53')
    else: print('Comparing N=30 with N=8 and N=30 with N=53 including lengths')

    lapse = 70128
    N = EU_Nodes_usage()
    link_dic = link_dict(N) # dictionary of links directly connected to each node
    
    N = EU_Nodes_regions()
    nLinks = np.zeros(96)
    link_dic_regions = link_dict(N,nLinks)

    N = EU_Nodes_superRegions()
    nlinks = np.zeros(12)
    link_dic_superRegions = link_dict(N,nlinks)

    for scheme in schemes:
        bars2(scheme)
        bars3(scheme)

if (('sensitivity' in task) and ('compare' in task) and ('all' in task)):
    if not length: print('Comparing N=53, N=30 and N=8 networks')
    else: print('Comparing N=53, N=30 and N=8 networks with lengths')
    for scheme in schemes:
        lapse = 70128
        N = EU_Nodes_usage()
        link_dic = link_dict(N) # dictionary of links directly connected to each node
        
        N = EU_Nodes_regions()
        nLinks = np.zeros(96)
        link_dic_regions = link_dict(N,nLinks)

        N = EU_Nodes_superRegions()
        nlinks = np.zeros(12)
        link_dic_superRegions = link_dict(N,nlinks)

        bars4(scheme)
