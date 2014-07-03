#! /usr/bin/env python
from pylab import *
from scipy import *
import numpy as np
from res_classes import *
#from classes import *
import networkx as nx
from res_tools import get_q#, ISO2LONG, biggestpair
import decimal
from EUgrid import *
from new_linkcolour_algorithm import get_neighbours
from matplotlib.colors import LinearSegmentedColormap
 

colwidth = (3.425)
dcolwidth = (2*3.425+0.236) 
cdict = {'red':  ((0.0, 1.0, 1.0),
                  (0.5, 0.0, 0.50),
                  (1.0, 0.0, 0.0)),
          'green': ((0.0, 0.0, 0.0),
                    (0.5, 1.0, 0.0),
                    (1.0, 0.0, 0.0)),
          'blue': ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                   (1.0, 1.0, 1.0))}
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
colors_countries = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers.
#sugar_colours = ['#F8CA00','#BD1550','#E97F02','#490A3D','#8A9B0F','#134B7C']
sugar_colours = ['#F8CA00','#AD0550','#FF9F22','#490A3D','#9AAB0F','#003B6C']
order = [24,11,10,22,3,9,17,13,4,23,12,0,14,16,26,25,15,6,18,2,8,7,19,21,20,5,1]
rolor= ( (19)/255. , (75)/255. , (124)/255.)
rrolor= ( 0.69,0.12,0.07 )
grolor= ( 0.6,0.67,0.06)
joloj= ( (236)/255. , (180)/255. , (131)/255.)
crolor=sugar_colours[4]

##################### Experiment #####################
colwidth = (3.425)
dcolwidth = (2*3.425+0.236)

blue = '#134b7c'
yellow = '#f8ca00'
orange = '#e97f02'
brown = '#876310'
green = '#4a8e05'
red = '#ae1215'
purple = '#4f0a3d'
darkred= '#4f1215'
pink = '#bd157d'
lightpink = '#d89bc2'
aqua = '#37a688'
darkblue = '#09233b'
lightblue = '#8dc1e0'
grayblue = '#4a7fa2'

blue_cycle = [darkblue, blue, grayblue, lightblue]

color_cycle = [blue, red, orange, purple, green, pink, lightblue, darkred, yellow]

au_cdict = {'red': ((0.0,int(yellow[1:3],16)/255.0,int(yellow[1:3],16)/255.0),
(0.5,int(green[1:3],16)/255.0,int(green[1:3],16)/255.0),
(1.0,int(blue[1:3],16)/255.0,int(blue[1:3],16)/255.0)),
'green': ((0.0,int(yellow[3:5],16)/255.0,int(yellow[3:5],16)/255.0),
(0.5,int(green[3:5],16)/255.0,int(green[3:5],16)/255.0),
(1.0,int(blue[3:5],16)/255.0,int(blue[3:5],16)/255.0)),
'blue': ((0.0,int(yellow[5:7],16)/255.0,int(yellow[5:7],16)/255.0),
(0.5,int(green[5:7],16)/255.0,int(green[5:7],16)/255.0),
(1.0,int(blue[5:7],16)/255.0,int(blue[5:7],16)/255.0))}

au_cmap = LinearSegmentedColormap('au_cmap',au_cdict,256)
	
	####################################################################

def nlargest(X,n=2):
    a=max(X)
    b=min(X)
    for i in X:
        if i>b and i<a:
            b=i
    return b
    
def draw_node_links_rel_network(n0,linksrel,maxi=0,ports="ex_and_im",new=False,mode="linear",N=None,title=None,alph=None,copper=True, modified=True,lapse=None): 
## links is the new_node_links_rel file or node_links_rel or variants of these.
#n0 is the countrynumber and the function draws the relative usage of the different links compared to the countrys total usage of links.
#Mode can be "linear", "square", "random" or "capped".
	close()
	if alph==None:
		alph="hetero"
	if title==None:
		title="powermix"
	if N==None:
		N=EU_Nodes()
	if lapse==None:
	    lapse=N[0].nhours
	    
	G=nx.Graph()
	nodelist=[]
	colours=[]
	alphas=np.zeros(len(N))
	alphas_remove=[]
	alpha=0
	tot_im_and_ex=0
	tot_im=0
	tot_ex=0
	tot_load=0
	for n in N:
		G.add_node(str(n.label))
		nodelist.append(str(n.label))
		
	
	if ports == "im":
	    if new:
		for t in range(lapse):
		    alpha+=(N[n0].power_mix[:,t]/(sum(N[n0].power_mix[:,t])+0.000000000001*0.000000000001)*N[n0].load[t])
		tot_load+=N[n0].mean*lapse
		print alpha
		alphas=alpha/tot_load
		alphas[n0]=0
		alphas=alphas/sum(alphas)
		for n in N:
		    if alphas[n.id]<=0.01:
			alphas_remove.append(n.label)
		print alphas
		print alphas_remove		
		s= sum(alphas)
		print s
	    else:
		if mode =="random":
		    for t in range(lapse):
			alpha+=(N[n0].power_mix[:,t]/(sum(N[n0].power_mix[:,t])+0.000000000001*0.000000000001)*N[n0].imports[t])	    
			tot_im+=N[n0].imports[t]
		    print alpha
		    alphas=alpha/tot_im
		    alphas[n0]=0
		    alphas=alphas/sum(alphas)
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		    print alphas
		    print alphas_remove		
		    s= sum(alphas)
		    print s
		else:
		    for t in range(lapse):
			alpha+=(N[n0].power_mix[:,t]/(sum(N[n0].power_mix[:,t])+0.000000000001*0.000000000001)*N[n0].get_import()[t])	    
			tot_im+=N[n0].get_import()[t]
		    print alpha
		    alphas=alpha/tot_im
		    alphas[n0]=0
		    alphas=alphas/sum(alphas)
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		    print alphas
		    print alphas_remove		
		    s= sum(alphas)
		    print s
		
	elif ports == "ex":
	    if new:
		for t in range(lapse):
		    alpha+=(N[n0].power_mix_ex[:,t]/(sum(N[n0].power_mix_ex[:,t])+0.000000000001*0.000000000001)*(N[n0].get_RES()[t]+N[n0].balancing[t])-N[n0].curtailment[t])	    
		    tot_ex+=N[n0].get_RES()[t]+N[n0].balancing[t]-N[n0].curtailment[t]
		print alpha
		alphas=alpha/tot_ex
		alphas[n0]=0
		alphas=alphas/sum(alphas)
		for n in N:
		    if alphas[n.id]<=0.01:
			alphas_remove.append(n.label)
		print alphas
		print alphas_remove		
		s= sum(alphas)
		print s
	    
	    else:
		if mode == "random":
		    for t in range(lapse):
			alpha+=N[n0].power_mix_ex[:,t]/(sum(N[n0].power_mix_ex[:,t])+0.000000000001*0.000000000001)*N[n0].exports[t]	    
			tot_ex+=N[n0].exports[t]
		    print alpha
		    alphas=alpha/tot_ex
		    alphas[n0]=0
		    alphas=alphas/sum(alphas)
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		    print alphas
		    print alphas_remove		
		    s= sum(alphas)
		    print s
		else:
		    for t in range(lapse):
			alpha+=N[n0].power_mix_ex[:,t]/(sum(N[n0].power_mix_ex[:,t])+0.000000000001*0.000000000001)*N[n0].get_export()[t]	    
			tot_ex+=N[n0].get_export()[t]
		    print alpha
		    alphas=alpha/tot_ex
		    alphas[n0]=0
		    alphas=alphas/sum(alphas)
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		    print alphas
		    print alphas_remove		
		    s= sum(alphas)
		    print s
	elif ports == "ex_and_im":
	    if new:
		for t in range(lapse):
		    alpha+=(N[n0].power_mix[:,t]/(sum(N[n0].power_mix[:,t])+0.000000000001*0.000000000001)*N[n0].load[t]+N[n0].power_mix_ex[:,t]/(sum(N[n0].power_mix_ex[:,t])+0.000000000001*0.000000000001)*(N[n0].get_RES()[t]+N[n0].balancing[t]-N[n0].curtailment[t]))
		    tot_im_and_ex+=N[n0].load[t]-N[n0].curtailment[t]+N[n0].get_RES()[t]+N[n0].balancing[t]
				
		#for n in N:
		    #alphas[n.id]=(alpha[n.id])
		print alpha
		print tot_im_and_ex
		alphas=alpha/tot_im_and_ex
		alphas[n0]=0
		alphas=alphas/sum(alphas)
		for n in N:
		    if alphas[n.id]<=0.01:
			alphas_remove.append(n.label)
		print alphas
		print alphas_remove		
		s= sum(alphas)
		print s			
		u=alphas[22]
		print u
	    else:
		if mode == "random":
		    for t in range(lapse):
			alpha+=N[n0].power_mix[:,t]/(sum(N[n0].power_mix[:,t])+0.00000000000000000001)*N[n0].imports[t]+N[n0].power_mix_ex[:,t]/(sum(N[n0].power_mix_ex[:,t])+0.000000000000000000001)*N[n0].exports[t]
			tot_im_and_ex+=N[n0].imports[t]+N[n0].exports[t]
		    print alpha
		    print "mojn"
		    print tot_im_and_ex
		    alphas=alpha/tot_im_and_ex
		    print "yo"
		    alphas[n0]=0
		    alphas=alphas/sum(alphas)
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		    print alphas
		    print alphas_remove		
		    s= sum(alphas)
		    print s			
		    u=alphas[22]
		    print u
		else:
		    for t in range(lapse):
			alpha+=N[n0].power_mix[:,t]/(sum(N[n0].power_mix[:,t])+0.00000000000000000000001)*N[n0].get_import()[t]+N[n0].power_mix_ex[:,t]/(sum(N[n0].power_mix_ex[:,t])+0.000000000001*0.000000000001)*N[n0].get_export()[t]
			tot_im_and_ex+=N[n0].get_import()[t]+N[n0].get_export()[t]
				    
		    #for n in N:
			#alphas[n.id]=(alpha[n.id])
		    print alpha
		    print tot_im_and_ex
		    alphas=alpha/tot_im_and_ex
		    alphas[n0]=0
		    alphas=alphas/sum(alphas)
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		    print alphas
		    print alphas_remove		
		    s= sum(alphas)
		    print s			
		    u=alphas[22]
		    print u
##################### Experiment #####################
	colwidth = (3.425)
	dcolwidth = (2*3.425+0.236)

	blue = '#134b7c'
	yellow = '#f8ca00'
	orange = '#e97f02'
	brown = '#876310'
	green = '#4a8e05'
	red = '#ae1215'
	purple = '#4f0a3d'
	darkred= '#4f1215'
	pink = '#bd157d'
	lightpink = '#d89bc2'
	aqua = '#37a688'
	darkblue = '#09233b'
	lightblue = '#8dc1e0'
	grayblue = '#4a7fa2'

	blue_cycle = [darkblue, blue, grayblue, lightblue]

	color_cycle = [blue, red, orange, purple, green, pink, lightblue, darkred, yellow]

	au_cdict = {'red': ((0.0,int(yellow[1:3],16)/255.0,int(yellow[1:3],16)/255.0),
	(0.5,int(green[1:3],16)/255.0,int(green[1:3],16)/255.0),
	(1.0,int(blue[1:3],16)/255.0,int(blue[1:3],16)/255.0)),
	'green': ((0.0,int(yellow[3:5],16)/255.0,int(yellow[3:5],16)/255.0),
	(0.5,int(green[3:5],16)/255.0,int(green[3:5],16)/255.0),
	(1.0,int(blue[3:5],16)/255.0,int(blue[3:5],16)/255.0)),
	'blue': ((0.0,int(yellow[5:7],16)/255.0,int(yellow[5:7],16)/255.0),
	(0.5,int(green[5:7],16)/255.0,int(green[5:7],16)/255.0),
	(1.0,int(blue[5:7],16)/255.0,int(blue[5:7],16)/255.0))}

	au_cmap = LinearSegmentedColormap('au_cmap',au_cdict,256)
########################################	
	testcmap=LinearSegmentedColormap.from_list("yellow_white_blue",[(1,1,0),(1,1,1),(0,0,1)])
	#cmap=mpl.cm.YlBu
	norm=mpl.colors.Normalize(vmin=0,vmax=0.5) #The biggest user 

	node_c=[]
    #ax4.arrow(0.15,0.35,0.0,-0.2,fc='k',ec='k',head_width=0.05,head_length=0.1)
	for n in N:
		node_c.append(au_cmap(  alphas[n.id]/0.5    ))
	node_c[n0] = (0,0,0,1)
	
	K,h,ListF=AtoKh_old(N)    

	x=len(ListF)
	
	m=max(linksrel[n0*x:(n0+1)*x])
	#print m
	p=m**(0.1)
	if maxi!=0:
	    #p=maxi**(0.1)
	    p=p+maxi**(0.1)-p
	
	#print p
	#if maxim !=0:
	    #p=maxim**(0.1)
	
	    
	#print p
	k=decimal.Decimal(str(p))
	k=round(k,2)
	
	for l in ListF:
		w=linksrel[n0*x+l[2]]		
		G.add_edge(l[0], l[1] , weight= w)
		
  
	pos=nx.spring_layout(G)

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
	if len(N)>27:
		pos['EE']=[1.0,0.94]
		pos['LV']=[0.95,0.83]
		pos['LT']=[0.87,0.72]
	if len(N)>30:
		pos['NS']=[0.325,.95]
		pos['NA']=[0.275,0.0]
	Graf=G
	if modified:	    
	    
	    summe=[d["weight"] for (u,v,d) in G.edges(data=True)]
	    #print G.edges()
	    #print summe
	    summe=sum(summe)
	    #print summe
	    counter=-1
	    for n in N:		
		idle_weights=[]		
		
		idle_weights=[d["weight"] for (u,v,d) in G.edges(data=True) if u==n.label or v==n.label]
		#print len(idle_weights)
		#if idle_weights[0]<=1.35:
		    #print "yo"
		if all(idle_weights<=p) and alphas[n.id]<=0.01:		    
		    G.remove_node(str(n.label))
		    del pos[str(n.label)]
		    nodelist.remove(str(n.label))
		    counter+=1
		    node_c.pop(n.id-counter)
	print nodelist
		    

	fig = figure(dpi=100,figsize=(13,7))

	ax1= fig.add_axes([-0.025,0.235,1.00,0.8]) #For displaying graph    
	nx.draw_networkx_nodes(G,pos,node_size=600,nodelist=nodelist,node_color=node_c,facecolor=(1,1,1))

#    nx.draw_networkx_edges(G,pos,edgelist=e1,width=1.5,edgecolor=(0.175,0.175,0.175))
	e0=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=p]	
	e1=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p and d['weight']<=p**2]
	e2=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**2 and d['weight']<=p**3]
	e3=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**3 and d['weight']<=p**4]
	e4=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**4 and d['weight']<=p**5]
	e5=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**5 and d['weight']<=p**6]
	e6=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**6 and d['weight']<=p**7]
	e7=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**7 and d['weight']<=p**8]
	e8=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**8 and d['weight']<=p**9]
	e9=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**9 and d['weight']<p**10]
	e10=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>=p**10-0.001]
	e100=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=p and u not in alphas_remove and v not in alphas_remove]
#    ax1.text(-0.05,1.05,"(a)",fontsize=12)
	if modified:
	    nx.draw_networkx_edges(G,pos,edgelist=e100,width=1.0,edge_color='c',alpha=1.0,style='dotted')
	else:
	    nx.draw_networkx_edges(G,pos,edgelist=e0,width=1.0,edge_color='c',alpha=1.0,style='dotted')
	nx.draw_networkx_edges(G,pos,edgelist=e1,width=0.5,edge_color='c',alpha=0.6,style='dashdot')
	nx.draw_networkx_edges(G,pos,edgelist=e2,width=1.0,edge_color='c',alpha=0.7,style='dashed')
	nx.draw_networkx_edges(G,pos,edgelist=e3,width=2.0,edge_color='c',alpha=0.8)
	nx.draw_networkx_edges(G,pos,edgelist=e4,width=3.0,edge_color=orange,alpha=0.9)
	nx.draw_networkx_edges(G,pos,edgelist=e5,width=3.0,edge_color='m',alpha=1.0)
	nx.draw_networkx_edges(G,pos,edgelist=e6,width=3.5,edge_color=aqua,alpha=.7)
	nx.draw_networkx_edges(G,pos,edgelist=e7,width=4.0,edge_color='g',alpha=.8)
	nx.draw_networkx_edges(G,pos,edgelist=e8,width=4.5,edge_color='b',alpha=0.9)
	nx.draw_networkx_edges(G,pos,edgelist=e9,width=5.0,edge_color='k',alpha=1.0)
	nx.draw_networkx_edges(G,pos,edgelist=e10,width=5.0,edge_color='r',alpha=1.0)
	nx.draw_networkx_labels(G,pos,font_size=13,font_color='w',font_family='sans-serif')	

#    nx.draw_networkx_labels(G,pos,font_size=8,font_color='k',font_family='sans-serif')

	ax1.axis('off')
	
	ax4= fig.add_axes([-0.075,0.075,1.5,.15]) #For displaying graph
	ax4.vlines(0.05*1.05+0.025,0.6,1.0,linewidth=1.0,color='c',alpha=1.0,linestyles='dotted')
	ax4.vlines(0.10*1.05+0.025,0.6,1.0,linewidth=0.5,color='c',alpha=0.6,linestyle='dashdot')
	ax4.vlines(0.15*1.05+0.025,0.6,1.0,linewidth=1.0,color='c',alpha=0.7,linestyle='dashed')
	ax4.vlines(0.20*1.05+0.025,0.6,1.0,linewidth=2.0,color='c',alpha=0.8)
	ax4.vlines(0.25*1.05+0.025,0.6,1.0,linewidth=3.0,color=orange,alpha=0.9)
	ax4.vlines(0.30*1.05+0.025,0.6,1.0,linewidth=3.0,color='m',alpha=1.0)
	ax4.vlines(0.35*1.05+0.025,0.6,1.0,linewidth=3.5,color=aqua,alpha=0.7)
	ax4.vlines(0.40*1.05+0.025,0.6,1.0,linewidth=4.0,color='g',alpha=0.8)
	ax4.vlines(0.45*1.05+0.025,0.6,1.0,linewidth=4.5,color='b',alpha=0.9)
	ax4.vlines(0.50*1.05+0.025,0.6,1.0,linewidth=5.0,color='k',alpha=1.0)
	ax4.vlines(0.55*1.05+0.025,0.6,1.0,linewidth=5.5,color='r',alpha=1.0)
	ax4.text(0.05*1.05+0.01,0.5,"$\leq$"+str(k)+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.10*1.05+0.01,0.5,"$\leq$"+str(round(p**2,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.15*1.05+0.01,0.5,"$\leq$"+str(round(p**3,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.20*1.05+0.01,0.5,"$\leq$"+str(round(p**4,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.25*1.05+0.01,0.5,"$\leq$"+str(round(p**5,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.30*1.05+0.01,0.5,"$\leq$"+str(round(p**6,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.35*1.05+0.01,0.5,"$\leq$"+str(round(p**7,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.40*1.05+0.01,0.5,"$\leq$"+str(round(p**8,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.45*1.05+0.01,0.5,"$\leq$"+str(round(p**9,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.50*1.05+0.01,0.5,r"<"+str(round(p**10,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.55*1.05+0.01,0.5,r"="+str(round(p**10,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.axis([0.0,1.0,0.0,1.2])
	ax4.axis('off')
	
	ax2= fig.add_axes([0.86,0.3,0.05,0.6]) #displaying colorbarlegend
	
	cbl=mpl.colorbar.ColorbarBase(ax2,au_cmap,norm,orientation= 'vertical' )
	#ax2.axis("off")
	alphas=np.array(alphas)
	if new:
	    if copper:
		if modified:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,copper, modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_linear_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
			np.save("./figures/new_linear_trade_alphas_relative_to_countrytrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,square,copper,modified,alpha="+str(alph), fontsize=9)
			np.save("./figures/new_square_trade_alphas_relative_to_countrytrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			savefig("./figures/new_squared_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		    elif mode == "random":
			np.save("./figures/new_random_trade_alphas_relative_to_countrytrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,random,copper,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_random_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			np.save("./figures/new_capped_trade_alphas_relative_to_countrytrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,copper,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		else:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_linear_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,square,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_squared_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,random,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_random_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
	    else:
		if modified:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_linear_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,square,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_squared_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,random,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_random_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		else:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_linear_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,square,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_squared_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,random,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_random_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
	else:
	    if copper:
		if modified:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,copper,modified,alpha="+str(alph), fontsize=9)
			np.save("./figures/old_linear_trade_alphas_relative_to_countrytrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			savefig("./figures/old_linear_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old,square,copper,modified,alpha="+str(alph), fontsize=9)
			np.save("./figures/old_squared_trade_alphas_relative_to_countrytrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			savefig("./figures/old_squared_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old random,copper,modified,alpha="+str(alph), fontsize=9)	
			np.save("./figures/old_random_trade_alphas_relative_to_countrytrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			savefig("./figures/old_random_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old capped,copper,modified,alpha="+str(alph), fontsize=9)	
			np.save("./figures/old_capped_trade_alphas_relative_to_countrytrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			savefig("./figures/old_capped_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		else:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_linear_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old,square,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_squared_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old random,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_random_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old capped,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_capped_link_usage_relative_to_countryflow_copper"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    
	    else:
		if modified:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_linear_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old,square,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_squared_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old random,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_random_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")	
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old capped,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_capped_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")	
		else:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_linear_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old,square,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_squared_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old random,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_random_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")	
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the countryflow for "+str(N[n0].label)+","+str(ports)+"ports,old capped,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_capped_link_usage_relative_to_countryflow_constr"+str(N[n0].id)+"_"+str(ports)+"ports_alpha="+str(alph)+".pdf")	
	


	
		
	
	show() # 
	return Graf

def draw_node_links_network(n0,links,linksflow,new=False,mode="linear",ports="ex_and_im",N=None,F=None,title=None,alph=None,copper=True,modified=True,lapse=None): 
## links is the new_links_mix file or links_mix or variants of these, linksflow is the links_flow file or variants of this, n0 is the countrynumber
# and it draws the countrys relative usage of the different links compared to the different flows over these links.
	#set new=True if you use new_links_mix
	close()
	if alph==None:
		alph="hetero"
	if title==None:
		title=="powermix"
	if N==None:
		N=EU_Nodes()
	if lapse==None:
	    lapse=N[0].nhours
	G=nx.Graph()
	nodelist=[]
	colours=[]
	alpha=np.zeros(len(N))
	tot_im_and_ex=np.zeros(len(N))
	alphas_remove=[]	
	for n in N:
		G.add_node(str(n.label))
		nodelist.append(str(n.label))
		n.power_mix_real=au.dcopy(n.power_mix)
		n.power_mix_ex_real=au.dcopy(n.power_mix_ex)
		
	

	
	K,h,ListF=AtoKh_old(N)
    
	x=len(N)
	m=0
	for i in range(len(ListF)):
		
		j=links[n0+i*x]/linksflow[i]*100
		if j>m:
			m=j
	p=m**(0.1)
	
	k=decimal.Decimal(str(p))
	k=round(k,2)
	
	for l in ListF:
		w=links[x*l[2]+n0]/linksflow[l[2]]*100	
		G.add_edge(l[0], l[1] , weight= w)
	##################################################### nodes usage #########################################
	if ports == "im":
	    if new:
		for t in range(lapse):
		    for n in N:
			alpha[n.id]+=(N[n.id].power_mix[n0,t]/(sum(N[n.id].power_mix[:,t])+0.000000000001*0.000000000001)*N[n.id].load[t]) #finds each nodes imports from each other node
		for n in N:
		    tot_load[n.id]=N[n.id].mean*lapse
		print alpha		
		alphas=alpha/tot_load		#normalises the imports from each node relative to its total imports	
		for n in N:
		    if alphas[n.id]<=0.01:
			alphas_remove.append(n.label) #removes the nodes with less than 1% 
				
		s= sum(alphas)
		print s #used to check, that they sum to 1
	    else:
		if mode =="random":
		    for t in range(lapse):
			for n in N:
			    alpha+=(N[n.id].power_mix[n0,t]/(sum(n.power_mix[:,t])+0.000000000001*0.000000000001)*n.imports[t])	    
		    for n in N:
			tot_im[n.id]=sum(n.imports)
		    
		    alphas=alpha/tot_im
		
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		   	
		    s= sum(alphas)
		    
		else:
		    for t in range(lapse):
			alpha+=(N[n0].power_mix[:,t]/(sum(N[n0].power_mix[:,t])+0.000000000001*0.000000000001)*N[n0].get_import()[t])	    
			tot_im+=N[n0].get_import()[t]
		    
		    alphas=alpha/tot_im
		    alphas[n0]=0
		    alphas=alphas/sum(alphas)
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		
		    s= sum(alphas)
		
		
	elif ports == "ex":
	    if new:
		for t in range(lapse):
		    alpha+=(N[n0].power_mix_ex[:,t]/(sum(N[n0].power_mix_ex[:,t])+0.000000000001*0.000000000001)*(N[n0].get_RES()[t]+N[n0].balancing[t])-N[n0].curtailment[t])	    
		    tot_ex+=N[n0].get_RES()[t]+N[n0].balancing[t]-N[n0].curtailment[t]
		
		alphas=alpha/tot_ex
		alphas[n0]=0
		alphas=alphas/sum(alphas)
		for n in N:
		    if alphas[n.id]<=0.01:
			alphas_remove.append(n.label)
		
		s= sum(alphas)
	    
	    
	    else:
		if mode == "random":
		    for t in range(lapse):
			alpha+=N[n0].power_mix_ex[:,t]/(sum(N[n0].power_mix_ex[:,t])+0.000000000001*0.000000000001)*N[n0].exports[t]	    
			tot_ex+=N[n0].exports[t]
		    
		    alphas=alpha/tot_ex
		    alphas[n0]=0
		    alphas=alphas/sum(alphas)
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		    	
		    s= sum(alphas)
		    
		else:
		    for t in range(lapse):
			alpha+=N[n0].power_mix_ex[:,t]/(sum(N[n0].power_mix_ex[:,t])+0.000000000001*0.000000000001)*N[n0].get_export()[t]	    
			tot_ex+=N[n0].get_export()[t]
		    
		    alphas=alpha/tot_ex
		    alphas[n0]=0
		    alphas=alphas/sum(alphas)
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		    	
		    s= sum(alphas)
		
	elif ports == "ex_and_im":
	    ############################ we find how much our node im- and exports with the different countrys relative the these countrys total im- and export from other countries ###################################
	    if new:
		
		for t in range(lapse):
		    if mod(t,100)==0 and t>0: 
			print "\r",round(100.0*(t/float(lapse)),2),"%",
			sys.stdout.flush()
		    

		    for n in N:
									
			n.power_mix[n.id,t]=0
			n.power_mix_ex[n.id,t]=0
			alpha[n.id]+=n.power_mix[n0,t]/(sum(n.power_mix[:,t])+0.000000000001*0.000000000001)*n.load[t]+n.power_mix_ex[n0,t]/(sum(n.power_mix_ex[:,t])+0.000000000001*0.000000000001)*(n.get_RES()[t]+n.balancing[t]-n.curtailment[t])
		for n in N:
		    tot_im_and_ex[n.id]=n.load.sum()+(n.get_RES()+n.balancing-n.curtailment).sum()	
		#for n in N:
		    #alphas[n.id]=(alpha[n.id]
		alphas=alpha/tot_im_and_ex
		
		for n in N:
		    if alphas[n.id]<=0.01:
			alphas_remove.append(n.label)
				
		s= sum(alphas)
			    
	    
	    
	    else:
		
		if mode == "random":
		    for t in range(lapse):
			for n in N:
			    n.power_mix[n.id,t]=0
			    n.power_mix_ex[n.id,t]=0
			    alpha[n.id]+=n.power_mix[n0,t]/(sum(n.power_mix[:,t])+0.00000000000000000001)*n.imports[t]+n.power_mix_ex[n0,t]/(sum(n.power_mix_ex[:,t])+0.000000000000000000001)*n.exports[t]
		    for n in N:
			tot_im_and_ex[n.id]=n.imports.sum()+n.exports.sum()
		    
		    alphas=alpha/tot_im_and_ex
		      
		    
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		    
		    s= sum(alphas)
		    print s
		    		
		    
		else:
		    
		    for t in range(lapse):
			if mod(t,100)==0 and t>0: 
			    print "\r",round(100.0*(t/float(lapse)),2),"%",
			    sys.stdout.flush()
			for n in N:
				n.power_mix[n.id,t]=0	
				n.power_mix_ex[n.id,t]=0	
				alpha[n.id]+=N[n.id].power_mix[n0,t]/(sum(N[n.id].power_mix[:,t])+0.00000000000000000000001)*N[n.id].get_import()[t]+N[n.id].power_mix_ex[n0,t]/(sum(N[n.id].power_mix_ex[:,t])+0.000000000001*0.000000000001)*N[n.id].get_export()[t]
		    for n in N:
			tot_im_and_ex[n.id]=N[n.id].get_import().sum()+N[n.id].get_export().sum()
				    
		    #for n in N:
			#alphas[n.id]=(alpha[n.id])
		    
		    alphas=alpha/tot_im_and_ex	
			
		    
		    for n in N:
			if alphas[n.id]<=0.01:
			    alphas_remove.append(n.label)
		   
		    	
		    s= sum(alphas)
		    print s			
		  
	for n in N:
	    
	    n.power_mix=n.power_mix_real
	    n.power_mix_ex=n.power_mix_ex_real
	###########################################################################################
		    
	
		

  
	pos=nx.spring_layout(G)


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
	if len(N)>27:
		pos['EE']=[1.0,0.94]
		pos['LV']=[0.95,0.83]
		pos['LT']=[0.87,0.72]
	if len(N)>30:
		pos['NS']=[0.325,.95]
		pos['NA']=[0.275,0.0]
	
	norm=mpl.colors.Normalize(vmin=0,vmax=0.5) #The biggest user 

	node_c=[]
    #ax4.arrow(0.15,0.35,0.0,-0.2,fc='k',ec='k',head_width=0.05,head_length=0.1)
	for n in N:
	    node_c.append(au_cmap(  alphas[n.id]/0.5    ))
	node_c[n0] = (0,0,0,1)

	if modified:
	    idle_nodes=[]
	    idle_weights=[]
	    counter=-1
	    
	    for n in N:
		idle_weights=[]
		 
		
		idle_weights=[d["weight"] for (u,v,d) in G.edges(data=True) if u==n.label or v==n.label] 
		if all(idle_weights<=p):
		    G.remove_node(str(n.label))
		    del pos[str(n.label)]
		    nodelist.remove(str(n.label))
		    counter+=1	    
		    node_c.pop(n.id-counter)	    
	       
	
	fig = figure(dpi=100,figsize=(13,7))
	

	ax1= fig.add_axes([-0.025,0.235,1.00,0.8]) #For displaying graph    
	
	nx.draw_networkx_nodes(G,pos,node_size=600,nodelist=nodelist,node_color=node_c,facecolor=(1,1,1))

	e0=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=p]
	e1=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p and d['weight']<=p**2]
	e2=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**2 and d['weight']<=p**3]
	e3=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**3 and d['weight']<=p**4]
	e4=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**4 and d['weight']<=p**5]
	e5=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**5 and d['weight']<=p**6]
	e6=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**6 and d['weight']<=p**7]
	e7=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**7 and d['weight']<=p**8]
	e8=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**8 and d['weight']<=p**9]
	e9=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**9]
	e10=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>=p**10-0.001]
	e100=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=p and u not in alphas_remove and v not in alphas_remove]

	if modified:
	    nx.draw_networkx_edges(G,pos,edgelist=e100,width=1.0,edge_color='c',alpha=1.0,style='dotted')
	else:
	    nx.draw_networkx_edges(G,pos,edgelist=e0,width=1.0,edge_color='c',alpha=1.0,style='dotted')
	nx.draw_networkx_edges(G,pos,edgelist=e1,width=0.5,edge_color='c',alpha=0.6,style='dashdot')
	nx.draw_networkx_edges(G,pos,edgelist=e2,width=1.0,edge_color='c',alpha=0.7,style='dashed')
	nx.draw_networkx_edges(G,pos,edgelist=e3,width=2.0,edge_color='c',alpha=0.8)
	nx.draw_networkx_edges(G,pos,edgelist=e4,width=3.0,edge_color=orange,alpha=0.9)
	nx.draw_networkx_edges(G,pos,edgelist=e5,width=3.0,edge_color='m',alpha=1.0)
	nx.draw_networkx_edges(G,pos,edgelist=e6,width=3.5,edge_color=aqua,alpha=.7)
	nx.draw_networkx_edges(G,pos,edgelist=e7,width=4.0,edge_color='g',alpha=.8)
	nx.draw_networkx_edges(G,pos,edgelist=e8,width=4.5,edge_color='b',alpha=0.9)
	nx.draw_networkx_edges(G,pos,edgelist=e9,width=5.0,edge_color='k',alpha=1.0)
	nx.draw_networkx_edges(G,pos,edgelist=e10,width=5.0,edge_color='r',alpha=1.0)
	nx.draw_networkx_labels(G,pos,font_size=13,font_color='w',font_family='sans-serif')	
	
	
	
	
	
	
	ax1.axis('off') 

	
	
	ax4= fig.add_axes([-0.075,0.075,1.5,.15]) #For displaying graph
	ax4.vlines(0.05*1.05+0.025,0.6,1.0,linewidth=1.0,color='c',alpha=1.0,linestyles='dotted')
	ax4.vlines(0.10*1.05+0.025,0.6,1.0,linewidth=0.5,color='c',alpha=0.6,linestyle='dashdot')
	ax4.vlines(0.15*1.05+0.025,0.6,1.0,linewidth=1.0,color='c',alpha=0.7,linestyle='dashed')
	ax4.vlines(0.20*1.05+0.025,0.6,1.0,linewidth=2.0,color='c',alpha=0.8)
	ax4.vlines(0.25*1.05+0.025,0.6,1.0,linewidth=3.0,color=orange,alpha=0.9)
	ax4.vlines(0.30*1.05+0.025,0.6,1.0,linewidth=3.0,color='m',alpha=1.0)
	ax4.vlines(0.35*1.05+0.025,0.6,1.0,linewidth=3.5,color=aqua,alpha=0.7)
	ax4.vlines(0.40*1.05+0.025,0.6,1.0,linewidth=4.0,color='g',alpha=0.8)
	ax4.vlines(0.45*1.05+0.025,0.6,1.0,linewidth=4.5,color='b',alpha=0.9)
	ax4.vlines(0.50*1.05+0.025,0.6,1.0,linewidth=5.0,color='k',alpha=1.0)
	ax4.vlines(0.55*1.05+0.025,0.6,1.0,linewidth=5.5,color='r',alpha=1.0)
	ax4.text(0.05*1.05+0.01,0.5,"$\leq$"+str(k)+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.10*1.05+0.01,0.5,"$\leq$"+str(round(p**2,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.15*1.05+0.01,0.5,"$\leq$"+str(round(p**3,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.20*1.05+0.01,0.5,"$\leq$"+str(round(p**4,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.25*1.05+0.01,0.5,"$\leq$"+str(round(p**5,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.30*1.05+0.01,0.5,"$\leq$"+str(round(p**6,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.35*1.05+0.01,0.5,"$\leq$"+str(round(p**7,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.40*1.05+0.01,0.5,"$\leq$"+str(round(p**8,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.45*1.05+0.01,0.5,"$\leq$"+str(round(p**9,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.50*1.05+0.01,0.5,r"<"+str(round(p**10,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.55*1.05+0.01,0.5,r"="+str(round(p**10,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.axis([0.0,1.0,0.0,1.2])
	
	ax4.axis("off")
	
	ax2= fig.add_axes([0.86,0.3,0.05,0.6]) #displaying colorbarlegend
	
	cbl=mpl.colorbar.ColorbarBase(ax2,au_cmap,norm,orientation= 'vertical' )
	
	
	alphas=np.array(alphas)
	if new:
	    if copper:
		if modified:
		    if mode =="linear":
			savefig("./figures/new_linear_trade_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,copper, modified,alpha="+str(alph), fontsize=9)
			np.save("./figures/new_linear_trade_alphas_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			#savefig("./figures/new_linear_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			savefig("./figures/new_squared"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,square,copper,modified,alpha="+str(alph), fontsize=9)
			np.save("./figures/new_squared_trade_alphas_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			#savefig("./figures/new_squared_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,random,copper,modified,alpha="+str(alph), fontsize=9)	
			np.save("./figures/new_random_trade_alphas_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			savefig("./figures/new_random_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,copper,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_trade_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
			np.save("./figures/new_capped_trade_alphas_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			#savefig("./figures/new_capped_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		else:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_linear_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,square,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_squared_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,random,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_random_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
	    else:
		if modified:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_linear_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,square,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_squared_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,random,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_random_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		else:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_linear_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,square,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_squared_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,random,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_random_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
	else:
	    if copper:
		if modified:
		    if mode =="linear":
			np.save("./figures/old_linear_trade_alphas_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			savefig("./figures/old_linear_trade_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,copper,modified,alpha="+str(alph), fontsize=9)
			#savefig("./figures/old_linear_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			savefig("./figures/old_squared_trade_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
			np.save("./figures/old_squared_trade_alphas_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,square,copper,modified,alpha="+str(alph), fontsize=9)
			#savefig("./figures/old_squared_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			np.save("./figures/old_random_trade_alphas_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			savefig("./figures/old_random_trade_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old random,copper,modified,alpha="+str(alph), fontsize=9)	
			#savefig("./figures/old_random_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			savefig("./figures/old_capped_trade_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph)+".pdf")
			np.save("./figures/old_capped_trade_alphas_relative_to_other_countriestrade_copper"+str(n0)+"_"+str(ports)+"ports_modified_alpha="+str(alph),alphas)
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old capped,copper,modified,alpha="+str(alph), fontsize=9)	
			#savefig("./figures/old_capped_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		else:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_linear_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+str(ports)+",ports,old,square,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_squared_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old random,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_random_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old capped,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_capped_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    
	    else:
		if modified:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_linear_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,square,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_squared_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old random,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_random_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")	
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old capped,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_capped_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")	
		else:
		    if mode =="linear":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_linear_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,square,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_squared_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old random,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_random_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")	
		    elif mode == "capped":
			#ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old capped,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_capped_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")	
	ax4.axis('off')
	

	show() # display
def draw_links_network(n0,links,linksflow,new=False,mode="linear",ports="ex and im",N=None,F=None,title=None,alph=None,copper=True,modified=False): 
## links is the new_links_mix file or links_mix or variants of these, linksflow is the links_flow file or variants of this, n0 is the countrynumber
# and it draws the countrys relative usage of the different links compared to the different flows over these links.
	#set new=True if you use new_links_mix
	close()
	if alph==None:
		alph="hetero"
	if title==None:
		title=="powermix"
	if N==None:
		N=EU_Nodes()
	G=nx.Graph()
	nodelist=[]
	colours=[]
	for n in N:
		G.add_node(str(n.label))
		nodelist.append(str(n.label))
		
	

	
	K,h,ListF=AtoKh_old(N)
    
	x=len(N)
	m=0
	for i in range(len(ListF)):
		
		j=links[n0+i*x]/linksflow[i]*100
		if j>m:
			m=j
	p=m**(0.1)
	print p	
	k=decimal.Decimal(str(p))
	k=round(k,2)
	print m
	for l in ListF:
		w=links[x*l[2]+n0]/linksflow[l[2]]*100	
		G.add_edge(l[0], l[1] , weight= w)
		

  
	pos=nx.spring_layout(G)


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
	if len(N)>27:
		pos['EE']=[1.0,0.94]
		pos['LV']=[0.95,0.83]
		pos['LT']=[0.87,0.72]
	if len(N)>30:
		pos['NS']=[0.325,.95]
		pos['NA']=[0.275,0.0]

	if modified:
	    idle_nodes=[]
	    idle_weights=[]
	    
	    for n in N:
		idle_weights=[]
		 
		
		idle_weights=[d["weight"] for (u,v,d) in G.edges(data=True) if u==n.label or v==n.label] 
		if all(idle_weights<=p):
		    G.remove_node(str(n.label))
		    del pos[str(n.label)]
		    nodelist.remove(str(n.label))	    	    
	       
	
	fig = figure(dpi=100,figsize=(13,7))
	

	ax1= fig.add_axes([-0.025,0.235,1.00,0.8]) #For displaying graph    
	nx.draw_networkx_nodes(G,pos,node_size=600,nodelist=nodelist,node_color=crolor,facecolor=(1,1,1))

	e0=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=p]
	e1=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p and d['weight']<=p**2]
	e2=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**2 and d['weight']<=p**3]
	e3=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**3 and d['weight']<=p**4]
	e4=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**4 and d['weight']<=p**5]
	e5=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**5 and d['weight']<=p**6]
	e6=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**6 and d['weight']<=p**7]
	e7=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**7 and d['weight']<=p**8]
	e8=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**8 and d['weight']<=p**9]
	e9=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>p**9]

	if modified:
	    nx.draw_networkx_edges(G,pos,edgelist=e0,width=1.0,edge_color='k',alpha=0.0,style='dotted')
	else:
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
	
	
	
	
	
	ax1.axis('off') 

#    nx.draw_networkx_labels(G,pos,font_size=8,font_color='k',font_family='sans-serif')
	

	ax4= fig.add_axes([-0.075,0.075,1.5,.15]) #For displaying graph
	ax4.vlines(0.06*1.05+0.025,0.6,1.0,linewidth=1.0,color='k',alpha=1.0,linestyles='dotted')
	ax4.vlines(0.12*1.05+0.025,0.6,1.0,linewidth=0.5,color='k',alpha=0.7,linestyle='dashed')
	ax4.vlines(0.18*1.05+0.025,0.6,1.0,linewidth=1.0,color='k',alpha=0.8,linestyle='dashed')
	ax4.vlines(0.24*1.05+0.025,0.6,1.0,linewidth=2.0,color='k',alpha=0.9,linestyle='dashed')
	ax4.vlines(0.30*1.05+0.025,0.6,1.0,linewidth=3.0,color='k',alpha=1.0,linestyle='dashed')
	ax4.vlines(0.36*1.05+0.025,0.6,1.0,linewidth=3.0,color='k',alpha=0.6)
	ax4.vlines(0.42*1.05+0.025,0.6,1.0,linewidth=3.5,color='k',alpha=0.7)
	ax4.vlines(0.48*1.05+0.025,0.6,1.0,linewidth=4.0,color='k',alpha=0.8)
	ax4.vlines(0.54*1.05+0.025,0.6,1.0,linewidth=4.5,color='k',alpha=0.9)
	ax4.vlines(0.60*1.05+0.025,0.6,1.0,linewidth=5.0,color='k',alpha=1.0)
	ax4.text(0.06*1.05+0.01,0.5,"$\leq$"+str(k)+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.12*1.05+0.01,0.5,"$\leq$"+str(round(p**2,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.18*1.05+0.01,0.5,"$\leq$"+str(round(p**3,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.24*1.05+0.01,0.5,"$\leq$"+str(round(p**4,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.30*1.05+0.01,0.5,"$\leq$"+str(round(p**5,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.36*1.05+0.01,0.5,"$\leq$"+str(round(p**6,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.42*1.05+0.01,0.5,"$\leq$"+str(round(p**7,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.48*1.05+0.01,0.5,"$\leq$"+str(round(p**8,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.54*1.05+0.01,0.5,"$\leq$"+str(round(p**9,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.text(0.60*1.05+0.01,0.5,"$\leq$"+str(round(p**10,2))+"$\%$",fontsize=9,rotation=-60)
	ax4.axis([0.0,1.0,0.0,1.2])
	
	if new:
	    if copper:
		if modified:
		    if mode =="linear":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,copper, modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_linear_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,square,copper,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_squared_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,random,copper,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_random_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,copper,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		else:
		    if mode =="linear":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_linear_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,square,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_squared_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,random,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_random_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
	    else:
		if modified:
		    if mode =="linear":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_linear_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,square,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_squared_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,random,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_random_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		else:
		    if mode =="linear":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,linear,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_linear_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,square,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/new_squared_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,random,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_random_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,new,capped,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/new_capped_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
	else:
	    if copper:
		if modified:
		    if mode =="linear":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,copper,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_linear_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,square,copper,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_squared_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old random,copper,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_random_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old capped,copper,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_capped_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		else:
		    if mode =="linear":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_linear_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+str(ports)+",ports,old,square,copper,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_squared_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old random,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_random_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "capped":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old capped,copper,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_capped_link_usage_relative_to_linkflow_copper"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    
	    else:
		if modified:
		    if mode =="linear":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_linear_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,square,constr,modified,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_squared_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old random,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_random_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")	
		    elif mode == "capped":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old capped,constr,modified,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_capped_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,modified,alpha="+str(alph)+".pdf")	
		else:
		    if mode =="linear":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,linear,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_linear_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "square":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old,square,constr,alpha="+str(alph), fontsize=9)
			savefig("./figures/old_squared_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")
		    elif mode == "random":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old random,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_random_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")	
		    elif mode == "capped":
			ax4.text(0.10,1.5, "Linkusage rel to the flow over the link for "+str(N[n0].label)+","+str(ports)+"ports,old capped,constr,alpha="+str(alph), fontsize=9)	
			savefig("./figures/old_capped_link_usage_relative_to_linkflow_constr"+str(N[n0].id)+","+str(ports)+"ports,alpha="+str(alph)+".pdf")	
	ax4.axis('off')
	

	show() # display

def draw_linkmix_network(l0,links,ports="ex and im",new=False,mode="linear",N=None,F=None,title=None,alph=None,copper=True): 
# l0 is the linknumber you want to look at, links is either links_mix.npy, new_links_mix.npy or variants of these.
# Draws each countrys share of the usage of the link. If you use new_links_mix remember to set new=True. 
	close()
	if alph==None:
		alph="hetero"
	if title==None:
		title=="linkmix"
	if N==None:
		N=EU_Nodes()
	G=nx.Graph()
	nodelist=[]
	colours=[]
	o=len(N)
	
	for n in N:
		G.add_node(str(n.label))
		nodelist.append(str(n.label))

	alpha=sum(links[o*l0:o*l0+30]) 
	print alpha
	

	
	alphas=[links[o*l0+x]/alpha*2 for x in range(len(N)) ] # We multiply by 2, because the colorbar is normed to 0.5 instead of 1. All our values in alphas are between 0 and 1 and the colors are in the range 0 to 1, so since the values in links are between 0 and 0.5 as max (no country can be using more than 50% of a link), it will now fit with the colourbar.
	print sum(alphas)
	

	#cmap=mpl.cm.PiYG
	au_cmap = LinearSegmentedColormap('au_cmap',au_cdict,256)
	norm=mpl.colors.Normalize(vmin=0,vmax=0.5) #The biggest user 

	node_c=[]
    
	for n in N:
		node_c.append(au_cmap(  alphas[n.id]    ))

	a,b,ListF=AtoKh_old(N)
	f,g,h,j,e=au.AtoKh(N)
	
    
	for l in ListF:
		w=4500	
		if l[2]==l0:
			G.add_edge(l[0],l[1],weight=13000)
		else:
			G.add_edge(l[0], l[1] , weight= w)
  
	pos=nx.spring_layout(G)


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
	if len(N)>27:
		pos['EE']=[1.0,0.94]
		pos['LV']=[0.95,0.83]
		pos['LT']=[0.87,0.72]
	if len(N)>30:
		pos['NS']=[0.325,.95]
		pos['NA']=[0.275,0.0]


	#fig = figure(dpi=100,figsize=(1.7*colwidth,1.7*colwidth*0.75))
	fig = figure(dpi=100,figsize=(13,7))

	ax1= fig.add_axes([-0.025,0.04,0.93,1.0]) #For displaying graph    
	nx.draw_networkx_nodes(G,pos,node_size=600,nodelist=nodelist,node_color=node_c,facecolor=(1,1,1))

	e0=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=700]
	e1=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>700 and d['weight']<=1200]
	e2=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1200 and d['weight']<=1800]
	e3=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1800 and d['weight']<=2400]
	e4=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>2400 and d['weight']<=3300]
	e5=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>3300 and d['weight']<=4000]
	e6=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>4000 and d['weight']<=5500]
	e7=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>5500 and d['weight']<=8000]
	e8=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>8000 and d['weight']<=12000]
	e9=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>12000]

	nx.draw_networkx_edges(G,pos,edgelist=e0,width=1.0,edge_color='k',alpha=1.0,style='dotted')
	nx.draw_networkx_edges(G,pos,edgelist=e1,width=0.5,edge_color='k',alpha=0.7,style='dashed')
	nx.draw_networkx_edges(G,pos,edgelist=e2,width=1.0,edge_color='k',alpha=0.8,style='dashed')
	nx.draw_networkx_edges(G,pos,edgelist=e3,width=2.0,edge_color='k',alpha=.9,style='dashed')
	nx.draw_networkx_edges(G,pos,edgelist=e4,width=3.0,edge_color='k',alpha=1.0,style='dashed')
	nx.draw_networkx_edges(G,pos,edgelist=e5,width=3.0,edge_color='k',alpha=0.6)
	nx.draw_networkx_edges(G,pos,edgelist=e6,width=3.5,edge_color='k',alpha=.2)
	nx.draw_networkx_edges(G,pos,edgelist=e7,width=4.0,edge_color='k',alpha=.8)
	nx.draw_networkx_edges(G,pos,edgelist=e8,width=4.5,edge_color='k',alpha=0.9)
	nx.draw_networkx_edges(G,pos,edgelist=e9,width=5.0,edge_color='k',alpha=1.0)
	nx.draw_networkx_labels(G,pos,font_size=13,font_color='k',font_family='sans-serif')
	
		
	ax1.axis('off') 

#    nx.draw_networkx_labels(G,pos,font_size=8,font_color='k',font_family='sans-serif')
	
	ax1.axis('off')
	ax2= fig.add_axes([0.86,0.05,0.05,0.92]) #displaying colorbarlegend
	cbl=mpl.colorbar.ColorbarBase(ax2,au_cmap,norm,orientation= 'vertical' )

	ax4= fig.add_axes([-0.075,0.075,1.5,.15]) #For displaying graph (the link)
	
	ax4.axis([0.0,1.0,0.0,1.2])
	ax4.axis('off')

    #plt.tight_layout
	
	if new:
	    if copper:
		if mode == "linear":
		    savefig("./figures/linkwork_new_linear_copper"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="square":
		    savefig("./figures/linkwork_new_square_copper"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="random":
		    savefig("./figures/linkwork_new_random_copper"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")	
	    else:
		if mode == "linear":
		    savefig("./figures/linkwork_new_linear_constr"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="square":
		    savefig("./figures/linkwork_new_square_constr"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="random":
		    savefig("./figures/linkwork_new_random_constr"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")		
	
	    
	else:
	    if copper:
		if mode == "linear":
		    savefig("./figures/linkwork_old_linear_copper"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="square":
		    savefig("./figures/linkwork_old_square_copper"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="random":
		    savefig("./figures/linkwork_old_random_copper"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")	
	    else:
		if mode == "linear":
		    savefig("./figures/linkwork_old_linear_constr"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="square":
		    savefig("./figures/linkwork_old_square_constr"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="random":
		    savefig("./figures/linkwork_old_random_constr"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")			
	
	show()
def draw_rel_linkmix_network(l0,links,ports="ex and im",new=False,mode="linear",N=None,F=None,title=None,alph=None, copper=True): 
# l0 is the linknumber you want to look at, 
## links is the new_node_links_rel file or node_links_rel or variants of these,
# Draws each countrys share of the usage of the link compared to the countrys overall usage of the network. 
#If you use new_node_links_rel remember to set new=True. 

######## In progress ############################################
	close()
	if alph==None:
		alph="hetero"
	if title==None:
		title=="linkmix"
	if N==None:
		N=EU_Nodes()
	G=nx.Graph()
	nodelist=[]
	colours=[]	
	
	for n in N:
		G.add_node(str(n.label))
		nodelist.append(str(n.label))
		
	K,h,ListF=AtoKh_old(N)    

	x=len(ListF)
	
	
	alphas=[links[n.id*x+l0]/100 for n in N]
	m=max(alphas)
	print m
	#print alphas
	
	
	cmap=mpl.cm.RdYlGn
	#norm=autoscale(alphas)
	norm=mpl.colors.Normalize(vmin=0,vmax=1)

	node_c=[]
    #ax4.arrow(0.15,0.35,0.0,-0.2,fc='k',ec='k',head_width=0.05,head_length=0.1)
	for n in N:
		node_c.append(cmap(  alphas[n.id]    ))
	
	f,g,h,j,e=au.AtoKh(EU_Nodes())
    
	for l in ListF:
		w=4500	
		if l[2]==l0:
			G.add_edge(l[0],l[1],weight=13000)
		else:
			G.add_edge(l[0], l[1] , weight= w)
  
	pos=nx.spring_layout(G)


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
	if len(N)>27:
		pos['EE']=[1.0,0.94]
		pos['LV']=[0.95,0.83]
		pos['LT']=[0.87,0.72]
	if len(N)>30:
		pos['NS']=[0.325,.95]
		pos['NA']=[0.275,0.0]


	fig = figure(dpi=100,figsize=(13,7))

	ax1= fig.add_axes([-0.025,0.04,0.93,1.0]) #For displaying graph    
	nx.draw_networkx_nodes(G,pos,node_size=600,nodelist=nodelist,node_color=node_c,facecolor=(1,1,1))
#    e1=[(u,v) for (u,v,d) in G.edges(data=True)]
#    nx.draw_networkx_edges(G,pos,edgelist=e1,width=1.5,edgecolor=(0.175,0.175,0.175))
	e0=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=700]
	e1=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>700 and d['weight']<=1200]
	e2=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1200 and d['weight']<=1800]
	e3=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1800 and d['weight']<=2400]
	e4=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>2400 and d['weight']<=3300]
	e5=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>3300 and d['weight']<=4000]
	e6=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>4000 and d['weight']<=5500]
	e7=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>5500 and d['weight']<=8000]
	e8=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>8000 and d['weight']<=12000]
	e9=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>12000]
#    ax1.text(-0.05,1.05,"(a)",fontsize=12)
	nx.draw_networkx_edges(G,pos,edgelist=e0,width=1.0,edge_color='k',alpha=1.0,style='dotted')
	nx.draw_networkx_edges(G,pos,edgelist=e1,width=0.5,edge_color='k',alpha=0.7,style='dashed')
	nx.draw_networkx_edges(G,pos,edgelist=e2,width=1.0,edge_color='k',alpha=0.8,style='dashed')
	nx.draw_networkx_edges(G,pos,edgelist=e3,width=2.0,edge_color='k',alpha=.9,style='dashed')
	nx.draw_networkx_edges(G,pos,edgelist=e4,width=3.0,edge_color='k',alpha=1.0,style='dashed')
	nx.draw_networkx_edges(G,pos,edgelist=e5,width=3.0,edge_color='k',alpha=0.6)
	nx.draw_networkx_edges(G,pos,edgelist=e6,width=3.5,edge_color='k',alpha=.2)
	nx.draw_networkx_edges(G,pos,edgelist=e7,width=4.0,edge_color='k',alpha=.8)
	nx.draw_networkx_edges(G,pos,edgelist=e8,width=4.5,edge_color='k',alpha=0.9)
	nx.draw_networkx_edges(G,pos,edgelist=e9,width=5.0,edge_color='k',alpha=1.0)
	nx.draw_networkx_labels(G,pos,font_size=13,font_color='k',font_family='sans-serif')
	
	
	if new:
	    if copper:
		if mode =="linear":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",new,linear,copper,alpha="+str(alph)+","+str(ports)+"ports", fontsize=9)
		elif mode == "square":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",new,square,copper,alpha="+str(alph)+","+str(ports)+"ports", fontsize=9)
		elif mode == "random":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",new,random,copper,alpha="+str(alph)+","+str(ports)+"ports", fontsize=9)	
	    else:
		if mode =="linear":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",new,linear,constr,alpha="+str(alph)+","+str(ports)+"ports", fontsize=9)
		elif mode == "square":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",new,square,constr,alpha="+str(alph)+","+str(ports)+"ports", fontsize=9)
		elif mode == "random":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",new,random,constr,alpha="+str(alph)+","+str(ports)+"ports", fontsize=9)	
	else:
	    if copper:
		if mode =="linear":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",old,linear,copper,alpha="+str(alph)+","+str(ports)+"ports", fontsize=9)
		elif mode == "square":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",old,square,copper,alpha="+str(alph)+","+str(ports)+"ports", fontsize=9)
		elif mode == "random":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",old,random,copper,alpha="+str(alph)+","+str(ports)+"ports", fontsize=9)	
	    else:
		if mode =="linear":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",old,linear,constr,alpha="+str(alph)+","+str(ports)+"ports", fontsize=9)
		elif mode == "square":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",old,square,constr,alpha="+str(alph)+","+str(ports)+"ports", fontsize=9)
		elif mode == "random":
		    ax1.text(0.03,0.02, "Each countrys linkusage relative to countryflow"+str(e[l0][0:1])+",old,random,constr,alph a="+str(alph)+","+str(ports)+"ports", fontsize=9)	
	
	ax1.axis('off') 

#    nx.draw_networkx_labels(G,pos,font_size=8,font_color='k',font_family='sans-serif')
	
	ax1.axis('off')
	ax2= fig.add_axes([0.86,0.05,0.05,0.92]) #displaying colorbarlegend
	cbl=mpl.colorbar.ColorbarBase(ax2,cmap,norm,orientation= 'vertical' )

	ax4= fig.add_axes([-0.075,0.075,1.5,.15]) #For displaying graph (the link)
	
	ax4.axis([0.0,1.0,0.0,1.2])
	ax4.axis('off')

    #plt.tight_layout
	
	if new:
	    if copper:
		if mode == "linear":
		    savefig("./figures/linkwork_new_linear_copper_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="square":
		    savefig("./figures/linkwork_new_square_copper_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="random":
		    savefig("./figures/linkwork_new_random_copper_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
	    else:
		if mode == "linear":
		    savefig("./figures/linkwork_new_linear_constr_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="square":
		    savefig("./figures/linkwork_new_square_constr_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="random":
		    savefig("./figures/linkwork_new_random_constr_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")	
	
	else:
	    if copper:
		if mode == "linear":
		    savefig("./figures/linkwork_old_linear_copper_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="square":
		    savefig("./figures/linkwork_old_square_copper_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="random":
		    savefig("./figures/linkwork_old_random_copper_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
	    else:
		if mode == "linear":
		    savefig("./figures/linkwork_old_linear_constr_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="square":
		    savefig("./figures/linkwork_old_square_constr_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
		elif mode=="random":
		    savefig("./figures/linkwork_old_random_constr_rel"+str(l0)+"alpha="+str(alph)+","+str(ports)+"ports.pdf")
	
	show()



def histogram_of_linkflow(l,F,n=None,links=None, ports="ex and im",alph=None,copper=True,new=False,N=None,mode="linear"):
#links should be a variant of new_linear_copper_link_mix_export_all_alpha=, l is the linknumber and n the countrynumber we want to look at.
    close()
    if alph == None:
	alph="hetero"
    if N == None:
	N=EU_Nodes()
	
    f,g,h,j,e=au.AtoKh(EU_Nodes())
    fig=plt.figure(figsize=(18,10))
    counter=0
    linksnew=[]
   # for t in range(70128):
	#if links[l][n][t]>1:
	    #linksnew.append(links[l][n][t])
	
	    #counter+=1
    #print counter

    plt.ion()
    plt.show()
    #plt.hist(linksnew,normed=1,bins=np.linspace(min(linksnew),max(linksnew),100))
    plt.hist(F[l],normed=1,bins=np.linspace(min(F[l]),max(F[l]),100))	
    ax=plt.gca()
    labels = ax.get_xticks().tolist()
    labels[0]="1"
    ax.set_xticklabels(labels)
    

    if new:
	if copper:
	    if mode == "linear":
		plt.title("histogram of the " +str(ports)+ "ports over link " +str(e[l][0:1])+" for "+str(N[n].label)+",new,linear,copper,alpha="+str(alph), fontsize=15)
		savefig("./figures/histogram_new_linear_copper"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")
	    elif mode=="square":
		plt.title("histogram of the " +str(ports)+ "ports over link " +str(e[l][0:1])+" for "+str(N[n].label)+",new,square,copper,alpha="+str(alph), fontsize=15)
		savefig("./figures/histogram_new_square_copper"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")
	    elif mode=="random":
		plt.title("histogram of the " +str(ports)+ "ports over link " +str(e[l][0:1])+" for "+str(N[n].label)+",new,random,copper,alpha="+str(alph), fontsize=15)
		savefig("./figures/histogram_new_random_copper"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")
	else:
	    if mode == "linear":
		plt.title("histogram of the " +str(ports)+ "ports over link " +str(e[l][0:1])+" for "+str(N[n].label)+",new,linear,constr,alpha="+str(alph), fontsize=15)
		savefig("./figures/histogram_new_linear_constr"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")
	    elif mode=="square":
		plt.title("histogram of the " +str(ports)+ "ports over link " +str(e[l][0:1])+" for "+str(N[n].label)+",new,square,constr,alpha="+str(alph), fontsize=15)
		savefig("./figures/histogram_new_square_constr"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")
	    elif mode=="random":
		plt.title("histogram of the " +str(ports)+ "ports over link " +str(e[l][0:1])+" for "+str(N[n].label)+",new,random,constr,alpha="+str(alph), fontsize=15)
		savefig("./figures/histogram_new_random_constr"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")	
	    
    else:
	if copper:
	    if mode == "linear":
		plt.title("Histogram of the " +str(ports)+ "ports over link " + str(e[l][0:1])+",old,linear,copper,alpha="+str(alph), fontsize=15)
		print("mojn")
		savefig("./figures/histogram_old_linear_copper"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")
	    elif mode=="square":
		plt.title("histogram of the " +str(ports)+ "ports over link " +str(e[l][0:1])+" for "+str(N[n].label)+",old,square,copper,alpha="+str(alph), fontsize=15)
		savefig("./figures/histogram_old_square_copper"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")
	    elif mode=="random":
		plt.title("histogram of the " +str(ports)+ "ports over link " +str(e[l][0:1])+" for "+str(N[n].label)+",old,random,copper,alpha="+str(alph), fontsize=15)
		savefig("./figures/histogram_old_random_copper"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")
	else:
	    if mode == "linear":
		plt.title("histogram of the " +str(ports)+ "ports over link " +str(e[l][0:1])+" for "+str(N[n].label)+",old,linear,constr,alpha="+str(alph), fontsize=15)
		savefig("./figures/histogram_old_linear_constr"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")
	    elif mode=="square":
		plt.title("histogram of the " +str(ports)+ "ports over link " +str(e[l][0:1])+" for "+str(N[n].label)+",old,square,constr,alpha="+str(alph), fontsize=15)
		savefig("./figures/histogram_old_square_constr"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")
	    elif mode=="random":
		plt.title("histogram of the " +str(ports)+ "ports over link " +str(e[l][0:1])+" for "+str(N[n].label)+",old,random,constr,alpha="+str(alph), fontsize=15)
		savefig("./figures/histogram_old_random_constr"+str(l)+"_"+str(n)+"alpha="+str(alph)+"_"+str(ports)+"ports.pdf")
plt.show()
	
