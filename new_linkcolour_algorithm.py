import csv
#import time
import sys
from pylab import *
from scipy.optimize import brentq
from scipy.optimize import fmin
#import datetime as dt
from localutils import *
from figutils import AtoKh_old
from time import time
from drawnet import draw_powmix_network
from EUgrid import *

def get_neighbours(n,matr=np.genfromtxt("./settings/eadmat.txt")):
    """ matr should be = np.genfromtxt("./settings/eadmat.txt")"""
    neighbours=[]    
    for j in range(len(matr[n])):
        if matr[n][j]>0 and n<>j:
            neighbours.append(j)
    return neighbours

def get_links(node_id,matr=np.genfromtxt("./settings/eadmat.txt")):
	""" matr should be = np.genfromtxt("./settings/eadmat.txt")"""
	links=[]
	linkno=-1
	for j in range(len(matr[node_id])):
		for i in range(len(matr[node_id])):
			if matr[i][j]>0 and i>j:
				linkno+=1
				if node_id==i: links.append([linkno,-1])
				if node_id==j: links.append([linkno,1])
	return links
    
def get_link_direction(n,N): #getting the first country-node in the link-label as start and the second as end
	a,b,c,d,e=au.AtoKh(EU_Nodes()) 
	
	for k in N:
		if str(k.label)==e[n][0][0:2]: start_node=k
	for k in N:
		if str(k.label)==e[n][0][6:8]: end_node=k
	
	return [start_node,end_node]

def track_new_square_link_usage_total(N,F,lapse=None,alph=None): #tracks the usage of each link for each country in the shared-balancing-scenario,alph sets the value of homogenous alphas
	
	if lapse== None:
		lapse=N[0].nhours   
	#if alph!=None:
		#file=open('new_square_link_usage_total_alpha='+str(alph)+'.dat','w+')
	
	#file = open('new_square_link_usage_total.dat', 'w+')
	
	matr=np.genfromtxt(N.pathadmat)
       
	a,b,c,d,e=au.AtoKh(EU_Nodes())     
	f,g,list_F=au.AtoKh_old(N)
	
	boxplot=np.zeros((4,len(e)))
	boksplotlabel=[]
	links=[]
	linksflow=[]
	cost=[]
	 		
	
	for n in range(len(e)):
		
		total_flow_link=0 
		total_link_mix=np.zeros((len(N)))   		
		link_mix_export=np.zeros((len(N),lapse))		
		link_mix_import=np.zeros((len(N),lapse))	
		total_link_mix_percentage=np.zeros((len(N)))
		#getting the first country-node in the link-label as start and the second as end
		start=get_link_direction(n,N)[0]
		end=get_link_direction(n,N)[1]
		# determining the links connected to start
		start.links=get_links(start.id,matr) 

		for t in range(lapse):
		## Update progress bar.
			if mod(t,100)==0 and t>0: 
				print "\r",round(100.0*(t/float(lapse)),2),"%",
				sys.stdout.flush()
			
    
			for l in start.links:
			    if end.label== list_F[l[0]][0] or end.label == list_F[l[0]][1]:
			#if flow indicates export from start_node			
				if F[l[0],t]*l[1]>0:
					if l[1]==1:
						sending_node_label=list_F[l[0]][0]
						receiving_node_label=list_F[l[0]][1]
					elif l[1]==-1:
						sending_node_label=list_F[l[0]][1] 
						receiving_node_label=list_F[l[0]][0]
					else:
						print "Warning (234nsd23): Link direction unknown!"
			#if flow indicates export from end_node
				else:				
					if l[1]==-1:
						sending_node_label=list_F[l[0]][0]
						receiving_node_label=list_F[l[0]][1]
					elif l[1]==1:
						sending_node_label=list_F[l[0]][1] 
						receiving_node_label=list_F[l[0]][0]
					else:
						print "Warning (234nsd23): Link direction unknown!"
				#finding the nodes corresponding to the sending and receiving nodes
				for k in N:
							
					if k.label==sending_node_label: 
						sending_node=k
					elif k.label==receiving_node_label:
						receiving_node=k
				# if sending_node and receiving_node are at the ends of the link we are looking at, create arrays with the flows being imported and exported from and to each country					
				if (sending_node==start and receiving_node==end) or (sending_node==end and receiving_node==start):
					link_mix_import[:,t]=receiving_node.power_mix_ex[:,t]*abs(F[l[0],t])/(sum(receiving_node.power_mix_ex[:,t])+0.000000000001)
					link_mix_export[:,t]=sending_node.power_mix[:,t]*abs(F[l[0],t])/(sum(sending_node.power_mix[:,t])+0.000000000001)
					
				#adding the flow over the link	for each timestep									
					total_flow_link+=abs(F[l[0],t])
				#adding the arrays with import and export to the array with the total_link_mix
					total_link_mix[:] += (link_mix_export[:,t]+link_mix_import[:,t])/2	
									
		# calculating each countrys usage of the link in percentages of the total usage of the link			   	    
		total_link_mix_percentage[:]=total_link_mix[:]/total_flow_link*100		
		print sum(total_link_mix_percentage[:])		
		links=np.append(links,total_link_mix[:])
		linksflow=np.append(linksflow,total_flow_link)			
		cost=np.append(cost,total_link_mix_percentage)
		
		# saving information to human-readable file
		#for k in N:			
			#x="country", str(k.label), "has a total usage of link", e[n], "of" ,str(total_link_mix_percentage[k.id]), "% after", str(t), "hours"
			#file.write(str(x) + '\n')			
			
		#file.write('\n')
		# sorting the total_link_mix_percentage, such that the higher percentages are first
		a=sorted(total_link_mix_percentage,reverse=True)
		#for the 4 highest values we save the usage for each link and the corresponding countries
		for i in range(4):
			boxplot[i,n]=a[i]
			for h in N:
				
				if a[i]==total_link_mix_percentage[h.id]:
					boksplotlabel=np.append(boksplotlabel,str(h.label))
	#if alph==None:
		np.save('linkcolouring/new_square_boxplot',boxplot) 
		np.save('linkcolouring/new_square_boksplotlabel',boksplotlabel)
		np.save('linkcolouring/new_square_links_mix',links)
		np.save('linkcolouring/new_square_links_flow',linksflow)
		np.save('linkcolouring/new_square_cost',cost)
	#else:
		np.save('linkcolouring/new_square_boxplot_alpha='+str(alph),boxplot) 
		np.save('linkcolouring/new_square_boksplotlabel_alpha='+str(alph),boksplotlabel)
		np.save('linkcolouring/new_square_links_mix_alpha='+str(alph),links)
		np.save('linkcolouring/new_square_links_flow_alpha='+str(alph),linksflow)
		np.save('linkcolouring/new_square_cost_alpha='+str(alph),cost)
	
	
	return boxplot,boksplotlabel

def track_link_usage_total(N,F,new=False,lapse=None,alph=None,mode="linear",copper=True): 
#tracks the usage of each link for each country, 
#alph sets the value of homogenous alphas.
#Mode can be "linear", "square" or "random".
#set new=True if new model is used.
	if alph==None:
	    alph="hetero"
	if lapse== None:
		lapse=N[0].nhours 
	if new:
	    if mode=="linear":
		
		if copper:
		    file=open('new_linear_copper_link_usage_total_alpha='+str(alph)+'.dat','w+')
		
	
		file = open('new_linear_constr_link_usage_total_alpha='+str(alph)+'.dat', 'w+')
	    elif mode=="square":
		
		if copper:
		    file=open('new_square_copper_link_usage_total_alpha='+str(alph)+'.dat','w+')
		
	
		file=open('new_square_constr_link_usage_total_alpha='+str(alph)+'.dat','w+')
	    
	    elif mode=="random":
		
		if copper:
		    file=open('new_random_copper_link_usage_total_alpha='+str(alph)+'.dat','w+')
		
	
		file=open('new_random_constr_link_usage_total_alpha='+str(alph)+'.dat','w+')
	else:
	    if mode=="linear":
		
		if copper:
		    file=open('old_linear_copper_link_usage_total_alpha='+str(alph)+'.dat','w+')
		
	
		file=open('old_linear_constr_link_usage_total_alpha='+str(alph)+'.dat','w+')
		
	    elif mode=="square":
		
		if copper:
		    file=open('old_square_copper_link_usage_total_alpha='+str(alph)+'.dat','w+')
		
	
		file=open('old_random_constr_link_usage_total_alpha='+str(alph)+'.dat','w+')
	    
	    elif mode=="random":
		
		if copper:
		    file=open('old_random_copper_link_usage_total_alpha='+str(alph)+'.dat','w+')
		
	
		file=open('old_random_constr_link_usage_total_alpha='+str(alph)+'.dat','w+')
	
	matr=np.genfromtxt(N.pathadmat)
       
	a,b,c,d,e=au.AtoKh(EU_Nodes())     
	f,g,list_F=au.AtoKh_old(N)
	
	boxplot=np.zeros((4,len(e)))
	boksplotlabel=[]
	links=[]
	linksflow=[]
	cost=[]
	links_ex=[]
	links_im=[]
	link_mix_export_all=[]
	link_mix_import_all=[]
	counter=0 	
	link_mix_export_variance_all=[]
	link_mix_import_variance_all=[]
	maxlex=[]
	maxlim=[]
	minlex=[]
	minlim=[]
	meanlex=[]
	meanlim=[]
	quant99lex=[]
	quant99lim=[]
		
	
	for n in range(len(e)):
	    
		counter+=1	
		link_mix_export_variance=np.zeros(len(N))
		link_mix_import_variance=np.zeros(len(N))	
		total_flow_link=0 
		total_export_mix=np.zeros((len(N)))
		total_import_mix=np.zeros((len(N)))
		total_link_mix=np.zeros((len(N)))   		
		link_mix_export=np.zeros((len(N),lapse))		
		link_mix_import=np.zeros((len(N),lapse))	
		total_link_mix_percentage=np.zeros((len(N)))
		#getting the first country-node in the link-label as start and the second as end
		start=get_link_direction(n,N)[0]
		end=get_link_direction(n,N)[1]
		# determining the links connected to start
		start.links=get_links(start.id,matr) 

		for t in range(lapse):
		## Update progress bar.
			if mod(t,100)==0 and t>0: 
				print "\r",round(100.0*(t/float(lapse)),2),"%",
				sys.stdout.flush()
			
    
			for l in start.links:
			#if flow indicates export from start_node
			    if end.label== list_F[l[0]][0] or end.label == list_F[l[0]][1]:
				if F[l[0],t]*l[1]>0:
					if l[1]==1:
						sending_node_label=list_F[l[0]][0]
						receiving_node_label=list_F[l[0]][1]
					elif l[1]==-1:
						sending_node_label=list_F[l[0]][1] 
						receiving_node_label=list_F[l[0]][0]
					else:
						print "Warning (234nsd23): Link direction unknown!"
			#if flow indicates export from end_node
				else:				
					if l[1]==-1:
						sending_node_label=list_F[l[0]][0]
						receiving_node_label=list_F[l[0]][1]
					elif l[1]==1:
						sending_node_label=list_F[l[0]][1] 
						receiving_node_label=list_F[l[0]][0]
					else:
						print "Warning (234nsd23): Link direction unknown!"
				#finding the nodes corresponding to the sending and receiving nodes
				for k in N:
							
					if k.label==sending_node_label: 
						sending_node=k
					elif k.label==receiving_node_label:
						receiving_node=k
				# if sending_node and receiving_node are at the ends of the link we are looking at, create arrays with the flows being imported and exported from and to each country					
				if (sending_node==start and receiving_node==end) or (sending_node==end and receiving_node==start):
					link_mix_import[:,t]=receiving_node.power_mix_ex[:,t]*abs(F[l[0],t])/(sum(receiving_node.power_mix_ex[:,t])+0.000000000001)
					link_mix_export[:,t]=sending_node.power_mix[:,t]*abs(F[l[0],t])/(sum(sending_node.power_mix[:,t])+0.000000000001)
					
				#adding the flow over the link	for each timestep									
					total_flow_link+=abs(F[l[0],t])
				#adding the arrays with import and export to the array with the total_link_mix
					total_link_mix[:] += (link_mix_export[:,t]+link_mix_import[:,t])/2	
					total_export_mix[:] += link_mix_export[:,t]
					total_import_mix[:] += link_mix_import[:,t]
									
		# calculating each countrys usage of the link in percentages of the total usage of the link		   	    
		total_link_mix_percentage[:]=total_link_mix[:]/total_flow_link*100		
		#print sum(total_link_mix_percentage[:])		
		links=np.append(links,total_link_mix[:])
		linksflow=np.append(linksflow,total_flow_link)		
		cost=np.append(cost,total_link_mix_percentage)
		links_ex=np.append(links_ex,total_export_mix[:])
		links_im=np.append(links_im,total_import_mix[:])
		link_mix_export_all.append(link_mix_export)
		link_mix_import_all.append(link_mix_import)
		
		
		for l in range(len(N)):
		    m=np.mean(link_mix_export[l])
		    o=np.mean(link_mix_import[l])
		    link_mix_export_variance[l]=(sum(((link_mix_export[l,t]-m)**2 for t in range(lapse)))/lapse)
		    link_mix_import_variance[l]=(sum(((link_mix_import[l,t]-o)**2 for t in range(lapse)))/lapse)
		    maxlex.append(max(link_mix_export_all[n][l]))
		    maxlim.append(max(link_mix_import_all[n][l]))
		    minlex.append(min(link_mix_export_all[n][l]))
		    minlim.append(min(link_mix_import_all[n][l]))
		    meanlex.append(np.mean(link_mix_export_all[n][l]))
		    meanlim.append(np.mean(link_mix_import_all[n][l]))	
		    quant99lex.append(np.percentile(link_mix_export_all[n][l],99))
		    quant99lim.append(np.percentile(link_mix_import_all[n][l],99))	    
		link_mix_export_variance_all.append(link_mix_export_variance)
		link_mix_import_variance_all.append(link_mix_import_variance)
		print sum(total_link_mix_percentage[:])		
		print counter
		
		
		# saving information to human-readable file
		for k in N:			
			x="country", str(k.label), "has a total usage of link", e[n], "of" ,str(total_link_mix_percentage[k.id]), "% after", str(t), "hours"
			file.write(str(x) + '\n')			
			
		file.write('\n')
		# sorting the total_link_mix_percentage, such that the higher percentages are first
		a=sorted(total_link_mix_percentage,reverse=True)
		#for the 4 highest values we save the usage for each link and the corresponding countries
		for i in range(4):
			boxplot[i,n]=a[i]
			for h in N:
				
				if a[i]==total_link_mix_percentage[h.id]:
				    boksplotlabel=np.append(boksplotlabel,str(h.label))
				    
	maxlex=np.array(maxlex).reshape((50,30))
	maxlim=np.array(maxlim).reshape((50,30))
	minlex=np.array(minlex).reshape((50,30))
	minlim=np.array(minlim).reshape((50,30))
	meanlex=np.array(meanlex).reshape((50,30))
	meanlim=np.array(meanlim).reshape((50,30))
	quant99lex=np.array(quant99lex).reshape((50,30))
	quant99lim=np.array(quant99lim).reshape((50,30))
	
	if new:
	    if copper:
	    				
		if mode=="linear":
			np.save('linkcolouring/new_linear_copper_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/new_linear_copper_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/new_linear_copper_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/new_linear_copper_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/new_linear_copper_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/new_linear_copper_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/new_linear_copper_links_im_alpha='+str(alph),links_im)		   
			np.save('linkcolouring/new_linear_copper_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/new_linear_copper_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/new_linear_copper_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/new_linear_copper_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/new_linear_copper_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/new_linear_copper_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/new_linear_copper_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/new_linear_copper_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/new_linear_copper_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/new_linear_copper_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/new_linear_copper_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/new_linear_copper_quant99lim_alpha='+str(alph),quant99lim)
		elif mode=="square":
			np.save('linkcolouring/new_square_copper_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/new_square_copper_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/new_square_copper_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/new_square_copper_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/new_square_copper_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/new_square_copper_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/new_square_copper_links_im_alpha='+str(alph),links_im)		   
			np.save('linkcolouring/new_square_copper_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/new_square_copper_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/new_square_copper_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/new_square_copper_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/new_square_copper_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/new_square_copper_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/new_square_copper_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/new_square_copper_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/new_square_copper_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/new_square_copper_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/new_square_copper_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/new_square_copper_quant99lim_alpha='+str(alph),quant99lim)

	    
		elif mode=="random":
			np.save('linkcolouring/new_random_copper_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/new_random_copper_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/new_random_copper_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/new_random_copper_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/new_random_copper_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/new_random_copper_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/new_random_copper_links_im_alpha='+str(alph),links_im)		    
			np.save('linkcolouring/new_random_copper_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/new_random_copper_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/new_random_copper_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/new_random_copper_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/new_random_copper_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/new_random_copper_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/new_random_copper_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/new_random_copper_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/new_random_copper_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/new_random_copper_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/new_random_copper_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/new_random_copper_quant99lim_alpha='+str(alph),quant99lim)
	    else:
		if mode=="linear":
			np.save('linkcolouring/new_linear_constr_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/new_linear_constr_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/new_linear_constr_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/new_linear_constr_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/new_linear_constr_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/new_linear_constr_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/new_linear_constr_links_im_alpha='+str(alph),links_im)		   
			np.save('linkcolouring/new_linear_constr_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/new_linear_constr_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/new_linear_constr_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/new_linear_constr_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/new_linear_constr_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/new_linear_constr_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/new_linear_constr_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/new_linear_constr_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/new_linear_constr_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/new_linear_constr_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/new_linear_constr_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/new_linear_constr_quant99lim_alpha='+str(alph),quant99lim)
		elif mode=="square":
			np.save('linkcolouring/new_square_constr_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/new_square_constr_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/new_square_constr_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/new_square_constr_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/new_square_constr_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/new_square_constr_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/new_square_constr_links_im_alpha='+str(alph),links_im)		   
			np.save('linkcolouring/new_square_constr_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/new_square_constr_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/new_square_constr_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/new_square_constr_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/new_square_constr_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/new_square_constr_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/new_square_constr_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/new_square_constr_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/new_square_constr_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/new_square_constr_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/new_square_constr_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/new_square_constr_quant99lim_alpha='+str(alph),quant99lim)

	    
		elif mode=="random":
			np.save('linkcolouring/new_random_constr_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/new_random_constr_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/new_random_constr_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/new_random_constr_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/new_random_constr_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/new_random_constr_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/new_random_constr_links_im_alpha='+str(alph),links_im)		    
			np.save('linkcolouring/new_random_constr_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/new_random_constr_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/new_random_constr_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/new_random_constr_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/new_random_constr_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/new_random_constr_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/new_random_constr_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/new_random_constr_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/new_random_constr_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/new_random_constr_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/new_random_constr_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/new_random_constr_quant99lim_alpha='+str(alph),quant99lim)
		

	else:
	    if copper:
		if mode=="linear":
		    
			np.save('linkcolouring/old_linear_copper_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/old_linear_copper_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/old_linear_copper_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/old_linear_copper_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/old_linear_copper_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/old_linear_copper_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/old_linear_copper_links_im_alpha='+str(alph),links_im)		
			np.save('linkcolouring/old_linear_copper_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/old_linear_copper_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/old_linear_copper_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/old_linear_copper_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/old_linear_copper_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/old_linear_copper_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/old_linear_copper_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/old_linear_copper_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/old_linear_copper_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/old_linear_copper_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/old_linear_copper_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/old_linear_copper_quant99lim_alpha='+str(alph),quant99lim)
			
		    
		elif mode=="square":
		    
			np.save('linkcolouring/old_square_copper_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/old_square_copper_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/old_square_copper_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/old_square_copper_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/old_square_copper_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/old_square_copper_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/old_square_copper_links_im_alpha='+str(alph),links_im)		   
			np.save('linkcolouring/old_square_copper_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/old_square_copper_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/old_square_copper_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/old_square_copper_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/old_square_copper_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/old_square_copper_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/old_square_copper_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/old_square_copper_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/old_square_copper_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/old_square_copper_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/old_square_copper_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/old_square_copper_quant99lim_alpha='+str(alph),quant99lim)

	    
		elif mode=="random":
		    
			np.save('linkcolouring/old_random_copper_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/old_random_copper_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/old_random_copper_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/old_random_copper_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/old_random_copper_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/old_random_copper_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/old_random_copper_links_im_alpha='+str(alph),links_im)		   
			np.save('linkcolouring/old_random_copper_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/old_random_copper_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/old_random_copper_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/old_random_copper_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/old_random_copper_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/old_random_copper_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/old_random_copper_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/old_random_copper_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/old_random_copper_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/old_random_copper_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/old_random_copper_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/old_random_copper_quant99lim_alpha='+str(alph),quant99lim)
	    else:
		if mode=="linear":
		    
			np.save('linkcolouring/old_linear_constr_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/old_linear_constr_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/old_linear_constr_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/old_linear_constr_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/old_linear_constr_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/old_linear_constr_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/old_linear_constr_links_im_alpha='+str(alph),links_im)		
			np.save('linkcolouring/old_linear_constr_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/old_linear_constr_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/old_linear_constr_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/old_linear_constr_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/old_linear_constr_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/old_linear_constr_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/old_linear_constr_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/old_linear_constr_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/old_linear_constr_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/old_linear_constr_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/old_linear_constr_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/old_linear_constr_quant99lim_alpha='+str(alph),quant99lim)
		    
		elif mode=="square":
		    
			np.save('linkcolouring/old_square_constr_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/old_square_constr_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/old_square_constr_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/old_square_constr_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/old_square_constr_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/old_square_constr_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/old_square_constr_links_im_alpha='+str(alph),links_im)		   
			np.save('linkcolouring/old_square_constr_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/old_square_constr_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/old_square_constr_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/old_square_constr_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/old_square_constr_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/old_square_constr_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/old_square_constr_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/old_square_constr_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/old_square_constr_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/old_square_constr_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/old_square_constr_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/old_square_constr_quant99lim_alpha='+str(alph),quant99lim)

	    
		elif mode=="random":
		    
			np.save('linkcolouring/old_random_constr_boxplot_alpha='+str(alph),boxplot) 
			np.save('linkcolouring/old_random_constr_boksplotlabel_alpha='+str(alph),boksplotlabel)
			np.save('linkcolouring/old_random_constr_links_mix_alpha='+str(alph),links)
			np.save('linkcolouring/old_random_constr_links_flow_alpha='+str(alph),linksflow)
			np.save('linkcolouring/old_random_constr_cost_alpha='+str(alph),cost)
			np.save('linkcolouring/old_random_constr_links_ex_alpha='+str(alph),links_ex)
			np.save('linkcolouring/old_random_constr_links_im_alpha='+str(alph),links_im)		   
			np.save('linkcolouring/old_random_constr_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
			np.save('linkcolouring/old_random_constr_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			np.save('linkcolouring/old_random_constr_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
			np.save('linkcolouring/old_random_constr_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
			np.save('linkcolouring/old_random_constr_maxlex_alpha='+str(alph),maxlex)
			np.save('linkcolouring/old_random_constr_maxlim_alpha='+str(alph),maxlim)
			np.save('linkcolouring/old_random_constr_minlex_alpha='+str(alph),minlex)
			np.save('linkcolouring/old_random_constr_minlim_alpha='+str(alph),minlim)
			np.save('linkcolouring/old_random_constr_meanlex_alpha='+str(alph),meanlex)
			np.save('linkcolouring/old_random_constr_meanlim_alpha='+str(alph),meanlim)
			np.save('linkcolouring/old_random_constr_quant99lex_alpha='+str(alph),quant99lex)
			np.save('linkcolouring/old_random_constr_quant99lim_alpha='+str(alph),quant99lim)

		
	
	return boxplot,boksplotlabel
	
def boksplot(a,b,alph=None,new=False,mode="linear",copper=True): # draws a columndiagram displaying the 4 biggest users and their usage of each link. 
#Insert new_boxplot and new_boksplotlabel or variants of these.
#Mode can be "linear", "square" or "random".
	if alph==None:
	    alph="hetero"	
	c=np.zeros((50))
	f,g,h,j,e=au.AtoKh(EU_Nodes())
	width=1.0
	n=np.arange(len(e))
	E = [i[0] for i in e]
	fig=plt.figure(figsize=(18,10))
	for i in range(50):
		c[i]=100
		plt.text(i+width/8,c[i]-a[0,i]+0.4,b[i*4],color='w')
		plt.text(i+width/8.,c[i]-a[0,i]-a[1,i]+0.4,b[1+i*4],color='w')
		plt.text(i+width/8.,c[i]-a[0,i]-a[1,i]-a[2,i]+0.4,b[2+i*4],color='w')
		plt.text(i+width/8.,c[i]-a[0,i]-a[1,i]-a[2,i]-a[3,i]+0.4,b[3+i*4],color='w')
	plt.ion()
	plt.show()
	p1=plt.bar(range(50),c,width,color=au.blue)
	p2=plt.bar(range(50),c-a[0,:],width,color=au.red)
	p3=plt.bar(range(50),c-a[0,:]-a[1,:],width,color=au.gree)
	p4=plt.bar(range(50),c-a[0,:]-a[1,:]-a[2,:],width,color=au.oran)
	p5=plt.bar(range(50),c-a[0,:]-a[1,:]-a[2,:]-a[3,:],width,color=au.purp)
	plt.xticks(n+width/2.,E,rotation=90)
	plt.yticks(np.arange(0,101,5))
	subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
	plt.legend( (p1[0], p2[0],p3[0],p4[0],p5[0]), ('Primary user', 'Secondary user','Tertiary user','Quaternary user', 'The rest'), 'upper center',fancybox=True,prop={'size':10})
	
	
	if new:
	    if copper:
		if mode =="linear":
		    plt.text(0.03,92, "Linkusg rel to linkflow,new,linear,copper,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_linear_copper_columndiagram_alpha="+str(alph)+".pdf")
		elif mode == "square":
		    plt.text(0.03,92, "Linkusg rel to linkflow,new,square,copper,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_square_copper_columndiagram_alpha="+str(alph)+".pdf")
		elif mode == "random":
		    plt.text(0.03,92, "Linkusg rel to linkflow,new,random,copper,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_random_copper_columndiagram_alpha="+str(alph)+".pdf")	
	    else:
		if mode =="linear":
		    plt.text(0.03,92, "Linkusg rel to linkflow,new,linear,constr,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_linear_constr_columndiagram_alpha="+str(alph)+".pdf")
		elif mode == "square":
		    plt.text(0.03,92, "Linkusg rel to linkflow,new,square,constr,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_square_constr_columndiagram_alpha="+str(alph)+".pdf")
		elif mode == "random":
		    plt.text(0.03,92, "Linkusg rel to linkflow,new,random,constr,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_random_constr_columndiagram_alpha="+str(alph)+".pdf")
	else:
	    if copper:
		if mode =="linear":
		    plt.text(0.03,92, "Linkusg rel to linkflow,old,linear,copper,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_linear_copper_columndiagram_alpha="+str(alph)+".pdf")
		elif mode == "square":
		    plt.text(0.03,92, "Linkusg rel to linkflow,old,square,copper,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_square_copper_columndiagram_alpha="+str(alph)+".pdf")
		elif mode == "random":
		    plt.text(0.03,92, "Linkusg rel to linkflow,old,random,copper,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_random_copper_columndiagram_alpha="+str(alph)+".pdf")
	    else:
		if mode =="linear":
		    plt.text(0.03,92, "Linkusg rel to linkflow,old,linear,constr,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_linear_constr_columndiagram_alpha="+str(alph)+".pdf")
		elif mode == "square":
		    plt.text(0.03,92, "Linkusg rel to linkflow,old,square,constr,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_square_constr_columndiagram_alpha="+str(alph)+".pdf")
		elif mode == "random":
		    plt.text(0.03,92, "Linkusg rel to linkflow,old,random,constr,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_random_constr_columndiagram_alpha="+str(alph)+".pdf")
		
	
	
		
	    

def new_square_boksplot(a,b,alph=None): # draws a columndiagram displaying the 4 biggest users and their usage of each link. 
#Insert new_boxplot and new_boksplotlabel or variants of these.

		
	c=np.zeros((50))
	f,g,h,j,e=au.AtoKh(EU_Nodes())
	width=1.0
	n=np.arange(len(e))
	E = [i[0] for i in e]
	for i in range(50):
		c[i]=100
		plt.text(i+width/8,c[i]-a[0,i]+0.4,b[i*4],color='w')
		plt.text(i+width/8.,c[i]-a[0,i]-a[1,i]+0.4,b[1+i*4],color='w')
		plt.text(i+width/8.,c[i]-a[0,i]-a[1,i]-a[2,i]+0.4,b[2+i*4],color='w')
		plt.text(i+width/8.,c[i]-a[0,i]-a[1,i]-a[2,i]-a[3,i]+0.4,b[3+i*4],color='w')
	plt.ion()
	plt.show()
	p1=plt.bar(range(50),c,width,color=au.blue)
	p2=plt.bar(range(50),c-a[0,:],width,color=au.red)
	p3=plt.bar(range(50),c-a[0,:]-a[1,:],width,color=au.gree)
	p4=plt.bar(range(50),c-a[0,:]-a[1,:]-a[2,:],width,color=au.oran)
	p5=plt.bar(range(50),c-a[0,:]-a[1,:]-a[2,:]-a[3,:],width,color=au.purp)
	plt.xticks(n+width/2.,E,rotation=90)
	plt.yticks(np.arange(0,101,5))
	subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=None)
	plt.legend( (p1[0], p2[0],p3[0],p4[0],p5[0]), ('Primary user', 'Secondary user','Tertiary user','Quaternary user', 'The rest'), 'upper center',fancybox=True,prop={'size':10})
	if alph!=None:		
		savefig("./figures/new_square_columndiagram_alpha="+str(alph)+".png")
	else:
		savefig("./figures/new_square_columndiagram.png")
	
	
def track_new_exports(N,F,admat='./settings/eadmat.txt',lapse=None):
	
	matr=np.genfromtxt(admat)

	if lapse==None:
		lapse=N[0].nhours

	a,b,list_F=au.AtoKh_old(N)
	start=time()
	for n in N:
		n.links = get_links(n.id,matr)
		n.power_mix_ex = np.zeros((len(N),lapse))
    
	for t in range(lapse):
		R=np.zeros(( len(N),len(N) ))
        
        ## Update progress bar.
		if mod(t,100)==0 and t>0: 
			print "\r",round(100.0*(t/float(lapse)),2),"%",
			sys.stdout.flush()
        
		for n in N:
            
            ## Add the nodes own export/sink strength.
			n.power_mix_ex[n.id,t] = n.load[t]
            
			for l in n.links:
                
                ## If the flow direction indicates export from node n on link l.
				if F[l[0],t]*l[1] > 0:
                    
                    ## Determine end of link using the link direction.
					if l[1]==-1:
						friend_label = list_F[l[0]][0]
					elif l[1]==1:
						friend_label= list_F[l[0]][1]
					else:
						print "Warning (234nsd23): Link direction unknown!"
                    
					for k in N:
						if k.label==friend_label: friend_id=k.id  ## lazy and messy way of getting id
                    # Holger Bech Nielsen
					n.power_mix_ex[friend_id,t] = abs(F[l[0],t])
                    
            ### So now we have the own+direct export vector
            ### we add it to the big matrix
			R[n.id]=n.power_mix_ex[:,t]
        
        ### Now we have the pre-done matrix
		target=np.arange(len(N))
		done=[]
        
		while len(done)<len(N): #All is not done yet.
			for n in target:
				if n not in done:  ### if it hasn't been done
					takers=[]
					for i in target:
						if R[n,i]>1. and n!=i: takers.append(i)
						#if n!=i: contributors.append(i)
					takers.sort()
                    ## Appears not to be used: contr=np.array(contributors)
                    
					if (np.in1d(takers,done).all()) or (len(takers)==0): ### check if it is doable all it's takers are done or takers is empty
						for c in takers:                           ### then, for all "done" takers, do the thingy
							export_from_ntoc = R[n,c]*1.0
							R[n,c] = 0
							
							R[n,:]   += R[c,:]*export_from_ntoc/sum(R[c,:])
							
						done.append(n)
						done.sort()
                        
		R[np.isnan(R)] = 0.0 #This should not be needed, but R says it is.
		 
		for n in N:
			n.power_mix_ex[:,t] = R[n.id,:]   
			
			          

	return N    

def track_new_imports(N,F,admat='./settings/eadmat.txt',lapse=None):
	
	matr=np.genfromtxt(admat)

	if lapse==None:
		lapse=N[0].nhours

	a,b,list_F=au.AtoKh_old(N)
	start=time()
	for n in N:
		n.links = get_links(n.id,matr)
		n.power_mix = np.zeros((len(N),lapse))
    
	for t in range(lapse):
		R=np.zeros(( len(N),len(N) ))
        
        ## Update progress bar.
		if mod(t,100)==0 and t>0: 
			print "\r",round(100.0*(t/float(lapse)),2),"%",
			sys.stdout.flush()
        
		for n in N:
            
            ## Add the nodes own contribution/source strength.
			n.power_mix[n.id,t] = n.get_RES()[t]+n.balancing[t]#-n.curtailment[t]
            
			for l in n.links:
                
                ## If the flow direction indicates import to node n on link l.
				if F[l[0],t]*l[1] < 0:
                    
                    ## Determine end of link using the link direction.
					if l[1]==-1:
						friend_label = list_F[l[0]][0]
					elif l[1]==1:
						friend_label= list_F[l[0]][1]
					else:
						print "Warning (234nsd23): Link direction unknown!"
                    
					for k in N:
						if k.label==friend_label: friend_id=k.id  ## lazy and messy way of getting id
                    # Holger Bech Nielsen
					n.power_mix[friend_id,t] = abs(F[l[0],t])
                    
            ### So now we have the own+direct input vector
            ### we add it to the big matrix
			R[n.id]=n.power_mix[:,t]
        
        ### Now we have the pre-done matrix
		target=np.arange(len(N))
		done=[]
        
		while len(done)<len(N): #All is not done yet.
			for n in target:
				if n not in done:  ### if it hasn't been done
					contributors=[]
					for i in target:
						if R[n,i]>1. and n!=i: contributors.append(i)
						
					contributors.sort()
                    ## Appears not to be used: contr=np.array(contributors)
                    
					if (np.in1d(contributors,done).all()) or (len(contributors)==0): ### check if it is doable all it's contributors are done or contributors is empty
						for c in contributors:                           ### then, for all "done" contributors, do the thingy
							supply_from_cton = R[n,c]*1.0
							R[n,c] = 0
							
							R[n,:]   += R[c,:]*supply_from_cton/sum(R[c,:])
							
						done.append(n)
						done.sort()
                        
		R[np.isnan(R)] = 0.0 #This should not be needed, but R says it is.
		 
		for n in N:
			n.power_mix[:,t] = R[n.id,:]   
			
			          

	return N
