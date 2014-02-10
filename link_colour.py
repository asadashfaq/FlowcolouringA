# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:28:29 2013

@author: anders
"""

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
from new_linkcolour_algorithm import *
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
	#print matr
	#print len(matr[29])
	for j in range(len(matr[node_id])):
		for i in range(len(matr[node_id])):
			if matr[i][j]>0 and i>j:
				linkno+=1
				if node_id==i: links.append([linkno,-1])
				if node_id==j: links.append([linkno,1])
	return links
    
def get_link_direction(n,N):
	a,b,c,d,e=au.AtoKh(EU_Nodes()) 
	
	for k in N:
		if str(k.label)==e[n][0][0:2]: start_node=k
	for k in N:
		if str(k.label)==e[n][0][6:8]: end_node=k
	
	return [start_node,end_node]
	
		
def track_node_links_usage(links,linksflow,N=None,ports="ex+im",new=False,mode="linear",alph=None,copper = True):
# as links use links_mix or variants of this, as linksflow use links_flow or variants of this, set new=True if you use the new model
#gives the linkusage for each country both relative to its total usage of links and to the links total flow. 
#Mode can be "linear", "square" or "random".
	if alph== None:
	    alph="hetero"
	if N==None:
		N=EU_Nodes()
	x=len(N)
	
	a,b,c,d,e=au.AtoKh(EU_Nodes())
	total_link_usage_node=[]
	total_link_usage_node_rel=[]
	total_usages=[]
	#cost=[]
	
	for n in N:
		total_usage=0
		n.link_usages=[]
		n.link_usages_normalized=[]
		
		for j in range (len(e)):
			n.link_usages=np.append(n.link_usages,links[n.id+x*j])
			total_usage+=links[n.id+x*j]
			k=linksflow[j]
			#cost=np.append(cost,links[n.id+x*j]/k)
		n.link_usages_normalized=np.append(n.link_usages_normalized,n.link_usages[:]/total_usage*100)
		total_link_usage_node_rel=np.append(total_link_usage_node_rel,n.link_usages_normalized)
		total_link_usage_node=np.append(total_link_usage_node,n.link_usages)
		total_usages.append(total_usage)
		
		
		print n.link_usages_normalized
		print sum(n.link_usages_normalized)
		#print cost
	if copper:
		if mode=="square":
			if new:
				np.save('linkcolouring/new_squared_copper_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel) 
				np.save('linkcolouring/new_squared_copper_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
				#np.save('linkcolouring/new_squared_cost_links_alpha='+str(alph),cost)
			else:
				np.save('linkcolouring/old_node_squared_copper_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel)
				np.save('linkcolouring/old_node_squared_copper_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
		elif mode == "linear":
			if new:
				np.save('linkcolouring/new_linear_copper_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel) 
				np.save('linkcolouring/new_linear_copper_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
				#np.save('linkcolouring/new_cost_links_alpha='+str(alph),cost)
			else:
				np.save('linkcolouring/old_node_linear_copper_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel)
				np.save('linkcolouring/old_node_linear_copper_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
				np.save('linkcolouring/old_node_linear_copper_links_total_usage_'+str(ports)+'_alpha='+str(alph),total_usages)
		elif mode == "random":
			if new:
				np.save('linkcolouring/new_random_copper_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel) 
				np.save('linkcolouring/new_random_copper_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
				#np.save('linkcolouring/new_cost_links_alpha='+str(alph),cost)
			else:
				np.save('linkcolouring/old_random_copper_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel)
				np.save('linkcolouring/old_random_copper_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
	else:
		if mode=="square":
			if new:
				np.save('linkcolouring/new_squared_constr_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel) 
				np.save('linkcolouring/new_squared_constr_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
				#np.save('linkcolouring/new_squared_cost_links_alpha='+str(alph),cost)
			else:
				np.save('linkcolouring/old_node_squared_constr_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel)
				np.save('linkcolouring/old_node_squared_constr_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
		elif mode == "linear":
			if new:
				np.save('linkcolouring/new_linear_constr_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel) 
				np.save('linkcolouring/new_linear_constr_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
				#np.save('linkcolouring/new_cost_links_alpha='+str(alph),cost)
			else:
				np.save('linkcolouring/old_node_linear_constr_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel)
				np.save('linkcolouring/old_node_linear_constr_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
		elif mode == "random":
			if new:
				np.save('linkcolouring/new_random_constr_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel) 
				np.save('linkcolouring/new_random_constr_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
				#np.save('linkcolouring/new_cost_links_alpha='+str(alph),cost)
			else:
				np.save('linkcolouring/old_random_constr_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel)
				np.save('linkcolouring/old_random_constr_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
			
	
	
	
def track_link_usage_total_link(N,F,n,lapse=None):
	
	if lapse== None:
		lapse=N[0].nhours   
	file = open('link_usage.dat', 'w+')
	matr=np.genfromtxt(N.pathadmat)       
	a,b,c,d,e=au.AtoKh(EU_Nodes())  
	f,g,list_F=au.AtoKh_old(N)	
	 
	total_flow_link=0 
	total_link_mix=np.zeros((len(N)))  		
	link_mix_export=np.zeros((len(N),lapse))		
	link_mix_import=np.zeros((len(N),lapse))	
	total_link_mix_percentage=np.zeros((len(N)))

	start=get_link_direction(n,N)[0]
	end=get_link_direction(n,N)[1]
	start.links=get_links(start.id,matr)
	
	for t in range(lapse):
	## Update progress bar.
		if mod(t,100)==0 and t>0: 
			print "\r",round(100.0*(t/float(lapse)),2),"%",
			sys.stdout.flush()
		
		 
    
		for l in start.links:
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
				
			for k in N:
						
				if k.label==sending_node_label: 
					sending_node=k
				elif k.label==receiving_node_label:
					receiving_node=k
									
			if (sending_node==start and receiving_node==end) or (sending_node==end and receiving_node==start):
				link_mix_import[:,t]=receiving_node.power_mix_ex[:,t]*abs(F[l[0],t])/(sum(receiving_node.power_mix_ex[:,t])+0.000000000001)
				link_mix_export[:,t]=sending_node.power_mix[:,t]*abs(F[l[0],t])/(sum(sending_node.power_mix[:,t])+0.000000000001)
					
														
				total_flow_link+=abs(F[l[0],t])
				total_link_mix[:] += (link_mix_export[:,t]+link_mix_import[:,t])/2					
				   	    
	total_link_mix_percentage[:]=total_link_mix[:]/total_flow_link*100		
	print sum(total_link_mix_percentage[:])
	
	
	for k in N:
		k.linkmix=[]
		k.linkmix=np.append(k.linkmix,total_link_mix_percentage[k.id])
		x="country", str(k.label), "has a total usage of link",str(e[n]) , "of" ,str(total_link_mix_percentage[k.id]), "% after", str(t), "hours"
		file.write(str(x) + '\n')
	
	#b=[]
	#for l in range(len(N)):
		#b=np.append(b,N[l].linkmix)
		
	#c=sorted(b,reverse=True)
	#d=np.zeros(len(N))
	#for l in range(len(N)):
		#for m in range(len(N)):
			#if c[l]==b[m]:
				#d[l]=m
	#return c,d
	return N

		
def track_link_usage_total_old(N,F,lapse=None):
	
	if lapse== None:
		lapse=N[0].nhours   
	
	file = open('link_usage_total.dat', 'w+')
	
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
		link_mix_import=np.zeros((len(N),lapse))		
		link_mix_export=np.zeros((len(N),lapse))	
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
		cost=np.append(cost,total_link_mix_percentage)
		print sum(total_link_mix_percentage[:])		
		links=np.append(links,total_link_mix[:])
		print sum(total_link_mix[:])	
		print total_flow_link
		linksflow=np.append(linksflow,total_flow_link)		
		
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
	np.save('linkcolouring/cost',cost)
	np.save('linkcolouring/boxplot',boxplot) 
	np.save('linkcolouring/boksplotlabel',boksplotlabel)
	np.save('linkcolouring/links_mix',links)
	np.save('linkcolouring/links_flow',linksflow)
	
	
	return boxplot,boksplotlabel
	
def boksplot_old(a,b): #insert boxplot and boksplot.label
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
	savefig("./figures/columndiagram.pdf")
	
def track_link_usage_imports(N,F,n,lapse=None):
	
	if lapse== None:
		lapse=N[0].nhours    
    
	matr=np.genfromtxt(N.pathadmat)
       
	a,b,c,d,e=au.AtoKh(EU_Nodes())  
    
	f,g,list_F=au.AtoKh_old(N)
	
   
	total_flow_link=0    
	total_link_mix_import=np.zeros((len(N)))
	link_mix_import=np.zeros((len(N),lapse))
	total_link_mix_percentage_import=np.zeros((len(N)))
          
	start=get_link_direction(n,N)[0]
	end=get_link_direction(n,N)[1]			
    
	start.links=get_links(start.id,matr)
    
    
       
	for t in range(lapse):
		for l in start.links:
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
            #if flow indicates import to end_node
				else:
				
					if l[1]==-1:
						sending_node_label=list_F[l[0]][0]
						receiving_node_label=list_F[l[0]][1]
					elif l[1]==1:
						sending_node_label=list_F[l[0]][1] 
						receiving_node_label=list_F[l[0]][0]
					else:
						print "Warning (234nsd23): Link direction unknown!"
					
				for k in N:
								
					if k.label==sending_node_label: 
						sending_node=k
						#sending_id_node=k
					elif k.label==receiving_node_label:
						receiving_node=k
						#receiving_id_node=k
			
			#if (sending_id==get_link_direction(n,N)[0] and receiving_id==get_link_direction(n,N)[1]) or (sending_id==get_link_direction(n,N)[1] and receiving_id==get_link_direction(n,N)[0]):	
				if (sending_node==start and receiving_node==end) or (sending_node==end and receiving_node==start):
					link_mix_import[:,t]=receiving_node.power_mix_ex[:,t]*abs(F[l[0],t])/(sum(receiving_node.power_mix_ex[:,t])+0.000000000001)
					
					total_link_mix_import[:] += link_mix_import[:,t]						
					total_flow_link+=abs(F[l[0],t])			
											
           
	total_link_mix_percentage_import[:]=total_link_mix_import[:]/total_flow_link*100
	print sum(total_link_mix_percentage_import[:])
	
	for k in N:
		print "country", k.label ,"has a total usage of link", e[n] ,"of", total_link_mix_percentage_import[k.id] ,"% after ",t," hours"
							
	
	
def track_link_usage_exports(N,F,n,lapse=None):
	
	if lapse== None:
		lapse=N[0].nhours    
    
	matr=np.genfromtxt(N.pathadmat)
       
	a,b,c,d,e=au.AtoKh(EU_Nodes())  
    
	f,g,list_F=au.AtoKh_old(N)
	
   
	total_flow_link=0    
	total_link_mix_export=np.zeros((len(N)))
	link_mix_export=np.zeros((len(N),lapse))
	total_link_mix_percentage_export=np.zeros((len(N)))
    
        
	start=get_link_direction(n,N)[0]
	end=get_link_direction(n,N)[1]
				
    
    #start.links=get_links(get_link_direction(n,N)[0],matr)
	start.links=get_links(start.id,matr)
    
    
       
	for t in range(lapse):
		for l in start.links:
			#if flow indicates import to start_node
			#if abs(F[l[0],t])>27:
				if F[l[0],t]*l[1]<0:
					if l[1]==-1:
						sending_node_label=list_F[l[0]][0]
						receiving_node_label=list_F[l[0]][1]
					elif l[1]==1:
						sending_node_label=list_F[l[0]][1] 
						receiving_node_label=list_F[l[0]][0]
					else:
						print "Warning (234nsd23): Link direction unknown!"
            #if flow indicates import to end_node
				else:
				
					if l[1]==1:
						sending_node_label=list_F[l[0]][0]
						receiving_node_label=list_F[l[0]][1]
					elif l[1]==-1:
						sending_node_label=list_F[l[0]][1] 
						receiving_node_label=list_F[l[0]][0]
					else:
						print "Warning (234nsd23): Link direction unknown!"
					
				for k in N:
								
					if k.label==sending_node_label: 
						sending_node=k
						
					elif k.label==receiving_node_label:
						receiving_node=k
						
			
			#if (sending_id==get_link_direction(n,N)[0] and receiving_id==get_link_direction(n,N)[1]) or (sending_id==get_link_direction(n,N)[1] and receiving_id==get_link_direction(n,N)[0]):	
				if (sending_node==start and receiving_node==end) or (sending_node==end and receiving_node==start):
					link_mix_export[:,t]=sending_node.power_mix[:,t]*abs(F[l[0],t])/(sum(sending_node.power_mix[:,t])+0.000000000001)
					
					total_link_mix_export[:] += link_mix_export[:,t]					
					total_flow_link+=abs(F[l[0],t])											
				
   
        
	total_link_mix_percentage_export[:]=total_link_mix_export[:]/total_flow_link*100
	print sum(total_link_mix_percentage_export[:])
		
    
	for k in N:
		print "country", k.label ,"has a total usage of link", e[n] ,"of", total_link_mix_percentage_export[k.id] ,"% after ",t," hours"
				
	
                   
 
    
def track_exports(N,F,admat='./settings/eadmat.txt',lapse=None):
	
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
			n.power_mix_ex[n.id,t] = n.get_import()[t]
            
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

def track_imports(N,F,admat='./settings/eadmat.txt',lapse=None):
	
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
			n.power_mix[n.id,t] = n.get_export()[t]
            
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

def storage():
    balancing=[]
    storages=[1,2,4,6,12,24]
    for s in storages:
        N=EU_Nodes(load_filename="Case_S_Beta_"+str(s)+".npz")
        balancing.append(sum(sum(n.balancing for n in N)))
    return balancing,storages

def make_table():
    Unc=get_quant(quant=1)/1000
    Q99=get_quant()/1000
    Int=0.4*get_quant()/1000
    a,b,c=AtoKh_old(Nodes())
    Pre=b[0:88]/1000
    for i in c:
        x=i[2]
        if mod(x,2)==0:
            print " & ",round(Pre[2*x],2)," & ",round(Int[2*x],2)," & ",round(Q99[2*x],2)," & ",round(Unc[2*x],2),"\\\\"
            print "\multirow{-2}{*}{",i[0] ,r" $\rightleftarrows$ ",i[1] ,"} & ",round(Pre[2*x+1],2)," & ",round(Int[2*x+1],2)," & ",round(Q99[2*x+1],2)," & ",round(Unc[2*x+1],2),r"\\"
        if mod(x,2)==1:
            print r"\rowcoor{rgrey} & ",round(Pre[2*x],2)," & ",round(Int[2*x],2)," & ",round(Q99[2*x],2)," & ",round(Unc[2*x],2),"\\\\"
            print r"\rowcolor{rgrey} \multirow{-2}{*}{",i[0] ,r" $\rightleftarrows$ ",i[1] ,"} & ",round(Pre[2*x+1],2)," & ",round(Int[2*x+1],2)," & ",round(Q99[2*x+1],2)," & ",round(Unc[2*x+1],2),r"\\"


shortfiles=['ISET_country_AT.npz', 'ISET_country_FI.npz', 'ISET_country_NL.npz', 'ISET_country_BA.npz', 'ISET_country_FR.npz', 'ISET_country_NO.npz', 'ISET_country_BE.npz', 'ISET_country_GB.npz', 'ISET_country_PL.npz', 'ISET_country_BG.npz', 'ISET_country_GR.npz', 'ISET_country_PT.npz', 'ISET_country_CH.npz', 'ISET_country_HR.npz', 'ISET_country_RO.npz', 'ISET_country_CZ.npz', 'ISET_country_HU.npz', 'ISET_country_RS.npz', 'ISET_country_DE.npz', 'ISET_country_IE.npz', 'ISET_country_SE.npz', 'ISET_country_DK.npz', 'ISET_country_IT.npz', 'ISET_country_SI.npz', 'ISET_country_ES.npz', 'ISET_country_LU.npz', 'ISET_country_SK.npz','ISET_country_LU.npz', 'ISET_country_SK.npz', 'ISET_country_SK.npz']
longfiles=['ISET_country_AT.npz', 'ISET_country_FI.npz', 'ISET_country_NL.npz', 'ISET_country_BA.npz', 'ISET_country_FR.npz', 'ISET_country_NO.npz', 'ISET_country_BE.npz', 'ISET_country_GB.npz', 'ISET_country_PL.npz', 'ISET_country_BG.npz', 'ISET_country_GR.npz', 'ISET_country_PT.npz', 'ISET_country_CH.npz', 'ISET_country_HR.npz', 'ISET_country_RO.npz', 'ISET_country_CZ.npz', 'ISET_country_HU.npz', 'ISET_country_RS.npz', 'ISET_country_DE.npz', 'ISET_country_IE.npz', 'ISET_country_SE.npz', 'ISET_country_DK.npz', 'ISET_country_IT.npz', 'ISET_country_SI.npz', 'ISET_country_ES.npz', 'ISET_country_LU.npz', 'ISET_country_SK.npz','ISET_country_LU.npz', 'ISET_country_SK.npz','ISET_country_LU.npz', 'ISET_country_SK.npz', 'ISET_country_SK.npz']
