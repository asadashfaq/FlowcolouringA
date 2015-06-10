#! /usr/bin/env/python
import csv
import sys
from pylab import *
from time import time
from EUgrid import *
import networkx as nx
import link_namer as ln

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
	a,b,c,d,e=au.AtoKh(N) 
	
	for k in N:
		if str(k.label)==e[n][0][0:2]: start_node=k
	for k in N:
		if str(k.label)==e[n][0][6:8]: end_node=k
	
	return [start_node,end_node]


def track_link_usage_total(N,F,new=False,lapse=None,alph=None,mode="linear",copper=True,sensitivity=False,constrained=False,logistic=False,gamma=False,heterogen=False):
#tracks the usage of each link for each country, 
#alph sets the value of homogenous alphas.
#Mode can be "linear", "square", "random" or "capped".
#set new=True if new model is used.
	if alph==None:
	    alph="hetero"
	if lapse== None:
		lapse=N[0].nhours 
	if new:
	    if copper:
		file=open('new_'+str(mode)+'_copper_link_usage_total_alpha='+str(alph)+'.dat','w+')
		
	
	    else:
		file = open('new_'+str(mode)+'_constr_link_usage_total_alpha='+str(alph)+'.dat', 'w+')
	    
	else:
	
		
	    if copper:
		file=open('old_'+str(mode)+'_copper_link_usage_total_alpha='+str(alph)+'.dat','w+')
		
	
	    else:
		 file=open('old_'+str(mode)+'_constr_link_usage_total_alpha='+str(alph)+'.dat','w+')
		
	# Check what kind of network is being used for naming savefiles see line 273
	if len(N) == 30:
		network = ''
	elif len(N) == 8:
		network = 'superRegions'
	else:
		network = 'regions'


	matr=np.genfromtxt(N.pathadmat)
       
	a,b,c,d,e=au.AtoKh(N)     
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
		linkDirection = ln.link_direction(n,N)
		start = linkDirection[0]
		end = linkDirection[1]
        
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
			
	    
		# sorting the total_link_mix_percentage, such that the higher percentages are first
		a=sorted(total_link_mix_percentage,reverse=True)
		#for the 4 highest values we save the usage for each link and the corresponding countries
		print a
		for i in range(4):
			boxplot[i,n]=a[i]
			for h in N:
				
				if a[i]==total_link_mix_percentage[h.id]:
				    boksplotlabel=np.append(boksplotlabel,str(h.label))
				    
	numNodes = len(N)
	numLinks = len(F)
	maxlex=np.array(maxlex).reshape((numLinks,numNodes))
	maxlim=np.array(maxlim).reshape((numLinks,numNodes))
	minlex=np.array(minlex).reshape((numLinks,numNodes))
	minlim=np.array(minlim).reshape((numLinks,numNodes))
	meanlex=np.array(meanlex).reshape((numLinks,numNodes))
	meanlim=np.array(meanlim).reshape((numLinks,numNodes))
	quant99lex=np.array(quant99lex).reshape((numLinks,numNodes))
	quant99lim=np.array(quant99lim).reshape((numLinks,numNodes))
	
	if new:
	    if copper:
		    np.save('linkcolouring/new_'+str(mode)+'_copper_boxplot_alpha='+str(alph),boxplot) 
		    np.save('linkcolouring/new_'+str(mode)+'_copper_boksplotlabel_alpha='+str(alph),boksplotlabel)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_links_mix_alpha='+str(alph),links)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_links_flow_alpha='+str(alph),linksflow)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_cost_alpha='+str(alph),cost)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_links_ex_alpha='+str(alph),links_ex)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_links_im_alpha='+str(alph),links_im)		   
		    np.save('linkcolouring/new_'+str(mode)+'_copper_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_maxlex_alpha='+str(alph),maxlex)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_maxlim_alpha='+str(alph),maxlim)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_minlex_alpha='+str(alph),minlex)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_minlim_alpha='+str(alph),minlim)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_meanlex_alpha='+str(alph),meanlex)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_meanlim_alpha='+str(alph),meanlim)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_quant99lex_alpha='+str(alph),quant99lex)
		    np.save('linkcolouring/new_'+str(mode)+'_copper_quant99lim_alpha='+str(alph),quant99lim)
		
	    else:
		    np.save('linkcolouring/new_'+str(mode)+'_constr_boxplot_alpha='+str(alph),boxplot) 
		    np.save('linkcolouring/new_'+str(mode)+'_constr_boksplotlabel_alpha='+str(alph),boksplotlabel)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_links_mix_alpha='+str(alph),links)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_links_flow_alpha='+str(alph),linksflow)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_cost_alpha='+str(alph),cost)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_links_ex_alpha='+str(alph),links_ex)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_links_im_alpha='+str(alph),links_im)		   
		    np.save('linkcolouring/new_'+str(mode)+'_constr_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_maxlex_alpha='+str(alph),maxlex)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_maxlim_alpha='+str(alph),maxlim)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_minlex_alpha='+str(alph),minlex)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_minlim_alpha='+str(alph),minlim)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_meanlex_alpha='+str(alph),meanlex)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_meanlim_alpha='+str(alph),meanlim)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_quant99lex_alpha='+str(alph),quant99lex)
		    np.save('linkcolouring/new_'+str(mode)+'_constr_quant99lim_alpha='+str(alph),quant99lim)
		
		

	else:
	    if copper:
			if sensitivity:
				np.save('./sensitivity/linkcolouring/'+network+'-old_'+str(mode)+'_copper_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
				np.save('./sensitivity/linkcolouring/'+network+'-old_'+str(mode)+'_copper_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
			elif constrained:
				np.save('./constrained/linkcolouring/'+str(mode)+'_link_mix_import_b_'+str(constrained),link_mix_import_all)
				np.save('./constrained/linkcolouring/'+str(mode)+'_link_mix_export_b_'+str(constrained),link_mix_export_all)
			elif logistic:
				np.save('./results/logistic/linkcolouring/'+str(mode)+'-'+str(logistic)+'_link_mix_import',link_mix_import_all)
				np.save('./results/logistic/linkcolouring/'+str(mode)+'-'+str(logistic)+'_link_mix_export',link_mix_export_all)
			elif gamma:
				np.save('./results/logistic/linkcolouring/'+str(mode)+'-g-'+str(gamma)+'_link_mix_import',link_mix_import_all)
				np.save('./results/logistic/linkcolouring/'+str(mode)+'-g-'+str(gamma)+'_link_mix_export',link_mix_export_all)
			elif type(heterogen) != bool:
				np.save('./linkcolouring/heterogen/'+str(mode)+'-b-'+str(heterogen)+'_link_mix_import',link_mix_import_all)
				np.save('./linkcolouring/heterogen/'+str(mode)+'-b-'+str(heterogen)+'_link_mix_export',link_mix_export_all)
			else:
				np.save('linkcolouring/old_'+str(mode)+'_copper_boxplot_alpha='+str(alph),boxplot) 
				np.save('linkcolouring/old_'+str(mode)+'_copper_boksplotlabel_alpha='+str(alph),boksplotlabel)
				np.save('linkcolouring/old_'+str(mode)+'_copper_links_mix_alpha='+str(alph),links)
				np.save('linkcolouring/old_'+str(mode)+'_copper_links_flow_alpha='+str(alph),linksflow)
				np.save('linkcolouring/old_'+str(mode)+'_copper_cost_alpha='+str(alph),cost)
				np.save('linkcolouring/old_'+str(mode)+'_copper_links_ex_alpha='+str(alph),links_ex)
				np.save('linkcolouring/old_'+str(mode)+'_copper_links_im_alpha='+str(alph),links_im)		
				np.save('linkcolouring/old_'+str(mode)+'_copper_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
				np.save('linkcolouring/old_'+str(mode)+'_copper_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
				np.save('linkcolouring/old_'+str(mode)+'_copper_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
				np.save('linkcolouring/old_'+str(mode)+'_copper_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
				np.save('linkcolouring/old_'+str(mode)+'_copper_maxlex_alpha='+str(alph),maxlex)
				np.save('linkcolouring/old_'+str(mode)+'_copper_maxlim_alpha='+str(alph),maxlim)
				np.save('linkcolouring/old_'+str(mode)+'_copper_minlex_alpha='+str(alph),minlex)
				np.save('linkcolouring/old_'+str(mode)+'_copper_minlim_alpha='+str(alph),minlim)
				np.save('linkcolouring/old_'+str(mode)+'_copper_meanlex_alpha='+str(alph),meanlex)
				np.save('linkcolouring/old_'+str(mode)+'_copper_meanlim_alpha='+str(alph),meanlim)
				np.save('linkcolouring/old_'+str(mode)+'_copper_quant99lex_alpha='+str(alph),quant99lex)
				np.save('linkcolouring/old_'+str(mode)+'_copper_quant99lim_alpha='+str(alph),quant99lim)
			
		
	    else:
		
		    
		    np.save('linkcolouring/old_'+str(mode)+'_constr_boxplot_alpha='+str(alph),boxplot) 
		    np.save('linkcolouring/old_'+str(mode)+'_constr_boksplotlabel_alpha='+str(alph),boksplotlabel)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_links_mix_alpha='+str(alph),links)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_links_flow_alpha='+str(alph),linksflow)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_cost_alpha='+str(alph),cost)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_links_ex_alpha='+str(alph),links_ex)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_links_im_alpha='+str(alph),links_im)		
		    np.save('linkcolouring/old_'+str(mode)+'_constr_link_mix_import_all_alpha='+str(alph),link_mix_import_all)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_link_mix_export_all_alpha='+str(alph),link_mix_export_all)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_link_mix_export_variance_all_alpha='+str(alph),link_mix_export_variance_all)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_link_mix_import_variance_all_alpha='+str(alph),link_mix_import_variance_all)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_maxlex_alpha='+str(alph),maxlex)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_maxlim_alpha='+str(alph),maxlim)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_minlex_alpha='+str(alph),minlex)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_minlim_alpha='+str(alph),minlim)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_meanlex_alpha='+str(alph),meanlex)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_meanlim_alpha='+str(alph),meanlim)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_quant99lex_alpha='+str(alph),quant99lex)
		    np.save('linkcolouring/old_'+str(mode)+'_constr_quant99lim_alpha='+str(alph),quant99lim)
		    
		

		
	
	return boxplot,boksplotlabel
	
def boksplot(a,b,alph=None,new=False,mode="linear",copper=True): # draws a columndiagram displaying the 4 biggest users and their usage of each link. 
#Insert new_boxplot and new_boksplotlabel or variants of these.
#a is an array representing the portion of usage of the different links for the 4 biggest users at each link.
#b is an array holding the corresponding labels of these different users.
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
		
	
	
		
	    

def costplot(links,linksflow,N,F,Graf=None,alph=None,new=False,mode="linear",copper=True,ports="ex_and_im",load=True): # prints different columndiagrams and saves values in arrays
#Mode can be "linear", "square", "random" or "capped".
#links is the links_mix or variations of that
#set load=True for a graph that sorts countrys by their average load.
	if N==None:
	    N=EU_Nodes()
	if alph==None:
	    alph="hetero"		
	d,f,g,r,e=au.AtoKh(EU_Nodes())	
	if new:
	    file=open('new_'+str(mode)+'_copper_costplot_alpha='+str(alph)+'.dat','w+')
	else:
	    file=open('old_'+str(mode)+'_copper_costplot_alpha='+str(alph)+'.dat','w+')
	
	capacities_today=[]
	capacities_year=[]
	width=1.0	
	gyell="#CCFF00"
	E = [n.label for n in N]
	D = [n.label for n in N]
	C = [n.label for n in N]
	fig=plt.figure(figsize=(13,7))
	blue_div_yellow_kms=np.zeros(len(N))
	blue_div_yellow_cost=np.zeros(len(N))
	blue_div_yellow=np.zeros(len(N))
	red_div_yellow_kms=np.zeros(len(N))
	red_div_yellow_cost=np.zeros(len(N))
	red_div_yellow=np.zeros(len(N))
	green_div_yellow_kms=np.zeros(len(N))
	green_div_yellow_cost=np.zeros(len(N))
	green_div_yellow=np.zeros(len(N))
	blue=np.zeros(len(N))
	Blue_GWKms=np.zeros(len(N))
	Green_GWkms=np.zeros(len(N))
	Blue_cost=np.zeros(len(N))
	Green_cost=np.zeros(len(N))
	green=np.zeros(len(N))
	c=np.zeros(len(N))
	c_kms=np.zeros(len(N))
	c_kms_avg=np.zeros(len(N))
	c_cost=np.zeros(len(N))
	c_cost_final=np.zeros(len(N))
	c_kms_avg_final=np.zeros(len(N))
	d=np.zeros(len(N))
	d_cost=np.zeros(len(N))
	cost=np.zeros(len(N))
	Graford=np.zeros(len(N))
	t=np.zeros(len(N))
	t_kms=np.zeros(len(N))
	t_cost=np.zeros(len(N))	
	today_usage=np.zeros(len(N))
	usage_year=np.zeros(len(N))
	usage_year_kms=np.zeros(len(N))
	usage_year_cost=np.zeros(len(N))
	today_usage_kms=np.zeros(len(N))
	total_usage=np.zeros(len(N))
	total_usage_year=np.zeros(len(N))
	total_usage_kms_year=np.zeros(len(N))
	total_usage_cost_year=np.zeros(len(N))
	total_usage_kms=np.zeros(len(N))
	total_usage_cost=np.zeros(len(N))
	today_usage_cost=np.zeros(len(N))
	cost_rev=np.zeros(len(N))
	usage_rel_to_share_of_load=np.zeros(len(N))
	usage_rel_to_share_of_load_kms=np.zeros(len(N))
	usage_rel_to_share_of_load_cost=np.zeros(len(N))
	j=np.arange(len(N))
	total_mean_load=0
	total_capacity=0
	total_capacity_kms=0
	total_capacities_2050_cost=0
	today_cap=0
	today_cap_kms=0
	today_cap_cost=0
	counter=0
	loads=[]
	tot_usage=0
	tot_usage_rel=0
	link_kms=[686,250,217,523,766,278,397,82,914,175,358,577,289,201,262,342,436,879,1109,1053,418,483,464,520,520,808,534,524,298,328,1052,503,753,691,305,371,117,642,446,281,291,317,163,812,356,591,523,491,277,263]
	link_prizes=[400,400,400,400,400,400,1500,1500,1500,400,1500,400,400,400,400,1500,400,400,400,400,400,1500,1500,400,400,1500,400,400,400,400,1500,400,400,400,400,400,400,400,400,400,400,400,400,1500,400,400,1500,400,400,400]
	
	link_kms_average=sum(link_kms)/len(link_kms)
	
	
	for i in range(len(e)):
	    capacities_today.append(max(r[i*2:i*2+2]))
	for l in range(len(F)):
	    capacities_year.append(max(abs(np.percentile(F[l],1)),abs(np.percentile(F[l],99))))
	    total_capacity+=max(abs(np.percentile(F[l],1)),abs(np.percentile(F[l],99)))
	    total_capacity_kms+=max(abs(np.percentile(F[l],1)),abs(np.percentile(F[l],99)))*link_kms[l]
	    if link_prizes[l]==1500: #If DC then a converter is needed at every end of the cable/HV-line
		total_capacities_2050_cost+=max(abs(np.percentile(F[l],1)),abs(np.percentile(F[l],99)))*150000*2+max(abs(np.percentile(F[l],1)),abs(np.percentile(F[l],99)))*link_kms[l]*link_prizes[l]
	    else: 
		total_capacities_2050_cost+=max(abs(np.percentile(F[l],1)),abs(np.percentile(F[l],99)))*link_kms[l]*link_prizes[l]
	Capacities_year=np.array(capacities_year)
	for n in N:
	    total_mean_load+=n.mean
	    loads.append(n.mean)
	loadssort=sorted(loads,reverse=True)  	   
	
		
	    
	for n in N:
	    cap=0	    
	    attached_links=get_links(n.id)
	    today_cap=0
	    today_cap_kms=0
	    today_cap_cost=0
	    cap_year=0
	    cap_year_kms=0
	    cap_cost_year=0
	    for l in attached_links:
		cap_year+=capacities_year[l[0]]
		cap_year_kms+=capacities_year[l[0]]*link_kms[l[0]]
		today_cap+=capacities_today[l[0]]
		today_cap_kms+=capacities_today[l[0]]*link_kms[l[0]]
		
		if link_prizes[l[0]]==1500: #If link is DC and two converter stations are needed
		    today_cap_cost+=capacities_today[l[0]]*link_kms[l[0]]*link_prizes[l[0]]+capacities_today[l[0]]*150000*2
		    cap_cost_year+=capacities_year[l[0]]*link_kms[l[0]]*link_prizes[l[0]]+capacities_year[l[0]]*150000*2
		else:    
		    today_cap_cost+=capacities_today[l[0]]*link_kms[l[0]]*link_prizes[l[0]]
		    cap_cost_year+=capacities_year[l[0]]*link_kms[l[0]]*link_prizes[l[0]]
	
	    today_usage[n.id]=today_cap*0.5 #red in "MW"
	    today_usage_kms[n.id]=today_cap_kms*0.5 #red in "MWKm"
	    today_usage_cost[n.id]=today_cap_cost*0.5 #red in Euros
	    usage_year[n.id]=cap_year*0.5 #green in "MW"
	    usage_year_kms[n.id]=cap_year_kms*0.5 #green in "MWKm"
	    usage_year_cost[n.id]=cap_cost_year*0.5 #green in Euros
	  
	    
	
	for n in N:	    
	    for l in range(len(F)):
		c[n.id]+=links[n.id+30*l]/linksflow[l]*max(abs(np.percentile(F[l],1)),abs(np.percentile(F[l],99)))
		c_kms[n.id]+=links[n.id+30*l]/linksflow[l]*max(abs(np.percentile(F[l],1)),abs(np.percentile(F[l],99)))*link_kms[l]		
	    t[n.id]=n.mean*total_capacity/total_mean_load #yellow in "MW"
	    t_kms[n.id]=n.mean*total_capacity_kms/total_mean_load #yellow in "MWKm"
	    t_cost[n.id]=n.mean*total_capacities_2050_cost/total_mean_load #yellow in "MW"
	    c_cost[n.id]=c_kms[n.id]/total_capacity_kms*total_capacities_2050_cost #blue in Euros
	    c[n.id]=c[n.id]/1000  #blue in "GW"
	    c_kms[n.id]=c_kms[n.id]/1000 #blue in "GW"
	    c_kms_avg[n.id]=c_kms[n.id]/link_kms_average
	    
	    
	k=sorted(c,reverse=True)
	k_kms=sorted(c_kms,reverse=True)
	k_cost=sorted(c_cost,reverse=True)
	
	
	if load:
	    for h in N:
		for n in N:
		    
		    if n.mean==loadssort[h.id]:		    
			E[h.id]=n.label
					    
			usage_rel_to_share_of_load[h.id]=t[n.id]/1000 #yellow in "GW"
			blue[h.id]=c[n.id]
			green[h.id]=usage_year[n.id]/1000
		    
			total_usage[h.id]=today_usage[n.id]/1000
			red_div_yellow[h.id]=total_usage[h.id]/usage_rel_to_share_of_load[h.id]
			
			total_usage_year[h.id]=usage_year[n.id]/1000
		
	    
			
		tot_usage+=k[h.id]
		tot_usage_rel+=t[h.id]
	    red_div_yellow[25]=0
	else:
	    for h in N:
		for n in N:		    
			if n.label==Grafsorted[h.id][0]:
			    E[h.id]=n.label	
			    
			    d[h.id]=c[n.id]/sum(abs(n.mismatch))*1000000000
			    usage_rel_to_share_of_load[h.id]=t[n.id]/1000
			   
			    total_usage[h.id]=today_usage[n.id]/1000
			    
			    total_usage_year[h.id]=usage_year[n.id]/1000
			    
		tot_usage+=k[h.id]
		tot_usage_rel+=t[h.id]
	   
		
	plt.ion()
	plt.show()	
	p1=plt.bar(np.arange(30),red_div_yellow,width=0.7,color=au.oran)
	plt.xlim(j[0]-0.5)	
	plt.xticks(j+0.4,E,rotation=90)
			
	
	if new:
	    if copper:
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    file.write("red_div_yellow_GW-new"+str(E)+str(usage_rel_to_share_of_load)+str(total_usage)+str(total_usage_year)+"red_div_yellow"+str(red_div_yellow) + '\n')
		    np.savez("./figures/new_GW_costplots/GW-numbers"+str(mode)+str(alph)+str(ports),labels=E, Blue=blue,Green=green,Yellow=usage_rel_to_share_of_load)
		    savefig("./figures/new_GW_costplots/red_yellow/new_"+str(mode)+"_copper_costdiagram_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_GW_costplots/red_yellow/new_"+str(mode)+"_constr_costdiagram_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	else:
	    if copper:
		    file.write("red_div_yellow_GW-old"+str(E)+str(usage_rel_to_share_of_load)+str(total_usage)+str(total_usage_year)+"red_div_yellow"+str(red_div_yellow) + '\n')
		    np.savez("./figures/old_GW_costplots/GW-numbers"+str(mode)+str(alph)+str(ports),labels=E, Blue=blue,Green=green,Yellow=usage_rel_to_share_of_load)
		    savefig("./figures/old_GW_costplots/red_yellow/old_"+str(mode)+"_copper_costdiagram_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_GW_costplots/red_yellow/old_"+str(mode)+"_constr_costdiagram_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	close()
	blue_div_green=blue/green
	blue_div_green_avg=round(blue_div_green.sum()/(len(N)),2)
	
	plt.ion()
	plt.show()	
	p1=plt.bar(np.arange(30),blue_div_green,width=0.7,color=au.cafe)
	plt.axhline(blue_div_green_avg,linestyle="dashed",color="k")
	plt.annotate(str(blue_div_green_avg),(-0.5,blue_div_green_avg),(1,blue_div_green_avg+0.5),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"))
	plt.ylabel('I-factor')
	plt.xlim(j[0]-0.5)	
	plt.xticks(j+0.4,E,rotation=90)
	
	if new:
	    if copper:
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    file.write("blue_div_green_GW-new"+str(E)+"blue_div_green"+str(blue_div_green) + '\n')
		    
		    savefig("./figures/new_GW_costplots/blue_div_green/new_"+str(mode)+"_copper_costdiagram_blue_div_green"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_GW_costplots/blue_div_green/new_"+str(mode)+"_constr_costdiagram_blue_div_green"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	else:
	    if copper:
		    file.write("red_div_yellow_GW-old"+str(E)+"blue_div_green"+str(blue_div_green) + '\n')
		    
		    savefig("./figures/old_GW_costplots/blue_div_green/old_"+str(mode)+"_copper_costdiagram_blue_div_green"+str(ports)+"ports_alpha="+str(alph)+".pdf")		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_GW_costplots/blue_div_green/old_"+str(mode)+"_constr_costdiagram_blue_div_green"+str(ports)+"ports_alpha="+str(alph)+".pdf")
	close()
	
	for h in N:
	    for n in N:
		if  n.label==E[h.id]:
		   
		    usage_rel_to_share_of_load_kms[h.id]=t_kms[n.id]/1000 #yellow in GWKms
		    Blue_GWKms[h.id]=c_kms[n.id]
		    Green_GWkms[h.id]=usage_year_kms[n.id]/1000
		    total_usage_kms[h.id]=today_usage_kms[n.id]/1000
		    red_div_yellow_kms[h.id]=total_usage_kms[h.id]/usage_rel_to_share_of_load_kms[h.id]

	blue_div_green_kms=Blue_GWKms/Green_GWkms
	blue_div_green_kms_avg=round(blue_div_green_kms.sum()/(len(N)),2)
	
	p1=plt.bar(np.arange(30),blue_div_green_kms,width=0.7,color=au.oran)
	plt.axhline(blue_div_green_kms_avg,linestyle="dashed",color="k")
	plt.annotate(str(blue_div_green_kms_avg),(-0.5,blue_div_green_kms_avg),(1,blue_div_green_kms_avg+0.5),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"))
	plt.ylabel('I-factor')
	plt.xlim(j[0]-0.5)	
	plt.xticks(j+0.4,E,rotation=90)
	
	
	if new:
	    if copper:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    file.write("blue_div_green_GWKm-new"+str(E)+"blue_div_green"+str(blue_div_green_kms) + '\n')
		    np.savez("./figures/new_GWKMs_costplots/GW-numbers"+str(mode)+str(alph)+str(ports),labels=E, Blue=Blue_GWKms,Green=Green_GWkms,Yellow=usage_rel_to_share_of_load_kms)
		    savefig("./figures/new_GWKMs_costplots/blue_div_green/new_"+str(mode)+"_copper_costdiagram_kms_blue_div_green"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_GWKMs_costplots/blue_div_green/new_"+str(mode)+"_constr_costdiagram_kms_blue_div_green"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	else:
	    if copper:
		    file.write("blue_div_green_GWKm-old"+str(E)+"blue_div_green"+str(blue_div_green) + '\n')
		    np.savez("./figures/old_GWKMs_costplots/GW-numbers"+str(mode)+str(alph)+str(ports),labels=E, Blue=Blue_GWKms,Green=Green_GWkms,Yellow=usage_rel_to_share_of_load_kms)
		    savefig("./figures/old_GWKMs_costplots/blue_div_green/old_"+str(mode)+"_copper_costdiagram_kms_blue_div_green"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_GWKMs_costplots/blue_div_green/old_"+str(mode)+"_constr_costdiagram_kms_blue_div_green"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	close()
	
	p1=plt.bar(np.arange(30),red_div_yellow_kms,width=0.7,color=au.oran)
	plt.xlim(j[0]-0.5)	
	plt.xticks(j+0.4,E,rotation=90)
	
	if new:
	    if copper:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    file.write("red_div_yellow_GWKm-new"+str(E)+str(usage_rel_to_share_of_load_kms)+str(total_usage_kms)+"red_div_yellow"+str(red_div_yellow_kms) + '\n')
		    
		    savefig("./figures/new_GWKMs_costplots/red_yellow/new_"+str(mode)+"_copper_costdiagram_kms_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_GWKMs_costplots/red_yellow/new_"+str(mode)+"_constr_costdiagram_kms_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	else:
	    if copper:
		    file.write("red_div_yellow_GWKm-old"+str(E)+str(usage_rel_to_share_of_load)+str(total_usage)+str(total_usage_year)+"red_div_yellow"+str(red_div_yellow) + '\n')
		   
		    savefig("./figures/old_GWKMs_costplots/red_yellow/old_"+str(mode)+"_copper_costdiagram_kms_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_GWKMs_costplots/red_yellow/old_"+str(mode)+"_constr_costdiagram_kms_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
	
	close()   
	    
	for h in N:
	    for n in N:
		if n.mean==loadssort[h.id]:
		    E[h.id]=n.label
		    usage_rel_to_share_of_load_cost[h.id]=t_cost[n.id] #yellow in Euros
		    total_usage_cost[h.id]=today_usage_cost[n.id] #red in Euros
		    Blue_cost[h.id]=c_cost[n.id]
		    Green_cost[h.id]=usage_year_cost[n.id]
		    #total_usage_cost_year[h.id]=usage_year_cost[n.id]
		    red_div_yellow_cost[h.id]=total_usage_cost[h.id]/usage_rel_to_share_of_load_cost[h.id]
		
	
	p1=plt.bar(np.arange(30),red_div_yellow_cost,width=0.7,color=au.oran)
	plt.xlim(j[0]-0.5)	
	plt.xticks(j+0.4,E,rotation=90)
	
	if new:
	    if copper:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    file.write("red_div_yellow_Euro-new"+str(E)+str(usage_rel_to_share_of_load_cost)+str(total_usage_cost)+"red_div_yellow"+str(red_div_yellow_cost) + '\n')
		    savefig("./figures/new_TEuro_costplots/red_yellow/new_"+str(mode)+"copper_costdiagram_cost_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		    np.save("./figures/new_TEuro_costplots/New_costs"+str(mode)+str(ports)+"alpha="+str(alph),Blue_cost)
		    np.save("./figures/new_TEuro_costplots/New_costs_old_order"+str(mode)+str(ports)+"alpha="+str(alph),c_cost)
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_TEuro_costplots/red_yellow/new_"+str(mode)+"_constr_costdiagram_cost_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	else:
	    if copper:
		    np.save("./figures/old_TEuro_costplots/Old_costs"+str(mode)+str(ports)+"alpha="+str(alph),Blue_cost)
		    np.save("./figures/old_TEuro_costplots/Old_costs_old_order"+str(mode)+str(ports)+"alpha="+str(alph),c_cost)
		    file.write("red_div_yellow_GWKm-old"+str(E)+str(usage_rel_to_share_of_load)+str(total_usage)+str(total_usage_year)+"red_div_yellow"+str(red_div_yellow) + '\n')
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    #plt.text(0.03,0.03, "old linear", fontsize=15, color='w')
		    savefig("./figures/old_TEuro_costplots/red_yellow/old_"+str(mode)+"_copper_costdiagram_cost_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_TEuro_costplots/red_yellow/old_"+str(mode)+"_constr_costdiagram_cost_red_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
	close()
	for h in N:
	    for n in N:
		if  n.label==E[h.id]:		    
		    usage_rel_to_share_of_load_kms[h.id]=t_kms[n.id]/1000
		    #total_usage_kms[h.id]=today_usage_kms[n.id]/1000
		    total_usage_kms_year[h.id]=usage_year_kms[n.id]/1000
		    green_div_yellow_kms[h.id]=total_usage_kms_year[h.id]/usage_rel_to_share_of_load_kms[h.id]
		   
    	green_div_yellow_kms_avg=round(green_div_yellow_kms.sum()/(len(N)),2)
	
	p1=plt.bar(np.arange(30),green_div_yellow_kms,width=0.7,color=gyell)
	plt.axhline(green_div_yellow_kms_avg,linestyle="dashed",color="k")
	plt.annotate(str(green_div_yellow_kms_avg),(-0.5,green_div_yellow_kms_avg),(1,green_div_yellow_kms_avg+0.5),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"))
	plt.ylabel('I-factor')
	plt.xlim(j[0]-0.5)	
	plt.xticks(j+0.4,E,rotation=90)
	
	if new:
	    if copper:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    file.write("green_div_yellow_GWKm-new"+str(E)+str(usage_rel_to_share_of_load_kms)+str(total_usage_kms_year)+"green_div_yellow"+str(green_div_yellow_kms) + '\n')
		    savefig("./figures/new_GWKMs_costplots/green_yellow/new_"+str(mode)+"_copper_costdiagram_kms_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_GWKMs_costplots/green_yellow/new_"+str(mode)+"_constr_costdiagram_kms_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	else:
	    if copper:
		    file.write("green_div_yellow_GWKm-old"+str(E)+str(usage_rel_to_share_of_load_kms)+str(total_usage_kms_year)+"green_div_yellow"+str(green_div_yellow_kms) + '\n')
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    #plt.text(0.03,0.03, "old linear", fontsize=15, color='w')
		    savefig("./figures/old_GWKMs_costplots/green_yellow/old_"+str(mode)+"_copper_costdiagram_kms_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_GWKMs_costplots/green_yellow/old_"+str(mode)+"_constr_costdiagram_kms_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	close()
	for h in N:
	    for n in N:
		if  n.label==E[h.id]:		    
		    usage_rel_to_share_of_load[h.id]=t[n.id]/1000
		    total_usage_year[h.id]=usage_year[n.id]/1000
		    green_div_yellow[h.id]=total_usage_year[h.id]/usage_rel_to_share_of_load[h.id]	
		    
	green_div_yellow_avg=round(green_div_yellow.sum()/(len(N)),2)
	
	p1=plt.bar(np.arange(30),green_div_yellow,width=0.7,color=gyell)
	plt.axhline(green_div_yellow_avg,linestyle="dashed",color="k")
	plt.annotate(str(green_div_yellow_avg),(-0.5,green_div_yellow_avg),(1,green_div_yellow_avg+0.5),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"))
	plt.ylabel('I-factor')
	plt.xlim(j[0]-0.5)	
	plt.xticks(j+0.4,E,rotation=90)
	
	if new:
	    if copper:
		    file.write("green_div_yellow_GWKm-new"+str(E)+str(usage_rel_to_share_of_load)+str(total_usage_year)+"green_div_yellow"+str(green_div_yellow) + '\n')
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    
		    savefig("./figures/new_GW_costplots/green_yellow/new_"+str(mode)+"_copper_costdiagram_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_GW_costplots/green_yellow/new_"+str(mode)+"_constr_costdiagram_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	else:
	    if copper:
		    file.write("green_div_yellow_GWKm-old"+str(E)+str(usage_rel_to_share_of_load)+str(total_usage_year)+"green_div_yellow"+str(green_div_yellow) + '\n')
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    #plt.text(0.03,0.03, "old linear", fontsize=15, color='w')
		    savefig("./figures/old_GW_costplots/green_yellow/old_"+str(mode)+"_copper_costdiagram_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
    
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_GW_costplots/green_yellow/old_"+str(mode)+"_constr_costdiagram_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	close()
	
	    
	for h in N:
	    for n in N:
		if n.mean==loadssort[h.id]:
		    E[h.id]==n.label
		    usage_rel_to_share_of_load_cost[h.id]=t_cost[n.id]
		    #total_usage_cost[h.id]=today_usage_cost[n.id]
		    total_usage_cost_year[h.id]=usage_year_cost[n.id]
		    green_div_yellow_cost[h.id]=total_usage_cost_year[h.id]/usage_rel_to_share_of_load_cost[h.id]
    
	p1=plt.bar(np.arange(30),green_div_yellow,width=0.7,color=gyell)
	plt.xlim(j[0]-0.5)	
	plt.xticks(j+0.4,E,rotation=90)
	
	if new:
	    if copper:
		    file.write("green_div_yellow_cost-new"+str(E)+str(usage_rel_to_share_of_load_cost)+str(total_usage_cost_year)+"green_div_yellow"+str(green_div_yellow_cost) + '\n')
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    
		    savefig("./figures/new_TEuro_costplots/green_yellow/new_"+str(mode)+"linear_copper_costdiagram_cost_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
			
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_TEuro_costplots/green_yellow/new_"+str(mode)+"_constr_costdiagram_cost_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	else:
	    if copper:
		    file.write("green_div_yellow_cost-old"+str(E)+str(usage_rel_to_share_of_load_cost)+str(total_usage_cost_year)+"green_div_yellow"+str(green_div_yellow_cost) + '\n')
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    #plt.text(0.03,0.03, "old linear", fontsize=15, color='w')
		    savefig("./figures/old_TEuro_costplots/green_yellow/old_"+str(mode)+"_copper_costdiagram_cost_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_TEuro_costplots/green_yellow/old_"+str(mode)+"_constr_costdiagram_cost_green_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
	close()
	for h in N:
	    for n in N:
		if  n.label==E[h.id]:
		    blue_div_yellow[h.id]=c[n.id]/t[n.id]*1000
		    
	blue_div_yellow_avg=round(blue_div_yellow.sum()/(len(N)),2)	    		    
	   
	p1=plt.bar(np.arange(30),blue_div_yellow,width=0.7,color=au.kelp)
	plt.axhline(blue_div_yellow_avg,linestyle="dashed",color="k")
	plt.annotate(str(blue_div_yellow_avg),(-0.5,blue_div_yellow_avg),(1,blue_div_yellow_avg+0.5),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"))
	plt.ylabel('I-factor')
	plt.xlim(j[0]-0.5)	
	plt.xticks(j+0.4,E,rotation=90)
	
	if new:
	    if copper:
		    file.write("blue_div_yellow-new"+str(E)+str(blue_div_yellow) + '\n')
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    
		    savefig("./figures/new_GW_costplots/blue_div_yellow/new_"+str(mode)+"_copper_costdiagram_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_GW_costplots/blue_div_yellow/new_"+str(mode)+"_constr_costdiagram_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	else:
	    if copper:
		    file.write("blue_div_yellow-old"+str(E)+str(blue_div_yellow) + '\n')
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    #plt.text(0.03,0.03, "old linear", fontsize=15, color='w')
		    savefig("./figures/old_GW_costplots/blue_div_yellow/old_"+str(mode)+"_copper_costdiagram_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_GW_costplots/blue_div_yellow/old_"+str(mode)+"_constr_costdiagram_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
	close()
	
	for h in N:
	    for n in N:
		if  n.label==E[h.id]:
		    blue_div_yellow_cost[h.id]=c_cost[n.id]/t_cost[n.id]
		    
		    		    
	   
	p1=plt.bar(np.arange(30),blue_div_yellow_cost,width=0.7,color=au.kelp)
	
	plt.xlim(j[0]-0.5)	
	plt.xticks(j+0.4,E,rotation=90)
	
	if new:
	    if copper:
		    file.write("blue_div_yellow_Euro-new"+str(E)+str(blue_div_yellow_cost) + '\n')
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    
		    savefig("./figures/new_TEuro_costplots/blue_div_yellow/new_"+str(mode)+"_copper_costdiagram_cost_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_TEuro_costplots/blue_div_yellow/new_"+str(mode)+"_constr_costdiagram_cost_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	else:
	    if copper:
		    file.write("blue_div_yellow_Euro-old"+str(E)+str(blue_div_yellow_cost) + '\n')
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    #plt.text(0.03,0.03, "old linear", fontsize=15, color='w')
		    savefig("./figures/old_TEuro_costplots/blue_div_yellow/old_"+str(mode)+"_copper_costdiagram_cost_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_TEuro_costplots/blue_div_yellow/old_"+str(mode)+"_constr_costdiagram_cost_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	close()
	for h in N:
	    for n in N:
		if  n.label==E[h.id]:
		    blue_div_yellow_kms[h.id]=c_kms[n.id]/t_kms[n.id]*1000
		    
		    		    
	blue_div_yellow_kms_avg=round(blue_div_yellow_kms.sum()/(len(N)),2)
	   
	p1=plt.bar(np.arange(30),blue_div_yellow_kms,width=0.7,color=au.kelp)
	plt.axhline(blue_div_yellow_kms_avg,linestyle="dashed",color="k")
	plt.annotate(str(blue_div_yellow_kms_avg),(-0.5,blue_div_yellow_kms_avg),(1,blue_div_yellow_kms_avg+0.5),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"))
	plt.ylabel('I-factor')
	plt.xlim(j[0]-0.5)	
	plt.xticks(j+0.4,E,rotation=90)
	
	if new:
	    if copper:
		    np.save("new_capacities-year"+str(mode),Capacities_year)
		    file.write("blue_div_yellow_kms-new"+str(E)+str(blue_div_yellow_kms) + '\n')
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    
		    savefig("./figures/new_GWKMs_costplots/blue_div_yellow/new_"+str(mode)+"_copper_costdiagram_kms_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,new,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/new_GWKMs_costplots/blue_div_yellow/new_"+str(mode)+"_constr_costdiagram_kms_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	else:
	    if copper:
		    #np.save("old_capacities-year"+str(mode),capacities-year)
		    file.write("blue_div_yellow_kms-old"+str(E)+str(blue_div_yellow_kms)  + '\n' )
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,copper,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    #plt.text(0.03,0.03, "old linear", fontsize=15, color='w')
		    savefig("./figures/old_GWKMs_costplots/blue_div_yellow/old_"+str(mode)+"_copper_costdiagram_kms_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	    else:
		
		    #plt.text(0.03,0.03, "cost rel to meanload,old,linear,constr,"+str(ports)+"ports,alpha="+str(alph), fontsize=15, color='w')
		    savefig("./figures/old_GWKMs_costplots/blue_div_yellow/old_"+str(mode)+"_constr_costdiagram_kms_blue_div_yell"+str(ports)+"ports_alpha="+str(alph)+".pdf")
		
	close()
	
	


def track_new_flows(N,F,admat='./settings/eadmat.txt',lapse=None,mode=None): # tracks the flows (downstream and upstream powermixes) for the unprocessed case
	
	matr=np.genfromtxt(admat)

	if lapse==None:
		lapse=N[0].nhours

	a,b,list_F=au.AtoKh_old(N)
	start=time()
	total_power_mixes=[]
	total_power_mixes_ex=[]
	for n in N:
		n.links = get_links(n.id,matr)
		n.power_mix = np.zeros((len(N),lapse))
		n.total_power_mix=np.zeros((len(N)))
		n.total_power_mix_ex=np.zeros((len(N)))
		
		n.power_mix_ex=np.zeros((len(N),lapse))
		n.total_power_mix_ex=np.zeros((len(N)))
    
	for t in range(lapse):
		R=np.zeros(( len(N),len(N) ))
		R_ex=np.zeros(( len(N),len(N) ))
        
        ## Update progress bar.
		if mod(t,100)==0 and t>0: 
			print "\r",round(100.0*(t/float(lapse)),2),"%",
			sys.stdout.flush()
        
		for n in N:
            
            ## Add the nodes own contribution/taking source/sink strength.
			n.power_mix[n.id,t] = n.get_RES()[t]+n.balancing[t]-n.curtailment[t]
			n.power_mix_ex[n.id,t] = n.load[t]		            
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
				
	#If the flow indicated export from node n on link l			
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
				    
                    
            ### So now we have the own+direct input vector and the own+direct export vector
            ### we add it to the big matrix	
			
            ### we add it to the big matrix
			R[n.id]=n.power_mix[:,t]
			R_ex[n.id]=n.power_mix_ex[:,t]
			
       ############################################### inflow ##########################################
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
				
	
	##################################### outflow powermix ##################################    
		target_ex=np.arange(len(N))
		done_ex=[]
        
		while len(done_ex)<len(N): #All is not done yet.
			for n in target_ex:
				if n not in done_ex:  ### if it hasn't been done
					takers=[]
					for i in target_ex:
						if R_ex[n,i]>1. and n!=i: takers.append(i)
						#if n!=i: contributors.append(i)
					takers.sort()
                    ## Appears not to be used: contr=np.array(contributors)
                    
					if (np.in1d(takers,done_ex).all()) or (len(takers)==0): ### check if it is doable all it's takers are done or takers is empty
						for c in takers:                           ### then, for all "done" takers, do the thingy
							export_from_ntoc = R_ex[n,c]*1.0
							R_ex[n,c] = 0
							
							R_ex[n,:]   += R_ex[c,:]*export_from_ntoc/sum(R_ex[c,:])
							
						done_ex.append(n)
						done_ex.sort()
                        
		R_ex[np.isnan(R_ex)] = 0.0 #This should not be needed, but R says it is.
		 
		for n in N:
			n.power_mix_ex[:,t] = R_ex[n.id,:]
			n.power_mix[:,t] = R[n.id,:]   
			n.total_power_mix[:]+=n.power_mix[:,t]
			n.total_power_mix_ex[:]+=n.power_mix_ex[:,t]
			
	for n in N: 
				
		total_power_mixes.append(n.total_power_mix)
		total_power_mixes_ex.append(n.total_power_mix_ex)
	
	return N#, total_power_mixes, total_power_mixes_ex
	
	
	  

def track_new_exports(N,F,admat='./settings/eadmat.txt',lapse=None): # tracks the flows ( upstream powermixes) for the unprocessed case
	
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
