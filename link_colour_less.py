#! /usr/bin/env/python
import sys
from pylab import *
import datetime as dt
from time import time
from EUgrid import *
from aurespf.tools import AtoKh_old

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
	
		
def track_node_links_usage(links,N=None,ports="ex_and_im",new=False,mode="linear",alph=None,copper = True):
# as links use links_mix or variants of this, set new=True if you use the new model
#gives the linkusage for each country relative to its total usage of links. The first 50 entries in 
#Mode can be "linear", "square", "random" or "capped".
	if alph== None:
	    alph="hetero"
	if N==None:
		N=EU_Nodes()
	x=len(N)
	
	a,b,c,d,e=au.AtoKh(N)
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
			#k=linksflow[j]
			#cost=np.append(cost,links[n.id+x*j]/k)
		n.link_usages_normalized=np.append(n.link_usages_normalized,n.link_usages[:]/total_usage*100) # calculates node n's usage of the links compared to node n's total usage of the links		
		total_link_usage_node_rel=np.append(total_link_usage_node_rel,n.link_usages_normalized) # stores them in an array with the 50 first entries storing the relative linkusages for node 0 and the 50 next for node 1 and so on.
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
				np.save('linkcolouring/old_random_copper_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel)
				np.save('linkcolouring/old_random_copper_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
		elif mode == "capped":
			if new:
				np.save('linkcolouring/new_capped_copper_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel) 
				np.save('linkcolouring/new_capped_copper_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
				#np.save('linkcolouring/new_cost_links_alpha='+str(alph),cost)
			else:
				np.save('linkcolouring/old_capped_copper_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel)
				np.save('linkcolouring/old_capped_copper_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
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
				np.save('linkcolouring/old_random_constr_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel)
				np.save('linkcolouring/old_random_constr_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
		elif mode == "capped":
			if new:
				np.save('linkcolouring/new_capped_constr_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel) 
				np.save('linkcolouring/new_capped_constr_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
				#np.save('linkcolouring/new_cost_links_alpha='+str(alph),cost)
			else:
				np.save('linkcolouring/old_capped_constr_node_links_rel_'+str(ports)+'_alpha='+str(alph),total_link_usage_node_rel)
				np.save('linkcolouring/old_capped_constr_node_links_'+str(ports)+'_alpha='+str(alph),total_link_usage_node)
			
			
	
	
	


def track_flows(N,F,admat='./settings/eadmat.txt',mode=None,lapse=None): #tracks the flows (downstream and upstream powermixes) for the processed case.
	
	matr=np.genfromtxt(admat)

	if lapse==None:
		lapse=N[0].nhours

	a,b,list_F=au.AtoKh_old(N)
	start=time()
	total_power_mixes=[]
	total_power_mix_exes=[]
	for n in N:
		n.links = get_links(n.id,matr)
		n.power_mix = np.zeros((len(N),lapse))
		n.total_power_mix=np.zeros((len(N)))
		n.total_power_mix_ex=np.zeros((len(N)))
		n.power_mix_ex=np.zeros((len(N),lapse))
    
	for t in range(lapse):
		R=np.zeros(( len(N),len(N) ))
        
        ## Update progress bar.
		if mod(t,100)==0 and t>0: 
			print "\r",round(100.0*(t/float(lapse)),2),"%",
			sys.stdout.flush()
        
		for n in N:
            
            ## Add the nodes own contribution/source strength.
		    if mode =="random":
				n.power_mix[n.id,t] = n.exports[t]
				n.power_mix_ex[n.id,t]=n.imports[t]
                #print n.power_mix[n.id,t]
			
		    else:
				n.power_mix[n.id,t]=n.get_export()[t]
				n.power_mix_ex[n.id,t]=n.get_import()[t]
            
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
                    
					if (np.in1d(contributors,done).all()) or (len(contributors)==0): ### check if it is doable all its contributors are done or contributors is empty
						for c in contributors:                           ### then, for all "done" contributors, do the thingy
							supply_from_cton = R[n,c]*1.0
							
							R[n,c] = 0
							
							R[n,:]   += R[c,:]*supply_from_cton/sum(R[c,:])
							
						done.append(n)
						done.sort()
                        
		R[np.isnan(R)] = 0.0 #This should not be needed, but R says it is.
		 
		for n in N:
			n.power_mix[:,t] = R[n.id,:]   
			n.total_power_mix[:]+=n.power_mix[:,t]	
		total_power_mixes.append(n.total_power_mix)
			
############################# outflows ##################################			
		
	
		    
	for t in range(lapse):
		R=np.zeros(( len(N),len(N) ))
        
        ## Update progress bar.
		if mod(t,100)==0 and t>0: 
			print "\r",round(100.0*(t/float(lapse)),2),"%",
			sys.stdout.flush()
        
		for n in N:
            
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
			n.total_power_mix_ex[:]+=n.power_mix_ex[:,t]
	for n in N:
	    total_power_mix_exes.append(n.total_power_mix_ex)
	
			          

	return N, total_power_mixes
def track_exports(N,F,admat='./settings/eadmat.txt',mode=None,lapse=None):
	
	matr=np.genfromtxt(admat)

	if lapse==None:
		lapse=N[0].nhours

	a,b,list_F=au.AtoKh_old(N)
	start=time()
	total_power_mix_exes=[]
	for n in N:
		n.links = get_links(n.id,matr)
		n.power_mix_ex = np.zeros((len(N),lapse))
		n.total_power_mix_ex = np.zeros((len(N)))
		total_power_mix_exes=[]
    
	for t in range(lapse):
		R=np.zeros(( len(N),len(N) ))
        
        ## Update progress bar.
		if mod(t,100)==0 and t>0: 
			print "\r",round(100.0*(t/float(lapse)),2),"%",
			sys.stdout.flush()
        
		for n in N:
		    if mode =="random":			
			
			n.power_mix_ex[n.id,t]=n.imports[t]
			#print n.power_mix[n.id,t]
			
		    else:
			
			n.power_mix_ex[n.id,t]=n.get_import()[t]
            ## Add the nodes own export/sink strength.
	    
            
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
			n.total_power_mix_ex[:]+=n.power_mix_ex[:,t]
			
	for n in N:
	    total_power_mix_exes.append(n.total_power_mix_ex)
			
			          

	return N, total_power_mix_exes

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
