#! /usr/bin/env python
from time import time
import sys, os
import gurobipy as gb
from classes import *
from tools import *
import math
import random

def AtoKh(N):
    """ For better understanding, go open admat.txt file in the settings folder. 
        This function transforms an adjacency table A into the incidence matrix K 
        of size N x L with K[n,l] if link 'l' starts at node and -1 if it ends 
        there. Also returns 'h', which holds the actual transmission limits.

        In this version, K is returned as row and column indices and values, and 
        used to build the problem in gurobi. The final entry returns a list of 
        links with names, which is very useful when you get lost in the numbers."""

    Ad=np.genfromtxt(N.pathadmat,dtype='d')
    L=0
    listFlows=[]
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    L+=1
    K_values=[]
    K_column_indices=[]
    K_row_indices=[]
    T_caps=np.zeros(L*2)
    L=0
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    K_values.extend([1,-1])
                    K_column_indices.extend([L,L])
                    K_row_indices.extend([j,i])
                    T_caps[2*L]=Ad[i,j]
                    T_caps[2*L+1]=Ad[j,i]
                    listFlows.append([str(N[j].label)+" to " +str(N[i].label), L])
                    L+=1   
    return K_row_indices,K_column_indices,K_values,T_caps, listFlows

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################



def print_status(start=0,t=1,l=1000000,relaxed=0,init = False, msg="power flows"):
    a=time()-start
    b=l*a*1.0
    if init: 
        (a,b)=(0.0,0.0)
        print "Building problem for", msg,"\n\n\n\n"
    print "\033[4A\r\033[0mRelaxations performed:    \033[34m%2.0f" %(relaxed)
    print "\033[0mProgress:    \033[34m%2.2f" %(t*100.0/l)+"%"
    print "\033[0mElapsed time \033[34m%2.0i:%02.0i" % (a/60,np.mod(a,60))
    print "\033[0mETA \033[34m%2.0i:%02.0i " % ((b/t - a)/60,np.mod(b/t-a,60)),
    print "(%2.0i:%02.0i)          " % ((b/t)/60,np.mod(b/t,60))
    sys.stdout.flush()

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################


def build_network(N, copper = 0, h0 = None, b=1):

    K_row,K_col,K_val,H,LF=AtoKh(N)
    Nlinks=len(K_row)/2
    Nnodes=len(N)
    if (h0 != None):
        H=h0
    h_neg=b*-H[1:Nlinks*2:2]
    h_pos=b*H[0:Nlinks*2:2]

    network = gb.Model()    

    f_names=['f'+str(i+1) for i in range(Nlinks)]
    s_names=['s'+str(i+1) for i in range(Nnodes)]
    b_names=['b'+str(i+1) for i in range(Nnodes)]
    c_names=['c'+str(i+1) for i in range(Nnodes)]
    names = f_names + s_names + b_names + c_names

    #upper and lower bounds
    h_lower=h_neg
    h_upper=h_pos
    s_upper=np.ones(Nnodes)*1e9
    s_lower=np.ones(Nnodes)
    b_upper=np.ones(Nnodes)*1e9
    b_lower=np.zeros(Nnodes)
    c_upper=np.ones(Nnodes)*1e9
    c_lower=np.zeros(Nnodes)
    u_bound=np.concatenate((h_upper,s_upper,b_upper,c_upper))
    l_bound=np.concatenate((h_lower,s_lower,b_lower,c_lower))

    if not copper:
        for i in range(len(names)):
            network.addVar(lb=l_bound[i],ub=u_bound[i],name=names[i])
    if copper:
        for i in np.arange(0,Nlinks,1):
            network.addVar(lb = -gb.GRB.INFINITY, ub = gb.GRB.INFINITY, name=names[i])
        for i in np.arange(Nlinks,Nlinks+3*Nnodes,1):
            network.addVar(lb=l_bound[i],ub=u_bound[i],name=names[i])
   

    network.update()

    a_vars = network.getVars()
### For reference in the rest of the code
###############
#   f_vars = a_vars[:Nlinks]
#   s_vars = a_vars[Nlinks:Nlinks+Nnodes]
#   b_vars = a_vars[Nlinks+Nnodes:-Nnodes]
#   c_vars = a_vars[-Nnodes:]
########################################
    for n in N:  ## Restricts non-load nodes to not balance
        if np.average(n.load) <= 1e-10: 
            a_vars[Nlinks+Nnodes:-Nnodes][n.id].ub = 0.0 

    network.update()

################### This part essentially builds the X/Y variables, as 
################### greater than the negative/positive part of D-KF+S

    for n in range(Nnodes):
        ind=[]
        val=[]
        for i in range(len(K_row)):
            if K_row[i]==n:
                ind.append('f'+str(K_col[i]+1))
                val.append(K_val[i])
        ind.append('s'+str(n+1))
        val.append(-1.0)
        ind.append('b'+str(n+1))
        ind.append('c'+str(n+1))
        val.append(-1.0)   ## kf-s-b < D  -->  x > -d+kf-s
        val.append( 1.0)   ## kf-s+c > D  -->  y > d-kf+s
        var=[]
        for i in ind:
            var.append(network.getVarByName(i))
        network.addConstr(lhs=gb.LinExpr(val,var),sense='=',rhs=1e6)

######################### Setting default values for the solver
###############################################################

    network.setParam("OutputFlag",0)
    network.setParam("FeasibilityTol",1e-4)          
    network.update()
    return network, names, a_vars

###############################################################
###############################################################
###############################################################
###############################################################


def solve(N, copper=0, h0=None, b=1.0, lapse=None, mode = "linear", msg="power flows", verbose = 0):
    """ This main function builds the network problem, gathers variables and runs 
        the time series. It compiles the results and returns them as an Nodes 
        object and a Nlinks x Nhours flow matrix.

        copper = 1 signals an unconstrained entwork. In reality, just one with very 
        high interconnector capacities (1.0 x 10^6). Default is 0.

        h0 can receive a customized interconnector capacity vector of size 
        (2*Nlinks). Else it uses the capacities saved in the admat matrix, part 
        of the Nodes object. Default is None.

        b is a linear multiplication factor for whatever interconnector capacities,
        either the default ones or the ones defined by h0. Default is 1.

        lapse defines the timespan to run the simulation, starting at 0. Default is 
        None and reverts to Nhours.

        mode = [linear, square, hybrid, random]
    """

###############################################################
################ Loading Default values #######################
###############################################################
    if not verbose: print str(msg)
    if lapse == None:
        lapse = np.arange(0,N[0].nhours,1)
    if type(lapse) == int:
        lapse = np.arange(0,lapse,1)

    K_row,K_col,K_val,H,LF=AtoKh(N)
    Nlinks=len(K_row)/2
    Nnodes=len(N)
    l = N[0].nhours
    M = np.zeros((Nnodes,Nlinks))
    for i in range(len(K_row)):
        M[K_row[i]][K_col[i]]=K_val[i] # M is adjacency matrix 'K'

    ADJ = (M, K_row, K_col, K_val)

    F = np.zeros((Nlinks,l)) #Flows
    S = np.zeros((Nnodes,l)) #Storage
    B = np.zeros((Nnodes,l)) #Balancing
    C = np.zeros((Nnodes,l)) #Curtailment

###############################################################
################# Setting up the model ########################
###############################################################

    network, names, a_vars = build_network(N, copper, h0, b)
    mean_n=[node.mean for node in N] 
###############################################################
##### Initial balancing constraints for random's case #########
###############################################################
    if mode =="random":
		
	for n in N: 
	    n.mismatchete=np.zeros(len(lapse))
	   
	    n.temp_mismatch=np.zeros(len(lapse))
	    
	    
	for t in lapse:
	    
	    sinks=[]
	    sources=[]
	    sur_or_def=np.zeros(len(lapse))	    
	    sur_or_def[t]=sum(n.mismatch[t] for n in N) # tells us whether the global european network has a total surplus or deficit of RES in the given hour.
	    for n in N:	    	
		if n.mismatch[t]>0:
		    sources.append(n)
		elif n.mismatch[t]<0:
		    sinks.append(n)
		    
		si=len(sinks)
		so=len(sources)			    
		n.temp_mismatch[t]=n.mismatch[t]
	    intervalssi=np.linspace(0,1,si+1,endpoint=True)
	    intervalsso=np.linspace(0,1,so+1,endpoint=True)
	    if sur_or_def[t]<=0: #global deficit
			
		while not all(n.temp_mismatch[t]==0 for n in sources):	    
		    rsi=random.random()
		    rso=random.random()
		    for i in range(si):
			if rsi >=intervalssi[i] and rsi<intervalssi[i+1]:												
			    if (sinks[i].temp_mismatch[t] <> 0):
				for o in range(so):
				    if rso >=intervalsso[o] and rso<intervalsso[o+1]:
					 if (sources[o].temp_mismatch[t] <> 0):	
					    trade=min(sources[o].temp_mismatch[t],-sinks[i].temp_mismatch[t],1000)					    
					    sources[o].temp_mismatch[t]-=trade
					    sinks[i].temp_mismatch[t]+=trade					   
					        
				  
		    
	    
	    elif sur_or_def[t]>0: #global surplus		
		while not all(n.temp_mismatch[t] == 0 for n in sinks):	    
		    rsi=random.random()
		    rso=random.random()
		    for i in range(si):			
			    if rsi >=intervalssi[i] and rsi<intervalssi[i+1]:	
				if (sinks[i].temp_mismatch[t] <> 0):										    
				    for o in range(so):						
					if rso >=intervalsso[o] and rso<intervalsso[o+1]:
					    if (sources[o].temp_mismatch[t] <> 0):
						trade=min(sources[o].temp_mismatch[t],-sinks[i].temp_mismatch[t],1000)					    
						sources[o].temp_mismatch[t]-=trade
						sinks[i].temp_mismatch[t]+=trade
		
	    
	    for n in N:		
		n.mismatchete[t]=n.mismatch[t]-n.temp_mismatch[t] #the mismatch that will actually be dealt with
		n.imports[t]=-min(n.mismatchete[t],0) #the part of the negative mismatches that will be covered by trading
		n.exports[t]=max(n.mismatchete[t],0) #the part of the positive mismatches that will be exported by trading.
		       
	
	    
	
				
		
			
				

###############################################################
###### Initial balancing constraints for capped's case #########
###############################################################

    if mode == "capped":
	weights = mean_n/sum(mean_n)
        Global_B = np.max( -sum(n.mismatch for n in N) )
        weights *= Global_B*0.25
        for x in range(Nnodes):
            a_vars[Nlinks+Nnodes:-Nnodes][x].ub = float(weights[x])
    network.update()
    


###############################################################
################# Run the time series  ########################
###############################################################

    start=time()
    relaxed = 0
    if verbose: print_status(init = True)
    for t in lapse:

        if mode not in ["hybrid", "capped"]:
            solution, r = _solve_flows_(network, N, t, ADJ, mode, mean_n)
            relaxed += r

        elif mode == "hybrid":
            solution, r = _hybrid_solve_(network, N, t, ADJ, mean_n)
            relaxed += r

        elif mode == "capped":
            solution, r = _capped_solve_(network, N, t, ADJ, mean_n)
            relaxed +=r
        ###### save hourly solutions in long-term matrices
        F[:,t]=solution[0:Nlinks]
        S[:,t]=solution[Nlinks:Nlinks+Nnodes]
        B[:,t]=solution[Nlinks+Nnodes:Nlinks+Nnodes*2]
        C[:,t]=solution[-Nnodes:]
	
	if mode == "random":
	    F[:,t]=solution[0:Nlinks]
	    S[:,t]=solution[Nlinks:Nlinks+Nnodes]
	    B[:,t]=solution[Nlinks+Nnodes:Nlinks+Nnodes*2]#+bal_mismatch[t] 
	    C[:,t]=solution[-Nnodes:]#+curtail_mismatch[t] 
	    
	    
        for i in N:
            i.ESS.charge(-S[i.id,t]) # The sign of S is opposite to the sign 
                                     # convention used inside the storage class.
            i.balancing[t] = B[i.id,t]
            i.curtailment[t] = C[i.id,t]
            i.storage_discharge[t] = S[i.id,t]
        if t>0 and verbose:
            print_status(start=start,t=t,l=l,relaxed=relaxed)
    if verbose: print_status(start=start,t=t,l=l,relaxed=relaxed)

    for i in N:
        i._update_()
        i.solved = True
    T=time()-start
    
    print "Calculation of "+msg+" took %2.0f:%02.0f" % (T/60.0-np.mod(t/60.0,1),np.mod(T,60))
    #if mode=="random":
	#return dcopy(N),dcopy(F)
    #else:
    return dcopy(N),dcopy(F)


###############################################################
###############################################################
###############################################################
###############################################################


def _solve_flows_(network,N,t, ADJ, mode, mean_n):
    """ Non-interactive function, receives a network, a timestep, a bunch of other 
        stuff :D and solves both steps of the optimisation process.
    """

    (M , K_row, K_col, K_val) = ADJ
    funny_nodes=0
    relaxed = 0
###### Build variable names ######
##################################
    Nlinks=len(K_row)/2
    Nnodes=len(N)

    a_vars = network.getVars()
    f_vars = a_vars[:Nlinks]
    s_vars = a_vars[Nlinks:Nlinks+Nnodes]
    b_vars = a_vars[Nlinks+Nnodes:-Nnodes]
    c_vars = a_vars[-Nnodes:]

#### Get mismatch & ESS power ####
##################################

    Delta=[i.mismatch[t] for i in N]
    if mode == "random":
	Delta=[n.mismatchete[t] for n in N]
	#Delta=[n.bongo[t] for n in N]
    s_upper=np.array(N.get_P_discharge(t)) ## NOTE: This may be very wrong!
    s_lower=np.array(-N.get_P_charge(t)) ## NOTE: This may be very wrong!
    for n in range(Nnodes):
        s_vars[n].ub=float(s_upper[n])
        s_vars[n].lb=float(s_lower[n])


#### Set deltas for X-Y Constraints ####
########################################
    r=0
    for bc_constraint in network.getConstrs():
        bc_constraint.setAttr("rhs",float(Delta[r]))
        r+=1

    network.update()
    
    #if mode == "random":
	#p=0
	#for bc_constraint in network.getConstrs():
	    #bc_constraint.setAttr("rhs",float(Delta_temp_random[p]))
	    #p+=1
    network.update()
    
#### Set objective factors based on mode ####
#############################################

    if mode in ["linear","capped","random"]:
        b1 = np.ones(Nnodes)
        c1 = np.ones(Nnodes)
        b2 = np.zeros(Nnodes)
        c2 = np.zeros(Nnodes)

    if mode in ["square"]: ## Shares curtailment! 
        b1 = np.zeros(Nnodes)
        c1 = np.zeros(Nnodes)
        b2 = np.ones(Nnodes)/mean_n
        c2 = np.zeros(Nnodes)/mean_n
    
    linear_factors = np.concatenate((b1 , c1))
    square_factors = np.concatenate((b2 , c2))
    
    


##### Define Step-1 objective #####
###################################
    step1_objective=gb.QuadExpr()
    step1_objective.addTerms(linear_factors, b_vars+c_vars)
    step1_objective.addTerms(square_factors, b_vars+c_vars, b_vars+c_vars)
    network.setObjective(expr=step1_objective, sense=1)


##### Solve Step-1 objective #########
######################################
    network.update()
    network.optimize()
    v=[]
    
    for i in network.getVars():
        v.append(i.x)
    BC_opt=network.objVal
    
##### Add constraints for Step-2 #####
######################################
    if mode in ["square"]:
        step2_cs = 0
        b = v[Nlinks+Nnodes:-Nnodes]
        #c = v[-Nnodes:]
        for n in range(Nnodes):
            b_vars[n].ub = b[n]*1.00001
            b_vars[n].lb = b[n]*0.99999
            #c_vars[n].ub = c[n]*1.00001
            #c_vars[n].lb = c[n]*0.99999
    else:
        step2_cs = 1
        step2_constraint=gb.LinExpr()
        step2_constraint.addTerms(linear_factors, b_vars+c_vars)
        network.addConstr(lhs=step2_constraint, sense="<", rhs=BC_opt, name="Step 2 constraint")
    network.update()

###### Set Step-2 Objective ######
##################################
    step2_objective=gb.QuadExpr()
    step2_objective.addTerms(np.ones(Nlinks), f_vars, f_vars)
    network.setObjective(expr=step2_objective,sense=1)


###### Solve Step 2 ######
##########################

    network.update()
    v=[]
    try:
        v=[]
        network.optimize()
        for i in network.getVars():
            v.append(i.x)
    except gb.GurobiError:
        #~ print "second step error"
        relaxed += 1
        v=[]
        networkRelaxed=network.copy()
        networkRelaxed.feasRelaxS(0,0,0,1)
        networkRelaxed.optimize()
        for i in networkRelaxed.getVars():
            v.append(i.x)
    
	
    

###### Clean-up ######
##########################
    for c in range(step2_cs):
        network.remove(network.getConstrs()[-1])
        network.update()
    if mode in ["square"]:
        for n in range(Nnodes):
            b_vars[n].ub = gb.GRB.INFINITY
            b_vars[n].lb = 0.0
            c_vars[n].ub = gb.GRB.INFINITY
            c_vars[n].lb = 0.0
    return v, relaxed


def _hybrid_solve_(network, N, t, ADJ, mean_n):
    mode = "linear"
    relaxed = 0
    solution, r = _solve_flows_(network, N, t, ADJ, mode, mean_n)
    relaxed+= r
    curt_const = solution[-len(N):]  #This is the Y value of solution
    for n in range(len(N)):
        network.getVars()[-len(N):][n].ub = curt_const[n]
    network.update()
    mode = "square"
    solution, r = _solve_flows_(network, N, t, ADJ, mode, mean_n)
    relaxed +=r
    mode = "hybrid"
    for n in range(len(N)):
        network.getVars()[-len(N):][n].ub = gb.GRB.INFINITY
    network.update()
    return solution, relaxed

def _capped_solve_(network, N, t, ADJ, mean_n):

    (M , K_row, K_col, K_val) = ADJ
    Nlinks=len(K_row)/2
    Nnodes=len(N)
    b_vars = network.getVars()[Nlinks+Nnodes:-Nnodes]
    relaxed = 0
    converged = 0
    nonconv = []
    preweight=[]
    
    network.update()
    while not converged:
        try:
            
            mode = "capped"
            solution, r =_solve_flows_(network, N, t, ADJ, mode, mean_n)
            converged = 1

        except gb.GurobiError:
            network.computeIIS()
            for x in range(Nnodes):
                if b_vars[x].getAttr("IISUB"):
                    nonconv.append(x)
                    preweight.append(b_vars[x].ub)
                    b_vars[x].ub += 1.0
                    b_vars[x].ub *= 1.001
            for const in network.getConstrs()[Nnodes:]:
                network.remove(const)

            network.update()
            relaxed += 1

    for n in range(len(nonconv)):
        b_vars[nonconv[n]].ub=preweight[n]
        a = dcopy(float(solution[Nlinks+Nnodes:-Nnodes][nonconv[n]]))
        if a>preweight[n]:
            b_vars[nonconv[n]].ub=a
    network.update()
    if relaxed>=1e9:
        print "After time = ",t
        raw_input()
        for n in nonconv:
            print N[n].id, N[n].label, N[n].mismatch[t], solution[Nlinks+Nnodes:-Nnodes][n], b_vars[n].ub

    
    return solution, relaxed
    
	
