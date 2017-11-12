# -*- coding: utf-8 -*-
########################################################
####### Genetic Scheduler 490; GS490 
####### 
####### John Hokkanen (johnhokkanen2014@u.northwestern.edu)
####### 
####### Provided on the basis that any bugs are reported to the author
####### and that you share revisions or enhancements with the class.
####### If you don't want to share your insights, then don't use the code.
#
# You should change the following four lines to make them appropriate to your system:
# Line 30 os.chdir('D:/Content/John/Northwestern/Classes/Algorithms/CaseStudy')
# Line 31 execfile('D:/Content/John/Northwestern/Classes/Algorithms/Code/fitness_function2.py')
# Line 54 pathtoOrderData = 'D:/Content/John/Northwestern/Classes/Algorithms/CaseStudy/OrderData.csv'
# Line 55 pathtoDistanceData = 'D:/Content/John/Northwestern/Classes/Algorithms/CaseStudy/FromToDist.csv'
# Set the gen_count to a small # or interrupt when desired.  Upon interruption, then run line this 
# code to get the output: fitness_function(routes[TopItem,],mydata,mydistances,zero,1)
# See the sample xls file to look at how to post-process the output.
#
# You probably will want to change the truck count or let it run longer to see how that changes results.
#####################################################
############## PSEUDO CODE ##########################
# Read in all data files and initialize system
# Create the initial population of random routes
# Reallocate routes to give high-volume (i.e., >2400 ft3) packages an entire truck (the system will then optimize the reallocation of packages in other trucks to fill the truck).
# Create the breeding population if proportional breeding enabled
# Loop for required generations or until fitness criteria achieved
# 	Create a child population of the same size as the current population:
# 		Always clone the most fit member of the parent population and make that the last child in the child population.
# 		If proportional breeding, randomly select parents based on proportional fitness.
# 		If parallel individual breeding, then each parent acts as the base for one child.
# 	Apply the mutation rules to each child in the population
# 		Create a unique mutation profile from the various mutation functions in  Appendix A for the individual
# 		Apply the mutations
# 	Calculate the fitness of the individual
# 	Determine if the parent plan should be replaced with the new individual’s plan
# 		If the child plan is more fit than the parent plan:
# 			Substitute the child plan as the new plan
# 			Reset the spawn counter
# 			If the new plan is the top plan, save it as the top item for printing purposes
# 		If search expansion is enabled:
# 			If the child plan = parent plan:
# 				Substitute the child plan as the new plan
# 				Reset the spawn counter
# 			Else if the spawn counter > specified threshold:
# 				If the difference between the two plans is less than the linear amount specified by the spawn counter:
# 					Substitute the child plan as the new plan
# 					Reset the spawn counter
# 	Calculate the proportional fitness of each individual in the population and create the breeding population if proportional fitness is used
# 	Output any diagnostic information about the generation to the console
# Output the top item best plan information
########################################################
import sys
import os
import numpy as np
import pandas as pd
import random
import pickle
from copy import deepcopy
os.getcwd()
from past.builtins import execfile
os.chdir('C:/CaseStudy')
execfile('C:/CaseStudy/fitness_function2.py')
######################## INITIALIZE CRITICAL VARIABLES ###########################
random.seed(0)            #seed value for recurrent testing
weekday = 1               #SET WEEKDAY - 1 to 5 (Last: 4)
trucks = 5                #SET HOW MANY TRUCKS TO USE (Last: 7)
pop_size = 5              #SET HOW MANY SIMULTANEOUS SOLUTIONS TO PURSUE (Minimum: 2)
gen_count=2500            #SET HOW MANY GENERATIONS (RECOMMENDED: 100000, and manually break at desired fitness level)
proportionalbreederenabled=0   #Proportional/Parallel breeding; 1=proportional and 0=parallel solutions 
                               #What does this mean?  If proportional, subsequent generations are driven by
                               #proportional fitness.  If parallel, then each individual breeds and, if the 
                               #progeny is better, then it is replaced and so all n solutions are pursued 
                               #independently.  Must have a large sample size to use this, i.e., >=50.
disableallprofiles = 1         #if =1, then specific profiles will be disabled and default profiles are used 
MutationComplexity=5           ##Not yet implemented 1,2,3...to 9; 5 is uniform distribution, 10="V" distribution.
                               #   At 1, simpler mutations are emphasized. 9 emphasises complex transformations
ForceSpecificPressures=1       #Not yet implemented This will emphasize separation of edges and individual gene mutations
enablesearchexpansion = 1      #Enable progeny mutation heuristic.  Recommended.
enablesearchupward =1          #Enable inferior progeny heuristic.  Recommended.
searchexpthold=25              #Threshold for inferior progeny. 25 generations recommended.
searchexppenalty=.1            #Penalty*attemptsatprogeny, ex: 100 failures*.1 = accept a 10 point inferior child
display_searchexpansion=1      #turn on/off display
display_MutationInfo=1         #Display the Mutation sequence
display_EveryRecord=0          #If every generation record displayed. Useful with MutationInfo
#pathtoOrderData = 'D:/Content/John/Northwestern/Classes/Algorithms/CaseStudy/OrderDataExperimental.csv'
pathtoOrderData = 'C:/CaseStudy/OrderData.csv'
pathtoDistanceData = 'C:/CaseStudy/FromToDist.csv'
firstrunever=0                 #set this to 1 for your first run so plan files are created.
loadbestplanever=-2            #{-2,-1,0,1} If -2: disable, -1:start random, 0: load into plan1, 1:load into all plans.
#
#########################################################
#### CODE STARTS HERE
#########################################################
#pathtoDistanceData = 'D:/Content/John/Northwestern/Classes/Algorithms/CaseStudy/DistanceToFrom.csv' old file
OrderData=pd.read_csv(pathtoOrderData, sep=',')
Distances2=pd.read_csv(pathtoDistanceData, sep=',')
#\\\\\\\\\\\\\\\\ EXPERIMENTAL //////////////
pop_size_scaling=0
pop_size_finish = 30
pop_size_transition = 30  #generations when switch
#^^^^^^^^^^^^^^^^ EXPERIMENTAL ^^^^^^^^^^^^^^
zero = 0
dispctr=-1
PrevPlan =-1
TopPlanctr = 0
PreviousTopItem=0
PreviousTopPlan=-1
displayit=0
dispstring = ""
Mutationlist = list()
mutationfunctions = 10
edgefunctionstart = 4
searchexpdisplay=""

mydata = OrderData[OrderData.dayofweek==weekday]
mydistances = Distances2[Distances2.dayofweek==weekday]
orders = len(mydata)
mydataindex = list(range(0,orders))
mydata['index']=mydataindex
mydata.set_index(['index'], inplace = True)
mydata['genenum'] = mydata.index

multiples = mydata[mydata.multipleorder==1]
multiplescount = len(multiples)
triples= mydata[mydata.tripleorder==1]
triplescount = len(triples)
if triplescount==0:
    mut_consolidatetriple_prob = 0    #don't use if there are no triples.

#create initial generation
plans = np.zeros(pop_size,dtype=int)    #this holds the fitness value for each plan
prevplanlist = np.zeros(pop_size,dtype=int)    #this holds the fitness value for each plan, previously
planscounts = np.zeros(pop_size,dtype=int)    #this holds the counts for static plans
attemptsatoffspring = np.zeros(pop_size,dtype=int)    #this holds the fitness value for each plan
Breeder = np.zeros(pop_size,dtype=int)    #this holds the index of the Breeders
planfitnessproportion = np.zeros(pop_size,dtype=int)    #this holds the fitness value for each plan
keepplan = np.zeros(pop_size,dtype=int)  #to keep a particular plan in stock
routes=np.empty((pop_size,trucks,),dtype=object)  #makes array of lists for chromosomes for trucks
routefitness = np.zeros((pop_size,trucks),dtype=int)
#same as above but to store temporary offspring
newplans = np.zeros(pop_size,dtype=int)
newroutes=np.empty((pop_size,trucks,),dtype=object)
newroutefitness = np.zeros((pop_size,trucks),dtype=int)
MutateDist = np.zeros(10,dtype=int)   #Mutation profiles
bestroute=np.empty((pop_size,trucks,),dtype=object)  #makes array of lists for chromosomes for trucks
bestplanvalue = np.zeros((pop_size),dtype=int)
bestplanvalue[:]=1000000
bestplangeneration = np.zeros((pop_size),dtype=int)
bestrouteever=np.empty((5,10,),dtype=object)
bestplanever = np.zeros((5,10),dtype=int)  #5 days, 10 trucks
bestplanever[:] = 1000000
#
if firstrunever==1:   #if the first time running the program, create the bestroute files.
    afile = open(r'GS490bestplan.pkl', 'wb')
    pickle.dump(bestplanever, afile)
    afile.close()
    bfile = open(r'GS490bestroute.pkl', 'wb')
    pickle.dump(bestrouteever, bfile)
    bfile.close()

if loadbestplanever!=-2:    
    afile = open(r'GS490bestplan.pkl', 'rb')   #read the bestroutefiles
    bestplanever = pickle.load(afile)
    afile.close()
    bfile = open(r'GS490bestroute.pkl', 'rb')
    bestrouteever = pickle.load(bfile)
    bfile.close()

###################################################
#Begin pre-processing setup.
#Check for high-volume orders that exceed volume
listoftrucks =list(range(0,trucks))  #this will be list of trucks to store orders
highvoltrucks = mydata[mydata.cube>2400]
totalvolume = sum(mydata['cube'])
if totalvolume/trucks > 3200:
    print("DOT Violation: You must increase the number of trucks.")
    sys.exit()

#initialize the population
for i in range(0, pop_size):
    stringstart = 0 
    randomplan=random.sample(list(range(0,orders)),orders) 
    ordersize = int(orders/trucks)  # this could be random, but might generate defectives
    for j in range (0, trucks):
        stringend = ordersize+stringstart
        if j==trucks-1:
            stringend=orders
        routes[i,j] = randomplan[stringstart:stringend]
        stringstart = stringend

#Check for high volume orders and repack accordingly
    notthistruck = deepcopy(listoftrucks)
    for highvolorder in highvoltrucks['orderid']:
        #Find the truck with the high volume order and remove all the packages other than that one
        for rte1 in range(0,trucks):
            genenum=int(mydata['genenum'][mydata.orderid==highvolorder])
            indexval = routes[i,rte1].index(genenum) if genenum in routes[i,rte1] else -1
            if indexval > -1:    #found it
                break
        #Remove this truck, then move all remaining orders, then put it back.
        routes[i,rte1].remove(genenum)  #remove highvol order, we'll put it back later
        notthistruck.remove(rte1)          #remove this truck from the available list
        rte1list = deepcopy(routes[i,rte1])
        for movegene in rte1list:
            routes[i,rte1].remove(movegene)
            rte2 = notthistruck[random.randint(0,len(notthistruck)-1)]
            routes[i,rte2].insert(0,movegene)
        routes[i,rte1].insert(0,genenum)

#
# IF INSTRUCTED TO USE PAST PLANS
if loadbestplanever==0:                            #load only once into plan 0 
    routes[0,] = bestrouteever[weekday, trucks]
elif loadbestplanever==1:                          #load into every plan 
    for i in range(0, pop_size):
        routes[i,] = bestrouteever[weekday, trucks]

#Calculate proportion of total value
for i in range(0, pop_size):
    plans[i]=fitness_function(routes[i,],mydata,mydistances,zero,0)
totalfitnessinpop = plans.sum(axis=0)   #sum total fitness in population
planfitnessproportion = totalfitnessinpop/plans  #this is inverse since we are minimizing
prevplanlist = deepcopy(plans)
#####
generation_ctr=-1   #start the generation counter
##########################################################################
##########################################################################
##########################                ################################
##########################      START     ################################
##########################                ################################
##########################################################################
##########################################################################
#### If you interrupt the processing Select from here down and RUN SELECTION       

while generation_ctr<gen_count:
    generation_ctr=generation_ctr+1
    ProbabilityList = list()    #Generate proportional probability list
    for i in range(0, pop_size):
        sublist = [i]*int(planfitnessproportion[i]) #repeat this index in proportion to its overall weight
        ProbabilityList.extend(sublist)

    if pop_size_scaling==1 and generation_ctr==pop_size_transition:
        pop_size = pop_size_finish
        TopItem = 0
            
    for i in range(0, pop_size):     #Select breeders randomly based from their probability distrib created
        Breeder[i] = int(random.sample(ProbabilityList,1)[0])
        #last item will not be used as it is for the TopItem
           
   # Breeder.sort()   #Order the breeder list for convenience for display
    
    for i in range(0,pop_size):
        if proportionalbreederenabled==0:
            newroutes[i] = deepcopy(routes[i,])
        else:
            if i<pop_size-1:            
                newroutes[i] = deepcopy(routes[Breeder[i,]])    
            else: #last one has extra copy of TopPlan so that it is always preserved
                newroutes[i] = deepcopy(routes[i,])    
#######################################################################################
############################ SET UP THE MUTATION PROTOCOL ############################
#######################################################################################
###################################################
################ MUTATION PROFILES ##############
################
#  1) Create a stacked list of mutations based on all the probabilities and counts
#  2) Randomize or sort the stack
#  3) Execute the stack one at a time
#

        del Mutationlist[:]
        randval=random.random()
        if randval <.005:
            Mutationlist.insert(0,random.randint(1,edgefunctionstart-1)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
        elif randval <.1:
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
        elif randval <.2:
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
        elif randval <.3:
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
        elif randval <.4:
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,edgefunctionstart+3)) 
            Mutationlist.insert(0,random.randint(1,edgefunctionstart+3)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
        elif randval <.5:
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,edgefunctionstart+3)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,edgefunctionstart+3)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
        elif randval <.6:
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
        elif randval <.7:
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
        elif randval <.85:
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
        elif randval <.95:
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
        elif randval <.995:
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(1,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
            Mutationlist.insert(0,random.randint(edgefunctionstart,mutationfunctions)) 
        else:
            Mutationlist.insert(0,20) 
#################################################            
#################################################            
# Now, select the sequencing of the mutation
# Sorted places genes first, then edges
# Reverse sorted places edges first, then gene mutations
# Random is just random.
######
        randomval = random.random()      
        if randomval<=.4:                
            random.shuffle(Mutationlist)
        elif randomval<.7:
            Mutationlist.sort()
        else:
            Mutationlist.sort(reverse=True)
########################################################
################       MUTATION SCHEDULER
#### Though a random mutation pattern has been generated,  
#### it may be overwritten with mandated sequences
#### This is also where you would place generation-specific
#### operations.
######################################                                   
        randomval = random.random()      
        if randomval<=.001:                
            Mutationlist = [20,20,20,20,20,20,20,20]
        elif randomval<=.02:
            Mutationlist = [20,20,21,21,22]
        elif randomval<=.03:
            Mutationlist = [21,21,22]
#################################################################
#################################################################
###################  APPLY THE MUTATIONS  #######################
#################################################################
#################################################################                               
        for mutation in Mutationlist:       
            if mutation==1:  #mut_singleswap_sametruck_prob
                rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                rte1len = len(newroutes[i,rte1])-1 #must subtract 1 for randint
                loc1 = random.randint(0,rte1len)
                loc2 = random.randint(0,rte1len)
                firstval =newroutes[i,rte1][loc1]
                secondval=newroutes[i,rte1][loc2]
                newroutes[i,rte1][loc1]=secondval
                newroutes[i,rte1][loc2]=firstval

            elif mutation==2:  #mut_singleswap_differenttruck_prob
                rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                rte2 = random.randint(0,trucks-1)  #must subtract 1 for randint
                if rte1==rte2:                    #if picked same chromo, pick again
                    rte2 = random.randint(0,trucks-1) #must subtract 1 for randint                
                rte1len = len(newroutes[i,rte1])-1 #must subtract 1 for randint
                rte2len = len(newroutes[i,rte2])-1 #must subtract 1 for randint
                loc1 = random.randint(0,rte1len)
                loc2 = random.randint(0,rte2len)
                firstval = newroutes[i,rte1][loc1]
                newroutes[i,rte1][loc1]=newroutes[i,rte2][loc2]
                newroutes[i,rte2][loc2]=firstval           
                
            elif mutation==3:  #mut_singledeladd_prob
                rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                rte2 = random.randint(0,trucks-1)  #must subtract 1 for randint
                rte1len = len(newroutes[i,rte1])-1 #must subtract 1 for randint
                rte2len = len(newroutes[i,rte2])-1 #must subtract 1 for randint
                if rte1len>0:
                    loc1 = random.randint(0,rte1len)
                    firstval = newroutes[i,rte1].pop(loc1)
                    loc2 = random.randint(0,rte2len)
                    newroutes[i,rte2].insert(loc2,firstval)
                #print(i,rte1,loc1,rte2,loc2)                        
    
            elif mutation==4:   #mut_seqmove_sametruck_prob
                rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                rte1len = len(newroutes[i,rte1])-1 #must subtract 1 for randint
                if rte1len >4:
                    loc1 = random.randint(0,rte1len-1) #must allow for at least 1 character
                    loc2 = random.randint(loc1+1,rte1len) # end must be after start
                    savedvals = np.zeros((100),dtype=int)
                    savedindex = 0
                    #print(i,rte1,"loc1 and 2:",loc1,loc2)
                    for indx in range(loc1,loc2+1):
                        savedvals[savedindex] = newroutes[i,rte1].pop(loc1) 
                        savedindex = savedindex+1
                    #print(savedvals)
                    newrte1len = len(newroutes[i,rte1])-1 #must subtract 1 for randint
                    newloc1 = random.randint(0,rte1len-1)
                    for indx in range(0,savedindex):
                        newroutes[i,rte1].insert(newloc1+indx,savedvals[indx])
    
            elif mutation==5:   #mut_seqinvmove_sametruck_prob
                rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                rte1len = len(newroutes[i,rte1])-1 #must subtract 1 for randint
                if rte1len >4:
                    loc1 = random.randint(0,rte1len-1) #must allow for at least 1 character
                    loc2 = random.randint(loc1+1,rte1len) # end must be after start
                    savedvals = np.zeros((100),dtype=int)
                    savedindex = 0
                    for indx in range(loc1,loc2+1):
                        savedvals[savedindex] = newroutes[i,rte1].pop(loc1) 
                        savedindex = savedindex+1
                    newrte1len = len(newroutes[i,rte1])-1 #must subtract 1 for randint
                    newloc1 = random.randint(0,rte1len-1)
                    indx2 = 0
                    for indx in range(savedindex-1,-1,-1):
                        newroutes[i,rte1].insert(newloc1+indx2,savedvals[indx])
                        indx2=indx2+1
                    #print(i,"from:",rte1, loc1,"to:",rte2,loc2)
                    #print(savedvals)
    
            elif mutation==6:   #mut_seqmove_differenttruck_prob
                rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                rte1len = len(newroutes[i,rte1])-1 #must subtract 1 for randint
                if rte1len >4:
                    loc1 = random.randint(0,rte1len-1) #must allow for at least 1 character
                    loc2 = random.randint(loc1+1,rte1len) # end must be after start
                    savedvals = np.zeros((100),dtype=int)
                    savedindex = 0
                    if rte1len-(loc2+loc1)>3:
                        for indx in range(loc1,loc2+1):
                            savedvals[savedindex] = newroutes[i,rte1].pop(loc1) 
                            savedindex = savedindex+1
                        rte2 = random.randint(0,trucks-1)  #must subtract 1 for randint
                        if rte1==rte2:
                            rte2 = random.randint(0,trucks-1)  #must subtract 1 for randint                
                        rte2len = len(newroutes[i,rte2])-1 #must subtract 1 for randint
                        loc2 = random.randint(0,rte2len)
                        for indx in range(0,savedindex):
                            newroutes[i,rte2].insert(loc2+indx,savedvals[indx])
                    #print(i,"from:",rte1, loc1,"to:",rte2,loc2)
                    #print(savedvals)
    
            elif mutation==7:   #mut_seqinvmove_differenttruck_prob
                rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                rte1len = len(newroutes[i,rte1])-1 #must subtract 1 for randint
                if rte1len >4:
                    loc1 = random.randint(0,rte1len-1) #must allow for at least 1 character
                    loc2 = random.randint(loc1+1,rte1len) # end must be after start
                    savedvals = np.zeros((100),dtype=int)
                    savedindex = 0
                    if rte1len-(loc2+loc1)>3:
                        for indx in range(loc1,loc2+1):
                            savedvals[savedindex] = newroutes[i,rte1].pop(loc1) 
                            savedindex = savedindex+1
                        rte2 = random.randint(0,trucks-1)  #must subtract 1 for randint
                        if rte1==rte2:                 #try again
                            rte2 = random.randint(0,trucks-1)  #must subtract 1 for randint                
                        if rte1==rte2:                 #try again
                            rte2 = random.randint(0,trucks-1)  #must subtract 1 for randint                
                        rte2len = len(newroutes[i,rte2])-1 #must subtract 1 for randint
                        loc2 = random.randint(0,rte2len)
                        indx2 = 0
                        for indx in range(savedindex-1,-1,-1):
                            newroutes[i,rte2].insert(loc2+indx2,savedvals[indx])
                            indx2=indx2+1
                    #print(i,"from:",rte1, loc1,"to:",rte2,loc2)
                    #print(savedvals)

            elif mutation==8:   #move sequence from begin/end to diff truck begin/end
                #print(newroutes[i])
                if random.random()<.5:   #begin to begin
                    rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                    rte1len = len(newroutes[i,rte1]) 
                    movelen = random.randint(1,2) #either 2 long or 3 long sequence
                    notthistruck =list(range(0,trucks))
                    notthistruck.remove(rte1)          #remove selected orig from list
                    if rte1len >movelen+1:
                        loc1 = 0
                        savedvals = np.zeros((100),dtype=int)
                        rte2 = notthistruck[random.randint(0,len(notthistruck)-1)]  #get alternative truck
                        savedvals[0] = newroutes[i,rte1].pop(loc1) 
                        savedvals[1] = newroutes[i,rte1].pop(loc1) 
                        if movelen==2: savedvals[2] = newroutes[i,rte1].pop(loc1) #three stop sequence
                        if movelen==2: newroutes[i,rte2].insert(0,savedvals[2]) #keep the same order
                        newroutes[i,rte2].insert(0,savedvals[1])
                        newroutes[i,rte2].insert(0,savedvals[0])
                    #print("begin to begin",i,"Movelen:",movelen," From:",rte1,rte2)
                    #print(newroutes[i])
                else:    #end to end
                    #print(newroutes[i])
                    rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                    rte1len = len(newroutes[i,rte1])-1  #subtract 1 for end position
                    movelen = random.randint(1,2) #either 2 long or 3 long sequence
                    notthistruck =list(range(0,trucks))
                    notthistruck.remove(rte1)          #remove selected orig from list
                    if rte1len>movelen:
                        loc1 = 0
                        savedvals = np.zeros((100),dtype=int)
                        rte2 = notthistruck[random.randint(0,len(notthistruck)-1)]  #get alternative truck
                        if movelen==2: 
                            savedvals[2] = newroutes[i,rte1].pop(rte1len) #three stop sequence
                            savedvals[1] = newroutes[i,rte1].pop(rte1len-1) 
                            savedvals[0] = newroutes[i,rte1].pop(rte1len-2) 
                        else:
                            savedvals[1] = newroutes[i,rte1].pop(rte1len) 
                            savedvals[0] = newroutes[i,rte1].pop(rte1len-1) 
                        rte2len = len(newroutes[i,rte2])-1  #subtract 1 for end position
                        if movelen==2: newroutes[i,rte2].insert(rte2len,savedvals[2]) #keep the same order
                        newroutes[i,rte2].insert(rte2len+1,savedvals[1])
                        newroutes[i,rte2].insert(rte2len+1,savedvals[0])                
                    #print("End to end",i,"Movelen:",movelen," From:",rte1,rte2)
                    #print(newroutes[i])

            elif mutation==9:   #move seq from begin to inv end OR from end to inv begin
                #print(newroutes[i])
                if random.random()<.5:   #begin to inv end
                    rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                    rte1len = len(newroutes[i,rte1]) 
                    movelen = random.randint(1,2) #either 2 long or 3 long sequence
                    notthistruck =list(range(0,trucks))
                    notthistruck.remove(rte1)          #remove selected orig from list
                    if rte1len >movelen+1:
                        loc1 = 0
                        savedvals = np.zeros((100),dtype=int)
                        rte2 = notthistruck[random.randint(0,len(notthistruck)-1)]  #get alternative truck
                        rte2len=len(newroutes[i,rte2])-1  #get last position in new truck
                        savedvals[0] = newroutes[i,rte1].pop(loc1)
                        savedvals[1] = newroutes[i,rte1].pop(loc1) 
                        if movelen==2: savedvals[2] = newroutes[i,rte1].pop(loc1) #three stop sequence
                        newroutes[i,rte2].insert(rte2len+1,savedvals[0]) #inverse the sequence
                        newroutes[i,rte2].insert(rte2len+1,savedvals[1])
                        if movelen==2: newroutes[i,rte2].insert(rte2len+1,savedvals[2]) 
                    #print("Begin to End",i,"Movelen:",movelen," From:",rte1,rte2)
                    #print(newroutes[i])
                else:   #inv end to begin
                    rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                    rte1len = len(newroutes[i,rte1])-1  #subtract 1 for end position
                    movelen = random.randint(1,2) #either 2 long or 3 long sequence
                    notthistruck =list(range(0,trucks))
                    notthistruck.remove(rte1)          #remove selected orig from list
                    if rte1len>movelen:
                        loc1 = 0
                        savedvals = np.zeros((100),dtype=int)
                        rte2 = notthistruck[random.randint(0,len(notthistruck)-1)]  #get alternative truck
                        if movelen==2: 
                            savedvals[2] = newroutes[i,rte1].pop(rte1len) #three stop sequence
                            savedvals[1] = newroutes[i,rte1].pop(rte1len-1) 
                            savedvals[0] = newroutes[i,rte1].pop(rte1len-2) 
                        else:
                            savedvals[1] = newroutes[i,rte1].pop(rte1len) 
                            savedvals[0] = newroutes[i,rte1].pop(rte1len-1) 
                        newroutes[i,rte2].insert(0,savedvals[0])                
                        newroutes[i,rte2].insert(0,savedvals[1])
                        if movelen==2: newroutes[i,rte2].insert(0,savedvals[2]) #keep the same order
                    #print("End to begin",i,"Movelen:",movelen," From:",rte1,rte2)
                    #print(newroutes[i])


            elif mutation==10:   #swap beginning and inv end - diff trucks
                rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                rte1len = len(newroutes[i,rte1]) 
                movelen = random.randint(1,2) #either 2 long or 3 long sequence
                if rte1len<=movelen+1:
                    rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                    rte1len = len(newroutes[i,rte1]) 
                if rte1len>movelen+1:
                    notthistruck =list(range(0,trucks))
                    notthistruck.remove(rte1)          #remove selected orig from list
                    rte2 = notthistruck[random.randint(0,len(notthistruck)-1)]  #get alternative truck
                    rte2len = len(newroutes[i,rte2])-1  #subtract 1 for end position 
                    if rte2len<=movelen:
                        rte2 = notthistruck[random.randint(0,len(notthistruck)-1)]  #get alternative truck
                        rte2len = len(newroutes[i,rte2])-1  #subtract 1 for end position 
                    if rte2len>movelen:
                        loc1 = 0
                        savedvals = np.zeros((100),dtype=int)
                        savedvals[0] = newroutes[i,rte1].pop(loc1)
                        savedvals[1] = newroutes[i,rte1].pop(loc1)
                        if movelen==2: savedvals[2] = newroutes[i,rte1].pop(loc1) #three stop sequence

                        savedvals[10] = newroutes[i,rte2].pop(rte2len)  #inverse the order for convenience
                        savedvals[11] = newroutes[i,rte2].pop(rte2len-1)
                        if movelen==2: savedvals[12] = newroutes[i,rte2].pop(rte2len-2) #three stop sequence
 
                        if movelen==2: newroutes[i,rte2].insert(rte2len,savedvals[2])  #insert 1 into 2
                        newroutes[i,rte2].insert(rte2len,savedvals[1])
                        newroutes[i,rte2].insert(rte2len,savedvals[0])

                        if movelen==2: newroutes[i,rte1].insert(0,savedvals[12]) #keep the same order
                        newroutes[i,rte1].insert(0,savedvals[11])
                        newroutes[i,rte1].insert(0,savedvals[10])
                    #print(i,"Movelen:",movelen," From:",rte1,rte2)
                
            elif mutation==20:   #mut_getcloseone_prob
                rte1 = random.randint(0,trucks-1)  #must subtract 1 for randint
                rte1len = len(newroutes[i,rte1])-1 #must subtract 1 for randint
                if rte1len==0:
                    loc1=0
                else:               
                    loc1 = random.randint(0,rte1len-1) #must allow for at least 1 character
                genenum = newroutes[i,rte1][loc1]  #get the gene num at the destination
                myproxrecord=mydata[mydata.genenum==genenum]
                howcloserandom = random.random()
                if howcloserandom <=.5:    #weightings for proximity of item
                    newzipid=myproxrecord['closest1'].iloc[0]
                elif howcloserandom <=.8:
                    newzipid=myproxrecord['closest2'].iloc[0]
                elif howcloserandom <=.95:
                    newzipid=myproxrecord['closest3'].iloc[0]
                else:
                    newzipid=myproxrecord['closest4'].iloc[0]
                myproxdestrecord=mydata[mydata.tozipid==newzipid]
                newgenenum=myproxdestrecord['genenum'].iloc[0]
                #Now locate the close-by destination gene for extraction
                for rte2 in range(0,trucks):
                    index = newroutes[i,rte2].index(newgenenum) if newgenenum in newroutes[i,rte2] else -1
                    if index > -1:
                        break
                rte2len = len(newroutes[i,rte2])
                if rte2len > 1 and index>-1:
                    newroutes[i,rte2].remove(newgenenum)
                    beforeorafter=int(random.random()+.5)
                    loc2 = loc1+beforeorafter
                    newroutes[i,rte1].insert(loc2,newgenenum)
            # print(i,"new location",rte1, loc1,"from:",rte2,index,"\n")
    
            elif mutation==21:   #mut_consolidatemultiple_prob
                foundit=0
                randommultiple = random.randint(0,multiplescount-1) 
                mymultiple=multiples.iloc[randommultiple]
                multzipid=mymultiple['tozipid']
                multgenenum=mymultiple['genenum']
                nextmultiple=multiples[(multiples.tozipid==multzipid) & (multiples.genenum!=multgenenum)]
                nextmultgenenum=int(nextmultiple['genenum'].iloc[0])
                for rte1 in range(0,trucks):
                    index = newroutes[i,rte1].index(multgenenum) if multgenenum in newroutes[i,rte1] else -1
                    if index > -1:
                        rte1len = len(newroutes[i,rte1])
                        if rte1len>1:    #require more than one order
                            foundit=1                        
                            newroutes[i,rte1].remove(multgenenum)
                        else:
                            foundit=-1   #only 1 order, so denied request                     
    
                if foundit==1:
                    for rte1 in range(0,trucks):
                        index = newroutes[i,rte1].index(nextmultgenenum) if nextmultgenenum in newroutes[i,rte1] else -1
                        if index > -1:
                            newroutes[i,rte1].insert(index,multgenenum)
    
            elif mutation==22:   #mut_consolidatetriple_prob
                if triplescount>0:
                    newroutescopy = deepcopy(newroutes[i,])
                    randommultiple = random.randint(0,triplescount-1) 
                    mytriple=triples.iloc[randommultiple]
                    multzipid=mytriple['tozipid']
                    triplegenenum=mytriple['genenum']
                    nexttriple=triples[(triples.tozipid==multzipid) & (triples.genenum!=triplegenenum)]
                    nexttriplegenenum0=int(nexttriple['genenum'].iloc[0])
                    nexttriplegenenum1=int(nexttriple['genenum'].iloc[1])
                    for rte1 in range(0,trucks):   #delete the 2nd and 3rd ones first                    
                        index = newroutes[i,rte1].index(nexttriplegenenum0) if nexttriplegenenum0 in newroutes[i,rte1] else -1
                        if index > -1:
                            newroutes[i,rte1].remove(nexttriplegenenum0)
        
                        index = newroutes[i,rte1].index(nexttriplegenenum1) if nexttriplegenenum1 in newroutes[i,rte1] else -1
                        if index > -1:
                            newroutes[i,rte1].remove(nexttriplegenenum1)
    
                    triplefail =0
                    for rte1 in range(0,trucks): 
                        if len(newroutes[i,rte1])==0:
                            triplefail = 1
                            
                    if triplefail==1:
                        newroutes[i,]=deepcopy(newroutescopy)
                    else:
                        for rte1 in range(0,trucks):    #now add them back in.
                            newindex = newroutes[i,rte1].index(triplegenenum) if triplegenenum in newroutes[i,rte1] else -1
                            if newindex > -1:
                                newroutes[i,rte1].insert(index,nexttriplegenenum0)
                                newroutes[i,rte1].insert(index,nexttriplegenenum1)
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###########################  Evaluate the Progeny  ############################
###############################################################################
## 1) Calculate the fitness of the child's plan
## 2) Decide if child replaces parent
##    a) if child is more fit OR
##    b) if searchexpansion is enabled, then child may succeed parent anyway.
##                   
        newplans[i]=fitness_function(newroutes[i,],mydata,mydistances,generation_ctr,0)  #get fitness of the routes
        attemptsatoffspring[i] = attemptsatoffspring[i] +1
        if newplans[i] < plans[i]:              #compare offspring to parent
            plans[i] = deepcopy(newplans[i])
            routes[i,] = deepcopy(newroutes[i,])
            attemptsatoffspring[i] = 0
            if plans[i] < bestplanvalue[i]:
                bestplanvalue[i] = deepcopy(plans[i])
                bestroute[i,] = deepcopy(routes[i,])
                bestplangeneration[i] = generation_ctr
                if plans[i] < bestplanever[weekday, trucks]:
                    bestplanever[weekday, trucks] = deepcopy(plans[i])
                    bestrouteever[weekday, trucks] = deepcopy(routes[i,])
                    if loadbestplanever!=-2:    
                        afile = open(r'GS490bestplan.pkl', 'wb')
                        pickle.dump(bestplanever, afile)
                        afile.close()
                        bfile = open(r'GS490bestroute.pkl', 'wb')
                        pickle.dump(bestrouteever, bfile)
                        bfile.close()                    
            
        elif enablesearchexpansion == 1:          #proceed only if search expansion enabled
            plandifference = newplans[i]-plans[i]  #how much does child vary from parent
            if plandifference == 0:               #if child is just as good as parent, then use child
                plans[i] = deepcopy(newplans[i])
                routes[i,] = deepcopy(newroutes[i,])
                attemptsatoffspring[i] = 0                
                searchexpdisplay = searchexpdisplay + "•" + str(i)
            elif attemptsatoffspring[i]>searchexpthold and enablesearchupward ==1:  #if worse & search upward
                if plandifference < searchexppenalty*attemptsatoffspring[i]:
                    plans[i] = deepcopy(newplans[i])
                    routes[i,] = deepcopy(newroutes[i,])
                    attemptsatoffspring[i] = 0                
                    searchexpdisplay = searchexpdisplay + "•P:" + str(i)  #Display a P 
                    
###################### END OF POPULATION GENERATION LOOP ###################### 
###############################################################################                      
### Calculate relative fitness proportion
    totalfitnessinpop = plans.sum(axis=0)   #sum total fitness in population
    planfitnessproportion = totalfitnessinpop/plans  #this is inverse since we are minimizing
    TopItem = plans[0:pop_size].argmin(axis=0)
    TopPlan = plans[TopItem]
    if TopItem<pop_size-1:   #one of core plans
        routes[pop_size-1,]= deepcopy(routes[TopItem,])
    else:    #duplicate of top item
        routes[PreviousTopItem,]= deepcopy(routes[TopItem,])
        plans[PreviousTopItem]=plans[TopItem]
        TopPlan = plans[TopItem]
        TopItem = PreviousTopItem
    PreviousTopItem=TopItem
#####################################################################
###########################  DISPLAY STATUS  ########################        
#####################################################################   
    if (display_EveryRecord == 1) or (prevplanlist.sum()!=plans.sum()):   #changed plans or heavy logging.
        if display_MutationInfo==1:
            print("\n",generation_ctr,TopItem,str(plans),Mutationlist," ",end="")
        else:
            print("\n",generation_ctr,TopItem,str(plans)," ",end="")
        prevplanlist = deepcopy(plans)
    if display_searchexpansion==1:
        print(searchexpdisplay, end="")
        searchexpdisplay=""
    if PreviousTopPlan != TopPlan:  #write plans to file
        PreviousTopPlan = TopPlan
#####################################################################
########################  SAVING DATA TO FILE  ######################        
#####################################################################   
# Code not designed.  
#CHANGE THIS TO A FLAG THAT ENABLES PRINTING
      #  outstr = '\n' + str(generation_ctr) + '\t' + str(routes)
      #  f_r.write(outstr)
      #  outstr = '\nGen:' + str(generation_ctr) + '   TopItem:' + str(TopItem) + ':' + str(TopPlan) + str(plans)
      #  f_hs.write(outstr)
                    
    #if generation_ctr/50 ==int(generation_ctr/50):
       # a=fitness_function(routes[TopItem,],mydata,Distances,zero,1)
#       
#Code for exporting data structures.  Use pickle to serialize the data 
#import pickle
#afile = open(r'routes_Thursday7-25plans.pkl', 'wb')
#pickle.dump(routes, afile)
#afile.close()
#
#reload object from file
#file2 = open(r'C:\d.pkl', 'rb')
#new_d = pickle.load(file2)
#file2.close()
#####################################################################
#####################################################################
#####################################################################
### If you halt the execution, you can use the following lines to print the best solution:
###
print("Gen:",generation_ctr,"   MeanFitness:",totalfitnessinpop/pop_size,"   TopItem:",TopItem,"(",TopPlan,")")
myitem= [TopItem]
for i in myitem:
    a=fitness_function(routes[i,],mydata,mydistances,zero,1)
#####################################################################
#####################################################################
#####################################################################
#Done

