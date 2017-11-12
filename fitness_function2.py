def fitness_function(planinfo,mydata,distances,generations,onoff):
#planinfo is the n chromosome set of routes
#mydata is the database of order data for this day
#distances is a database with distance information
#generations is simply the generation counter for printing
#onoff is a flag as to whether console display is on or off, 1=print, 0=noprint
    ########## SET CRITICAL VARIABLES
    #################################
    overnightpenaltyenabled = 0
    overnightpenaltybase = 50     #minimum # of miles before an overnight stay is pursued
    #if penalty base above set to 100, then overnight stays will be drastically curtailed.
    #if penalty set to 0, then a reduction of 10 miles can result in an overnight.
    volumepenaltybase=100000      #This is a constraint, so penalty must exceed any other value
                                  #due to randomness at startup.  
    traveltimepenaltybase =500    #set to 500; opportunity for optimization
#########\/  \/  \/ Added 05.08.2017  \/  \/  \/  \/  \/
    permitelevenhourdrive = 1     #allow 11 hour drive, but add 10 hour break.
    elevenhourpenalty = 200      #This is when the 11 hour rule simply requires a 10 hour stop on the return drive.
    Noelevenhourdrivespenalty=10000 #This is when the 11 hour rule creates an unacceptable plan; set to 10000
#########/\  /\  /\  End of add /\  /\  /\  /\  /\  /\
    #################################
    ###### END OF CRITICAL VARIABLES
    #################################    
    trucks=len(planinfo)
    generationsfactor = int(generations)          #experimental idea to change fitness function by gen count
    ######## I think this is the wrong approach - a better approach is to pass the fitness value
    ######## and if the fitness value is high (i.e., not very fit), then alter the algorithm.
    ######## or could calculate the ratio of the fitness to the most fit.  If the ratio shows that we
    ######## have a bad plan, then zero out the penalties and let the routes mutate.  
    #####################################################################################
    routefitvalue = np.zeros((trucks),dtype=float)
    destgene = np.zeros((100,trucks),dtype=int)
    dest = np.zeros((100,trucks),dtype=int)
    volume = np.zeros((100,trucks),dtype=int)
    traveldist = np.zeros((100,trucks),dtype=int)
    travelmin = np.zeros((100,trucks),dtype=float)
    unloadmin = np.zeros((100,trucks),dtype=float)
    unloadminraw = np.zeros((100,trucks),dtype=float)
    cumtravelmin = np.zeros((100,trucks),dtype=float)
    cumunloadmin = np.zeros((100,trucks),dtype=float)
    cumdisttrav = np.zeros((100,trucks),dtype=float)
    cumvolume = np.zeros((100,trucks),dtype=float)
    fromzip = np.zeros((100,trucks),dtype=float)
    tozip = np.zeros((100,trucks),dtype=float)

    totalvalue = 0
    totaldist=0
    if onoff>0:
        print("\ntruck,stop,orderid,from,to,travelmin,cumtravel,volume,cumvolume,unloadmin,cumunload,miles,cummiles,overntmin",end="")    
       # print("truck\tstop\torderid\tfrom\tto\ttravelmin\tcumtravel\tvolume\tcumvolume\tunloadmin\tcumunload\tmiles\tcummiles")    

    for j in range (0,trucks):
        index=0
        for routeitem in (planinfo[j]):
            destgene[index,j] = routeitem
            myrecord=mydata[mydata.genenum==routeitem]

            thisorderid=myrecord['orderid'].iloc[0]
            thiszipid=myrecord['tozipid'].iloc[0]
            milesfromdc=myrecord['milesfromdc'].iloc[0]
            thisvolume=myrecord['cube'].iloc[0]
            minutesfromdc=myrecord['minutesfromdc'].iloc[0]
            thisunloadtime=myrecord['minutesunload'].iloc[0]
            thisunloadtimeraw=myrecord['minutesunloadraw'].iloc[0]

            dest[index,j] = thiszipid #get actual destination number
            volume[index,j] = thisvolume
            unloadmin[index,j] = thisunloadtime
            unloadminraw[index,j] = thisunloadtimeraw            
            if index==0:   #for first destination
                traveldist[index,j] = milesfromdc #use pre-computed distance from 1887
                travelmin[index,j] = minutesfromdc
                firsttravelmin = minutesfromdc
                previousid = thiszipid
                cumtravelmin[index,j]=travelmin[index,j]
                cumunloadmin[index,j]=unloadmin[index,j]
                cumdisttrav[index,j]=traveldist[index,j]
                fromzip[index,j]=-1
                tozip[index,j]=thiszipid
                cumvolume[index,j]=thisvolume
            if index!=0:
                if previousid==thiszipid:
                    thismiles=0    
                else:
                    myrecord2=distances[(distances.fromid==previousid) & (distances.toid==thiszipid) ]
                    thismiles=myrecord2['distance'].iloc[0]
                traveldist[index,j] = thismiles #use distance in database
                travelmin[index,j] = thismiles*1.5   #60/40 =1.5 or minutes per mile
                samedest = 0
                if dest[index,j]==dest[index-1,j]:
                    samedest=2
                    if j>1:
                        if dest[index,j]==dest[index-2,j]:
                            samedest=3                            
                if samedest==2:
                    sumunload = unloadminraw[index-1,j] + unloadminraw[index,j]
                    if sumunload<=30:
                        unloadmin[index,j] = 0
                    else:
                        unloadmin[index,j] = sumunload-unloadminraw[index,j]
                if samedest==3:
                    sumunload = unloadminraw[index-1,j] + unloadminraw[index,j] + unloadminraw[index-2,j]
                    if sumunload<=30:
                        unloadmin[index,j] = 0
                    else:
                        unloadmin[index,j] = sumunload-unloadminraw[index,j]
                    
                cumvolume[index,j]=thisvolume+cumvolume[index-1,j]
                cumtravelmin[index,j]=travelmin[index,j]+cumtravelmin[index-1,j]
                cumunloadmin[index,j]=unloadmin[index,j]+cumunloadmin[index-1,j]
                cumdisttrav[index,j]=traveldist[index,j]+cumdisttrav[index-1,j]
                fromzip[index,j]=previousid
                tozip[index,j]=thiszipid
                previousid = thiszipid
            if onoff==1:
               # print(j,"\t",index,"\t",thisorderid,"\t",int(fromzip[index,j]),"\t",int(tozip[index,j]),"\t",travelmin[index,j],"\t",cumtravelmin[index,j],"\t",thisvolume,"\t",cumvolume[index,j],"\t",unloadmin[index,j],"\t",cumunloadmin[index,j],"\t",traveldist[index,j],"\t",cumdisttrav[index,j])
                print("\n",j,",",index,",",thisorderid,",",int(fromzip[index,j]),",",int(tozip[index,j]),",",travelmin[index,j],",",cumtravelmin[index,j],",",thisvolume,",",cumvolume[index,j],",",unloadmin[index,j],",",cumunloadmin[index,j],",",traveldist[index,j],",",cumdisttrav[index,j],end="")
            index=index+1

        dest[index,j]=-1    #last destination is distribution center
        volume[index,j] = 0 #no volume to unload
        thisvolume = 0
        unloadmin[index,j] = 0 #no unload time
        traveldist[index,j] = milesfromdc #use previous record travel distance
        travelmin[index,j] = minutesfromdc #use previous minutes record
        lasttravelmin = travelmin[index,j]

        cumvolume[index,j]=thisvolume+cumvolume[index-1,j]
        cumtravelmin[index,j]=travelmin[index,j]+cumtravelmin[index-1,j]
        cumunloadmin[index,j]=unloadmin[index,j]+cumunloadmin[index-1,j]
        cumdisttrav[index,j]=traveldist[index,j]+cumdisttrav[index-1,j]
        fromzip[index,j]=previousid
        tozip[index,j]=-1
        if onoff>0:
            print("\n",j,",",index,",",thisorderid,",",int(fromzip[index,j]),",",int(tozip[index,j]),",",travelmin[index,j],",",cumtravelmin[index,j],",",thisvolume,",",cumvolume[index,j],",",unloadmin[index,j],",",cumunloadmin[index,j],",",traveldist[index,j],",",cumdisttrav[index,j],end="")
        
        travelminsum = travelmin.sum(axis=0)[j] 
        unloadminsum = unloadmin.sum(axis=0)[j] 
        traveldistsum = traveldist.sum(axis=0)[j] 
        volumesum = volume.sum(axis=0)[j] 
        middletravelmin = travelminsum - firsttravelmin - lasttravelmin
        totalmiddlemin = middletravelmin + unloadminsum
        completeroutemin = travelminsum + unloadminsum
        if firsttravelmin > 240:            
            totalmiddlemin = totalmiddlemin + (firsttravelmin-240)
            firsttravelmin = 240
        penaltymin = 0
        if totalmiddlemin > 600:
            penaltymin = (totalmiddlemin - 600)
        penaltyvolume = 0
        if volumesum>3200:
            penaltyvolume = (volumesum-3200)
        overnightmin = completeroutemin - (14*60)
        overnightpenalty=0

#### NEw code 05.08.2017 ####
# Find out how many travel hours within the 14 hour DOT workday.
        if overnightmin >0:
            daytimetravelmin = travelminsum - overnightmin
        else:
            daytimetravelmin = travelminsum
        if daytimetravelmin <= 11*60:    #restricted to 11 hours of travel in 14 hour day
            routefitvalue[j] = traveldistsum
        else:
            if permitelevenhourdrive==0:
                routefitvalue[j] = traveldistsum + Noelevenhourdrivespenalty
            else:
                minutesbeforelasttravel = travelminsum-lasttravelmin
                if minutesbeforelasttravel>11*60:   #Overage comes earlier, so bad route - must decline
                    routefitvalue[j] = traveldistsum + Noelevenhourdrivespenalty
                else:     #overage is in the last travel bit back
                    overageminutes = travelminsum-660
                    if overageminutes > overnightmin:  #is overage more from driving or more from 14 hours 
                        overnightmin = overageminutes  #if driving, then reset overnight
                    routefitvalue[j] = traveldistsum
##################End new code#######################
        if onoff>0:
            if(overnightmin>0): print(",",overnightmin,end="")
            print("\n",end="")
        
        if overnightmin > 0 and overnightpenaltyenabled==1:
            overnightpenalty = overnightpenaltybase + 2*overnightmin

        routefitvalue[j] = routefitvalue[j] + overnightpenalty

        if penaltyvolume > 0:                       #Must add in volume penalty
            routefitvalue[j] = routefitvalue[j] + penaltyvolume*2+volumepenaltybase
        if penaltymin > 0:
            routefitvalue[j] = routefitvalue[j] + penaltymin*2+traveltimepenaltybase
        

        totaldist =traveldistsum + totaldist       #This is just distance fitness value
        totalvalue = totalvalue + routefitvalue[j]  
        
    if onoff>0:
        print("")
        print("Fitness: ", totalvalue,"TotDist: ,",totaldist)
        print("")

    return totalvalue
