########################################################
########################################################
"""
Coding started 11/27/2017
by Steve Reedy
"""
########################################################
########################################################

#import modules section

import os #needed to move around filesystem
import csv #needed to read and write csv files
from cvxopt import matrix, solvers #matrix needed to create input matrices to solver, solver needed to optimize
import time #to check run times
import sys # to get arguments from command line, check python version
import glob #to get names of files in directories
import gc #to clean up memory
import re
import subprocess
from mipcl_py.mipshell.mipshell import *

#########################################################
#########################################################



#########################################################
#########################################################
"""There are three main data structures in this program.
permdata - a dictionary that contains (mostly) data that is used by each day 
           in the program; namely
              system conditions
              resource parameters
              under generation and over generation penalty factors
              load as requirements

daydata - a dictionary that contains the things needed for a particular day; namely
              60 day resource file
              cdr 12302
              cdr 12354
              output data

timestampdata - a dictionary that contains the things needed for a particular timestamp;
              P,q,G,h,A,b matrices
              optimization results
"""

def get_system_conditions(wdname):
    """
    System Conditions File Format
          0- timestamp
          1- repeated hour flag
          2- Reg up AS obl for that interval
          3- Reg down AS obl for that interval
          4- RRS AS obl for that interval
          5- NSRS AS obl for that interval
          6- Reg up AS price for that interval
          7- Reg down AS price for that interval
          8- RRS AS price for that interval
          9- NSRS AS price for that interval
          10-system lambda for that interval
          11-GTBD for that interval

    """
    os.chdir(wdname) #changes working directory to input working directory
    with open('system_conditions.csv', 'r') as csvfile: #open the system condition file
        filereader = csv.reader(csvfile)
        permdata = {}
        systemconditions = {}
        for row in filereader:
            if row[0] != 'SCED_TIMESTAMP':
                timestamp = (row[0],row[1])
                systemconditions[timestamp] = {}
                systemconditions[timestamp]['RegDownMW'] = float(row[2])
                systemconditions[timestamp]['RegUpMW'] = float(row[3])
                systemconditions[timestamp]['RRSMW'] = float(row[4])
                systemconditions[timestamp]['NSRSMW'] = float(row[5])
                systemconditions[timestamp]['RegDownPrice'] = float(row[6])
                systemconditions[timestamp]['RegUpPrice'] = float(row[7])
                systemconditions[timestamp]['RRSPrice'] = float(row[8])
                systemconditions[timestamp]['NSRSPrice'] = float(row[9])
                systemconditions[timestamp]['SystemLambda'] = float(row[11])
                systemconditions[timestamp]['GTBD'] = float(row[10])
                systemconditions[timestamp]['LoadNSRS'] = 0.0
                systemconditions[timestamp]['LoadRRS'] = 0.0
                systemconditions[timestamp]['LoadRegUp'] = 0.0
                systemconditions[timestamp]['OutNSRS'] = 0.0
                systemconditions[timestamp]['OutRRS'] = 0.0
                systemconditions[timestamp]['OutRegUp'] = 0.0
        permdata['SystemConditions'] = systemconditions
    return permdata

def get_quickstarts(wdname,permdata):
    os.chdir(wdname) #changes working directory to input working directory
    with open('quickstarts.csv', 'r') as csvfile: #open the system condition file
        filereader = csv.reader(csvfile)
        quickstarts = []
        for quickstart in filereader:
            quickstarts.append(quickstart[0])
    permdata['QuickStarts'] = quickstarts
    return permdata

def write_log_entry(entry,wdname,permdata):
    parentdir = os.path.dirname(wdname)
    resultsdir = os.path.join(parentdir,'nfrcresults')
    os.chdir(resultsdir)
    filename = 'log' + os.path.basename(wdname) + '.csv'
    if permdata['PythonVersion'] == 3:
        with open(filename, 'w',newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            spamwriter.writerow(entry)

    else:
        with open(filename, 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            spamwriter.writerow(entry)



def get_load_as(wdname,permdata):
    os.chdir(wdname) #changes working directory to input working directory
    filename = glob.glob('60d_Load*.csv')[0] #get name of load resource file in directory
    with open(filename, 'r') as csvfile: #open the system condition file
        filereader = csv.reader(csvfile)    
        for row in filereader:
            if row[0] != 'SCED Time Stamp':
                timestamp = (row[0],row[1])
                try:
                    permdata['SystemConditions'][timestamp]['LoadNSRS'] += float(row[8])
                    permdata['SystemConditions'][timestamp]['LoadRRS'] += float(row[7])
                    permdata['SystemConditions'][timestamp]['LoadRegUp'] += float(row[9])
                    if row[3] == 'OUTL':
                        permdata['SystemConditions'][timestamp]['OutNSRS'] += float(row[8])
                        permdata['SystemConditions'][timestamp]['OutRRS'] += float(row[7])
                        permdata['SystemConditions'][timestamp]['OutRegUp'] += float(row[9])
                except KeyError:
                    pass
    return permdata


def get_resource_parameters(wdname):
    """
    Resource Parameters File Format
    0- Resource
    1- Reg up AS obl Historical MWhs trailing twelve months
    2- Reg down AS obl Historical MWhs -trailing twelve months
    3- RRS AS obl Historical MWhs -trailing twelve months
    4- NSRS AS obl Historical MWhs -trailing twelve months
    5- Ramp Rate Up -not used, no data
    6- Ramp Rate Down -not used, no data
    7- NFRC -max in last year
    """
    rp = {} #temp variable
    os.chdir(wdname) #changes working directory to input working directory
    filename = glob.glob('resource_parameters*.csv')[0] #get name of resource pararmeters file in directory
    with open(filename, 'r') as csvfile: #open the resource parameter file
        filereader = csv.reader(csvfile)
        for row in filereader:
            if row[0] != 'RES_NAME':
                resource = row[0]
                rp[resource] = {}
                rp[resource]['RegUpHMWh'] = float(row[1])
                rp[resource]['RegDownHMWh'] = float(row[2])
                rp[resource]['RRSHMWh'] = float(row[3])
                rp[resource]['NSRSHMWh'] = float(row[4])
                rp[resource]['RampRateUp'] = float(row[5])
                rp[resource]['RampRateDown'] = float(row[6])
                rp[resource]['NFRC'] = float(row[7])
    return rp

def open_days(filename): #opens daylist file and gets list of all possible days to optimize
    with open(filename, 'r') as csvfile:
        import pdb; pdb.set_trace()
        filereader = csv.reader(csvfile)
        dayfolders  =[]
        count = 0
        for row in filereader:
            dayfolders.append((row[0],count))
            count += 1
    return dayfolders

def get_grd_file(wdname,resourceparameters,permdata,dayfolder): #gather the data in the 60day Sced GRD report
    themonth = int(dayfolder[4:6])
    theday = int(dayfolder[6:8])
    prequickstartstatus = themonth <= 2 or (themonth == 3 and theday <= 9)
    startload = time.time()
    quickstarts = permdata['QuickStarts']
    #for quickstart in quickstarts: print(quickstart) #troubleshooting
    daydata = {}
    daydata['ResourceParameters'] = resourceparameters
    daydata['Log'] = []
    daydata['Summary'] = [] #initialize the summary file data
    daydata['PriceCheck'] = [] #initialize the price check  file
    daydata['RRSCommit'] = {} #initialize RRSCommit info
    os.chdir(wdname)
    badstatuslist = ['OUT','EMR','OFF','OFFNS','Telemetered Resource Status'] #shouldn't be included in optimization
    daydata['TimeStamps'] = [] #initialize list of all timestamps for the dictionary
    daydata['60DayGRD'] = [] #initialize output file
    filename = glob.glob('60d_SCED_Gen*.csv')[0] #get name of 60 day sced report
    with open(filename, 'r') as csvfile: #open the 60 day SCED report
        offerfilereader = csv.reader(csvfile)
        for row in offerfilereader: #get a row of data
            if row[0] == 'SCED Time Stamp': #if header row, add new column headers
                row.append('SCED_COST')
                row.append('SCED_BP')
                row.append('SCED_LMP')
                row.append('MIP_RRS_COMMIT')
                row.append('COOPT_COST')
                row.append('COOPT_BP')
                row.append('COOPT_NSRS_AWARD')
                row.append('COOPT_RRS_AWARD')
                row.append('COOPT_REG_UP_AWARD')
                row.append('COOPT_LMP')
                row.append('COOPT_NSRS_PRICE')
                row.append('COOPT_RRS_PRICE')
                row.append('COOPT_REG_UP_PRICE')
                row.append('COOPT_HSL-LDL')
                row.append('COOPT_RRS_OFFER_CAP')
                row.append('COOPT_REG_UP_OFFER_CAP')
                row.append('COOPT_RAMP_UP')
                row.append('COOPT_NFRC')
                row.append('COOPT_HDL')
                row.append('E_PLUS_5_SEVENTHS_REG')
                row.append('RRS_PLUS_REG_PLUS_OFFER_CAP')
                row.append('NSRS_PLUS_RRS_PLUS_REG_OFFER_CAP')
            else: #otherwise, add the time stamp to the timestamp set
                timestamp = (row[0],row[1])
                if timestamp not in daydata['TimeStamps']:
                    daydata['TimeStamps'].append(timestamp) # daydata['TimeStamps'] is a list of all the timestamps in the day
                    daydata[timestamp] = {} # initialize dictionary with timestamp data
                resource = row[2]
                daydata[timestamp][resource] = {}
                if row[151] not in badstatuslist:
                    daydata[timestamp][resource]['OriginalRRSAward'] = float(row[156]) >= 0.1 
                    daydata[timestamp][resource]['NFRC'] = daydata['ResourceParameters'][resource]['NFRC']
                    if prequickstartstatus and row[2] in quickstarts and row[151] == 'ON' and float(row[148]) < 1.0: #check if this is a quickstart unit operating in quickstart mode before March 10th
                        row[151] = 'OFFQS'
                        #print(row) #for troubleshooting
            daydata['60DayGRD'].append(row) #build output file
    return daydata

def get_binding_constraint_data(wdname,daydata): #get cdr12302 data
    os.chdir(wdname)
    daydata['NetworkConstraints'] = []
    for filename in glob.glob('cdr.00012302.*.csv'):
        print(filename)
        daydata['Log'].append(filename)
        with open(filename, 'r') as csvfile: #open the binding constraint report
            filereader = csv.reader(csvfile)
            for row in filereader: #get a row of data
                timestamp = (row[0],row[1])
                if timestamp in daydata['TimeStamps'] : #check to see if relevant
                    constraintid = int(row[2])
                    maxsp = float(row[6])
                    limit = float(row[7])
                    contingency = row[4]
                    constraintname = row[3]
                    try:
                        daydata[timestamp]['NetworkConstraints'][constraintid] = [maxsp,limit,contingency,constraintname] #list of binding constraints per interval
                        daydata[timestamp]['NetworkConstraints']['IDList'].add(constraintid)
                    except KeyError:
                        daydata[timestamp]['NetworkConstraints'] = {}
                        daydata[timestamp]['NetworkConstraints'][constraintid] = [maxsp,limit,contingency,constraintname]
                        daydata[timestamp]['NetworkConstraints']['IDList'] = set()
                        daydata[timestamp]['NetworkConstraints']['IDList'].add(constraintid)
                    try:
                        daydata['NetworkConstraints'].append([timestamp,constraintid,maxsp,limit]) #masterlist of binding constraints
                    except KeyError:
                        daydata['NetworkConstraints'] = [[timestamp,constraintid,maxsp,limit]]
    return daydata

def get_shift_factors(wdname,daydata): # get cdr12354 data
    os.chdir(wdname)
    for constraint in daydata['NetworkConstraints']: #initialize shift factor dictionaries
        timestamp = constraint[0]
        constraintid = constraint[1]
        try:
            daydata[timestamp]['ShiftFactors'][constraintid] = {}
        except KeyError:
            daydata[timestamp]['ShiftFactors'] = {}
            daydata[timestamp]['ShiftFactors'][constraintid] = {}
    filenames = glob.glob('cdr.00012354.*.csv')
    for filename in filenames:
        with open(filename, 'r') as csvfile: #open the binding constraint report
            filereader = csv.reader(csvfile)
            for row in filereader:
                if row[0] != 'SCED_TIMESTAMP':
                    try:
                        timestamp = (row[0],row[1])
                        constraintid = int(row[2])
                        resourcename = row[5]
                        shiftfactor = row[6]
                        try:
                            daydata[timestamp]['ShiftFactors'][constraintid][resourcename] = float(shiftfactor)
                        except KeyError:
                            pass
                        if re.search('_CC\d+_\d+',resourcename): #If a combined cycle, also save a shift factor under train name in case configuration mismatch in data
                            cctrainname = re.sub('_\d+','',resourcename)  #strip off _CC# from end of resource name
                            try:
                                daydata[timestamp]['ShiftFactors'][constraintid][cctrainname] = float(shiftfactor) #save SF under CC train name
                                #print(timestamp,constraintid,cctrainname,shiftfactor)
                            except KeyError:
                                pass
                    except IndexError:
                        pass
    return daydata

def true_up_shift_factors(daydata,timestamp,timestampdata):
    dummyvariable = 0.0
    try: 
        timestampdata['ConstraintIDs'] = list(daydata[timestamp]['NetworkConstraints']['IDList'])
    except KeyError:
        timestampdata['ConstraintIDs'] = []
    constraintids = timestampdata['ConstraintIDs']
    print(constraintids)
    for constraintid in constraintids:
        for resource in timestampdata['Resources']:
            try:  
                dummyvariable = daydata[timestamp]['ShiftFactors'][constraintid][resource]
            except KeyError:  
                try: #check to see if this is the CC configuration mismatch issue
                    cctrainname = re.sub('_\d+','',resource) #strip off combined cycle configuration # from end of resource name
                    daydata[timestamp]['ShiftFactors'][constraintid][resource] = daydata[timestamp]['ShiftFactors'][constraintid][cctrainname]#look for SF under CC train name
                    print(cctrainname)
                except KeyError: #if not use 0
                    daydata[timestamp]['ShiftFactors'][constraintid][resource] = 0.0
    return daydata

def get_ramp_rate(wdname,daydata): # get the ramp rate data
    os.chdir(wdname)
    for filename in glob.glob('*SysCon.csv'):
        with open(filename, 'r') as csvfile: #open the binding constraint report
            filereader = csv.reader(csvfile)
            for row in filereader:
                if row[0] != 'SCED_TIMESTAMP':
                    timestamp = (row[0],row[1])
                    resource = row[2]
                    try:
                        daydata[timestamp]['RampRates'][resource] = float(row[3]) * 5.0
                    except KeyError:
                        try:
                            daydata[timestamp]['RampRates'] = {}
                            daydata[timestamp]['RampRates'][resource] = float(row[3]) * 5.0
                        except KeyError:
                            daydata[timestamp] = {}
                            daydata[timestamp]['RampRates'] = {}
                            daydata[timestamp]['RampRates'][resource] = float(row[3]) * 5.0         
    return daydata


def get_nfrc_data(wdname,daydata): # get the ramp rate data
    os.chdir(wdname)
    for filename in glob.glob('nfrc*.csv'):
        with open(filename, 'r') as csvfile: #open the binding constraint report
            filereader = csv.reader(csvfile)
            for row in filereader:
                if row[0] != 'SCED_TIMESTAMP':
                    timestamp = (row[0],row[1])
                    resource = row[2]
                    if daydata[timestamp][resource]['OriginalRRSAward']:
                        daydata[timestamp][resource]['NFRC'] = float(row[3])
                    """try:
                        daydata[timestamp]['NFRC'][resource] = max(float(row[3]), 0.0)
                    except KeyError:
                        try:
                            daydata[timestamp]['NFRC'] = {}
                            daydata[timestamp]['NFRC'][resource] = max(float(row[3]), 0.0)
                        except KeyError:
                            daydata[timestamp] = {}
                            daydata[timestamp]['NFRC'] = {}
                            daydata[timestamp]['NFRC'][resource] = max(float(row[3]), 0.0)"""
    return daydata

def get_powerbalance(permdata): # load undergen and overgen data
    permdata['UnderGen'] = [(5.0,250.0),(5.0,300.0),(10.0,400.0),(10.0,500.0),(10.0,1000.0),(10.0,2250.0),
                            (50.0,4500.0),(50.0,6000.0),(50.0,7500.0),(100000.0,9001.0)]
    permdata['OverGen'] = [(100000.0,250.0)]
    return permdata

def get_sced_segment_info(timestamp,daydata): #gather the SCED data for the timestamp
    pr = range(75,145,2) # iterable of the 35 different offer segments
    mw = range(74,144,2) # iterable of the 35 different offer segments
    badstatuslist = ['OUT','EMR','OFF','OFFNS','Telemetered Resource Status'] #shouldn't be included in optimization
    timestampdata = {}
    timestampdata['Resources'] = set()
    timestampdata['LDLSum'] = 0.0
    counter = 0
    for row in daydata['60DayGRD']:
        status = row[151]
        resource = row[2]
        if (row[0],row[1]) == timestamp and status not in badstatuslist: #check to see if the row contains data relevant to this timestamp
            if status == 'ONTEST':
                output = float(row[153])
                hdl = output
                ldl = output
            elif status == 'ONRR':
                hdl = 0.0
                ldl = 0.0
            else:
                hdl = float(row[147])
                ldl = float(row[150])
            if hdl < 0.05:
                hdl = 0
            timestampdata['Resources'].add(resource)
            timestampdata['LDLSum'] += ldl
            timestampdata[resource] = ldl
            for point in range (1,35,1): #for each segment, grab prold, prnew, mwold, mwnew (defining points of segment)
                mwnew = float(row[mw[point]])
                mwold = float(row[mw[point - 1]])
                prnew = float(row[pr[point]])
                prold = float(row[pr[point - 1]])
                if (mwnew > mwold and mwnew > ldl and mwold < hdl) : #check to se if this is a segment that is within HDL/LDL and should be considered/optimized
                    if mwold < ldl: #check to see if need to adjust MWOLD and PROLD for LDL
                        prold = prold + (ldl - mwold) * (prnew - prold) / (mwnew - mwold) #linearly interpolated price
                        mwold = ldl
                    if mwnew > hdl: #check to see if need to adjust mwnew and prnew for hdl
                        prnew = prnew + (hdl - mwnew) * (prnew - prold) / (mwnew - mwold) #linearly interpolated price
                        mwnew = hdl
                    if (mwnew > mwold and mwnew > ldl and mwold < hdl): #check again to see if segment should still be optimized
                                                                        #If so, load up segment dictionary with data
                        timestampdata[counter] = {} #initialize segment dictionarey
                        timestampdata[counter]['MWOLD'] = mwold
                        timestampdata[counter]['MWNEW'] = mwnew
                        timestampdata[counter]['PROLD'] = prold
                        timestampdata[counter]['PRNEW'] = prnew
                        timestampdata[counter]['ResourceName'] = resource
                        counter += 1
    timestampdata['Counter'] = counter
    print(counter)
    daydata['Log'].append(counter)
    timestampdata['Resources'] = list(timestampdata['Resources'])
    return timestampdata

def create_sced_matrices(timestamp,timestampdata,daydata,permdata): # create the matrices for the sced re run
    """ P is a (numvar + numconstraints + undergens + overgens) x (numvar + numconstraints + undergens + overgens) matrix with diagonal elements
          = price delta/mw delta for the corresponding segment for the first numvar rows/columns (0's for all other entries)
        q is a (numvar + numconstraints + undergens + overgens) x 1 matrix =
          for the first numvar rows represent the low level price for each segment (prold),
          the next numconstraint rows represent the violation amounts for each of the modeled network constraints  - the value for each is the Max Shadow Price of the constraint
          the next undergens rows are the steps of the undergeneration power balance penalty curve (up)
          the last overgens rows are the overgeneration powerbalance penalty curve
        G is a (numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens - 1) + numconstraints rows by numvar + numconstraints +  undergens + overgens columns matrix
          first numvar + numconstraints + undergens + overgens rows are a negative diagonal identity matrix representing low limits for each optimized segment and violation MWs and power balance segments
          second numvar + undergens - 1 + overgens - 1 rows are a positive diagonal identity matrix representing high limits (MW high - MW low for segments, step size for first 8 undergeneration power balance steps) plus a rectangular 0 matrix on right side
          last numconstraints rows are shift factors (segment to constraint) - one row for each modeled constraint
        h is a (numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens - 1) + numconstraints by one column matrix
          first (numvar + numconstraints + undergens + overgens) rows are 0
          second numvar rows are high limits (MW high - MW low) for each segment
          next undergens - 1 + overgens - 1 rows are undergeneration power balance step sizes
          last numconstraints are the mathlimit - LDL flows for each constraint
        A is a 1 row by numvar + numconstraints + undergens + overgens columns matrix.
          First numvar columns are 1.0,
          next numconstraints are 0.0,
          next undergens are 1.0,
          next overgens are -1.0
        b is a 1 by 1 matrix with a single entry of GTBD - sum of LDL energy"""
    startmatrix = time.time()
    timestampdata['UnderGens'] = len(permdata['UnderGen'])
    timestampdata['OverGens'] = len(permdata['OverGen'])
    timestampdata['NumVar'] = timestampdata['Counter']
    try:
        timestampdata['NumConstraints'] = len(daydata[timestamp]['NetworkConstraints']) - 1
    except KeyError:
        timestampdata['NumConstraints'] = 0
    try:
        timestampdata['ConstraintIDs'] = list(daydata[timestamp]['NetworkConstraints']['IDList'])
    except KeyError:
        timestampdata['ConstraintIDs'] = []
    numvar = timestampdata['NumVar']
    undergens = timestampdata['UnderGens']
    overgens = timestampdata['OverGens']
    numconstraints = timestampdata['NumConstraints']
    constraintids = timestampdata['ConstraintIDs']



########################################
######### Create P Matrix ##############
########################################

    timestampdata['PLists'] = []
    for counter in range(numvar): #Create P matrix for segment variable columns
        segment = timestampdata[counter]
        plist = [0.0] * (numvar + numconstraints + undergens + overgens)
        plist[counter] = (segment['PRNEW'] - segment['PROLD']) / (segment['MWNEW'] - segment['MWOLD'])
        timestampdata['PLists'].append(plist)
    for counter in range(numconstraints + undergens + overgens): #Create P matrix for network constraint, undergen and overgen columns
        plist = [0.0] * (numvar + numconstraints + undergens + overgens)
        timestampdata['PLists'].append(plist)
    timestampdata['P'] =   matrix(timestampdata['PLists'])

########################################
######### Create q Matrix ##############
########################################


    timestampdata['qList'] = (numvar + numconstraints + undergens + overgens) * [0.0] #Initialize q vector
    for counter in range(numvar): #linear cost for segements
        segment = timestampdata[counter]
        timestampdata['qList'][counter] = segment['PROLD']
    for constraint in range(numconstraints): #linear costs for network constraint violations (max shadow prices)
        try: constraintid = constraintids[constraint]
        except IndexError: 
            print(timestamp,constraint,numconstraints,constraintids)
            daydata['Log'].append((timestamp,constraint,numconstraints,constraintids))
        timestampdata['qList'][numvar + constraint] = daydata[timestamp]['NetworkConstraints'][constraintid][0] 
    for undergen in range(undergens): #undergen penalty factors
        timestampdata['qList'][numvar + numconstraints + undergen] = permdata['UnderGen'][undergen][1]
    for overgen in range(overgens): #overgen penalty factors
        timestampdata['qList'][numvar + numconstraints + undergens + overgen] = permdata['OverGen'][overgen][1]
    timestampdata['q'] =   matrix(timestampdata['qList'])

########################################
######### Create G Matrix ##############
########################################

    timestampdata['GLists'] = [] #initialize list of lists
    for counter in range(numvar): #Create G matrix for segment variables
        segment = timestampdata[counter]
        resource = segment['ResourceName']
        glist = [0.0] * ((numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens - 1) + numconstraints)
        glist[counter] = -1.0 #each segment energy award >= 0.0
        glist[counter + numvar + numconstraints + undergens + overgens] = 1.0 #each segment energy award <= segment length
        for constraint in range(numconstraints):
            try:
                constraintid = constraint + 1 #need to increase constraint number by one because constraint ids start at 1
                glist[constraint + numvar + numconstraints + undergens + overgens + numvar + undergens - 1 + overgens - 1 ] = \
                  daydata[timestamp]['ShiftFactors'][constraintid][resource]
            except KeyError:
                pass #not all resourcces have a shift factor to a given constraint.  If no shift factor exists in the data, use 0.0
        timestampdata['GLists'].append(glist)
    for constraint in range(numconstraints): #Create G Matrix for network limit violation variables
        glist = [0.0] * ((numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens -1) + numconstraints)
        glist[numvar + constraint] = -1.0
        glist[constraint + numvar + numconstraints + undergens + overgens + numvar + undergens - 1 + overgens - 1] = -1.0
        timestampdata['GLists'].append(glist)
    for undergen in range(undergens): #Create G matrix for undergen variable
        glist = [0.0] * ((numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens -1) + numconstraints)
        glist[numvar + numconstraints + undergen] = -1.0 #undergen >= 0
        if undergen <= undergens - 2:
            glist[2 * numvar + numconstraints + undergens + overgens + undergen] = 1.0 #first 8 undergen segments <= max
        timestampdata['GLists'].append(glist)
    for overgen in range(overgens): #create G matrix for overgen
        glist = [0.0] * ((numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens -1) + numconstraints)
        glist[numvar + numconstraints + undergens + overgen] = -1.0 #overgen >= 0
        if overgen <= overgens - 2:
            glist[numvar + numconstraints + undergens + overgens + undergens - 1 + overgen] = 1.0 #first 8 undergen segments <= max
        timestampdata['GLists'].append(glist)
    timestampdata['G'] =   matrix(timestampdata['GLists'])

########################################
######### Create h Matrix ##############
########################################

    timestampdata['hList'] = [0.0] * (numvar + numconstraints + undergens + overgens) #set inequality limits for >=0 stuff - segment variables, network limit violations, under/overgen
    for counter in range(numvar):
        segment = timestampdata[counter]
        timestampdata['hList'].append(segment['MWNEW'] - segment['MWOLD'])
    for undergen in range(undergens - 1):
        timestampdata['hList'].append(permdata['UnderGen'][undergen][0])
    for overgen in range(overgens - 1):
        timestampdata['hList'].append(permdata['OverGen'][overgen][0])
    for constraint in range(numconstraints):
        constraintid = constraintids[constraint]
        timestampdata['hList'].append(daydata[timestamp]['NetworkConstraints'][constraintid][1])
        for resource in timestampdata['Resources']: #remove the LDL flow from the constraint limit
            try:
                sf = daydata[timestamp]['ShiftFactors'][constraintid][resource]
                timestampdata['hList'][numvar + numconstraints + undergens + overgens + numvar + undergens - 1 + overgens - 1 + constraint] -= sf * timestampdata[resource]
            except KeyError:
                pass
    timestampdata['h'] =   matrix(timestampdata['hList'])

########################################
######### Create A Matrix ##############
########################################

    timestampdata['AList'] = numvar * [1.0]
    for constraint in range(numconstraints):
        timestampdata['AList'].append(0.0)
    for undergen in range(undergens):
        timestampdata['AList'].append(1.0)
    for overgen in range(overgens):
        timestampdata['AList'].append(-1.0)
    timestampdata['A'] =   matrix(timestampdata['AList'],(1,numvar + numconstraints + undergens + overgens))

########################################
######### Create b Matrix ##############
########################################

    timestampdata['b'] = matrix(permdata['SystemConditions'][timestamp]['GTBD'] - timestampdata['LDLSum'])
    return timestampdata

def re_run_sced(timestampdata): #re-run sced with data from 12302, 12354, and 60 day sced report
    P = timestampdata['P']
    q = timestampdata['q']
    G = timestampdata['G']
    h = timestampdata['h']
    A = timestampdata['A']
    b = timestampdata['b']
    timestampdata['Solution']=solvers.qp(P, q, G, h, A, b)
    return timestampdata

def consolidate_sced_results(timestamp,timestampdata,daydata,permdata): #Take the optimization results and consolidate into results
    sol = timestampdata['Solution']
    numvar = timestampdata['NumVar']
    undergens = timestampdata['UnderGens']
    overgens = timestampdata['OverGens']
    numconstraints = timestampdata['NumConstraints']
    constraintids = timestampdata['ConstraintIDs']
    optstatus = sol['status']
    if optstatus == 'optimal':
        timestampdata['Prices'] = {}
        timestampdata['Summary'] = []
        timestampdata['PriceCheck'] = []
        try:
            for constraint in range(numconstraints):
                constraintid = constraintids[constraint]
                daydata[timestamp]['NetworkConstraints'][constraintid].append(sol['z'][numvar + numconstraints + undergens + overgens + numvar + undergens - 1 + overgens - 1 + constraint])
        except KeyError:
            pass
        for resource in timestampdata['Resources']:
            timestampdata['Prices'][(timestamp,resource)] = -1*sol['y'][0]
            try:
                for constraint in range(numconstraints):
                    constraintid = constraintids[constraint]
                    try:
                        timestampdata['Prices'][(timestamp,resource)] -= daydata[timestamp]['ShiftFactors'][constraintid][resource] * \
    		       sol['z'][numvar + numconstraints + undergens + overgens + numvar + undergens - 1 + overgens - 1 + constraint]
                    except KeyError:
                        pass
            except KeyError:
                pass
        klugecost = {}
        productioncost = 0.0
        for resource in timestampdata['Resources']:
            klugecost[(timestamp,resource)] = 0.0
        for counter in range(numvar):
            resource = timestampdata[counter]['ResourceName']
            x = sol['x'][counter]
            P = timestampdata['P'][counter,counter]
            q = timestampdata['q'][counter]
            timestampdata[resource] += x
            klugecost[(timestamp,resource)] += 0.5 * P * x * x + q * x
            productioncost += 0.5 * P * x * x + q * x


        for row in daydata['60DayGRD']:
            grdtimestamp = (row[0],row[1])
            grdresource = row[2]
            if grdtimestamp == timestamp:
                try:
                    row.append(klugecost[(grdtimestamp,grdresource)])
                    row.append(timestampdata[grdresource])
                    row.append(timestampdata['Prices'][(grdtimestamp,grdresource)])
                except KeyError:
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
        timestampdata['Summary'].append(timestamp[0])
        timestampdata['PriceCheck'].append(timestamp[0])
        timestampdata['Summary'].append(timestamp[1])
        timestampdata['PriceCheck'].append(timestamp[1])
        timestampdata['Summary'].append(optstatus)
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['GTBD'])
        timestampdata['PriceCheck'].append(permdata['SystemConditions'][timestamp]['GTBD'])
        timestampdata['Summary'].append(timestampdata['LDLSum'])

        try: timestampdata['Summary'].append(-1*sol['y'][0])
        except KeyError: timestampdata['Summary'].append('NA')
        try: timestampdata['Summary'].append(sol['primal objective'])
        except KeyError: timestampdata['Summary'].append('NA')
        try: timestampdata['Summary'].append(productioncost)
        except KeyError: timestampdata['Summary'].append('NA')
        try: timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['GTBD']  * -1 * sol['y'][0])
        except KeyError: timestampdata['Summary'].append('NA')

        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['RegUpMW'])
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['RRSMW'])
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['NSRSMW'])
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['RegUpPrice'])
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['RRSPrice'])
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['NSRSPrice'])
    else:
        timestampdata['Prices'] = {}
        timestampdata['Summary'] = []
        timestampdata['PriceCheck'] = []
        for row in daydata['60DayGRD']:
            grdtimestamp = (row[0],row[1])
            grdresource = row[2]
            if grdtimestamp == timestamp:
                row.append('NotOptimal')
                row.append('NotOptimal')
                row.append('NotOptimal')
        try:
            for constraint in range(numconstraints):
                constraintid = constraintids[constraint]
        #try:
            #for constraintid in daydata[timestamp]['NetworkConstraints']:
                daydata[timestamp]['NetworkConstraints'][constraintid].append('Not Optimized')
        except KeyError:
            pass
        timestampdata['Summary'].append(timestamp[0])
        timestampdata['Summary'].append(timestamp[1])
        timestampdata['Summary'].append(optstatus)
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['GTBD'])
        timestampdata['Summary'].append(timestampdata['LDLSum'])
        timestampdata['Summary'].append('NA')
        timestampdata['Summary'].append('NA')
        timestampdata['Summary'].append('NA')
        timestampdata['Summary'].append('NA')
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['RegUpMW'])
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['RRSMW'])
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['NSRSMW'])
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['RegUpPrice'])
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['RRSPrice'])
        timestampdata['Summary'].append(permdata['SystemConditions'][timestamp]['NSRSPrice'])
    daydata['Summary'].append(timestampdata['Summary'])
    daydata['PriceCheck'].append(timestampdata['PriceCheck'])
    try: 
        daydata['WelfareToBeat'][timestamp] = sol['primal objective']
    except KeyError:
        daydata['WelfareToBeat'] = {}
        daydata['WelfareToBeat'][timestamp] = sol['primal objective']     
    return daydata

def get_mip_cooptimization_segment_info(timestamp,daydata,permdata,optimizationround): #get the data for real time cooptimization
    timestampdata = {}
    timestampdata['UnderGens'] = len(permdata['UnderGen'])
    timestampdata['OverGens'] = len(permdata['OverGen'])
    timestampdata['LDLSum'] = 0.0
    try:
        timestampdata['NumConstraints'] = len(daydata[timestamp]['NetworkConstraints']) - 1
    except KeyError:
        timestampdata['NumConstraints'] = 0
    try:
        timestampdata['ConstraintIDs'] = list(daydata[timestamp]['NetworkConstraints']['IDList'])
    except KeyError:
        timestampdata['ConstraintIDs'] = []
    counter = 0
    timestampdata['Resources'] = set()
    timestampdata['LDLResources'] = set()
    badstatuslist = ['OUT','EMR','OFF','OFFNS','Telemetered Resource Status'] #shouldn't be included in optimization
    pr = range(75,145,2) # iterable of the 35 different offer segments
    mw = range(74,144,2) # iterable of the 35 different offer segments
    for row in daydata['60DayGRD']: #get a row of data
        status = row[151]
        if status not in badstatuslist and (row[0],row[1]) == timestamp: #only look at stuff that should be optimized inthis go round
            resource = row[2]
############ Section that gets HDL, LDL, BP data for the row
            if status == 'ONTEST':
                output = float(row[153])
                hsl = output
                lsl = output
                ldl = output
                hasl = output
                lasl = output
                hdl = output
            elif status in ('STARTUP','SHUTDOWN','NA'):
                output = float(row[153])
                basepoint = float(row[152])
                hsl = basepoint
                lsl = basepoint
                ldl = basepoint
                hasl = basepoint
                lasl = basepoint
                hdl = basepoint
            elif status == 'ONRR':
                hdl = 0.0
                hsl = float(row[145])
                ldl = 0.0
                hasl = 0.0
                lasl = 0.0
                lsl = 0.0
                output = float(row[153])                
            else:
                hdl = float(row[147])
                hsl = max(float(row[145]),hdl)
                ldl = min(float(row[150]),hdl)
                hasl = max(float(row[146]),hdl)
                lasl = min(float(row[149]),hdl)
                lsl = min(float(row[148]),hdl)
                output = float(row[153])

            if optimizationround == 3:
                try:
                    rrsaward = float(row[156])
                    regaward = float(row[154])
                    try:
                        regdown = float(row[155])
                    except ValueError:
                        regdown = 0.0
                except IndexError:
                    rrsaward = 0.0
                    regaward = 0.0
                    regdown = 0.0


############## Correct HDL for situations that require correction
            if hdl < 0.05: hdl = 0 #accounting for observed SCED behavior
            if (hdl == hasl or status == 'ONREG') and status != 'ONRR' and status != 'ONTEST': #adjust HDL if HASL limited or ramp sharing
                try: 
                    hdl = min(daydata[timestamp]['RampRates'][resource], hsl)
                except KeyError:
                    try:
                        hdl = min(max(output + daydata['ResourceParameters'][resource]['RampRateUp'], hdl), hsl)
                        #print('Missing ramp rate up for:', resource, timestamp)
                        daydata['Log'].append(('Missing ramp rate up for:', resource, timestamp))
                    except KeyError:
                        hdl = min(hdl, hsl)
                        #print('missing permanent ramp rate data for ',resource)
                        daydata['Log'].append(('missing permanent ramp rate data for ',resource))
            timestampdata['LDLSum'] += ldl
            timestampdata['LDLResources'].add(resource)
            timestampdata[resource] = {}
            timestampdata[resource]['LDL'] = ldl
            timestampdata[resource]['HDL'] = hdl
            timestampdata[resource]['LSL'] = lsl
            timestampdata[resource]['HSL'] = hsl
            timestampdata[resource]['Output'] = output
            timestampdata[resource]['Status'] = status
            timestampdata[resource]['OrigNSRSRRSRegAward'] = float(row[157]) + float(row[156]) + float(row[154])
            timestampdata[resource]['OrigRRSRegAward'] = float(row[156]) + float(row[154])
            timestampdata[resource]['OrigRRSAward'] = float(row[156])
            timestampdata[resource]['OrigNSRSAward'] = float(row[157])
            timestampdata[resource]['OrigRegAward'] = float(row[154])
            if optimizationround == 3:
                timestampdata[resource]['PreviousAward'] = {}
                timestampdata[resource]['PreviousAward']['RRS'] = rrsaward
                timestampdata[resource]['PreviousAward']['RegUp'] = regaward
                timestampdata[resource]['PreviousAward']['RegDown'] = regdown

            for point in range (1,35,1): #for each segment, grab PROLD, PRNEW, MWOLD, MWNEW (defining points of segment)
                mwnew = float(row[mw[point]])
                mwold = float(row[mw[point - 1]])
                prnew = float(row[pr[point]])
                prold = float(row[pr[point - 1]])
                if (mwnew > mwold and mwnew > ldl and mwold < hsl) : #check to see if this is a segment that is within HSL/LDL and should be considered/optimized
                  #timestampdata[counter] = {} #initialize segment dictionary
                  if mwold < ldl: #check to see if need to adjust mwold and prold for LDL
                      prold = prold + (ldl - mwold) * (prnew - prold) / (mwnew - mwold) #linearly interpolated price
                      mwold = ldl
                  if mwnew > hsl: #check to see if need to adjust mwnew and prnew for HDL
                      prnew = prnew + (hsl - mwnew) * (prnew - prold) / (mwnew - mwold) #linearly interpolated price
                      mwnew = hsl # I think I can take this out and the constraints will handle this, but leaving language in in case
                  if (mwnew > mwold and mwnew > ldl and mwold < hsl) : #check to see if this is a segment that is within HSL/LDL and should be considered/optimized
#### Load up segment dictionary with data
                    if (mwnew >= mwold + 1.0 and prnew != prold):
                      tempmwold = mwold
                      tempmwnew = mwold + 1.0
                      tempprold = prold
                      while tempmwnew < mwnew:
                          tempprnew = tempprold + (prnew - prold)/(mwnew - mwold) * 1.0
                          #print(counter, tempmwold,tempmwnew,tempprold,tempprnew,mwold,mwnew,prold,prnew)
                          timestampdata[counter] = {}
                          timestampdata[counter]['MWOLD'] = tempmwold
                          timestampdata[counter]['MWNEW'] = tempmwnew
                          timestampdata[counter]['PROLD'] = tempprold
                          timestampdata[counter]['PRNEW'] = tempprnew
                          timestampdata[counter]['ResourceName'] = row[2]
                          timestampdata['Resources'].add(row[2])
                          counter += 1                          
                          tempmwold = tempmwnew
                          tempmwnew = tempmwold + 1.0
                          tempprold = tempprnew
                      #print(counter, tempmwold,tempmwnew,tempprold,tempprnew,mwold,mwnew,prold,prnew)  
                      timestampdata[counter] = {}  
                      timestampdata[counter]['MWOLD'] = tempmwold
                      timestampdata[counter]['MWNEW'] = mwnew
                      timestampdata[counter]['PROLD'] = tempprold
                      timestampdata[counter]['PRNEW'] = prnew
                      timestampdata[counter]['ResourceName'] = row[2]
                      timestampdata['Resources'].add(row[2])
                      counter += 1
                    else:
                      #print(counter, 'no subdivision', row[2],mwold,mwnew,prold,prnew)
                      timestampdata[counter] = {}
                      timestampdata[counter]['MWOLD'] = mwold
                      timestampdata[counter]['MWNEW'] = mwnew
                      timestampdata[counter]['PROLD'] = prold
                      timestampdata[counter]['PRNEW'] = prnew
                      timestampdata[counter]['ResourceName'] = row[2]
                      timestampdata['Resources'].add(row[2])
                      counter += 1

    timestampdata['Counter'] = counter
    timestampdata['Resources'] = list(timestampdata['Resources'])
    timestampdata['LDLResources'] = list(timestampdata['LDLResources'])
    return timestampdata

class Coopt(Problem):
    def model(self,timestamp,timestampdata,daydata,permdata,optimizationround,networklimits):
        """

        optimized variables

          eamt - energy - numvar long, per offer segment

          nsamt - nonspin amount - numresources long, per resource
          rrsamt - rrs amount - numresources long, per resource
          regamt - regup amount - numresources long, per resource
          undergen - energy penalty factor
          nsshort - non-spin shortage amount

          rrsshort - rrs-shortage amount
          regupshort - reg up shortage amount
          conviol - constraint violations
          rrscommit - whether a unit gets nfrc subtracted or not
        """   
        




        nfrc = []
        hsl = []
        ldl = []
        nsshortcost = 75.0
        rrsshortcost = 9001.0
        regupshortcost = []
        regupshortlimit = []
        undergens = len(permdata["UnderGen"])
        numvar = timestampdata['Counter']
        numresources = len(timestampdata['Resources']) 
        numconstraints = timestampdata['NumConstraints']
        constraintids = timestampdata['ConstraintIDs']
        numnfrc = 0
        print(numvar,numresources,numconstraints)
        for resource in range(numresources):
            resname = timestampdata['Resources'][resource]
            nfrc.append(daydata[timestamp][resname]['NFRC'])
            hsl.append(timestampdata[resname]['HSL'])
            ldl.append(timestampdata[resname]['LDL'])
        for n in range(undergens):
            regupshortcost.append(permdata["UnderGen"][n][1])
            if n <= undergens - 2: regupshortlimit.append(permdata["UnderGen"][n][0])

        conviolcost = []
        for constraintnum in range(numconstraints):
            constraintid = constraintids[constraintnum]
            conviolcost.append(daydata[timestamp]['NetworkConstraints'][constraintid][0])


        eoc = []
        for segment in range(numvar):
            eoc.append((timestampdata[segment]['PRNEW'] + timestampdata[segment]['PROLD'])/2)

##############################################
##### Declare vriables for optimization ######
##############################################
            
        self.eamt = eamt = VarVector([numvar],"eamt")
        self.nsamt = nsamt = VarVector([numresources],"nsamt")
        self.rrsamt = rrsamt = VarVector([numresources],"rrsamt")
        self.regamt = regamt = VarVector([numresources],"regamt")
        self.undergen = undergen = Var("undergen")
        self.nsshort = nsshort = Var("nsshort")
        self.rrsshort = rrsshort = Var("rrsshort")
        self.regupshort = regupshort = VarVector([undergens],"regupshort")
        self.conviol = conviol = VarVector([numconstraints],"conviol")
        self.rrscommit = rrscommit = VarVector([numresources],"rrscommit",BIN)

###########################################################################
############# Define Objective Function for optimization ##################
###########################################################################

        minimize( \
            sum_(regupshortcost[h] * regupshort[h] for h in range(undergens)) +\
            sum_(eoc[segment] * eamt[segment] for segment in range(numvar)) + \
            nsshortcost * nsshort +\
            rrsshortcost * rrsshort +\
            9001.0 * undergen +\
            sum_(conviolcost[i] * conviol[i] for i in range(numconstraints))\
        )
###########################################################################
############ Define Constraints for optimization ##########################
###########################################################################


        for segment in range(numvar): #segment mw award <= segment length

            eamt[segment] <= timestampdata[segment]['MWNEW'] - timestampdata[segment]['MWOLD']

        for resource in range(numresources): #total resource responsibility <= resource hsl
            nsamt[resource] + rrsamt[resource] + regamt[resource] +\
            sum_(eamt[segment] for segment in range(numvar) if timestampdata[segment]['ResourceName'] == timestampdata['Resources'][resource]) \
            <= hsl[resource] - ldl[resource] - nfrc[resource] * rrscommit[resource]


        for resource in range(numresources): #reg up <= ramp rate
            resname = timestampdata['Resources'][resource]
            if daydata['ResourceParameters'][resname]['RegUpHMWh'] and timestampdata[resname]['Status'] != 'OFFQS':
                regamt <= daydata[timestamp]['RampRates'][resname]
            else:
                regamt <= 0

        rrslimit = {}
        for resource in range(numresources): #set RRS max award
            resname = timestampdata['Resources'][resource]
            if daydata['ResourceParameters'][resname]['RRSHMWh'] and timestampdata[resname]['Status'] != 'OFFQS':
                if timestampdata[resname]['Status'] == 'ONRR':
                    rrslimit[resource] = (max(hsl[resource] - ldl[resource], timestampdata[resname]['OrigRRSAward'])) 
                elif timestampdata[resname]['Status'] not in ['OFFQS','FRRSUP']:
                    rrslimit[resource] = max(0.2 * hsl[resource], timestampdata[resname]['OrigRRSAward'])
                else:
                    rrslimit[resource] = 0.000
            else: rrslimit[resource] = 0.000

        for resource in range(numresources): #RRS award only if rrscommit(binary optimization variable) = 1 also RRS Award <= RRS Max
            rrsamt[resource] <= rrslimit[resource] * rrscommit[resource]


        for resource in range(numresources): #RRS + Reg award <= 2 * ramp
            resname = timestampdata['Resources'][resource]
            rrsamt[resource] + regamt[resource] <= max(2 * daydata[timestamp]['RampRates'][resname], timestampdata[resname]['OrigRRSAward'] +\
                                                                                 timestampdata[resname]['OrigRegAward'])

        for resource in range(numresources): #NSRS + RRS + Reg <= 6*ramp
            resname = timestampdata['Resources'][resource]
            nsamt[resource] + rrsamt[resource] + regamt[resource] <= max(6 * daydata[timestamp]['RampRates'][resname], \
                                          timestampdata[resname]['OrigRegAward'] + timestampdata[resname]['OrigRRSAward'] +\
                                          timestampdata[resname]['OrigNSRSAward'])
                                
        for resource in range(numresources): #5/7 Reg + dispatch up 
            resname = timestampdata['Resources'][resource]
            (5/7) * regamt[resource] + \
            sum_(eamt[segment] for segment in range(numvar) if timestampdata[segment]['ResourceName'] == timestampdata['Resources'][resource]) \
            - timestampdata[resname]['Output'] +  timestampdata[resname]['LDL']\
            <= daydata[timestamp]['RampRates'][resname]

        for segment in range(undergens - 1):
            regupshort[segment] <= regupshortlimit[segment]

        for constraintnum in range(numconstraints):
            constraintid = constraintids[constraintnum]
            sum_(daydata[timestamp]['ShiftFactors'][constraintid][timestampdata[segment]['ResourceName']]*eamt[segment] \
            for segment in range(numvar)) <= networklimits[constraintid] + conviol[constraintnum]


        for resource in range(numresources):
            resname = timestampdata['Resources'][resource]
            if timestampdata[resname]['Status'] == 'ONRR': 
                sum_(eamt[segment] for segment in range(numvar) if timestampdata[segment]['ResourceName'] == timestampdata['Resources'][resource]) +\
                nsamt[resource] + \
                regamt[resource] \
                <= 0.0

        for resource in range(numresources):
            resname = timestampdata['Resources'][resource]
            if daydata['ResourceParameters'][resname]['NSRSHMWh'] < 0.1:
                nsamt[resource] <= 0.0

        sum_(eamt[segment] for segment in range(numvar)) + undergen == \
        permdata['SystemConditions'][timestamp]['GTBD'] - timestampdata['LDLSum']
        sum_(nsamt[resource] for resource in range(numresources)) +  nsshort == permdata['SystemConditions'][timestamp]['NSRSMW']
        sum_(rrsamt[resource] for resource in range(numresources)) + rrsshort == permdata['SystemConditions'][timestamp]['RRSMW']
        sum_(regamt[resource] for resource in range(numresources)) + sum_(regupshort[h] for h in range(undergens)) == permdata['SystemConditions'][timestamp]['RegUpMW']

def mip_coopt(timestamp,timestampdata,daydata,permdata,optimizationround):
    """
    optimized variables
        eamt - energy - numvar long, per offer segment
        nsamt - nonspin amount - numresources long, per resource
        rrsamt - rrs amount - numresources long, per resource
        regamt - regup amount - numresources long, per resource
        undergen - energy penalty factor
        nsshort - non-spin shortage amount
        rrsshort - rrs-shortage amount
        regupshort - reg up shortage amount
        conviol - constraint violations
    """     
    print('NUMVAR ###########',timestampdata['Counter'])  
    numconstraints = timestampdata['NumConstraints']
    constraintids = timestampdata['ConstraintIDs']
    networklimits = {}
    for constraintnum in range(numconstraints):
        constraintid = constraintids[constraintnum]
        currentlimit = daydata[timestamp]['NetworkConstraints'][constraintid][1] #get Limit
        for resource in timestampdata['LDLResources']:
            try:
                sf = daydata[timestamp]['ShiftFactors'][constraintid][resource]
                ldl = timestampdata[resource]['LDL']
                currentlimit -= sf * ldl #subtract out LDL component
            except KeyError:
                pass
        networklimits[constraintid] = currentlimit 
    prob = Coopt("test1")
    prob.model(timestamp,timestampdata,daydata,permdata,optimizationround,networklimits)
    try:
        prob.optimize(False)
    except KeyboardInterrupt:
        prob.is_solution = False
    rrscommit = {}
    numresources = len(timestampdata['Resources'])
    for resnum in range(numresources):
        resource = timestampdata['Resources'][resnum]
        rrscommit[resource] = prob.rrscommit[resnum].val > 0.5
    if prob.is_solution == True:
        tempsol = 'optimal'
    else:
        tempsol = 'not optimal'
    daydata['RRSCommit'][timestamp] = rrscommit
    daydata['Summary'][-1].append(tempsol)
    daydata[timestamp]['MIPStatus'] = tempsol
    return daydata

def get_cooptimization_segment_info(timestamp,daydata,permdata,optimizationround): #get the data for real time cooptimization
    timestampdata = {}
    timestampdata['UnderGens'] = len(permdata['UnderGen'])
    timestampdata['OverGens'] = len(permdata['OverGen'])
    timestampdata['LDLSum'] = 0.0
    try:
        timestampdata['NumConstraints'] = len(daydata[timestamp]['NetworkConstraints']) - 1
    except KeyError:
        timestampdata['NumConstraints'] = 0
    try:
        timestampdata['ConstraintIDs'] = list(daydata[timestamp]['NetworkConstraints']['IDList'])
    except KeyError:
        timestampdata['ConstraintIDs'] = []
    counter = 0
    timestampdata['Resources'] = set()
    timestampdata['LDLResources'] = set()
    badstatuslist = ['OUT','EMR','OFF','OFFNS','Telemetered Resource Status'] #shouldn't be included in optimization
    pr = range(75,145,2) # iterable of the 35 different offer segments
    mw = range(74,144,2) # iterable of the 35 different offer segments
    for row in daydata['60DayGRD']: #get a row of data
        status = row[151]
        if status not in badstatuslist and (row[0],row[1]) == timestamp: #only look at stuff that should be optimized inthis go round
            resource = row[2]
############ Section that gets HDL, LDL, BP data for the row
            if status == 'ONTEST':
                output = float(row[153])
                hsl = output
                lsl = output
                ldl = output
                hasl = output
                lasl = output
                hdl = output
            elif status in ('STARTUP','SHUTDOWN','NA'):
                output = float(row[153])
                basepoint = float(row[152])
                hsl = basepoint
                lsl = basepoint
                ldl = basepoint
                hasl = basepoint
                lasl = basepoint
                hdl = basepoint
            elif status == 'ONRR':
                hdl = 0.0
                hsl = float(row[145])
                ldl = 0.0
                hasl = 0.0
                lasl = 0.0
                lsl = 0.0
                output = float(row[153])                
            else:
                hdl = float(row[147])
                hsl = max(float(row[145]),hdl)
                ldl = min(float(row[150]),hdl)
                hasl = max(float(row[146]),hdl)
                lasl = min(float(row[149]),hdl)
                lsl = min(float(row[148]),hdl)
                output = float(row[153])

            if optimizationround == 3:
                try:
                    rrsaward = float(row[156])
                    regaward = float(row[154])
                    try:
                        regdown = float(row[155])
                    except ValueError:
                        regdown = 0.0
                except IndexError:
                    rrsaward = 0.0
                    regaward = 0.0
                    regdown = 0.0


############## Correct HDL for situations that require correction
            if hdl < 0.05: hdl = 0 #accounting for observed SCED behavior
            if (hdl == hasl or status == 'ONREG') and status != 'ONRR' and status != 'ONTEST': #adjust HDL if HASL limited or ramp sharing
                try: 
                    hdl = min(daydata[timestamp]['RampRates'][resource], hsl)
                except KeyError:
                    try:
                        hdl = min(max(output + daydata['ResourceParameters'][resource]['RampRateUp'], hdl), hsl)

                        #print('Missing ramp rate up for:', resource, timestamp)
                        daydata['Log'].append(('Missing ramp rate up for:', resource, timestamp))
                    except KeyError:
                        hdl = min(hdl, hsl)
                        #print('missing permanent ramp rate data for ',resource)
                        daydata['Log'].append(('missing permanent ramp rate data for ',resource))
            timestampdata['LDLSum'] += ldl
            timestampdata['LDLResources'].add(resource)
            timestampdata[resource] = {}
            timestampdata[resource]['LDL'] = ldl
            timestampdata[resource]['HDL'] = hdl
            timestampdata[resource]['LSL'] = lsl
            timestampdata[resource]['HSL'] = hsl
            timestampdata[resource]['Output'] = output
            timestampdata[resource]['Status'] = status
            timestampdata[resource]['OrigNSRSRRSRegAward'] = float(row[157]) + float(row[156]) + float(row[154])
            timestampdata[resource]['OrigRRSRegAward'] = float(row[156]) + float(row[154])
            timestampdata[resource]['OrigRRSAward'] = float(row[156])
            timestampdata[resource]['OrigNSRSAward'] = float(row[157])
            timestampdata[resource]['OrigRegAward'] = float(row[154])
            if optimizationround == 3:
                timestampdata[resource]['PreviousAward'] = {}
                timestampdata[resource]['PreviousAward']['RRS'] = rrsaward
                timestampdata[resource]['PreviousAward']['RegUp'] = regaward
                timestampdata[resource]['PreviousAward']['RegDown'] = regdown

            for point in range (1,35,1): #for each segment, grab PROLD, PRNEW, MWOLD, MWNEW (defining points of segment)
                mwnew = float(row[mw[point]])
                mwold = float(row[mw[point - 1]])
                prnew = float(row[pr[point]])
                prold = float(row[pr[point - 1]])
                if (mwnew > mwold and mwnew > ldl and mwold < hsl) : #check to see if this is a segment that is within HSL/LDL and should be considered/optimized
                  timestampdata[counter] = {} #initialize segment dictionary
                  if mwold < ldl: #check to see if need to adjust mwold and prold for LDL
                      prold = prold + (ldl - mwold) * (prnew - prold) / (mwnew - mwold) #linearly interpolated price
                      mwold = ldl
                  if mwnew > hsl: #check to see if need to adjust mwnew and prnew for HDL
                      prnew = prnew + (hsl - mwnew) * (prnew - prold) / (mwnew - mwold) #linearly interpolated price
                      mwnew = hsl # I think I can take this out and the constraints will handle this, but leaving language in in case
                  if (mwnew > mwold and mwnew > ldl and mwold < hsl) : #check to see if this is a segment that is within HSL/LDL and should be considered/optimized
#### Load up segment dictionary with data
                      timestampdata[counter]['MWOLD'] = mwold
                      timestampdata[counter]['MWNEW'] = mwnew
                      timestampdata[counter]['PROLD'] = prold
                      timestampdata[counter]['PRNEW'] = prnew
                      timestampdata[counter]['ResourceName'] = row[2]
                      timestampdata['Resources'].add(row[2])
                      counter += 1
    timestampdata['Counter'] = counter
    timestampdata['Resources'] = list(timestampdata['Resources'])
    timestampdata['LDLResources'] = list(timestampdata['LDLResources'])
    return timestampdata

def create_cooptimization_matrices(timestamp,timestampdata,daydata,permdata,pfudge):
    """ P is a (numvar + 1 + numresources + 1 + numresources + 1 + numresources + undergens + numconstraints) ^ 2 with diagonal elements
          = price delta/mw delta for the corresponding segment for the first numvar rows/columns (all 0.000001's other)"""
    timestampdata['Trouble'] = []
    numvar = timestampdata['Counter']
    undergens = timestampdata['UnderGens']
    overgens = timestampdata['OverGens']
    numresources = len(timestampdata['Resources'])
    numconstraints = timestampdata['NumConstraints']
    constraintids = timestampdata['ConstraintIDs']
    timestampdata['PLists'] = []
    rrscommit = daydata['RRSCommit'][timestamp]
    for counter in range(numvar): #Create P matrix for segment variable columns
        segment = timestampdata[counter]
        plist = [0.0] * (numvar + 1 + numresources + 1 + numresources + 1 + numresources + undergens + numconstraints)
        try:
            plist[counter] = max((segment['PRNEW'] - segment['PROLD']) / (segment['MWNEW'] - segment['MWOLD']),pfudge) #making sure all diagonal P >0
        except ZeroDivisionError: #troubleshooting should never happen
            print(TimeStamp, segment['ResourceName'], segment['MWOLD'], segment['MWNEW'])
        timestampdata['PLists'].append(plist)
    for counter in range(1 + numresources + 1 + numresources + 1 + numresources + undergens + numconstraints): #Create P matrix for energy snuff, ns, ns ins, rrs, rrs ins, reg, reg ins (10) and constr viol columns
        plist = [0.0] * (numvar + 1 + numresources + 1 + numresources + 1 + numresources + undergens + numconstraints)
        plist[numvar + counter] = pfudge #making sure all diagonal P >0
        timestampdata['PLists'].append(plist)
    timestampdata['P'] =   matrix(timestampdata['PLists'])

    """q is a (numvar + 1 + numresources + 1 + numresources + 1 + numresources + undergens + numconstraints) x 1 array/vector =
       for the first numvar rows represent the low level price for each segment,
       the next row is the undergeneration penalty
       the next numresources rows are nonspin costs (0)
       the next row is the NonSpin insufficiency penalty ($75)
       the next numresources rows are the rrs costs (0)
       the next row is the rrs insufficiency penalty ($9001)
       the next numresources rows are the reg up costs ($0)
       the next undergen rows are the reg up insufficiency costs
       the next numconstraint rows represent the violation amounts for each of the modeled network constraints  - the value for each is the Max Shadow Price of the constraint
       """

    timestampdata['qList'] = (numvar + 1 + numresources + 1 + numresources + 1 + numresources + undergens + numconstraints) * [0.0] #Initialize q vector
    for counter in range(numvar): #linear cost for segements
        segment = timestampdata[counter]
        timestampdata['qList'][counter] = segment['PROLD']
    timestampdata['qList'][numvar] = 9001 #undergen penalty factor
    timestampdata['qList'][numvar + 1 + numresources] = 75 #nonspin insufficiency penalty factor
    timestampdata['qList'][numvar + 1 + numresources + 1 + numresources] = 9001 #rrs insufficiency penalty factor
    positioncount = 0
    for position in range((numvar + 1 + numresources + 1 + numresources + 1 + numresources),(numvar + 1 + numresources + 1 + numresources + 1 + numresources + undergens)):
        timestampdata['qList'][position] = permdata['UnderGen'][positioncount][1] #regulation up penalty factors (from undergen penalties)
        print("q",positioncount,position,timestampdata['qList'][position])
        timestampdata['Trouble'].append([positioncount,position,timestampdata['qList'][position]])
        positioncount += 1
    constraint = 0
    for constraintnum in range(numconstraints):
        constraintid = constraintids[constraintnum]
        timestampdata['qList'][numvar + 1 + numresources + 1 + numresources + 1 + numresources + undergens + constraintnum] \
                = daydata[timestamp]['NetworkConstraints'][constraintid][0]

    timestampdata['q'] =   matrix(timestampdata['qList'])

    """
    G is a (2 * numvar + 9 * numresources + 2 * numconstraints + 2 * undergens + 2) rows by (numvar + 1 + numresources + 1 + numresources + 1 + numresources + undergens + numconstraints) columns matrix
    first numvar rows are: energy (per bid segment) >= 0 MW  ## 0
    the next row is: energy insufficiency >=0 ## 1* numvar
    the next numresources rows are: NSRS (per resource) >= 0 ## 1 * numvar + 1 
    the next row is: NSRS insufficiency >= 0 1* numvar + 1 + 1 * numresources
    the next numresources rows are: RRS (per resource) >= 0 ## 1 * numvar + 2 + 1 * numresources
    the next row is RRS insufficiency >= 0 ## 1 * numvar + 2 + 2 * numresources
    the next numresources rows are: RegUp >= 0  ## 1 * numvar + 3 + 2 * numresources
    the next undergen rows are: RegUp insufficiency (per power balance penalty curve segments) >= 0 ## 1 * numvar + 3 + 3 * numresources
    the next numconstraints rows are: constraint violation >= 0 ## 1 * numvar + 3 + 3 * numresources + undergens
    the next numvar rows are: energy (per bid segment) <= bid segment MW  ## 1 * numvar + 3 + 3 * numresources + undergens + numconstraints
    the next numresources rows are: Energy + NSRS + RRS + RegUp <= HSL - LDL ## 2 * numvar + 3 + 3 * numresources + undergens + numconstraints
    the next numresources rows are: RegUp <= 5 minute Ramp Rate Up ## 2 * numvar + 3 + 4 * numresources + undergens + numconstraints
    the next numresources rows are: RRS <= 0.2 * HSL ## 2 * numvar + 3 + 5 * numresources + undergens + numconstraints
    the next numresources rows are: RRS + RegUp <= 10 minute ramp rate up) ## 2 * numvar + 3 + 6 * numresources + undergens + numconstraints
    the next numresources rows are: NSRS + RRS + RegUp <= 30 minute ramp rate up ## 2 * numvar + 3 + 7 * numresources + undergens + numconstraints
    the next numresources rows are: Energy + (5/7) * RegUp <= Output + 5 minute ramp rate up - LDL ## 2 * numvar + 3 + 8 * numresources + undergens + numconstraints
    the next (undergens - 1) rows are: Reg Insufficient segment energy <= segment MWs ## 2 * numvar + 3 + 9 * numresources + undergens + numconstraints
    the next numconstraint rows are: Sum of SF * Energy <= mathlimt - Sum of SF * LDLs ## 2 * numvar + 2 + 9 * numresources + 2 * undergens + numconstraints 
    the next numresources rows are: Sum of Energy + RegUp + NSRS <= LDL if status = 'ONRR' ## 2 * numvar + 2 + 9 * numresources + 2 * undergens + 2 * numconstraints
    the next numresources rows are: NSRS <= 0 if haven't sold NSRS before  ## 2 * numvar + 2 + 10 * numresources + 2 * undergens + 2 * numconstraints
    """
    Glists = [] #initialize list of lists
    ###############################################################
    ### create columns for each segment energy awards #############
    ###############################################################
    for column in range(numvar): 
        resource = timestampdata[column]['ResourceName']
        glist = (2 * numvar + 11 * numresources + 2 * numconstraints + 2 * undergens + 2) * [0.0]
        glist[column] = -1.0 #energy (per bid segment) >= 0 MW
        glist[numvar + 3 + 3 * numresources + undergens + numconstraints + column] = 1.0 #energy (per bid segment) <= bid segment MW
        resourcenum = 0
        for rowresource in timestampdata['Resources']:
            if rowresource == resource:
                glist[2 * numvar + 3 + 3 * numresources + undergens + numconstraints + resourcenum] = 1.0 # Energy + NSRS + RRS + RegUp <= HSL - LDL
                glist[2 * numvar + 3 + 8 * numresources + undergens + numconstraints + resourcenum] = 1.0 # Energy + (5/7) * RegUp <= Output + 5 minute ramp rate up - LDL
                glist[2 * numvar + 2 + 9 * numresources + 2 * undergens + 2 * numconstraints + resourcenum] = 1.0 # Energy + NSRS + RegUp <= LDL if 'ONRR'
            resourcenum += 1
        try:
            for constraintnum in range(numconstraints):
                constraintid = constraintids[constraintnum]
                try:
                    glist[2 * numvar + 9 * numresources + numconstraints + 2 * undergens + 2 + constraintnum] = daydata[timestamp]['ShiftFactors'][constraintid][resource]
                except KeyError:
                    pass
                constraintnum += 1
        except KeyError:
            pass
        Glists.append(glist)

    ###############################################################
    ### create column for energy insufficiency ####################
    ###############################################################
    glist = (2 * numvar + 11 * numresources + 2 * numconstraints + 2 * undergens + 2) * [0.0] #energy insufficiency column
    glist[numvar] = -1.0
    Glists.append(glist)

    ###############################################################
    ############# create columns for NSRS #########################
    ###############################################################

    for resourcenum in range(numresources): 
        glist = (2 * numvar + 11 * numresources + 2 * numconstraints + 2 * undergens + 2) * [0.0]
        glist[numvar + 1 + resourcenum] = -1.0 #NSRS (per resource) >= 0 
        glist[2 * numvar + 3 + 3 * numresources + undergens + numconstraints + resourcenum] = 1.0 # Energy + NSRS + RRS + RegUp <= HSL - LDL
        glist[2 * numvar + 3 + 7 * numresources + undergens + numconstraints + resourcenum] = 1.0 # NSRS + RRS + RegUp <= 30 minute ramp rate up
        glist[2 * numvar + 2 + 9 * numresources + 2 * undergens + 2 * numconstraints + resourcenum] = 1.0 # Energy + NSRS + RegUp <= LDL if 'ONRR'
        glist[2 * numvar + 2 + 10 * numresources + 2 * undergens + 2 * numconstraints + resourcenum] = 1.0 # NSRS <= 0 if not sold NSRS before
        Glists.append(glist)

    ###############################################################
    ############# create column for NSRS insufficency #############
    ###############################################################

    glist = (2 * numvar + 11 * numresources + 2 * numconstraints + 2 * undergens + 2) * [0.0] #NSRS insufficiency column
    glist[numvar + 1 + numresources] = -1.0
    Glists.append(glist)

    ###############################################################
    ############# create columns for RRS ##########################
    ###############################################################


    for resourcenum in range(numresources): 
        glist = (2 * numvar + 11 * numresources + 2 * numconstraints + 2 * undergens + 2) * [0.0]
        glist[numvar + numresources + 2 + resourcenum] = -1.0 # RRS  >= 0
        glist[2 * numvar + 3 + 3 * numresources + undergens + numconstraints + resourcenum] = 1.0 # Energy + NSRS + RRS + RegUp <= HSL - LDL
        glist[2 * numvar + 3 + 5 * numresources + undergens + numconstraints + resourcenum] = 1.0 # RRS <= 0.2 * HSL
        glist[2 * numvar + 3 + 6 * numresources + undergens + numconstraints + resourcenum] = 1.0 # RRS + RegUp <= 10 minute ramp rate up)
        glist[2 * numvar + 3 + 7 * numresources + undergens + numconstraints + resourcenum] = 1.0 # NSRS + RRS + RegUp <= 30 minute ramp rate up
        Glists.append(glist)

    ###############################################################
    ############# create column for RRS insufficency ##############
    ###############################################################

    glist = (2 * numvar + 11 * numresources + 2 * numconstraints + 2 * undergens + 2) * [0.0] #RRS insufficiency column
    glist[numvar + 2 + 2 * numresources] = -1.0
    Glists.append(glist)

    ###############################################################
    ############# create columns for RegUp ########################
    ###############################################################

    for resourcenum in range(numresources): 
        glist = (2 * numvar + 11 * numresources + 2 * numconstraints + 2 * undergens + 2) * [0.0]
        glist[numvar + 2 * numresources + 3 + resourcenum] = -1.0 # RegUp  >= 0
        glist[2 * numvar + 3 + 3 * numresources + undergens + numconstraints + resourcenum] = 1.0 # Energy + NSRS + RRS + RegUp <= HSL - LDL
        glist[2 * numvar + 3 + 4 * numresources + undergens + numconstraints + resourcenum] = 1.0 # RegUp <= 5 minute Ramp Rate Up
        glist[2 * numvar + 3 + 6 * numresources + undergens + numconstraints + resourcenum] = 1.0 # RRS + RegUp <= 10 minute ramp rate up)
        glist[2 * numvar + 3 + 7 * numresources + undergens + numconstraints + resourcenum] = 1.0 # NSRS + RRS + RegUp <= 30 minute ramp rate up
        glist[2 * numvar + 3 + 8 * numresources + undergens + numconstraints + resourcenum] = 5.0/7.0 # NSRS + RRS + RegUp <= 30 minute ramp rate up
        glist[2 * numvar + 2 + 9 * numresources + 2 * undergens + 2 * numconstraints + resourcenum] = 1.0 # Energy + NSRS + RegUp <= LDL if 'ONRR'
        Glists.append(glist)

    ###############################################################
    ########## create columns for RegUp insufficency ##############
    ###############################################################

    for column in range(undergens):
        glist = (2 * numvar + 11 * numresources + 2 * numconstraints + 2 * undergens + 2) * [0.0]
        glist[numvar + 3 * numresources + 3 + column] = -1.0 # RegUp insufficiency >= 0
        if column <= (undergens - 2):
            glist[2 * numvar + 3 + 9 * numresources + undergens + numconstraints + column] = 1.0 # powerbalance segment energy <= segment MWs
        Glists.append(glist)
        print("G",column,2 * numvar + 3 + 9 * numresources + undergens + numconstraints + column)

    ###############################################################
    ######## create columns for Constraint Violations #############
    ###############################################################

    for column in range(numconstraints):
        glist = (2 * numvar + 11 * numresources + 2 * numconstraints + 2 * undergens + 2) * [0.0]
        glist[numvar + 3 + 3 * numresources + undergens + column] = -1.0 # constraint violation >= 0
        glist[2 * numvar + 2 + 9 * numresources +  2 * undergens + numconstraints + column] = -1.0 # Sum of SF * Energy <= mathlimt - Sum of SF * LDLs 
        Glists.append(glist)

    timestampdata['G'] =   matrix(Glists)
    Glists = []
    print('finishing G')

    """
    h is a (2 * numvar + 9 * numresources + 2 * numconstraints + 2 * undergens + 2) rows by 1 column matrix/vector
    first numvar rows are: energy (per bid segment) >= 0 MW
    the next row is: energy insufficiency >=0
    the next numresources rows are: NSRS (per resource) >= 0 
    the next row is: NSRS insufficiency >= 0
    the next numresources rows are: RRS (per resource) >= 0
    the next row is RRS insufficiency >= 0
    the next numresources rows are: RegUp >= 0
    the next undergen rows are: RegUp insufficiency (per power balance penalty curve segments) >= 0
    the next numconstraints rows are: constraint violation >= 0
    the next numvar rows are: energy (per bid segment) <= bid segment MW
    the next numresources rows are: Energy + NSRS + RRS + RegUp <= HSL - LDL
    the next numresources rows are: RegUp <= 5 minute Ramp Rate Up
    the next numresources rows are: RRS <= 0.2 * HSL
    the next numresources rows are: RRS + RegUp <= 10 minute ramp rate up)
    the next numresources rows are: NSRS + RRS + RegUp <= 30 minute ramp rate up
    the next numresources rows are: Energy + (5/7) * RegUp <= Output + 5 minute ramp rate up - LDL
    the next (undergens - 1) rows are: Reg Insufficient segment energy <= segment MWs
    the next numconstraint rows are: Sum of SF * Energy <= mathlimt - Sum of SF * LDLs 
    the next numresources rows are: Sum of Energy + RegUp + NSRS <= LDL if status = 'ONRR'
    the next numresources rows are: NSRS <= 0 if haven't sold NSRS before
    """
    hlist = (2 * numvar + 11 * numresources + 2 * numconstraints + 2 * undergens + 2) * [0.0]
    for counter in range(numvar):
        hlist[numvar + 3 + 3 * numresources + undergens + numconstraints + counter] = timestampdata[counter]['MWNEW'] - timestampdata[counter]['MWOLD'] # Energy <= segment MW
    rescounter = 0
    for resource in timestampdata['Resources']:
        origaward = timestampdata[resource]['OrigRRSRegAward']
        orignsrsaward = timestampdata[resource]['OrigNSRSAward']
        origregaward = timestampdata[resource]['OrigRegAward']
        origrrsaward = timestampdata[resource]['OrigRRSAward']
        orignsrsrrsregupaward = timestampdata[resource]['OrigNSRSRRSRegAward']
        if daydata[timestamp]['MIPStatus'] == 'not optimal':
            rrscommit[resource] = origrrsaward >= 0.5 #if MIP commit optimization not optimal, just go with actual NFRC dispatch
        hdl = timestampdata[resource]['HDL']
        output = timestampdata[resource]['Output']
        status = timestampdata[resource]['Status'] 
        hasnfrc = daydata['ResourceParameters'][resource]['NFRC'] >= 0.05
        nfrc = daydata[timestamp][resource]['NFRC']
        try:
            ramp = daydata[timestamp]['RampRates'][resource]
        except KeyError:
            ramp = max(timestampdata[resource]['HDL'] - timestampdata[resource]['Output'], 0)
            print('missing ramp rate in G data ',resource,timestamp) #log
        try:
            regupmwh = daydata['ResourceParameters'][resource]['RegUpHMWh'] 
        except KeyError:
            regupmwh = 0.0
            print('missing regupmwh from daydata', resource) #log
        try:
            rrsmwh = daydata['ResourceParameters'][resource]['RRSHMWh'] 
        except KeyError:
            rrsmwh = 0.0
            print('missing rrsmwh from daydata', resource) #log

        try:
            nsrsmwh = daydata['ResourceParameters'][resource]['NSRSHMWh'] 
        except KeyError:
            nsrsmwh = 0.0
            print('missing nsrsmwh from permdata', resource) #log
        try:
            freqaward = timestampdata[resource]['PreviousAward']['RRS'] >= 0.1 
        except KeyError:
            freqaward = False
        if timestampdata[resource]['Status'] == 'ONTEST': 
            hlist[2 * numvar + 3 * numresources + 3 + undergens + numconstraints + rescounter] = 0.0 #no dispatchable headroom for ONTEST units
            timestampdata[resource]['NFRCUsed'] = 'NA'
        elif rrscommit[resource]:
            hlist[2 * numvar + 3 * numresources + 3 + undergens + numconstraints + rescounter] = \
              max(timestampdata[resource]['HSL'] - timestampdata[resource]['LDL'] - nfrc, 0.0) # Energy + RRS + RegUp <= HSL - LDL - NFRC if MIP committed RRS for resource
            timestampdata[resource]['NFRCUsed'] = nfrc
        else:
            hlist[2 * numvar + 3 * numresources + 3 + undergens + numconstraints + rescounter] = \
              max(timestampdata[resource]['HSL'] - timestampdata[resource]['LDL'], 0) # Energy + NSRS + RRS + RegUp <= HSL - LDL
            timestampdata[resource]['NFRCUsed'] = 0.0
        timestampdata[resource]['HSL-LDLUsed'] = hlist[2 * numvar + 3 * numresources + 3 + undergens + numconstraints + rescounter]

        if (regupmwh or origregaward) and timestampdata[resource]['Status']  not in  ['ONRR', 'OFFQS']:
            hlist[2 * numvar + 4 * numresources + 3 + undergens + numconstraints + rescounter] = ramp # RegUp <= 5 minute Ramp Rate Up
        timestampdata[resource]['RegUpOffCapUsed'] = hlist[2 * numvar + 4 * numresources + 3 + undergens + numconstraints + rescounter]

        if ((rrsmwh or origrrsaward) and (not(hasnfrc) or rrscommit[resource])  ):  #only award RRS if resource is qualified AND it either got an RRS/NFRC commit from MIP or has no NFRC -commit doesn't matter
            if timestampdata[resource]['Status'] == 'ONRR':
                hlist[2 * numvar + 5 * numresources + 3 + undergens + numconstraints + rescounter] = max(timestampdata[resource]['HSL'] , timestampdata[resource]['OrigRRSAward'])
            elif timestampdata[resource]['Status'] not in ['OFFQS','FRRSUP']:
                hlist[2 * numvar + 5 * numresources + 3 + undergens + numconstraints + rescounter] = max(0.2 * (timestampdata[resource]['HSL'] ), timestampdata[resource]['OrigRRSAward'])
            # RRS <= 0.2 * HSL unless OFFQS or FFRSUP or (optround 2 and the resource has nfrc) in which case RRS <= 0
        timestampdata[resource]['RRSOffCapUsed'] = hlist[2 * numvar + 5 * numresources + 3 + undergens + numconstraints + rescounter]

        hlist[2 * numvar + 6 * numresources + 3 + undergens + numconstraints + rescounter] = max(2 * ramp, origaward) # RRS + RegUp <= 10 minute ramp rate up
        timestampdata[resource]['RRSRegEnergyCapUsed'] = hlist[2 * numvar + 6 * numresources + 3 + undergens + numconstraints + rescounter]

        hlist[2 * numvar + 7 * numresources + 3 + undergens + numconstraints + rescounter] = max(6 * ramp, orignsrsrrsregupaward) # NSRS + RRS + RegUp <= 30 minute ramp rate up
        timestampdata[resource]['NSRSRRSRegEnergyCapUsed'] = hlist[2 * numvar + 7 * numresources + 3 + undergens + numconstraints + rescounter]

        #  Energy + (5/7) * RegUp <= Output + 5 minute ramp rate up - LDL
        hlist[2 * numvar + 8 * numresources + 3 + undergens + numconstraints + rescounter] =  max(timestampdata[resource]['Output'] + ramp - timestampdata[resource]['LDL'], 0.0)
        timestampdata[resource]['E + 5/7 Reg Used'] = hlist[2 * numvar + 8 * numresources + 3 + undergens + numconstraints + rescounter]


        #  If 'ONRR', Energy + NSRS + RegUp <= LDL, else <= 10,000 (arbitrarilty large number to ensure non-binding)
        if status != 'ONRR':
            hlist[2 * numvar + 9 * numresources + 2 + 2 * undergens + 2 * numconstraints + rescounter] =  10000
        # if resource has sold NSRS before, can sell now, otherwise can not.
        if nsrsmwh or orignsrsaward:
            hlist[2 * numvar + 10 * numresources + 2 + 2 * undergens + 2 * numconstraints + rescounter] =  10000            

        rescounter += 1

    for reginsuff in range(undergens - 1):
        hlist[2 * numvar + 9 * numresources + 3 + undergens + numconstraints + reginsuff] = permdata['UnderGen'][reginsuff][0]
        print("h",reginsuff,2 * numvar + 9 * numresources + 3 + undergens + numconstraints + reginsuff)
        # Reg Insufficient segment energy <= segment MWs

    for constraintnum in range(numconstraints):
        constraintid = constraintids[constraintnum]
        hlist[2 * numvar + 9 * numresources + 2 + 2 * undergens + numconstraints + constraintnum] = daydata[timestamp]['NetworkConstraints'][constraintid][1] #get Limit
        for resource in timestampdata['LDLResources']:
            try:
                sf = daydata[timestamp]['ShiftFactors'][constraintid][resource]
                ldl = timestampdata[resource]['LDL']
                hlist[2 * numvar + 9 * numresources + 2 + 2 * undergens + numconstraints + constraintnum] -= sf * ldl #subtract out LDL component
            except KeyError:
                pass

    timestampdata['h'] =   matrix(hlist)

    """
    A is a 4 rows by (numvar + 1 + numresources + 1 + numresources + 1 + numresources + undergens + numconstraints) columns matrix
    the first row is for the total GTBD constraint
    the second row is for the total NSRS constraint
    the third row is for the total RRS constraint
    the fourth row is for the total RegUp constraint
    """

    Alists = []
    for column in range(numvar + 1): #handle energy awards and energy insufficiency
        Alist = [1.0,0.0,0.0,0.0]
        Alists.append(Alist)

    for column in range(numresources + 1): #handle nonspin awards and nonspin insufficiency
        Alist = [0.0,1.0,0.0,0.0]
        Alists.append(Alist)

    for column in range(numresources + 1): #handle rrs awards and insufficiency
        Alist = [0.0,0.0,1.0,0.0]
        Alists.append(Alist)

    for column in range(numresources + undergens): #handle reg up awards and insufficiency
        Alist = [0.0,0.0,0.0,1.0]
        Alists.append(Alist)

    for column in range(numconstraints):
        Alist = [0.0,0.0,0.0,0.0]
        Alists.append(Alist)
    
    timestampdata['A'] = matrix(Alists)

    """
    b is a 4 row by 1 column vector/matrix
    the first row is GTBD - SUM of LDLs
    the second row is total NSRS needed
    the third row os total RRS needed
    the fourth row is total RegUp needed
    """

    syscon = permdata['SystemConditions'][timestamp]
    gtbd = syscon['GTBD']
    ldlsum = timestampdata['LDLSum']
    nsrs = syscon['NSRSMW']
    rrs = syscon['RRSMW']
    regup = syscon['RegUpMW']
    regdown = syscon['RegDownMW']
    timestampdata['b'] = matrix([gtbd - ldlsum, nsrs, rrs, regup])
    print([gtbd - ldlsum, nsrs, rrs, regup])

    return timestampdata

def cooptimize(timestampdata):

    P = timestampdata['P']
    q = timestampdata['q']
    G = timestampdata['G']
    h = timestampdata['h']
    A = timestampdata['A']
    b = timestampdata['b']
    timestampdata['Solution']=solvers.qp(P, q, G, h, A, b)

    return timestampdata

def consolidate_cooptimization_results(timestamp,timestampdata,daydata,permdata,pfudge):

    sol = timestampdata['Solution']
    optstatus = sol['status']
    numvar = timestampdata['Counter']
    undergens = timestampdata['UnderGens']
    overgens = timestampdata['OverGens']
    numconstraints = timestampdata['NumConstraints']
    constraintids = timestampdata['ConstraintIDs']
    numresources = len(timestampdata['Resources']) 
    summary = []
    pricecheck = []

    timestampdata['ASPrices'] = {} #initialize ASPrice dict
    timestampdata['ASAwards'] = {} #initialize ASAwards dict
    for troublepoint in timestampdata['Trouble']:
        print(troublepoint,sol['x'][troublepoint[1]],timestampdata['q'][troublepoint[1]])

    if optstatus == 'optimal':
        optstatus = optstatus + str(pfudge)
        try:
    #add shadow price to daydata[timestamp]['NetworkConstraints'][constraintid] list
            for constraintnum in range(numconstraints):
                constraintid = constraintids[constraintnum]
                mathflow = 0.0
                for resource in timestampdata['LDLResources']:
                    try:
                        mathflow += daydata[timestamp]['ShiftFactors'][constraintid][resource] * timestampdata[resource]['LDL']
                    except KeyError:
                        pass
                for counter in range(numvar):
                    try:
                        resource = timestampdata[counter]['ResourceName']
                        x = sol['x'][counter] #energy award
                        mathflow += daydata[timestamp]['ShiftFactors'][constraintid][resource] * x
                    except KeyError:
                        pass
                violatedmw = sol['x'][numvar + 1 + numresources + 1 + numresources + 1 + numresources + undergens  + constraintnum]
                if True:
                    daydata[timestamp]['NetworkConstraints'][constraintid].append(sol['z'][2 * numvar + 9 * numresources + numconstraints + 2 * undergens + 2 + constraintnum])
                    daydata[timestamp]['NetworkConstraints'][constraintid].append(mathflow)
                    daydata[timestamp]['NetworkConstraints'][constraintid].append(violatedmw)
        except KeyError:
            pass
        permdata['SystemConditions'][timestamp]['BPSum'] = 0.0
        for resource in timestampdata['LDLResources']: #get prices for each resource that was dispatched
            timestampdata[resource]['EnergyPrice'] = -1*sol['y'][0] #system lambda
            for constraintnum in range(numconstraints): #add in the congestion component
                constraintid = constraintids[constraintnum]
                try:
                    timestampdata[resource]['EnergyPrice'] -= daydata[timestamp]['ShiftFactors'][constraintid][resource] * \
                      daydata[timestamp]['NetworkConstraints'][constraintid][5]
                except (KeyError, IndexError):
                    pass
            timestampdata[resource]['EnergyAward'] = timestampdata[resource]['LDL']
            timestampdata[resource]['NSRSAward'] = 0.0
            timestampdata[resource]['RRSAward'] = 0.0
            timestampdata[resource]['RegUpAward'] = 0.0
            timestampdata[resource]['ProductionCost'] = 0.0
            permdata['SystemConditions'][timestamp]['BPSum'] += timestampdata[resource]['LDL']
        timestampdata['ASPrices']['NSRS'] = -1*sol['y'][1]  #get the AS prices
        timestampdata['ASPrices']['RRS'] = -1*sol['y'][2]
        timestampdata['ASPrices']['RegUp'] = -1*sol['y'][3]

        nsrs = 0.0
        rrs = 0.0
        reg = 0.0
        for resourcenum in range(numresources):
            nsrs += sol['x'][numvar + 1 + resourcenum]
            rrs += sol['x'][numvar + 1 + numresources + 1 + resourcenum]
            reg += sol['x'][numvar + 1 + numresources + 1 + numresources + 1 + resourcenum]
        timestampdata['ASAwards']['NSRS'] = nsrs
        timestampdata['ASAwards']['RRS'] = rrs
        timestampdata['ASAwards']['RegUp'] = reg
        timestampdata['UnitProductionCost'] = 0.0


        for counter in range(numvar):
            x = sol['x'][counter]
            P = timestampdata['P'][counter,counter]
            q = timestampdata['q'][counter]
            resource = timestampdata[counter]['ResourceName']
            timestampdata[resource]['EnergyAward'] += x
            timestampdata[resource]['ProductionCost'] += 0.5 * P * x * x + q * x
            timestampdata['UnitProductionCost'] += 0.5 * P * x * x + q * x
            permdata['SystemConditions'][timestamp]['BPSum'] += x #for price checking system conditions file
        for resourcenum in range(numresources):
            resource = timestampdata['Resources'][resourcenum]
            timestampdata[resource]['NSRSAward'] += sol['x'][numvar + 1 + resourcenum]
            timestampdata[resource]['RRSAward'] += sol['x'][numvar + 1 + numresources + 1 + resourcenum]
            timestampdata[resource]['RegUpAward'] += sol['x'][numvar + 1 + numresources + 1 + numresources + 1 + resourcenum]


        permdata['SystemConditions'][timestamp]['PBPCMW'] = sol['x'][numvar]
        permdata['SystemConditions'][timestamp]['SystemLambda'] = -1 * sol['y'][0]
    else:
        permdata['SystemConditions'][timestamp]['BPSum'] = 'NA' #for price checking system conditions file
        permdata['SystemConditions'][timestamp]['PBPCMW'] = 'NA'
        permdata['SystemConditions'][timestamp]['SystemLambda'] = 'NA'
        if True:
            for constraintnum in range(numconstraints):
                constraintid = constraintids[constraintnum]
                daydata[timestamp]['NetworkConstraints'][constraintid].append('NA')#for binding constraint file
                daydata[timestamp]['NetworkConstraints'][constraintid].append('NA')
                daydata[timestamp]['NetworkConstraints'][constraintid].append('NA')

    for row in daydata['60DayGRD']:
        grdtimestamp = (row[0],row[1])
        grdresource = (row[2])        
        if True:
            if grdtimestamp == timestamp:
                try:
                    hdl = timestampdata[grdresource]['HDL']
                    hsl = timestampdata[grdresource]['HSL']
                    ldl = timestampdata[grdresource]['LDL']
                    output = timestampdata[grdresource]['Output']
                    try:
                        rampup = daydata[timestamp]['RampRates'][grdresource]
                    except KeyError:
                        rampup = max((hdl-output),daydata['ResourceParameters'][grdresource]['RampRateUp']) #log
                    try:
                        rrsunitcommit = daydata['RRSCommit'][timestamp][grdresource] >= 0.1
                    except KeyError:
                        rrsunitcommit = False                      
                    try:
                        freqaward = timestampdata[grdresource]['PreviousAward']['RRS'] >= 0.1
                    except KeyError:
                        freqaward = False
                    try:
                        nfrc = timestampdata[grdresource]['NFRCUsed']
                    except KeyError:
                        nfrc = 0.0
                    nsrsoffcap = max(hsl-ldl-nfrc,0.0)
                    try:
                        rrsoffcap = timestampdata[grdresource]['RRSOffCapUsed']
                    except KeyError:
                        rrsoffcap = 0.0
                    try:
                        regoffcap = timestampdata[grdresource]['RegUpOffCapUsed']
                    except KeyError:
                        regoffcap = 0.0
                    try:
                        hslminusldl = timestampdata[grdresource]['HSL-LDLUsed']
                    except KeyError:
                        hslminusldl = 0.0
                    try:
                        hdlminusldl = timestampdata[grdresource]['E + 5/7 Reg Used']
                    except KeyError:
                        hdlminusldl = 0.0
                    try:
                        rrsregenoffcap = timestampdata[grdresource]['RRSRegEnergyCapUsed']
                    except KeyError:
                        rrsregenoffcap = 0.0
                    try:
                        nsrsrrsregenoffcap = timestampdata[grdresource]['NSRSRRSRegEnergyCapUsed']
                    except KeyError:
                        nsrsrrsregenoffcap = 0.0
                    row.append(rrsunitcommit)
                    row.append(timestampdata[grdresource]['ProductionCost'])
                    row.append(timestampdata[grdresource]['EnergyAward'])
                    row.append(timestampdata[grdresource]['NSRSAward'])
                    row.append(timestampdata[grdresource]['RRSAward'])
                    row.append(timestampdata[grdresource]['RegUpAward'])
                    row.append(timestampdata[grdresource]['EnergyPrice'])
                    row.append(timestampdata['ASPrices']['NSRS'])
                    row.append(timestampdata['ASPrices']['RRS'])
                    row.append(timestampdata['ASPrices']['RegUp'])
                    row.append(hslminusldl)
                    row.append(rrsoffcap)
                    row.append(regoffcap)
                    row.append(rampup)
                    row.append(nfrc)
                    row.append(hdl)
                    row.append(hdlminusldl)
                    row.append(rrsregenoffcap)
                    row.append(nsrsrrsregenoffcap)
                except KeyError:
                    row.append('FALSE')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')
                    row.append('NA')

    summary.append(optstatus)
    if True:
        try:    pricecheck.append(permdata['SystemConditions'][timestamp]['BPSum'])
        except KeyError:  pricecheck.append('NA')
        try:    pricecheck.append(permdata['SystemConditions'][timestamp]['PBPCMW'])
        except KeyError:  pricecheck.append('NA')
        try:    pricecheck.append(permdata['SystemConditions'][timestamp]['SystemLambda'])
        except KeyError:  pricecheck.append('NA')
        try:    pricecheck.append(timestampdata['ASAwards']['NSRS'])
        except KeyError:  pricecheck.append('NA')
        try:    pricecheck.append(timestampdata['ASAwards']['RRS'])
        except KeyError:  pricecheck.append('NA')
        try:    pricecheck.append(timestampdata['ASAwards']['RegUp'])
        except KeyError:  pricecheck.append('NA')
        try:    pricecheck.append(timestampdata['ASPrices']['NSRS'])
        except KeyError:  pricecheck.append('NA')
        try:    pricecheck.append(timestampdata['ASPrices']['RRS'])
        except KeyError:  pricecheck.append('NA')
        try:    pricecheck.append(timestampdata['ASPrices']['RegUp'])
        except KeyError:  pricecheck.append('NA')

        try:    summary.append(-1*sol['y'][0])
        except KeyError:  summary.append('NA')
        try:    summary.append(sol['primal objective'])
        except KeyError:  summary.append('NA')
        try:    summary.append(timestampdata['UnitProductionCost'])
        except KeyError:  summary.append('NA')
        try:    summary.append(permdata['SystemConditions'][timestamp]['GTBD'] *-1*sol['y'][0])
        except KeyError:  summary.append('NA')
        try:    summary.append(timestampdata['ASAwards']['RegUp'])
        except KeyError:  summary.append('NA')
        try:    summary.append(timestampdata['ASAwards']['RRS'])
        except KeyError:  summary.append('NA')
        try:    summary.append(timestampdata['ASAwards']['NSRS'])
        except KeyError:  summary.append('NA')
        try:    summary.append(timestampdata['ASPrices']['RegUp'])
        except KeyError:  summary.append('NA')
        try:    summary.append(timestampdata['ASPrices']['RRS'])
        except KeyError:  summary.append('NA')
        try:    summary.append(timestampdata['ASPrices']['NSRS'])
        except KeyError:  summary.append('NA')
        try:    summary.append(permdata['SystemConditions'][timestamp]['LoadNSRS'])
        except KeyError:  summary.append('NA')
        try:    summary.append(permdata['SystemConditions'][timestamp]['LoadRRS'])
        except KeyError:  summary.append('NA')
        try:    summary.append(permdata['SystemConditions'][timestamp]['LoadRegUp'])
        except KeyError:  summary.append('NA')
        try:    summary.append(permdata['SystemConditions'][timestamp]['OutNSRS'])
        except KeyError:  summary.append('NA')
        try:    summary.append(permdata['SystemConditions'][timestamp]['OutRRS'])
        except KeyError:  summary.append('NA')
        try:    summary.append(permdata['SystemConditions'][timestamp]['OutRegUp'])
        except KeyError:  summary.append('NA')
            
    for item in summary:
        daydata['Summary'][-1].append(item)

    for item in pricecheck:
        daydata['PriceCheck'][-1].append(item)

    return daydata

def save_output_for_day(wdname,daydata,permdata):
    ##################################
    ######### Save GRD file ##########
    ##################################
    parentdir = os.path.dirname(wdname)
    resultsdir = os.path.join(parentdir,'nfrcresults')
    os.chdir(resultsdir)
    filename = 'output' + os.path.basename(wdname) + '.csv'
    if permdata['PythonVersion'] == 3:
        with open(filename, 'w',newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            for datum in daydata['60DayGRD']:
                spamwriter.writerow(datum)
    else:
        with open(filename, 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            for datum in daydata['60DayGRD']:
                spamwriter.writerow(datum)
    if permdata['CloudStore']:
        origin = os.path.join(resultsdir,filename)
        destination = os.path.join(permdata['GSBucket'],filename)
        subprocess.check_call(['gsutil','cp','-r',origin,destination])
        subprocess.check_call(['rm',origin])
        

    ##################################            
    ####### Save summary file ########
    ##################################

    filename = 'summary' + os.path.basename(wdname) + '.csv'
    headerrow = ['SCED_TIMESTAMP','REPEATED_HOUR_FLAG','SCED_OPT_STATUS','GTBD','LDL_SUM','SECD_SYSTEM_LAMBDA','SCED_PRIMAL',
    'SCED_PRODUCTION_COST','SCED_ENERGY_COST','ORIG_REG_UP_MW','ORIG_RRS_MW','ORIG_NSRS_MW','ORIG_REG_UP_PRICE','ORIG_RRS_PRICE','ORIG_NSRS_PRICE',
    'MIP_OPT_STATUS','COOPT_OPT_STATUS','COOPT_SYSTEM_LAMBDA','COOPT_PRIMAL','COOPT_UNIT_PROD_COST','COOPT_ENERGY_COST',
    'COOPT_REG_UP_MW','COOPT_RRS_MW','COOPT_NSRS_MW','COOPT_REG_UP_PRICE','COOPT_RRS_PRICE','COOPT_NSRS_PRICE','LOAD_NSRS','LOAD_RRS','LOAD_REG_UP',
    'OUT_NSRS','OUT_RRS','OUT_REG_UP']
    if permdata['PythonVersion'] == 3:
        with open(filename, 'w', newline = '') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter = ',')
            spamwriter.writerow(headerrow)
            for row in daydata['Summary']:
                spamwriter.writerow(row)
    else:
        with open(filename, 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter = ',')
            spamwriter.writerow(headerrow)
            for row in daydata['Summary']:
                spamwriter.writerow(row)
    if permdata['CloudStore']:
        origin = os.path.join(resultsdir,filename)
        destination = os.path.join(permdata['GSBucket'],filename)
        subprocess.check_call(['gsutil','cp','-r',origin,destination])
        subprocess.check_call(['rm',origin])

    ##################################            
    ## Save binding constraint file ##
    ##################################

    filename = 'binding_constraints' + os.path.basename(wdname) + '.csv'    
    headerrow = ['SCED_TIMESTAMP','REPEATED_HOUR_FLAG','CONSTR_ID','CONTINGENCY_NAME',
                    'CONSTR_NAME','MATHLIMIT','MATHFLOW','VIOLATED_MW','SHADOW_PRICE','MAX_SHADOW_PRICE','SCED_SHADOW_PRICE']
    if permdata['PythonVersion'] == 3:
        with open(filename, 'w', newline = '') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter = ',')
            spamwriter.writerow(headerrow)
            for constraint in daydata['NetworkConstraints']:
                timestamp = constraint[0]
                scedtimestamp = timestamp[0]
                repeatedhourflag = timestamp[1]
                constraintid = constraint[1]
                value = daydata[timestamp]['NetworkConstraints'][constraintid]
                contingency = value[2]
                constraintname = value[3]
                mathlimit = constraint[3]
                mathflow = value[6]
                try:
                    if value[7] < 0.1: violatedmw = 0.0
                    else: violatedmw = value[7]
                except TypeError:
                    violatedmw = value[7]
                try:
                    if value[5] < 0.1: shadowprice = 0.0
                    else: shadowprice = value[5]
                except TypeError:
                    shadowprice = value[5]
                maxshadowprice = constraint[2]
                try:
                    if value[4] < 0.1: scedshadowprice = 0.0
                    else: scedshadowprice = value[4]
                except TypeError:
                    scedshadowprice = value[4]
                row = [scedtimestamp,repeatedhourflag,constraintid,contingency,
                       constraintname,mathlimit,mathflow,violatedmw,shadowprice,maxshadowprice,scedshadowprice]
                spamwriter.writerow(row)

    else:
        with open(filename, 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter = ',')
            spamwriter.writerow(headerrow)
            for constraint in daydata['NetworkConstraints']:
                timestamp = constraint[0]
                scedtimestamp = timestamp[0]
                repeatedhourflag = timestamp[1]
                constraintid = constraint[1]
                value = daydata[timestamp]['NetworkConstraints'][constraintid]
                contingency = value[2]
                constraintname = value[3]
                mathlimit = constraint[3]
                mathflow = value[6]
                try:
                    if value[7] < 0.1: violatedmw = 0.0
                    else: violatedmw = value[7]
                except TypeError:
                    violatedmw = value[7]
                try:
                    if value[5] < 0.1: shadowprice = 0.0
                    else: shadowprice = value[5]
                except TypeError:
                    shadowprice = value[5]
                maxshadowprice = constraint[2]
                try:
                    if value[4] < 0.1: scedshadowprice = 0.0
                    else: scedshadowprice = value[4]
                except TypeError:
                    scedshadowprice = value[4]
                row = [scedtimestamp,repeatedhourflag,constraintid,contingency,
                       constraintname,mathlimit,mathflow,violatedmw,shadowprice,maxshadowprice,scedshadowprice]
                spamwriter.writerow(row)
    if permdata['CloudStore']:
        origin = os.path.join(resultsdir,filename)
        destination = os.path.join(permdata['GSBucket'],filename)
        subprocess.check_call(['gsutil','cp','-r',origin,destination])
        subprocess.check_call(['rm',origin])

    ##################################            
    ###### Save price check file #####  -only for ercot price validation - commented out for production
    ##################################

    """filename = 'price_check_data' + os.path.basename(wdname) + '.csv'    
    headerrow = ['SCED_TIMESTAMP','REPEATED_HOUR_FLAG','GTBD','TOTAL_BP','PBPC_MW','SYSTEM LAMBDA','Cleared_NSRS',
                 'Cleared_RRS','Cleared_Reg_Up','NSRS_Price','RRS_Price','Reg_Up_Price']
    if permdata['PythonVersion'] == 3:
        with open(filename, 'w', newline = '') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter = ',')
            spamwriter.writerow(headerrow)
            for row in daydata['PriceCheck']:
                spamwriter.writerow(row)
    else:
        with open(filename, 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter = ',')
            spamwriter.writerow(headerrow)
            for row in daydata['PriceCheck']:
                spamwriter.writerow(row)
    if permdata['CloudStore']:
        origin = os.path.join(resultsdir,filename)
        destination = os.path.join(permdata['GSBucket'],filename)
        subprocess.check_call(['gsutil','cp','-r',origin,destination])
        subprocess.check_call(['rm',origin])"""
                                 
# Loop through a day's worth of optimizations
def run_day(wdname,permdata,dayfolder):
    pfudges = [0.0,0.001,0.1,0.2,0.3,0.4,0.5]
    firststart = time.time() #get start to measure run time for the day of optimization
    resourceparameters = get_resource_parameters(wdname) #get resource parameters (ramp, can they sell AS,etc.)
    daydata = get_grd_file(wdname,resourceparameters,permdata,dayfolder) #load up the 60 day sced generator report
    daydata = get_binding_constraint_data(wdname,daydata) #get cdr12302 information
    daydata = get_shift_factors(wdname,daydata) #get cdr12354 information
    daydata = get_ramp_rate(wdname,daydata) #get ramp rate data
    #daydata = get_resource_parameters(wdname,daydata) #get resource parameters (ramp, can they sell AS,etc.)
    daydata = get_nfrc_data(wdname,daydata) #get nfrc data
    for timestamp in daydata['TimeStamps']: #iterate through the timestamps and optimize
        starttimestamp = time.time() #get start time to mesure length of optimization of one sced interval.
        write_log_entry([starttimestamp -firststart,timestamp[0],'start optimization'],wdname,permdata)
        #sced re-run section
        timestampdata = get_sced_segment_info(timestamp,daydata)
        daydata = true_up_shift_factors(daydata,timestamp,timestampdata)
        timestampdata = create_sced_matrices(timestamp,timestampdata,daydata,permdata)
        timestampdata = re_run_sced(timestampdata)
        daydata = consolidate_sced_results(timestamp,timestampdata,daydata,permdata)
        welfaretobeat = daydata['WelfareToBeat'][timestamp]
        gc.collect()
        #first try cooptimization section
        optimizationround = 2
        timestampdata = get_mip_cooptimization_segment_info(timestamp,daydata,permdata,optimizationround)
        if timestamp not in [('07/10/2017 15:00:17','N'),('08/12/2017 11:00:17','N')]:
            daydata = mip_coopt(timestamp,timestampdata,daydata,permdata,optimizationround)
        else:
            daydata[timestamp]['MIPStatus'] = 'not optimal'
            daydata['Summary'][-1].append('not optimized')
            daydata['RRSCommit'][timestamp] = {}
            
        print(timestamp, 'day so far: ',time.time() - firststart)
        for pfudge in pfudges:
            print('Opt Round: ',optimizationround,' PFUDGE: ',pfudge)
            timestampdata = get_cooptimization_segment_info(timestamp,daydata,permdata,optimizationround)
            timestampdata = create_cooptimization_matrices(timestamp,timestampdata,daydata,permdata,pfudge)
            timestampdata = cooptimize(timestampdata)
            if timestampdata['Solution']['status'] == 'optimal' :
                daydata = consolidate_cooptimization_results(timestamp,timestampdata,daydata,permdata,pfudge)
                break
    save_output_for_day(wdname,daydata,permdata) #save output files for the day

def main():
    wdname = sys.argv[1] #get working directory where information is stored from command line entry
    try:
        startday = sys.argv[2] #get day to start analysis on
    except IndexError:
        startday = 'none'
    try:
        numofdays = int(sys.argv[3]) #get number of days to work
    except IndexError:
        numofdays = 'none'
    try:
        gsbucket = sys.argv[4]
    except IndexError:
        gsbucket = 'none'
        cloudstore = False
    if gsbucket == 'none': cloudstore = False
    else: cloudstore = True
    os.chdir(wdname) #change to correct data directory
    permdata = get_system_conditions(wdname) #get the system condition data for the year (GTBD, etc.)
    permdata = get_quickstarts(wdname,permdata) #get the list of units that operated with 'OFFQS' status in 2017
    permdata['CloudStore'] = cloudstore
    permdata['GSBucket'] = gsbucket
    permdata['PythonVersion'] = sys.version_info[0] #get python version (need to know when saving files)
    permdata = get_powerbalance(permdata) #get powerbalance penalty curve parameters
    daylistfilename = os.path.join(wdname, 'daylist.csv')
    days = open_days(daylistfilename) #get list of days to be analyzed
    analyze = False
    endnum = 1000000
    if startday == 'none':
        for dayfolder in days: #loop through day folders to be analyzed
            fulldayfolder = os.path.join(wdname,dayfolder[0])
            permdata = get_load_as(fulldayfolder,permdata)
            run_day(fulldayfolder, permdata, dayfolder) #cooptimize and compare a day a day
            gc.collect() #collect the garbage
    else:
        for dayfolder in days: #loop through day folders to be analyzed
            if dayfolder[0] == startday:
                analyze = True
                endnum = dayfolder[1] + numofdays
            if dayfolder[1] >= endnum:
                analyze = False
            if analyze:
                fulldayfolder = os.path.join(wdname,dayfolder[0])
                permdata = get_load_as(fulldayfolder,permdata)
                run_day(fulldayfolder, permdata, dayfolder[0]) #cooptimize and compare a day a day
                gc.collect() #collect the garbage
    if permdata['CloudStore']: subprocess.check_call(['sudo','poweroff'])

if __name__ == "__main__":
    main()
