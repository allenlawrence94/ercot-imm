
########################################################
########################################################
"""
Coding started 3/13/2019
by MAT&T
"""
########################################################
########################################################

from rtcoopt import *

def min_run_day(wdname,permdata,dayfolder):
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

def min_run(wdname, startday, numofdays):
    startday = 'none'
    gsbucket = 'none'
    cloudstore = False
    os.chdir(wdname)
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
    for dayfolder in days: #loop through day folders to be analyzed
        if dayfolder[0] == startday:
            analyze = True
            endnum = dayfolder[1] + numofdays
        if dayfolder[1] >= endnum:
            analyze = False
        if analyze:
            fulldayfolder = os.path.join(wdname,dayfolder[0])
            permdata = get_load_as(fulldayfolder,permdata)
            min_run_day(fulldayfolder, permdata, dayfolder[0]) #cooptimize and compare a day a day
            gc.collect() #collect the garbage
