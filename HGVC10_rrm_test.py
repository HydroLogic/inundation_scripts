'''
_________________________________________________________________________________________________

Script Name: Valley Bottom Classification (including Break-in-Slope (BiS), Q100 flood level, and
hillslope steepness calculators)
Description: This model classifies the valley bottoms and hillslopes at stream-segment scale in
    the following steps:
        1. Evaluates the 100 year flood using a reach averaged Manning's equation algorithm.
        2. Selects the points of greatest concavity in the valley bottom or "Break-in-Slope"
        3. Clips the BiS layer by the Q100 layer leaving a geomorphically and hydrologically
          defined valley bottom
        4. Evaluates the slope of adjacent hillsides and categorizes each side of every segment
        5. Combines all the valley bottom and hillslope information to classify the characteristics
          of the fluvial riparian zone for each stream segment
The following files are required:
        1. A shapefile of all stream segments to be analyzed (with stream lines divided at junction
          nodes and at changes in average-slope class)  NOTE: these must come from a DEM, not
            bluelines
        2. The filled-DEM that was the source of the stream delineation
        3. An estimated Manning's roughness value

Created By:   Daniel W Baker
              With the assistance of:
                  Erick Carlson
                  Brian Bledsoe
                  Thomas Flowe

Using modifications of:
      http://forums.esri.com/Thread.asp?c=93&f=1728&t=191496&g=1 (adapted for shapefile loops)
      http://arcscripts.esri.com/details.asp?dbid=15756 (adapted for stream extensions to split hillslopes)
      CreateFeaturesFromTextFile (Python script built into ArcGIS 9.3)
      
Updates
6/5/2013    - Started changeover from v9.3 to v10.0 as guided by
            http://www.geospatialtraining.com/blog/index.php/checklist-for-updating-python-scripts-from-arcgis-9-3-to-arcgis-10/
            -replaced gp. with arcpy.ASCII3DToFeatureClass_3d
6/6/2013    - Completed checklist (from geospatialtraining.com link above)
            - Had to check syntax of all Geoprocessing commands
                -Now save to a file (instead of pre-initate)
                -some need arcpy. others not, and
                -had to validate format (some must be set equal to variable, others still inc   luded output)
             - Progressed down through LINE 804
                 -PYTHON .split needs str()
6.7.2014    -Completed changeover from v9.3 to v10.0
            - Note: v10.0 handles geoprocessing functions from Spatial Analyst different from the other
            groups (basically inverting the process of creating a new dataset).
                - For ALL OTHERS:
                    (1) Initate the file
                    (2) Proceed function name by 'arcpy.' and be careful with capitialization and spelling
                    they both must be exact
                - Spatial Analyst functions (EucAllocation, CostDistance, etc)
                    (1) DO NOT initate file
                    (2) no need for 'arcpy.' preceeding function name
                    (3) set function equal to variable (no output files in parentheses of parameters)
                    (4) .save() the variable to the desired location
7.25.2013   -Folded in RSAC edits as sent to Dave Merritt 6/26/2013
             -Commented Edits with "RSAC" for searchability
7.29.2013   -Fixed Hill Slope categorization (due to previous elseif statement it would assume at least
                one of the hillslopes was moderate or steep)
4.7.2014    -Worked to incorporate changes from Dave Theobald including - note that all changes were in segment creation
                (hence not in this code)
4.17.2014   -Put all temp shapefiles in a folder within TEMP (/segs)
6.5.2014    -Created addtional local (TEMP/SEGS) files for error checking. Files include
                LINE    FILE
                740	    slope_pct_100
                1050	vw_final
                1058	Curvature
                1069	BiS
                1132	HG_final_sides
                1166    Distance from Stream
                1302    Hill Split
                1293	Hill_RT

6.6.2014    -Curvature Cutoff (for BiS) line 
                 -Previously took range of curvature and used cutoff of 20% of range above min. Changed to (for BiS)
                 to 30% above minimum value curvature (always neg). Curvatures are defined as negative when concave upward.
            -Removed cutline if statement (line 1300) as cutlines are now stable

6.9.2014    -Removed 6.5.2014 TEMP/SEGS files
            -Redefined extent of valley_sections to that of the DEM


__________________________________________________________________________________________________
'''

#print ' '

        
# ##############################################################
# A. General script setup

#Import standard library modules
import sys, os, csv, string, random, time, logging, flog, traceback, linecache 
import datetime

thentime = datetime.datetime.now()  # used to note start time of model run

import arcpy   # New for ArcGIS v.10
from arcpy import env
from arcpy.sa import *

arcpy.ResetEnvironments()

####################################################################################
#Check out necessary extensions
arcpy.CheckOutExtension("Spatial")
arcpy.CheckOutExtension("3D")

# Set OverWriteOutput to 1 to copy over existing outputs, 0 to through exceptions for existing output
arcpy.env.OverWriteOutput = True

# ##############################################################
# B. Establish input parameter values

base  = "C:/GIS"
root    = "srout"
folder = "A001"
userworkfolder = base + '/' + root
# userworkfolder = 'C:/GIS/SRouteMedbow/'

# Create userworkspace and userworkspace/temp directories
filenum = 1
dir_exists = False 
while not dir_exists:
    if filenum > 100:
        break
    try:
        userworkspace = userworkfolder + '/' + "B" + str(filenum).zfill(3)
        os.mkdir(userworkspace) # Create the workspace
        dir_exists = True
    except:
        filenum += 1
##        print "  Working Folder =" + str(filenum).zfill(3)
        pass

# Hard coded locations
##inDEM = 'C:/GIS/SRouteMedbow/filled_dem.img'                      # Filled DEM 
##valley_section = 'C:/GIS/SRouteMedbow/A009/val_segs_shp.shp'      # Set the input shape file
##da_km = 'C:/GIS/SRouteMedbow/da_km2.img'                          # Flow accumulation in km2
##strm_cells = 'C:/GIS/SRouteMedbow/strm_net.img'                   # Raster stream cells
##valley_block = 'C:/GIS/SRouteMedbow/A009/valley_block.shp'        # Raster watershed cells
##Q100_raster =  'C:/GIS/SRouteMedbow/q100_cms.tif'                 # Q100 raster (same grid size and extent as DEM)

# Locations using userworkspace and root
inDEM           = userworkfolder + '/' + folder + '/' + root + "_fdem"
valley_section  = userworkfolder + '/' + folder + '/' + root +  "_segs" + ".shp"
da_km           = userworkfolder + '/' + folder + '/' + root +  "_da_km"
strm_cells      = userworkfolder + '/' + folder + '/' + root +  "_strm"
valley_block    = userworkfolder + '/' + folder + '/' + root +  "_blks" + ".shp"
Q100_raster     = userworkfolder +  '/' + "q100_cms.tif"

start_ARCID = 0    # Set to zero to start at beginning !!KEEP AT ZERO, UNKONWN ERRORS WHEN STARTING ELSEWHERE
seg_max =10     # Maximum ARCID of segments to evaluate in this model run

alpha = 2.26        # coefficient - BF channel width coefficients (from Faustini/WEMAP)
beta = 0.31         # power - BF channel width coefficients (from Faustini/WEMAP)
hill_buff_dist = 250.0  # Distance beyond Hydro-Geo valley bottom to evaluate hill slopes
HSthresh_up = 0.70       # Upper Hillslope slope threshold
HSthresh_low = 0.30       # Lower hillslope slope threshold (i.e. 0.30 = 30%)
debris_RO = 15.0    # Input from user
glacial_min_elev = 2500 # Lower limit of glacial influence
analyst = 'Dan Baker'

Q_tol = 2.0         # Adjustment for certainty of Q100
mannings_n= 0.03    # Estimated Manning's roughness value
delete_temp = "NO"  # Delete files from /temp folder
diff_tol = 0.1         # Acceptible tolerance for % diff btw Q_est & Q_calc
flood_min = 0.5     # Minimum flood elevation (think of as vertical resolution of DEM)
iter_max = 4        # Num. of iterations w/ depth lower then flood_min before exiting Q100 calculations

# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
##userworkspace = sys.argv[1]        # Folder used to store data                            
##valley_section =sys.argv[2]         # Set the input shape file
##inDEM =         sys.argv[3]         # Filled DEM 
##da_km =         sys.argv[4]         # Flow accumulation in km2
##strm_cells =    sys.argv[5]         # Raster stream cells
##start_ARCID =   int(sys.argv[6])         # Set to zero to start at beginning
##seg_max =       int(sys.argv[7])         # Number of segments to evaluate in this model run
##Q100_raster =   sys.argv[8]        # Raster of Q100 stream cell values
##Q_tol =         float(sys.argv[9])    # Adjustment for certainty of Q100
##mannings_n=     float(sys.argv[10])  # Estimated Manning's roughness value
##alpha =         float(sys.argv[11])     # coefficient - BF channel width coefficients (from Faustini/WEMAP)
##beta =          float(sys.argv[12])     # exponent - BF channel width coefficients (from Faustini/WEMAP)
##hill_buff_dist =float(sys.argv[13])    # Distance beyond Hydro-Geo valley bottom to evaluate hill slopes
##HSthresh_up =        float(sys.argv[14])    # Upper hillslope slope threshold (i.e. 0.30 = 30%)
##HSthresh_low =        float(sys.argv[15])    # Lower Hillslope slope threshold
##debris_RO =     float(sys.argv[16])    # Input from user
##glacial_min_elev=float(sys.argv[17])   # Lower limit of glacial influence

# Open error log
log_file = userworkspace + "/ErrorLog.txt"
flog = open(log_file, 'w')
# logging.basicConfig(filename=log_file, level=logging.DEBUG)

arcpy.env.workspace = userworkspace
print '  Workspace is set to:', str(env.workspace)
flog.write("START TIME: "+ '{:%Y-%m-%d %H:%M:%S}'.format(thentime))

try:
    os.mkdir(userworkspace + '/temp') # Create the /TEMP folder
    os.mkdir(userworkspace + '/temp/seg') # Create the /TEMP/seg folder
except:
    print '  ERROR - establishing temp directory (may already exist)'
    pass
arcpy.AddMessage(arcpy.GetMessages(2))

flog.write( '\n' + 'ANALYST: ' + analyst+ '\n')
flog.write( 'INPUTS to model run'+ '\n')
flog.write( '  Workspace is set to:' + str(env.workspace)+ '\n')
flog.write( '  valley_section=' + valley_section+ '\n')
flog.write( '  inDEM =' +inDEM+ '\n')
flog.write( '  da_km =' + da_km+ '\n')
flog.write( '  strm_cells =' + strm_cells+ '\n')
flog.write( '  valley_block =' + valley_block+ '\n')
flog.write( '  Q100_raster =' + Q100_raster+ '\n')
flog.write( '  start_ARCID =' + str(start_ARCID)+ '\n')
flog.write( '  seg_max =' + str(seg_max)+ '\n')
flog.write( '  mannings_n =' + str(mannings_n)+ '\n')
flog.write( '  alpha =' + str(alpha)+ '\n')
flog.write( '  beta =' + str(beta)+ '\n')
flog.write( '  hill_buff_dist =' + str(hill_buff_dist)+ '\n')
flog.write( '  Q_tol =' + str(Q_tol)+ '\n')
flog.write( '  HSthresh_up =' + str(HSthresh_up)+ '\n')
flog.write( '  HSthresh_low =' + str(HSthresh_low)+ '\n')
flog.write( '  debris_RO =' + str(debris_RO)+ '\n')
flog.write( '  glacial_min_elev =' + str(glacial_min_elev)+ '\n')
flog.write ('OUTPUT from MODEL'+ '\n')

# Additional General settings for model run
arcpy.SnapRaster = inDEM   # Snap all created rasters to inDEM
inBasename = 'S'      # Basename of files (keep below 3 char)
tempCellSize = 10.0
arcpy.CellSize = tempCellSize

# Create Shapefile from inDEM extent
inDEM_sh = userworkspace + '/temp' + '/inDEM_sh' + '.shp'
reclassifyRanges = "0.000000 50000.000000 1"   # Set the reclassify ranges

inDEM_1 = Reclassify(inDEM, "Value", reclassifyRanges, "NODATA")
inDEM_1.save(userworkspace + '/temp' + '/inDEM_1')

arcpy.RasterToPolygon_conversion(inDEM_1, inDEM_sh, "NO_SIMPLIFY", "VALUE")

# Expand spatial reference of the valley_section to that of the DEM (spatialRef also used for cutlines)

### Reset Extent to full Extent of DEM
##arcpy.env.extent = inDEM_sh

dataset = arcpy.Describe(inDEM_sh)
tempExtent = dataset.Extent
##arcpy.Extent = tempExtent

    # ###########################################################################
    # C. Create channel raster and channel bankfull-width shapefile
##try:
try:
    # Calculate bankfull channel width raster
    strm_accum = ExtractByMask(da_km, strm_cells)
    strm_accum.save(userworkspace + '/temp' + '/strm_accum')  # Accumulation of stream cells (in km2)

    strm_power = Power(strm_accum, beta)
    strm_power.save(userworkspace + '/temp' + '/strm_power')  # Intermediate calculation step

    strm_wdth = Times(strm_power, alpha)
    strm_wdth.save(userworkspace + '/strm_wdth')    # Bankfull channel width - raster
  
except:
    arcpy.AddMessage(arcpy.GetMessages(2))
    print '  ERROR: Unsuccessful creating bankfull-width shapefile'
    print arcpy.GetMessages(2)

# #############################################################################
# D. Calculate Slope for entire DEM

all_slp_pct = Slope(inDEM, "PERCENT_RISE", "1")
all_slp_pct.save(userworkspace + '/temp' + '/all_slp_pct')  # slope in pct of entire DEM

all_slp_100 = Raster("temp/all_slp_pct") / 100.0
all_slp_100.save(userworkspace + '/temp' + '/all_slp_100')  # slope in decimal pct of entire DEM

# #############################################################################
# E. Extract individual stream segments from comprehensive stream segment shapefile  

# Define input shapefile and field of interst
inShapeFile = valley_section
inField = "ARCID"                       

spatialRef = arcpy.Describe(all_slp_100).spatialReference

# Create empty parameter dictionaries for later use
s_length_dict = {}         # Empty stream length dictionary
iter_ARCID_dict = {}        # dictionary to match up the loop number (or FID) with the ARCID
slope_class_dict = {}       # dictionary of slope classes for each stream segment

##try:
# Open a cursor on the input shape file attribute table
cursor1 = arcpy.SearchCursor(inShapeFile)
##    arcpy.GetCount(inShapeFile) = seg_max

# Open a file to track the names of the new shape files
##    newShapeList = open(userworkspace + '/temp' + '/Stream_Segment_list.txt','w')
row1 = cursor1.next()
#print 
print "Starting separation of stream segments..."
# print("  Input shape file: " + inShapeFile)
# print("  Field to be used for separating shapes: " + inField)

n = 0
val = 0

while row1:
##    if val+1 > seg_max:
##        print "  EARLY OUT: Met preset stream segment limit"
##        break
    val = row1.getValue(inField)
    slope_class = row1.getValue("GRID_CODE")
    val_s = str("%05d" % (val))
    select_exp = inField + '=' + str(val)
    select_exp = inField + '=' + str(val)
    outShapeFile = userworkspace + '/temp/seg' + '/S_' + inBasename + val_s + '.shp'
##        print("valley_section - New shape file to be created: " + outShapeFile)
    try:
        # Put the feature into a new shape file on it's own
        arcpy.Select_analysis(inShapeFile,outShapeFile,select_exp)

        # Write the new shapefile names and locations to a text file (currently not needed)
##        print "Stream segment -", outShapeFile, 'exported to', userworkspace
##            newShapeList.write(userworkspace + '/temp' + '/' + outShapeFile)
    except:
        print("  ERROR -   Could not create " + userworkspace + '/temp' + '/' + outShapeFile)
        break

    # Create a stream length dictionary (val_s, length in m) for later use
    temp_length = row1.getValue("SLength")
##        feat = row.shape
##        temp_length=feat.length
    s_length_dict[val_s]=temp_length

    # Create a dictionary of iteration_number vs. ARCID of segment (during the editing
    #   process some segments were deleted, thus there are fewer segments than the number
    #   of loops necessary to complete the watershed)
    iter_ARCID_dict[n]=val
##        print n, ",", iter_ARCID_dict[n]

    # Create a dictionary of stream slope classes
    slope_class_dict[val_s]=slope_class
    
    # Read the next record from the search cursor        
    row1 = cursor1.next()
    n += 1
    
# Delete the cursors
del row1
del cursor1

# Close the file listing the new shape files
##    newShapeList.close()

print '  FINISHED separation of', n , 'stream segments'    

##except:
##    arcpy.AddMessage(arcpy.GetMessages(2))
##    print arcpy.GetMessages(2)

# End of valley segment loop
# ***********************************************************************************

# #############################################################################
# F. Build cutlines to separate right and left hillslopes of stream segments

print 'Preparing cut lines for hillslopes of 1st order streams'
# Input/Output files
inputlines = valley_section
textfile = userworkspace + '/temp' + '/Segment_Extension.txt'
outputlines = userworkspace + '/temp' + '/valley_extension' +'.shp'
##ext_distance = 0.0
ext_distance = 20 * hill_buff_dist   # distance to extend stream segments

#print "  Extension distance = " + str(ext_distance)

try:
    os.remove(textfile)
except:
    pass

#Create a text file and write polylines to the first line.
f = open(textfile,'a')
thestring = "ARCID,E1S0,x1,y1,x2,y2\n"
f.writelines(thestring)
f.close()   

# Create search cursor
cursor6 = arcpy.SearchCursor(inputlines)
row6 = cursor6.next()

counter = 0
#start the row iteration
n = 0
val = 0
while row6:
    if val+1 > seg_max:
        print "  EARLY OUT: Met preset stream segment limit"
        break
    val = row6.getValue("ARCID")
##    print '  Working on entry', val
    # Create the geometry object
    feat = row6.Shape
    #get coordinate values as lists
    firstpoint1 = feat.firstPoint
    lastpoint1 = feat.lastPoint
    midpoint1 = feat.centroid
    #split the lists by the blank space between the coordinate pairs
    firstpoint = str(firstpoint1).split(" ")
    lastpoint = str(lastpoint1).split(" ")
    midpoint = str(midpoint1).split(" ")
    #get the x and y values as array positions 0 and 1, and convert them to floating point numbers from the native string literals
    startx = float(firstpoint[0])
    starty = float(firstpoint[1])
    endx = float(lastpoint[0])
    endy = float(lastpoint[1])
    midx = float(midpoint[0])
    midy = float(midpoint[1])

##    print '  Start x,y =', str("%d"%(startx)),',',str("%d"%(starty))
##    print '  Mid x,y =', str("%d"%(midx)),',',str("%d"%(midy))
##    print '  End x,y =', str("%d"%(endx)),',',str("%d"%(endy))

    # Add 'End' line    
    #if the line is horizontal or vertical the slope and negreciprocal will fail so do this instead.
    if endy==midy or endx==midx:
        if endy == midy:
            y1 = endy 
            y2 = endy 
            x1 = endx 
            if endx > midx:
                x2 = endx + ext_distance
            elif endx < midx:
                x2 = endx - ext_distance
        if endx==midx:
            y1 = endy 
            x1 = endx 
            x2 = endx 
            if endy > midy:
                y2 = endy + ext_distance
            elif endy < midy:
                y2 = endy - ext_distance
    
    else:
        #get the slope of the line
        m = ((endy - midy)/(endx - midx))
##        print '  Slope (end) =', m
        if endy > midy and endx > midx:
            y1 = endy
            y2 = endy + (((ext_distance**2.0)/(1 + ((1/m)**2.0)))**(0.5))
            x1 = endx
            x2 = endx + (((ext_distance**2.0)/(1 + (m**2.0)))**(0.5))
        elif endy < midy and endx > midx:
            y1 = endy
            y2 = endy - (((ext_distance**2.0)/(1 + ((1/m)**2.0)))**(0.5))
            x1 = endx
            x2 = endx + (((ext_distance**2.0)/(1 + (m**2.0)))**(0.5))
        elif endy < midy and endx < midx:
            y1 = endy
            y2 = endy - (((ext_distance**2.0)/(1 + ((1/m)**2.0)))**(0.5))
            x1 = endx
            x2 = endx - (((ext_distance**2.0)/(1 + (m**2.0)))**(0.5))
        elif endy > midy and endx < midx:
            y1 = endy
            y2 = endy + (((ext_distance**2.0)/(1 + ((1/m)**2.0)))**(0.5))
            x1 = endx
            x2 = endx - (((ext_distance**2.0)/(1 + (m**2.0)))**(0.5))
        else:
            print "  ERROR - somewhere in end line"
        
    f = open(textfile,'a')
#    thestring = str(val) + " 0\n" + "0 "+ str(x1)+" "+str(y1) + "\n" + "1 " + str(x2) + " " + str(y2) +"\n"
    thestring = str(val)+", 1, " + str(x1)+", "+str(y1)+", " + str(x2) + ", " + str(y2) +"\n"
    f.writelines(thestring)
    f.close()   
    del x1
    del x2
    del y1
    del y2

##    # Need to increase counter to it builds another text file    
##    counter = counter + 1
    
    # Add 'Start' line    
    #if the line is horizontal or vertical the slope and negreciprocal will fail so do this instead.
    if starty==midy or startx==midx:
        if starty == midy:
            y1 = starty 
            y2 = starty 
            x1 = startx 
            if startx > midx:
                x2 = startx + ext_distance
            elif startx < midx:
                x2 = startx - ext_distance

        if startx==midx:
            y1 = starty 
            x1 = startx 
            x2 = startx 
            if starty > midy:
                y2 = starty + ext_distance
            elif starty < midy:
                y2 = starty - ext_distance
    
    else:
        #get the slope of the line
        m = ((starty - midy)/(startx - midx))
##        print '  Slope (start) =', m
        if starty > midy and startx > midx:
            y1 = starty
            y2 = starty + (((ext_distance**2.0)/(1 + ((1/m)**2.0)))**(0.5))
            x1 = startx
            x2 = startx + (((ext_distance**2.0)/(1 + (m**2.0)))**(0.5))
        elif starty < midy and startx > midx:
            y1 = starty
            y2 = starty - (((ext_distance**2.0)/(1 + ((1/m)**2.0)))**(0.5))
            x1 = startx
            x2 = startx + (((ext_distance**2.0)/(1 + (m**2.0)))**(0.5))
        elif starty < midy and startx < midx:
            y1 = starty
            y2 = starty - (((ext_distance**2.0)/(1 + ((1/m)**2.0)))**(0.5))
            x1 = startx
            x2 = startx - (((ext_distance**2.0)/(1 + (m**2.0)))**(0.5))
        elif starty > midy and startx < midx:
            y1 = starty
            y2 = starty + (((ext_distance**2.0)/(1 + ((1/m)**2.0)))**(0.5))
            x1 = startx
            x2 = startx - (((ext_distance**2.0)/(1 + (m**2.0)))**(0.5))
        else:
            print "ERROR - somewhere in start line"
        
    f = open(textfile,'a')
##    thestring = str(val) + " 0\n" + "0 "+ str(x1)+" "+str(y1) + "\n" + "1 " + str(x2) + " " + str(y2) +"\n"
    thestring = str(val)+", 0, " + str(x1)+", "+str(y1)+", " + str(x2) + ", " + str(y2) +"\n"
    f.writelines(thestring)
    f.close()   
    del x1
    del x2
    del y1
    del y2
    
    n += 1
    row6 = cursor6.next()

del row6
del cursor6

#Write extension line points to a line shapefile
arcpy.XYToLine_management(textfile,outputlines,"x1","y1","x2","y2","GEODESIC","ARCID", spatialRef)
# 'spatialRef' as defined above for 'all_slp_100 ~Line 263'

print '  Finished creating shapefile from Extensions.'    

# ################################################################
# Merge extension lines with 1st order stream lines to form cut lines

hill_cut_merge = userworkspace + '/temp'+ '/hill_cut_merge' + '.shp'
hill_cut_final = userworkspace + '/temp'+ '/hill_cut_final' + '.shp'

try:
    arcpy.Merge_management([inputlines, outputlines], hill_cut_merge)
    arcpy.Dissolve_management(hill_cut_merge, hill_cut_final, "ARCID")
    
except:
    arcpy.AddMessage(arcpy.GetMessages(2))
    print arcpy.GetMessages(2)
    
# ###############################################################
# G. Extract individual cutlines from comprehensive cutline shapefile

inShapeFile = hill_cut_final
inField = "ARCID"                       

try:
    # Open a cursor on the input shape file attribute table
    cursor5 = arcpy.SearchCursor(inShapeFile)

    # Open a file to track the names of the new shape files
##    newShapeList = open(userworkspace + '/temp' + '/Stream_Segment_list.txt','w')
    row = cursor5.next()

    print "Starting separation of hillslope cut lines..."
    #print("  Input shape file: " + inShapeFile)
    #print("  Field to be used for separating shapes: " + inField)

    n = 0
    val = 0
    while row:
        if val+1 > seg_max:
            print "  EARLY OUT: Met preset stream segment limit"
            break
        val = row.getValue(inField)
        val_s = str("%05d" % (val))
        select_exp = inField + '=' + str(val)
        outShapeFile = userworkspace + '/temp/seg' + '/CL_' + inBasename + val_s + '.shp'
        # Put the feature into a new shape file on it's own
        try:
            arcpy.Select_analysis(inShapeFile,outShapeFile,select_exp)
        except:
            print("ERROR -   Could not create " + userworkspace + '/temp' + '/' + outShapeFile)
            break

        # Read the next record from the search cursor        
        row = cursor5.next()
        n += 1
        
    # Delete the cursors
    del cursor5
    del row
except:
    arcpy.AddMessage(arcpy.GetMessages(2))
    print arcpy.GetMessages(2)

# End of valley cut line 
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#print
print 'Begin determination the valley bottom for each valley block'

# #################################################################################
# H. Allocate sections of basin to each stream segment using Euclidean Allocation #Moved to valley segment creation

# May need these checks if sequentail folders aren't used
##try:
##    arcpy.delete(valley_bl_sh)
##except:
##    print("Euclidean Allocation " + valley_al_sh + " did not pre-exist.")

##try:
### Initiate parameters   
##inField = "GRIDCODE"
##valley_bl_sh = userworkspace + '/temp' + '/valley_bl_sh' + '.shp'
##valley_block = userworkspace + '/valley_block'+'.shp'
##
### Convert Watersheds from raster to polygon
##arcpy.RasterToPolygon_conversion(val_seg_ras, valley_bl_sh, "NO_SIMPLIFY", "VALUE") 
##arcpy.Dissolve_management(valley_bl_sh, valley_block, inField)
   
##except:
##    arcpy.AddMessage(arcpy.GetMessages(2))
##    print "ERROR: Failed to perform Euclidean Allocation"
##    print arcpy.GetMessages(2)
    
### &&&&&&&&&&&&&&&&&&&&&&&&
### TIME TEST START
##time_test = 'Test1'
##start_time = datetime.now()

# Trim the blocks by a reclassified/converted .shp of the DEM
##    try:
#Iniate shapefile

##arcpy.Clip_analysis(valley_block_sq, inDEM_sh, valley_block)    # Hydro-Geo valley bottom # deleted "+'.shp'"
##arcpy.AddField_management(inDEM_sh, "Shape_area", "DOUBLE")
##arcpy.CalculateField_management (inDEM_sh, "Shape_area", "!shape.area!", "PYTHON_9.3")

##except: #Added by DB 4/14/2014
##    arcpy.AddMessage(arcpy.GetMessages(2))
##    logging.debug ("  ARCID " + val_s + " Error: " + arcpy.GetMessages(2))
    
'''
 //////////////////////////////////////////////////////////////////////////
 Summary of the following steps to determine the valley bottom for
   each valley block
   I. Create individual valley blocks from larger 'valley_block' shapefile
   J. Extract DEM and hillslopes by valley block
   K. Calculate flood depth and extent of Q100+ using Manning's equation
   L. Proportionally expand area greater than Q100 for BiS analysis
   M. Creating lower limit for BiS analysis using multiple of channel bankfull width
   N. Combine lower and upper limits of fluvial valley bottom
   O. Calculate BiS layer for each block
   '''

# #######################################################################
# I. Extract individual valley blocks from comprehensive valley blocks shapefile 

# Define input shapefile and field of interest
inShapeFile = valley_block
inField = "GRIDCODE"             # Set the fields of interest

##try:
cursor2 = arcpy.SearchCursor(inShapeFile)
##newShapeList = open(userworkspace + '/temp' + '/Valley_block_list.txt','w')
row2 = cursor2.next()
#print("  Input shape file: " + inShapeFile)
#print("  Field to be used for separating shapes: " + inField)

step_ARCID = 0
d = 0
if start_ARCID > 1:
    print "  NOTICE - NOT starting at beginning, starting with ARCID", start_ARCID
while iter_ARCID_dict[d] < start_ARCID:
    row2 = cursor2.next()
##    print '    Current ARCID', iter_ARCID_dict[d], 'Start ARCID', start_ARCID
    d += 1
    
m=0     #   Starting value for segment limit
val = 0

midtime = datetime.datetime.now()  # used to note start time of model run        
while row2:
##        print "top loop"
##    try:
    if val +1 > seg_max: # m is the row number of the current segment
        print 
        print "  EARLY OUT: Met preset stream segment limit"
        break

    print '----------------------------'
    
    # Reset Extent to full Extent of DEM
    dataset = arcpy.Describe(inDEM)
    tempExtent = dataset.Extent
    arcpy.Extent = tempExtent
    arcpy.env.extent = inDEM_sh

    # Extract name of block from 'inField'
    val = row2.getValue(inField)\
    #open('last-segment.txt','w+').write(str(int(val + 1)))  #Added by RSAC, incorrect syntax (and I don't know what it does)
    val_s = str("%05d" % (val))
    select_exp = inField + '=' + str(val)

    print "#" + str(m+1) + " Valley Section for ARCID", val_s

    # Define name for new individual valley block (toggle permanent name to keep)
    outShapeFile = userworkspace + '/temp/seg' + '/B_' + inBasename + val_s + '.shp'
##        outShapeFile = userworkspace + '/temp' + '/block' + '.shp'
    stream_segment = userworkspace + '/temp/seg' + '/S_' + inBasename + val_s + '.shp'  # referencing previously established files, can't make intermediate
    block_buff = userworkspace + '/temp' + '/block_buff.shp'
    # next copy the selected shape into a new shape file
##        try:
    arcpy.Select_analysis(inShapeFile,outShapeFile,select_exp)
##        newShapeList.write(userworkspace + '/' + outShapeFile)
    arcpy.Buffer_analysis(outShapeFile, block_buff, "10.0", "FULL", "FLAT", "NONE")
        #print "    Valley Block exported"
##        except:
##            print("  ERROR - COULD NOT CREATE " + outShapeFile)

    # Reset Extent to Extent of outShapeFile
    arcpy.env.extent = block_buff #Added and removed 6/11/2014

##        dataset = arcpy.Describe(outShapeFile)
##        tempExtent = dataset.Extent
##        arcpy.Extent = tempExtent

# #######################################################################
#   J. Extract DEM and hillslopes by valley block

    outDEM = ExtractByMask(inDEM, outShapeFile)
    outDEM.save(userworkspace + '/temp' + '/outDEM')

    # Extract decimal slope by valley block    
    slope_pct_100 = ExtractByMask(all_slp_100, outShapeFile)
    slope_pct_100.save(userworkspace + '/temp' + '/slope_pct_100')
##        slope_pct_100.save(userworkspace + '/temp/seg' + '/SP1_' + inBasename + val_s) ## DB: 6/9/2014 saves unique version

# #######################################################################
#   K. Calculate flood depth and extent of Q100+ using Manning's equation
#       (Used for upper/outer extent of valley bottom)

    print '  Calculating Q100 for stream block...'  

    # Retrieve  stream length from s_length_dict    
    try:
        s_length = float(s_length_dict[val_s])
##            print '    Stream Length =', str(s_length)[:8]
    except:
        print '  ERROR - Could not retreive valley length for ', val_s

    # Pull Q100 value for each segment from Q100_raster
    Q100_stream = ExtractByMask(Q100_raster, stream_segment)
    Q100_stream.save(userworkspace + '/temp' + '/Q100_stream')

    Q100val = arcpy.GetRasterProperties_management(Q100_stream, "Maximum") # RSCAC changed "Mean" to "Maximum"
    Q100 = float(Q100val.getOutput(0))
    
    #print '    Length =', str(s_length)[:8],' Q100 =', str(Q100)[:6]
    if Q100 == 0.0:
        print "    ERROR reading Q100 raster, Q100 set to 20 m3/s"
        Q100 = 20.0        

    # Initialize Values needed for this section
    n = 0                   # Reset counter for number of iterations of Q calc
    q = 0                   # Rest counter of number of iterations with d below flood_min
    Q_calc = 0.0            # Reset Q_calc from previous segment
    Volume = 0.0
    Q_est = Q100 * Q_tol    # Use Q100 * Q_tol for flood estimate
##    flood_depth = 0.108*(Q_est**0.428)      # Starting flood depth in meters (based on regression)
    flood_depth = 0.208*(Q_est**0.428)      # Starting flood depth in meters (based on regression)
                                            #   If starting depths are consistenly too high or low
                                            #   (thus model is taking too many loops to converge)
                                            #   adjust this equation.  Feel free to fit to local data

##    print "    Starting Flood Depth ", str(flood_depth)[:5]
    Q_diff = abs((Q_est - Q_calc) / Q_est)  # % difference btw estimated and calculated Q100

    #Initialize datasets calculated in this step
    file_loc = userworkspace + '/temp/seg' + '/SV_' + val_s + '.txt'

    # Extract stream DEM and calculate slope using 'rise over run'
    s_dem = ExtractByMask(outDEM, stream_segment)
    s_dem.save(userworkspace + '/temp' + '/s_dem')
    elev_min_result = arcpy.GetRasterProperties_management(s_dem, "Minimum")
    elev_max_result = arcpy.GetRasterProperties_management(s_dem, "Maximum")
    elev_min = float(elev_min_result.getOutput(0))
    elev_max = float(elev_max_result.getOutput(0))
    slope = (elev_max - elev_min) / s_length
    
    #print '    elev_min =', str(elev_min)[:6], 'elev_max', str(elev_max)[:6], 'slope =', str(slope)[:6]
        
    # Set slope equal to 0.1% if DEM shows zero slope for segment   
    if slope == 0:
        slope = 0.001
        print "    NOTICE: DEM Slope = 0, therefore assigned slope of 0.1% (0.001)"
    flagForContinue = False

    # Iterate Manning's equation to converge on Q_est
    try:
        while (Q_diff > diff_tol) and q <= iter_max:
##                print "loop here"
    ##        print "Q_diff = " + str(Q_diff)

            # Limit number of iterations to converge on Q
            if n >= 12:
                print "  EARLY OUT:  Met iteration limit"
                break
            elif q > iter_max:
                n = 11  # must be n-1 iteration limit above
                flood_depth = flood_min
                print "  NOTICE: Flood Depth below vertical accuracy limits of DEM."
                #print "  Q100 depth set to " + str(flood_min)

            # Provide escape if flood raster is only a single cell wide (hence all
            #   cells in it have a value of zero, and the volume is zero)
            p = 0
    ##            class zero_volume(Exception):
    ##                pass
            while Volume == 0.0:
                if p == 0:
                    pass
                elif p > 3:
                    print "    ERROR - Unable to find depth that gives non-zero flood volume"
                    flood_depth = flood_min
                    Q_calc = 0.0
                    break  # need to find a way to pass the remainder of Manning's loop
    ##                    raise zero_volume
                else:
                    print "    CAUTION - Envoking Zero Volume Loop, step", p, 'D =', str(flood_depth)[:5]
                    flood_depth = flood_depth * 1.5

                # Build flood raster (Norman's flood simulator)
                arcpy.Extent = "MAXOF"
                flood_raster = CostDistance(stream_segment, slope_pct_100, flood_depth)
                #flood_raster.save(userworkspace + '/temp/seg' '/F_' + inBasename + val_s) # Don't need flood raster saved for each
                flood_raster.save(userworkspace + '/temp' '/Flood_r')

                # Extract flood geometry with Surface Volume
                arcpy.SurfaceVolume_3d(flood_raster, file_loc, "ABOVE", "0")

                # Read output from Surface Area Tool
                FO = open(file_loc)
                SA = csv.DictReader(FO)
        ##        SA = csv.DictReader(open(file_loc)) 
                for row in SA:
                    Area_2D = float(row[' Area_2D'])
                    Area_3D = float(row[' Area_3D'])
                    Volume = float(row[' Volume'])
        ##        print "   ", str("%d"%(Area_2D)), str("%d"%(Area_3D)), str("%d"%(Volume)), file_loc

                # Close SurfaceVolume output file
                FO.close() 
                p +=1

            # Preliminary calculations for Manning's Eq   
            flood_vol = (Area_2D * flood_depth)- Volume 
            xc_area = flood_vol / s_length
            wtd_perimeter = Area_3D / s_length
            hyd_radius = xc_area / wtd_perimeter
            Q100width = Area_2D / s_length
        
            # Manning's equation calculations and comparison with Q_est
            Q_calc = (1 / mannings_n) * (hyd_radius**(0.66666666)) * xc_area * (slope**(0.5))
            Q_diff = abs((Q_est - Q_calc) / Q_est)

            # Adjust flood depth for next iteration   
            if n < 5:
                Divisor = 2
            else:
                Divisor *= 1.4
            Q_compare = ((Q_calc - Q_est) / (Q_calc)) # Used to adust flood_depth in next step
            flood_depth_old = flood_depth   # Store previous flood depth for reporting purposes
            flood_depth = ((1+(Q_compare/-Divisor))*flood_depth_old)     # Bisection type convergance equation  
                
            #print '    Loop',n,': Q=',str("%f"%(Q_calc)),'d=',str("%f"% (flood_depth_old)),"d_new=", str("%f"%(flood_depth))

            # Limit new flood depth to no more than 3x more or less than the previous step
            if (flood_depth / flood_depth_old ) > 3.0:
                flood_depth_temp = flood_depth
                flood_depth = (flood_depth_temp**(0.5))
                #print '      ADJUSTED: New depth of ',str(flood_depth_temp)[:5],'too great, adjusted to',str(flood_depth)[:5]

            elif (flood_depth / flood_depth_old ) < 0.3:
                flood_depth_temp = flood_depth
    ##                flood_depth = (flood_depth_old/3)
                flood_depth = flood_depth_temp**(2.0)
                #print '    ADJUSTED: New depth',str(flood_depth_temp)[:5],'too shallow, adjusted to', str(flood_depth)[:5]
                
            # Increment for next step of convergance loop
            n += 1
            if flood_depth < flood_min:
                q += 1

##    except zero_volume:
##        pass
    except Exception,e:  #RSAC added 'Exception,e'
        print "    ERROR - Could not calculate Mannings Equation"
        print '    Final for',val_s,': Q_calc=', str(Q_calc)[:5], "(", str(Q_est)[:5] ,"), Q_diff=", str(Q_diff)[:6]
        arcpy.AddMessage(arcpy.GetMessages(2))
        print arcpy.GetMessages(2)
        flagForContinue = True  #Added by RSAC
        raise e  #Added by RSAC

    if flood_depth_old < flood_min:
        print "    NOTICE: Q100 convergance on depth below vertical tolerance of DEM. depth(Q100) set to " + str(flood_min)[:3] + "m"
        flood_depth = flood_min
        flood_depth_old = flood_min
        # Build flood raster (Norman's flood simulator)
        arcpy.Delete_management(flood_raster)    
        flood_raster = CostDistance(stream_segment, slope_pct_100, flood_depth)
        ###flood_raster.save(userworkspace + '/temp/seg' '/F_' + inBasename + val_s) #Don't need flood raster saved for each
        flood_raster.save(userworkspace + '/temp' '/Flood_r')

        # Extract flood geometry with Surface Volume
        try:
            os.remove(file_loc)
        except:
            print "    Could not delete ", file_loc
        arcpy.SurfaceVolume_3d(flood_raster, file_loc, "ABOVE", "0")

        # Read output from Surface Area Tool
        FO = open(file_loc)
        SA = csv.DictReader(FO)
##        SA = csv.DictReader(open(file_loc)) 
        for row in SA:
            Area_2D = float(row[' Area_2D'])
            Area_3D = float(row[' Area_3D'])
            Volume = float(row[' Volume'])
##        print "   ", str("%d"%(Area_2D)), str("%d"%(Area_3D)), str("%d"%(Volume)), file_loc

        # Close SurfaceVolume output file
        FO.close() 
            
        # Preliminary calculations for Manning's Eq   
        flood_vol = (Area_2D * flood_depth)- Volume 
        xc_area = flood_vol / s_length
        wtd_perimeter = Area_3D / s_length
        hyd_radius = xc_area / wtd_perimeter
        Q100width = Area_2D / s_length
    
        # Manning's equation calculations and comparison with Q_est
        Q_calc = (1 / mannings_n) * (hyd_radius**(0.66666666)) * xc_area * (slope**(0.5))
            
        #print '    Loop',n,': Q=',str("%f"%(Q_calc)),'d=',str("%f"% (flood_depth))

    # Report final values for 
    print '    Final for',val_s,': Q_calc=', str(Q_calc)[:5], "Q_diff=", str(Q_diff)[:5], "Width =", str(Q100width)[:4]

##    # &&&&&&&&&&&&&&&&&&&&&&&&
##    # TIME TEST START
##    time_test = 'Test1'
##    start_time = datetime.now()
    
    # Create Shapefile for upper limits of valley width (based upon Q100)
    reclassifyRanges = "0.000000 30.000000 1"   # Set the reclassify ranges
    vw_upper = Reclassify(flood_raster, "Value", reclassifyRanges, "NODATA")
    vw_upper.save(userworkspace + '/temp' + '/vw_upper')

    # ####################################################################
    # L. Proportionally expand area greater than Q100 for BiS analysis 

    # Initialize files/values for this step        
    h_mult = 2.0    # Depth multiplier of Q100 depth for BiS analysis
                    #   i.e. if 'h_mult = 2.0', all area within 2x the depth of Q100
                    #   will be condidered for BiS analysis
    UL_Depth = flood_depth * h_mult
    #print '    Upper Limit (UL) Depth = ',str(UL_Depth)[:6]
    
# 7/26/2013 DB: Encountering an error in the creation of UL_Flood_ra for some segments (and apparently some model
# runs.
    try:
        UL_flood_ra = CostDistance(stream_segment, slope_pct_100, UL_Depth)    
        UL_flood_ra.save(userworkspace + '/temp' + '/UL_flood_ra')
        #print '    Flood raster created...'

    except:
        UL_flood_ra = vw_upper
        print '    ERROR in BiS analysis, used Q100(exact) for upper limit instead of 2x' # DB should 

# Reclassify the flood_raster to a single value to define the outer limit of the valley edge
    try:
        reclassifyRanges = "0.000000 30.000000 1"   # Set the reclassify ranges
        UL_flood_1 = Reclassify(UL_flood_ra, "Value", reclassifyRanges, "NODATA")
        UL_flood_1.save(userworkspace + '/temp' + '/UL_flood_1')
        #print '    Single value flood raster created...'
    except:
        print '    ERROR reclassifying UL_flood_ra'
##    # &&&&&&&&&&&&&&&&&&&&&&&&
##    # TIME TEST FINISH
##    finish_time = datetime.now()
##    delta_time = (finish_time)-(start_time)
##    print time_test, 'run time',str(delta_time)[:-7] 

# //////////////////////////////////////////////////////////////////////////

    # #########################################################################
    # M. Creating lower limit for BiS analysis using multiple of channel bankfull width

    print '  Performing BiS analysis...'
    
    #Initialize files     
    vw_ll = userworkspace + '/temp' + '/vw_ll' + '.shp'
    vw_lower_nd = userworkspace + '/temp' + '/vw_lower_nd'
    buff_width = 0.0    # Re-set buffer width to zero
    BF_mult = 1.0   # Multplier for BF width to set lower possible limit of valley edge
                    #   i.e. If 'BF_mult = 2.0' twice the modeled width of the bankful
                    #   channel will be excluded from BiS analysis

    # Set extent for upper and lower limits to the allocation block for this segment  
    dataset = arcpy.Describe(outShapeFile)
    tempExtent = dataset.Extent
    arcpy.Extent = tempExtent
##    print "tempExtent= ", str(tempExtent)

    # Extract and buffer lower limit of stream
    strm_block = ExtractByMask(strm_wdth, outShapeFile)
    strm_block.save(userworkspace + '/temp' + '/strm_block')

    BF_width_result = arcpy.GetRasterProperties_management (strm_block, "Mean")
    BF_width = float(BF_width_result.getOutput(0))

##    # Create a stream length dictionary (val_s, length in m) for later use
##    s_BFwidth_dict[val_s]=BF_width

    # Calculate buffer around stream segments equal to BF_width *     
    buff_width = BF_mult * BF_width
    #print '    Inner buffer width = ' + str(buff_width)[:5]
    arcpy.Buffer_analysis(stream_segment, vw_ll, buff_width, "FULL", "FLAT", "NONE")
    arcpy.FeatureToRaster_conversion(vw_ll, 'ARCID', vw_lower_nd, tempCellSize)
    
    # ##########################################################################3
    # N. Combine lower and upper limits of fluvial valley bottom 

    # Identify cells of lower limit with 'No Data' (lower limit cells will be =0
    vw_lower = IsNull(vw_lower_nd)
    vw_lower.save(userworkspace + '/temp' + '/vw_lower')

    # Add upper limit
    vw_plus = Plus(UL_flood_1, vw_lower)
    vw_plus.save(userworkspace + '/temp' + '/vw_plus')

    # Possible valley  bottom cells will be equal have 'Value = 2'
    #   (1=outside of lower limit + 1=inside of upper)
    vw_final = Con(vw_plus, "1", "", 'Value = 2')   
    vw_final.save(userworkspace + '/temp' + '/vw_final')
##        vw_final.save(userworkspace + '/temp/seg' + '/VW_' + inBasename + val_s) ## DB: 6/9/2014 saves unique version

    # #######################################################################
    # O. Calculate BiS layer for each block
##    try:
    #Initialize files     
    curvature = userworkspace + '/temp' + '/curvature'
##        curvature = userworkspace + '/temp/seg' + '/Crv_' + inBasename + val_s ## DB: 6/9/2014 saves unique version

    # Initialize final outputs for geomorphic, hydrologic, and hydro-geomorphic valley bottoms
    # Toggle to make permanent 
    HG_final = userworkspace + '/temp' + '/HG_final' + '.shp'
    H_final_many = userworkspace + '/temp' + '/H_final_many' + '.shp'
    H_final = userworkspace + '/temp' + '/H_final' + '.shp'
    G_final = userworkspace + '/temp' + '/G_final' + '.shp'
    G_final_temp = userworkspace + '/temp' + '/G_final_temp' + '.shp'
    
    BiS_mult = 2.0      # Mult if you want > or < a single STD above min BiS value !! NOT CURRENTLY USED
    BiS_stat_name = 'Mean' # Can't do 'median' on ArcGIS 9.2 or with floating point DEM's

    #Calculate curvature and extract by possible valley bottom    
    arcpy.Curvature_3d(outDEM, curvature)   # Could use slp-slp times reclassed curv, but that would require more steps
    BiS_surf = ExtractByMask(curvature, vw_final)
    BiS_surf.save(userworkspace + '/temp' + '/BiS_surf')
##        BiS_surf.save(userworkspace + '/temp/seg' + '/BiS_' + inBasename + val_s)## DB: 6/9/2014 saves unique version

    # Determine the minimum values for the BiS_surf to extract
    BiS_Max_result = arcpy.GetRasterProperties_management (BiS_surf, "Maximum")
    BiS_Min_result = arcpy.GetRasterProperties_management (BiS_surf, "Minimum")
    BiS_Max = float(BiS_Max_result.getOutput(0))
    BiS_Min = float(BiS_Min_result.getOutput(0))


##    if (1.2*(BiS_mult * BiS_STD))>(-BiS_Min): # In rare cases when BiS_mult > 1 BiS_UL can be less then BiS_Min and therefore no points are selected
##        BiS_UL = BiS_Mean - (BiS_STD)
##    else:
##        BiS_UL = BiS_Mean - (BiS_mult * BiS_STD)
##        BiS_UL = (( BiS_Max - BiS_Min ) * 0.2) + BiS_Min # 20% greater than BiS_Min DB: justification for this?
    BiS_UL = (( -BiS_Min ) * 0.3) + BiS_Min # 20% greater than BiS_Min DB: justification for this? edited 6/6/2014
    
    BiS_exp = 'Value < ' + str(BiS_UL)
    print '    Bis_UL =', str(BiS_UL)[:5]
##        print '  BiS_exp ', BiS_exp
    BiS_final = Con(BiS_surf, "1", "", BiS_exp)
    BiS_final.save(userworkspace + '/temp' + '/BiS_final')

    # Extract values of flood depth as determined by BiS using Zonal Statistics
    BiS_uDepth = ZonalStatistics(BiS_final, "VALUE", UL_flood_ra, BiS_stat_name, "DATA")
    BiS_uDepth.save(userworkspace + '/temp' + '/BiS_uDepth')
    BiS_stat_result = arcpy.GetRasterProperties_management (BiS_uDepth, "Mean")
    BiS_stat = float(BiS_stat_result.getOutput(0))
    
    # Assign lower limit of BiS above channel to 1.0 meter   
    if BiS_stat <= 1.0:
        BiS_stat = 1.0

    # Flood block to final BiS_stat value, then reclass and convert to .shp
    BiS_btmR = CostDistance(stream_segment, slope_pct_100, BiS_stat)
    BiS_btmR.save(userworkspace + '/temp' + '/BiS_btmR')

    BiS_1 = Reclassify(BiS_btmR, "Value", reclassifyRanges, "DATA")
    BiS_1.save(userworkspace + '/temp' + '/BiS_1')

    inField = "GRIDCODE"
    arcpy.RasterToPolygon_conversion(BiS_1, G_final_temp, "NO_SIMPLIFY")
    arcpy.Dissolve_management(G_final_temp, G_final, inField)    # Keep this Geomorphic valley bottom and compile to a single shapefile

    # Convert Q100 to polygon and clip edges of BiS valley bottom     
    arcpy.RasterToPolygon_conversion(vw_upper, H_final_many, "NO_SIMPLIFY")
    arcpy.Dissolve_management(H_final_many,H_final,"GRIDCODE")
    arcpy.Clip_analysis(G_final, H_final, HG_final)    #Hydro-Geo valley bottom

    # ##################################################################
    # Accurately calculate the three different valley bottom widths (HG, G, and H)
    # NOTE: this section contains two techniques for calculating valley width,
    #   both work, but the 'edge' technique (converting the valley edge to a line
    #   and then sampling a Euclidean distance to the stream raster 
    #   with the lines to get a mean value) works better for sinuous streams
    #   than the 'area' technique (dividing the valley bottom area by the stream
    #   length to get width).  Thus the 'area' technique is commented out.
    
    # Initiate fields for this step
    block_minus = userworkspace + '/temp' + '/block_minus' + '.shp'

    HG_final_line = userworkspace + '/temp' + '/HG_final_line' + '.shp'
    HG_final_sep = userworkspace + '/temp' + '/HG_final_sep' + '.shp'
    HG_final_sides = userworkspace + '/temp' + '/HG_final_sides' + '.shp'
##        HG_final_sides = userworkspace + '/temp/seg' + '/HGs_' + inBasename + val_s + '.shp' ## DB: 6/9/2014 saves unique version

    H_final_line = userworkspace + '/temp' + '/H_final_line' + '.shp'
    H_final_sides = userworkspace + '/temp' + '/H_final_sides' + '.shp'
    H_final_sep = userworkspace + '/temp' + '/H_final_sep' + '.shp'

    G_final_line = userworkspace + '/temp' + '/G_final_line' + '.shp'
    G_final_sides = userworkspace + '/temp' + '/G_final_sides' + '.shp'
    G_final_sep = userworkspace + '/temp' + '/G_final_sep' + '.shp'

    # Back buffer the stream segment block 
    arcpy.Buffer_analysis(outShapeFile, block_minus, "-5 Meters", "" , "FLAT", "NONE")

    # Set extent for upper and lower limits to the allocation block for this segment  
    dataset = arcpy.Describe(outShapeFile)
    tempExtent = dataset.Extent
    arcpy.Extent = tempExtent

    # Create raster of distance from the stream
    DistFromStr = EucDistance(stream_segment) 
##        DistFromStr.save(userworkspace + '/temp' + '/DistFromStr')
    DistFromStr.save(userworkspace + '/temp/seg' + '/DFS_' + inBasename + val_s)

    # Calculate average width for HG_final (HG Valley Bottom) 'edge' technique
    arcpy.PolygonToLine_management(HG_final, HG_final_line)
    arcpy.Clip_analysis(HG_final_line, block_minus, HG_final_sep)
    arcpy.Dissolve_management(HG_final_sep, HG_final_sides, "FID")
    HG_width_r = ZonalStatistics(HG_final_sides, "ID", DistFromStr, "MEAN", "DATA") #RSAC changed "FID" to "ID"
    HG_width_r.save(userworkspace + '/temp' + '/HG_width_r')

    HG_width1_result = arcpy.GetRasterProperties_management (HG_width_r, "Mean")
    HG_width1 = float(HG_width1_result.getOutput(0))
    HG_width = 2.0*(HG_width1)  # Must multiply by 2.0 as raster distance is from only one side to the stream
    #print '    HG_width =', str(HG_width)[0:6]

    # 'area' technique for HG_width
##    cursor3 = arcpy.SearchCursor(HG_final)
##    row3 = cursor3.next()
##    feat = row3.shape
##    HG_area = feat.area
##    HG_width = HG_area / float(s_length_dict[val_s])
##    del cursor3
##    del row3

##    print 'Area, Width =',temp_area, HG_width
##    print '***Width_edges, Width_area =', HG_width1, HG_width 

    # Calculate Valley/BF width ratio
    V_BF_ratio = HG_width / BF_width

    # Calculate average width for Q100 (Hydro Valley Bottom)'edge' technique
    arcpy.PolygonToLine_management(H_final, H_final_line)
    arcpy.Clip_analysis(H_final_line, block_minus, H_final_sep)
    arcpy.Dissolve_management(H_final_sep, H_final_sides, "FID")
    H_width_r = ZonalStatistics(H_final_sides, "ID", DistFromStr, "MEAN", "DATA") #RSAC changed "FID" to "ID"
    H_width_r.save(userworkspace + '/temp' + '/H_width_r')

    #Q100_width = 2.0*(arcpy.GetRasterProperties_management (H_width_r, "Mean"))  # Must multiply by 2.0 as raster distance is from only one side to the stream

    Q100_width1_result = arcpy.GetRasterProperties_management (H_width_r, "Mean")
    Q100_width1 = float(Q100_width1_result.getOutput(0))
    Q100_width = 2.0*(Q100_width1)  # Must multiply by 2.0 as raster distance is from only one side to the stream
    #print '    Q100 width =', str(Q100_width1)[:5], '(with multiplier)'
     
    # 'area' technique for Q100_width
##    cursor7 = arcpy.SearchCursor(H_final)
##    row7 = cursor7.next()
##    feat = row7.shape
##    Q100_area = feat.area
##    Q100_width = Q100_area / float(s_length_dict[val_s])
##    del cursor7
##    del row7
##    print 'Area, Width =',temp_area, Q100_width

    # Calculate average width for BiS (Geomorphic Valley Bottom)'edge' technique
    arcpy.PolygonToLine_management(G_final, G_final_line)
    arcpy.Clip_analysis(G_final_line, block_minus, G_final_sep)
    arcpy.Dissolve_management(G_final_sep, G_final_sides, "FID")
    G_width_r = ZonalStatistics(G_final_sides, "ID", DistFromStr, "MEAN", "DATA") #RSAC changed "FID" to "ID"
    G_width_r.save(userworkspace + '/temp' + '/G_width_r')
    
##    BiS_width = 2.0*(arcpy.GetRasterProperties_management (G_width_r, "Mean"))  # Must multiply by 2.0 as raster distance is from only one side to the stream

    BiS_width1_result = arcpy.GetRasterProperties_management (G_width_r, "Mean")
    BiS_width1 = float(BiS_width1_result.getOutput(0))
    BiS_width = 2.0*(BiS_width1)  # Must multiply by 2.0 as raster distance is from only one side to the stream
    #print '    BiS width =', str(BiS_width)[:5]

    # 'area' technique for BiS_width
##    cursor8 = arcpy.SearchCursor(G_final)
##    row8 = cursor8.next()
##    feat = row8.shape
##    BiS_area = feat.area
##    BiS_width = BiS_area / float(s_length_dict[val_s])
##    del cursor8
##    del row8
##    print 'Area, Width =',temp_area, BiS_width
    
##    except:
##        print "    ERROR - Could not calculate BiS for this block"
##        arcpy.AddMessage(arcpy.GetMessages(2))
##        print arcpy.GetMessages(2)


    print "    FINISHED - BiS analysis"

## //////////////////////////////////////////////////////////////////////////
## Summary of the following steps to classify the hillslope steepness for
##     each valley block
##   P.	Buffer valley bottoms to determine areas of hillslopes to analyze
##   Q.	Cut hillslopes by extension lines
##   R.	Separate split hillslopes into right and left
##   S.	Analyze hillslope 'steepness'

    # ########################################################################
    # P.	Buffer valley bottoms to determine areas of hillslopes to analyze

    print '  Classifying hill slopes...'    

    try:
        hill_buff = userworkspace + '/temp' + '/hill_buff' + '.shp'
        hill_buff_d = userworkspace + '/temp' + '/hill_buff_d' + '.shp'
        hill_noVB = userworkspace + '/temp' + '/hill_noVB' + '.shp'
        hill_clip = userworkspace + '/temp' + '/hill_clip' + '.shp'

        hill_buff_str = str(hill_buff_dist) + " Meters"
        print '    Hill_buff_str =', hill_buff_str
    
        arcpy.Buffer_analysis(HG_final, hill_buff_d, hill_buff_str, "FULL", "FLAT", "NONE") # leave inField for actual model
        arcpy.Dissolve_management(hill_buff_d, hill_buff, "GRIDCODE")
        arcpy.Erase_analysis(hill_buff,HG_final,hill_noVB)
        arcpy.Clip_analysis(hill_noVB, outShapeFile, hill_clip)
    except:
        arcpy.AddMessage(arcpy.GetMessages(2))
        print arcpy.GetMessages(2)

    # ############################################################################
    # Q.	Cut hillslopes by extension lines
    # The hillslopes are cut by extensions of the channel segments.  This is done
    #   particularly for 1st order streams, but also for higher order streams that
    #   may have connected hillslopes at the upstream or downstream ends.
    #print "    Building 1st order cutlines..."

##    try:
    cutline = userworkspace + '/temp/seg' + '/CL_' + inBasename + val_s + '.shp'  # Referencing previously established file
    cut_clip = userworkspace + '/temp'+ '/cut_clip' + '.shp'
    cut_buff = userworkspace + '/temp'+ '/cut_buff' + '.shp'
    hill_union = userworkspace + '/temp'+ '/hill_union' + '.shp'
    hill_erase = userworkspace + '/temp'+ '/hill_erase' + '.shp'
    
##        hill_split = userworkspace + '/temp'+ '/hill_split' + '.shp'
    hill_split = userworkspace + '/temp/seg' + '/HS_' + inBasename + val_s + '.shp'
    
##        if os.path.exists(cutline):   
    arcpy.Clip_analysis(cutline, outShapeFile, cut_clip)
    arcpy.Buffer_analysis(cut_clip, cut_buff, "1 Meter", "FULL", "FLAT", "NONE") #This is the source of many of the errors
    #union_str = hill_clip + ";" + cut_buff 
    union_str = '\"%(hill)s\"; \"%(cut)s\"' % {"hill":hill_clip,"cut":cut_buff}# RSAC replaced above line with this one        
    arcpy.Union_analysis(union_str, hill_union,"ONLY_FID")
    arcpy.MultipartToSinglepart_management(hill_union, hill_erase)
    arcpy.Erase_analysis(hill_erase,cut_buff,hill_split)
##        else:
##            arcpy.MultipartToSinglepart_management(hill_clip, hill_split)
##            print'    NOTE: No existing cutline'

    # #################################################################
    # R. Separate split hillslopes into right and left
    hill_right = userworkspace + '/temp'+ '/hill_right' + '.shp'
##        hill_right = userworkspace + '/temp/seg' + '/HR_' + inBasename + val_s + '.shp'## DB: 6/9/2014 saves unique version
    hill_left = userworkspace + '/temp'+ '/hill_left' + '.shp'
##        hill_left = userworkspace + '/temp/seg' + '/HL_' + inBasename + val_s + '.shp' ## DB: 6/9/2014 saves unique version
    
    # Add area to each entry
    arcpy.AddField_management(hill_split, "Poly_Area", "Float")
    arcpy.CalculateField_management(hill_split, "POLY_AREA", "float(!SHAPE.AREA!)", "PYTHON")

    cursor10 = arcpy.SearchCursor(hill_split,"","","","POLY_AREA D")
    row10 = cursor10.next()
    n_cat = arcpy.GetCount_management(hill_split)
    print '   ', n_cat, 'hillslope area(s), classifying steepness...'

    # Largest hillslope are is always labeled "Right"
    if n_cat == 0:
        pass
    elif n_cat >= 1:
        FID_num = row10.getValue("FID")   
        R_L_exp = "FID=" + str(FID_num)
        arcpy.Select_analysis(hill_split,hill_right,R_L_exp)
        row10 = cursor10.next()

    # If more than one hillslope are present the second is labeled "Left"      
    if n_cat <= 1: #RSAC added pass, changed ">" to "<="
        pass
    elif n_cat >1: # added by RSAC
        try:
            FID_num = row10.getValue("FID")   
            R_L_exp = "FID=" + str(FID_num)
            arcpy.Select_analysis(hill_split,hill_left,R_L_exp)
        except:
            pass
    else:
        print "    CAUTION - Only one hillslope present for segment", val_s
    del cursor10
    del row10

##    except:
##        arcpy.AddMessage(arcpy.GetMessages(2))
##        print arcpy.GetMessages(2)

    # ########################################################################
    # S. Analyze hillslope 'steepness'
    arcpy.env.overwriteOutput = True
    
    s = 0
    r = 0
    while s <=1:
##          ---------------------
##            if s == 0:
##                hill_side = hill_right
##                side = '"R"'
##                hill_sl_cat = userworkspace + '/temp' + '/hill_R_cat' + '.shp'
##            elif n_cat > 1:
##                hill_side = hill_left
##                side = '"L"'
##                hill_sl_cat = userworkspace + '/temp' + '/hill_L_cat' + '.shp'
##          ---------------------
        if n_cat == 0:  # RSAC replaced block ---- above with below
            print "    CAUTION - No hillslopes present"
            side = "NA"
        elif s == 0:
##            print "    Classifying hillslope steepness RIGHT..."
            hill_side = hill_right
            side = '"R"'
            hill_sl_cat = userworkspace + '/temp' + '/hill_R_cat' + '.shp'
        elif n_cat > 1:
##            print "    Classifying hillslope steepness LEFT..."
            hill_side = hill_left
            side = '"L"'
            hill_sl_cat = userworkspace + '/temp' + '/hill_L_cat' + '.shp'
##          ---------------------
        else:
            break

        hill_sl_sh = userworkspace + '/temp' + '/hill_sl_sh' + '.shp'
        hill_cat = userworkspace + '/temp' + '/hill_cat' + '.shp'
        
        HS_up_plus = HSthresh_up + 0.0001
        HS_low_plus = HSthresh_low + 0.0001

        if n_cat ==0:
            report_stat = 99
        else:
            hill_sl_pct = ExtractByMask(slope_pct_100, hill_side)
            hill_sl_pct.save(userworkspace + '/temp' + '/hill_sl_pct')

            reclassifxyRanges = "0.00 "+str(HSthresh_low)+" 1; "+str(HS_low_plus)+" "+str(HSthresh_up)+" 2; "+str(HS_up_plus)+" 1000.0 3"
            hill_sl_cl = Reclassify(hill_sl_pct, "Value", reclassifxyRanges, "NODATA")
            hill_sl_cl.save(userworkspace + '/temp'+ '/hill_sl_cl')
            arcpy.RasterToPolygon_conversion(hill_sl_cl, hill_sl_sh, "NO_SIMPLIFY")
            arcpy.Dissolve_management(hill_sl_sh, hill_sl_cat, "GRIDCODE")

            # Associate hillslopes with current segment
            arcpy.AddField_management(hill_sl_cat, "ARCID", "Short")
            arcpy.CalculateField_management (hill_sl_cat, "ARCID", val)

            # Assign "R" or "L"         
            arcpy.AddField_management(hill_sl_cat, "R_OR_L", "TEXT","", "", "5")
            arcpy.CalculateField_management (hill_sl_cat, "R_OR_L", side)

            # Add area to each entry
            arcpy.AddField_management(hill_sl_cat, "Poly_Area", "Float")
            arcpy.CalculateField_management(hill_sl_cat, "POLY_AREA", "float(!SHAPE.AREA!)", "PYTHON")

            # Zero out stats from previous segment   
            area1 = 0.0
            area2 = 0.0
            area3 = 0.0
            area_tot = 0.0
            prop1 = 0.0
            prop2 = 0.0
            prop3 = 0.0
            report_stat = 0

        # Calculate proportions for each category and determine a the final slope category
            cursor4 = arcpy.SearchCursor(hill_sl_cat)
            row4 = cursor4.next()

            while row4:
                P_area = row4.getValue("Poly_Area")
                grid = row4.getValue("GRIDCODE")
                if grid ==1:
                    area1 = P_area
                elif grid ==2:
                    area2 = P_area
                elif grid == 3:
                    area3 = P_area
                else:
                    print "    ERROR - Incorrect slope class reported"
                row4 = cursor4.next()
            del cursor4
            del row4
            del P_area

            area_tot = area1 + area2 + area3
            prop1 = area1 / area_tot
            prop2 = area2 / area_tot
            prop3 = area3 / area_tot
##            print '    area1',side, str(area1)[:9], str(prop1)[:5]
##            print '    area2',side, str(area2)[:9], str(prop2)[:5]
##            print '    area3',side, str(area3)[:9], str(prop3)[:5]
##            print '    Total',side, str(area_tot)[:9], str(prop1+prop2+prop3)[:5]

    # Limits of proportion of hillslopes in each category for classification
            thresh1 = 0.75      # % < HSthresh_low required for  'low' hillslope 
            thresh3 = 0.25      # % > HSthresh_up required for 'high' hillslope
            if prop1 >= thresh1:
                report_stat = 1
            elif prop3 >= thresh3:
                report_stat = 3
            else:
                report_stat = 2
    ##        print '    report_stat=', report_stat

            arcpy.AddField_management(hill_sl_cat, "HS_Cat", "SHORT")
            arcpy.CalculateField_management (hill_sl_cat, "HS_Cat", report_stat)

            # Pass values to report_stat (1st time through is R, 2nd time L)
            if s == 0:
                cat_right = report_stat
                
            elif n_cat > 1:
                cat_left = report_stat                    
            else:
                cat_left = 0  #  Value of "cat_left = 0" indicates that only a single hillslope was present
            s += 1
    #print "      Cat_right", cat_right, "and Cat_left", cat_left

##            except:
##                print "    ERROR - Did not complete hill slope steepness"
##                arcpy.AddMessage(arcpy.GetMessages(2))
##                print arcpy.GetMessages(2)
##                s += 1
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
    # #############################################################################
    # T. Calculate hillslope coupling statistic

    sideR=1
    sideL=1

    if cat_right == 1:   # If cat_right is 2 or 3, sideR must be 1
        sideR = 0

    if cat_left == 1:  # If cat_left is 2 or 3, sideL must be 1
        sideL = 0

    if cat_left == 0:  # If cat_left is 0, only cat_right exists so **assume** both sides are the same
        sideL = sideR

    coup_stat = ((sideR + sideL) * debris_RO) / (HG_width - BF_width)
    #print "    Hillslope coupling statistic =", str(coup_stat)[:4]
##    print "    Finished hill slope classification."

    # #############################################################################
    # U. Calculate Valley Classification

    valley_class = 0
    valley_name = '"UNC"'
    
    slope_class = slope_class_dict[val_s]    
    if slope_class == 3:
        if coup_stat >= 0.75:
            valley_class = 1        # High energy coupled
            valley_name = '"HEC"'
        else:
            valley_class = 2        # High energy open
            valley_name = '"HEO"'
    elif slope_class == 2:
        if V_BF_ratio <= 7:
            valley_class = 3        # Moderate energy confined
            valley_name = '"MEC"'
        else:
            valley_class = 4        # Medium energy open
            valley_name = '"MEO"'
    elif slope_class == 1:
        valley_class = 8            # Low energy floodplain
        valley_name = '"LEF"'
    else:
        print "    ERROR - Incorrect slope classification!"

    if V_BF_ratio < 3 and cat_right == 3 and cat_left ==3:
        valley_class = 5            # Canyon
        valley_name = '"CAN"'
    elif V_BF_ratio >= 3 and cat_right == 3 and cat_left ==3:
        valley_class = 6            # Gorge
        valley_name = '"GOR"'
    elif slope_class < 3 and HG_width >= 100 and elev_min > glacial_min_elev:
        valley_class = 7
        valley_name = '"GLA"'

    print '  FINAL HGVC for ARCID ',val_s,'= #' + str(valley_class), "(", str(valley_name[1:-1]),")"
    try:
        flog.write ("  ARCID " + val_s + " Complete, HGVC = #" + str(valley_class)+ "(" + str(valley_name[1:-1]) + ")"+ '\n')
    except:
        print '    ERROR: writing to log file'
        
    # ########################################################################
    # V. Build output shapefiles of all valley bottoms

    print '  Appending results to previous:'
##    print'     Building output Shapefiles for Valley Bottoms'

    # Reset Extent to full Extent of DEM
    dataset = arcpy.Describe(inDEM)
    tempExtent = dataset.Extent
    arcpy.Extent = tempExtent

    # ________________________________________________________________
    # HYDRO-GEO (BiS clippled by Q100): Remove unused fields and add field corresponding to segment number
    FieldName = "ARCID"
    FieldName2 = "S_Length"
    FieldName3 = "BF_Width"
    FieldName6 = "Grad_Mean"
    FieldName7 = "HG_V_Width"
    FieldName8 = "V_BF_Ratio"
    FieldName4 = "Coupling"
    FieldName5 = "HS_L_Cat"
    FieldName9 = "HS_R_Cat"
    FieldName10 = "Min_Elev"
    FieldName11 = "Val_Class"
    FieldName12 = "Val_Cl_Abv"
    FieldName1 = "Slope_Cl"
    try:
        arcpy.AddField_management(HG_final, FieldName, "Short")
        arcpy.CalculateField_management (HG_final, FieldName, val)

        arcpy.DeleteField_management(HG_final, "ID; GRIDCODE")
        
        arcpy.AddField_management(HG_final, FieldName2, "Long", "8")
        arcpy.CalculateField_management (HG_final, FieldName2, float(s_length_dict[val_s]))
        
        arcpy.AddField_management(HG_final, FieldName3, "Float")
        arcpy.CalculateField_management (HG_final, FieldName3, float(BF_width))

        arcpy.AddField_management(HG_final, FieldName6, "Float")
        arcpy.CalculateField_management (HG_final, FieldName6, float(slope))

        arcpy.AddField_management(HG_final, FieldName7, "Float")
        arcpy.CalculateField_management (HG_final, FieldName7, float(HG_width))

        arcpy.AddField_management(HG_final, FieldName8, "Float")
        arcpy.CalculateField_management (HG_final, FieldName8, float(V_BF_ratio))        

        arcpy.AddField_management(HG_final, FieldName4, "Float")
        arcpy.CalculateField_management (HG_final, FieldName4, float(coup_stat))

        arcpy.AddField_management(HG_final, FieldName5, "Short")
        arcpy.CalculateField_management (HG_final, FieldName5, cat_left)

        arcpy.AddField_management(HG_final, FieldName9, "Short")
        arcpy.CalculateField_management (HG_final, FieldName9, cat_right)

        arcpy.AddField_management(HG_final, FieldName10, "Short")
        arcpy.CalculateField_management (HG_final, FieldName10, elev_min)

        arcpy.AddField_management(HG_final, FieldName1, "Short")
        arcpy.CalculateField_management (HG_final, FieldName1, slope_class)

        arcpy.AddField_management(HG_final, FieldName11, "Short")
        arcpy.CalculateField_management (HG_final, FieldName11, valley_class)

        arcpy.AddField_management(HG_final, FieldName12, "TEXT", "", "", "5")
        arcpy.CalculateField_management (HG_final, FieldName12, valley_name)

    except:
        print "  ERROR updating fields - HYDRO-GEO", arcpy.GetMessages(2)

    # ________________________________________________________________
    # HYDRO (Q100): Remove unused fields and add field corresponding to segment number
    FieldName = "ARCID"
    FieldName2 = "S_Length"
    FieldName3 = "BF_Width"
    FieldName4 = "Grad_Mean"
    FieldName5 = "Depth_Q100"
    FieldName6 = "Width_Q100"
    FieldName1 = "Q100_calc"
    FieldName7 = "XC_Area"
    FieldName8 = "Hyd_Radius"
    
    try:
        arcpy.AddField_management(H_final, FieldName, "short", "6")
        arcpy.CalculateField_management (H_final, FieldName, val)

        arcpy.DeleteField_management(H_final, "ID; GRIDCODE")
        
        arcpy.AddField_management(H_final, FieldName2, "long")
        arcpy.CalculateField_management (H_final, FieldName2, float(s_length_dict[val_s]))
        
        arcpy.AddField_management(H_final, FieldName3, "Float")
        arcpy.CalculateField_management (H_final, FieldName3, float(BF_width))

        arcpy.AddField_management(H_final, FieldName4, "Float")
        arcpy.CalculateField_management (H_final, FieldName4, float(slope))

        arcpy.AddField_management(H_final, FieldName7, "Float")
        arcpy.CalculateField_management (H_final, FieldName7, float(xc_area))

        arcpy.AddField_management(H_final, FieldName8, "Float")
        arcpy.CalculateField_management (H_final, FieldName8, float(hyd_radius))        

        arcpy.AddField_management(H_final, FieldName1, "Float")
        arcpy.CalculateField_management (H_final, FieldName1, float(Q_calc))

        arcpy.AddField_management(H_final, FieldName5, "Float")
        arcpy.CalculateField_management (H_final, FieldName5, float(flood_depth_old))

        arcpy.AddField_management(H_final, FieldName6, "Float")
        arcpy.CalculateField_management (H_final, FieldName6, float(Q100_width))
    except:
        print "  ERROR updating fields - HYDRO", arcpy.GetMessages(2)

    # ________________________________________________________________
    # GEOMORPHOLOGIC (BiS): Remove unused fields and add field corresponding to segment number
    FieldName = "ARCID"
    FieldName2 = "S_Length"
    FieldName3 = "Depth_BiS"
    FieldName4 = "Width_BiS"
    
##    try:
    arcpy.AddField_management(G_final, FieldName, "short", "6")
    arcpy.CalculateField_management (G_final, FieldName, val)

    arcpy.DeleteField_management(G_final, "ID; GRIDCODE")
    
    arcpy.AddField_management(G_final, FieldName2, "LONG", "8")
    arcpy.CalculateField_management (G_final, FieldName2, float(s_length_dict[val_s]))

    arcpy.AddField_management(G_final, FieldName3, "Float")
    arcpy.CalculateField_management (G_final, FieldName3, float(BiS_stat))

    arcpy.AddField_management(G_final, FieldName4, "Float")
    arcpy.CalculateField_management (G_final, FieldName4, float(BiS_width))

    
##    except:
##        print "ERROR updating fields - HYDRO-GEO", arcpy.GetMessages(2)

    # ________________________________________________________________
    # Append the Hydro-Geomorphic valley bottoms into a single .shp
    VB_HydGeo = userworkspace + '/ValleyBottom_HydroGeo' + '.shp'
    try:
        arcpy.Append_management(HG_final, VB_HydGeo, "NO_TEST")
    except:
        arcpy.CopyFeatures_management(HG_final, VB_HydGeo) 
       
    # ________________________________________________________________
    # Append the Hydrologic (Q100) valley bottoms into a single .shp
    VB_Hyd = userworkspace + '/ValleyBottom_Hydro_Q100' + '.shp'
    try:
        arcpy.Append_management(H_final, VB_Hyd, "NO_TEST")
    except:
        arcpy.CopyFeatures_management(H_final, VB_Hyd) 

    # ________________________________________________________________
    # Append the Geomorphic (BiS) valley bottoms into a single .shp
    VB_Geo = userworkspace + '/ValleyBottom_Geo_BiS' + '.shp'
    try:
        arcpy.Append_management(G_final, VB_Geo, "NO_TEST")
    except:
        arcpy.CopyFeatures_management(G_final, VB_Geo) 
    print '    Finished appending to Valley Bottoms'

    # ________________________________________________________________
    # Append the hillslope classifications into a single .shp
    HS_cat = userworkspace + '/Hillslope_categories' + '.shp'
    # Left First
    hill_sl_cat = userworkspace + '/temp' + '/hill_L_cat' + '.shp'
    try:
        arcpy.Append_management(hill_sl_cat, HS_cat, "NO_TEST")
    except:
        arcpy.CopyFeatures_management(hill_sl_cat, HS_cat)
    # Right Second
    hill_sl_cat = userworkspace + '/temp' + '/hill_R_cat' + '.shp'
    try:
        arcpy.Append_management(hill_sl_cat, HS_cat, "NO_TEST")
    except:
        arcpy.CopyFeatures_management(hill_sl_cat, HS_cat)

    print '    Finished appending to Hillslopes'   

##    except: #release with BIG try:except:
##        msgs = arcpy.GetMessages(0)
##        flog.write("  ARCID " + val_s +  " Error: " + msgs + '\n')
##        failed_seg = userworkspace + '/temp/seg' + '/S_' + inBasename + val_s + '.shp'
##        failed_seg_lib = userworkspace + '/Failed_Segs' + '.shp'
##        failed_bl = userworkspace + '/temp/seg' + '/B_' + inBasename + val_s + '.shp'
##        failed_bl_lib = userworkspace + '/Failed_Blks' + '.shp'
##        try:
##            arcpy.Append_management(failed_seg, failed_seg_lib, "NO_TEST")
##            arcpy.Append_management(failed_bl, failed_bl_lib, "NO_TEST")
##        except:
##            arcpy.CopyFeatures_management(failed_seg, failed_seg_lib)
##            arcpy.CopyFeatures_management(failed_bl, failed_bl_lib)
##
##        pass  #Added by RSAC, presumably to pass a failed segment 
##    
##    # Delete un-needed shapefiles referenced from earlier in code
##    arcpy.delete(stream_segment)
##    arcpy.delete(cutline)

# Read the next record from the search cursor     
    m += 1
    row2 = cursor2.next()

# Delete the cursors
del cursor2
del row2
# Close the file listing the new shape files
##newShapeList.close()
#print

# End of valley block loop (it's a long one!)
# ***********************************************************************************

print '__________________________________________________________________ '

# ####################################################################################
# Delete all contents of 'Temp' directory

##TempFolder = userworkspace + '/temp'
##
##try:
##    #
##    class no_delete(Exception):
##        pass
##
##    if delete_temp == "NO":
##        raise no_delete
##
##    print "  Deleting temporary files"
##    # Delete all Shapefiles in TempFolder
##    try:
##        arcpy.Workspace = TempFolder
##        shapes = arcpy.ListFeatureClasses()
##        shp = shapes.next()
##        while shp:
##            arcpy.delete_management(shp)
##            shp = shapes.next()
##        del shp, shapes
##    except:
##        print '  ERROR - Could not delete \temp Shapefiles'
##        arcpy.AddMessage(arcpy.GetMessages(2))
##        print arcpy.GetMessages(2)
##        if shp:
##            del shp
##        if shapes:
##            del shapes
##
##    # Delete all Rasters in TempFolder
##    try:
##        arcpy.Workspace = TempFolder
##        rasters = arcpy.ListRasters("", "GRID")
##        rast = rasters.next()
##        while rast:
##            arcpy.delete(rast)
##            rast = rasters.next()
##        del rast, rasters
##    except:
##        print "  ERROR - Could not delete \temp Rasters"
##        arcpy.AddMessage(arcpy.GetMessages(2))
##        if rast:
##            del rast
##        if rasters:
##            del rasters        
##
##    ### Delete Arc Info folder
##    ##try:
##    ##    os.remove(tempfolder, '/info/arc.dir')
##    ##    os.rmdir(tempfolder, '/info')
##    ##except:
##    ##    print 'ERROR - Could not delete txt files'
##    ##    arcpy.AddMessage(arcpy.GetMessages(2))
##    ##    print arcpy.GetMessages(2)
##        
##    # Delete remaining text files
##    for the_file in os.listdir(TempFolder):
##        file_path = TempFolder + '/' + the_file
##    ##    print file_path
##        try:
##            if not the_file == 'info':
##                os.remove(file_path)
##        except:
##            print
##            print '  ERROR - Could not delete txt files'
##            arcpy.AddMessage(arcpy.GetMessages(2))
##            print arcpy.GetMessages(2)
##
##except no_delete:
##    print "  No deletion selected - DONE."
##    
##    ### Delete \Temp directory # DB: Can't currently delete because I can't get rid of temp\info\arc.dir
##    ##try:
##    ##    os.rmdir(userworkspace + '\temp')
##    ##except:
##    ##    print 'ERROR - Could not delete temp directory'
##    ##    arcpy.AddMessage(arcpy.GetMessages(2))
##    ##    print arcpy.GetMessages(2)

# ##############################################################################################
#Calculate Time Elapsed in Model

try:
    three = 2+1
finally:
    nowtime = datetime.datetime.now()
    diff = nowtime - thentime
    diff2 = nowtime - midtime
    print 'Model run ',str(diff)[:-7], 'for', m, 'seg, (', int(diff.seconds / m), 'sec per segment)'
    print 'After pre-process ',str(diff2)[:-7], '(', int(diff2.seconds / m), 'sec per segment)'
    flog.write( 'COMPLETE: Model run time '+str(diff)[:-7]+ ' for '+ str(m) + ' segments = ('+ str(int(diff.seconds / m))+ ' sec per segment)'+ '\n')
    flog.write( 'After pre-process '+str(diff2)[:-7]+ '('+ str(int(diff2.seconds / m))+ 'sec per segment)')
# Close Error Log
flog.close()

print '__________________________________________________________________ '
del arcpy

import sys
sys.exit(0)
    
##except:
##    arcpy.AddMessage(arcpy.GetMessages(2))
##    print arcpy.GetMessages(2)
##del arcpy




