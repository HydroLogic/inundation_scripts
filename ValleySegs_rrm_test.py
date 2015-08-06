
#Import standard library modules
import sys, os, csv, string, time
from datetime import datetime, date, time

thentime = datetime.now()  # used to note start time of model run
print 'Start Time is:'; thentime

# Import arcpy module
import arcpy   # New for ArcGIS v.10
from arcpy import env
from arcpy.sa import *

print '  Set up environment...'
# Set environment settings
arcpy.ResetEnvironments()

# Check out any necessary licenses
arcpy.CheckOutExtension("Spatial")

# Set OverWriteOutput to 1 to copy over existing outputs, 0 to through exceptions for existing output
arcpy.env.OverWriteOutput = True

# File names and locations
root    = "rmorrison"
folder  = "C:/GIS"
userworkfolder = folder + '/' + root 
##dem_     = userworkfolder + '/' + root + "_dem.img"
dem = 'C:/GIS/srout/srout_dem.img' # Note slashed have to be changed from \ to /
##fdem = 'C:/GIS/nrout/A008/nrout_fdem' 
##fdir = 'C:/GIS/nrout/A008/temp/nrout_fdir' 
##facc = 'C:/GIS/nrout/A008/temp/nrout_facc' 
##strm = 'C:/GIS/nrout/A008/nrout_strm'
##strm_slp = 
##strm_link_img
##link_slp_raw 

# Create workspace folder (Iteratively numbered)
filenum = 1
dir_exists = False 
while not dir_exists:
    if filenum > 500:
        break
    try:
        userworkspace = userworkfolder + '/' + "A" + str(filenum).zfill(3)
        os.mkdir(userworkspace) # Create the workspace
        dir_exists = True
    except:
        filenum += 1
##        print "  Working Folder =" + str(filenum).zfill(3)
        pass
try:
    os.mkdir(userworkspace + '/temp') # Create the /TEMP folder
except:
    print '  ERROR - establishing temp directory (may already exist)'
    arcpy.AddMessage(arcpy.GetMessages(2))
    pass

arcpy.env.workspace = userworkspace

print '  Workspace is set to:', str(env.workspace)

# Set output and temp names and file locations
fdem_    = userworkspace + '/' + root + "_fdem"
fdir_    = userworkspace + '/temp/' + root + "_fdir"
facc_    = userworkspace + '/temp/' + root + "_facc"
strm_    = userworkspace + '/' + root + "_strm"
da_km_   = userworkspace + '/' + root + "_da_km"
da_mi_   = userworkspace + '/temp/' + root + "_da_mi"
segs_   = userworkspace + '/' + root + "_segs" + ".shp"
blks_   = userworkspace + '/' + root + "_blks" + ".shp"
Minimum_Mapping_Unit__cells_ = "\"COUNT\" > 30"
DA_Threshold_Eq = "VALUE > 30000" 

print 'dem =', dem
print 'fdem =', fdem_
print 'fdir =',fdir_
print 'facc =',facc_
print 'strm =',strm_
print da_km_  
print da_mi_   
print segs_   
print blks_

# Set Geoprocessing environments
arcpy.env.extent = dem
env.workspace = userworkfolder

# ############################################################################
# This section of code creates the necessary input files for the HGVC script

# Process: Fill
print ' Fill'
fdem = Fill(dem, "")
fdem.save(fdem_)  # filled DEM

# Process: Flow Direction
print ' Dir'
fdir = FlowDirection(fdem, "NORMAL")
fdir.save(fdir_)

# Process: Flow Accumulation
print ' Acc'
facc = FlowAccumulation(fdir, "", "FLOAT")
facc.save(facc_)

# Process: Con
print ' DA_Threshold'
strm = Con(facc, "1", "", DA_Threshold_Eq)
strm.save(strm_)

# Process: Divide
print ' Acc km2'
da_km = Float(Raster(facc_) / 10000.0)
da_km.save(da_km_)

# Process: Divide (2)
print ' Acc mi'
da_mi = Float(Raster(facc_) / 25899.8811)
da_mi.save(da_mi_)

# Process: Slope
print ' Slope'
tempEnvironment0 = arcpy.env.mask
arcpy.env.mask = strm_
strm_slp = Slope(fdem, "PERCENT_RISE", "1")
strm_slp.save(userworkspace + '/temp' + '/strm_slp')    # Channel slope
##strm_slp_ = '"%s"' % strm_slp

# Process: Stream Link - Create individual stream links separated by nodes @ junctions
print ' StreamLink'
strm_link_img = StreamLink(strm, fdir)
strm_link_img.save(userworkspace + '/temp' + '/strm_link_img')
               
# Process: Zonal Statistics
print ' zonal stat'
link_slp_raw = ZonalStatistics(strm_link_img, "Value", strm_slp, "MEAN", "DATA")
link_slp_raw.save(userworkspace + '/temp' + '/link_slp_raw')
##link_slp_raw_ =  "'" + link_slp_raw + "'" 

# Process: Raster Calculator - Calc average of local slope and link slope to smooth out local variation
print ' Raster1'
strm_slp_   = userworkspace + '/temp/' + 'strm_slp'
link_slp_raw_   = userworkspace + '/temp/' + 'link_slp_raw'

strm_slp_mean = Float((Raster(strm_slp_) + Raster(link_slp_raw_)) /2.0)
strm_slp_mean.save(userworkspace + '/temp' + '/strm_slp_mean')

# Process: Zonal Statistics (2)
print ' zonal stat2'
seg_slp_mean = ZonalStatistics(strm_link_img, "Value", strm_slp_mean, "MEAN", "DATA")
seg_slp_mean.save(userworkspace + '/temp' + '/seg_slp_mean')

# Process: Reclassify 
seg_slp_cls = Reclassify(seg_slp_mean, "Value", "0 0.10000000000000001 1;0.10000000000000001 3 2;3 10000 3", "DATA")
seg_slp_cls.save(userworkspace + '/temp' + '/seg_slp_cls')
seg_slp_cls.save(userworkspace + '/temp' + '/seg_slp_cls')

# Process: Region Group - Group segments by value
print ' RegionGroup2'
val_segs_r = RegionGroup(seg_slp_cls, "EIGHT", "WITHIN", "ADD_LINK", "")
val_segs_r.save(userworkspace + '/temp' + '/val_segs_r')

# Process: Raster to Polyline - Create polyline of stream segments
print ' Raster to Polyline'
val_segs_shpA = (userworkspace + '/temp' + '/val_segs_shpA.shp')
arcpy.RasterToPolyline_conversion(val_segs_r, val_segs_shpA, "ZERO", "20", "SIMPLIFY", "LINK")

# Copy so I have a record of the original Raster to Polyline
val_segs_shp = (segs_)
arcpy.CopyFeatures_management(val_segs_shpA, val_segs_shp)

# Process: Surface Length
print ' Surface Length'
arcpy.AddField_management(val_segs_shp, "SLength", "FLOAT")
arcpy.CalculateField_management (val_segs_shp, "SLength", "!shape.length@meters!", "PYTHON_9.3")

# Process: Remove all short segments
print ' Delete short segments'
try:
    arcpy.Delete_management("val_segs_tbl") # Delete table if it currently exists
except:
    print '  Note: val_segs_tbl does not exist yet'
arcpy.MakeFeatureLayer_management(val_segs_shp, "val_segs_tbl")
arcpy.SelectLayerByAttribute_management ("val_segs_tbl", "NEW_SELECTION", "\"SLength\" <= 50.0")
arcpy.DeleteFeatures_management("val_segs_tbl")
#val_segs_shp2.save(userworkspace + '/val_segs_tbl')

# Convert final valley segments back to raster for watershed delineation
print ' Polyline to Raster'
val_seg_ras = (userworkspace + '/temp' + '/val_segs_ras')
arcpy.PolylineToRaster_conversion(val_segs_shp, "ARCID", val_seg_ras,"", "", 10.0)

# Process: Watersheds
print ' Watersheds'
arcpy.env.mask = fdir_
val_seg_ws = Watershed(fdir, val_seg_ras)
val_seg_ws.save(userworkspace + '/temp' +  '/val_seg_ws')

# Initiate parameters   
inField = "GRIDCODE"
valley_bl_sh = userworkspace + '/temp' +  '/valley_bl_sh' + '.shp'
valley_block = blks_

# Convert Watersheds from raster to polygon
arcpy.RasterToPolygon_conversion(val_seg_ws, valley_bl_sh, "NO_SIMPLIFY", "VALUE") 
arcpy.Dissolve_management(valley_bl_sh, valley_block, inField)



print '("_______________________________________________________________")'

    
# ##############################################################################################
#Calculate Time Elapsed in Model
nowtime = datetime.now()
diff = nowtime - thentime
print 'COMPLETE: Model run time',str(diff)[:-7]

try:
    three = 2+1
finally:
    nowtime = datetime.now()
    diff = nowtime - thentime
    print 'Model run ',str(diff)[:-7]
    
sys.exit(0)