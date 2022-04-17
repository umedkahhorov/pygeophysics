######################## SEGY ##################################
import segyio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shutil import copyfile
import utm

##1
def read_segy(filename,coors=True,static=False):
    """
    read segy file using segyio, coors=True returns XYZ coors
    filename    segy file (test.segy)
    

    Output:
    data        data array
    hdr         dataframe pandas, sou XYZ, rec XYZ, offsets
                if coors=True 

    """
    with segyio.open(filename,'r',ignore_geometry=True) as segyfile:

        segyfile.mmap()
        if coors==True:
            # Extract XY coordinates for all traces
            sourceX = segyfile.attributes(segyio.TraceField.SourceX)[:]
            sourceY = segyfile.attributes(segyio.TraceField.SourceY)[:]
            sou_depth = segyfile.attributes(segyio.TraceField.SourceSurfaceElevation)[:]

            groupX = segyfile.attributes(segyio.TraceField.GroupX)[:]
            groupY = segyfile.attributes(segyio.TraceField.GroupY)[:]
            rec_depth = segyfile.attributes(segyio.TraceField.ReceiverGroupElevation)[:]
            offsets = segyfile.attributes(segyio.TraceField.offset)[:]
            # read static 
            if static==True:
                tot_stat = segyfile.attributes(segyio.TraceField.TotalStaticApplied)[:]

            sc_2 = segyfile.attributes(segyio.TraceField.ElevationScalar)[:]
            sc_2 = sc_2[0]
            sc_1 = segyfile.attributes(segyio.TraceField.SourceGroupScalar)[:]
            sc_1 = sc_1[0]

            if sc_1 < 0:
                sc_1 = abs(sc_1)

            if sc_2 < 0:
                sc_2 = abs(sc_2)
            
            if sc_1==0 or sc_2==0:
                sc_1 = 1
                sc_2 = 1

            hdr = [sourceX/sc_1,sourceY/sc_1,sou_depth/sc_2,groupX/sc_1,groupY/sc_1,rec_depth/sc_2,offsets]
            
            hdr = pd.DataFrame(data=np.array(hdr).T,columns=['sourceX','sourceY','sou_depth','groupX','groupY','rec_depth',
                                                       'offsets'])

        data = segyfile.trace.raw[:]
        data = np.squeeze(data)

        if coors==True:
            if static==True:
                hdr['tot_stat'] = tot_stat
            return data,hdr
        else:
            return data

##2
def recip_segy(filename):
    """
    read segy file using segyio and switch XYZ coors (Reciprocity)

    Output:
    segy file       copied original and changed source and receivers coordinates
    
    Example:
    recip_segy('Data/24_data_wind.sgy') --> Data/24_data_wind_recip.sgy
    """ 
    
    # make e new filename
    _ = filename.split('.')
    new_filenname = _[0] + '_recip.' + _[-1]
    
    # copy file to store results
    copyfile(filename,new_filenname)

    with segyio.open(new_filenname,'r+',ignore_geometry=True) as segyfile:

        segyfile.mmap()

        # Extract XY coordinates for all traces
        sourceX = segyfile.attributes(segyio.TraceField.SourceX)[:]
        sourceY = segyfile.attributes(segyio.TraceField.SourceY)[:]
        sou_depth = segyfile.attributes(segyio.TraceField.SourceSurfaceElevation)[:]

        groupX = segyfile.attributes(segyio.TraceField.GroupX)[:]
        groupY = segyfile.attributes(segyio.TraceField.GroupY)[:]
        rec_depth = segyfile.attributes(segyio.TraceField.ReceiverGroupElevation)[:]
        # write coordinates back
        
        for i in range(len(sourceX)):
            segyfile.header[i] = {segyio.TraceField.SourceX:(groupX[i])}
            segyfile.header[i] = {segyio.TraceField.SourceY:(groupY[i])}
            segyfile.header[i] = {segyio.TraceField.SourceSurfaceElevation:(rec_depth[i])}

            segyfile.header[i] = {segyio.TraceField.GroupX:(sourceX[i])}
            segyfile.header[i] = {segyio.TraceField.GroupY:(sourceY[i])}
            segyfile.header[i] = {segyio.TraceField.ReceiverGroupElevation:(sou_depth[i])}
##3
def deg_to_utm(filename,copy=False):
    """Read sgy file using segyio and change coordinates
    kr0917_A3obs_s001-1.sgy - Linename_SourceNumber-Component.sgy 

    """
    latlon = lambda latitude_y,longitude_x: utm.from_latlon(latitude_y,longitude_x)
    
    shot_index = int(filename.split('_')[-1].split('-')[0].replace('s',""))
    if copy:
        _ = filename.split('.')
        new_filename = _[0] + '_deg.' + _[-1]
        #copy file to store results
        copyfile(filename,new_filename)

    with segyio.open(filename,'r+',ignore_geometry=True) as segyfile:
        segyfile.mmap()
        # Extract source x,y
        sourceX = segyfile.attributes(segyio.TraceField.SourceX)[:]
        sourceY = segyfile.attributes(segyio.TraceField.SourceY)[:]
        # Extract receivers x,y
        groupX = segyfile.attributes(segyio.TraceField.GroupX)[:]
        groupY = segyfile.attributes(segyio.TraceField.GroupY)[:]
        # Convert deg. to lat/long
        sourceX = sourceX/3600/1000
        sourceY = sourceY/3600/1000
        groupX  = groupX/3600/1000
        groupY  = groupY/3600/1000
        # Convert to UTM
        s_utm = latlon(sourceY,sourceX)
        r_utm = latlon(groupY,groupX)
        sx    = np.around(s_utm[0][:],2)*100
        sy    = np.around(s_utm[1][:],2)*100
        rx    = np.around(r_utm[0][:],2)*100
        ry    = np.around(r_utm[1][:],2)*100
        sx = sx.astype(int)
        sy = sy.astype(int)
        rx = rx.astype(int)
        ry = ry.astype(int)

        traces = list(range(1,len(sx)+1))
        # Write back chanched headers
        segyfile.header = {segyio.TraceField.ShotPoint:shot_index}
        segyfile.header = {segyio.TraceField.SourceGroupScalar:-100}
        
        for i in range(len(sourceX)):
            segyfile.header[i] = {segyio.TraceField.SourceX:(sx[i])}
            segyfile.header[i] = {segyio.TraceField.SourceY:(sy[i])}
            segyfile.header[i] = {segyio.TraceField.GroupX:(rx[i])}
            segyfile.header[i] = {segyio.TraceField.GroupY:(ry[i])}
            segyfile.header[i] = {segyio.TraceField.TRACE_SEQUENCE_LINE:(traces[i])}
        print ('Shot', shot_index, 'deg_to_utm done')
##4
def rotate(filename,radian=-0.16574, sx_min =-470591.6737, sy_mean= 4315246.5733,copy=False):
    """
    Rotate xy after copy segy files. Rotate xy along x -axis, eliminate y-axis 3D --> 2D
    kr0917_A3obs_s001-1.sgy - Linename_SourceNumber-Component.sgy
    copy=True saves original file with "orig"
   vp_grad_muted.rss """
    def rotate_xy(xy,radian):
        c, s = np.cos(radian), np.sin(radian)
        R = np.array([[c, s], [-s, c]])
        xyhat = np.matmul(R,xy)
        x = xyhat.T[:,0]
        y = xyhat.T[:,1]
        return x,y
    
    shot_index = int(filename.split('_')[-1].split('-')[0].replace('s',""))
    
    if copy:
        _ = filename.split('.')
        filenname = _[0] + '_rot.'+ _[-1]

    with segyio.open(filename,'r+',ignore_geometry=True) as segyfile:
        segyfile.mmap()
        # Extract XY coordinates for all traces
        sourceX = segyfile.attributes(segyio.TraceField.SourceX)[:]
        sourceY = segyfile.attributes(segyio.TraceField.SourceY)[:]
        groupX = segyfile.attributes(segyio.TraceField.GroupX)[:]
        groupY = segyfile.attributes(segyio.TraceField.GroupY)[:]
        
        # scales
        sc_coors = segyfile.attributes(segyio.TraceField.SourceGroupScalar)[:]
        sc_coors = abs(sc_coors)
        
        xy_s = [sourceX/sc_coors,sourceY/sc_coors]
        sx,sy = rotate_xy(xy_s,radian)
        xy_r = [groupX/sc_coors,groupY/sc_coors]
        rx,ry = rotate_xy(xy_r,radian)
        
        sx = np.round((sx - sx_min),2)*sc_coors
        rx = np.round((rx - sx_min),2)*sc_coors
        
        sy = np.round((sy - sy_mean),2)*sc_coors
        ry = np.round((ry - sy_mean),2)*sc_coors
        
        sx = sx.astype(int)
        sy = sy.astype(int)
        rx = rx.astype(int)
        ry = ry.astype(int)

        for i in range(len(sx)):
            segyfile.header[i] = {segyio.TraceField.SourceX:(sx[i])}
            segyfile.header[i] = {segyio.TraceField.SourceY:(sy[i])}
            segyfile.header[i] = {segyio.TraceField.GroupX:(rx[i])}
            segyfile.header[i] = {segyio.TraceField.GroupY:(ry[i])}
        print ('Shot', shot_index, 'rotated')
##4 
def write_segy(weigths,filename,copy=False):
    """
    write segy file not chaning headers
    segy data --> replace with new data  
    """
    if copy:
        _ = filename.split('.')
        new_filename = _[0] + 'weight.' + _[-1]
        #copy file to store results
        copyfile(filename,new_filename)
        filename = new_filename

    with segyio.open(filename,'r+',ignore_geometry=True) as segyfile:
        segyfile.mmap()
        # write data
        segyfile.trace.raw[:] = weigths
##5
def read_model(filename):
    """
    read models and ignore headers
    
    Output:
    data - traces, times
    """
    with segyio.open(filename,ignore_geometry=True) as segyfile:
        segyfile.mmap()
        data = segyfile.trace.raw[:]
    data = np.squeeze(data)
    return data
