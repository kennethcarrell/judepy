# -*- coding: utf-8 -*-
"""
rev 2 Dec20-2023   add Dons idea to subtract jupiter from earlier in night
lots of work on improvinng centroid of jupiter by using interpolation between pixels

rev1 12-2-23 clean up variables and change to float arithmetic in jupiter centroid
Jup2StarsBFNov.py   Nov 30 2023  by Bill Fisher
This has extensive addins to Jup3StarsH1Dv1.py  Sept 16, 2023  
to 1) centroid Jupiter 
   2) remove glow of Jupiter   using  1) subtraction of two regions on either side of star 2 
                               or     2) one region 180 degree from star two   
                                      both relative jupiter center 
                               then a method is used to remove the remaing background gredient from the star2 ROI
                               
   The code also has a bit of code so a jupiter image can be added to a non conjunction night data set  for testing 
   


Jup3StarsH1Dv1.py  Sept 16, 2023

Opens a series of FITS files, subtracts a Masterdark file, divides by a 
MasterFlat file, then calculates centroids using two methods. 

An Inner radius defines the centroid measurement area, and a set of middle
and outer radii defines the background annulus region.

Centroids are calculated first using Howell's method (similar to MaxIm DL), 
which is a center-of-mass centroid. It first subtracts the mean background 
level, then subtracts 2.874*background RMS. The center-of-mass is calculated 
from this clipped signal. It is fairly independent of measurement radii.

The program also calculates the centroids using photutils.centroids routines
of 1D gauss (fit to the projected values on x and y). This requires an odd value
for the BOX DIAMETER, nominally = 11. The stars should fit into this box.
The centroids are averaged between Howell's C-of-Mass and the 1D value.

This version requires an initial file Jup3StarJ.txt to indicate FITS file dimension,
approximate centroid positions from the first image in the series, number of 
image files in the series, and the BOX dimension for the centroids. 
If Jupiter is not in this image series, then simply use one of the 
star locations instead of Jupiter, and ignore those centroids.

The location of the input parameter file Jup3StarJ.txt is in Line 59.
File series is specified in Line 60. The images should be named consecutively.
If the first file number is not 1, change firstfile in Line 61.
If the image files do not end in .FITS, change that in Line 62.
Verify the number of digits in the file numbering is correct in line 63

Requires dark.fits (which already has nominal 100 counts subtracted from it,
to prevent zeroes in the image files.) Location is in Line 64.

Requires flat.fits, normalized to approximately 1.0, to improve pixel-pixel gain.
If a flat image is not available, comment out that line in the code, or create
an artificial image with every pixel = 1.000. Location is in Line 65.
The final output file name is specified in Line 66.

If a star is saturated, a warning will be printed in the output Console. 
The saturation value is in Line 68.
Progress is reported in the output Console for every nth file, where n is found
in Line 69. This can be any number, but is meant to avaoid hiding the warnings.

The ring diameters for the estimated positions is pretty broad, and are listed
in Lines 70-72 and 76-78. These can probably stay fixed.
The ring diameters for the final centroids are found in Lines 73-75 and 79-81.

"""
#### Set up matplotlib and use a nicer set of plot parameters
from pyexpat.model import XML_CTYPE_MIXED
import numpy as np
from math import sqrt,exp    #10-30-23BF needed for jup shine correction
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.io import fits
from photutils.centroids import (centroid_1dg,centroid_com, centroid_sources)
#from photutils.detection import IRAFStarFinder


#some helper functions

def jupmodel(radius): 
    jupRadius=50
    #jupPolRadius=48  # for later possible ellipse model
    #jupAngle=92
    dark=1200
    pixAdu=0
    if(radius < 50):pixAdu=65500
    else:pixAdu= dark+ 65500*exp(-(radius-jupRadius)/4.5)
    if( pixAdu>65500):pixAdu=65500
    return(int(pixAdu))


def centroid_JupbyEdge(data, startx, starty):  # for later possible expansion to ellipse model  , angJup, sizejup):
    """
    Calculate the center of jupiter like    
    Parameters
    ----------
    data : `~numpy.ndarray`
        The input n-dimensional array. The image should be a
        background-subtracted cutout image containing a single source.
    Returns
    -------
    centroid : `~numpy.ndarray`
        The coordinates of the centroid in pixel order (e.g., ``(x, y)``      
    """
    # preserve input data - which should be a small cutout image
    data = data.copy()
    #print("on entry    at int x,y= ",startx,starty)
    # get a rough estimate
        
    pixN=starty+250.
    pixS=starty-250.
    pixW=startx+250.
    pixE=startx-250.
    
    pixNa=starty+250
    pixSa=starty-250
    pixWa=startx+250
    pixEa=startx-250

    i=0
    bDone= 0
    while (bDone <1):
        if(data[starty+i][startx]< 18000 ):
            pixNa =starty+i
            rv0=float(data[starty+i][startx])
            rvm=float(data[starty+i-1][startx])
            pixN= float(starty+i) + (18000.-rv0)/(rvm-rv0)
            bDone=1
        i=i+1
        if i>200: bDone=1
    i=0
    bDone= 0
    while (bDone < 1):
        if(data[starty-i][startx]< 18000 ):
            rv0=float(data[starty-i][startx])
            rvp=float(data[starty-i+1][startx])
            pixS= float(starty-i) - (18000.-rv0)/(rvp-rv0)
            pixSa=starty-i
            bDone=1
        i=i+1
        if i>200: bDone=1        
    i=0
    bDone= 0
    while( bDone <1):
        if(data[starty][startx+i]< 18000 ):
            pixWa=startx+i
            rv0=float(data[starty][startx+i])
            rvm=float(data[starty][startx+i-1])
            pixW= float(startx+i) + (18000.-rv0)/(rvm-rv0)
            bDone=1
        i=i+1
        if i>200: bDone=1
    i=0
    bDone= 0
    while (bDone <1):
        if(data[starty][startx-i]< 18000 ):
            rv0=float(data[starty][startx-i])
            rvp=float(data[starty][startx-i+1])
            pixE= float(startx-i) - (18000.-rv0)/(rvp-rv0)
            pixEa=startx-i
            bDone=1
        i=i+1
        if i>200: bDone=1    
    starty= int((pixNa+pixSa)/2)  
    startx= int((pixWa+pixEa)/2)  
    rstarty= (pixN+pixS)/2.  
    rstartx= (pixW+pixE)/2. 

    # 
    #print("after step 1    at int x,y= ",startx,starty)
    ##  60 pix away should be clear of sat jupiter  
    #estval= float(jupmodel(60))
    #valYp= float(data[starty + 60][startx])
    #valYn= float(data[starty - 60][startx])
    #valXp= float(data[starty][startx + 60])
    #valXn= float(data[starty][startx - 60])
      
    #ychange= 2*(valYp-valYn)/estval    
    #xchange= 2*(valXp-valXn)/estval   
    #newY= startx+ ychange
    #newX= startx+ xchange
    #print("after step 2    at x,y= ",newX,newY)      
    #startx=int(newX)
    #starty=int(newY)
    startx=int(rstartx+.5)
    starty=int(rstarty+.5)
    #move  a little closer and redo with last estimate
    ##  53 pix away should be clear of sat jupiter 
    estval= float(jupmodel(53))
    valYp= float(data[starty + 53][startx])
    valYn= float(data[starty - 53][startx])
    valXp= float(data[starty][startx + 53])
    valXn= float(data[starty][startx - 53])
       
    ychange= 2*(valYp-valYn)/estval    
    xchange= 2*(valXp-valXn)/estval
       
    besty= float(starty) + ychange
    bestx= float(startx) + xchange  
    
    #besty= rstarty
    #bestx= rstartx
    #print("final step 3    at x,y= ",bestx,besty)
    return np.array([bestx,besty])




def subtractGradient(image, XcEst,YcEst, jBacsize):   #removes gradient for a region jbacsize around XcEst,YcEst in image
    halfjBacsize=int(jBacsize/2)
    #now compute the remaining slope in the background and subtract it 
    sumN=0   #for sum of north row
    sumE=0   #for sum of east row
    sumS=0   #for sum of south row
    sumW=0   #for sum of West row
    newXE= int(-halfjBacsize+XcEst) #west,east
    newXW= int(halfjBacsize+XcEst-1) #west,east
    for j in np.arange( 0 ,jBacsize,1):
        newY=int(j+YcEst-halfjBacsize) #west
        sumW= sumW + image[newY][newXW]
        sumE= sumE + image[newY][newXE]
        newYS= int(-halfjBacsize+YcEst) #west,east
        newYN= int(halfjBacsize+YcEst-1) #west,east   
    for i in np.arange( 0 ,jBacsize,1):
        newX=int(i+XcEst-halfjBacsize) #west
        sumN= sumN + image[newYN][newX]
        sumS= sumS + image[newYS][newX]   
    aveN=sumN/jBacsize
    aveE=sumE/jBacsize
    aveS=sumS/jBacsize
    aveW=sumW/jBacsize
    slopeY= (aveN-aveS)/jBacsize
    slopeX= (aveW-aveE)/jBacsize
    for j in np.arange( 0 ,jBacsize,1): 
        newY=int(j+YcEst-halfjBacsize)
        for i in np.arange( 0 ,jBacsize,1):
            newX=int(i-halfjBacsize+XcEst)
                 
            image[newY][newX] = image[newY][newX] +  (jBacsize-i)*slopeX +(jBacsize-j)*slopeY
    k=-1
    if k==0 :   # -1 ignores print  
        hduj = fits.PrimaryHDU(np.uint16(image))
        hdulj = fits.HDUList([hduj])
        hdulj.writeto("C:\\jupiterd\\JupTest12-6gradcorInFunc.fits",overwrite=True)
    return()


def roughCentroidFwhm(data, XcEst,YcEst, roisize):
    halfroisize = int(roisize/2)
    largePixVal=0
    bestX=XcEst
    bestY= YcEst
    #find brightest pixel
    for j in np.arange( 0 ,roisize,1): 
        newY=int(j+YcEst-halfroisize)
        for i in np.arange( 0 ,roisize,1):
            newX=int(i-halfroisize+XcEst)
            if(data[newY][newX] > largePixVal ):
                largePixVal=data[newY][newX]
                bestX=newX
                bestY= newY
    #now find halfmax pixels
    pY=bestY
    nY=bestY
    pX=bestX
    nX=bestX
    curY=bestY
    curX=bestX
    i=0
    half=np.uint16(largePixVal/2)
    bDone= 0
    while (bDone <1):
        if(data[curY+i][curX]< half ):
            pY=curY+i
            bDone=1
        i=i+1
        if i>(halfroisize+2) : bDone=1
    i=0
    bDone= 0
    while (bDone < 1):
        if(data[curY-i][curX]< half ):
            nY=curY-i
            bDone=1
        i=i+1
        if i> (halfroisize+2) : bDone=1        
    i=0
    while( bDone <1):
        if(data[curY][curX+i]< half ):
            pX=curX+i
            bDone=1
        i=i+1
        if i>(halfroisize+2): bDone=1
    i=0
    while (bDone <1):
        if(data[curY][curX-i]< half ):
            nX=curX-i
            bDone=1
        i=i+1
        if i>(halfroisize+2): bDone=1  
    bestX= (pX -nX)/2.0 + curX 
    bestY= (pY -nY)/2.0 + curY
    
    fwhm= ((pX-nX) + (pY-nY))/2.
    
    return np.array([bestX,bestY,fwhm])   

def interpolatePixValue(data, xpos, ypos):
    nx= int(xpos-.000001)
    ny= int(ypos-.000001)
    delx = xpos-float(nx)
    dely = ypos-float(ny)
    value= data[ny][nx]*(1.0-delx)*(1-dely)  + data[ny+1][nx]*(1.-delx)*(dely) + data[ny][nx+1]*(delx)*(1.-dely) + data[ny+1][nx+1]*(delx)*(dely)
    return value    

def align2Jup(baseImage, subtractImage, basejX,basejY,offsetX,offsetY,floor):
    boxsize =140
    halfboxsize=70
    testBox= np.zeros((boxsize,boxsize))
    nbasejX=int(basejX+.5)
    nbasejY=int(basejY+.5)

    bestXoffset=offsetX
    bestYoffect=offsetY
    bestavepix=65000.
    rsum=0.0 
    xinc=0.0
    yinc=0.0   
  
    #while nDone<30:
    for m in np.arange(0,7,1):
        yinc= 0.1*float(m) - 0.4 
        for j in np.arange( 0 ,boxsize,1):
                yimagepix= nbasejY- halfboxsize +j 
                ryjup=basejY+ float(j-halfboxsize) +  offsetY  + yinc              
                for i in np.arange(0 ,boxsize,1):
                    ximagepix= nbasejX- halfboxsize +i                 
                    rxjup=basejX+ + float(i-halfboxsize) + offsetX  + xinc  
                    jupval=interpolatePixValue(subtractImage, rxjup, ryjup)
                    testBox[j][i] = baseImage[yimagepix][ximagepix] -jupval + 200 #+ darkAve - jupval
                    if testBox[j][i] >200 :rsum=rsum+ testBox[j][i]
                    else :testBox[j][i]=199 
                #endi
            #endj
        avepix=rsum/(float(boxsize**2))                             
        if(avepix < bestavepix):
            bestavepix=avepix           
            best_xinc=xinc
            best_yinc=yinc
    yinc=best_yinc       
    rsum=0.0
    for l in np.arange(0,9,1):
        xinc= 0.1*float(l) - 0.5
        rsum=0.0
        for j in np.arange( 0 ,boxsize,1):
            yimagepix= nbasejY- halfboxsize +j 
            ryjup=basejY+ float(j-halfboxsize) +  offsetY  + yinc              
            for i in np.arange(0 ,boxsize,1):
                ximagepix= nbasejX- halfboxsize +i                 
                rxjup=basejX+ + float(i-halfboxsize) + offsetX  + xinc  
                jupval=interpolatePixValue(subtractImage, rxjup, ryjup)
                testBox[j][i] = baseImage[yimagepix][ximagepix] -jupval + 200 #+ darkAve - jupval
                if testBox[j][i] >200 :rsum=rsum+ testBox[j][i]
                else :testBox[j][i]=199 
             #endi
        #endj
        avepix=rsum/(float(boxsize**2))
                              
        if(avepix < bestavepix):
                bestavepix=avepix
                #bestXoffset= offsetX+ xinc
                #bestYoffset= offsetY+ yinc
                best_xinc=xinc
                best_yinc=yinc
    #endl 
    #       
    #endm
    bestX= offsetX +best_xinc
    bestY= offsetY +best_yinc
    k=0
    if k== -1 :   # -1 ignores print
        for j in np.arange( 0 ,boxsize,1):
                yimagepix= nbasejY- halfboxsize +j 
                ryjup=basejY+ float(j-halfboxsize) +  offsetY  + best_yinc              
                for i in np.arange(0 ,boxsize,1):
                    ximagepix= nbasejX- halfboxsize +i                 
                    rxjup=basejX+ + float(i-halfboxsize) + offsetX  + best_xinc  
                    jupval=interpolatePixValue(subtractImage, rxjup, ryjup)
                    testBox[j][i] = baseImage[yimagepix][ximagepix] -jupval + 200 #+ darkAve - jupval
                    if testBox[j][i] >200 :rsum=rsum+ testBox[j][i]
                    else :testBox[j][i]=199 
        hduje = fits.PrimaryHDU(np.uint16(testBox))
        hdulje = fits.HDUList([hduje])
        hdulje.writeto("C:\\jupiterd\\testbox.fits",overwrite=True)
        
    if k==-1:
        baseImage[nbasejY][nbasejX]=5000   # mark center pix
        hduje = fits.PrimaryHDU(np.uint16(baseImage))
        hdulje = fits.HDUList([hduje])
        hdulje.writeto("C:\\jupiterd\\baseJup.fits",overwrite=True)
        
        subtractImage[nbasejY+ int(offsetY)][nbasejX+int(offsetX)]=10000   # mark center pix
        hduje = fits.PrimaryHDU(np.uint16(subtractImage))
        hdulje = fits.HDUList([hduje])
        hdulje.writeto("C:\\jupiterd\\Jupsubtract.fits",overwrite=True)

    return np.array([bestX,bestY , avepix])


#end helper functions



#begin main program 
################  file names listed here
#parameter_file = "C:/Users/sales/Desktop/Jup3StarJ.txt"
parameter_file = "C:/jupiterd/Jup3StarJB.txt"
#file_series ="F:\\jupDef\\Don\\TestImages\\ED127 StarJupiter 64arcsec\\L-"
#file_series ="F:\\jupDef\\Don\\TestImages\\ED127 StarOnly\\L-"
file_series ="F:\\jupDef\\Don\\TestImages\\RC8 StarJupiter 48arcsec\\L-"
#"F:\\jupDef\\2023-11-07\\sao93016\\21_27_04\\sao93016_"

firstfile = 1              # the first file might not start with 00001
fileext="fits"             # file extension
Numberformat = '{:05}'     # number of digits in the file numbering list
#dark = "C:\\Users\\sales\\Desktop\\2023-09-25\\MasterDark.fits"

file_output = "3StarsNovBFrev2.csv" 
################ parameters listed here
sat = 60000                # value at which the image is considered saturated
printnumber = 10            # number of files to pass before printing file number
ring_inner_estimate = 13   # ring_inner radius for first estimates nom 15
ring_middle_estimate = 18  # ring_middle radius for first estimates nom 25
ring_outer_estimate = 20   # ring_outer radius for first estimates nom 35
ring_inner_final = 7       # ring_inner radius for final calculations nom 7 10
ring_middle_final = 14     # ring_middle radius for final calculations nom 14
ring_outer_final = 18      # ring_outer radius for final calculations nom 18 20
###################################

#jupiter cleanup parameters
######################################################################################
useJupCorrect = 4       #  1- 2ang  2- 180  3- both   4 subtract earlier jupiter
jupECX=1103.0   # xposition for earlier jupiter image  (needs to be within 10 pixels of center)
jupECY=635.8    # xposition
jupEarlier= "F:\\jupDef\\Don\\TestImages\\RC8 StarJupiter 64arcsec\\L-00001-C.fits"
useGradientCorrect= 1   # use a gradient cleanup at end of options above
useAlign2jup=0  # this does a search to align two jupiter images as part of option 4  it takes 1-2 sec per image  default is off
#######################################################################################


########### programming starts here

#10-30-23BF create a box to hold the ave jupershine background from cw  and ccw from the target star 
jBacsize= 2*ring_outer_estimate+4
imageJBacBox= np.zeros((jBacsize,jBacsize))
darkAve=1300.0     #true dark - 100 pedestal
#end BF addin    #10-30-23BF needed for jup shine correction

##dark_image = fits.getdata(dark, ext=0)
## add flat field file here
##flat_image = fits.getdata(flat, ext=0)
##### read centroid parameters from file
CP = np.loadtxt(parameter_file, delimiter=",",max_rows = 6)
XX =  int(CP[0][0]); YY =  int(CP[0][1])
X1Est=int(CP[1][0]); Y1Est=int(CP[1][1]); 
X2Est=int(CP[2][0]); Y2Est=int(CP[2][1]);
X3Est=int(CP[3][0]); Y3Est=int(CP[3][1]);
JXEst=int(CP[4][0]); JYEst=int(CP[4][1]); 
number_of_files=int(CP[5][0]); box=int(CP[5][1]); #box==odd, diam(1D centroid)
print ("number of input files =", number_of_files)
clip_multiplier = 2.874   # Howell's method ratio 2.874

#BF this is code used to set up image of jupiter used for testing shine removal
#1hr30m to 1hr before
#sepStartDelX= 0.0   #pixel
#sepStartDelY= -110.8  #pixel
#sepEndDelX= 16.5   #pixel
#sepEndDelY= -100.6  #pixel

#30 min to closest 
sepStartDelX= 33.1   #pixel
sepStartDelY= -93.6  #pixel
sepEndDelX= 49.6   #pixel
sepEndDelY= -83.4  #pixel

incX=sepEndDelX-sepStartDelX
incY=sepEndDelY-sepStartDelY
addJup=0     # this is used to add a simulated jupiter for testing jupiter glow removal methods 
if(addJup>0):
    jupglow= "F:\\jupDef\\2023-11-07\\jupiter\\21_20_06\\jupiter11-7crop500_05.fits"
    #"C:\\jupiterD\\jupglow\\11-7\\jup11-7-23_60-1s_crop.fit"   
    imageJup= fits.getdata(jupglow,ext=0)    
    jShineSize=500
    halfjShineSize=250    
    for j in np.arange( 0 ,jShineSize,1):
            for i in np.arange( 0 ,jShineSize,1):
                if(imageJup[j][i] >64000):imageJup[j][i]=50000     # to avoid overflow artifacts
                
    #hduj = fits.PrimaryHDU(np.uint16(imageJup))
    #hdulj = fits.HDUList([hduj])
    #hdulj.writeto("C:\\jupiterd\\JupTest2Crop.fits",overwrite=True)     
#bf end jupiter  addin image setup

#12-20=23 setup file if subtract earlier jupiter is used
#jupECX=1103.0  moved earlier
#jupECY=635.8
#sameDayJup=1   
if (useJupCorrect==4 ):    #12-20  change from samedayJup >0 to this to match correction
    imageJupEarlier= fits.getdata(jupEarlier,ext=0)
    xdimJupE = np.size(imageJupEarlier,1)
    ydimJupE = np.size(imageJupEarlier,0)
    k=-1
    if k== 0 :   # -1 ignores print
        hduje = fits.PrimaryHDU(np.uint16(imageJupEarlier))
        hdulje = fits.HDUList([hduje])
        hdulje.writeto("C:\\jupiterd\\JupEarlier.fits",overwrite=True)
        
    jupECX, jupECY =centroid_JupbyEdge(imageJupEarlier, int(jupECX),int(jupECY))


######## make subarray masks, only need one time, never changes
estimate_backgnd_mask = np.zeros((2*ring_outer_estimate, 2*ring_outer_estimate))
for j in np.arange(0, 2*ring_outer_estimate, 1):
    for i in np.arange(0, 2*ring_outer_estimate, 1):
        if ((i-ring_outer_estimate)**2 + (j-ring_outer_estimate)**2)  \
            < (ring_outer_estimate)**2 and \
            ((i-ring_outer_estimate)**2 + (j-ring_outer_estimate)**2)  \
            > (ring_middle_estimate)**2:
                estimate_backgnd_mask[i][j] = 1
final_backgnd_mask = np.zeros((2*ring_outer_final, 2*ring_outer_final))
for j in np.arange(0, 2*ring_outer_final, 1):
    for i in np.arange(0, 2*ring_outer_final, 1):
        if ((i-ring_outer_final)**2 + (j-ring_outer_final)**2)  \
            < (ring_outer_final)**2 and \
            ((i-ring_outer_final)**2 + (j-ring_outer_final)**2)  \
            > (ring_middle_final)**2:
                final_backgnd_mask[i][j] = 1
estimate_centroid_mask = np.zeros((2*ring_inner_estimate, 2*ring_inner_estimate))
for j in np.arange(0, 2*ring_inner_estimate, 1):
    for i in np.arange(0, 2*ring_inner_estimate, 1):
        if ((i-ring_inner_estimate)**2 + (j-ring_inner_estimate)**2)  \
            < (ring_inner_estimate)**2:
                estimate_centroid_mask[i][j] = 1   
final_centroid_mask = np.zeros((2*ring_inner_final, 2*ring_inner_final))
for j in np.arange(0, 2*ring_inner_final, 1):
    for i in np.arange(0, 2*ring_inner_final, 1):
        if ((i-ring_inner_final)**2 + (j-ring_inner_final)**2)  \
            < (ring_inner_final)**2:
                final_centroid_mask[i][j] = 1              
#plt.imshow(final_centroid_mask, cmap='gray')
## repeat for Jupiter rings
#bFjupiter masks no longer needed

#######start the loop to measure the images
Results = open(file_output, 'a')
for jj in np.arange(0, number_of_files, 1):  # typically 5 to 5000 loops
    image = np.zeros((YY, XX))         # clear the image each loop
    current_file = file_series + Numberformat.format(jj+firstfile) + "-A." +fileext
    #current_file = file_series + Numberformat.format(jj+firstfile) + "." +fileext
    hdul = fits.open(current_file)
    image = fits.getdata(current_file, ext=0)
    
    rX2Est,rY2Est,fwhm= roughCentroidFwhm(image, X2Est,Y2Est, jBacsize)
    X2Est=int(rX2Est)
    Y2Est= int(rY2Est)
#BF for now we are not doing dark or flat field calibreation   
##  subtract dark image one time at start
##    image = image - dark_image
#      image=image-1100    # fix pedestal instead of dark         
##  add flat field division here
##    image = image/flat_image


    #BF 11-3 code to add jupiter for testing
    jupGlowY = JYEst
    jupGlowX = JXEst
   
    #addJup=1 this is define earlier in file setup for this feature
    if(addJup>0):
        curSepDelX= int(sepStartDelX + (jj*incX)/number_of_files )
        curSepDelY= int(sepStartDelY + (jj*incY)/number_of_files )
        jupGlowX=X2Est-curSepDelX
        jupGlowY=Y2Est - curSepDelY
        for j in np.arange( 0 ,jShineSize,1):
            newY=int(j+Y2Est-halfjShineSize -curSepDelY)
            for i in np.arange( 0 ,jShineSize,1):
                newX=int(i-halfjShineSize+X2Est- curSepDelX)
                image[newY][newX] =  image[newY][newX] + imageJup[j][i]
    if jj== -1:   # -1 ignores print
        hduj = fits.PrimaryHDU(np.uint16(image))
        hdulj = fits.HDUList([hduj])
        hdulj.writeto("C:\\jupiterd\\JupTestEarlyJupAdded.fits",overwrite=True)
    #end code to add jupiter for testing 

    #BF version with two roi shine fix location of jupiter from prior estimate
    #useJupCorrect = 4       #  1- 2ang  2- 180  3- both   4 subtract earlier jupiter   moved to start of program 
    if(useJupCorrect > 0  and useJupCorrect <4):
        jupCX2, jupCY2 =centroid_JupbyEdge(image, int(jupGlowX),int(jupGlowY))   #for possible furture expansion to ellipse model of jupiter  add ,95,101)
        JXfinal=jupCX2
        JYfinal=jupCY2
        
        jupGlowX=jupCX2
        jupGlowY=jupCY2
        sepsq= (rY2Est-jupGlowY)**2 +  (rX2Est-jupGlowX)**2 
        sep= sqrt(sepsq)
     
        theta= np.arctan2((jupGlowY-rY2Est)/sep,(jupGlowX-rX2Est)/sep)   
        thetaD= np.degrees( theta)
        shiftRad= (jBacsize +14 )/sep   #add sub  38 degrees .662 rad
    
        imageJBacBox = np.zeros((YY, XX))
        halfjBacsize= int(jBacsize/2)
        for j in np.arange( 0 ,jBacsize,1):
            pixY=int(j+Y2Est-halfjBacsize)       
            for i in np.arange( 0 ,jBacsize,1):
                pixX=int(i-halfjBacsize+X2Est)           
                pixSep= sqrt( (pixY-jupGlowY)**2 +  (pixX-jupGlowX)**2) 
                sepAsec=0.491*pixSep
                pixTheta= np.arctan2((jupGlowY-pixY)/pixSep,(jupGlowX-pixX)/pixSep)            
                pixTheataD= np.degrees(pixTheta)
                
                if(useJupCorrect == 1  or useJupCorrect == 3):
                    #single pixel version
                    pixXcw=  int(jupGlowX- pixSep*np.cos(pixTheta+ shiftRad) )
                    pixYcw=  int(jupGlowY- pixSep*np.sin(pixTheta+ shiftRad) )
                    pixXccw= int(jupGlowX- pixSep*np.cos(pixTheta- shiftRad) )
                    pixYccw= int(jupGlowY- pixSep*np.sin(pixTheta- shiftRad) )
                    #imageJBacBox[j][i]= (image[pixYcw][pixXcw] + image[pixYccw][pixXccw])/2
                    #interpolated pixel version
                    valcw= interpolatePixValue(image, (jupGlowX- pixSep*np.cos(pixTheta+ shiftRad)) , (jupGlowY- pixSep*np.sin(pixTheta+ shiftRad)) )
                    valccw= interpolatePixValue( image, (jupGlowX- pixSep*np.cos(pixTheta- shiftRad)), (jupGlowY- pixSep*np.sin(pixTheta- shiftRad))  )
                    imageJBacBox[j][i]=(valcw +valccw)/2.0
                    if(jj==-1):                    
                        image[pixYccw][pixXccw]=8000+10*j+i    # mark where samples taken diagnostic only can mess up correction since target pixel may be used more than once
                        image[pixYcw][pixXcw]=9000+10*j+i

                if(useJupCorrect == 2  or useJupCorrect == 3):
                    pixX180= int(jupGlowX- pixSep*np.cos(pixTheta -3.14159) )
                    pixY180= int(jupGlowY- pixSep*np.sin(pixTheta -3.14159) )
                    val180= interpolatePixValue(image,jupGlowX- pixSep*np.cos(pixTheta -3.14159), jupGlowY- pixSep*np.sin(pixTheta -3.14159)  )
                    
                    #imageJBacBox[j][i]=  image[pixY180][pixX180]
                    imageJBacBox[j][i]= val180
                    #image[pixY180][pixX180]=10000+10*j+i  # mark where samples taken diagnostic only
                if(useJupCorrect==3): 
                    #imageJBacBox[j][i]= (image[pixYcw][pixXcw] + image[pixYccw][pixXccw] +image[pixY180][pixX180])/3
                    imageJBacBox[j][i]= (valcw +valccw +val180)/3.0
        if jj==-1 :   # -1 ignores print  
            hduj = fits.PrimaryHDU(np.uint16(image))
            hdulj = fits.HDUList([hduj])
            hdulj.writeto("C:\\jupiterd\\JupTest12-11DonSimRC8_jupEcor.fits",overwrite=True)                       
 #                 
        #image = np.zeros((YY, XX))
        #if(jj==30): print ("before pix x y val= 1348,398 ", image[398][1348])
        for j in np.arange( 0 ,jBacsize,1):
            newY=int(j+Y2Est-halfjBacsize)
            for i in np.arange( 0 ,jBacsize,1):
                newX=int(i-halfjBacsize+X2Est)
                 
                image[newY][newX] = 200.0 + darkAve + image[newY][newX] - imageJBacBox[j][i]   # 200. is pedesstall to avoid underflow
                image[j][i]=imageJBacBox[j][i]-darkAve -200.  # diagnostic shows the correction in the corner of the image
                if(image[newY][newX] >60000):
                    print("sat pix frame",jj ,"x,y",newX,newY,"val= ",image[newY][newX] )
                if(image[newY][newX] <0):
                    print("neg pix frame",jj ,"x,y",newX,newY,"val= ",image[newY][newX] )   
        #if(jj==30): print ("after pix x y val= 1348,398 ", image[398][1348] )
        
        if useGradientCorrect == 1 : subtractGradient(image,X2Est,Y2Est,jBacsize)
        #subtractGradient(image,X1Est,Y1Est,jBacsize)
        #subtractGradient(image,X3Est,Y3Est,jBacsize)
                
        if jj==0 :   # -1 ignores print  
            hduj = fits.PrimaryHDU(np.uint16(image))
            hdulj = fits.HDUList([hduj])
            hdulj.writeto("C:\\jupiterd\\JupTest12-11DonSimRC8_30.fits",overwrite=True)
        
        
        
     #useJupCorrect=4      # subtract supiter from at least an hour earlier method of shine reduction 
    if (useJupCorrect == 4):   
        if jj== -1:   # -1 ignores print
            hduj = fits.PrimaryHDU(np.uint16(image))
            hdulj = fits.HDUList([hduj])
            hdulj.writeto("C:\\jupiterd\\JupTestEarlyJupAdded.fits",overwrite=True)
            
        #testX,testY= centroid_JupbyEdge(imageJupEarlier,int(jupECX), int(jupECY) )         

        jupXimage,jupYimage= centroid_JupbyEdge(image,int(jupGlowX), int(jupGlowY) ) #start by finding centroid of  jupiter in current jupiter
        JXfinal=jupXimage
        JYfinal=jupYimage
        roisize=jBacsize # 300 #jBacsize
        halfroisize=int(roisize/2)
        rdifX=jupECX-jupXimage
        rdifY=jupECY-jupYimage
        sep= sqrt(rdifX**2 + rdifY**2)
        difX= int(jupECX+.5) -int(jupXimage+.5)
        difY= int(jupECY+.5) -int(jupYimage+.5)
        #halfjBacsize= int(jBacsize/2) 
        #imageJupEarlier2 = imageJupEarlier.copy()

        
        if useAlign2jup== 1 :
            xoffsetout,yoffsetout, avepixout = align2Jup(image, imageJupEarlier, jupXimage,jupYimage,rdifX,rdifY,darkAve)
            rdifx=xoffsetout
            rdify=yoffsetout
# now subtract the earlier jupiter image from current image
#          
        for j in np.arange( 0 ,roisize,1):
            yimagepix= Y2Est- halfroisize +j 
            yjup=yimagepix+difY  
            ryjup=rdifY+yimagepix #+0.09 #+0.095
            for i in np.arange(0 ,roisize,1):
                ximagepix= X2Est- halfroisize +i 
                xjup=ximagepix + difX 
                rxjup=rdifX+ximagepix # +1.5 #+2.1  
                jupval=interpolatePixValue(imageJupEarlier, rxjup, ryjup)
                image[yimagepix][ximagepix] =  image[yimagepix][ximagepix]  + darkAve - jupval # add dark ave back image so no negatives
                
                #image[j][i]= imageJupEarlier[yjup][xjup]  # correction in corner
                
                              
                
        #if jj==0 :   # -1 ignores print  
        #    hduje = fits.PrimaryHDU(np.uint16(imageJupEarlier))
        #    hdulje = fits.HDUList([hduje])
        #    hdulje.writeto("C:\\jupiterd\\JupEarlier0.fits",overwrite=True)

        if useGradientCorrect == 1 :subtractGradient(image,X2Est,Y2Est,roisize)
        if jj==-1:
            current_file_out = "C:\\jupiterd\\out2\\DonSim_Ed_64-48_out" + Numberformat.format(jj+firstfile)  + "." + fileext
            hduj = fits.PrimaryHDU(np.uint16(image))
            hdulj = fits.HDUList([hduj])
            hdulj.writeto(current_file_out,overwrite=True) 

        if jj==0 :   # -1 ignores print  
            hduj = fits.PrimaryHDU(np.uint16(image))
            hdulj = fits.HDUList([hduj])
            hdulj.writeto("C:\\jupiterd\\JupTest12-13_DonRC8-48-64_0.fits",overwrite=True)
        if jj==399 :   # -1 ignores print  
            hduj = fits.PrimaryHDU(np.uint16(image))
            hdulj = fits.HDUList([hduj])
            hdulj.writeto("C:\\jupiterd\\JupTest12-11DonED_48-64_399.fits",overwrite=True)   
            
    if(useJupCorrect ==0 ):   # needed for jupiter centroid print 
        JXfinal=JXEst
        JYfinal=JYEst
            
            ########## initialize the background arrays
    B1_estimate = np.zeros((2*ring_outer_estimate, 2*ring_outer_estimate))
    B2_estimate = np.zeros((2*ring_outer_estimate, 2*ring_outer_estimate))
    B3_estimate = np.zeros((2*ring_outer_estimate, 2*ring_outer_estimate))
   
    ####### form the estimate background subarrays
    for j in np.arange(0, 2*ring_outer_estimate, 1):
        for i in np.arange(0, 2*ring_outer_estimate, 1):
            B1_estimate[j][i] = image [j+Y1Est-ring_outer_estimate][i+X1Est-ring_outer_estimate]
            B2_estimate[j][i] = image [j+Y2Est-ring_outer_estimate][i+X2Est-ring_outer_estimate]
            B3_estimate[j][i] = image [j+Y3Est-ring_outer_estimate][i+X3Est-ring_outer_estimate]    
    #### mask the estimate background subarrays
    B1Mask_estimate=B1_estimate * estimate_backgnd_mask
    B2Mask_estimate=B2_estimate * estimate_backgnd_mask
    B3Mask_estimate=B3_estimate * estimate_backgnd_mask   
    ##### assign NaN outside mask
    ##### this is needed to get correct RMS and mean by ignoring the NaN elements
    for j in np.arange(0, 2*ring_outer_estimate, 1):
        for i in np.arange(0, 2*ring_outer_estimate, 1):
            if B1Mask_estimate[j][i] == 0 : B1Mask_estimate[j][i] = "NaN"
            if B2Mask_estimate[j][i] == 0 : B2Mask_estimate[j][i] = "NaN"
            if B3Mask_estimate[j][i] == 0 : B3Mask_estimate[j][i] = "NaN"
    ##### calculate the standard deviations and means
    B1_estimate_RMS =  np.nanstd(B1Mask_estimate)       
    B2_estimate_RMS =  np.nanstd(B2Mask_estimate)       
    B3_estimate_RMS =  np.nanstd(B3Mask_estimate)       
        
    B1_estimate_mean = np.nanmean(B1Mask_estimate) 
    B2_estimate_mean = np.nanmean(B2Mask_estimate)
    B3_estimate_mean = np.nanmean(B3Mask_estimate)
    
    ##### for the centroids, use a smaller, inner array
    C1_estimate = np.zeros((2*ring_inner_estimate, 2*ring_inner_estimate))
    C2_estimate = np.zeros((2*ring_inner_estimate, 2*ring_inner_estimate))
    C3_estimate = np.zeros((2*ring_inner_estimate, 2*ring_inner_estimate))
    
    ####### form the estimate centroid subarrays
    for j in np.arange(0, 2*ring_inner_estimate, 1):
        for i in np.arange(0, 2*ring_inner_estimate, 1):
            C1_estimate[j][i] = image [j+Y1Est-ring_inner_estimate][i+X1Est-ring_inner_estimate]
            C2_estimate[j][i] = image [j+Y2Est-ring_inner_estimate][i+X2Est-ring_inner_estimate]
            C3_estimate[j][i] = image [j+Y3Est-ring_inner_estimate][i+X3Est-ring_inner_estimate]
    
    #plt.imshow(C3_estimate, cmap='gray')
    #hdum = fits.PrimaryHDU(np.uint16(C2_estimate))
    #hdulm = fits.HDUList([hdum])
    #hdulm.writeto("C:\\jupiterd\\C2estimate.fits",overwrite=True)
  
    
    ##### check if any of the stars are saturated
    if (np.max(C1_estimate)>sat): print ("Star1 is saturated in file", jj+1)
    if (np.max(C2_estimate)>sat): print ("Star2 is saturated in file", jj+1)
    if (np.max(C3_estimate)>sat): print ("Star3 is saturated in file", jj+1)  
    ##### subtract background mean and clip_multiplier*RMS
    C1_estimate_centroid = C1_estimate-B1_estimate_mean-clip_multiplier*B1_estimate_RMS
    C2_estimate_centroid = C2_estimate-B2_estimate_mean-clip_multiplier*B2_estimate_RMS
    C3_estimate_centroid = C3_estimate-B3_estimate_mean-clip_multiplier*B3_estimate_RMS
   
    #### mask the estimate centroid subarrays, putting zeros in the corners
    C1Mask_estimate=C1_estimate_centroid * estimate_centroid_mask
    C2Mask_estimate=C2_estimate_centroid * estimate_centroid_mask
    C3Mask_estimate=C3_estimate_centroid * estimate_centroid_mask
   
    ##### set negative values in the clipped sub-arrays to zero
    for j in np.arange(0, 2*ring_inner_estimate, 1):
        for i in np.arange(0, 2*ring_inner_estimate, 1):
            if C1Mask_estimate[j][i] < 0: C1Mask_estimate[j][i] = 0
            if C2Mask_estimate[j][i] < 0: C2Mask_estimate[j][i] = 0
            if C3Mask_estimate[j][i] < 0: C3Mask_estimate[j][i] = 0
    
    ##### initialize the weighting arrays, using 1-D arrays
    C1X_estimate_weights = np.zeros(2*ring_inner_estimate)
    C1Y_estimate_weights = np.zeros(2*ring_inner_estimate)
    C2X_estimate_weights = np.zeros(2*ring_inner_estimate)
    C2Y_estimate_weights = np.zeros(2*ring_inner_estimate)
    C3X_estimate_weights = np.zeros(2*ring_inner_estimate)
    C3Y_estimate_weights = np.zeros(2*ring_inner_estimate)
   
    ##### calculate the weighting arrays
    for i in np.arange(0, 2*ring_inner_estimate, 1):
        C1X_estimate_weights[i] = i+X1Est-ring_inner_estimate
        C2X_estimate_weights[i] = i+X2Est-ring_inner_estimate
        C3X_estimate_weights[i] = i+X3Est-ring_inner_estimate
        C1Y_estimate_weights[i] = i+Y1Est-ring_inner_estimate
        C2Y_estimate_weights[i] = i+Y2Est-ring_inner_estimate
        C3Y_estimate_weights[i] = i+Y3Est-ring_inner_estimate
    
##### sum over all rows and columns to get the total signal
    C1sum_estimate = np.sum(C1Mask_estimate)
    C2sum_estimate = np.sum(C2Mask_estimate)
    C3sum_estimate = np.sum(C3Mask_estimate)
    
    ##### sum the rows and columns to start the moment calculations
    C1Xsum_estimate = np.sum(C1Mask_estimate, axis=0)
    C1Ysum_estimate = np.sum(C1Mask_estimate, axis=1)
    C2Xsum_estimate = np.sum(C2Mask_estimate, axis=0)
    C2Ysum_estimate = np.sum(C2Mask_estimate, axis=1)
    C3Xsum_estimate = np.sum(C3Mask_estimate, axis=0)
    C3Ysum_estimate = np.sum(C3Mask_estimate, axis=1)
    
    ##### calculate the moments to get the centroids
    Centroid1X_estimate = np.sum(C1Xsum_estimate*C1X_estimate_weights)/C1sum_estimate
    Centroid1Y_estimate = np.sum(C1Ysum_estimate*C1Y_estimate_weights)/C1sum_estimate
    Centroid2X_estimate = np.sum(C2Xsum_estimate*C2X_estimate_weights)/C2sum_estimate
    Centroid2Y_estimate = np.sum(C2Ysum_estimate*C2Y_estimate_weights)/C2sum_estimate
    Centroid3X_estimate = np.sum(C3Xsum_estimate*C3X_estimate_weights)/C3sum_estimate
    Centroid3Y_estimate = np.sum(C3Ysum_estimate*C3Y_estimate_weights)/C3sum_estimate
   
    ##### calculate the new array centroids based on estimate arrays
    X1final = int(round(Centroid1X_estimate))
    Y1final = int(round(Centroid1Y_estimate))
    X2final = int(round(Centroid2X_estimate))
    Y2final = int(round(Centroid2Y_estimate))
    X3final = int(round(Centroid3X_estimate))
    Y3final = int(round(Centroid3Y_estimate))
    
    ##### repeat previous calculations, replacing "estimate" with "final"
    ########## initialize the background arrays
    B1_final = np.zeros((2*ring_outer_final, 2*ring_outer_final))
    B2_final = np.zeros((2*ring_outer_final, 2*ring_outer_final))
    B3_final = np.zeros((2*ring_outer_final, 2*ring_outer_final))
   
    ####### form the final background subarrays
    for j in np.arange(0, 2*ring_outer_final, 1):
        for i in np.arange(0, 2*ring_outer_final, 1):
            B1_final[j][i] = image [j+Y1final-ring_outer_final][i+X1final-ring_outer_final]
            B2_final[j][i] = image [j+Y2final-ring_outer_final][i+X2final-ring_outer_final]
            B3_final[j][i] = image [j+Y3final-ring_outer_final][i+X3final-ring_outer_final]
##            print(B3_final, ring_outer_final)
   
    #### mask the estimate background subarrays
    B1Mask_final=B1_final * final_backgnd_mask
    B2Mask_final=B2_final * final_backgnd_mask
    B3Mask_final=B3_final * final_backgnd_mask
   
    ##### assign NaN outside mask
    ##### this is needed to get correct STD and mean by ignoring the NaN elements
    for j in np.arange(0, 2*ring_outer_final, 1):
        for i in np.arange(0, 2*ring_outer_final, 1):
            if B1Mask_final[j][i] == 0 : B1Mask_final[j][i] = "NaN"
            if B2Mask_final[j][i] == 0 : B2Mask_final[j][i] = "NaN"
            if B3Mask_final[j][i] == 0 : B3Mask_final[j][i] = "NaN"
   
    ##### calculate the standard deviations and means
    B1_final_RMS =  np.nanstd(B1Mask_final)       
    B2_final_RMS =  np.nanstd(B2Mask_final)       
    B3_final_RMS =  np.nanstd(B3Mask_final)       
       
    B1_final_mean = np.nanmean(B1Mask_final)
    B2_final_mean = np.nanmean(B2Mask_final)  
    B3_final_mean = np.nanmean(B3Mask_final)  
      
    ##### for the centroids, use a smaller, inner array
    C1_final = np.zeros((2*ring_inner_final, 2*ring_inner_final))
    C2_final = np.zeros((2*ring_inner_final, 2*ring_inner_final))
    C3_final = np.zeros((2*ring_inner_final, 2*ring_inner_final))
   
    ####### form the final centroid subarrays
    for j in np.arange(0, 2*ring_inner_final, 1):
        for i in np.arange(0, 2*ring_inner_final, 1):
            C1_final[j][i] = image [j+Y1final-ring_inner_final][i+X1final-ring_inner_final]
            C2_final[j][i] = image [j+Y2final-ring_inner_final][i+X2final-ring_inner_final]
            C3_final[j][i] = image [j+Y3final-ring_inner_final][i+X3final-ring_inner_final]
   
    ##### subtract mean to calculate SNRs, then mask
    C1_forSNR =  (C1_final-B1_final_mean) * final_centroid_mask 
    C2_forSNR =  (C2_final-B2_final_mean) * final_centroid_mask  
    C3_forSNR =  (C3_final-B3_final_mean) * final_centroid_mask  
   
## insert photutils here to calculate the centroids via 1-D gaussian    
    C1DX1, C1DY1 = centroid_sources(image, X1final, Y1final, box_size=box, centroid_func=centroid_1dg)                                       
    C1DX2, C1DY2 = centroid_sources(image, X2final, Y2final, box_size=box, centroid_func=centroid_1dg)                                       
    C1DX3, C1DY3 = centroid_sources(image, X3final, Y3final, box_size=box, centroid_func=centroid_1dg) 
    ##### subtract clip_multiplier*RMS from B1-BJ     
    C1_final_centroid = (C1_forSNR-clip_multiplier*B1_final_RMS)*final_centroid_mask
    C2_final_centroid = (C2_forSNR-clip_multiplier*B2_final_RMS)*final_centroid_mask
    C3_final_centroid = (C3_forSNR-clip_multiplier*B3_final_RMS)*final_centroid_mask
    
    ##### set negative values in the clipped sub-arrays to zero
    for j in np.arange(0, 2*ring_inner_final, 1):
        for i in np.arange(0, 2*ring_inner_final, 1):
            if C1_final_centroid[j][i] < 0: C1_final_centroid[j][i] = 0
            if C2_final_centroid[j][i] < 0: C2_final_centroid[j][i] = 0
            if C3_final_centroid[j][i] < 0: C3_final_centroid[j][i] = 0
    
    ##### initialize the weighting arrays, using 1-D arrays
    C1X_final_weights = np.zeros(2*ring_inner_final)
    C1Y_final_weights = np.zeros(2*ring_inner_final)
    C2X_final_weights = np.zeros(2*ring_inner_final)
    C2Y_final_weights = np.zeros(2*ring_inner_final)
    C3X_final_weights = np.zeros(2*ring_inner_final)
    C3Y_final_weights = np.zeros(2*ring_inner_final)
    
    ##### calculate the weighting arrays
    for i in np.arange(0, 2*ring_inner_final, 1):
        C1X_final_weights[i] = i+X1final-ring_inner_final
        C2X_final_weights[i] = i+X2final-ring_inner_final
        C3X_final_weights[i] = i+X3final-ring_inner_final
        C1Y_final_weights[i] = i+Y1final-ring_inner_final
        C2Y_final_weights[i] = i+Y2final-ring_inner_final
        C3Y_final_weights[i] = i+Y3final-ring_inner_final
    
    ##### sum over all rows and columns to get the total signal
    C1sum_final = np.sum(C1_final_centroid)
    C2sum_final = np.sum(C2_final_centroid)
    C3sum_final = np.sum(C3_final_centroid)
    
    ##### sum the rows and columns to start the moment calculations
    C1Xsum_final = np.sum(C1_final_centroid, axis=0)
    C1Ysum_final = np.sum(C1_final_centroid, axis=1)
    C2Xsum_final = np.sum(C2_final_centroid, axis=0)
    C2Ysum_final = np.sum(C2_final_centroid, axis=1)
    C3Xsum_final = np.sum(C3_final_centroid, axis=0)
    C3Ysum_final = np.sum(C3_final_centroid, axis=1)
   
    ##### calculate the moments to get the centroids
    Centroid1X_final = np.sum(C1Xsum_final*C1X_final_weights)/C1sum_final
    Centroid1Y_final = np.sum(C1Ysum_final*C1Y_final_weights)/C1sum_final
    Centroid2X_final = np.sum(C2Xsum_final*C2X_final_weights)/C2sum_final
    Centroid2Y_final = np.sum(C2Ysum_final*C2Y_final_weights)/C2sum_final
    Centroid3X_final = np.sum(C3Xsum_final*C3X_final_weights)/C3sum_final
    Centroid3Y_final = np.sum(C3Ysum_final*C3Y_final_weights)/C3sum_final
    
    ##### calculate the SNRs
    SNR1final = np.sum(C1_forSNR)/(B1_final_RMS*1.77*ring_inner_final)
    SNR2final = np.sum(C2_forSNR)/(B2_final_RMS*1.77*ring_inner_final)
    SNR3final = np.sum(C3_forSNR)/(B3_final_RMS*1.77*ring_inner_final)
    ######## save centroid results and SNR values
## now average the values from Howell's COM and the 1Dgaussian  
##
    Results.write("{:s},,{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:f},,{:f},{:f},{:f},{:f},{:f},{:f},{:f}\n"\
                  .format(hdul[0].header['DATE-AVG'],Centroid1X_final,Centroid1Y_final,\
                   Centroid2X_final,Centroid2Y_final,\
                   Centroid3X_final,Centroid3Y_final,\
                   JXfinal,JYfinal,\
                   int(SNR1final),int(SNR2final),int(SNR3final),\
                   C1DX1[0],C1DY1[0],C1DX2[0],C1DY2[0],C1DX3[0],C1DY3[0],box   ))      
    #### use the final coordinates as starting point for the next average file series  
    X1Est = X1final; Y1Est = Y1final; 
    X2Est = X2final; Y2Est = Y2final;
    X3Est = X3final; Y3Est = Y3final; 
    JXEst = JXfinal; JYEst = JYfinal; 
    testjj = jj+1-(int((jj+1)/printnumber))*printnumber;
    if (testjj < 0.5): print("done with output file number", jj+1)      

##### done with loops                                                                           
print("done with", number_of_files, "output files")        
Results.close()  
 
#BF main code end 


#extra code  not ready to delete yet         
##plt.imshow(image, cmap='gray')



# 
# this was first version of gradient subtract code before I moved it to a function        
        #now compute the remaining slope in the background and subtract it 
        #sumN=0   #for sum of north row
        #sumE=0   #for sum of east row
        #sumS=0   #for sum of south row
        #sumW=0   #for sum of West row
        #newXE= int(-halfjBacsize+X2Est) #west,east
        #newXW= int(halfjBacsize+X2Est-1) #west,east
        #for j in np.arange( 0 ,jBacsize,1):
        #    newY=int(j+Y2Est-halfjBacsize) #west
        #    sumW= sumW + image[newY][newXW]
        #    sumE= sumE + image[newY][newXE]
        #newYS= int(-halfjBacsize+Y2Est) #west,east
        #newYN= int(halfjBacsize+Y2Est-1) #west,east   
        #for i in np.arange( 0 ,jBacsize,1):
        #    newX=int(i+X2Est-halfjBacsize) #west
        #    sumN= sumN + image[newYN][newX]
        #    sumS= sumS + image[newYS][newX]   
        #aveN=sumN/jBacsize
        #aveE=sumE/jBacsize
        #aveS=sumS/jBacsize
        #aveW=sumW/jBacsize
        #slopeY= (aveN-aveS)/jBacsize
        #slopeX= (aveW-aveE)/jBacsize
        #for j in np.arange( 0 ,jBacsize,1): 
        #    newY=int(j+Y2Est-halfjBacsize)
        #    for i in np.arange( 0 ,jBacsize,1):
        #        newX=int(i-halfjBacsize+X2Est)
        #         
        #       image[newY][newX] = image[newY][newX] +  (jBacsize-i)*slopeX +(jBacsize-j)*slopeY


'''
test interpolation 
testdata= np.zeros((2,3))
testdata[0][0] = 10
testdata[0][1] = 20
testdata[0][2] = 30
testdata[1][0] = 1010
testdata[1][1] = 1020
testdata[1][2] = 1030
out= interpolatePixValue(testdata, .0, .0)
'''