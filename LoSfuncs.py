import requests, json
import numpy as np
from ftplib import FTP
from osgeo import ogr
from osgeo import osr
import os.path
import zipfile
import time
from math import floor, ceil, atan2, pi, sqrt
import pickle
import simplekml
from osgeo import gdal
import urllib.request

def loadMapsWCS(bbox, dim, mapType):
    if os.path.isfile('wcsMap.tif'):
        os.remove('wcsMap.tif')

    #dtmMap = np.memmap('Kort/terrainMap',dtype='float32',mode='w+',shape=(dim[0],dim[1]))
    mapFile = 'wcsMap.tif'

    #BBOX = "BBOX="+str(858275)+","+str(6108438)+","+str(897160)+","+str(6144049)
    BBOX = "BBOX="+str(bbox[0])+","+str(bbox[1])+","+str(bbox[2])+","+str(bbox[3])

    DIM = "WIDTH="+str(dim[0])+"&HEIGHT="+str(dim[1])
    if mapType == 't':
        coverage = "dhm_terraen"
    elif mapType == 's':
        coverage = "dhm_overflade"
    url = "http://services.kortforsyningen.dk/dhm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=1.0.0&COVERAGE="+coverage+"&CRS=EPSG:25832&"+BBOX+"&"+DIM+"&FORMAT=image/gtiff&login=ASN&password=ASN"
    print(url)
    urllib.request.urlretrieve(url,mapFile)
    
    m = gdal.Open(mapFile)
    dtmMap = np.array(m.GetRasterBand(1).ReadAsArray())
    
    return dtmMap

def loadMapsHybrid(nRange, mRange, mapType):
    inputSize = 2500
    scaledSize = inputSize

    if mapType == 't':
        prefix = 'Kort/DTM_1km_'
        
    elif mapType == 's':
        prefix = 'Kort/DSM_1km_'

    Map = np.zeros([scaledSize*len(mRange),scaledSize*(len(nRange))])
    
    for n in nRange:
        #print(n)
        for m in mRange:
             
            nStart = (n-nRange[0])*inputSize
            mStart = (mRange[-1]-m)*inputSize

            nEnd = nStart+inputSize
            mEnd = mStart+inputSize
       
            fileName = prefix+str(m)+'_'+str(n)+'.tif'
            if os.path.isfile(fileName):
                dt = gdal.Open(fileName)
                Map[mStart:mEnd,nStart:nEnd] = np.array(dt.GetRasterBand(1).ReadAsArray())
                dt=None
            else:
                #BBOX = "BBOX="+str(858275)+","+str(6108438)+","+str(897160)+","+str(6144049)
                bbox = [n*1000,m*1000,(n+1)*1000,(m+1)*1000]
                print(bbox)
                Map[mStart:mEnd,nStart:nEnd] = loadMapsWCS(bbox, [2500,2500], mapType)

    return np.fliplr(np.transpose(Map))

def getMapExtents(BasePosition,windowSize):
    source = osr.SpatialReference()
    source.ImportFromEPSG(4326) # Coordinate system of DAWA results
  
    target = osr.SpatialReference()
    target.ImportFromEPSG(25832) # Coordinate system of maps
    
    transform = osr.CoordinateTransformation(source, target)
         
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(BasePosition.get('long'),BasePosition.get('lat'))
    point.Transform(transform)
            
    nMin=int(floor((point.GetX()-BasePosition.get('radius'))/1e3)) #
    nMax=int(floor((point.GetX()+BasePosition.get('radius'))/1e3)) # Determine relevant map files to load
    mMin=int(floor((point.GetY()-BasePosition.get('radius'))/1e3)) #
    mMax=int(floor((point.GetY()+BasePosition.get('radius'))/1e3)) #

    mapFiles = [[n,m] for n in range(nMin,nMax+1,windowSize) for m in range(mMin,mMax+1,windowSize)]
    return mapFiles

def saveKML(name, coordinateList, addressList, LoSList):
    kmlC = simplekml.Kml()
    kmlUC = simplekml.Kml()

    for m in range(0,len(LoSList)):
        addStr = addressList[m].get('addvejnavn')+' '+addressList[m].get('husnr')
        
        if (LoSList[m]):
            kmlC.newpoint(name=addStr, description="med LoS",coords=[(coordinateList[m][0],coordinateList[m][1])])  # lon, lat optional height
            #kmlC.newpoint(description="Daekket", coords=[(coordinateList[m][0],coordinateList[m][1])])  # lon, lat optional height
        else:
            kmlUC.newpoint(name=addStr, description="uden LoS",coords=[(coordinateList[m][0],coordinateList[m][1])])  # lon, lat optional height

    kmlC.save(name+'_Google_Earth_med_LoS.kml')
    kmlUC.save(name+'_Google_Earth_uden_LoS.kml')           
    return 0

def saveData(name, addressList, LoSList, BW, coveredBase):
    with open('indstillinger','rb') as f:
        name, basePosition, logon, mapRange = pickle.load(f)
        
    covered = open(name+'_med_LoS.csv', 'wb')
    notCovered = open(name+'_uden_LoS.csv', 'wb')

    keys = ['name','long','lat','radius','hbase', 'hclient', 'freq', 'download', 'upload', 'thetamin', 'thetamax']
    for n in range(0,len(basePosition)):
        BS = []
        for key in keys:
            BS.append(str(basePosition[n].get(key)))
            
        s = str(BS)
        covered.write(bytes(s, 'latin-1'))
        notCovered.write(bytes(s, 'latin-1'))


    covLevel = str(len([n for n in LoSList if n==1]))+'/'+str(len(LoSList))
    covered.write(bytes(covLevel, 'latin-1'))
    notCovered.write(bytes(covLevel, 'latin-1'))
    
    s='\nid;kvh;kommunenavn;kommunekode;postnr;vejnavn;vejkode;husnr;bogstav;download_tek_mbits;upload_tek_mbits;download_udt_privat_mbits;upload_udt_privat_mbits;download_udt_erhverv_mbits;upload_udt_erhverv_mbits;base(r)\n'
    covered.write(bytes(s, 'latin-1'))
    notCovered.write(bytes(s, 'latin-1'))
    
    
    for m in range(0,len(LoSList)):
        coveredBasesStr = ''
        for i3, cb in enumerate(coveredBase[m]):
            if i3 < len(coveredBase[m])-1:
                coveredBasesStr += str(cb)+', '
            else:
                coveredBasesStr += str(cb)

        addStr = addressList[m].get('uid')+';;'+addressList[m].get('kommunenavn')+';'+addressList[m].get('kommunekode')+';'\
        +addressList[m].get('postnr')+';'+addressList[m].get('addvejnavn')+';'+addressList[m].get('vejkode')+';'\
        +addressList[m].get('husnr')+';;'+str(BW[m][0])+';'+str(BW[m][1])+';'+str(BW[m][0])+';'+str(BW[m][1])+';'\
        +str(BW[m][0])+';'+str(BW[m][1])+';'+coveredBasesStr+'\n'
        if LoSList[m]:
            covered.write(bytes(addStr, 'latin-1'))
        else:
            notCovered.write(bytes(addStr, 'latin-1'))
            
            
    covered.close()
    notCovered.close()
    return 0

def findAddresses(BasePosition):
    DAWA_URL = 'http://dawa.aws.dk/adgangsadresser?cirkel='+str(BasePosition.get('long'))+','+str(BasePosition.get('lat'))+','+str(BasePosition.get('radius'))+'&format=geojson'
    print(DAWA_URL)
    r = requests.get(DAWA_URL)
    
    try:
        features = json.loads(r.text)['features']
    except KeyError:
        features = []
    coordinates = []
    uid = []
    address = []
    for i in features:
        coordinates.append(i['geometry']['coordinates'])
        uid.append(i['properties']['id'])
        address.append({
            'uid':i['properties']['id'],
            'kommunenavn':i['properties']['kommunenavn'],
            'kommunekode':i['properties']['kommunekode'],
            'postnr':i['properties']['postnr'],
            'addvejnavn':i['properties']['vejnavn'],
            'vejkode':i['properties']['vejkode'],
            'husnr':i['properties']['husnr']
            })

    return coordinates, address

def getBuildingCoords(address_uid, numElems=None):
    DAWA_URL = 'http://dawa.aws.dk/bygninger?adgangsadresseid='+str(address_uid)+'&format=geojson'
    idx = 0
    while idx < 5:
        try:
            r = requests.get(DAWA_URL, timeout=10)
            idx = 5
        except:
            print('BBR points error')
            print(str(address_uid))
            r = {}
            time.sleep(0.5)
            idx += 1
        
    try:
        features = json.loads(r.text)['features']
    except KeyError:
        features = []
        
    if len(features) == 0:
        return features
    
    arr1 = features[0]
    arr2 = arr1['geometry']
    arr3 = arr2['coordinates']
    arr4 = arr3[0]
    arr = arr4[0:-1]
    
    if numElems == None:
        return arr
    else:
        return shrink_list(arr, numElems) 

def shrink_list(arr, numElems):
    if len(arr)>numElems:
        idx = np.round(np.linspace(0, len(arr) - 1, numElems)).astype(int)
        return [arr[i] for i in idx]
    else:
        return arr

def downloadMaps(self, logon, BasePosition):
    print('Downloading maps')
    self.statusSig.sig.emit("Downloader kort")
    source = osr.SpatialReference()
    source.ImportFromEPSG(4326) # Coordinate system of DAWA results 
    target = osr.SpatialReference()
    target.ImportFromEPSG(25832) # Coordinate system of maps
    transform = osr.CoordinateTransformation(source, target)
        
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(BasePosition.get('long'),BasePosition.get('lat'))
    point.Transform(transform)    
    nZipMin=int(floor((point.GetX()-BasePosition.get('radius'))/1e4)) #
    nZipMax=int(floor((point.GetX()+BasePosition.get('radius'))/1e4)) # Determine relevant zip files to download
    mZipMin=int(floor((point.GetY()-BasePosition.get('radius'))/1e4)) #
    mZipMax=int(floor((point.GetY()+BasePosition.get('radius'))/1e4)) #
    nZipRange = list(range(nZipMin,nZipMax+1))
    mZipRange = list(range(mZipMin,mZipMax+1))

    loggedOn = 0
    
    if not os.path.isdir("Kort"):
        os.mkdir("Kort")

    for n in nZipRange:
        for m in mZipRange:
            terrainZipName = 'DTM_'+str(m)+'_'+str(n)+'_TIF_UTM32-ETRS89.zip'
            self.statusSig.sig.emit("Downloader kort - "+terrainZipName)
            if os.path.isfile('Kort/'+terrainZipName):
                Download = (time.time()-os.stat('Kort/'+terrainZipName).st_mtime)>(60*60*24*365)

            else:
                if not loggedOn:
                    ftp = FTP('ftp.kortforsyningen.dk',user=logon[0], passwd = logon[1])
                ftp.cwd('/dhm_danmarks_hoejdemodel/DTM')
                Download = (terrainZipName in ftp.nlst())

            if Download:
                # Download
                ftp.retrbinary('RETR '+terrainZipName, open('Kort/'+terrainZipName, 'wb').write,102400)
                # Unzip
                zf = zipfile.ZipFile('Kort/'+terrainZipName)
                zf.extractall('Kort')
                zf.close()
                os.remove('Kort/'+terrainZipName)
                open('Kort/'+terrainZipName, 'a').close()


            surfaceZipName = 'DSM_'+str(m)+'_'+str(n)+'_TIF_UTM32-ETRS89.zip'
                  
            if os.path.isfile('Kort/'+surfaceZipName):
                Download = (time.time()-os.stat('Kort/'+surfaceZipName).st_mtime)>(60*60*24*365)

            else:
                if not loggedOn:
                    ftp = FTP('ftp.kortforsyningen.dk',user=logon[0], passwd = logon[1])
                ftp.cwd('/dhm_danmarks_hoejdemodel/DSM')
                Download = (surfaceZipName in ftp.nlst())

            if Download:
                
                # Download
                ftp.retrbinary('RETR '+surfaceZipName, open('Kort/'+surfaceZipName, 'wb').write,102400)
                # Unzip
                zf = zipfile.ZipFile('Kort/'+surfaceZipName)
                zf.extractall('Kort')
                zf.close()
                os.remove('Kort/'+surfaceZipName)
                open('Kort/'+surfaceZipName, 'a').close()

    if loggedOn:
        ftp.quit()
    return 0

def findCoveredAddresses(BasePosition, inCoordinates, inAddress, clientHeights=None, BPHeight=None):
    # This function finds the addresses within range of the base position and calculates their line of sight.
    #
    # Base positions are given as a dictionary with the keys:
    #'name'
    #'long', 'lat' (WGS84)
    #'radius' (Maximum coverage radius, m)
    #'hbase' (Antenna height, m)
    #'hclient' (Client antenna height over ground, m)
    #'rhclient' (Client antenna height over roof, m)
    #'download' 
    #'upload'
    #'freq' (MHz)
    #'thetamin', 'thetamax' (For directional antennas; angles clockwise from North)
    # 
    # For each 1 km x 1 km it calculates all the sight lines that pass through, saves a png-image of the map, and emits a plot signal also containing
    # sample data points to display one of the sight lines with 1. Fresnel zone.
    # It also emits status signals with text strings to be displayed to the user.
    # in the end, it returns a list of address coordinates, address data, and LoS-status (0 or 1).
    #
    # The function is run once for each base station. Then, the coverage data should be aggregated and exported:
    #
    #coveredCoordinates = []
    #coveredUid = []
    #coveredAddress = []
    #coveredLoS = []
    #
    #for i in range(0,len(self.addressList)):
    #  if (self.addressList[i].get('uid') not in coveredUid):
    #     coveredCoordinates.append(coordinatesList[i])
    #     coveredUid.append(self.addressList[i].get('uid'))
    #     coveredAddress.append(self.addressList[i])
    #     coveredLoS.append(self.LoSList[i])
    #     
    #  elif (self.LoSList[i]):
    #      #coveredUid[coveredUid.index(addressList[i].get('uid'))] = LoSList[i]
    #      coveredLoS[coveredUid.index(self.addressList[i].get('uid'))] = self.LoSList[i]
    #
    #saveData(self.name,coveredAddress, coveredLoS, bwList)
    #saveKML(self.name,coveredCoordinates, coveredAddress, coveredLoS, self.basePosition[0])
    print('Calculating for '+str(len(inAddress))+' addresses')
    # Set parameters
    windowSize = 1                      # Size of calculation window (in no. of map tiles)
    inputSize = 2500                    # Side length of map tiles (in pixels)
    resolution = 0.4                    # Map resolution in m
    Kfactor = 1                         # 1 is a conservative value.
    Rearth = 6371e3*Kfactor             # Effective earth radius
    #WantedClearance = 1                 # Fraction of first Fresnel zone that needs to be clear.
    source = osr.SpatialReference()
    source.ImportFromEPSG(4326)         # Coordinate system of DAWA results
    target = osr.SpatialReference()
    target.ImportFromEPSG(25832)        # Coordinate system of maps
    
    # Find the necessary map file range
    mapFiles = getMapExtents(BasePosition,windowSize)
  
    # Show the map and sight lines
    nMin = min([n for [n,_] in mapFiles])
    nMax = max([n for [n,_] in mapFiles])
    mMin = min([m for [n,m] in mapFiles])
    mMax = max([m for [n,m] in mapFiles])
    
    imgDim = [nMax-nMin+1,mMax-mMin+1]
    imgDim = [(n+n%windowSize)*100 for n in imgDim]     
    mapMat = np.zeros(imgDim)
    posMat = np.zeros(imgDim)
    slMat = np.zeros(imgDim)
    
    wl = 3e2/BasePosition.get('freq') # c/f (0.06 m @ 5 GHz)
    
    transform = osr.CoordinateTransformation(source, target)
    
    basePos = ogr.Geometry(ogr.wkbPoint)
    basePos.AddPoint(BasePosition.get('long'),BasePosition.get('lat'))
    basePos.Transform(transform)
    BPCoords=[basePos.GetX()/1e3,basePos.GetY()/1e3]
                                           
    thetamin = BasePosition.get('thetamin')
    thetamax = BasePosition.get('thetamax')

    #####################################################################################################################################################
    # Calculate relevant parameters for each address, for later use
    coordinates = []
    CPCoords = []
    address = []
    angle = []
    if clientHeights == None:
        clientHeights = []
        clientHeights_stat = True
    else:
        clientHeights_stat = False
    totDist = []
    LoS = []
    
    for n in range(0,len(inAddress)):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(inCoordinates[n][0],inCoordinates[n][1])
        point.Transform(transform)
        inAngle = (pi+atan2(basePos.GetX()-point.GetX(),basePos.GetY()-point.GetY()))*180/pi
        
        distance = sqrt((basePos.GetX()-point.GetX())**2+(basePos.GetY()-point.GetY())**2)

        if (inAngle>thetamin)&(inAngle<thetamax) & (distance<=BasePosition.get('radius')):
            coordinates.append(inCoordinates[n])
            address.append(inAddress[n])            
            angle.append(inAngle)

            CPCoords.append([point.GetX()/1e3,point.GetY()/1e3])
            if clientHeights_stat == True:
                clientHeights.append(-999)
            totDist.append(sqrt((basePos.GetX()-point.GetX())**2+(basePos.GetY()-point.GetY())**2))
            LoS.append(1)
                 
    #####################################################################################################################################################
    terrainMap=np.zeros([inputSize,inputSize]) # Initialize array

    sSquare = 5
    maxDiff = 5
    heights=[]
    #####################################################################################################################################################
    
    if len(coordinates):
        
        #####################################################################################################################################################
        # Find client and base heights
        if clientHeights_stat == True:
            for [n,m] in mapFiles:
                CPInside = any([(N>=n) & (N<=n+windowSize) & (M>=m) & (M<=m+windowSize) for [N,M] in CPCoords])
                BPInside = any([(N>=n) & (N<=n+windowSize) & (M>=m) & (M<=m+windowSize) for [N,M] in [BPCoords]])
                
                if CPInside | BPInside: # Only load file if there are addresses or base stations inside...
                    print("Loaded "+ str([n,m]))
                    terrainMap = loadMapsHybrid(range(n,n+windowSize), range(m,m+windowSize), 't')
                    surfMap = loadMapsHybrid(range(n,n+windowSize), range(m,m+windowSize), 's')
    
                    for nClient in range(0,len(address)):
                        if ((CPCoords[nClient][0]>n) & (CPCoords[nClient][0]<n+windowSize) & (CPCoords[nClient][1]>m) & (CPCoords[nClient][1]<m+windowSize)):
                            CPMatrixCoords = [floor((CPCoords[nClient][0]-n)*inputSize), floor((CPCoords[nClient][1]-m)*inputSize)]
                            tHeight = terrainMap[CPMatrixCoords[0],CPMatrixCoords[1]] + BasePosition.get('hclient')
                            sHeight = surfMap[CPMatrixCoords[0],CPMatrixCoords[1]]
                            posMat[floor((CPCoords[nClient][0]-nMin)*100), floor((CPCoords[nClient][1]-mMin)*100)] = 255
    
                            heights=[]
                            if BasePosition.get('rhclient')>0:
                                NRange = range(max([0,CPMatrixCoords[0]-sSquare]),min([inputSize,CPMatrixCoords[0]+sSquare]))
                                MRange = range(max([0,CPMatrixCoords[1]-sSquare]),min([inputSize,CPMatrixCoords[1]+sSquare]))
                                heights.extend([surfMap[N,M]+BasePosition.get('rhclient') for N in NRange for M in MRange if surfMap[N,M] <= (sHeight+maxDiff)])
                                                        
                            heights.append(tHeight)
                            maxH = max(heights)
                            
                            clientHeights[nClient] = maxH
                       
                    if BPInside:
                        BPMatrixCoords = [floor((BPCoords[0]-n)*inputSize), floor((BPCoords[1]-m)*inputSize)]
                        BPHeight = terrainMap[BPMatrixCoords[0],BPMatrixCoords[1]] + BasePosition.get('hbase')
                        posMat[floor((BPCoords[0]-nMin)*100), floor((BPCoords[1]-mMin)*100)] = 255
    
                    terrainMap = None
                    surfMap = None
                    maxH = None
                    NRange = None
                    MRange = None
                    sHeight = None
                    tHeight = None
                    CPMatrixCoords = None
                    CPInside = None
                    BPInside = None
                    
        else:
            terrainMap = None
            surfMap = None
            maxH = None
            NRange = None
            MRange = None
            sHeight = None
            tHeight = None
            CPMatrixCoords = None
            CPInside = None
            BPInside = None
        #####################################################################################################################################################
        
        #####################################################################################################################################################
        # Do the LoS-calculations subsection by subsection
        hCentre = []
        hFresnel = []
        hTerrain = []
        sortedMapFiles=[]
        for offset in range(0,len(mapFiles)):
            sortedMapFiles.extend([[n,m] for [n,m] in mapFiles if ((n==nMin+offset)|(n==nMax-offset)|(m==mMin++offset)|(m==mMax-offset))&([n,m] not in sortedMapFiles)])
        
        for [n,m] in sortedMapFiles:
            dist = []
            hFresnel = []
            hTerrain = []
            hCentre = []
            sMapLoaded = 0 # Never load the same subsection twice
            
            cornerPoints = [[n,m],[n+windowSize,m],[n+windowSize,m+windowSize],[n,m+windowSize]]

            cornerMatrixCoords = [[floor((corner[0]-n)*inputSize), floor((corner[1]-m)*inputSize)] for corner in cornerPoints]

            BPMatrixCoords = [floor((BPCoords[0]-n)*inputSize), floor((BPCoords[1]-m)*inputSize)]

            for nClient in [n for n in range(0,len(address)) if LoS[n]]:
                CPMatrixCoords = [floor((CPCoords[nClient][0]-n)*inputSize), floor((CPCoords[nClient][1]-m)*inputSize)]
                    
                CPInside = any([(N>=n) & (N<n+windowSize) & (M>=m) & (M<m+windowSize) for [N,M] in [CPCoords[nClient]]])
                BPInside = any([(N>=n) & (N<n+windowSize) & (M>=m) & (M<m+windowSize) for [N,M] in [BPCoords]])

                #intersects = 0

                intPoints = []
                for corner in range(-1,3):
                    x00 = BPMatrixCoords[0]
                    y00 = BPMatrixCoords[1]

                    x10 = cornerMatrixCoords[corner][0]
                    y10 = cornerMatrixCoords[corner][1]

                    x01 = CPMatrixCoords[0]-BPMatrixCoords[0]
                    y01 = CPMatrixCoords[1]-BPMatrixCoords[1]

                    x11 = cornerMatrixCoords[corner+1][0]-cornerMatrixCoords[corner][0]
                    y11 = cornerMatrixCoords[corner+1][1]-cornerMatrixCoords[corner][1]

                    d = x11*y01 - x01*y11

                    if d!=0:
                        s = (1/d) * ((x00-x10)*y01-(y00-y10)*x01)
                        t = -(1/d) * (-(x00-x10)*y11+(y00-y10)*x11)
                    else:
                        s=-1
                        t=-1

                    if (s>=0) & (s<=1) & (t>=0) & (t<=1):
                        intPoints.append([x10+s*x11, y10+s*y11])
                
                if len(intPoints)>2:
                    print("Error: Sight line intercepts bounding box at more than 2 points!")
                    
                if (len(intPoints)>1) | CPInside | BPInside: # Load the map section if a sight line crosses its border, or if the client station is inside.
                    if not sMapLoaded:
                        surfaceMap = loadMapsHybrid(range(n,n+windowSize), range(m,m+windowSize), 's')
                        
                        print("Loaded "+ str([n,m]))
                        sMapLoaded = 1
                        mapMat[(floor(n-nMin)*100):(floor(n-nMin)*100+100*windowSize), (floor(m-mMin)*100):(floor(m-mMin)*100+100*windowSize)] = surfaceMap[::25,::25]

                    if len(intPoints)==2:
                        x0, y0 = intPoints[0]
                        x1, y1 = intPoints[1]

                    if CPInside & BPInside:
                        [x0, y0] = BPMatrixCoords
                        [x1, y1] = CPMatrixCoords

                    elif CPInside:
                        x0, y0 = intPoints[0]
                        [x1, y1] = CPMatrixCoords
                        
                    elif BPInside:
                        [x0, y0] = BPMatrixCoords
                        x1, y1 = intPoints[0]
                        
                    num = sqrt((x0-x1)**2+(y0-y1)**2) # Choose an appropriate number of nearest-neighbour sampling points
                    x, y = np.linspace(x0, x1, num), np.linspace(y0, y1, num)
                    x=x.astype(np.int)
                    y=y.astype(np.int)

                    xn = [x[n] for n in range(0,len(x)) if (  (x[n]>0) & (x[n]<(windowSize*inputSize)) & (y[n]>0) & (y[n]<(windowSize*inputSize))   )]
                    yn = [y[n] for n in range(0,len(x)) if (  (x[n]>0) & (x[n]<(windowSize*inputSize)) & (y[n]>0) & (y[n]<(windowSize*inputSize))   )]

                    x = xn
                    y = yn

                    zi = surfaceMap[x, y]

                    dist = [sqrt(   (x[n]-CPMatrixCoords[0])**2 + (y[n]-CPMatrixCoords[1])**2)  *resolution for n in range(0,len(x))]
                    dist = [dist[n] for n in range(0,len(dist)) if ((totDist[nClient]>dist[n])&(dist[n]>0))]
                    
                    earthCurvature = [Rearth*(np.cos((distance-totDist[nClient]/2)/Rearth)-np.cos(totDist[nClient]/(2*Rearth))) for distance in dist]
                    earthCurvature = [ec-Rearth*(np.cos((-totDist[nClient]/2)/Rearth)-np.cos(totDist[nClient]/(2*Rearth))) for ec in earthCurvature]
                    
                    hTerrain = [zi[n]+earthCurvature[n] for n in range(0,len(dist))]
                    
                    hCentre = [clientHeights[nClient]+(BPHeight-clientHeights[nClient])*distance/totDist[nClient] for distance in dist]

                    hFresnel = [hCentre[n] - sqrt(  wl*dist[n]*(totDist[nClient]-dist[n])/totDist[nClient]  ) for n in range(0,len(dist))]

                    xp, yp = np.linspace(x0/25+(n-nMin)*100, x1/25+(n-nMin)*100, ceil(num/25)), np.linspace(y0/25+(m-mMin)*100, y1/25+(m-mMin)*100, ceil(num/25))
                    xp=xp.astype(np.int)
                    yp=yp.astype(np.int)
                    xp1 = [xp[n] for n in range(0,len(xp)) if (  (xp[n]>0) & (xp[n]<len(posMat)) & (yp[n]>0) & (yp[n]<len(posMat))   )]
                    yp1 = [yp[n] for n in range(0,len(xp)) if (  (xp[n]>0) & (xp[n]<len(posMat)) & (yp[n]>0) & (yp[n]<len(posMat))   )]
                    slMat[xp1,yp1] = np.ones(len(xp1))
                    if any([hFresnel[n]<hTerrain[n] for n in range(0,len(hFresnel))]):
                        LoS[nClient] = 0
                        if ((dist[0]<2.0)&(hFresnel[0]<hTerrain[0]))|((dist[-1]<2.0)&(hFresnel[-1]<hTerrain[-1])):
                            print("Warning: Client below terrain!")
                
    print('Clients processed:')
    print(len(address))
    
    return coordinates, address, LoS, clientHeights, BPHeight


if __name__=='__main__':
    
    arr = [1,2,3,4,5,6,7]
    a = shrink_list(arr, 6)
    print(a)