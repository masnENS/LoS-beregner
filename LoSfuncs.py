import requests, json
import numpy as np
from ftplib import FTP
from osgeo import ogr
from osgeo import osr
import os.path
import zipfile
import time
from math import floor, atan2, pi, sqrt
import pickle
import simplekml
from osgeo import gdal
import urllib.request
from _collections import OrderedDict
import socket
import AllMapsfuncs as allmap
#import pydevd


def loadMapsWCS(bbox, dim, mapType, basename):
    mapFile = 'tmp_tifs/'+basename+'.tif'
    
    if os.path.isfile(mapFile):
        m = gdal.Open(mapFile)
        dtmMap = np.array(m.GetRasterBand(1).ReadAsArray())
        print('Loaded from tmp_tifs: ',basename)
        return dtmMap
    
    BBOX = "BBOX="+str(bbox[0])+","+str(bbox[1])+","+str(bbox[2])+","+str(bbox[3])
    DIM = "WIDTH="+str(dim[0])+"&HEIGHT="+str(dim[1])
    
    if mapType == 't':
        coverage = "dhm_terraen"
    elif mapType == 's':
        coverage = "dhm_overflade"
    
    url = "http://services.kortforsyningen.dk/dhm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=1.0.0&COVERAGE="+coverage+"&CRS=EPSG:25832&"+BBOX+"&"+DIM+"&FORMAT=image/gtiff&login=ASN&password=ASN"
    print(url)
    
    getre = False
    while getre == False:
        try:
            urllib.request.urlretrieve(url,mapFile)
            getre = True
        except:
            conn = False
            while conn == False:
                time.sleep(10)
                conn = check_internet()
                if conn == False:
                    print('No internet connection!')
                else:
                    print('Connection established')
    
    m = gdal.Open(mapFile)
    dtmMap = np.array(m.GetRasterBand(1).ReadAsArray())
    
    return dtmMap

def loadMapsHybrid(nRange, mRange, mapType, basename):
    inputSize = 2500
    scaledSize = inputSize

    if mapType == 't':
        prefix = 'AllMaps_DTM_npz/DTM_1km_'
        
    elif mapType == 's':
        prefix = 'AllMaps_DSM_npz/DSM_1km_'

    Map = np.zeros([scaledSize*len(mRange),scaledSize*(len(nRange))])
    Map = Map.astype('float32')
    
    for n in nRange:
        #print(n)
        for m in mRange:
             
            nStart = (n-nRange[0])*inputSize
            mStart = (mRange[-1]-m)*inputSize

            nEnd = nStart+inputSize
            mEnd = mStart+inputSize
       
            fileName = prefix+str(m)+'_'+str(n)+'.npz'
            if os.path.isfile(fileName):
                #dt = gdal.Open(fileName)
                #Map[mStart:mEnd,nStart:nEnd] = np.array(dt.GetRasterBand(1).ReadAsArray(),dtype='float32')
                #dt=None
                with np.load(fileName) as data:
                    Map[mStart:mEnd,nStart:nEnd] = data['a'].astype('float32')/10.
                
            else:
                #BBOX = "BBOX="+str(858275)+","+str(6108438)+","+str(897160)+","+str(6144049)
                bbox = [n*1000,m*1000,(n+1)*1000,(m+1)*1000]
                print('Filename: ',fileName)
                print('Downloaded temp map: ',(n,m))
                Map[mStart:mEnd,nStart:nEnd] = loadMapsWCS(bbox, [2500,2500], mapType, basename+'_'+str((n,m)))

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
    
    idx = 0
    while idx < 5:
        try:
            #r = requests.get(DAWA_URL)
            with requests.Session() as s:
                r = s.get(DAWA_URL, timeout=10)
            
            #r = requests.get(DAWA_URL, timeout=10)
            
            try:
                r_txt = r.text
                json.loads(r_txt)
                idx = 5
            except:
                #print('Error in loading JSON: ',address_uid)
                r_txt = {}
                r = None
                time.sleep(np.random.uniform(3,10))
                idx += 1
                if idx == 5:
                    print('Error in loading JSON for BasePosition: ',BasePosition)
            
        except:
            conn = False
            while conn == False:
                time.sleep(10)
                conn = check_internet()
                if conn == False:
                    print('No internet connection!')
                else:
                    print('Connection established')
    
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

def getBuildingCoords(address_uid):
    numElems = 6
    DAWA_URL = 'http://dawa.aws.dk/bygninger?adgangsadresseid='+str(address_uid)+'&format=geojson'
    idx = 0
    while idx < 5:
        try:
            with requests.Session() as s:
                r = s.get(DAWA_URL, timeout=10)
            
            #r = requests.get(DAWA_URL, timeout=10)
            
            try:
                r_txt = r.text
                json.loads(r_txt)
                idx = 5
            except:
                #print('Error in loading JSON: ',address_uid)
                r_txt = {}
                r = None
                time.sleep(np.random.uniform(3,10))
                idx += 1
                if idx == 5:
                    print('Error in loading JSON: ',address_uid)
        except:
            print('BBR points error')
            print(str(address_uid))
            conn = False
            while conn == False:
                time.sleep(10)
                conn = check_internet()
                if conn == False:
                    print('No internet connection!')
                else:
                    print('Connection established')       
    
    try:
        features = json.loads(r_txt)['features']
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
    self.statusSig.sig.emit("Downloader kort: Igang")
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
    
    #print('nMin: ',nZipMin)
    #print('nMax: ',nZipMax)
    #print('mMin: ',mZipMin)
    #print('nMax: ',mZipMax)
    
    loggedOn = 0
    
    if not os.path.isdir("Kort"):
        os.mkdir("Kort")

    for n in nZipRange:
        for m in mZipRange:
            terrainZipName = 'DTM_'+str(m)+'_'+str(n)+'_TIF_UTM32-ETRS89.zip'
            
            if os.path.isfile('AllMaps_DTM/'+terrainZipName):
                Download = (time.time()-os.stat('AllMaps_DTM/'+terrainZipName).st_mtime)>(60*60*24*365)

            else:
                print('Logon FTP')
                if not loggedOn:
                    ftp = FTP('ftp.kortforsyningen.dk',user=logon[0], passwd = logon[1])
                ftp.cwd('/dhm_danmarks_hoejdemodel/DTM')
                Download = (terrainZipName in ftp.nlst())
                print('Down: ',Download)

            if Download:
                # Download
                print('Downloading T: ',(n,m))
                ftp.retrbinary('RETR '+terrainZipName, open('Kort/'+terrainZipName, 'wb').write,102400)
                # Unzip
                zf = zipfile.ZipFile('Kort/'+terrainZipName)
                zf.extractall('Kort')
                zf.close()
                os.remove('Kort/'+terrainZipName)
                open('Kort/'+terrainZipName, 'a').close()


            surfaceZipName = 'DSM_'+str(m)+'_'+str(n)+'_TIF_UTM32-ETRS89.zip'
                  
            if os.path.isfile('AllMaps_DSM/'+surfaceZipName):
                Download = (time.time()-os.stat('AllMaps_DSM/'+surfaceZipName).st_mtime)>(60*60*24*365)

            else:
                print('Logon FTP2')
                if not loggedOn:
                    ftp = FTP('ftp.kortforsyningen.dk',user=logon[0], passwd = logon[1])
                ftp.cwd('/dhm_danmarks_hoejdemodel/DSM')
                Download = (surfaceZipName in ftp.nlst())
                print('Down2: ',Download)

            if Download:
                print('Downloading T: ',(n,m))
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
    
    self.statusSig.sig.emit("Downloader kort: Færdig")    
    
    return 0

def findCoveredAddresses(BasePosition, inCoordinates, inAddress, clientHeights=None,  BPHeight=None):
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
    
    #pydevd.settrace(suspend=False, trace_only_current_thread=True)    
    
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
    #mapFiles = getMapExtents(BasePosition,windowSize)
  
    # Show the map and sight lines
    #nMin = min([n for [n,_] in mapFiles])
    #nMax = max([n for [n,_] in mapFiles])
    #mMin = min([m for [n,m] in mapFiles])
    #mMax = max([m for [n,m] in mapFiles])
    
    wl = 3e2/BasePosition.get('freq') # c/f (0.06 m @ 5 GHz)
    
    transform = osr.CoordinateTransformation(source, target)
    
    basePos = ogr.Geometry(ogr.wkbPoint)
    basePos.AddPoint(BasePosition.get('long'),BasePosition.get('lat'))
    basePos.Transform(transform)
    BPCoords=[basePos.GetX()/1e3,basePos.GetY()/1e3]
    
    #print('Base coords: ',BPCoords)
                                  
    thetamin = BasePosition.get('thetamin')
    thetamax = BasePosition.get('thetamax')

    #####################################################################################################################################################
    # Calculate relevant parameters for each address, for later use
    coordinates = []
    CPCoords = []
    address = []
    angle = []
    
    
    
    if BPHeight == None:
        clientHeights = {}
        clientHeights_stat = True
    else:
        clientHeights_stat = False
    totDist = []
    LoS = []
    
    for n, ia in enumerate(inAddress):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(inCoordinates[n][0],inCoordinates[n][1])
        point.Transform(transform)
        inAngle = (pi+atan2(basePos.GetX()-point.GetX(),basePos.GetY()-point.GetY()))*180/pi
        
        distance = sqrt((basePos.GetX()-point.GetX())**2+(basePos.GetY()-point.GetY())**2)
        
        
        insideAngle = False
        if distance <= BasePosition.get('radius'):
            if thetamax > thetamin:
                if inAngle > thetamin and inAngle < thetamax:
                    insideAngle = True
            
            if thetamin > thetamax:                
                if inAngle < thetamax and inAngle >= 0:
                    insideAngle = True
                elif inAngle > thetamin:
                    insideAngle = True
        
        if insideAngle:
        #if (inAngle>=thetamin)&(inAngle<thetamax) & (distance<=BasePosition.get('radius')):
            coordinates.append(inCoordinates[n])
            address.append(ia) ####
            angle.append(inAngle)

            CPCoords.append([point.GetX()/1e3,point.GetY()/1e3])
            
            if clientHeights_stat == True:
                clientHeights[ia['uid']] = None
            totDist.append(sqrt((basePos.GetX()-point.GetX())**2+(basePos.GetY()-point.GetY())**2))
            LoS.append(1)
                 
    #####################################################################################################################################################
    terrainMap=np.zeros([inputSize,inputSize]) # Initialize array

    sSquare = 5
    maxDiff = 5
    #heights=[]
    #####################################################################################################################################################
    
    mapLoadCounter = 0
    
    if len(coordinates):
        
        #####################################################################################################################################################
        # Find client and base heights
        if clientHeights_stat == True:
            n = floor(BPCoords[0])
            m = floor(BPCoords[1])
            terrainMap = loadMapsHybrid(range(n,n+windowSize), range(m,m+windowSize), 't', BasePosition['name'])
            BPMatrixCoords = [floor((BPCoords[0]-n)*inputSize), floor((BPCoords[1]-m)*inputSize)]
            BPHeight = terrainMap[BPMatrixCoords[0],BPMatrixCoords[1]] + BasePosition.get('hbase')
            
            terrainMap = None
            
        #####################################################################################################################################################
        
        #####################################################################################################################################################
        # Do the LoS-calculations subsection by subsection
        
        map_dict = findSquares(BPCoords, CPCoords, address)
        sortedMapFiles = sortDictDist(map_dict, BPCoords)
        
        addr_dict = {}        
        for idx0,addr in enumerate(address):
            addr_dict[addr['uid']] = idx0 
        
        
        
        for [n,m] in sortedMapFiles:
            
            dist = []
            hFresnel = []
            hTerrain = []
            hCentre = []
            #sMapLoaded = 0 # Never load the same subsection twice
                        
            uids = sortedMapFiles[(n,m)]
            
            tloaded = 0
            sloaded = 0
            
            for uid in uids:
                nClient = addr_dict[uid]
                
                if not LoS[nClient]:
                    continue
                
                if sloaded == 0:
                    cornerPoints = [[n,m],[n+windowSize,m],[n+windowSize,m+windowSize],[n,m+windowSize]]
                    cornerMatrixCoords = [[floor((corner[0]-n)*inputSize), floor((corner[1]-m)*inputSize)] for corner in cornerPoints]
                    BPMatrixCoords = [floor((BPCoords[0]-n)*inputSize), floor((BPCoords[1]-m)*inputSize)]  
                    
                    surfaceMap = loadMapsHybrid(range(n,n+windowSize), range(m,m+windowSize), 's', BasePosition['name'])
                    #print('Map section: ',(n,m))
                    mapLoadCounter += 1
                    sloaded = 1
                    
                CPMatrixCoords = [floor((CPCoords[nClient][0]-n)*inputSize), floor((CPCoords[nClient][1]-m)*inputSize)]
                
                CPInside = any([(N>=n) & (N<n+windowSize) & (M>=m) & (M<m+windowSize) for [N,M] in [CPCoords[nClient]]])
                BPInside = any([(N>=n) & (N<n+windowSize) & (M>=m) & (M<m+windowSize) for [N,M] in [BPCoords]])
                
                if (CPInside==True): #and (clientHeights[uid]==None) and (clientHeights_stat == True):
                    #CPMatrixCoords = [floor((CPCoords[nClient][0]-n)*inputSize), floor((CPCoords[nClient][1]-m)*inputSize)]
                    #print('TRUE: ',nClient)
                    
                    if tloaded == 0:
                        terrainMap = loadMapsHybrid(range(n,n+windowSize), range(m,m+windowSize), 't', BasePosition['name'])
                        tloaded = 1
                    
                    tHeight = terrainMap[CPMatrixCoords[0],CPMatrixCoords[1]] + BasePosition.get('hclient')
                    sHeight = surfaceMap[CPMatrixCoords[0],CPMatrixCoords[1]]
                    #posMat[floor((CPCoords[nClient][0]-nMin)*100), floor((CPCoords[nClient][1]-mMin)*100)] = 255

                    heights=[]
                    if BasePosition.get('rhclient')>0:
                        NRange = range(max([0,CPMatrixCoords[0]-sSquare]),min([inputSize,CPMatrixCoords[0]+sSquare]))
                        MRange = range(max([0,CPMatrixCoords[1]-sSquare]),min([inputSize,CPMatrixCoords[1]+sSquare]))
                        heights.extend([surfaceMap[N,M]+BasePosition.get('rhclient') for N in NRange for M in MRange if surfaceMap[N,M] <= (sHeight+maxDiff)])
                                                
                    heights.append(tHeight)
                    maxH = max(heights)
                    
                    clientHeights[uid] = maxH
                
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
                    
                    num = int(np.round(sqrt((x0-x1)**2+(y0-y1)**2))) # Choose an appropriate number of nearest-neighbour sampling points
                    x = np.linspace(x0, x1, num, dtype=np.int)
                    y = np.linspace(y0, y1, num, dtype=np.int)
                    
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
                    hCentre = [clientHeights[uid]+(BPHeight-clientHeights[uid])*distance/totDist[nClient] for distance in dist]
                    hFresnel = [hCentre[n] - sqrt(  wl*dist[n]*(totDist[nClient]-dist[n])/totDist[nClient]  ) for n in range(0,len(dist))]
                    
                                        
                    if any([hFresnel[n]<hTerrain[n] for n in range(0,len(hFresnel))]):
                        LoS[nClient] = 0
                    
    print('Addresses with LoS: ',np.sum(LoS))
    print('Map loaded counter: ',mapLoadCounter)
    #if clientHeights_stat == False:
    #    print('a')
    print('Clients processed:')
    print(len(address))
    
    
    return coordinates, address, LoS, clientHeights, BPHeight

def findSquares(BPCoords,CPCoords,inAddress):
    map_dict = {}
    BPMatrixCoords = (int(BPCoords[0]), int(BPCoords[1]))
    
    for idx,CPCoord in enumerate(CPCoords):    
        CPMatrixCoords = ( int(CPCoord[0]), int(CPCoord[1]) )
        
        uid = inAddress[idx]['uid']
        ####################
        
        if (floor(CPCoord[0]),floor(CPCoord[1])) not in map_dict:
            map_dict[(floor(CPCoord[0]),floor(CPCoord[1]))] = []
        
        map_dict[(floor(CPCoord[0]),floor(CPCoord[1]))].append(uid)
        
        angle = np.arctan2(CPCoord[1]-BPCoords[1],CPCoord[0]-BPCoords[0])
        
        y_ax = np.linspace(CPMatrixCoords[1],BPMatrixCoords[1],int(np.abs(BPMatrixCoords[1]-CPMatrixCoords[1]))+1)
        x_ip = np.repeat(np.floor( (y_ax-BPCoords[1])*(1/np.tan(angle)) + BPCoords[0] ),2)
        
        y_ax = np.repeat(y_ax,2)-[1,0]*int(len(y_ax))
        
        x_ax = np.linspace(CPMatrixCoords[0],BPMatrixCoords[0],int(np.abs(BPMatrixCoords[0]-CPMatrixCoords[0]))+1)
        y_ip = np.repeat(np.floor( (x_ax-BPCoords[0])*(np.tan(angle)) + BPCoords[1] ),2)
        
        x_ax = np.repeat(x_ax,2)-[1,0]*int(len(x_ax))
        
        for idx0,mt in enumerate(y_ax):
            n = int(x_ip[idx0])
            m = int(mt)
            pmx = CPCoord[0]-BPCoords[0]
            pmy = CPCoord[1]-BPCoords[1]
            
            if (pmy >= 0 and (CPMatrixCoords[1] >= m >= BPMatrixCoords[1])) or (pmy < 0 and (CPMatrixCoords[1] <= m <= BPMatrixCoords[1])):
                if (pmx >= 0 and (CPMatrixCoords[0] >= n >= BPMatrixCoords[0])) or (pmx < 0 and (CPMatrixCoords[0] <= n <= BPMatrixCoords[0])):
                    if (n,m) not in map_dict:
                        map_dict[(n,m)] = []
                    
                    if uid not in map_dict[(n,m)]: 
                        map_dict[(n,m)] += [uid]
        
        for idx0,nt in enumerate(x_ax):
            n = int(nt)
            m = int(y_ip[idx0])
            pmx = CPCoord[0]-BPCoords[0]
            pmy = CPCoord[1]-BPCoords[1]
            
            if (pmx >= 0 and (CPMatrixCoords[0] >= n >= BPMatrixCoords[0])) or (pmx < 0 and (CPMatrixCoords[0] <= n <= BPMatrixCoords[0])):
                if (pmy >= 0 and (CPMatrixCoords[1] >= m >= BPMatrixCoords[1])) or (pmy < 0 and (CPMatrixCoords[1] <= m <= BPMatrixCoords[1])):
                    if (n,m) not in map_dict:
                        map_dict[(n,m)] = []
                
                    if uid not in map_dict[(n,m)]: 
                        map_dict[(n,m)] += [uid]
    return map_dict

def sortDictDist(map_dict, BPCoords):
    dist_dict = {}
    for key in map_dict:
        
        dist = []
        for p in [(0,0),(1,0),(0,1),(1,1)]:
            dist += [sqrt( (key[0]+p[0]-BPCoords[0])**2 + (key[1]+p[1]-BPCoords[1])**2 )]
        
        dist_dict[key] = max(dist)
    
    sort_dict = OrderedDict(sorted(dist_dict.items(), key=lambda t: t[1],reverse=True))
    
    for key in map_dict:
        sort_dict[key] = map_dict[key]
        
    return sort_dict

def check_internet(host="8.8.8.8", port=53, timeout=3):
    """
    Host: 8.8.8.8 (google-public-dns-a.google.com)
    OpenPort: 53/tcp
    Service: domain (DNS/TCP)
    """
    try:
        socket.setdefaulttimeout(timeout)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect((host, port))
        return True
    except socket.error as ex:
        print(ex)
        return False

if __name__=='__main__':
    
    #arr = [1,2,3,4,5,6,7]
    #a = shrink_list(arr, 6)
    #print(a)
    
    #print(getBuildingCoords('0a3f5084-533e-32b8-e044-0003ba298018'))
    
    #DAWA_URL = 'http://services.kortforsyningen.dk/dhm?REQUEST=GetCoverage&SERVICE=WCS&VERSION=1.0.0&COVERAGE=dhm_overflade&CRS=EPSG:25832&BBOX=635000,6138000,636000,6139000&WIDTH=2500&HEIGHT=2500&FORMAT=image/gtiff&login=ASN&password=ASN'
    
    #r = requests.get(DAWA_URL)
    
    #print(r.text)
    print(check_internet())