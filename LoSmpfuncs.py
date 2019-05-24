import LoSfuncs as lf

class findCoveredMP():
    def __init__(self):
        print('Multiprocessing module loaded')
        
    def findBase(self, BP):
        inCoordinates, inAddress = lf.findAddresses(BP)
        
        coordinates, address, LoS, clientHeights, BPHeight = lf.findCoveredAddresses(BP, inCoordinates, inAddress)
        
        inCoordinatesT3 = list(range(0,len(inCoordinates)))
    
        LoST_tf = list(range(0,len(LoS)))
        
        for i5, add0 in enumerate(address):
            if LoS[i5] == 1:
                LoST_tf[i5] = True
                inCoordinatesT3[i5] = []
            else:
                LoST_tf[i5] = False
                inCoordinatesT3[i5] = lf.getBuildingCoords(address_uid=add0['uid'], numElems=6)
        
        while False in LoST_tf:
            inCoordT = []
            inAddressT = []
            clientHeightsT = []
            
            for i5, LoST in enumerate(LoS):
                if LoST == 0:
                    if len(inCoordinatesT3[i5])==0:
                        LoST_tf[i5] = None
                        
                    else:
                        tt = inCoordinatesT3[i5].pop(0)
                        inCoordT += [tt]
                        inAddressT += [address[i5]]
                        clientHeightsT += [clientHeights[i5]]
            
            coordinatesT2, addressT2, LoST2, _, _ = lf.findCoveredAddresses(BP, inCoordT, inAddressT, clientHeightsT, BPHeight)
            
            for i5, add0 in enumerate(address):
                for i6, add1 in enumerate(addressT2):
                    if add0['uid']==add1['uid']:
                        if LoST2[i6] == 1:
                            LoST_tf[i5] = 'Found'
                        
                            LoS[i5] = LoST2.copy()[i6]
                            coordinates[i5] = coordinatesT2.copy()[i6]
                            
        
        BW = [[BP.get('download'),BP.get('upload')] for _ in LoS]
        
        data_dic= {
            'name':BP.get('name'),
            'address':address.copy(),
            'LoS':LoS.copy(),
            'coordinates':coordinates.copy(),
            'BW':BW.copy(),
            }
        
        return data_dic
        