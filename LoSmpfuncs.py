import LoSfuncs as lf
from multiprocessing.dummy import Pool as ThreadPool
#import pydevd
import timtoc
import time
import pickle
import os

class findCoveredMP():
    def __init__(self):
        print('Multiprocessing module loaded')
        
    def findBase(self, BP):
        tstart = time.time()
        
        inCoordinates, inAddress = lf.findAddresses(BP)
        
        #inCoordinates = [[11.29590515, 55.50009737]]
        #inAddress = [{'uid': '0a3f5083-6e7d-32b8-e044-0003ba298018', 'kommunenavn': 'Kalundborg', 'kommunekode': '0326', 'postnr': '4270', 'addvejnavn': 'Tovesvej', 'vejkode': '1840', 'husnr': '41'}]

        with timtoc.Timer('First'):
            coordinates, address, LoS, clientHeights, BPHeight = lf.findCoveredAddresses(BP, inCoordinates, inAddress)
        
        get_uids = []
        for idx0, sight in enumerate(LoS):
            if sight == 0:
                get_uids += [address[idx0]['uid']]
        
        with timtoc.Timer('Get BBR points'):
            pool_url = ThreadPool(16)
            get_uids_res = pool_url.map(lf.getBuildingCoords, get_uids)
            pool_url.close()
            pool_url.join()
        
            inCoordinatesT3 = list(range(0,len(inCoordinates)))
            LoST_tf = list(range(0,len(LoS)))
        
        
            pc = 0
            for i5, add0 in enumerate(address):
                if LoS[i5] == 1:
                    LoST_tf[i5] = True
                    inCoordinatesT3[i5] = []
                else:
                    #if (pc == 1000):
                    #    print('Progress report: ',i5)
                    #    pc = 0
                    LoST_tf[i5] = False
                    #inCoordinatesT3[i5] = lf.getBuildingCoords(address_uid=add0['uid'])
                    inCoordinatesT3[i5] = get_uids_res[pc]
                    pc += 1
        
        #pool_url = ThreadPool(4)
        #results = pool_url.map(lf.getBuildingCoords, urls)
        #pool_url.close()
        #pool_url.join()
        
        
        while False in LoST_tf:
            inCoordT = []
            inAddressT = []
            clientHeightsT = {}
            
            for i5, LoST in enumerate(LoS):
                if LoST == 0:
                    if len(inCoordinatesT3[i5])==0:
                        LoST_tf[i5] = None
                        
                    else:
                        tt = inCoordinatesT3[i5].pop(0)
                        inCoordT += [tt]
                        inAddressT += [address[i5]]
                        clientHeightsT[address[i5]['uid']] = clientHeights[address[i5]['uid']]
            
            #pydevd.settrace(suspend=False, trace_only_current_thread=True)
            #
            #tstart = time.time()
            if len(inCoordT)!=0:
                with timtoc.Timer('Second'):
                    coordinatesT2, addressT2, LoST2, _, _ = lf.findCoveredAddresses(BP, inCoordT, inAddressT, clientHeightsT, BPHeight)
                #print('[Second] Elapsed: %s' % ((time.time() - tstart)/len(addressT2)))
            
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
        
        if os.path.exists('pickles/'+BP['name']+'.pickle'):
            os.remove('pickles/'+BP['name']+'.pickle')
        
        with open('pickles/'+BP['name']+'.pickle', 'wb') as handle:
            pickle.dump(data_dic, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        print('[MP] Elapsed: %s' % (time.time() - tstart))
        return data_dic
        