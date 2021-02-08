import hashlib
import os
from ftplib import FTP
import zipfile
#from PIL import Image
import numpy as np
from osgeo import gdal
from multiprocessing.dummy import Pool as ThreadPool
from pathlib import Path
from multiprocessing import Pool

def round_up(n, decimals=0):
    multiplier = 10 ** decimals
    return np.ceil(n * multiplier) / multiplier

def check_md5(folder, fname):
    # Get original md5 from reference file
    with open(folder+'/'+fname+'.md5','r') as f:
        original_md5 = f.readline().strip()
    
    # Open,close, read file and calculate MD5 on its contents 
    with open(folder+'/'+fname+'.tif','rb') as file_to_check:
        # read contents of the file
        data = file_to_check.read()    
        # pipe contents of the file through
        md5_returned = hashlib.md5(data).hexdigest()
    
    # Finally compare original MD5 with freshly calculated
    if original_md5 == md5_returned:
        print('MD5 verified.')
        return True
    else:
        print('MD5 verification failed!.')
        return False

def check_md5s(folder, fnames):
    md5sgood = []
    for fname in fnames:
        md5sgood += [check_md5(folder, fname)]
        
    if False in md5sgood:
        return False
    else:
        return True

def get_file_names(folder, ftype):
    tifnames = []
    for file in os.listdir(folder):
        if file.endswith('.'+str(ftype)):
            tifnames += [file[:-(len(ftype)+1)]]
    return tifnames

class sdfeKort():
    def __init__(self):
        print('Attempting to log on to FTP server')
        self.ftp = FTP('ftp.kortforsyningen.dk',user='ens_dahj', passwd = 'EnergiKEFM19')
        print('Logon succeeded')
    
    def get_DTMzipnames(self):
        self.ftp.cwd('/dhm_danmarks_hoejdemodel/DTM')
        DTMrange = self.ftp.nlst()
        
        for i2, name in enumerate(DTMrange):
            if len(name) > 4:
                if name[-4::] != '.zip':
                    DTMrange.pop(i2)
            else:
                DTMrange.pop(i2)
        return DTMrange
        
    def get_DSMzipnames(self):
        self.ftp.cwd('/dhm_danmarks_hoejdemodel/DSM')
        DSMrange = self.ftp.nlst()
        
        for i2, name in enumerate(DSMrange):
            if len(name) > 4:
                if name[-4::] != '.zip':
                    DSMrange.pop(i2)
            else:
                DSMrange.pop(i2)
        return DSMrange
        
    def get_DSMfiles(self, fnames):
        
        #print('Attempting to log on to FTP server')
        ftp = FTP('ftp.kortforsyningen.dk',user='ens_dahj', passwd = 'EnergiKEFM19')
        #print('Logon succeeded')
        typ = 'S'
        
        folder = 'AllMaps_DSM'
        #save_folder = 'AllMapsNpz'
        options = {}
        options['unzip'] = True
        options['zip_delete'] = True
        #options['numpy'] = True
        
        if type(fnames) is str:
            fnames = [fnames]
        
        if not os.path.isdir(folder):
            os.makedirs(folder)
        
        ftp.cwd('/dhm_danmarks_hoejdemodel/D'+typ+'M')
        for fname in fnames:
            
            if fname[0:3] != ('D'+typ+'M'):
                print('Ikke rigtig: ',fname)
                continue
            
            if os.path.exists(folder+'/'+fname):
                
                #print('Size: ',Path(folder+'/'+fname).stat().st_size)
                
                
                if Path(folder+'/'+fname).stat().st_size < 1:
                
                    continue
            
            #print('Nope: ',fname)
            
            print('Downloading ' + fname)
            
            
            ftp.retrbinary('RETR '+fname, open(folder+'/'+fname, 'wb').write,102400)
            # Unzip
            if options.get('unzip') == True:
                zf = zipfile.ZipFile(folder+'/'+fname)
                zf.extractall(folder)
                zf.close()
                
                if options.get('zip_delete') == True:
                    os.remove(folder+'/'+fname)
                    open(folder+'/'+fname, 'a').close()
            
            # Convert to numpy
            #if options.get('numpy') == True:
            #    conv_tif_np16int(folder, save_folder, options.get('tm_delete'))
            
        ftp.close()
                
def conv_tif_np16int(tname):
    load_folder = 'AllMaps_DSM'
    save_folder = 'AllMaps_DSM_npz'
    tm_delete = False
    print('Converting: ',tname)
    
    dt = gdal.Open(load_folder+'/'+tname+'.tif')
    ima_np = np.array(dt.GetRasterBand(1).ReadAsArray(),dtype='float32')
    
    #ima = Image.open(load_folder+'/'+tname+'.tif')
    
    #print(ima)
    
    #ima_np = np.array(ima)
    
    #print(ima_np)
    
    ima_np16 = np.array( round_up(ima_np, decimals=1)*10 , dtype='int16')
    
    np.savez_compressed(save_folder+'/'+tname+'.npz', a=ima_np16)
    
    if tm_delete == True:
        os.remove(load_folder+'/'+tname+'.md5')
        os.remove(load_folder+'/'+tname+'.tif')
        
    return None

if __name__ == "__main__":
    
    '''
    n = np.array([1.14,2.14,3.14,4.14,5.14,6.4])
    a = round_up(n, decimals=1)
    print(a)
    b = np.array( round_up(a, decimals=1)*10 , np.dtype='int16')
    print(b)
    print(b/10)
    '''
    
    '''
    tifs = get_tif_names('Kort')
    print(check_md5s('Kort', tifs))
    check_md5('Kort','DSM_1km_6059_670')
    '''

    #options = {}
    #options['unzip'] = True
    #options['zip_delete'] = True
    #options['numpy'] = True
    ##options['tm_delete'] = True
    
    #pool_url = ThreadPool(12)
    #pool_url.map(sk.get_DTMfiles, znames)
    #pool_url.close()
    #pool_url.join()
    
    sk = sdfeKort()
    znames = sk.get_DSMzipnames()
    sk.get_DSMfiles(znames)
    
    
    #a = np.array([[1.15,2,3],[4,5,6]])
    #print(np.shape(a))
    #b = round_up(a,1)*10
    
    #print(b)
    
    #conv_tif_np16int('AllMaps','AllMapsNpz')
    
    
    
    load_folder = 'AllMaps_DSM'
    tnames = get_file_names(load_folder, ftype='tif')
    
    pool = Pool()
    
    data_out = pool.map(conv_tif_np16int, tnames)
    
    pool.close()
    pool.join()
    
    
    
