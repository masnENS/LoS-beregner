import multiprocessing as mp
import sys
import os.path
import pickle
from PyQt5.QtCore import pyqtSignal, QObject, QThread, Qt
from PyQt5.QtWidgets import QMainWindow, QGridLayout, QWidget, QLabel, QLineEdit, QPushButton, QTableWidget, QTableWidgetItem, QMessageBox, QFileDialog, QFormLayout, QApplication
from csv import reader, writer
from os.path import basename
import time
import LoSfuncs as lf
import LoSmpfuncs as lfmp
from multiprocessing import Pool
import datetime

#C:\Python34\Scripts\pyinstaller.exe H:\Fast_trdls\Python\07-07-2017\LoSCalc.py --onefile

def except_hook(cls, exception, traceback):
    sys.__excepthook__(cls, exception, traceback)

class finishedSignal(QObject): # Send signal to GUI when calculations have finished
    sig = pyqtSignal(str)

class PLSignal(QObject): # Send updated plot data to GUI during calculations
    sig = pyqtSignal(list,list)

class statusSignal(QObject): # Send signal to GUI to update status messages during calculations
    sig = pyqtSignal(str)

class calculate(QThread): # The actual LoS-calculations
    def __init__(self, basePosition, inputOK, parent=None):
            QThread.__init__(self, parent)
            self.FSig = finishedSignal()
            self.PLSig = PLSignal()
            self.statusSig = statusSignal()

            self.basePosition = basePosition
            self.inputOK = inputOK
            
    def setParams(self, basePosition, inputOK, name):
        self.basePosition = basePosition
        self.inputOK = inputOK
        self.name = name
        
    def run(self):
        if self.inputOK:  # Don't run if the input validation failed
            #######################################################################################################
            # Uncomment the next lines if the program should automatically download maps
            logon = ['',''] # ['user','password'], need user and password from kortforsyningen
            for BP in self.basePosition:
                lf.downloadMaps(self, logon, BP)
            #######################################################################################################
            
            self.statusSig.sig.emit('Beregner - vent venligst')
            start_time = time.time()
            multibase = lfmp.findCoveredMP()
            cpucount = mp.cpu_count()
            self.statusSig.sig.emit('Beregner - vent venligst - bruger '+str(cpucount)+' kerner')
            pool = Pool()
            
            data_out = pool.map(multibase.findBase, self.basePosition)
            pool.close()
            pool.join()
            t_out = time.time() - start_time
            
            #######################################################################################################
            # Organize covered addresses after multiprocessing
            coordinatesList = [] # List of client position coordinates
            addressList = [] # List of client position addresses
            LoSList = [] # LoS for the client position (0 or 1)
            bwList = [] # Bandwidth available to the address
            baseList = []
            
            for base in data_out:
                for i4, add in enumerate(base['address']):
                    addressList += [add]
                    coordinatesList += [base['coordinates'][i4]]
                    LoSList += [base['LoS'][i4]]
                    bwList += [base['BW'][i4]]
                    baseList += [base['name']]
            
            #######################################################################################################
            # Deduplicate covered addresses
            coveredCoordinates = []
            coveredUid = []
            coveredAddress = []
            coveredLoS = []
            coveredBase = []
            coveredBW = []

            for i in range(0,len(addressList)):
                if (addressList[i].get('uid') not in coveredUid):
                    coveredCoordinates.append(coordinatesList[i])
                    coveredUid.append(addressList[i].get('uid'))
                    coveredAddress.append(addressList[i])
                    coveredLoS.append(LoSList[i])
                    coveredBase.append([baseList[i]])
                    coveredBW.append(bwList[i])
                    
                    cvold = coveredBase[coveredUid.index(addressList[i].get('uid'))]
                    if baseList[i] not in cvold:
                        print('regiruesgbuire')
                        coveredBase[coveredUid.index(addressList[i].get('uid'))] += [baseList[i]]
                 
                elif (LoSList[i]):
                    #coveredUid[coveredUid.index(addressList[i].get('uid'))] = LoSList[i]
                    coveredLoS[coveredUid.index(addressList[i].get('uid'))] = LoSList[i]

                if (addressList[i].get('uid') in coveredUid):
                    if coveredBW[coveredUid.index(addressList[i].get('uid'))][0] < bwList[i][0]:
                        coveredBW[coveredUid.index(addressList[i].get('uid'))] = bwList[i]
                    elif coveredBW[coveredUid.index(addressList[i].get('uid'))][0] == bwList[i][0]:
                        if coveredBW[coveredUid.index(addressList[i].get('uid'))][1] < bwList[i][1]:
                            coveredBW[coveredUid.index(addressList[i].get('uid'))] = bwList[i]
                
                cvold = coveredBase[coveredUid.index(addressList[i].get('uid'))]  
                if baseList[i] not in cvold:
                    coveredBase[coveredUid.index(addressList[i].get('uid'))] += [baseList[i]]
                    
            print("Calculations took %s seconds ---" % (t_out))
            
            lf.saveData(self.name,coveredAddress, coveredLoS, coveredBW, coveredBase)
            lf.saveKML(self.name,coveredCoordinates, coveredAddress, coveredLoS)
            
            #######################################################################################################
            
        self.FSig.sig.emit('OK') # Calculations are done!

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self,parent)

        # Load settings from "indstillinger"-file, if it exists
        if os.path.isfile('indstillinger'):
            with open('indstillinger','rb') as f:
                self.name, self.basePosition, self.logon, self.mapRange = pickle.load(f)

        # Otherwise, use default settings     
        else:
            self.name = ''
            self.basePosition = []
            self.logon = ['ASN','ASN']
            self.mapRange = [0,0,1,0]

            with open('indstillinger', 'wb') as f:  # Python 3: open(..., 'wb') # Create a settings file
                pickle.dump([self.name, self.basePosition, self.logon, self.mapRange], f,0)

        # Create a simple GUI layout
        self.mainWin = QWidget()
        self.layout = QGridLayout()
        self.layout.setColumnStretch(0, 1)
        self.layout.setColumnStretch (1, 1)
        self.layout.setColumnMinimumWidth (0, 100)
        self.layout.setColumnMinimumWidth (1, 100)

        # Input field used for name of output files
        self.lNavn = QLabel()
        self.lNavn.setText("Virksomhedsnavn:")
        self.lNavn.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.layout.addWidget(self.lNavn,0,0)

        self.leNavn = QLineEdit()
        self.leNavn.setText(str(self.name))
        self.layout.addWidget(self.leNavn,0,1)

        # Load and save base stations from/to csv files
        self.loadBtn = QPushButton("Hent fra csv")
        self.loadBtn.clicked.connect(self.importBasePositions)
        self.layout.addWidget(self.loadBtn,1,0)

        self.saveBtn = QPushButton("Gem som csv")
        self.saveBtn.clicked.connect(self.saveSettings) 
        self.saveBtn.clicked.connect(self.exportBasePositions)
        self.layout.addWidget(self.saveBtn,1,1)

        # The input table contains all base stations to be calculated
        self.inputTable = QTableWidget()
        self.inputTable.setRowCount(1000)
        #self.keys = ['name','long','lat','radius','hbase', 'hclient','rhclient', 'download', 'upload', 'freq', 'thetamin', 'thetamax']
        self.keys = ['name','lat','long','radius','hbase', 'download', 'upload', 'freq', 'thetamin', 'thetamax', 'hclient','rhclient']
        #self.titles = ['Navn',  'Maks. radius [m]', 'Senderens \n højde [m]', 'Modtagerens \n højde (jord) [m]'
        #              , 'Modtagerens \n højde (tag) [m]' ,'Download \n [Mbit/s]','Upload \n [Mbit/s]','Frekvens \n [MHz]','Min. vinkel','Maks. vinkel']
        self.titles = ['Navn', 'Breddegrad\n [WGS84]', 'Længdegrad\n [WGS84]', 'Maks. radius [m]', 'Senderens \n højde [m]','Download \n [Mbit/s]','Upload \n [Mbit/s]','Frekvens \n [MHz]','Min. vinkel','Maks. vinkel', 'Modtagerens \n højde (jord) [m]', 'Modtagerens \n højde (tag) [m]'] 
        
        #self.defaultValues = ['','','','','','','','','','5800','0','360']
        self.defaultValues = ['','','','','','','','5800','0','360','10','2']# Usually, base stations operate at 5.8 GHz and 360 degrees
        
        self.baseLength = len(self.keys)
        self.nreq = 6
        self.inputTable.setColumnCount(self.baseLength)

        self.inputTable.setHorizontalHeaderLabels(self.titles)

        # Set the values in the input table, either from the settings file, or from default values
        for n in range(0,1000):
            for m in range(0,self.baseLength):
                try:
                    self.inputTable.setItem(n, m, QTableWidgetItem(str(self.basePosition[n].get(self.keys[m]))))
                
                except IndexError:
                    self.inputTable.setItem(n, m, QTableWidgetItem(self.defaultValues[m]))

        self.layout.addWidget(self.inputTable,2,0,1,0)

        self.inputOK = 1 # Initially, input is assumed to be OK...

        # Push to calculate!
        self.calcBtn = QPushButton("Beregn")
        self.calcThread = calculate(self.basePosition, self.inputOK)
        self.calcBtn.clicked.connect(self.saveSettings) # Save settings and input before calculating
        
        self.calcThread.FSig.sig.connect(self.calcThreadComplete) # Let the GUI know when calculations are done
        self.calcThread.PLSig.sig.connect(self.plotUpdate) # Send updated maps and graphs to the GUI
        self.calcThread.statusSig.sig.connect(self.statusUpdate) # Update the status message on the calc button
        
        self.calcBtn.clicked.connect(self.doCalculations) # Now we're ready to calculate!
        
        self.layout.addWidget(self.calcBtn,3,0,1,0)

        #self.image = QImage(480, 480, QImage.Format_ARGB32) # Displays the maps with sight lines
        #intial_color = qRgb(0, 0, 0)
        #self.image.fill(qRgb(0,0,0))
        #self.image_label = QLabel(" ")
        #self.image_label.setPixmap(QPixmap.fromImage(self.image))
        #self.layout.addWidget(self.image_label,4,0,1,1)

        #self.plot = pg.PlotWidget() # Displays surface heights and Fresnel zones
        #self.layout.addWidget(self.plot,4,0,1,1)
        
        self.layout.setRowMinimumHeight(2,300)
        self.layout.setRowMinimumHeight(4,100)
        self.layout.setRowStretch(2,100)
        
        self.centralwidget = QWidget(self)
        self.setCentralWidget(self.centralwidget)
        
        self.centralwidget.setLayout(self.layout)
        # End main window

    ###########################################################################################################################################################
    # Convert base station info from input table
    def convertInputTable(self):
        text = [[str(self.inputTable.item(m,n).text()) for n in range(0,self.baseLength)] for m in range(0,1000)]
       
        self.basePosition = [] 
        self.inputOK = 1
        for n in range(0,1000):
            inputBase = text[n]
            basePos = {}
            
            if any(inputBase[0:self.nreq])&(self.inputOK > 0):
                try:
                    basePos[self.keys[0]] = str(inputBase[0])
                    for i in range(1,self.baseLength):
                        basePos[self.keys[i]] = float(inputBase[i].replace(',','.'))

                    if (basePos['long'] > 8)&(basePos['long'] < 16)&(basePos['lat'] > 54)&(basePos['lat'] < 58):
                        self.basePosition.append(basePos)
                    
                    else:  
                        inputError = QMessageBox()
                        inputError.setText(str('Fejl i linje '+str(n+1)+', koordinater er ikke i Danmark'))
                        inputError.exec_();
                        self.inputOK = 0
                  
                except ValueError:
                    inputError = QMessageBox()
                    inputError.setText(str('Fejl i linje '+str(n+1)))
                    inputError.exec_();
                    self.inputOK = 0
        if len(self.basePosition)==0:
            self.inputOK = 0
            inputError = QMessageBox()
            inputError.setText('Fejl: Intet input')
            inputError.exec_();
            
    # End conversion
    ###########################################################################################################################################################

    ###########################################################################################################################################################
    # Base position import
    def importBasePositions(self):
        leFile = QLineEdit()
        leFile.setText(QFileDialog.getOpenFileName(self, 'Åben fil', '',"CSV (*.csv *.CSV)")[0])
       
        if leFile.text():
            self.leNavn.setText(os.path.splitext(basename(leFile.text()))[0])
            with open(leFile.text(),'rt', encoding='latin-1') as csvfile:
                inputReader = reader(csvfile, delimiter=';', quotechar = '|')
                n=-1
                for row in inputReader:
                    if n>=0:
                        for m in range(0,self.baseLength):
                            try:
                                self.inputTable.setItem(n, m, QTableWidgetItem(str(row[m])))
                            except IndexError:
                                self.inputTable.setItem(n, m, QTableWidgetItem(self.defaultValues[m]))
                    n = n+1
            for n2 in range(n,1000):
                for m in range(0,self.baseLength):
                    self.inputTable.setItem(n2, m, QTableWidgetItem(self.defaultValues[m]))
                    
    # End base position import
    ###########################################################################################################################################################

    ###########################################################################################################################################################
    # Base position export
    def exportBasePositions(self):
        leFile = QLineEdit()
        #leFile.setText(QFileDialog.getSaveFileName(self, 'Åben fil', '',"CSV (*.csv *.CSV)")[0])
        leFile.setText(self.name+'_'+datetime.datetime.now().strftime('%d-%b-%Y-%H%M')+'.csv') #'%-d-%b-%Y-%-H:%M'
        
        if leFile.text():
            with open(leFile.text(),'w', newline='') as csvfile:
                baseWriter = writer(csvfile, delimiter=';', quotechar = '|')

                baseWriter.writerow([x for x in ['Navn','Breddegrad','Længdegrad','Maks. Radius (m)','Senderens højde(m)','Download (Mbit/s)','Upload (Mbit/s)','Frekvens (Mhz)','Min. Vinkel','Maks. Vinkel']])
                for n in range(0,len(self.basePosition)):
                    BP = [self.basePosition[n].get(i) for i in self.keys]
                    if any(BP[0:]):
                        baseWriter.writerow([str(x).replace('.', ',') for x in BP])
                        
    # End base position export
    ###########################################################################################################################################################

    ###########################################################################################################################################################
    # Terms and conditions
    def terms(self):
        self.termsWin = QWidget()
        tlayout = QFormLayout()
        
        lTerms = QLabel()
        lTerms.setWordWrap(True)
        lTerms.setText("Ved at anvende softwaren accepterer brugeren denne aftale, jf. nedenfor. \n \
        Energistyrelsen har udviklet denne software med det formål at hjælpe udbydere af faste trådløse forbindelser med at vurdere deres dækning. \n \
        Energistyrelsen har tilstræbt, at softwaren - på grundlag af de oplysninger, brugeren selv indtaster - så præcist som muligt vurderer, hvilke adresser der vil kunne dækkes fra en given basestation. \n \
        Energistyrelsen fraskriver sig ethvert ansvar for tab eller skade, – såvel direkte som indirekte – der måtte opstå som følge af brugen af softwaren. \n \
        Energistyrelsen understreger, at softwarens teoretisk beregnede resultater er vejledende, \n \
        og at den reelle dækning af faste trådløse forbindelser blandt andet vil afhænge af nøjagtigheden af de benyttede terrænkort samt af de praktiske muligheder for placering af klientudstyr. \n \
        Energistyrelsen anbefaler derfor, at brugeren sørger for at verificere alle beregningsresultater.")
        tlayout.addRow(lTerms)
        
        okBtn = QPushButton("Accepter")
        tlayout.addRow(okBtn)
        okBtn.clicked.connect(self.termsWin.close)
        
        self.termsWin.setWindowTitle("Betingelser")
        self.termsWin.setWindowModality(Qt.ApplicationModal)
        self.termsWin.setLayout(tlayout)
        
        #sys.exit(self.terms.exec_())
    # End terms and conditions
    ###########################################################################################################################################################

    ###########################################################################################################################################################
    # Save settings
    def saveSettings(self):
        self.convertInputTable()
        self.name = self.leNavn.text()
       
        with open('indstillinger', 'wb') as f:  # Python 3: open(..., 'wb')
            pickle.dump([self.name, self.basePosition, self.logon, self.mapRange], f)
    # End save settings
    ###########################################################################################################################################################

    def doCalculations(self):
        if not self.calcThread.isRunning():
            self.calcThread.exiting=False
            self.calcThread.setParams(self.basePosition, self.inputOK, self.leNavn.text())
            self.calcThread.start()
            self.calcBtn.setText('Beregner - vent venligst')
            self.calcBtn.setEnabled(False)
            self.loadBtn.setEnabled(False)
            self.saveBtn.setEnabled(False)
            self.inputTable.setEnabled(False)
            self.leNavn.setEnabled(False)

    def calcThreadComplete(self):
        self.calcBtn.setText('Beregning færdig. Tryk for at beregne igen.')
        self.calcBtn.setEnabled(True)
        self.loadBtn.setEnabled(True)
        self.saveBtn.setEnabled(True)
        self.inputTable.setEnabled(True)
        self.leNavn.setEnabled(True)
        
    def plotUpdate(self,imgMat, heights):
        colours = ['w','w','g']

        self.plot.clear()
        for n in range(0,len(heights)):
            self.plot.plot(heights[n][0],heights[n][1],pen=colours[n])

        #pixmap = QPixmap()
        #pixmap.load(os.getcwd() + "\img.png")
        
        #if not pixmap.isNull():
                #print('Load nyt billede:')
                #self.image_label.setPixmap( (pixmap).scaledToWidth(480) )

    def statusUpdate(self, text):
        self.calcBtn.setText(text)
                
if __name__=='__main__':
    mp.freeze_support()
    app = QApplication(sys.argv)
    window = MainWindow()
    window.setWindowTitle("LoS-beregner version 1.0")
    
    window.resize(960, 800)
    window.show()

    sys.excepthook = except_hook
    app.exec_()
    
    
