# ============================================================================================= #
# ============================================================================================= #
"""
Class developped by Alan Loh to handle Fermi/LAT analyses.

Requirements:
--> work on sappcfermi
--> download the FT1 and FT2 files somewhere.
--> make sure that in your .bashrc (or else) you have the variables FERMI_DIR and LATEXTDIR properly defined, for e.g.: 
export FERMI_DIR=/dsm/saplxglast/glast/sas/Binary/ScienceTools-RELEASE-11-05-02
export LATEXTDIR=/dsm/saplxglast/home/aloh/ALAN/Extended_archive_v15/Templates
--> have your python path pointing to my script directory:
export PYTHONPATH=/dsm/saplxglast/home/aloh/ALAN/Scripts:$PYTHONPATH
"""

import os, sys, glob, shutil, subprocess, time, itertools
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
#import gt_apps as apps # Fermi 

# For plots
from matplotlib import rc
rc('font',**{'family':'serif', 'serif':['Computer Modern']})
rc('text', usetex=True)
import matplotlib as mpl
mpl.use('Agg') # get rid of the 'DISPLAY' error when running under a distant server
import matplotlib.pylab as plt
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams["text.latex.preamble"].append(r'\usepackage[dvips]{graphicx}\usepackage{amssymb}\usepackage{amsmath}')
import matplotlib.ticker as mtick


# ============================================================================================= #
# ============================================================================================= #
# ============================================================================================= #
class algamma(object):
    def __init__(self):
        # Observation files
        self.ft1   = None
        self.ft2   = None
        self.ephem = None
        self.image = None

        # Observation properties
        self.tstart   = None
        self.tstop    = None
        self.metstart = None
        self.metstop  = None
        self.emin     = None
        self.emax     = None
        self.nevents  = None
        self.passver  = None 
        self.evclass  = None
        self.evtype   = None
        self.ra       = None
        self.dec      = None
        self.rad      = None
        self.zmax     = 90
        self.irf      = None
        self.frac     = None
        self.tsmin    = 25

        # Source models / analysis
        self.diffModel = None
        self.model     = None
        self.lcmodel   = None  
        self.srcmodel  = 'powerlaw'
        self.mode      = 'binned' # /or unbinned
        self.fermicat  = None

        # Paths
        self.datapath = None
        self.workpath = None
        self.suffix   = ''
        self.diffpath = os.path.join(os.environ['FERMI_DIR'], 'diffuseModels/v5r0')
        self.extpath  = os.environ['LATEXTDIR']

        # Plots
        self.csys  = 'CEL' # coordinate system 'CEL' / 'GAL'
        self.imwid = 100 # pixels
        self.binsz = 0.2 # deg / pix
        self.color = 'black'
        self.lblue = '#79a6d2'
        self.loran = '#ffcc66'

        # General
        self.clobber = False

    # ============================================================================================= #
    #                                     Attributes Checking                                       #
    # ============================================================================================= #
    @property
    def ft1(self):
        return self._ft1
    @ft1.setter
    def ft1(self, f):
        if f is None:
            self._ft1 = None
        else:
            if not os.path.isfile(f):
                print("\t=== FT1 file '{}' not found. ===".format(f))
                self._ft1 = None
            else:
                self._ft1 = os.path.abspath(f)
                self._readFT1()
                self._initNames()

    @property
    def ft2(self):
        return self._ft2
    @ft2.setter
    def ft2(self, f):
        if f is None:
            self._ft2 = None
        else:
            if not os.path.isfile(f):
                print("\t=== FT2 file '{}' not found. ===".format(f))
                self._ft2 = None
            else:
                self._ft2 = os.path.abspath(f)

    @property
    def datapath(self):
        return self._datapath
    @datapath.setter
    def datapath(self, d):
        if d is None:
            self._datapath = None
        else:
            self._datapath = os.path.abspath(d)

    @property
    def workpath(self):
        return self._workpath
    @workpath.setter
    def workpath(self, w):
        if w is None:
            self._workpath = None
        else:
            self._workpath = os.path.abspath(w)
            self._initNames()
            if not os.path.isdir(self._workpath): 
                os.makedirs(self._workpath)

    @property
    def model(self):
        return self._model
    @model.setter
    def model(self, m):
        if m is None:
            self._model = None
        else:
            if os.path.isfile(m):
                # the XML model has already been created, no need for a second one
                self._model = m
            else:
                # create the model using make3FGLxml
                self._model = m
                if 
                self.makeModel()
                
    @property
    def suffix(self):
        return self._suffix
    @suffix.setter
    def suffix(self, s):
        if s is '':
            self._suffix = ''
        elif isinstance(s, str):
            self._suffix = s
            self._initNames()
        else:
            self._suffix = ''

    # ============================================================================================= #
    #                                           Analysis                                            #
    # ============================================================================================= #
    def makeModel(self):
        """ Make a XML model based on info inside the FT1 file
        """

        # Get the script
        modelScript = os.path.join(self.datapath, 'make3FGLxml.py')
        if not os.path.isfile(modelScript):
            # download it
            print("\t=== Downloading make3FGLxml.py ===")
            os.system('wget https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/make3FGLxml.py -O {}'.format(modelScript))

        # Create the model using Tyrel's script
        galModel = os.path.join(self.diffpath, 'gll_iem_v06.fits')
        isoModel = os.path.join(self.diffpath, 'iso_'+self.irf+'_v06.txt')
        if (not os.path.isfile(galModel)) or (not os.path.isfile(isoModel)):
            print("\t=== Unable to find the diffuse models, check the variable '$FERMI_DIR' ===")
            return
        if not os.path.isdir(self.extpath):
            print("\t=== Unable to find models of extended sources, check the variable '$LATEXTDIR' ===")
            return
        if not os.path.isfile(self.fermicat):
            # download it
            print("\t=== Downloading 3FGL catalog ===")
            os.system('wget https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4yr_catalog/gll_psc_v16.fit -O {}'.format(self.fermicat))

        os.popen("python {} {} {} -o {} -G {} -g 'gll_iem_v06'\
            -I {} -i 'iso_source_v06' -e {} -r 5 -R 10 -ER 10\
            -s 9 -m False -GIF False".format(modelScript, self.fermicat,
            self.ft1, self.model, galModel, isoModel, self.extpath))

        # Add the target to the model
        tmpName = self.model + '.tmp'
        rfil = open(self.model, 'r')
        wfil = open(tmpName, 'w')
        # Copy the XML to the temporary model
        wfil.writelines([l for l in rfil.readlines() if not l=='</source_library>']) # copy everything but the last line
        wfil.write(' <source ROI_Center_Distance="0.00" name="TARGET" type="PointSource">\n')
        wfil.write('  <spectrum type="PowerLaw2">\n')
        wfil.write('   <parameter free="1" max="1000" min="1e-05" name="Integral" scale="1e-08" value="0.3591824258"/>\n')
        wfil.write('   <parameter free="1" max="1" min="-5" name="Index" scale="1" value="-2.7"/>\n')
        wfil.write('   <parameter free="0" max="1000000" min="20" name="LowerLimit" scale="1" value="100"/>\n')
        wfil.write('<parameter free="0" max="1000000" min="20" name="UpperLimit" scale="1" value="100000"/>\n')
        wfil.write('  </spectrum>\n')
        wfil.write('  <spatialModel type="SkyDirFunction">\n')
        wfil.write('   <parameter free="0" max="360.0" min="-360.0" name="RA" scale="1.0" value="'+str(self.ra)+'"/>\n')
        wfil.write('   <parameter free="0" max="360.0" min="-360.0" name="DEC" scale="1.0" value="'+str(self.dec)+'"/>\n')
        wfil.write('  </spatialModel>\n')
        wfil.write(' </source>\n')
        wfil.write('</source_library>\n')
        rfil.close()
        wfil.close()

        os.remove(self.model)
        os.rename(tmpName, self.model)
        
        print("\t=== Source model {} added ===".format(self.model))
        return
    def getSrc(self):
        """ Get the source list contained in an xml file
        """
        xml    = open(self.model, 'r')
        keywd1 = ['RA', 'DEC', 'PointSource']
        ra  = []
        dec = []
        nam = []
        sep = []
        target = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree, frame='icrs')        
        for line in xml :
            if keywd1[0] in line:
                ra.append( float(line.split('"')[-2]) )
            if keywd1[1] in line:
                dec.append( float(line.split('"')[-2]) )
                s = SkyCoord(ra=ra[-1]*u.degree, dec=dec[-1]*u.degree, frame='icrs')
                sep.append(target.separation(s).deg)
            if keywd1[2] in line:
                nam.append( line.split('"')[3].split()[-1] ) # no '3FGL'
        xml.close()

        if self.csys == 'GAL':
            srcPos  = SkyCoord(np.array(ra)*u.degree, np.array(dec)*u.degree, frame='icrs')
            ra, dec = srcPos.galactic.l.deg, srcPos.galactic.b.deg

        srcs = Table([ra, dec, nam, sep], names=('RA', 'DEC', 'Name', 'Separation'))
        return srcs
    def getResult(self, param='Flux', src='TARGET'):
        """ Get the result from a likelihood fit

        Parameters
        ----------
        param: string
            Parameter to look for ('Flux', 'TS value', 'Index', 'upplim') 
        """
        lines = open(self.outgtlike, 'r').readlines()

        if param == 'upplim':
            value = -1
            error = 0 # default value
            for line in lines:
                if 'Upper limit' in line: 
                    value = line.split()[-2]
        else:
            CodeString = ''
            for line in lines:
                if not 'Upper limit' in line:
                    CodeString += line[:-1]

            MyData = eval(CodeString)
            Values = MyData[src][param].split()
            value  = float(Values[0])
            try: 
                error = float(Values[2])
            except:
                error = 0
        
        return value, error
    def rmGt(self):
        """ Remove all subsidiary files
        """
        gtfiles = [self.outselect, self.outmktime, self.outltcube,
            self.outbincub, self.outbinmap, self.outbinexp, 
            self.outexpmap, self.outsrcmap, 
            os.path.join(self.workpath, 'SrcList_cntspec'+self.suffix+'.fits'),
            os.path.join(self.workpath, 'SrcList_cntspec'+self.suffix+'.log')]
        for f in gtfiles:
            if os.path.isfile(f):
                os.remove(f)
        return
    def pulsEphem(self):
        """ Add a column PULSE_PHASE to a FT1 file
        """

        hduMain = fits.open(self.ft1)

        # --------------------------------------------------------------------------------------------- #
        # Split the FT1 file every 4000 events
        noEv   = 0
        deltEv = 5000
        count  = 0
        wfil   = open(os.path.dirname(self.ft1) + os.path.basename('tmpFT1.lis'), 'w')
        while noEv <= self.nevents:
            hduCols = []
            for colname, form, uni in zip(hduMain['EVENTS'].columns.names, hduMain['EVENTS'].columns.formats, hduMain['EVENTS'].columns.units):
                hduCols.append( fits.Column(name=colname, array=hduMain['EVENTS'].data[colname][noEv:noEv+deltEv], format=form, unit=uni) )
            # Updte the tstart and tstop in the header in order for tempo2 to work...
            hduMain['EVENTS'].header.set('TSTART', hduMain['EVENTS'].data['TIME'][noEv:noEv+deltEv][0])
            hduMain['EVENTS'].header.set('TSTOP', hduMain['EVENTS'].data['TIME'][noEv:noEv+deltEv][-1])
            newHDU  = fits.BinTableHDU.from_columns(hduCols, name='EVENTS', header=hduMain['EVENTS'].header) 
            hdulist = fits.HDUList([hduMain['PRIMARY'], newHDU, hduMain['GTI']])
            tmpName = os.path.dirname(self.ft1)+os.path.basename('tempFT1_'+str(count)+'.fits')
            hdulist.writeto(tmpName, clobber=True)
            wfil.write(tmpName + '\n')
            noEv  += deltEv
            count += 1
        if noEv != self.nevents:
            hduCols = []
            noEv -= deltEv
            for colname, form, uni in zip(hduMain['EVENTS'].columns.names, hduMain['EVENTS'].columns.formats, hduMain['EVENTS'].columns.units):
                hduCols.append( fits.Column(name=colname, array=hduMain['EVENTS'].data[colname][noEv:self.nevents], format=form, unit=uni) )
            hduMain['EVENTS'].header.set('TSTART', hduMain['EVENTS'].data['TIME'][noEv:self.nevents][0])
            hduMain['EVENTS'].header.set('TSTOP', hduMain['EVENTS'].data['TIME'][noEv:self.nevents][-1])
            newHDU  = fits.BinTableHDU.from_columns(hduCols, name='EVENTS', header=hduMain['EVENTS'].header)
            hdulist = fits.HDUList([hduMain['PRIMARY'], newHDU, hduMain['GTI']])
            tmpName = os.path.dirname(self.ft1)+os.path.basename('tempFT1_'+str(count)+'.fits')
            hdulist.writeto(tmpName, clobber=True)
            wfil.write(tmpName + '\n')
        wfil.close()

        hduMain.close()

        # --------------------------------------------------------------------------------------------- #
        # Run tempo2 for each piece of the FT1
        rfil = open(os.path.dirname(self.ft1) + 'tmpFT1.lis', 'r')
        percent = 0
        nbFiles = sum(1 for line in open(os.path.dirname(self.ft1) + 'tmpFT1.lis', 'r'))
        count   = 0
        for tmpFil in rfil:
            # Print a progression bar every 5%
            if ( count / np.floor(nbFiles) * 100 ) >= percent:
                self._progressBar(percent, printEvery=5)
                percent += 5
            with open(os.devnull, 'wb') as devnull:
                subprocess.check_call(['/dsm/fermi/fermifast/glast/tempo2-2013.9.1/tempo2',
                    '-gr', 'fermi', '-ft1', tmpFil[:-1], '-ft2', self.ft2, '-f', self.ephem,
                    '-phase'], stdout=devnull, stderr=subprocess.STDOUT)
            count += 1
        # Replace the old ft1 by the new one with the PULSE_PHASE column
        #os.remove()
        self._gtSelect(data = os.path.dirname(self.ft1) + os.path.basename('tmpFT1.lis'))




        #self.nevents
        #J2032+4127_54683_57791_chol_pos.par
        #os.popen("tempo2 -gr fermi -ft1 {} -ft2 {} -f {} -phase".format(self.ft1, self.ft2, self.ephem))
    def pulsCut(self, phcut):
        """ Copy a FT1 whithout the phases of a nearby pulsar
        """
        if not os.path.isfile(self.outmktime):
            ("\t=== Make sure that self.workpath points towards a valid directory ===")
            return

        if 'PULSE_PHASE' not in fits.getdata(self.outmktime).columns.names:
            print("\t=== FT1 file does not have a 'PULSE_PHASE' column, run tempo2 first ===")
            return

        frac   = np.sum(phcut[1::2]) - np.sum(phcut[0::2])
        outfil = self.outmktime[:-5]+'_NoPulse_'+str(frac)+'.fits'
        cutCmd = "(PULSE_PHASE >= {} && PULSE_PHASE <= {})".format(phcut[0], phcut[1])
        if len(phcut) > 2:
            for i in np.arange(2, len(phcut), 2):
                cutCmd += " || (PULSE_PHASE >= {} && PULSE_PHASE <= {})".format(phcut[i], phcut[i+1])
        
        os.popen('fcopy "{}[{}]" {}'.format(self.outmktime, cutCmd, outfil))

        print("\t=== File '{}' created ===".format(outfil))
        return
    def diffRsp(self, n=2):
        """ Add the diffuse component response to a FT1 file
        """ 
        #  /dsm/saplxglast/home/aloh/ALAN/CygX-3_FermiAnalysis/MODEL/Model_GTDIFFRSP.xml 

        # --------------------------------------------------------------------------------------------- #
        # Create a directory to store the little FT1 pieces
        self.workpath = os.path.join(self.datapath, 'GTDIFFRSP', 'Runs')
        if self.diffModel is None:
            print('\t=== A diffuse model file needs to be provided --> self.diffModel ===')
            return
        tmpMETstart = self.metstart
        tmpMETstop  = self.metstop

        # --------------------------------------------------------------------------------------------- #
        # Split the calculation
        tWidth = 2e6
        self.metstart = tmpMETstart
        self.metstop  = tmpMETstart + tWidth

        count = 0
        while self.metstop < tmpMETstop + tWidth:
            filesRunning = glob.glob( os.path.join(self.workpath, 'tmp_*.py') )
            while len(filesRunning) == n: 
                # Limit the number of parallel runs
                time.sleep(60)
                filesRunning = glob.glob( os.path.join(self.workpath, 'tmp_*.py') )

            fil = os.path.join( self.workpath, 'tmp_'+str(count)+'.py' )
            tmp = open(fil, 'w')
            tmp.write("import algamma; import os; a=algamma.algamma(); a.ft1='{}';\
                a.ft2='{}'; a.diffModel='{}'; a.metstart={}; a.metstop={}; a.suffix='_Chk_{}';\
                a.workpath='{}'; a._initNames(); a._gtSelect(); a._gtMktime();\
                a._gtDiffrsp(); os.remove(a.outselect); os.remove('{}')"
                .format(os.path.join(self.datapath, self.ft1),
                os.path.join(self.datapath, self.ft2), self.diffModel,
                self.metstart, self.metstop, count, self.workpath, fil))
            # Launch the file
            os.popen("nohup python {} &".format(fil))
            tmp.close()
            
            count += 1
            self.metstart += tWidth
            self.metstop  += tWidth

        self._mergeDiffrsp()
        return
    def binAnalysis(self):
        """ Pipeline to do a binned Fermi analysis
        """
        self.mode = 'binned'
        # --------------------------------------------------------------------------------------------- #
        # Make sure that another working directory is selected
        if self.workpath == self.datapath:
            print("\t=== Variable 'self.workpath' is equal to 'self.datapath', provide another ===")
            return
        else:
            if os.path.isfile(self.outgtlike):
                print("\t=== Directory {} already contains a complete analysis, remove the .dat file ===".format(self.workpath))
                return
            else:
                pass
            print("\t=== Binned analysis will be computed in '{}' ===".format(self.workpath))

        # --------------------------------------------------------------------------------------------- #
        # Create a temporary python script and launch the Science Tools
        fil = os.path.join(self.workpath, 'tmp_BinnedAnalysis'+self.suffix+'.py')
        tmp = open(fil, 'w')
        tmp.write("import algamma; import os; a=algamma.algamma(); a.ft1='{}';\
            a.ft2='{}'; a.metstart={}; a.metstop={}; a.emin={}; a.emax={}; a.suffix='{}';\
            a.workpath='{}'; a._gtSelect(); a._gtMktime();\
            a._gtLtcube(); a._gtBincube(); a._gtExpmap(); a._gtSrcmap();\
            a._gtLike(); os.remove('{}')".format(self.ft1, self.ft2, 
            self.metstart, self.metstop, self.emin, self.emax,
            self.suffix, self.workpath, fil))
        # Launch the file
        os.popen("nohup python {} &".format(fil))
        tmp.close()

        return
    def unbinAnalysis(self):
        """ Pipeline to do an unbinned Fermi analysis
        """
        self.mode = 'unbinned'
        # --------------------------------------------------------------------------------------------- #
        # Make sure that another working directory is selected
        if self.workpath == self.datapath:
            print("\t=== Variable 'self.workpath' is equal to 'self.datapath', provide another ===")
            return
        else:
            if os.path.isfile(self.outgtlike):
                print("\t=== Directory {} already contains a complete analysis, remove the .dat file ===".format(self.workpath))
                return
            else:
                pass
            print("\t=== Unbinned analysis will be computed in '{}' ===".format(self.workpath))

        # --------------------------------------------------------------------------------------------- #
        # Create a temporary python script and launch the Science Tools
        fil = os.path.join(self.workpath, 'tmp_UnbinnedAnalysis'+self.suffix+'.py')
        tmp = open(fil, 'w')
        tmp.write("import algamma; import os; a=algamma.algamma(); a.ft1='{}';\
            a.ft2='{}'; a.metstart={}; a.metstop={}; a.emin={}; a.emax={}; a.suffix='{}';\
            a.workpath='{}'; a._gtSelect(); a._gtMktime(); a._gtLtcube(); a._gtExpmap();\
            a._gtLike(); os.remove('{}')".format(self.ft1, self.ft2, 
            self.metstart, self.metstop, self.emin, self.emax,
            self.suffix, self.workpath, fil))
        # Launch the file
        os.popen("nohup python {} &".format(fil))
        tmp.close()

        return
    def compareModel(self):
        """ Outline the major differences from the 3FGL model and the fitted sources
        """

        # --------------------------------------------------------------------------------------------- #
        # Store the Model parameters
        lines   = open(self.model, 'r').readlines()
        MyModel = {}
        for line in lines:
            if ('<source' in line) & ('name=' in line):
                srcNam = line.split('"')[3]
                MyModel[ srcNam ] = {}
            elif ('<parameter' in line) & ('free="1"' in line):
                parNam = line.split('"')[7]
                parVal = float(line.split('"')[11])
                MyModel[ srcNam ][ parNam ] = parVal
            else:
                pass

        # --------------------------------------------------------------------------------------------- #
        # Store the fitted results
        lines = open(self.outgtlike, 'r').readlines()
        CodeString = ''
        for line in lines:
            if not 'Upper limit' in line:
                CodeString += line[:-1]
        MyData = eval(CodeString) # create a dictionnary

        # --------------------------------------------------------------------------------------------- #
        # Compare
        for key in MyData.keys():
            if 'TS value' in MyData[key].keys():
                # The source has been fitted
                print("--- {} ---".format(key))
                for k in MyModel[key].keys():
                    difference = 100* (MyModel[key][k] - float(MyData[key][k].split()[0])) / MyModel[key][k] 
                    print("{0:s} differs by {1:.2f} per cent".format(k, difference))
        return
    def doLCURVE(self, dt, shift=None, n=5):
        """ Light curve computation

        Parameters
        ----------
        dt: float
            Integration time (in hours)
        shift: float
            Starting time shift between two consecutive bins
        n: integer
            Number of parallel analyses (warning...)
        
        Example
        ----------
            import algamma
            a = algamma.algamma()
            a.ft1 = 'FT1_filtered_Chk_0.fits'
            a.ft2 = 'CygX3_FT2.fits'
            a.workpath = '../Lightcurves/LCTEST'
            a.outmodel = '../Binned/FirstAttempt/OutModel.xml'
            a.lightCurve(dt=24, n=2)
        """

        # --------------------------------------------------------------------------------------------- #
        # Make sure that another working directory is selected
        if self.workpath == self.datapath:
            print("\t=== Variable 'self.workpath' is equal to 'self.datapath', provide another ===")
            return
        # Check that FT1 has diffuse columns
        header = fits.getheader(self.ft1, ext=1)
        if header['DIFRSP0'] == 'NONE':
            print("\t=== Selected FT1 has no diffuse columns, run self.diffRsp() first ===")
            return
        if shift is None:
            shift = dt
        if not os.path.isfile(self.outmodel):
            print("\t=== Unable to find source model '{}' --> check self.outmodel ===".format(self.outmodel))
            return
        if (self.lcmodel is None) or (not os.path.isfile(self.lcmodel)):
            print("\t=== Creating a source model for the light curve ===")
            self.lcmodel = os.path.join(self.workpath, 'LCModel.xml') 
            rfil  = open(self.outmodel, 'r')
            wfil  = open(self.lcmodel, 'w')
            isSrc = False
            for line in rfil:
                if 'TARGET' in line: isSrc = True
                if 'LowerLimit'  in line: isSrc = False # let norm + index free
                if isSrc:
                    wfil.write(line)
                else:
                    # fix every parameter
                    wfil.write(line.replace('free="1"', 'free="0"'))
            rfil.close()
            wfil.close()
        # Convert dt and shift in seconds
        dt    *= 3600
        shift *= 3600

        # --------------------------------------------------------------------------------------------- #
        # Initilialize the result storage FITS file
        fitsNnam = os.path.join(self.workpath, 'LCresults.fits')
        if not os.path.isfile(fitsNnam):
            prihdr  = fits.Header() 
            prihdr.set('USER-CRE', 'A.Loh', 'File generator')
            prihdr.set('FT1', self.ft1, 'Data file')
            prihdr.set('FT2', self.ft2, 'Spacecraft file')
            prihdr.set('START', self.metstart, 'Start time (MET)')
            prihdr.set('STOP', self.metstop, 'Stop time (MET)')
            prihdr.set('RA', self.ra, 'RA (deg)')
            prihdr.set('DEC', self.dec, 'Dec (deg)')
            prihdr.set('DT', dt, 'Dec (deg)')
            prihdr.set('SHIFT', shift, 'Dec (deg)')
            hdu    = fits.PrimaryHDU(header=prihdr) # hdu table that will be filled
            hdu.writeto(fitsNnam)
            ct       = []
            tstart   = []
            tstop    = []
            mjd      = []
            mjderr   = []
            flux     = []
            fluxerr  = []
            index    = []
            indexerr = []
            ts       = []
            upperlim = []
            status   = []
            hdus = fits.open(fitsNnam)
        else:
            hdus     = fits.open(fitsNnam)
            ct       = hdus['LIGHTCURVE'].data['count']
            tstart   = hdus['LIGHTCURVE'].data['tstart']
            tstop    = hdus['LIGHTCURVE'].data['tstop']
            mjd      = hdus['LIGHTCURVE'].data['mjd']
            mjderr   = hdus['LIGHTCURVE'].data['mjderr']
            flux     = hdus['LIGHTCURVE'].data['flux']
            fluxerr  = hdus['LIGHTCURVE'].data['fluxerr']
            index    = hdus['LIGHTCURVE'].data['index']
            indexerr = hdus['LIGHTCURVE'].data['indexerr']
            ts       = hdus['LIGHTCURVE'].data['ts']
            upperlim = hdus['LIGHTCURVE'].data['upperlim']
            status   = hdus['LIGHTCURVE'].data['status']
            hdus.remove( hdus['LIGHTCURVE'] ) # it will be updated
            pass

        # --------------------------------------------------------------------------------------------- #
        # Split the calculation and create all the required files
        if len(ct) == 0:
            # Start from scratch
            tmpMETstart = self.metstart
            tmpMETstop  = self.metstop
            count       = 0
        else:
            # Start from the end of previous calculation
            tmpMETstart = tstop[-1]
            tmpMETstop  = self.metstop
            count       = ct[-1] + 1
        self.metstart = tmpMETstart
        self.metstop  = tmpMETstart + dt
        while self.metstop < tmpMETstop + dt:
            ct.append(count)
            tstart.append(self.metstart)
            tstop.append(self.metstop)
            mjd.append( self._MET2MJD( (self.metstop+self.metstart)/2. ) )
            mjderr.append( dt/86400./2. ) # in days
            flux.append(-1)
            fluxerr.append(-1)
            index.append(-1)
            indexerr.append(-1)
            ts.append(-1)
            upperlim.append(-1)
            status.append('pending')
            fil = os.path.join( self.workpath, 'tmp_'+str(count)+'.py' )
            if not os.path.isfile(fil):
                tmp = open(fil, 'w')
                tmp.write("import algamma; import os; from astropy.io import fits; a=algamma.algamma();\
                    a.ft1='{}'; a.ft2='{}'; a.mode='unbinned'; a.model='{}'; a.metstart={};\
                    a.metstop={}; a.suffix='_{}'; a.workpath='{}'; a._gtSelect(); a._gtMktime();\
                    a._gtLtcube(); a._gtExpmap(); a._gtLike(); a._uppLim()\nif os.path.isfile(a.outgtlike):\n\
                    \ta.rmGt(); a._lcResults(); os.remove('{}')"
                    .format(self.ft1, self.ft2, self.lcmodel, self.metstart, self.metstop,
                    count, self.workpath, fil))
                tmp.close()
            count += 1
            self.metstart += shift
            self.metstop  += shift
        # Create the fits file
        c1  = fits.Column(name='count', array=ct, format='I')
        c2  = fits.Column(name='tstart', array=tstart, format='D')
        c3  = fits.Column(name='tstop', array=tstop, format='D')
        c4  = fits.Column(name='mjd', array=mjd, format='D')
        c5  = fits.Column(name='mjderr', array=mjderr, format='D')
        c6  = fits.Column(name='flux', array=flux, format='D')
        c7  = fits.Column(name='fluxerr', array=fluxerr, format='D')
        c8  = fits.Column(name='index', array=flux, format='D')
        c9  = fits.Column(name='indexerr', array=fluxerr, format='D')
        c10 = fits.Column(name='ts', array=ts, format='D')
        c11 = fits.Column(name='upperlim', array=upperlim, format='D')
        c12 = fits.Column(name='status', array=status, format='10A')
        hdus.append(fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12], name='LIGHTCURVE'))
        hdus.writeto(fitsNnam, clobber=True)
        hdus.close()

        # --------------------------------------------------------------------------------------------- #
        # Run the analysis
        while fits.getdata(fitsNnam, ext=1)['status'].count('pending').sum() > 0:
            nbRun = fits.getdata(fitsNnam, ext=1)['status'].count('running').sum()
            while nbRun == n: 
                # Limit the number of parallel runs
                time.sleep(10)
                try:
                    nbRun = fits.getdata(fitsNnam, ext=1)['status'].count('running').sum()
                except:
                    time.sleep(30)
                    nbRun = fits.getdata(fitsNnam, ext=1)['status'].count('running').sum()
            ctRun = fits.getdata(fitsNnam, ext=1)['count'][fits.getdata(fitsNnam, ext=1)['status'] == 'pending']

            # Launch the file
            hdu = fits.open(fitsNnam)
            hdu['LIGHTCURVE'].data['status'][ctRun[0]] = 'running'
            hdu.writeto(fitsNnam, clobber=True)
            hdu.close()
            fil = os.path.join( self.workpath, 'tmp_'+str(ctRun[0])+'.py' )
            os.popen("nohup python {} &".format(fil))

        self.lCurve() # plot
        return
    def doTSMAP(self):
        """ Compute and plot a residual TS map
        """
        # --------------------------------------------------------------------------------------------- #
        # Make sure that another working directory is selected
        if self.workpath == self.datapath:
            print("\t=== Variable 'self.workpath' is equal to 'self.datapath', provide another ===")
            return
        else:
            if os.path.isfile(self.outtsmap):
                print("\t=== Directory {} already contains a TSmap ===".format(self.workpath))
                self.tsMap()
                return
            else:
                pass
            print("\t=== TS map will be computed in '{}' ===".format(self.workpath))

        # --------------------------------------------------------------------------------------------- #
        # Create a temporary python script and launch the Science Tools
        fil = os.path.join(self.workpath, 'tmp_TSmap.py')
        tmp = open(fil, 'w')
        tmp.write("import algamma; import os; a=algamma.algamma(); a.ft1='{}'; a.ft2='{}';\
            a.suffix='{}'; a.mode='{}'; a.csys='{}'; a.workpath='{}'; a._gtTSmap(); os.remove('{}')".format(self.ft1,
            self.ft2, self.suffix, self.mode, self.csys, self.workpath, fil))
        # Launch the file
        os.popen("nohup python {} &".format(fil))
        tmp.close()

        return
    def doSPEC(self, nbBins=4, n=2):
        """ Launch a series of likelihood analyses in order to derive a spectrum

        Parameters
        ----------
        nbBins: int
            Number of energy bins between self.emin and self.emax
        n: int
            Number of parrallel computations
        """    
        # --------------------------------------------------------------------------------------------- #
        # Determine the energy edges for the calculation 
        eEdges = np.logspace( np.log10(self.emin), np.log10(self.emax), nbBins+1)
        tmpEmin, tmpEmax = self.emin, self.emax

        # --------------------------------------------------------------------------------------------- #
        # Initilialize the result storage FITS file
        fitsNnam = os.path.join(self.workpath, 'SPECresults.fits')
        prihdr   = fits.Header() 
        prihdr.set('USER-CRE', 'A.Loh', 'File generator')
        prihdr.set('FT1', self.ft1, 'Data file')
        prihdr.set('FT2', self.ft2, 'Spacecraft file')
        prihdr.set('START', self.metstart, 'Start time (MET)')
        prihdr.set('STOP', self.metstop, 'Stop time (MET)')
        prihdr.set('EMIN', tmpEmin, 'Min Energy (MeV)')
        prihdr.set('EMAX', tmpEmax, 'Max energy (MeV)')
        prihdr.set('RA', self.ra, 'RA (deg)')
        prihdr.set('DEC', self.dec, 'Dec (deg)')
        prihdr.set('DT', dt, 'Dec (deg)')
        prihdr.set('SHIFT', shift, 'Dec (deg)')
        hdu      = fits.PrimaryHDU(header=prihdr) # hdu table that will be filled
        hdu.writeto(fitsNnam)
        ct       = []
        emin     = []
        emax     = []
        flux     = []
        fluxerr  = []
        index    = []
        indexerr = []
        ts       = []
        upperlim = []
        status   = []
        hdus = fits.open(fitsNnam)

        # --------------------------------------------------------------------------------------------- #
        # Loop through the energy bins and launch a likelihood fit
        # Split the calculation and create all the required files
        for count in xrange(nbBins):
            self.emin   = eEdges[count]
            self.emax   = eEdges[count + 1]
            self.suffix = '_'+str(count)
            ct.append(count)
            emin.append(self.emin)
            emax.append(self.emax)
            flux.append(-1)
            fluxerr.append(-1)
            index.append(-1)
            indexerr.append(-1)
            ts.append(-1)
            upperlim.append(-1)
            status.append('pending')

            fil = os.path.join( self.workpath, 'tmp_'+str(count)+'.py' )
            tmp = open(fil, 'w')
            if self.mode == 'binned':
                tmp.write("import algamma; import os; a=algamma.algamma(); a.ft1='{}';\
                    a.ft2='{}'; a.emin={}; a.emax={}; a.suffix='_{}';\
                    a.workpath='{}'; a._gtSelect(); a._gtMktime();\
                    a._gtLtcube(); a._gtBincube(); a._gtExpmap(); a._gtSrcmap();\
                    a._gtLike(); a._uppLim()\nif os.path.isfile(a.outgtlike):\n\
                    \ta.rmGt(); a._specResults(); os.remove('{}')".format(self.ft1, self.ft2, 
                    self.emin, self.emax, count, self.workpath, fil))
            elif self.mode == 'unbinned':
                tmp.write("import algamma; import os; from astropy.io import fits; a=algamma.algamma();\
                    a.ft1='{}'; a.ft2='{}'; a.mode='unbinned'; a.model='{}'; a.emin={};\
                    a.emax={}; a.suffix='_{}'; a.workpath='{}'; a._gtSelect(); a._gtMktime();\
                    a._gtLtcube(); a._gtExpmap(); a._gtLike(); a._uppLim()\nif os.path.isfile(a.outgtlike):\n\
                    \ta.rmGt(); a._specResults(); os.remove('{}')"
                    .format(self.ft1, self.ft2, self.lcmodel, self.emin, self.emax,
                    count, self.workpath, fil))
            else:
                pass
            tmp.close()

        # --------------------------------------------------------------------------------------------- #
        # Create the fits file
        c1  = fits.Column(name='count', array=ct, format='I')
        c2  = fits.Column(name='emin', array=tstart, format='D')
        c3  = fits.Column(name='emax', array=emax, format='D')
        c4  = fits.Column(name='flux', array=flux, format='D')
        c5  = fits.Column(name='fluxerr', array=fluxerr, format='D')
        c6  = fits.Column(name='index', array=flux, format='D')
        c7  = fits.Column(name='indexerr', array=fluxerr, format='D')
        c8 = fits.Column(name='ts', array=ts, format='D')
        c9 = fits.Column(name='upperlim', array=upperlim, format='D')
        c10 = fits.Column(name='status', array=status, format='10A')
        hdus.append(fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10], name='SPECTRUM'))
        hdus.writeto(fitsNnam, clobber=True)
        hdus.close()

        # --------------------------------------------------------------------------------------------- #
        # Run the analysis
        while fits.getdata(fitsNnam, ext=1)['status'].count('pending').sum() > 0:
            nbRun = fits.getdata(fitsNnam, ext=1)['status'].count('running').sum()
            while nbRun == n: 
                # Limit the number of parallel runs
                time.sleep(10)
                try:
                    nbRun = fits.getdata(fitsNnam, ext=1)['status'].count('running').sum()
                except:
                    time.sleep(30)
                    nbRun = fits.getdata(fitsNnam, ext=1)['status'].count('running').sum()
            ctRun = fits.getdata(fitsNnam, ext=1)['count'][fits.getdata(fitsNnam, ext=1)['status'] == 'pending']

            # Launch the file
            hdu = fits.open(fitsNnam)
            hdu['SPECTRUM'].data['status'][ctRun[0]] = 'running'
            hdu.writeto(fitsNnam, clobber=True)
            hdu.close()
            fil = os.path.join( self.workpath, 'tmp_'+str(ctRun[0])+'.py' )
            os.popen("nohup python {} &".format(fil))

        return     


    # ============================================================================================= #
    #                                             Plot                                              #
    # ============================================================================================= #
    def pulsPhase(self, ra=None, dec=None, shaded=[]):
        """ Plot an histogram of counts around 1deg of a pulsar
        """

        # --------------------------------------------------------------------------------------------- #
        # Photon selection
        if (self.image is None) or (not os.path.isfile(self.image)):
            print("\t=== Variable self.image needs to point towards a file ===")
            return
        data   = fits.getdata(self.image)
        pulsar = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
        events = SkyCoord(ra=data['RA']*u.degree, dec=data['DEC']*u.degree, frame='icrs') 
        sep    = pulsar.separation(events).deg

        # --------------------------------------------------------------------------------------------- #
        # Plot
        phaseplt         = FermiPlot(savepath='', xsize=8.5, ysize=5)
        phaseplt.figname = os.path.join( os.path.dirname(self.image), 'PulsePhase.pdf')
        phaseplt.xlabel  = r'Phase'
        phaseplt.ylabel  = r'Number of events'
        phaseplt.color   = self.lblue
        phaseplt.shadecol= self.loran
        phaseplt.xmin    = 0.
        phaseplt.xmax    = 1.
        phaseplt.histo(data['PULSE_PHASE'][sep <= 1.], bins=100, shaded=shaded)
        phaseplt.save()

        print("\t=== Figure '{}' created ===".format( os.path.abspath(phaseplt.figname) ))
        return
    def apertLC(self, dt=86400, spec=-2.7, **kwargs):
        """ Plot an aperture photometry light curve
        """

        # --------------------------------------------------------------------------------------------- #
        # Compute and plot aperture photometry light curve
        if not os.path.isfile(self.outapert):
            self._gtAperture(dt=dt, spec=spec)

        aperTab = Table.read(self.outapert)
        timeMJD = self._MET2MJD(aperTab['TIME'])
        tref    = int(np.floor( timeMJD[0] / 100.0)) * 100 # round to lowest hundred

        deltaT = timeMJD[-1]+(aperTab['TIMEDEL'][-1]/86400.)/2. - (timeMJD[0]-(aperTab['TIMEDEL'][0]/86400.)/2.)

        lcplt         = FermiPlot(savepath='', xsize=8.5, ysize=5)
        lcplt.figname = os.path.join(self.workpath, 'AperturePhotometryLC.pdf')
        lcplt.xlabel  = r'Time (MJD $-$ {})'.format(tref)
        lcplt.ylabel  = r'Flux (ph\,cm$^{-2}$\,s$^{-1}$)'
        lcplt.mksize  = 2
        lcplt.scinot  = True
        lcplt.xmin    = timeMJD[0]-(aperTab['TIMEDEL'][0]/86400.)/2.-tref   - 0.05*deltaT
        lcplt.xmax    = timeMJD[-1]+(aperTab['TIMEDEL'][-1]/86400.)/2.-tref + 0.05*deltaT
        lcplt.color   = self.color
        lcplt.plot(x=timeMJD - tref, y=aperTab['RATE'], 
            yerr=aperTab['RATE_ERROR'], xerr=(aperTab['TIMEDEL']/86400.)/2., **kwargs)
        lcplt.save()

        print("\t=== Figure '{}' created ===".format(lcplt.figname))


        # --------------------------------------------------------------------------------------------- #
        # Compute the periodogram
        import scipy.signal as signal

        pmax = timeMJD[-1] - timeMJD[0]
        pmin = aperTab['TIMEDEL'][-1]/86400
        
        freq    = np.linspace(1./pmax, 1./pmin, 10000)
        flux    = np.nan_to_num(aperTab['RATE'].byteswap().newbyteorder())
        # Filter the zeros
        indOK   = flux > 0. 
        time    = timeMJD[indOK]
        flux    = flux[indOK]
        error   = aperTab['RATE_ERROR'][indOK]
        # weighted mean (https://fermi.gsfc.nasa.gov/science/mtgs/workshops/da2010_india/Robin_Corbet_aperture_photometry.pdf)
        wmean   = np.sum(flux/(error**2.)) / np.sum(1./(error**2.)) 
        weighct = (flux - wmean) / (error**2.)
        pgram   = signal.lombscargle(time-time[0], weighct, freq)

        #orbitalmodulation = 0.1 # day.period-1
        #time = np.linspace(0.01, 10, 1000) 
        #flux = np.cos(2.*np.pi*1./orbitalmodulation*time)
        #freq=np.linspace(1./0.01*2*np.pi, 1./10*2*np.pi, 10000)
        #pgram = signal.lombscargle(time, flux, freq)


        normval = time.shape[0]
        power   = np.sqrt( 4*(pgram/normval) )

        pgmplt         = FermiPlot(savepath='', xsize=8.5, ysize=5)
        pgmplt.figname = os.path.join(self.workpath, 'Periodogram.pdf')
        pgmplt.xlabel  = r'Period (days)'
        pgmplt.ylabel  = r'Power'
        pgmplt.mksize  = 0
        pgmplt.lstyle  = '-'
        pgmplt.scinot  = True
        pgmplt.color   = self.color
        pgmplt.xmode   = 'log'
        pgmplt.plot(x=1./freq, y=power, **kwargs)
        pgmplt.save()
        
        return
    def residLike(self):
        """ Plot the residuals after a likelihood analysis
        """

        # --------------------------------------------------------------------------------------------- #
        # Compute the residuals
        if self.csys == 'GAL':
            # Redo some file computations with this coordinate system
            self.outbinexp = os.path.join(self.workpath, 'BinExpMapGAL'+self.suffix+'.fits')
            self.outbincub = os.path.join(self.workpath, 'BinCubeGAL'+self.suffix+'.fits')
            self.outsrcmap = os.path.join(self.workpath, 'SrcMapsGAL'+self.suffix+'.fits')
            self.outresid  = os.path.join(self.workpath, 'ResidGAL'+self.suffix+'.fits')
            self.outresig  = os.path.join(self.workpath, 'ResSigmaGAL'+self.suffix+'.fits')

            self._gtExpmap()
            self._gtBincube()
            self._gtSrcmap()
        else:
            # Nothing to add
            pass
    
        self._gtBinmap()
        self._gtModel()
        # Create the residual count map (count_map - model_map)
        if not os.path.isfile(self.outresid):
            os.popen("farith {} {} {} ops=SUB".format(self.outbinmap, self.outgtmod,
                self.outresid))
        # Create the sigma-residual map (residual_map/sqrt(model_map))
        if not os.path.isfile(self.outresig):
            os.popen("ftpixcalc {} '(a-b)/sqrt(b)' a={} b={}".format(self.outresig,
                self.outbinmap, self.outgtmod))

        # --------------------------------------------------------------------------------------------- #
        # Get the sources to overplot
        srcs = self.getSrc()
        srcs = srcs[(srcs['Separation'] <= 3.) & ([not i.endswith('c') for i in srcs['Name']])]
        # Plot the residuals
        resplt1           = FermiMap()
        resplt1.savepath  = self.workpath
        resplt1.image     = self.outresig
        resplt1.figname   = 'ResSigma.pdf'
        dmin, dmax        = np.abs(resplt1.datamin), resplt1.datamax
        resplt1.datamin   = - min(dmin, dmax)
        resplt1.datamax   = + min(dmin, dmax)
        resplt1.cbarlabel = r'Residual $\sigma$/pixel'
        resplt1.mapSky()
        resplt1.srcSky(srcs['RA'], srcs['DEC'], srcs['Name'])
        resplt1.save()
        print("\t=== Figure '{}' created ===".format( os.path.join(resplt1.savepath, resplt1.figname) ))

        resplt2           = FermiMap()
        resplt2.savepath  = self.workpath
        resplt2.image     = self.outresid
        resplt2.figname   = 'Residuals.pdf'
        dmin, dmax        = np.abs(resplt2.datamin), resplt2.datamax
        resplt2.datamin   = - min(dmin, dmax)
        resplt2.datamax   = + min(dmin, dmax)
        resplt2.cbarlabel = r'Residual counts/pixel'
        resplt2.mapSky()
        resplt2.srcSky(srcs['RA'], srcs['DEC'], srcs['Name'])
        resplt2.save()
        print("\t=== Figure '{}' created ===".format( os.path.join(resplt2.savepath, resplt2.figname) ))
        return
    def ctMap(self, showSrc=False):
        """ Function to plot a count map
        """
        self._gtBinmap()

        mapplt           = FermiMap()
        mapplt.savepath  = self.workpath
        mapplt.image     = self.outbinmap
        mapplt.figname   = 'CMAP.pdf'
        mapplt.cbarlabel = r'Counts'
        mapplt.mapSky()
        if showSrc:
            srcs = self.getSrc()
            srcs = srcs[(srcs['Separation'] <= 3.) & ([not i.endswith('c') for i in srcs['Name']])]
            mapplt.srcSky(srcs['RA'], srcs['DEC'], srcs['Name'])
        mapplt.save()

        print("\t=== Figure '{}' created ===".format( os.path.join(mapplt.savepath, mapplt.figname) ))
        return
    def tsMap(self):
        """ Plot a TS map
        """
        mapplt           = FermiMap()
        mapplt.savepath  = self.workpath
        mapplt.image     = self.outtsmap
        mapplt.figname   = 'TSMAP.pdf'
        mapplt.cbarlabel = r'TS'
        mapplt.mapSky()
        if showSrc:
            srcs = self.getSrc()
            srcs = srcs[(srcs['Separation'] <= 3.) & ([not i.endswith('c') for i in srcs['Name']])]
            mapplt.srcSky(srcs['RA'], srcs['DEC'], srcs['Name'])
        mapplt.save()

        print("\t=== Figure '{}' created ===".format( os.path.join(mapplt.savepath, mapplt.figname) ))
        return
    def lCurve(self):
        """ Plot a light curve
        """ 

        # --------------------------------------------------------------------------------------------- #
        # Read data
        fitsNnam = os.path.join(self.workpath, 'LCresults.fits')
        lcTab    = Table.read(fitsNnam)
        if (self.tstart is not None) and (self.tstop is not None):
            lcTab = lcTab[ (self.tstart <= lcTab['mjd']) & (lcTab['mjd'] <= self.tstop)]
        lcTab = lcTab[lcTab['flux'] != -1.] # avoid undone analyses

        timeMJD = lcTab['mjd']
        tref    = int(np.floor( timeMJD[0] / 100.0)) * 100 # round to lowest hundred
        timeMJD -= tref
        ts      = lcTab['ts']
        detect  = lcTab['ts'] >= self.tsmin
        undet   = lcTab['ts'] < self.tsmin
        flux    = lcTab['flux'][detect]
        fluxerr = lcTab['fluxerr'][detect]
        upperl  = lcTab['upperlim'][undet]
        upperl[upperl == -1.] = 0. # for when it failed
        scale   = 10**int(np.floor(np.log10(  np.mean(  np.concatenate( (flux, upperl), axis=0) )  )))   

        # --------------------------------------------------------------------------------------------- #
        # Plot
        lcplt         = FermiPlot(savepath='', xsize=8.5, ysize=6)
        lcplt.figname = os.path.join(self.workpath, 'LightCurve.pdf')
        lcplt.xlabel  = r'Time (MJD $-$ {})'.format(tref)
        lcplt.ylabel  = [r'Flux ($10^{%d}$ ph\,cm$^{-2}$\,s$^{-1}$)'%(int(np.log10(scale))), r'TS']
        lcplt.hline   = [None, self.tsmin]
        deltaY = max(np.concatenate((flux+fluxerr, upperl), axis=0)) - min(np.concatenate((flux-fluxerr, upperl), axis=0))
        lcplt.ymin    = [(min(np.concatenate((flux-fluxerr, upperl-upperl*0.1), axis=0)) - 0.05*deltaY) / scale, min(ts) - 0.05*(max(ts)-min(ts))]
        lcplt.ymax    = [(max(np.concatenate((flux+fluxerr, upperl), axis=0)) + 0.05*deltaY) / scale, max(ts) + 0.05*(max(ts)-min(ts))]
        deltaX = (timeMJD[-1] + lcTab['mjderr'][-1]) - (timeMJD[0] - lcTab['mjderr'][0])      
        lcplt.xmin    = timeMJD[0] - lcTab['mjderr'][0] - 0.05*deltaX
        lcplt.xmax    = timeMJD[-1] + lcTab['mjderr'][-1] + 0.05*deltaX
        lcplt.fill    = [item for sublist in zip( timeMJD[detect]-lcTab['mjderr'][detect], timeMJD[detect]+lcTab['mjderr'][detect]  ) for item in sublist]
        lcplt.shadecol= self.loran 
        if len(flux) == 0:
            lcplt.mksize  = [2, 2]
            lcplt.ymode   = ['linear', 'linear']
            lcplt.color   = ['gray', 'black']
            lcplt.prop    = [3, 1]
            lcplt.limit   = [True, False]
            lcplt.multiplot(x = [ timeMJD[undet], timeMJD ],
                y    = [ upperl/scale, ts ],
                xerr = [ lcTab['mjderr'][undet], lcTab['mjderr']],
                yerr = [ upperl/scale*0.1, None])
        else:
            lcplt.mksize  = [2, 2, 2]
            lcplt.ymode   = ['linear', 'linear', 'linear']
            lcplt.color   = ['gray', 'black', 'black']
            lcplt.prop    = [3, 1]
            lcplt.limit   = [[True, False], False]
            lcplt.multiplot(x = [ [timeMJD[undet], timeMJD[detect]], timeMJD ],
                y    = [ [upperl/scale, flux/scale], ts ],
                xerr = [ [lcTab['mjderr'][undet], lcTab['mjderr'][detect]], lcTab['mjderr']],
                yerr = [ [upperl/scale*0.1, fluxerr/scale], None])
        lcplt.save()

        print("\t=== Figure '{}' created ===".format(lcplt.figname)) 
        return  
    def indFlux(self):
        """ Plot the power-law index versus the gamma-ray flux
        """
        # --------------------------------------------------------------------------------------------- #
        # Read data
        fitsNnam = os.path.join(self.workpath, 'LCresults.fits')
        lcTab    = Table.read(fitsNnam)
        if (self.tstart is not None) and (self.tstop is not None):
            lcTab = lcTab[ (self.tstart <= lcTab['mjd']) & (lcTab['mjd'] <= self.tstop)]
        lcTab = lcTab[lcTab['flux'] != -1.] # avoid undone analyses

        detect   = lcTab['ts'] >= self.tsmin
        flux     = lcTab['flux'][detect]
        fluxerr  = lcTab['fluxerr'][detect]
        index    = lcTab['index'][detect]
        indexerr = lcTab['indexerr'][detect]
        scale    = 10**int(np.floor(np.log10( np.mean(flux) )))

        # --------------------------------------------------------------------------------------------- #
        # Plot
        indplt         = FermiPlot(savepath='', xsize=8.5, ysize=6)
        indplt.figname = os.path.join(self.workpath, 'IndvsFlux.pdf')
        indplt.xlabel  = r'Flux ($10^{%d}$ ph\,cm$^{-2}$\,s$^{-1}$)'%(int(np.log10(scale)))
        indplt.ylabel  = r'Index'
        indplt.mksize  = 2
        indplt.color   = 'black'
        indplt.plot(x=flux/scale, xerr=fluxerr/scale, y=index, yerr=indexerr)
        indplt.save()

        print("\t=== Figure '{}' created ===".format(indplt.figname)) 
        return
    def spec(self):
        """ Plot a spectrum
        """
        fitsNnam = os.path.join(self.workpath, 'SPECresults.fits')


    # ============================================================================================= #
    #                                           Cyg X-3                                             #
    # ============================================================================================= #
    def cygx3IndFlux(self):
        """ Plot the power-law index versus the gamma-ray flux
        """
        # --------------------------------------------------------------------------------------------- #
        # Read data
        fitsNnam = os.path.join(self.workpath, 'LCresults.fits')
        lcTab    = Table.read(fitsNnam)
        detect   = lcTab['ts'] >= self.tsmin
        lcTab    = lcTab[detect]   

        ind08      = (lcTab['mjd'] > 54700) & (lcTab['mjd'] < 54900) 
        flux08     = lcTab['flux'][ind08]
        fluxerr08  = lcTab['fluxerr'][ind08]
        index08    = lcTab['index'][ind08]
        indexerr08 = lcTab['indexerr'][ind08]

        ind09      = (lcTab['mjd'] > 54900) & (lcTab['mjd'] < 55100) 
        flux09     = lcTab['flux'][ind09]
        fluxerr09  = lcTab['fluxerr'][ind09]
        index09    = lcTab['index'][ind09]
        indexerr09 = lcTab['indexerr'][ind09]

        scale = 10**int(np.floor(np.log10( np.mean(  np.concatenate( (flux08, flux09), axis=0) ) ))) 

        # --------------------------------------------------------------------------------------------- #
        # Plot
        indplt         = FermiPlot(savepath='', xsize=8.5, ysize=6)
        indplt.figname = os.path.join(self.workpath, 'IndvsFlux.pdf')
        indplt.xlabel  = r'Flux ($10^{%d}$ ph\,cm$^{-2}$\,s$^{-1}$)'%(int(np.log10(scale)))
        indplt.ylabel  = r'Index'
        indplt.mksize  = 2
        indplt.color   = self.lblue
        indplt.label   = r'2008'
        indplt.plot(x=flux08/scale, xerr=fluxerr08/scale, y=index08, yerr=indexerr08)
        indplt.color   = self.loran
        indplt.label   = r'2009'
        indplt.plot(x=flux09/scale, xerr=fluxerr09/scale, y=index09, yerr=indexerr09)
        indplt.save()

        print("\t=== Figure '{}' created ===".format(indplt.figname)) 
        return
    def cygx3MWLC(self):
        """ Plot a multi-wavelength light curve for Cyg X-3
        """
        # --------------------------------------------------------------------------------------------- #
        # Fermi data
        fitsNnam = os.path.join(self.workpath, 'LCresults.fits')
        lcTab    = Table.read(fitsNnam)
        if (self.tstart is not None) and (self.tstop is not None):
            lcTab = lcTab[ (self.tstart <= lcTab['mjd']) & (lcTab['mjd'] <= self.tstop)]
        lcTab = lcTab[lcTab['flux'] != -1.] # avoid undone analyses

        timeMJD = lcTab['mjd']
        tref    = int(np.floor( timeMJD[0] / 100.0)) * 100 # round to lowest hundred
        timeMJD -= tref
        ts      = lcTab['ts']
        detect  = lcTab['ts'] >= self.tsmin
        undet   = lcTab['ts'] < self.tsmin
        flux    = lcTab['flux'][detect]
        fluxerr = lcTab['fluxerr'][detect]
        upperl  = lcTab['upperlim'][undet]
        upperl[upperl == -1.] = 0. # for when it failed
        scale   = 10**int(np.floor(np.log10(  np.mean(  np.concatenate( (flux, upperl), axis=0) )  )))

        # --------------------------------------------------------------------------------------------- #
        # X-ray data
        batFile  = os.path.join(self.workpath, 'CygX3_BAT.fits')
        #maxiFile = os.path.join(self.workpath, 'CygX3_MAXI.csv')
        maxiFile = os.path.join(self.workpath, 'CygX3_MAXI.dat')
        asmFile  = os.path.join(self.workpath, 'CygX3_ASM.fits')
        if not os.path.isfile(batFile) or self.clobber:
            os.system('wget http://swift.gsfc.nasa.gov/results/transients/CygX-3.lc.fits -O {}'.format(batFile))
        batTab  = Table.read(batFile)
        if not os.path.isfile(maxiFile) or self.clobber:
            #os.system('wget http://www.maxi.jaxa.jp/obs/agn_etc/data/J2032+409/J2032+409.txt -O {}'.format(maxiFile))
            os.system('wget http://134.160.243.77/star_data/J2032+409/J2032+409_g_lc_1day_all.dat -O {}'.format(maxiFile))
        maxiTab = Table.read(maxiFile, format='ascii') #t, f2-20, e2-20, f2-4, e2-4, f4-10, e4-10, f10-20, e10-20
        if not os.path.isfile(asmFile) or self.clobber:
            os.system('wget https://www.dropbox.com/s/65qrhi1oyifvjfn/CygX3_ASM.fits?dl=0 -O {}'.format(asmFile))
        asmTab  = Table.read(asmFile) #mjd, 1.3-12.2 keV, 1.3-3.0 keV, 3.0-5.0 keV, and 5.0-12.2 keV
        asmTab  = asmTab[asmTab['col1'] > 54500]

        # --------------------------------------------------------------------------------------------- #
        # Radio data
        amiFile  = os.path.join(self.workpath, 'CygX3_AMI.fits')
        ovroFile = os.path.join(self.workpath, 'CygX3_OVRO.fits')
        if not os.path.isfile(amiFile) or self.clobber:
            os.system('wget https://www.dropbox.com/s/bz9xbdbq6hbrant/AMI_2008_14.fits?dl=0 -O {}'.format(amiFile))
        amiTab = Table.read(amiFile)
        if not os.path.isfile(ovroFile) or self.clobber:
            os.system('wget https://www.dropbox.com/s/rs7xlztd66j6fej/CygX3_OVRO.fits?dl=0 -O {}'.format(ovroFile))
        ovroTab = Table.read(ovroFile)
        ovroOff = - 0.124

        # --------------------------------------------------------------------------------------------- #
        # Plot
        lcplt         = FermiPlot(savepath='', xsize=8.5, ysize=17)
        lcplt.figname = os.path.join(self.workpath, 'CygX3_MWLC.pdf')
        lcplt.xlabel  = r'Time (MJD $-$ {})'.format(tref)
        lcplt.ylabel  = [r'Flux density (Jy)', r'Count rate', r'Rate (cm$^{-2}$\,s$^{-1}$)', r'Flux ($10^{%d}$ ph\,cm$^{-2}$\,s$^{-1}$)'%(int(np.log10(scale))), r'TS']
        lcplt.label   = [r'AMI', r'OVRO', r'ISS/MAXI ($\times 30$ ct\,cm$^{-2}$\,s$^{-1}$)', r'RXTE/ASM (ct\,s$^{-1}$)', r'\textit{Swift}/BAT', None, r'\textit{Fermi}/LAT', None]
        lcplt.hline   = [None, None, None, None, self.tsmin]

        deltaY = max(np.concatenate((flux+fluxerr, upperl), axis=0)) - min(np.concatenate((flux-fluxerr, upperl), axis=0))
        lcplt.ymin    = [5.e-2,
            None,
            -0.01,
            (min(np.concatenate((flux-fluxerr, upperl-upperl*0.1), axis=0)) - 0.05*deltaY) / scale,
            min(ts) - 0.05*(max(ts)-min(ts))]
        lcplt.ymax    = [3.e1,
            None,
            0.08,
            (max(np.concatenate((flux+fluxerr, upperl), axis=0)) + 0.05*deltaY) / scale, 
            max(ts) + 0.05*(max(ts)-min(ts))]
        deltaX = (timeMJD[-1] + lcTab['mjderr'][-1]) - (timeMJD[0] - lcTab['mjderr'][0])      
        lcplt.xmin    = timeMJD[0] - lcTab['mjderr'][0] - 0.05*deltaX
        lcplt.xmax    = timeMJD[-1] + lcTab['mjderr'][-1] + 0.05*deltaX

        lcplt.fill    = [item for sublist in zip( timeMJD[detect]-lcTab['mjderr'][detect], timeMJD[detect]+lcTab['mjderr'][detect]  ) for item in sublist]
        lcplt.shadecol= self.loran 

        lcplt.mksize  = [1, 1, 1, 1, 1, 2, 2, 2]
        lcplt.ymode   = ['log', 'log', 'linear', 'linear', 'linear', 'linear', 'linear', 'linear']
        lcplt.color   = ['black', self.lblue, 'black', self.lblue, 'black', 'gray', 'black', 'black']
        lcplt.prop    = [3, 3, 3, 3, 1]
        lcplt.limit   = [[False, False], [False, False], False, [True, False], False]
        lcplt.multiplot(x = [ [amiTab['MJD']-tref, ovroTab['mjd']-tref],
            [maxiTab['col1']-tref, asmTab['col1']-tref],
            batTab['TIME']+np.ones(len(batTab['TIME']))*0.5-tref,
            [timeMJD[undet], timeMJD[detect]],
            timeMJD ],
            y    = [ [amiTab['Jy'], ovroTab['flux']+ovroOff],
            [maxiTab['col4']*30, asmTab['col6']], 
            batTab['RATE'],
            [upperl/scale, flux/scale],
            ts ],
            xerr = [ [None, None],
            [None, None],
            np.ones(len(batTab['TIME']))*0.5,
            [lcTab['mjderr'][undet], lcTab['mjderr'][detect]],
            lcTab['mjderr']],
            yerr = [ [None, None],
            [maxiTab['col5']*30, asmTab['col7']],
            batTab['ERROR'],
            [upperl/scale*0.1, fluxerr/scale],
            None])
        lcplt.save()

        print("\t=== Figure '{}' created ===".format(lcplt.figname)) 
        return 


    # ============================================================================================= #
    #                                      Internal functions                                       #
    # ============================================================================================= #
    def _readFT1(self):
        """ Read a FT1 file and extract informations from header each time the ft1 attribute is updated.
        """

        mainHead = fits.getheader(self._ft1, ext=0)
        dataHead = fits.getheader(self._ft1, ext=1)
        irfsPass = {'P8R2': 'P8R2_SOURCE_V6'}

        # --------------------------------------------------------------------------------------------- #
        # Fill the class attributes with informations from the data file
        if os.path.dirname(self.ft1) == '':
            # Get the absolute name to the current path
            self.datapath = os.getcwd()
            self.workpath = os.getcwd() # by default it's the same
        else:
            self.datapath = os.path.dirname(self.ft1)
            self.workpath = os.getcwd() #os.path.dirname(self.ft1)
        if 'NoPulse' in os.path.basename(self.ft1):
            self.frac = float(os.path.basename(self.ft1).split('_')[-1].replace('.fits', '')) 
        self.tstart   = Time(mainHead['DATE-OBS'], format='isot', scale='utc')
        self.tstop    = Time(mainHead['DATE-END'], format='isot', scale='utc')
        self.metstart = mainHead['TSTART']
        self.metstop  = mainHead['TSTOP']
        for i in dataHead.keys():
            if 'DSTYP' in i:
                if 'ENERGY' in dataHead[i]:
                    self.emin = float(dataHead['DSVAL'+i[-1]].split(':')[0])
                    self.emax = float(dataHead['DSVAL'+i[-1]].split(':')[1])
                    break
        self.nevents  = dataHead['NAXIS2']
        self.passver  = dataHead['PASS_VER']
        try:
            self.irf  = irfsPass[self.passver]
        except:
            print("\t=== self.irf needs to be filled manually===")
        for i in dataHead.keys():
            if isinstance(dataHead[i], str):
                if 'BIT_MASK(EVENT_CLASS' in dataHead[i]:
                    self.evclass  = dataHead['DSTYP'+i[-1]].split(',')[1]
                    break
        for i in dataHead.keys():
            if isinstance(dataHead[i], str):
                if 'BIT_MASK(EVENT_TYPE' in dataHead[i]:
                    self.evtype = dataHead['DSTYP'+i[-1]].split(',')[1]
                    break
            else:
                self.evtype = None
        firstFound   = False
        for i in dataHead.keys():
            if (dataHead[i] == 'POS(RA,DEC)') and (not firstFound):
                pointingInfo = dataHead['DSVAL'+i[-1]].split('(')[1].split(')')[0].split(',')
                firstFound   = True 
            elif (dataHead[i] == 'POS(RA,DEC)') and (firstFound):
                # The FT1 has two positions informations (classic...)
                # Need to remove the second one 
                print("\t=== Multiple central postions found, removing the secondary ===")
                hdus = fits.open(self._ft1)
                hdus['EVENTS'].header.remove('DSTYP' + i[-1])
                hdus['EVENTS'].header.remove('DSUNI' + i[-1])
                hdus['EVENTS'].header.remove('DSVAL' + i[-1])
                # Rename the remaining keywords
                for j in dataHead.keys():
                    if j[:5] in ['DSVAL', 'DSUNI', 'DSTYP', 'DSREF']:
                        if int(j[-1]) > int(i[-1]):
                            hdus['EVENTS'].header.rename_keyword(j, j[:5]+str(int(j[-1])-1), force=True)
                hdus['EVENTS'].header.set('NDSKEYS', hdus['EVENTS'].header['NDSKEYS']-1)
                hdus.writeto(self._ft1, clobber=True) # replace the existing FT1
                break
            else:
                pass
        self.ra  = float(pointingInfo[0])
        self.dec = float(pointingInfo[1])
        self.rad = float(pointingInfo[2])
        if 'FT1_filtered' not in os.path.basename(self.ft1):
            # It's not at least a processed filtered FT1 file
            self.fermicat = os.path.join(self.datapath, 'gll_psc_v16.fit')
            self.model    = os.path.join(self.datapath, os.path.basename(self.ft1[:-5])+'_Model.xml')
        return
    def _initNames(self):
        """ Initialize attributes based on the current data
        """
        self.outselect = os.path.join(self.workpath, 'FT1_selected'+self.suffix+'.fits')
        self.outmktime = os.path.join(self.workpath, 'FT1_filtered'+self.suffix+'.fits')
        self.outltcube = os.path.join(self.workpath, 'LtCube'+self.suffix+'.fits')
        self.outbincub = os.path.join(self.workpath, 'BinCube'+self.suffix+'.fits')
        self.outbinmap = os.path.join(self.workpath, 'CMAP'+self.suffix+'.fits')
        self.outbinexp = os.path.join(self.workpath, 'BinExpMap'+self.suffix+'.fits')
        self.outexpmap = os.path.join(self.workpath, 'ExpMap'+self.suffix+'.fits')
        self.outsrcmap = os.path.join(self.workpath, 'SrcMaps'+self.suffix+'.fits')
        self.outgtlike = os.path.join(self.workpath, 'Results'+self.suffix+'.dat')
        self.outmodel  = os.path.join(self.workpath, 'OutModel'+self.suffix+'.xml')
        self.outapert  = os.path.join(self.workpath, 'LC_ApPhoto'+self.suffix+'.fits')
        self.outgtmod  = os.path.join(self.workpath, 'GtModel'+self.suffix+'.fits')
        self.outresid  = os.path.join(self.workpath, 'Resid'+self.suffix+'.fits')
        self.outresig  = os.path.join(self.workpath, 'ResSigma'+self.suffix+'.fits')
        self.outtsmap  = os.path.join(self.workpath, 'TSMmap'+self.suffix+'.fits')
        return
        # self.outfind    = self.dir + self.src + '_FindSrc'+self.suffix+'.txt'
    def _progressBar(self, percent, printEvery=10):
        """ Print a progression bar every 10%
        """
        floor = int(percent)
        sys.stdout.write('\r' * (floor + 9))
        sys.stdout.write('[')
        sys.stdout.write('=' * (floor/printEvery))
        sys.stdout.write('>] {:02.2f}%'.format(percent))
        sys.stdout.flush()
    def _MET2MJD(self, met):
        """ Convert Fermi's Mission Elapsed Time into Modified Julian Date
        """ 
        return (met - 10676868.60)/86400. + 52033.575
    def _mergeDiffrsp(self):
        """ Merge all the FT1 files with diffuse columns in them into one single file
        """
        # --------------------------------------------------------------------------------------------- #
        # Create a file listing all the FT1 files
        ft1Files = np.array(glob.glob(os.path.join(self.workpath, '*_Chk_*.fits')))
        tmpFT1   = os.path.join(self.workpath, 'ft1.lis')
        chkNb    = np.array([int(os.path.basename(ff).split('_')[-1][:-5]) for ff in ft1Files])
        sortInd  = np.argsort(chkNb)
        wfil     = open( tmpFT1, 'w')
        for f in ft1Files[sortInd]:
            wfil.write(f + '\n')
        wfil.close()

        # --------------------------------------------------------------------------------------------- #
        # Merging everything
        self.outselect = os.path.join(self.workpath, 'FT1_Diffuse'+self.suffix+'.fits')
        self._gtSelect(data=tmpFT1)
        os.remove(tmpFT1)

        return
    def _lcResults(self):
        """ Update the result FITS file with LC measurement
        """
        fitsNnam = os.path.join(self.workpath, 'LCresults.fits')
        noCt     = int( self.suffix.split('_')[-1] ) 
        
        hdu = fits.open(fitsNnam)
        hdu['LIGHTCURVE'].data['flux'][noCt] = self.getResult(param='Flux')[0]
        hdu['LIGHTCURVE'].data['fluxerr'][noCt] = self.getResult(param='Flux')[1]
        hdu['LIGHTCURVE'].data['index'][noCt] = self.getResult(param='Index')[0]
        hdu['LIGHTCURVE'].data['indexerr'][noCt] = self.getResult(param='Index')[1]
        hdu['LIGHTCURVE'].data['ts'][noCt] = self.getResult(param='TS value')[0]
        hdu['LIGHTCURVE'].data['upperlim'][noCt] = self.getResult(param='upplim')[0] 
        hdu['LIGHTCURVE'].data['status'][noCt] = 'done'
        hdu.writeto(fitsNnam, clobber=True)
        hdu.close()
    def _specResults(self):
        """ Update the result FITS file with Spectrum measurement
        """
        fitsNnam = os.path.join(self.workpath, 'SPECresults.fits')
        noCt     = int( self.suffix.split('_')[-1] ) 
        
        hdu = fits.open(fitsNnam)
        hdu['SPECTRUM'].data['flux'][noCt] = self.getResult(param='Flux')[0]
        hdu['SPECTRUM'].data['fluxerr'][noCt] = self.getResult(param='Flux')[1]
        hdu['SPECTRUM'].data['index'][noCt] = self.getResult(param='Index')[0]
        hdu['SPECTRUM'].data['indexerr'][noCt] = self.getResult(param='Index')[1]
        hdu['SPECTRUM'].data['ts'][noCt] = self.getResult(param='TS value')[0]
        hdu['SPECTRUM'].data['upperlim'][noCt] = self.getResult(param='upplim')[0] 
        hdu['SPECTRUM'].data['status'][noCt] = 'done'
        hdu.writeto(fitsNnam, clobber=True)
        hdu.close()
    def _gtSelect(self, data=None):
        if data is None:
            data = self.ft1
        if (os.path.isfile(self.outselect) or os.path.isfile(self.outmktime)) and not self.clobber:
            print("\t=== gtselect already performed. ===")
            return
        else:
            os.popen("gtselect infile={} zmax={} ra={} dec={} rad={} emin={} emax={}\
                evclass={} evtype={} outfile={} tmin={} tmax={}".format(data, self.zmax,
                self.ra, self.dec, self.rad, self.emin, self.emax, self.evclass, self.evtype,
                self.outselect, self.metstart, self.metstop))
            return
    def _gtMktime(self):
        """ Select GTIs
        """
        if os.path.isfile(self.outmktime):
            print("\t=== '{}' already exists ===".format(self.outmktime))
            return
        else:
            if not os.path.isfile(self.outselect):
                self._gtSelect()

        os.popen("gtmktime evfile={} scfile={} outfile={} filter='(DATA_QUAL>0)&&(LAT_CONFIG==1)'\
            roicut='no'".format(self.outselect, self.ft2, self.outmktime))
        return
    def _gtDiffrsp(self):
        """ Add the diffuse colmuns to an event file
        """
        if not os.path.isfile(self.outmktime):
            self._gtMktime()
        os.popen("gtdiffrsp evfile={} scfile={} srcmdl={} irfs={} convert=Y chatter=4\
            evtype={}".format(self.outmktime, self.ft2, self.diffModel, self.irf, self.evtype))
        return
    def _gtBincube(self):
        """ Create a 3-D (binned) counts map
        """
        if os.path.isfile(self.outbincub):
            print("\t=== '{}' already exists ===".format(self.outbincub))
            return
        else:
            if not os.path.isfile(self.outmktime):
                self._gtMktime()

        # Image width must be comprised within the acceptance cone
        imWidth = int( np.floor(self.rad* 2**(0.5)) )
        imWipix = int(imWidth / self.binsz)
        # 10 energy bins per decade
        enBins  = int( np.ceil( (np.log10(self.emax) - np.log10(self.emin)) * 10 ) )

        # Coordinate system
        if self.csys == 'GAL':
            center_icrs = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree, frame='icrs')
            self.ra  = center_icrs.galactic.l.deg
            self.dec = center_icrs.galactic.b.deg
 
        #os.popen("gtbin evfile={} scfile={} outfile={} algorithm=CCUBE ebinalg=LOG\
        #    emin={} emax={} enumbins={} nxpix={} nypix={} binsz={} coordsys={} xref={}\
        #    yref={} axisrot=0 proj=AIT".format(self.outmktime, self.ft2, self.outbincub,
        #    self.emin, self.emax, enBins, imWipix, imWipix, self.binsz, self.csys, self.ra, self.dec))

        # Don't know why but os.popen() doesn't work for this one
        cmd = ['gtbin', 'evfile='+self.outmktime, 'scfile='+self.ft2, 'outfile='+self.outbincub,
            'algorithm=CCUBE', 'ebinalg=LOG', 'emin='+str(self.emin), 'emax='+str(self.emax),
            'enumbins='+str(enBins), 'nxpix='+str(imWipix), 'nypix='+str(imWipix), 'binsz='+str(self.binsz),
            'coordsys='+self.csys, 'xref='+str(self.ra), 'yref='+str(self.dec), 'axisrot=0', 'proj=AIT'] 
        with open( os.path.join(self.workpath, 'gtbin_out.txt'), 'w') as out:
            subprocess.call(cmd, stdout=out) 

        if self.csys == 'GAL':
            self.ra  = center_icrs.galactic.ra.deg
            self.dec = center_icrs.galactic.dec.deg
        return
    def _gtBinmap(self):
        """ Make a counts map from the event data
        """
        if os.path.isfile(self.outbinmap) and (not self.clobber):
            print("\t=== '{}' already exists ===".format(self.outbinmap))
            return
        else:
            if not os.path.isfile(self.outmktime):
                self._gtMktime()

        # Image width must be comprised within the acceptance cone
        imWidth = int( np.floor(self.rad* 2**(0.5)) ) # deg
        imWipix = int(imWidth / self.binsz)

        # Coordinate system
        if self.csys == 'GAL':
            center_icrs = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree, frame='icrs')
            self.ra  = center_icrs.galactic.l.deg
            self.dec = center_icrs.galactic.b.deg

        os.popen("gtbin evfile={} scfile=none outfile={} algorithm=CMAP emin={}\
            emax={} nxpix={} nypix={} binsz={} coordsys={} xref={} yref={} axisrot=0\
            proj=AIT".format(self.outmktime, self.outbinmap, self.emin, self.emax,
            imWipix, imWipix, self.binsz, self.csys, self.ra, self.dec))

        if self.csys == 'GAL':
            self.ra  = center_icrs.ra.deg
            self.dec = center_icrs.dec.deg
        return
    def _gtLtcube(self):
        """ Compute livetimes and exposure
        """
        if os.path.isfile(self.outltcube):
            print("\t=== '{}' already exists ===".format(self.outltcube))
            return
        else:
            if not os.path.isfile(self.outmktime):
                self._gtMktime()

        os.popen("gtltcube evfile={} scfile={} sctable='SC_DATA' outfile={} dcostheta=0.025\
            binsz=1.0 tmin=0.0 tmax=0.0 zmax={}".format(self.outmktime, self.ft2, self.outltcube,
            self.zmax))
        return
    def _gtExpmap(self):
        """ Compute exposure map
        """
        if self.mode   == 'binned':
            if os.path.isfile(self.outbinexp):
                print("\t=== '{}' already exists ===".format(self.outbinexp))
                return
            else:
                if not os.path.isfile(self.outltcube):
                    self._gtLtcube()

            # Image width must encompass sources up tp 10deg from ROI + 10 additional deg
            imWidth = self.rad + 30
            imWipix = int(np.ceil( (imWidth/self.binsz) / 100.0)) * 100 # round to next hundred
            # 10 energy bins per decade
            enBins  = int( np.ceil( (np.log10(self.emax) - np.log10(self.emin)) * 10 ) )

            # Coordinate system
            if self.csys == 'GAL':
                center_icrs = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree, frame='icrs')
                self.ra  = center_icrs.galactic.l.deg
                self.dec = center_icrs.galactic.b.deg

            os.popen("gtexpcube2 infile={} cmap=none outfile={} evtype={} irfs={} nxpix={}\
                nypix={} binsz={} coordsys={} xref={} yref={} axisrot=0 proj=AIT emin={}\
                emax={} enumbins={}".format(self.outltcube, self.outbinexp, self.evtype, self.irf,
                imWipix, imWipix, self.binsz, self.csys, self.ra, self.dec, self.emin, self.emax,
                enBins))

            if self.csys == 'GAL':
                self.ra  = center_icrs.ra.deg
                self.dec = center_icrs.dec.deg
            return
        elif self.mode == 'unbinned':
            if os.path.isfile(self.outexpmap):
                print("\t=== '{}' already exists ===".format(self.outexpmap))
                return
            else:
                if not os.path.isfile(self.outltcube):
                    self._gtLtcube()
            # Image width must encompass sources up tp 10deg from ROI + 10 additional deg
            imWidth = int(self.rad + 20)
            imWipix = int(imWidth/0.5) # half deg pix
            # 5 energy bins per decade
            enBins  = int( np.ceil( (np.log10(self.emax) - np.log10(self.emin)) * 5 ) )

            os.popen("gtexpmap evfile={} expcube={} outfile={} evtable='EVENTS' scfile={}\
                irfs={} evtype={} srcrad={} nlong={} nlat={} nenergies={}".format(self.outmktime,
                    self.outltcube, self.outexpmap, self.ft2, self.irf, self.evtype, imWidth,
                    imWipix, imWipix, enBins))
            return
        else:
            print("\t=== Mode '{}' unrecognized ('binned'/'unbinned') ===".format(self.mode))
            return
            
        # Correcting exposure in case of nopulse analysis
        if self.frac is not None:
            for fil in [self.outexpmap, self.outbinexp]:
                tmpName = os.path.join(self.workpath, 'tmpexpmap.fits')
                os.popen("fcarith {} {} {} MUL".format(fil, self.frac, tmpName))
                os.popen("fappend {}+1 {}".format(fil, tmpName))
                os.popen("fappend {}+2 {}".format(fil, tmpName))
            os.remove(fil)
            os.rename(tmpName, fil)

        return
    def _gtSrcmap(self):
        """ Compute source map
        """
        if os.path.isfile(self.outsrcmap):
            print("\t=== '{}' already exists ===".format(self.outsrcmap))
            return
        else:
            if not os.path.isfile(self.outbincub):
                self._gtBincube()
            if not os.path.isfile(self.outbinexp):
                self._gtExpmap()

        os.popen("gtsrcmaps expcube={} cmap={} scfile={} srcmdl={} bexpmap={}\
            outfile={} irf={} ptsrc=no".format(self.outltcube, self.outbincub,
            self.ft2, self.model, self.outbinexp, self.outsrcmap, self.irf))
        return
    def _gtLike(self):
        """ Performs unbinned or binned likelihood analysis of LAT data
        """
        if os.path.isfile(self.outgtlike):
            print("\t=== '{}' already exists ===".format(self.outgtlike))
            return

        specFile = os.path.join(self.workpath, 'SrcList_cntspec'+self.suffix+'.fits') 
        specLog  = os.path.join(self.workpath, 'SrcList_cntspec'+self.suffix+'.log')    
        if self.mode   == 'binned':
            if not os.path.isfile(self.outsrcmap):
                self._gtSrcmap()
            os.popen("gtlike refit=no plot=no statistic=BINNED optimizer=NewMinuit ftol=1e-2 sfile={}\
                cmap={} expcube={} bexpmap={} srcmdl={} irfs={} results={} specfile={} > {}"
                .format(self.outmodel, self.outsrcmap, self.outltcube, self.outbinexp, 
                self.model, self.irf, self.outgtlike, specFile, specLog))
            return
        elif self.mode == 'unbinned':
            if not os.path.isfile(self.outexpmap):
                self._gtExpmap()
            os.popen("gtlike expcube={} evfile={} expmap={} sfile={} statistic=UNBINNED\
                optimizer=NewMinuit ftol=1e-2 scfile={} tsmin=yes srcmdl={} irfs={} results={}\
                specfile={} > {}".format(self.outltcube, self.outmktime, self.outexpmap,
                self.outmodel, self.ft2, self.model, self.irf, self.outgtlike, specFile, specLog))
            return
        else:
            print("\t=== Mode '{}' unrecognized ('binned'/'unbinned') ===".format(self.mode))
            return
    def _uppLim(self):
        """ Return the 95 per cent upper limit on the flux in units of ph/cm2/s.
        """
        if self.getResult(param='TS value')[0] >= self.tsmin:
            print("\t=== TS value {} is above TSmin {}, no need to compute an upperlimit  ==="
                .format(self.getResult(param='TS value')[0], self.tsmin))
            return

        from UpperLimits import UpperLimits
        import UnbinnedAnalysis as UA
    
        like = UA.unbinnedAnalysis(evfile=self.outmktime, scfile=self.ft2, expmap=self.outexpmap,
                expcube=self.outltcube, irfs=self.irf, optimizer="NewMinuit", srcmdl=self.model)
        like.fit(0)
        ul = UpperLimits(like)

        try:
            upp, norm=ul['TARGET'].bayesianUL(emin=self.emin, emax=self.emax, cl=0.95) 
        except:
            upp = -1
        wf = open(self.outgtlike, 'a')
        wf.write("\nUpper limit on source 'TARGET': {} ph/cm2/s.".format(upp))
        wf.close()
        return
    def _gtModel(self):
        """ Create a model map of the region based on the fit parameters
        """
        if os.path.isfile(self.outgtmod):
            print("\t=== '{}' already exists ===".format(self.outgtmod))
            return
        else:
            if not os.path.isfile(self.outsrcmap):
                self._gtSrcmap()
            if not os.path.isfile(self.outmodel):
                self.gtLike()

        os.popen("gtmodel srcmaps={} srcmdl={} outfile={} irfs={} expcube={}\
            bexpmap={}".format(self.outsrcmap, self.outmodel, self.outgtmod,
            self.irf, self.outltcube, self.outbinexp))
        return
    def _gtTSmap(self):
        """ Create a residual TS map
        """
        if os.path.isfile(self.outtsmap):
            # Already exists
            return

        if self.csys == 'GAL':
            center_icrs = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree, frame='icrs')
            self.ra  = center_icrs.galactic.l.deg
            self.dec = center_icrs.galactic.b.deg

        model = os.path.join(self.workpath, 'TSmapModel.xml') 
        rfil  = open(self.outmodel, 'r')
        wfil  = open(model, 'w')
        isSrc = False
        isDif = False
        for line in rfil:
            if (isSrc) and ('<source name' in line):
                # Arrived to a new source, restart copying
                isSrc = False
            if (isDif) and ('<source name' in line) and ('PointSource' in line):
                isDif = False
            if 'TARGET' in line:
                isSrc = True
            if ('<source name="gll_iem_v06"' in line) or ('<source name="iso_source_v06"' in line): 
                isDif = True
            
            if isSrc:
                # Do not copy the Target model to make it appear in the TS map
                pass
            else:
                if isDif:
                    # Leave Diffuse model normalizations free
                    wfil.write(line)
                else:
                    # Make sur the gtlike output source model has all source parameters fixed
                    wfil.write(line.replace('free="1"', 'free="0"'))
        rfil.close()
        wfil.close()

        # Launch the gttsmap tool 
        if self.mode == 'binned':
            os.popen("gttsmap evfile={} scfile={} bexpmap={} expcube={} cmap={} srcmdl={}\
                outfile={} evtype={} irfs=CALDB optimizer=NewMinuit statistic=BINNED ftol=1e-2\
                coordsys={} proj=AIT nxpix={} nypix={} binsz={} xref={} yref={}".format(self.outmktime,
                self.ft2, self.outbinexp, self.outltcube, self.outbincub, model, self.outtsmap, self.evtype,
                self.csys, self.imwid, self.imwid, self.binsz, self.ra, self.dec))
        elif self.mode == 'unbinned':
            os.popen("gttsmap evfile={} scfile={} expmap={} expcube={} srcmdl={}\
                outfile={} evtype={} irfs=CALDB optimizer=NewMinuit statistic=UNBINNED ftol=1e-2\
                coordsys={} proj=AIT nxpix={} nypix={} binsz={} xref={} yref={}".format(self.outmktime,
                self.ft2, self.outexpmap, self.outltcube, model, self.outtsmap, self.evtype,
                self.csys, self.imwid, self.imwid, self.binsz, self.ra, self.dec))
        else:
            return

        if self.csys == 'GAL':
            self.ra  = center_icrs.ra.deg
            self.dec = center_icrs.dec.deg
        return
    def _gtAperture(self, dt, spec):
        """ Compute an aperture photometry lightcurve
        """
        self.rad = 1
        if os.path.isfile(self.outapert):
            print("\t=== '{}' already exists ===".format(self.outapert))
            return

        self._gtSelect()
        self._gtMktime()
        # Use gtbin in 'Light curve mode':
        os.popen("gtbin algorithm=LC evfile={} outfile={} scfile={}\
            tbinalg=LIN tstart={} tstop={} dtime={}".format(self.outmktime,
            self.outapert, self.ft2, self.metstart, self.metstop, dt))
        # Add the exposure column
        os.popen("gtexposure infile={} scfile={} irfs={} srcmdl='none' specin={}"
            .format(self.outapert, self.ft2, self.irf, spec))
        # Correct for solar system barycenter
        os.popen("gtbary evfile={} outfile={} scfile={} ra={} dec={} tcorrect=BARY"
            .format(self.outapert, self.outapert, self.ft2, self.rad, self.dec))
        # Convert and compute errors
        os.popen("ftcalc {} {} RATE 'counts/exposure' clobber=YES".format(self.outapert, self.outapert))
        os.popen("ftcalc {} {} RATE_ERROR 'error/exposure' clobber=YES".format(self.outapert, self.outapert))

        return
        

# ============================================================================================= #
# ============================================================================================= #
# ============================================================================================= #
class FermiPlot(object):
    """ A class to do fancy plots
    """
    def __init__(self, savepath=None, xsize=8, ysize=8):
        def cm2inch(size):
            """ Return the conversion in inches
            """
            return size / 2.54

        # --------------------------------------------------------------------------------------------- #
        # Attributes / figure properties
        self.savepath  = savepath
        if self.savepath is None: 
            self.savepath = os.getcwd() + '/' # By default, save in the current directory
        self.figname   = 'MyFigure.pdf'
        self.ymode     = 'linear' # or 'log' 
        self.xmode     = 'linear' 
        self.ftsize1   = 9
        self.ftsize2   = 6
        self.xlabel    = 'X Label'
        self.ylabel    = 'Y Label'
        self.xsize     = xsize # in cm
        self.ysize     = ysize
        self.color     = 'black'
        self.color2    = 'black'
        self.lstyle    = ''
        self.mksize    = 4 
        self.mktype    = 'o'
        self.fill      = None
        self.shadecol  = 'gray'
        self.mag       = False
        self.labelx    = -0.1
        self.lwdth     = 1
        # Attributes / plot properties
        self.xmin      = None
        self.xmax      = None
        self.ymin      = None
        self.ymax      = None
        self.alpha     = 1.0
        self.limit     = False
        self.yerr      = None
        self.xerr      = None
        self.shaded    = []
        self.hline     = None
        self.label     = None
        self.loc       = "best"
        self.ncol      = 1
        self.fillbelow = False
        self.scinot    = False
        self.raster    = False
        self.mkfill    = True
        self.prop      = []

        # --------------------------------------------------------------------------------------------- #
        # Initialization of plotting routines
        self.f, self.ax = plt.subplots(figsize=(cm2inch(self.xsize), cm2inch(self.ysize)))

    # ============================================================================================= #
    #                                         Normal plot                                           #
    # ============================================================================================= #
    def plot(self, x, y, **kwargs):
        # --------------------------------------------------------------------------------------------- #
        # Attributes
        self._evalKwargs(kwargs)

        plt.yscale(self.ymode)
        plt.xscale(self.xmode)

        # --------------------------------------------------------------------------------------------- #
        # Plot upper limits (down-arrows)
        if self.limit: 
            self.ax.errorbar(x, y, xerr=self.xerr, yerr=[self.yerr, np.zeros(len(self.yerr))],
                fmt='none', ecolor=self.color, elinewidth=0.5, alpha=self.alpha,
                capsize=0, barsabove=False, lolims=False, uplims=False, xlolims=False,
                xuplims=False, errorevery=1, capthick=None, rasterized=self.raster)
            self.ax.plot(x, y - self.yerr, marker='v', color=self.color, alpha=self.alpha,
                markersize=self.mksize, linestyle='', markeredgecolor=self.color, rasterized=self.raster) 
        # Normal plot
        else:
            # Marker empty or not
            if self.mkfill:
                mkfacecol = self.color
            else: 
                mkfacecol = 'none' 
            # Plot
            self.ax.errorbar(x, y, yerr=self.yerr, xerr=self.xerr, fmt=self.mktype, ecolor=self.color,
                elinewidth=0.5, capsize=0, linestyle=self.lstyle, markerfacecolor=mkfacecol,
                markeredgecolor=self.color, color=self.color, markersize=self.mksize, label=self.label,
                alpha=self.alpha, barsabove=False, errorevery=1, capthick=None, linewidth=self.lwdth,
                rasterized=self.raster)
            # Label of the plot and legend box
            if self.label is not None:
                self.ax.legend(loc=self.loc, prop={'size':self.ftsize2}, frameon=True, numpoints=1,
                ncol=self.ncol)

        # --------------------------------------------------------------------------------------------- #
        # Add features
        if self.hline is not None:
            # Draw a horizontal line:
            self.ax.axhline(y=self.hline, color=self.color2, linestyle=':') 
        if len(self.shaded) != 0:
            # Draw shaded region between two x points:
            for i in range(len(self.shaded)/2) :
                self.ax.axvspan(self.shaded[i*2], self.shaded[i*2+1], facecolor=self.shadecol,
                edgecolor='none', linewidth=0., zorder=-10, alpha=0.5, rasterized=self.raster)
        if self.fillbelow:
            # Colorize the area below the plotted curve
            ymin, ymax = self.ax.get_ylim()
            
        self._plotDisplay()
    def vline(self, value, zorder=1):
        """ Draw a vertical line
        """
        self.ax.axvline(x=value, color=self.color, linestyle=self.lstyle) 
    def plotBetween(self, x, yBottom, yTop, **kwargs):
        """ Fill an area
        """
        self._evalKwargs(kwargs)

        self.ax.fill_between(x, yBottom, yTop, facecolor=self.color, edgecolor='none', alpha=0.5,
            rasterized=self.raster, zorder=-10)

        self._plotDisplay()
    def multiplot(self, x, y, **kwargs):
        """ Draw several plots with the same X axis.

        Parameters
        ----------
        self.prop: list
            e.g.: [2,1] --> first subplot's height will be twice the second's 
        """

        # --------------------------------------------------------------------------------------------- #
        # Attributes
        self._evalKwargs(kwargs)
        # Remove the previous and create the new framework
        plt.delaxes(self.ax)
        count = 0
        colcount = 0
        # Get the min and max values of the X-axis
        xmin = []
        xmax = []
        for i in range( len(x) - 1):
            if hasattr(x[i][0], "__len__"):
                for j in range( len(x[i]) - 1):
                    xmin.append( min(x[i][j]) )
                    xmax.append( max(x[i][j]) )
            else:
                xmin.append( min(x[i]) )
                xmax.append( max(x[i]) )
        if self.xmin is not None:
            xmin = [self.xmin]
        if self.xmax is not None:
            xmax = [self.xmax]
        deltaX = max(xmax) - min(xmin)
        xmin   = min(xmin) - 0.05*deltaX
        xmax   = max(xmax) + 0.05*deltaX

        # --------------------------------------------------------------------------------------------- #
        # Iterate over the number of subplots 
        for nSP in range( len(self.prop) ):
            # --------------------------------------------------------------------------------------------- #
            # Initialize the subplot properties
            self.ax = plt.subplot2grid( (sum(self.prop), 1), (count, 0), rowspan=self.prop[nSP])
            count   += self.prop[nSP] # Keep track of the size of the plot
            # Extract the errors if any are given
            if self.yerr is not None:
                yerrSP = self.yerr[nSP]
            if self.xerr is not None:
                xerrSP = self.xerr[nSP] 
            # Set the y-axis and x-axis scales
            try:
                ymode  = self.ymode[colcount]
            except:
                ymode  = self.ymode
            self.ax.set_yscale(ymode)
            self.ax.set_xscale(self.xmode)

            # --------------------------------------------------------------------------------------------- #
            # Iterate over the different curves to plot in the same subplot
            if hasattr(y[nSP][0], "__len__"):
                for nCurv in range( len(y[nSP]) ):
                    # Read the plot properties
                    try: color     = self.color[colcount]
                    except: color  = self.color
                    try: mksize    = self.mksize[colcount]
                    except: mksize = self.mksize
                    try: alpha     = self.alpha[colcount]
                    except: alpha  = self.alpha
                    try: ncol      = self.ncol[colcount]
                    except: ncol   = self.ncol
                    try: loc       = self.loc[colcount]
                    except: loc    = self.loc
                    try: legend    = self.label[colcount]
                    except: legend = self.label 
                    try: lstyle    = self.lstyle[colcount]
                    except: lstyle = self.lstyle
                    try: mktype    = self.mktype[colcount]
                    except : mktype= self.mktype

                    # Extract the errors if any are given
                    if (self.yerr is not None) and (hasattr(self.yerr[nSP][nCurv], "__len__")):
                        yerrnCurv = self.yerr[nSP][nCurv]
                    else:
                        yerrnCurv = None
                    if (self.xerr is not None) and (hasattr(self.xerr[nSP][nCurv], "__len__")):
                        xerrnCurv = self.xerr[nSP][nCurv] 
                    else:
                        xerrnCurv = None

                    # Plot limits as down-arraows
                    if (self.limit is not None) and (self.limit[nSP][nCurv]):
                        self.ax.errorbar(x[nSP][nCurv], y[nSP][nCurv], xerr=xerrnCurv, 
                            yerr=[yerrnCurv, np.zeros( len(yerrnCurv) )], fmt='none', 
                            ecolor=color, elinewidth=0.5, alpha=alpha, capsize=0, 
                            barsabove=False, lolims=False, uplims=False, xlolims=False, 
                            xuplims=False, errorevery=1, capthick=None, zorder=nCurv, legend=None)
                        self.ax.plot(x[nSP][nCurv], y[nSP][nCurv]-yerrnCurv, marker='v',
                            color=color, alpha=alpha, markersize=mksize, linestyle='',
                            markeredgecolor=color, zorder=nCurv)
                    # Fill an area between y[nSP][0][0] and y[nSP][0][1]
                    #elif hasattr(y[nSP][nCurv], "__len__"):
                    #    self.ax.fill_between(x[nSP][nCurv], y[nSP][nCurv][0], y[nSP][nCurv][1], facecolor=self.color, edgecolor='none', alpha=0.5,
                    #        rasterized=self.raster, zorder=-10)
                    # Plot a 'normal' curve
                    else:
                        if (legend is not None) and (legend != 'None') :
                            graph = self.ax.errorbar(x[nSP][nCurv], y[nSP][nCurv], yerr=yerrnCurv, 
                                xerr=xerrnCurv, fmt=mktype, ecolor=color, elinewidth=0.5, capsize=0,
                                linestyle=lstyle, markerfacecolor=color, markeredgecolor=color, 
                                color=color, markersize=mksize, label=legend, linewidth=self.lwdth, 
                                barsabove=False, errorevery=1, capthick=None, alpha=alpha, zorder=nCurv)
                            # Handling of the labels of the curves
                            handles, labels         = self.ax.get_legend_handles_labels()
                            handle_list, label_list = [], []
                            for k in xrange( len(labels) ):
                                if labels[k] in self.label:
                                    handle_list.append(handles[k])
                                    label_list.append(labels[k])
                            self.ax.legend(handle_list, label_list, loc="best", prop={'size':self.ftsize2},
                                frameon=True, numpoints=1, ncol=ncol, handletextpad=0.1)
                        else:
                            graph = self.ax.errorbar(x[nSP][nCurv], y[nSP][nCurv], yerr=yerrnCurv,
                                xerr=xerrnCurv, fmt=mktype, ecolor=color, elinewidth=0.5, capsize=0,
                                linestyle=lstyle, markerfacecolor=color, markeredgecolor=color, 
                                color=color, markersize=mksize, alpha=alpha, linewidth=self.lwdth,
                                barsabove=False, errorevery=1, capthick=None, zorder=nCurv)
                    colcount += 1
            # --------------------------------------------------------------------------------------------- #
            # There is only one curve per subplot
            else:
                # Read the plot properties
                try: color     = self.color[colcount]
                except: color  = self.color
                try: mksize    = self.mksize[colcount]
                except: mksize = self.mksize
                try: alpha     = self.alpha[colcount]
                except: alpha  = self.alpha
                try: ncol      = self.ncol[colcount]
                except: ncol   = self.ncol
                try: loc       = self.loc[colcount]
                except: loc    = self.loc
                try: legend    = self.label[colcount]
                except: legend = self.label 
                try: lstyle    = self.lstyle[colcount]
                except: lstyle = self.lstyle
                try: mktype    = self.mktype[colcount]
                except : mktype= self.mktype

                # Extract the errors if any are given
                if (self.yerr is not None) and (hasattr(self.yerr[nSP], "__len__")):
                    yerrSP = self.yerr[nSP]
                else:
                    yerrSP = None
                if (self.xerr is not None) and (hasattr(self.xerr[nSP], "__len__")):
                    xerrSP = self.xerr[nSP] 
                else:
                    xerrSP = None
                # Plot
                if (self.limit is not None) and (self.limit[nSP]):
                    self.ax.errorbar(x[nSP], y[nSP], xerr=xerrSP, 
                        yerr=[yerrSP, np.zeros( len(yerrSP) )], fmt='none', 
                        ecolor=color, elinewidth=0.5, alpha=alpha, capsize=0, 
                        barsabove=False, lolims=False, uplims=False, xlolims=False, 
                        xuplims=False, errorevery=1, capthick=None, legend=None)
                    self.ax.plot(x[nSP], y[nSP]-yerrSP, marker='v',
                        color=color, alpha=alpha, markersize=mksize, linestyle='',
                        markeredgecolor=color)
                else:
                    self.ax.errorbar(x[nSP], y[nSP], yerr=yerrSP, xerr=xerrSP, fmt=mktype, ecolor=color,
                        elinewidth=0.5, capsize=0, linestyle=lstyle, markerfacecolor=color, 
                        markeredgecolor=color, markersize=mksize, label=legend, alpha=alpha, color=color,
                        barsabove=False, errorevery=1, capthick=None)
                colcount += 1
                if legend is not None:
                    # Handling of the labels of the curves
                    self.ax.legend(loc="best", prop={'size':self.ftsize2}, frameon=True, numpoints=1,
                        ncol=ncol, handletextpad=0.1)
                    handles, labels = self.ax.get_legend_handles_labels()
                    handle_list, label_list = [], []
                    for k in xrange(len(labels)):
                        if labels[k] in self.label:
                            handle_list.append(handles[k])
                            label_list.append(labels[k])
                    self.ax.legend(handle_list, label_list, loc="best", prop={'size':self.ftsize2}, 
                        frameon=True, numpoints=1, ncol=ncol, handletextpad=0.1)

            # --------------------------------------------------------------------------------------------- #
            # Make pretty each subplot

            # Shift the x-label
            self.ax.yaxis.set_label_coords(self.labelx, 0.5)
            # Set the y-label for each subplot
            self.ax.set_ylabel(self.ylabel[nSP], fontsize=self.ftsize1, multialignment='center')
            self._plotDisplay()

            # Dimensions
            self.ax.set_xlim(xmin, xmax) # Every subplot has the same x-axis 
            ymin, ymax = self.ax.get_ylim()
            try: ymin  = self.ymin[nSP]
            except: pass
            try: ymax  = self.ymax[nSP]
            except: pass
            self.ax.set_ylim(ymin, ymax) 

            # Draw a horizontal line
            if (self.hline is not None) and (self.hline[nSP] is not None):
                # Multiple h-line to draw
                self.ax.axhline(y=self.hline[nSP], color='black', linestyle=':')
            # Fill an area
            if self.fill is not None:
                #self.ax.fill_between(x[nSP][nCurv], y[nSP][nCurv][0], y[nSP][nCurv][1], facecolor=self.color, edgecolor='none', alpha=0.5,
                #    rasterized=self.raster, zorder=-10)
                for k in range(len(self.fill)/2):
                    self.ax.axvspan(self.fill[k*2], self.fill[k*2+1], facecolor=self.shadecol, 
                        edgecolor="none", linewidth=0., zorder=-10, alpha=0.5)
            # For all upper subplot, remove the last ticks
            if nSP != len(self.prop)-1:
                plt.setp(self.ax.get_xticklabels(), visible=False)
                self.ax.set_xlabel('')
                ymincheck, ymaxcheck=self.ax.get_ylim()
                if ymaxcheck > ymincheck:
                    self.ax.get_yticklabels()[0].set_visible(False)
                else: # in case of a revert y axis...
                    self.ax.get_yticklabels()[-1].set_visible(False)

        self.f.subplots_adjust(hspace=0)
    def histo(self, data, bins=10, shaded=[], **kwargs):
        self._evalKwargs(kwargs)

        plt.hist(data, bins=bins, facecolor=self.color, edgecolor=self.color)
        if len(shaded)!=0:
            for i in range(len(shaded)/2) :
                self.ax.axvspan(shaded[i*2], shaded[i*2+1], facecolor=self.shadecol, edgecolor="none", linewidth=0., alpha=0.5)
                self.ax.axvline(x=shaded[i*2], color='black', linestyle=':') 
                self.ax.axvline(x=shaded[i*2+1], color='black', linestyle=':') 
        self._plotDisplay()
    def figtext(self, x, y, s='text', fontsize=9, halign='left'):
        """ Print some text at the relative location (x, y)
        """
        plt.figtext(x=x, y=y, s=s, fontsize=self.ftsize1, horizontalalignment=halign)
    def text(self, x, y, s='text', halign='left', valign='bottom', rot='horizontal'):
        """ Put some text
        """
        self.ax.text(x=x, y=y, s=s, fontsize=self.ftsize2, horizontalalignment=halign,
            verticalalignment=valign, rotation=rot)
    def save(self):
        """ Save the figure
        """
        self.f.savefig(self.savepath + self.figname, dpi=300, facecolor='none', edgecolor='none', 
            transparent=True, bbox_inches='tight')
        plt.close(self.f)


    # ============================================================================================= #
    #                                      Internal functions                                       #
    # ============================================================================================= #
    def _evalKwargs(self, kwargs):
        """ Evaluate all keyword arguments called inside a function
        """

        for key, value in kwargs.iteritems():
            if   key == 'xmin':      self.xmin      = value
            elif key == 'xmax':      self.xmax      = value
            elif key == 'ymin':      self.ymin      = value
            elif key == 'ymax':      self.ymax      = value
            elif key == 'raster':    self.raster    = value
            elif key == 'yerr':      self.yerr      = value
            elif key == 'xerr':      self.xerr      = value
            elif key == 'limit':     self.limit     = value # boolean
            elif key == 'shaded':    self.shaded    = value # list
            elif key == 'hline':     self.hline     = value
            elif key == 'label':     self.label     = value
            elif key == 'fillbelow': self.fillbelow = value # boolean
            elif key == 'alpha':     self.alpha     = value
            elif key == 'scinot':    self.scinot    = value
            elif key == 'prop':      self.prop      = value
            else:
                print("\t=== Key '{}' unknown ===".format(key))
    def _plotDisplay(self):
        """ Default plotting settings
        """

        # Axis parameters
        self.ax.set_xlabel(self.xlabel, fontsize=self.ftsize1)
        if isinstance(self.ylabel, basestring):
            self.ax.set_ylabel(self.ylabel, fontsize=self.ftsize1)
        plt.setp(self.ax.get_xticklabels(), rotation='horizontal', fontsize=self.ftsize2)
        plt.setp(self.ax.get_yticklabels(), rotation='horizontal', fontsize=self.ftsize2)
        self.ax.minorticks_on()
        # Scientific notation or not
        if self.scinot:
            def _scienNot_(x, y):
                if x == 0.:
                    return r'0.0'
                else:
                    return r'$%1.1f \! \times \! 10^{%1.0f}$'%( x/(10**np.floor(np.log10(x))), np.floor(np.log10(x)) )
            self.ax.yaxis.set_major_formatter(mtick.FuncFormatter( _scienNot_ ))
        # Plot limits
        if not hasattr(self.xmin, "__len__") and not hasattr(self.xmax, "__len__"):
            xmin, xmax = self.ax.get_xlim()
            if self.xmin is not None:
                xmin = self.xmin
            if self.xmax is not None:
                xmax = self.xmax
            self.ax.set_xlim(xmin, xmax)
        if not hasattr(self.ymin, "__len__") and not hasattr(self.ymax, "__len__"):
            ymin, ymax = self.ax.get_ylim()
            if self.ymin is not None:
                ymin = self.ymin
            if self.ymax is not None:
                ymax = self.ymax
            
            self.ax.set_ylim(ymin, ymax)


# ============================================================================================= #
# ============================================================================================= #
# ============================================================================================= #
class FermiMap(object):
    """ A class to plot fancy maps
    """
    def __init__(self, xsize=7.2, ysize=7.2):
        # Figure properties
        self.savepath  = None
        self.xsize     = xsize
        self.ysize     = ysize
        self.figname   = 'mymap.pdf'
        self.ftsize1   = 9
        self.ftsize2   = 6
        self.ftsize3   = 5
        self.color     = 'black'
        self.colormap  = 'Spectral_r'
        self.cbarlabel = r'cbarlabel'
        self.xlabel    = None
        self.ylabel    = None
        self.imwidth   = None
        self.imheight  = None
        self.datamin   = None
        self.datamax   = None
        self.csys      = None
        self.obsdate   = None

        # Plotting routine
        self.f = plt.figure()
        self.f.set_size_inches(self.xsize, self.ysize)
        self.gc       = None
        self.image    = None
        self.grid     = True
        self.contlist = None
        self.scl      = 1.
        self.scale    = 'linear' # 'log', 'sqrt', 'arcsinh', 'power'

    # ============================================================================================= #
    #                                     Attributes Checking                                       #
    # ============================================================================================= #
    @property
    def savepath(self):
        return self._savepath
    @savepath.setter
    def savepath(self, s):
        if s is None:
            self._savepath = os.path.abspath(os.getcwd())
        else:
            self._savepath = os.path.abspath(s)
            if not os.path.isdir(self._savepath):
                os.makedirs(self._savepath)

    @property
    def xsize(self):
        return self._xsize
    @xsize.setter
    def xsize(self, x):
        def cm2inch(size):
            """ Return the conversion in inches
            """
            return size / 2.54
        if x is None:
            self._xsize = cm2inch(8)
        else:
            self._xsize = cm2inch(x)

    @property
    def ysize(self):
        return self._ysize
    @ysize.setter
    def ysize(self, y):
        def cm2inch(size):
            """ Return the conversion in inches
            """
            return size / 2.54
        if y is None:
            self._ysize = cm2inch(8)
        else:
            self._ysize = cm2inch(y)

    @property
    def image(self):
        return self._image
    @image.setter
    def image(self, im):
        if im is None:
            self._image = None
        else:
            if os.path.isfile(im):
                self._image = im
                data   = fits.getdata(self._image, ext=0)
                header = fits.getheader(self._image, ext=0)
                self.obsdate = Time(header['DATE-OBS'], format="isot", scale="utc")
                if 'GLON' in header['CTYPE1']:
                    self.csys = 'GAL'
                elif 'RA' in header['CTYPE1']:
                    self.csys = 'CEL'
                else:
                    self.csys = 'CEL'
                    print("\t=== Coordinate in fits header '{}' not understood ===".format(header['CTYPE1']))
                self.imwidth  = header['NAXIS1']
                self.imheight = header['NAXIS2']
                self.datamin  = data.min()
                self.datamax  = data.max()
            else:
                print("\t=== Unable to find the image '{}' ===".format(im))
            

    # ============================================================================================= #
    #                                          Image plot                                           #
    # ============================================================================================= #
    def mapSky(self):
        """ Plot a nice map
        """
        import aplpy

        # Plot with aplpy
        self.gc = aplpy.FITSFigure(self.image, figure=self.f, 
            dimensions=[0,1], slices=[0,0], subplot=[0.1, 0.9, 0.9, 0.9])
        
        # Coordinate Grid
        if self.grid:
            self.gc.add_grid()
            self.gc.grid.set_color(self.color)
            self.gc.grid.set_alpha(0.3)
            self.gc.grid.set_linewidth(0.2)

        self._colorBar()
        self._plotDisplay()
    def circle(self, center, rad):
        """ Overplot a circle on the map 
        """
        self.gc.show_circles(center[0], center[1], rad, facecolor='none', edgecolor=self.color, linewidth=0.5)
    def srcSky_v0(self, ra, dec, names):
        """ Add sources names onto a map
        """
        import networkx as nx

        G = nx.DiGraph()
        data_nodes = []
        init_pos   = {}
        for rai, deci, name in zip(ra, dec, names):
            data_str = 'data_{0}'.format(name)
            G.add_node(data_str)
            G.add_node(name)
            G.add_edge(name, data_str)
            data_nodes.append(data_str)
            init_pos[data_str] = (rai, deci)
            init_pos[name]     = (rai, deci)

        pos = nx.spring_layout(G, pos=init_pos, fixed=data_nodes, k=0.0147) #0.126

        pos_after  = np.vstack([pos[d] for d in data_nodes])
        pos_before = np.vstack([init_pos[d] for d in data_nodes])
        scale, shift_x = np.polyfit(pos_after[:,0], pos_before[:,0], 1)
        scale, shift_y = np.polyfit(pos_after[:,1], pos_before[:,1], 1)
        shift = np.array([shift_x, shift_y])
        for key, val in pos.items():
            pos[key] = (val*scale) + shift

        for name, data_str in G.edges():
            self.gc.add_label(pos[name][0], pos[name][1], name,
                size=self.ftsize3, color=self.color)
            self.gc.show_lines([np.array([[pos[name][0], pos[data_str][0]], [pos[name][1], pos[data_str][1]]])],
                color=self.color, linewidth=0.2)
        return
    def srcSky_v1(self, ra, dec, names):
        """ Add sources names onto a map
        """
        srcs = Table([ra, dec, names], names=('ra', 'dec', 'name'))
        raSort  = np.argsort(ra)

        # Distribute labels left and right of the figure
        srcRight = srcs[raSort][:np.floor(len(srcs)/2)]
        decSort  = np.argsort(srcRight['dec'])
        pixDy = self.imheight / (len(srcRight)+1)
        pixY  = self.imheight
        pixX  = 0.77*self.imwidth
        for src in srcRight[list(reversed(decSort))]:
            pixY -= pixDy
            raLabel, decLabel = self.gc.pixel2world(pixX, pixY)
            if src['name'] == 'TARGET':
                src['name'] = raw_input("\n--> Name of the target: ")
            self.gc.add_label(raLabel, decLabel, src['name'],
                size=self.ftsize3, color=self.color, horizontalalignment='left')
            self.gc.show_lines([np.array([[raLabel, src['ra']], [decLabel, src['dec']]])],
                color=self.color, linewidth=0.2)

        srcLeft = srcs[raSort][np.floor(len(srcs)/2):]
        decSort = np.argsort(srcLeft['dec'])
        pixDy = self.imheight / (len(srcLeft)+1)
        pixY  = self.imheight
        pixX  = 0.23*self.imwidth
        for src in srcLeft[list(reversed(decSort))]:
            pixY -= pixDy
            raLabel, decLabel = self.gc.pixel2world(pixX, pixY)
            if src['name'] == 'TARGET':
                src['name'] = raw_input("\n--> Name of the target: ")
            self.gc.add_label(raLabel, decLabel, src['name'],
                size=self.ftsize3, color=self.color, horizontalalignment='right')
            self.gc.show_lines([np.array([[raLabel, src['ra']], [decLabel, src['dec']]])],
                color=self.color, linewidth=0.2)
        return
    def srcSky(self, ra, dec, names):
        """ Add sources names onto a map
        """
        srcs = Table([ra, dec, names], names=('ra', 'dec', 'name'))
        srcX, srcY = self.gc.world2pixel(srcs['ra'], srcs['dec'])

        # Sources on the left
        leftInd = srcX <= self.imwidth/2.
        srcLeft = srcs[leftInd]
        sortYin = np.argsort(srcY[leftInd])
        pixDy   = self.imheight / (len(srcLeft)+1)
        pixY    = self.imheight
        pixX    = 0.26*self.imwidth
        for src in srcLeft[list(reversed(sortYin))]:
            pixY -= pixDy
            raLabel, decLabel = self.gc.pixel2world(pixX, pixY)
            if src['name'] == 'TARGET':
                src['name'] = raw_input("\n--> Name of the target: ")
            self.gc.add_label(raLabel, decLabel, src['name'],
                size=self.ftsize3, color=self.color, horizontalalignment='right')
            self.gc.show_lines([np.array([[raLabel, src['ra']], [decLabel, src['dec']]])],
                color=self.color, linewidth=0.2)

        # Sources on the right
        rightInd = srcX > self.imwidth/2.
        srcRight = srcs[rightInd]
        sortYin  = np.argsort(srcY[rightInd])
        pixDy    = self.imheight / (len(srcRight)+1)
        pixY     = self.imheight
        pixX     = 0.74*self.imwidth
        for src in srcRight[list(reversed(sortYin))]:
            pixY -= pixDy
            raLabel, decLabel = self.gc.pixel2world(pixX, pixY)
            if src['name'] == 'TARGET':
                src['name'] = raw_input("\n--> Name of the target: ")
            self.gc.add_label(raLabel, decLabel, src['name'],
                size=self.ftsize3, color=self.color, horizontalalignment='left')
            self.gc.show_lines([np.array([[raLabel, src['ra']], [decLabel, src['dec']]])],
                color=self.color, linewidth=0.2)
        return
    def save(self):
        self.f.savefig(os.path.join(self.savepath, self.figname), dpi=200, facecolor='None', edgecolor='None', transparent=True, bbox_inches='tight')
        plt.close(self.f)
        #plt.close(self.gc)

    # ============================================================================================= #
    #                                      Internal functions                                       #
    # ============================================================================================= #
    def _plotDisplay(self):
        """ Default settings for map plotting
        """
        self.gc.tick_labels.set_xformat('ddd')
        self.gc.tick_labels.set_yformat('ddd')
        if self.csys == 'GAL':
            if self.xlabel is None: self.xlabel = r'Galactic longitude $l$ $(^{\circ})$'
            if self.ylabel is None: self.ylabel = r'Galactic latitude $b$ $(^{\circ})$'
        else:
            if self.xlabel is None: self.xlabel = r'RA (J2000)'
            if self.ylabel is None: self.ylabel = r'Dec (J2000)'
        self.gc.axis_labels.set_xtext(self.xlabel)
        self.gc.axis_labels.set_ytext(self.ylabel)
        self.gc.set_axis_labels_font(size=self.ftsize1)
        self.gc.tick_labels.set_font(size=self.ftsize2) # <====== perhaps a string here?
        self.gc.ticks.set_color('black')
    def _colorBar(self):
        """ Add a color bar to the figure
        """
        self.gc.show_colorscale(cmap=self.colormap, vmin=self.datamin, vmax=self.datamax, stretch=self.scale)        
        axcb = self.f.add_axes([0.1, 0.70, 0.9, 0.025])
        if self.scale   == 'linear': 
            normcb = mpl.colors.Normalize(vmin=self.datamin/self.scl, vmax=self.datamax/self.scl)
        elif self.scale == 'log':
            normcb = mpl.colors.LogNorm(vmin=self.datamin/self.scl, vmax=self.datamax/self.scl)
        else:
            print("\t=== Normalization not implemented for '{}' ===".format(self.scale))

        cbar = mpl.colorbar.ColorbarBase(axcb, cmap=self.colormap, norm=normcb, orientation='horizontal')
        cbar.solids.set_edgecolor('face')
        cbar.set_label(self.cbarlabel, fontsize=self.ftsize1, horizontalalignment='center')
        cbar.ax.tick_params(labelsize=self.ftsize2)







