#!/usr/bin/env python

__author__ = "Cyril Tasse, and Alan Loh"
#__copyright__ = "Copyright 2007, The Cogent Project"
__credits__ = ["Cyril Tasse", "Alan Loh"]
#__license__ = "GPL"
__version__ = "1.0.1"
#__maintainer__ = "Rob Knight"
#__email__ = "rob@spot.colorado.edu"
#__status__ = "Production"

# =========================================================================
#                                TO DESCRIPTION
#   Blablabla
#   
#   Example:
#   python ms2dynspec.py --msfile /usr/data/ --data CORRECTED --model PREDICT --sols killms.npz --srclist SRCPOS.txt --rad 10
# =========================================================================

# =========================================================================
#                                TO DO
# - convertSrclist: only keep (RA, Dec) wich are within the field -- DONE
# - convertSrclist: add some other ~random positions on wich to compute dynamic spectra for comparison
# =========================================================================

import sys
import os
import argparse
from pyrap.tables import table
from astropy.time import Time
from astropy import units as uni
from astropy.io import fits
from astropy import coordinates as coord
from astropy import constants as const
import numpy as np
import glob, os
import pylab


# =========================================================================
# =========================================================================
def angSep(ra1, dec1, ra2, dec2):
    """ Find the angular separation of two sources (ra# dec# in deg) in deg
        (Stolen from the LOFAR scripts), works --> compared with astropy (A. Loh)
    """
    b = (np.pi / 2) - np.radians(dec1)
    c = (np.pi / 2) - np.radians(dec2)
    temp = (np.cos(b) * np.cos(c)) + (np.sin(b) * np.sin(c) * np.cos(np.radians(ra1 - ra2)))
    if abs(temp) > 1.0:
        temp = 1.0 * cmp(temp, 0)
    return np.degrees(np.arccos(temp))

def convertSrclist(srclist, ms, fieldradius=10):
    """ Check that the source list provided is in the correct format
        return a dictionnary
    """
    tField  = table("%s::FIELD"%ms)
    ra0, dec0 = np.degrees(tField.getcol("PHASE_DIR").ravel()) # field center
    srcdict = {}
    for line in open(srclist):
        words = line.split()
        srcname = ' '.join(words[:-2])
        srcra   = float(words[-1])
        srcdec  = float(words[-2])
        if angSep(ra0, dec0, srcra, srcdec) <= fieldradius:
            # Only keep if the source is within 'fieldradius' from the image center
            srcdict[ srcname ] = {'RA': srcra, 'Dec': srcdec}
    tField.close()
    return srcdict

def checkArgs(args):
    """ Check that the arguments are correct
    """
    if (not os.path.isfile(args.msfile)) and (not os.path.exists(args.msfile)): # TO DO: if it is a directory where ae stored MS files...
        raise Exception("\t=== MS file/directory '{}'' not found ===".format(args.msfile))
    else:
        if args.msfile[-2:].upper() == 'MS':
            t = table(args.msfile, ack=False)
            args.msfile = [args.msfile]
        else:
            args.msfile = sorted(glob.glob("%s*"%args.msfile))
            t = table(args.msfile[0], ack=False)

    if args.data not in t.colnames():
        raise Exception("\t=== DATA column '{}' not found ===".format(args.data))

    if args.model not in t.colnames():
        raise Exception("\t=== MODEL column '{}' not found ===".format(args.model))

    if not os.path.isfile(args.sols):
        raise Exception("\t=== Solution file '{}'' not found ===".format(args.sols))
    else:
        if args.sols[-3:] != 'npz':
            raise Exception("\t=== Solution '{}' should be a .npz file ===".format(args.sols))

    if args.srclist is not None:
        if not os.path.isfile(args.srclist):
             raise Exception("\t=== Source list file '{}'' not found ===".format(args.srclist))
    else:
        args.srclist = convertSrclist(srclist=args.srclist, ms=args.msfile[0], fieldradius=args.rad)

    t.close()
    return # Everything's ok

# =========================================================================
# =========================================================================
class DynSpecMS():
    def __init__(self, ListMSName, ColName="DATA", ColName="PREDICT_KMS", UVRange=[1.,1000.], Sols='killMS.ALL_SOLS_1km_KAFCA.p.sols.npz'):
        self.ListMSName = ListMSName
        self.ColName    = ColName
        self.ModelName  = ModelName
        self.OutName    = self.ListMSName[0].split("/")[-1].split("_")[0]
        self.UVRange    = UVRange
        self.Sols       = Sols
        self.ReadMSInfos()

    def ReadMSInfos(self):
        print "Reading data"
        DicoMSInfos = {}
        t0  = table(self.ListMSName[0], ack=False)
        tf0 = table("%s::SPECTRAL_WINDOW"%self.ListMSName[0], ack=False)
        self.ChanWidth = tf0.getcol("CHAN_WIDTH").ravel()[0]
        tf0.close()
        self.times = np.unique(t0.getcol("TIME"))
        t0.close()
        
        for iMS, MSName in enumerate(sorted(self.ListMSName)):
            print MSName
            try:
                t = table(MSName, ack=False)
            except:
                DicoMSInfos[iMS] = {"Readable": False}
                print "problem reading %s"%MSName
                continue

            # Extract Jones matrices that will be appliedto the visibilities
            try:
                JonesSols = np.load("%s/%s"%(MSName, self.Sols))
                JonesSols = self.NormJones(JonesSols) # Normalize Jones matrices
            except:
                print "Problem reading solutions %s/%s"%(MSName, self.Sols)
                JonesSols = None

            tf = table("%s::SPECTRAL_WINDOW"%MSName, ack=False)
            ThisTimes = np.unique(t.getcol("TIME"))
            if not np.allclose(ThisTimes, self.times):
                raise ValueError("should have the same times")
            DicoMSInfos[iMS] = {"MSName": MSName,
                            "ChanFreq":   tf.getcol("CHAN_FREQ").ravel(),  # Hz
                            "ChanWidth":  tf.getcol("CHAN_WIDTH").ravel(), # Hz
                            "times":      ThisTimes,
                            "startTime":  Time(ThisTimes[0]/(24*3600.), format='mjd', scale='utc').isot,
                            "stopTime":   Time(ThisTimes[-1]/(24*3600.), format='mjd', scale='utc').isot,
                            "deltaTime":  (ThisTimes[-1] - ThisTimes[0])/3600., # h
                            "Readable":   True,
                            "JonesSols":  JonesSols}
            if DicoMSInfos[iMS]["ChanWidth"][0] != self.ChanWidth:
                raise ValueError("should have the same chan width")

        t.close()
        tf.close()
        self.nMS         = len(DicoMSInfos)
        self.DicoMSInfos = DicoMSInfos
        self.FreqsAll    = [DicoMSInfos[iMS]["ChanFreq"] for iMS in DicoMSInfos.keys() if DicoMSInfos[iMS]["Readable"]]
        self.Freq_minmax = np.min(self.FreqsAll), np.max(self.FreqsAll)
        self.NTimes      = self.times.size
        f0, f1           = self.Freq_minmax
        self.NChan       = int((f1 - f0)/self.ChanWidth) + 1
        self.IGrid       = np.zeros((self.NChan, self.NTimes, 4), np.complex128)

        return
        
        
    def StackAll(self, ra=None, dec=None):   
        print "Stacking data"
        
        f0, _ = self.Freq_minmax
        self.ra  = ra
        self.dec = dec

        self.WeightsGrid = np.zeros_like(self.IGrid)
        DicoMSInfos      = self.DicoMSInfos
        uv0, uv1         = np.array(self.UVRange) * 1000

        for iMS, MSName in enumerate(sorted(self.ListMSName)):
            print "%i/%i"%(iMS, self.nMS)
            if not DicoMSInfos[iMS]["Readable"]: continue

            t      = table(MSName, ack=False)
            tField = table("%s::FIELD"%MSName)
            self.ra0, self.dec0 = tField.getcol("PHASE_DIR").ravel() # radians!
            phaseCent = coord.SkyCoord(self.ra0, self.dec0, unit=uni.rad)
            if (self.ra is not None) and (self.dec is not None):
                phaseCent = coord.SkyCoord(self.ra, self.dec, unit=uni.deg) # not phase center but still...
            data   = t.getcol(self.ColName) - t.getcol(self.ModelName) # Residuals
            #data   = t.getcol(self.ModelName) # Model
            #data   = t.getcol(self.ColName) # data

            flag   = t.getcol("FLAG")
            times  = t.getcol("TIME")
            A0, A1 = t.getcol("ANTENNA1"), t.getcol("ANTENNA2")

            u, v, w = t.getcol("UVW").T
            d = np.sqrt(u**2 + v**2 + w**2)
            indUV = np.where( (d<uv0)|(d>uv1) )[0]
            flag[indUV, :, :] = 1 # flag according to UV selection
            data[flag] = 0 # put down to zeros flagged visibilities

            ich0 = int( (DicoMSInfos[iMS]["ChanFreq"][0] - f0)/self.ChanWidth )
            nch  = DicoMSInfos[iMS]["ChanFreq"].size

            # ------------------------------------------------------------------------------------------------
            #                               Pixel direction to sum visibilities
            # ------------------------------------------------------------------------------------------------
            if ra is None and dec is None:
                # considering central pixel for the dynamic spectrum computation
                ra  = np.degrees(self.ra0)
                dec = np.degrees(self.dec0)
                self.ra  = ra
                self.dec = dec
                
            # Considering another position than the phase center
            u0 = u.reshape( (-1, 1, 1) )
            v0 = v.reshape( (-1, 1, 1) )
            w0 = w.reshape( (-1, 1, 1) )
            f  = DicoMSInfos[iMS]["ChanFreq"].reshape( (1, nch, 1) )
            l, m = self.radec2lm(np.radians(ra), np.radians(dec))
            n = np.sqrt(1. - l**2. - m**2.)
            k = np.exp( -2.*np.pi*1j* f/const.c.value *(u0*l + v0*m + w0*(n-1)) ) # Phasing term
            # ------------------------------------------------------------------------------------------------

            for iTime in range(self.NTimes):
                indRow = np.where(times==self.times[iTime])[0]
                f   = flag[indRow, :, :]
                d   = data[indRow, :, :]
                A0s = A0[indRow]
                A1s = A1[indRow]
                kk  = k[indRow, :, :]

                # ------------------------------------------------------------------------------------------------
                #                                          Apply Jones matrices
                # ------------------------------------------------------------------------------------------------
                if DicoMSInfos[iMS]["JonesSols"] is not None:
                    self.J = DicoMSInfos[iMS]["JonesSols"]['G']
                    t0 = DicoMSInfos[iMS]["JonesSols"]['t0']
                    t1 = DicoMSInfos[iMS]["JonesSols"]['t1']
                    # Time slot for the solution
                    iTJones = np.where( t0 <= self.times[iTime] )[0][-1]
                    # Facet used
                    posFacets = coord.SkyCoord(DicoMSInfos[iMS]["JonesSols"]['ra'], DicoMSInfos[iMS]["JonesSols"]['dec'], unit=uni.rad)
                    iDJones = phaseCent.separation(posFacets).deg.argmin()                
                    # construct corrected visibilities
                    J0 = self.J[iTJones, 0, A0s, iDJones, 0, 0]
                    J1 = self.J[iTJones, 0, A1s, iDJones, 0, 0]
                    J0 = J0.reshape((-1, 1, 1))*np.ones((1, nch, 1))
                    J1 = J1.reshape((-1, 1, 1))*np.ones((1, nch, 1))
                    dcorr = J0.conj() * d * J1
                # ------------------------------------------------------------------------------------------------

                #ds=np.sum(d*kk, axis=0) # without Jones
                ds = np.sum(dcorr * kk, axis=0) # with Jones
                ws = np.sum(1-f, axis=0)
                self.IGrid[ich0:ich0+nch, iTime, :]       = ds
                self.WeightsGrid[ich0:ich0+nch, iTime, :] = ws

            # Fill properties
            self.tStart = DicoMSInfos[iMS]["startTime"]
            self.tStop  = DicoMSInfos[iMS]["stopTime"] 
            self.fMin   = self.Freq_minmax[0]
            self.fMax   = self.Freq_minmax[1]

        W = self.WeightsGrid
        G = self.IGrid
        W[W == 0] = 1
        Gn = G/W 

        GOut=np.zeros_like(self.IGrid)
        GOut[:, :, 0] = Gn[:, :, 0] + Gn[:, :, 3]
        GOut[:, :, 1] = Gn[:, :, 0] - Gn[:, :, 3]
        GOut[:, :, 2] = Gn[:, :, 1] + Gn[:, :, 2]
        GOut[:, :, 3] = 1j*(Gn[:, :, 1] - Gn[:, :, 2])
        self.GOut = GOut

    def NormJones(self, jMatrices):
        """ Normalisation of Jones Matrices
        """
        t0   = jMatrices['Sols']['t0']
        t1   = jMatrices['Sols']['t1']
        sols = jMatrices['Sols']['G']
        ra   = jMatrices['ClusterCat']['ra']
        dec  = jMatrices['ClusterCat']['dec']
        for iTime in xrange(sols.shape[0]):
            for iFreq in xrange(sols.shape[1]):
                for iAnt in xrange(sols.shape[2]):
                    for iDir in xrange(sols.shape[3]):
                        # 2D matrix [[xx, xy], [yx, yy]]
                        sols[iTime, iFreq, iAnt, iDir] /= np.linalg.norm(sols[iTime, iFreq, iAnt, iDir])
        return {'G':sols, 't0':t0, 't1':t1, 'ra':ra, 'dec':dec}     
    def WriteFits(self):
        """ Store the dynamic spectrum in a FITS file
        """
        fitsname = "DynSpec_{}_ra{}_dec{}.fits".format(self.OutName, self.ra, self.dec)

        # Create the fits file
        prihdr  = fits.Header() 
        prihdr.set('DATE-CRE', Time.now().iso.split()[0], 'Date of file generation')
        prihdr.set('CHAN-WID', self.ChanWidth, 'Frequency channel width')
        prihdr.set('FRQ-MIN', self.fMin, 'Minimal frequency')
        prihdr.set('FRQ-MAX', self.fMax, 'Maximal frequency')
        prihdr.set('OBS-STAR', self.tStart, 'Observation start date')
        prihdr.set('OBS-STOP', self.tStop, 'Observation end date')
        hdus    = fits.PrimaryHDU(header=prihdr) # hdu table that will be filled
        hdus.writeto(fitsname, clobber=True)

        label = ["I", "Q", "U", "V"]
        hdus = fits.open(fitsname)
        for ipol in range(4):
            Gn = self.GOut[:, :, ipol].real
            hdr   = fits.Header()
            hdr.set('STOKES', label[ipol], '')
            hdus.append( fits.ImageHDU(data=Gn, header=hdr, name=label[ipol]) )
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(fitsname, clobber=True)
        hdulist.close()
        print("\t=== Dynamic spectrum '{}' written ===".format(fitsname))
        return

    def PlotSpec(self):
        fig = pylab.figure(1, figsize=(15,8))
        label = ["I", "Q", "U", "V"]
        pylab.clf()
        for ipol in range(4):
            Gn = self.GOut[:, :, ipol].T.real
            v0, v1 = np.percentile(Gn.ravel(), [10., 90.])
            pylab.subplot(2, 2, ipol+1)
            pylab.imshow(Gn, interpolation="nearest", aspect="auto", vmin=3*v0, vmax=10*v1)
            pylab.title(label[ipol])
            pylab.colorbar()
            pylab.ylabel("Time bin")
            pylab.xlabel("Freq bin")
        pylab.tight_layout()
        pylab.draw()
        pylab.show(False)
        fig.savefig("DynSpec_{}_ra{}_dec{}.png".format(self.OutName, self.ra, self.dec))
        return Gn


    def radec2lm(self, ra, dec):
        # ra and dec must be in radians
        l = np.cos(dec) * np.sin(ra - self.ra0)
        m = np.sin(dec) * np.cos(self.dec0) - np.cos(dec) * np.sin(self.dec0) * np.cos(ra - self.ra0)
        return l, m
# =========================================================================
# =========================================================================


# =========================================================================
# =========================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--msfile", type=str, help="Name of MS file / directory", required=True)
    parser.add_argument("--data", type=str, default='CORRECTED', help="Name of DATA column", required=True)
    parser.add_argument("--model", type=str, help="Name of MODEL column", required=True)
    parser.add_argument("--sols", type=str, help="Jones solutions", required=True)
    parser.add_argument("--srclist", type=str, default=None, help="List of targets", required=True)
    parser.add_argument("--rad", type=float, default=10, help="Radius of the field", required=False)
    args = parser.parse_args()    
    checkArgs(args) # Verify that everything is correct before launching the dynamic spectrum computation

    D = DynSpecMS(ListMSName=args.msfile, ColName=args.data, ColName=args.model, Sols=args.sols)
    for key in args.srclist.keys():
        D.StackAll(ra=args.srclist[key]['RA'], dec=args.srclist[key]['Dec'])
        D.WriteFits()
        #D.PlotSpec 
# =========================================================================
# =========================================================================
