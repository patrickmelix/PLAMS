from scm.plams import Settings




### ==== DFT SETTINGS ==== ###
def DFT():
    s = Settings()  
    s.input.ams.task             = 'GeometryOptimization'
    s.input.adf.basis.type       = 'TZ2P'
    s.input.adf.basis.core       = 'None'
    s.input.adf.xc.hybrid        = 'B3LYP'
    s.input.adf.xc.Dispersion    = 'GRIMME3 BJDAMP'
    s.input.adf.Relativity.Level = 'None'
    s.input.adf.NumericalQuality = 'Good'
    s.input.adf.Symmetry         = 'NOSYM'
    s.input.ams.UseSymmetry      = 'No'
    return s


def DFT_oxidized():
    s = Settings()
    s.input.adf.Unrestricted     = 'Yes'
    s.input.adf.SpinPolarization = '1.0'
    s.input.ams.System.Charge    = '1.0'
    return s


def DFT_reduced():
    s = Settings()
    s.input.adf.Unrestricted      = 'Yes'
    s.input.adf.SpinPolarization  = '1.0'
    s.input.ams.System.Charge     = '-1.0'
    return s


### ==== DFTB SETTINGS ==== ###
def DFTB():
    s = Settings()
    s.input.ams.task   = 'GeometryOptimization'
    s.input.DFTB
    s.input.DFTB.Model = "GFN1-xTB" 
    return s


def DFTB_oxidized():
    s = Settings()
    s.input.ams.System.Charge   = '1.0'
    return s


def DFTB_reduced():
    s = Settings()
    s.input.ams.System.Charge    = '-1.0'
    return s


### ==== FREQ AND SOLVATION SETTINGS ==== ###
def frequencies():
    s = Settings()
    s.input.ams.properties.NormalModes   = 'Yes'
    s.input.ams.Properties.PESPointCharacter     = 'No'
    s.input.ams.NormalModes.ReScanFreqRange      = '-1000 0'
    s.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = -20
    return s


def COSMO():
    s = Settings() 
    s.input.adf.Solvation.Solv = "Name=Dichloromethane"
    return s