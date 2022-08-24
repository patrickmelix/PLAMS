import os
try:
    from rdkit import Chem
    has_rdkit = True
except ImportError:
    has_rdkit = False

from scm.utils.runsubprocess import RunSubprocess
from ...mol.molecule import Molecule
from ...tools.kftools import KFFile
from ..molecule.rdkit import to_rdmol

__all__ = ['PropertyPrediction']

class PropertyPrediction:
    properties = [
        'boilingpoint',
        'criticalpressure',
        'criticaltemp',
        'criticalvol',
        'dielectricconstant',
        'entropygas',
        'flashpoint',
        'gidealgas',
        'hcombust',
        'hformstd',
        'hfusion',
        'hidealgas',
        'hsublimation',
        'meltingpoint',
        'molarvol',
        'vaporpressure',
        'parachor',
        'solubilityparam',
        'tpt',
        'vdwarea',
        'vdwvol',
    ]

    units = {
        'boilingpoint' : 'K',
        'criticalpressure' : 'bar',
        'criticaltemp' : 'K',
        'criticalvol' : 'L/mol',
        'dielectricconstant' : '',
        'entropygas' : 'J/(mol K)',
        'flashpoint' : 'K',
        'gidealgas' : 'kJ/mol',
        'hcombust' : 'kJ/mol',
        'hformstd' : 'kJ/mol',
        'hfusion' : 'kJ/mol',
        'hidealgas' : 'kJ/mol',
        'hsublimation' : 'kJ/mol',
        'meltingpoint' : 'K',
        'molarvol' : 'L/mol',
        # density?!
        'vaporpressure' : 'bar',
        'parachor' : '',
        'solubilityparam' : '(cal/cm^3)^1/2',
        'tpt' : 'K', #triple point temperature
        'vdwarea' : 'angstrom^2',
        'vdwvol' : 'angstrom^3',
    }
    def __init__(self, molecule, temperature=298.15, n_temperature=10):
        """
            Usage:

            .. code-block:: python

                p = PropertyPrediction('CCCCO', temperature=340)
                print(f"SMILES: {p.smiles}, Temperature: {p.temperature} for the vaporpressure")
                for key in p.properties:
                    if key in p.results:
                        print(f"{key} {p.results[key]:.2f} {p.units[key]}")
                    else:
                        print(f"{key} not calculated")

            molecule: str or Molecule
                Can be given as a SMILES string or a Molecule

            temperature: float or 2-tuple of float (min, max)
                Temperature or temperature range. A temperature range only affects the property 'vaporpressure'

            n_temperature: int
                If Temperature is a 2-tuple range, number of temperature to sample in range

            Returns: dictionary of with property names as keys. The Units are given on https://www.scm.com/doc/COSMO-RS/Property_Prediction.html

            .. code-block:: none

                'boilingpoint',
                'criticalpressure',
                'criticaltemp',
                'criticalvol',
                'dielectricconstant',
                'entropygas',
                'flashpoint',
                'gidealgas',
                'hcombust',
                'hformstd',
                'hfusion',
                'hidealgas',
                'hsublimation',
                'meltingpoint',
                'molarvol',
                'vaporpressure',
                'parachor',
                'solubilityparam',
                'tpt',
                'vdwarea',
                'vdwvol',

            
        """


        executable = os.path.join( os.path.expandvars("$AMSBIN") , "prop_prediction" )
        if not os.path.isfile(executable):
            raise RuntimeError("ERROR: cannot find prop_prediction ... has amsbashrc been executed?")

        if isinstance(molecule, str):
            smiles = molecule
        elif isinstance(molecule, Molecule):
            if not has_rdkit:
                raise RuntimeError("Cannot pass in a Molecule to prop_prediction when RDkit is not installed")
            rdmol = to_rdmol(molecule)
            smiles = Chem.MolToSmiles(rdmol)
        else:
            raise TypeError("molecule must be str (SMILES) or Molecule, but got {}".format(type(molecule)))

        self.smiles = smiles
        self.temperature = temperature

        results_file      =  "tmp_results18954.crskf"
        while os.path.exists(results_file):
            results_file = "0"+results_file
            if '0'*1000 in results_file:
                raise RuntimeError("Too many temporary 0000*crskf files exist. Remove them, and do not call prop_prediction massively in parallel")

        subprocess_string = " --smiles '" + smiles + "' "
        try:
            subprocess_string += f' --temperature {temperature[0]} --temperature {temperature[1]} --n {n_temperature}'
        except TypeError:
            subprocess_string += f' --temperature {temperature}'

        RunSubprocess( executable + subprocess_string + " -o " + results_file )
        crskf = KFFile( results_file )

        results = {}
        for k in self.properties:
            try:
                results[k] = crskf.read( "PROPPREDICTION" , k )
                #ret[k]['unit'] = crskf.read( "PROPPREDICTION" , k+'(Units)' ).strip().rstrip('\x00')
            except KeyError:
                pass

        try:
            os.remove(results_file)
        except OSError:
            pass

        self.results = results


