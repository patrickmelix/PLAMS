import os
try:
    from rdkit import Chem
    has_rdkit = True
except ImportError:
    has_rdkit = False

try:
    from scm.utils.runsubprocess import RunSubprocess
    has_scm = True
except ImportError:
    has_scm = False

from ...mol.molecule import Molecule
from ...tools.kftools import KFFile
from ..molecule.rdkit import to_rdmol
import numpy as np

__all__ = ['PropertyPrediction', 'PropertyPredictionList']

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
        #'vaporpressure',
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
    def __init__(self, molecule, temperatures=298.15, n_temperatures=10):
        """
            Usage:

            .. code-block:: python

                p = PropertyPrediction('CCCCO', temperatures=340)
                print(f"SMILES: {p.smiles}, Temperatures: {p.temperatures} for the vaporpressure")
                for key in p.properties:
                    if key in p.results:
                        print(f"{key} {p.results[key]:.2f} {p.units[key]}")
                    else:
                        print(f"{key} not calculated")
                print(p.vaporpressures)

            molecule: str or Molecule
                Can be given as a SMILES string or a Molecule

            temperatures: float or 2-tuple of float (min, max)
                Temperature or temperature range. A temperature only affects self.vaporpressures

            n_temperatures: int
                If Temperature is a 2-tuple range, number of temperature to sample in range


            The Units are given on https://www.scm.com/doc/COSMO-RS/Property_Prediction.html. Keys in the dictionary:

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
        if not has_scm or not os.path.isfile(executable):
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

        results_file      =  "tmp_results18954.crskf"
        while os.path.exists(results_file):
            results_file = "0"+results_file
            if '0'*1000 in results_file:
                raise RuntimeError("Too many temporary 0000*crskf files exist. Remove them, and do not call prop_prediction massively in parallel")

        subprocess_string = " --smiles '" + smiles + "' "
        try:
            subprocess_string += f' --temperature {temperatures[0]} --temperature {temperatures[1]} --n {n_temperatures}'
            self.temperatures = np.linspace(temperatures[0], temperatures[1], n_temperatures)
        except TypeError:
            subprocess_string += f' --temperature {temperatures}'
            self.temperatures = np.array([temperatures]).flatten()

        RunSubprocess( executable + subprocess_string + " -o " + results_file )
        crskf = KFFile( results_file )

        # 'vaporpressure' is not part of self.properties
        self.vaporpressures = crskf.read( "PROPPREDICTION", 'vaporpressure')
        self.vaporpressures = np.array([self.vaporpressures]).flatten()

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

class PropertyPredictionList:
    def __init__(self, smiles_list, temperatures=298.15, n_temperatures=10):
        self.predictions = {}
        self.smiles_list = smiles_list
        for smiles in smiles_list:
            self.predictions[smiles] = PropertyPrediction(smiles, temperatures=temperatures, n_temperatures=n_temperatures)

    def write_csv(self, filename):
        """
        Writes a csv file with properties
        """

        with open(filename, 'w') as f:
            f.write(self.get_csv_string())

    def get_csv_string(self):
        """
            Returns a string with one line per SMILES string and a header line.

            The vaporpressure is not given here, for the vaporpressure use get_csv_vaporpressure_string.
        """
        # print header
        header = "SMILES"
        for k in PropertyPrediction.properties:
            if k == 'vaporpressure': 
                continue
            unit = PropertyPrediction.units[k]
            header += f",{k} [{unit}]"

        lines = header + '\n'

        for smiles in self.smiles_list:
            # run the property prediction 
            results = self.predictions[smiles].results 

            # print all properties in the same order as the header
            line = smiles
            for k in PropertyPrediction.properties:
                if k == 'vaporpressure':
                    continue
                value = results.get(k, None) 
                # if the value is exactly 0 then it is unreliable
                if value:
                    line += ",{:.6f}".format(value)
                else:
                    line += ",N/A"
            lines += line + '\n'

        return lines

    def write_vaporpressures_csv(self, filename):
        """ Writes the vapor pressures to a .csv file with one line per SMILES and one column per temperature """
        with open(filename, 'w') as f:
            f.write(self.get_vaporpressures_csv_string())

    def get_vaporpressures_csv_string(self):
        header = ['SMILES/T (K). Vapor pressures in {}.'.format(PropertyPrediction.units['vaporpressure'])]
        header += self.predictions[self.smiles_list[0]].temperatures.tolist()
        header = ','.join(str(x) for x in header)
        lines = header + '\n'
        for smiles in self.smiles_list:
            lines += smiles + ','
            lines += ','.join(str(x) for x in self.predictions[smiles].vaporpressures) + '\n'
        return lines



