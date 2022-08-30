#!/usr/bin/env amspython
from scm.plams import *

"""

This example uses the Property Prediction tool to estimate some properties for
a molecule represented by a SMILES string.

You can also use the PropertyPredictionList class to run predictions on a list of
SMILES strings and print the results to csv files.

The vapor pressures are handled separately from other properties vapor pressures are temperature-dependent.

The output is written to vaporpressures.csv and properties.csv

Run this script using
$AMSBIN/amspython PropertyPrediction.py

"""

def main():
    single_compound_example()
    multi_compound_example()

def single_compound_example():
    print("Single compound example for SMILES CCO")
    p = PropertyPrediction('CCO')
    print("Boiling point: {} {}".format(p.results['boilingpoint'], p.units['boilingpoint']))
    print("Available properties: {}".format(p.properties))
    # vapor pressures are handled separately
    print("Temperatures ({}): {}".format('K', p.temperatures))
    print("Vapor pressures ({}): {}".format(p.units['vaporpressure'], p.vaporpressures))

def multi_compound_example():
    smiles_list = [
        'CCO', 
        'CCOC', 
        'OCCCN', 
        'C', 
        'C1=CC=C(C=C1)COCC2=CC=CC=C2'
    ]
    pp = PropertyPredictionList(smiles_list, temperatures=(280,340), n_temperatures=8)
    print("Writing info about SMILES {} to properties.csv and vaporpressures.csv".format(smiles_list))
    pp.write_csv("properties.csv")
    pp.write_vaporpressures_csv("vaporpressures.csv")


if __name__ == '__main__':
    main()

