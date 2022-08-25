#!/usr/bin/env amspython
from scm.plams import *

"""

This example uses the Property Prediction tool to estimate some properties for
a list of molecules given as SMILES strings.

The output is written to table.csv

"""

def main():
    smiles_list = ['CCO', 'CCOC', 'OCCCN', 'C', 'C1=CC=C(C=C1)COCC2=CC=CC=C2']

    f = open("table.csv", "w")

    # print header
    header = "SMILES"
    for k in PropertyPrediction.properties:
        # the vapor pressure is temperature-dependent, so let's ignore it in this example
        if k == 'vaporpressure': 
            continue
        unit = PropertyPrediction.units[k]
        header += f",{k} [{unit}]"
    print(header, file=f)

    for smiles in smiles_list:
        # run the property prediction 
        results = PropertyPrediction(smiles).results

        # print all properties in the same order as the header
        line = smiles
        for k in PropertyPrediction.properties:
            if k == 'vaporpressure':
                continue
            value = results.get(k, None) 
            # if the value is exactly 0 then it is unreliable
            if value:
                line += ",{:.3f}".format(value)
            else:
                line += ",N/A"
        print(line, file=f)

    f.close()

if __name__ == '__main__':
    main()

