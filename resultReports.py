"""
Read the reports in .csv file and stores them in an array to be used
by the main module to be comparing the results

"""

import csv;

def Enum(*sequential, **named):
    Enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), Enums)

#CTReport = Enum('epidural', 'intraparenchymal', 'intraventricular', 'subarachnoid', 'subdural', 'any');

class Reports:

    reportsDict = dict();

    def __init__(self, reportsFilename):
        with open(reportsFilename, newline='') as csvfile:
            reader = csv.DictReader(csvfile, lineterminator='\n');
            for row in reader:
                if row['Label'] == '1':
                    self.reportsDict[row['ID'].rpartition('_')[0]] = True;
                
    def getReports(self):
        return self.reportsDict;
