import os
from time import localtime, strftime
from daetools.pyDAE.data_reporters import daeDelegateDataReporter, daeJSONFileDataReporter, daePandasDataReporter

def setupDataReporters(simulation):

    datareporter = daeDelegateDataReporter()

    dr = daeJSONFileDataReporter()

    datareporter.AddDataReporter(dr)

    json_filename = simulation.report_filename
    simName = os.path.basename(json_filename) + '_' + strftime(" [%d.%m.%Y %H:%M:%S]", localtime())
    dr.Connect(json_filename, simName)

    print('%s (%s): %s' % (
    dr.__class__.__name__, dr.ConnectString, 'connected' if dr.IsConnected() else 'NOT connected'))

    return datareporter


def setupPandaDataReporters(simulation):

    datareporter = daeDelegateDataReporter()

    dr = daePandasDataReporter()

    datareporter.AddDataReporter(dr)


    return datareporter