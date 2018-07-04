
# Before you have to run:
# python -m daetools.dae_plotter.plotter &

import argparse

from daetools.pyDAE import *
from daetools.pyDAE.data_reporters import *
from daetools_extended.daesimulation_extended import daeSimulationExtended


def read_data(args):
    json_data=open(args.input).read()
    data = json.loads(json_data)
    return data

def configure(args):
    cfg = daeGetConfig()
    cfg.SetBoolean('daetools.activity.printHeader', False)
    cfg.SetFloat('daetools.IDAS.MaxStep',args.MaxStep)
    cfg.SetFloat('daetools.IDAS.relativeTolerance',args.relative_tolerance)
    cfg.GetInteger('daetools.IDAS.MaxNumSteps',args.MaxNumSteps)
    return cfg


def save_reports(args):

    # Save the model report and the runtime model report
    xmlfile1 = '{0}.model.xml'.format(args.input, )
    xmlfile1_modified = '{0}.model2.xml'.format(args.input, )
    htmlfile1 = '{0}.model.hml'.format(args.input, )
    xslfile1 = "examples/dae-tools.xsl"

    simulation.m.SaveModelReport(xmlfile1)
    #convert_to_xhtml(xmlfile1, htmlfile1, xslfile1)
    inject_external(xmlfile1,xmlfile1_modified)

    xmlfile2 = '{0}.model-rt.xml'.format(args.input, )
    xmlfile2_modified = '{0}.model2-rt.xml'.format(args.input, )
    htmlfile2 = '{0}.model-rt.hml'.format(args.input, )
    xslfile2 = "examples/dae-tools-rt.xsl"
    inject_external(xmlfile2,xmlfile2_modified)

    simulation.m.SaveModelReport(xmlfile2)
    # convert_to_xhtml(xmlfile1, htmlfile2, xslfile2)


def inject_external(xmlfile, xmlfile_modified):

    origin_str = ('dae-tools.xsl', 'dae-tools.css')
    destination_str = ('https://s3-us-west-2.amazonaws.com/jaimenms/daetools/dae-tools.xsl',
                       'https://s3-us-west-2.amazonaws.com/jaimenms/daetools/dae-tools.css')

    try:
        with open(xmlfile) as f, open(xmlfile_modified, "w") as g:

            text = f.read()
            for str, str_modified in zip(origin_str, destination_str):
                text = text.replace(str, str_modified)
            g.write(text)
    except:
        print("File not generated")


def convert_to_xhtml(infile, outfile, xslfile):

    from lxml import etree
    xslt_doc = etree.parse(xslfile)
    xslt_transformer = etree.XSLT(xslt_doc)
    source_doc = etree.parse(infile)
    output_doc = xslt_transformer(source_doc)
    output_doc.write(outfile, pretty_print=True)



def get_reporter(args):

    if args.format == 'json':
        dr = daeJSONFileDataReporter()
    elif args.format == 'xml':
        dr = daeXMLFileDataReporter()
    elif args.format == 'mat':
        dr = daeMatlabMATFileDataReporter()
    elif args.format == 'xslx':
        dr = daeExcelFileDataReporter()
    elif args.format == 'vtk':
        dr = daeVTKDataReporter()
    elif args.format == 'csv':
        dr = daeCSVFileDataReporter()
    else:
        dr = None

    return dr


def get_name(args, data):

    if not args.name:
        simName = data['name']
    else:
        simName = args.name

    return simName


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Simulate a model according to DAETOOLS based on a json data. '
                                                 'It is necessary to have an openned dae_plotter thread before '
                                                 'executing. For that, please execute the following command:  '
                                                 'python -m daetools.dae_plotter.plotter &')
    parser.add_argument('input', help='Path of the json input file.')
    parser.add_argument('--format', default='gui', help='Format of output, where gui is actually '
                                                        'a graphical interface.',
                        choices=['gui', 'json', 'xml', 'mat','vtk','xlsx','csv'])
    parser.add_argument('--output', help='Path to the output file (not used if format is gui).')
    parser.add_argument('--name', help='Simulation name.')
    parser.add_argument('--reporting_interval', type=int, default= 3600, help='Reporting interval in seconds.')
    parser.add_argument('--time_horizon', type=int, default= 20*24*3600, help='Time horizon in seconds')
    parser.add_argument('--relative_tolerance', type=float, default= 1e-6, help='Relative tolerance for the integration '
                                                                                'method.')
    parser.add_argument('--MaxStep', type=int, default= 10., help='IDAS.MaxStep parameter.')
    parser.add_argument('--MaxNumSteps', type=int, default= 1000000, help='IDAS.MaxNumSteps parameter.')

    args = parser.parse_args()

    if args.format != 'gui' and not args.output:
        args.output = '{0}.output.{1}'.format(args.input,args.format)

    # Read data
    data = read_data(args)

    # Configure
    cfg = configure(args)

    # Reports
    dr = get_reporter(args)

    # Name
    simName = get_name(args, data)

    # Instantiate
    simulation = daeSimulationExtended(simName, data=data, set_reporting = True, reporting_interval = args.reporting_interval, time_horizon = args.time_horizon)

    # Gui Option
    if args.format == 'gui':

        qtApp = daeCreateQtApplication(sys.argv)
        simulator = daeSimulator(qtApp, simulation=simulation)
        simulator.exec_()

    else:

        dr.Connect(args.output, simName)

        solver = daeIDAS()
        solver.RelativeTolerance = args.relative_tolerance
        log = daePythonStdOutLog()

        # Initialize
        simulation.Initialize(solver, dr, log)
        print("Number of equations", simulation.NumberOfEquations)
        print("Number of variables", simulation.TotalNumberOfVariables)
        save_reports(args)

        # Solve at time = 0
        simulation.SolveInitial()

        # Run
        simulation.Run()

        # Clean up
        simulation.Finalize()