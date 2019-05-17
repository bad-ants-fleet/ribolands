import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def main(args):
    outname = args.name if args.name else '{}.pdf'.format(args.input_filename)

    if args.files:
        #import numpy as np
        #rs = np.random.RandomState(365)
        #values = rs.randn(365, 4).cumsum(axis=0)
        #dates = pd.date_range("1 1 2016", periods=365, freq="D")
        #data = pd.DataFrame(values, dates, columns=["A", "B", "C", "D"])
        #data = data.rolling(7).mean()
        #print(data)

        li = []
        for fi in args.files:
            print(fi)
            df = pd.read_csv(fi, header=None, 
                names=['length', 
                    'number_structures', 'hidden_graph_size', 'number_edges', 
                    'time_algorithm', 'time_simulations', 'time_total', 
                    'deleted_nodes', 'still_reachables', 'rejected_deletions', 
                    'treekin_default', 'treekin_expo', 'treekin_plusI', 'treekin_fail', 
                    'fp_expansion', 'fp_coarse_grain', 'fp_prune', 'fp_total', 'dG_min'],
                comment='#', delim_whitespace=True, float_precision='2')
            df['sequence'] = fi
            li.append(df)
        df = pd.concat(li, axis=0, ignore_index=True)
        sns.lmplot(x="length", y="number_structures", hue="sequence", 
                data=df, fit_reg=False)
        #sns.pairplot(df, vars=["length", "number_structures", "time_algorithm"], 
        #        hue="sequence")
        #sns.relplot(x="length", y="number_edges", hue="sequence", data=df);


        plt.savefig(outname, bbox_inches='tight')
        return
    else :
        df = pd.read_csv(args.input_filename, header=None, 
            names=['length', 
                'number_structures', 'number_edges', 'hidden_graph_size', 'hidden_number_edges', 
                'time_algorithm', 'time_simulations', 'time_total', 
                'deleted_nodes', 'still_reachables',
                'treekin_default', 'treekin_expo', 'treekin_plusI', 'treekin_fail', 
                'fp_expansion', 'fp_connected',  'fp_coarse_grain', 'fp_prune', 'fp_total', 'dG_min'],
            comment='#', delim_whitespace=True, float_precision='2')

    mi = None
    ma = None


    #fig = plt.figure(figsize(10,5))
    #ax = fig.add_subplot(2, 2, 2)
    #fig.get_axes().set_xlim(20, 30)

    disp = ['number_structures']
    disp.append('number_edges')
    if args.graph_data:
        disp.append('hidden_graph_size')
        disp.append('hidden_number_edges')

    if args.time_data:
        disp.append('time_algorithm')
        disp.append('time_simulations')
        #disp.append('time_total')

    if args.prune_data:
        disp.append('deleted_nodes')
        disp.append('still_reachables')

    if args.treekin_data:
        disp.append('treekin_default')
        disp.append('treekin_expo')
        disp.append('treekin_plusI')
        disp.append('treekin_fail')

    if args.findpath_data:
        disp.append('fp_expansion')
        disp.append('fp_connected')
        disp.append('fp_coarse_grain')
        disp.append('fp_prune')
        disp.append('fp_total')

    if args.dG_min:
        disp.append('dG_min')

    if mi is not None and ma is not None:
        df.iloc[mi:ma, :].plot(x='length', y=disp, subplots=True)
    else:
        df.plot(x='length', y=disp, subplots=True, legend=True)

    lenx = 15
    leny = 4 * len(disp)

    fig = plt.gcf()
    fig.set_size_inches(lenx, leny)

    plt.savefig(outname, bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='%(prog)s [options] filename')

    parser.add_argument("-v", "--verbose", action='count', default=0,
            help="""Track process by writing verbose output to STDOUT during
            calculations. Use --logfile if you want to see *just* verbose
            information via STDOUT.""")

    parser.add_argument('input_filename', default=None, nargs='?', metavar='<str>',
            help="Path to the input file.")

    parser.add_argument('--files', default=None, nargs='+', metavar='<str>',
            help="Input files.")

    parser.add_argument('--name', action='store', help="Graph output name.")

    parser.add_argument('-g', '--graph-data', action='store_true',
            help="Graph data.")

    parser.add_argument('-t', '--time-data', action='store_true',
            help="Runtime behavior.")

    parser.add_argument('-p', '--prune-data', action='store_true',
            help="Prune data.")

    parser.add_argument('-f', '--findpath-data', action='store_true',
            help="Findpath data.")

    parser.add_argument('-k', '--treekin-data', action='store_true',
            help="Treekin data.")

    parser.add_argument('-d', '--dG_min', action='store_true',
            help="Number of structures.")

    args = parser.parse_args()

    main(args)
