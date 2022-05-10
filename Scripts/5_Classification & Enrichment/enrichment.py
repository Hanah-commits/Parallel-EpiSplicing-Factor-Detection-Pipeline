import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
import pylab
import operator

epi = pd.read_csv('0_Files/features_epi.csv', delimiter='\t')
nonepi = pd.read_csv('0_Files/features_nonepi.csv', delimiter='\t')


def adjust_pvalue(df):
    pval_cols = df.columns.tolist()[1:]  # skipping gene-id column
    new_cols = []
    col_names = []
    for col in pval_cols:

        # get indices of null values
        na_idx = df[df[col] == 1.0].index.tolist()

        # adjust non-null p values
        pvals = df[col].values.tolist()
        pvals = [x for x in pvals if x != 1]
        adj_pval = p_adjust_bh(pvals).tolist()

        # insert null at original indices
        for idx in na_idx:
            adj_pval.insert(idx, None)

        new_cols.append(adj_pval)
        col_names.append(col)

    # adjusted p values as new df
    df1 = pd.DataFrame(columns=col_names)
    for i in range(len(new_cols)):
        df1[col_names[i]] = new_cols[i]

    return df1


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]



def enrichment (filename):
    epi = pd.read_csv(filename, delimiter='\t')

    epi.drop(['gene', 'flanks', 'mean_dpsi_per_lsv_junction', 'H3K27me3', 'H3K9me3', 'H3K27ac', 'H3K4me3'], axis = 1, inplace=True)

    epi = adjust_pvalue(epi)

    rbp_pval = {}
    sfs = ['HNRNPL', 'HNRNPH1', 'HNRNPK', 'SFPQ', 'HNRNPA1', 'HNRNPA2B1', 'PTBP1', 'HNRNPF', 'HNRNPH2', 'HNRNPM', 'FUS',
          'YBX1', 'PCBP1', 'HNRNPC', 'HuR', 'TARDBP', 'HNRNPU', 'PCBP2', 'SRSF9', 'SRSF1', 'SRSF7', 'U2AF2', 'SRSF10',
          'SRSF2', 'RALY', 'MBNL1', 'RBM4', 'PABPN1', 'RBM3', 'TIA1', 'KHDRBS1', 'RBM28', 'PABPC1', 'RBM5', 'SART3',
          'SNRNP70', 'FXR2', 'FXR1', 'ESRP2', 'HNRNPA1L2', 'ZNF638', 'SNRPA', 'RBM8A', 'FMR1', 'DAZAP1', 'RBM42',
          'ZCRB1', 'KHDRBS3', 'QKI', 'KHDRBS2', 'ZC3H10', 'RBM24', 'RBFOX1', 'BRUNOL5', 'BRUNOL4', 'BRUNOL6']

    for sf in sfs:
        rbp_pval[sf] = epi[sf][epi[sf] < 0.01].count() / len(epi)

    rbp_pval = dict( sorted(rbp_pval.items(), key=operator.itemgetter(1),reverse=True))

    return [*rbp_pval][:10]


def qqplot_interpolate(type, sfs, name):

    color = {
        'epi': 'g',
        'nonepi': 'r',
        'Epi&NonEpi': 'b'
    }

    if type == 0:

        fig = plt.figure(figsize=(10, 7))
        # fig.text(0.5, 0.04, '', ha='center')
        # fig.text(0.04, 0.5, 'common Y', va='center', rotation='vertical')
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        columns = 3
        rows = math.ceil(len(sfs)/3)
        i = 1
        for sf in sfs:
            test2 = epi[sf].values.tolist()
            test1 = nonepi[sf].values.tolist()

            # Calculate quantiles
            test1.sort()
            quantile_levels1 = np.arange(len(test1), dtype=float) / len(test1)

            test2.sort()
            quantile_levels2 = np.arange(len(test2), dtype=float) / len(test2)

            # Use the smaller set of quantile levels to create the plot(epigenes)
            quantile_levels = quantile_levels2

            # We already have the set of quantiles for the smaller data set
            quantiles2 = test2

            # We find the set of quantiles for the larger data set using linear interpolation
            quantiles1 = np.interp(quantile_levels, quantile_levels1, test1)

            # Add a reference line
            maxval = max(test1[-1], test2[-1])
            minval = min(test1[0], test2[0])

            fig.add_subplot(rows, columns, i)
            # Plot the quantiles to create the qq plot
            pylab.plot(quantiles1, quantiles2, 'o', markerfacecolor="None", markeredgecolor=color[name])
            pylab.plot([minval, maxval], [minval, maxval], 'k-')
            pylab.title(sf)

            i += 1

        # plt.show()
        plt.savefig('../Output_Files/enrichedRBP_'+name+'.png')

    else:

        test2 = epi['mean_dpsi_per_lsv_junction'].values.tolist()
        test1 = nonepi['mean_dpsi_per_lsv_junction'].values.tolist()

        # Calculate quantiles
        test1.sort()
        quantile_levels1 = np.arange(len(test1), dtype=float) / len(test1)

        test2.sort()
        quantile_levels2 = np.arange(len(test2), dtype=float) / len(test2)

        # Use the smaller set of quantile levels to create the plot(epigenes)
        quantile_levels = quantile_levels2

        # We already have the set of quantiles for the smaller data set
        quantiles2 = test2

        # We find the set of quantiles for the larger data set using linear interpolation
        quantiles1 = np.interp(quantile_levels, quantile_levels1, test1)

        # Plot the quantiles to create the qq plot
        pylab.plot(quantiles1, quantiles2, 'o')

        # Add a reference line
        maxval = max(test1[-1], test2[-1])
        minval = min(test1[0], test2[0])
        pylab.plot([minval, maxval], [minval, maxval], 'k-')
        pylab.title('dPSI distribution')
        # pylab.show()




if __name__ == "__main__":

    types = ['epi', 'nonepi']
    enriched = {}
    for t in types:
        most_impt = enrichment('0_Files/pvals_rbp'+t+'.csv')
        enriched[t] = most_impt

    inter = list(set(enriched[types[0]]) & set(enriched[types[1]]))
    qqplot_interpolate(0, inter, name='Epi&NonEpi')

    for k,v in enriched.items():
        enriched = [x for x in v if x not in inter]
        qqplot_interpolate(0, enriched, name=k)
