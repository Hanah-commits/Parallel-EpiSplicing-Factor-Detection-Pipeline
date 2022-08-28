import pandas as pd
import math
import sys
import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pylab
import operator
from scipy import stats


epi = pd.read_csv('0_Files/features_epi.csv', delimiter='\t')
nonepi = pd.read_csv('0_Files/features_nonepi.csv', delimiter='\t')

epi["label"] = "epi"
nonepi["label"] = "nonepi"

sfs = pd.read_csv('0_Files/impt_features.csv', delimiter='\t')['Unnamed: 0'].values.tolist()

# RBPs with no binding site in any flank
all_zero = []
for column in epi:  # iterates by-name
    if epi[column].isna().all() or (epi[column] == 0).all():
        all_zero.append(column)

for column in nonepi:  # iterates by-name
    if nonepi[column].isna().all() or (nonepi[column] == 0).all():
        all_zero.append(column)

sfs = [sf for sf in sfs if sf not in all_zero]

with open('paths.json') as f:
    d = json.load(f)

hms = d["Histone modifications"]


def adjust_pvalue(df):
    pval_cols = df.columns.tolist()
    new_cols = []
    col_names = []
    for col in pval_cols:

        # get indices of null values
        na_idx = df[df[col].isnull()].index.tolist()

        # adjust non-null p values
        pvals = df[col].values.tolist()
        pvals = [x for x in pvals if not math.isnan(x)]
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
    cols = ['gene', 'flanks', 'mean_dpsi_per_lsv_junction'] + hms
    epi.drop(cols, axis = 1, inplace=True)

    epi = adjust_pvalue(epi)

    rbp_pval = {}
    sfs = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'PTBP1', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']

    for sf in sfs:
        rbp_pval[sf] = epi[sf][epi[sf] < 0.01].count() / len(epi)

    rbp_pval = dict( sorted(rbp_pval.items(), key=operator.itemgetter(1),reverse=True))

    return [*rbp_pval][:10]


def qqplot_interpolate(type, sfs, name, output_dir):

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
        plt.savefig(output_dir+'enrichedRBP_'+name+'.png')

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


def significane_test(test):

    prob_dict = {}
    sfs = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'PTBP1', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']

    if test == 'binomial':
        epi_pvals = pd.read_csv('0_Files/pvals_rbpepi.csv', delimiter='\t')[sfs]
        epi_pvals = adjust_pvalue(epi_pvals)

        p =0.5 #probability of success in one trial
        n = len(epi) # number of trials
        for sf in sfs:
            k = epi_pvals[sf][epi_pvals[sf] < 0.01].count() #number of successes
            prob = stats.binom.cdf(k, n, p)
            prob_dict[sf] = prob

    elif test == 'welch':

        for sf in sfs:
            prob = stats.ttest_ind(epi[sf], nonepi[sf], equal_var = False)
            prob_dict[sf] = prob[1]
    
    elif test == 'mannwhitney':

        for sf in sfs:
            prob = stats.mannwhitneyu(epi[sf], nonepi[sf])
            prob_dict[sf] = prob[1]

    
    print(prob_dict)


if __name__ == "__main__":

    types = ['epi', 'nonepi']
    enriched = {}
    for t in types:
        most_impt = enrichment('0_Files/pvals_rbp'+t+'.csv')
        enriched[t] = most_impt

    inter = list(set(enriched[types[0]]) & set(enriched[types[1]]))
    qqplot_interpolate(0, inter, name='Epi&NonEpi', output_dir=sys.argv[1])

    for k,v in enriched.items():
        enriched = [x for x in v if x not in inter]
        qqplot_interpolate(0, enriched, name=k, output_dir=sys.argv[1])
