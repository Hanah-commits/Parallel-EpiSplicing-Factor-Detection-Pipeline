import pandas as pd
import math
import json
import numpy as np
from scipy import stats
from argparse import ArgumentParser

# Get the process name, use it in the output directory

p = ArgumentParser()
p.add_argument("output_dir")
p.add_argument("--process", "-p",
        help="The name of the process")

args = p.parse_args()
proc = args.process


epi = pd.read_csv(f'{proc}_0_Files/features_epi.csv', delimiter='\t')
nonepi = pd.read_csv(f'{proc}_0_Files/features_nonepi.csv', delimiter='\t')

epi["label"] = "epi"
nonepi["label"] = "nonepi"

sfs = pd.read_csv(f'{proc}_0_Files/impt_features.csv', delimiter='\t')['Unnamed: 0'].values.tolist()

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
    data = json.load(f)
d = data[proc]

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


def significane_test(test):

    prob_dict = {}
    # sfs = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']

    if test == 'binomial':
        epi_pvals = pd.read_csv(f'{proc}_0_Files/pvals_rbpepi.csv', delimiter='\t')[sfs]
        epi_pvals = adjust_pvalue(epi_pvals)

        p =0.5 #probability of success in one trial
        n = len(epi) # number of trials
        for sf in sfs:
            k = epi_pvals[sf][epi_pvals[sf] < 0.01].count() #number of successes
            prob = stats.binom.cdf(k, n, p)
            prob_dict[sf] = 1- prob

        for k,v in prob_dict.items():
            if (v <0.05):
                print(k) 

    elif test == 'welch':

        enriched_epi = []
        enriched_nonepi = []

        for sf in sfs:
            prob = stats.ttest_ind(epi[sf], nonepi[sf], equal_var = False, alternative='greater')
            prob_dict[sf] = prob[1]
            if prob[1] < 0.05 :
                enriched_epi.append(sf)

        for sf in sfs:
            prob = stats.ttest_ind(nonepi[sf], epi[sf], equal_var = False, alternative='greater')
            prob_dict[sf] = prob[1]
            if prob[1] < 0.05 :
                enriched_nonepi.append(sf)


    elif test == 'mannwhitney':

        enriched_epi = []
        enriched_nonepi = []

        for sf in sfs:
    
            prob = prob = stats.mannwhitneyu(epi[sf],nonepi[sf], alternative='greater')
            if prob[1] <= 0.05:
                enriched_epi.append(sf)

            prob = prob = stats.mannwhitneyu(nonepi[sf],epi[sf], alternative='greater')
            if prob[1] <= 0.05:
                enriched_nonepi.append(sf)


    elif test == 'mannwhitney2':

        for sf in sfs:
            prob = stats.mannwhitneyu(epi[sf], nonepi[sf])
            if prob[1] <= 0.05:
                prob_dict[sf] = prob[1]

        enriched_epi = []
        enriched_nonepi = []
        for k,v in prob_dict.items():
            nonzero_epi = epi[epi[k]!=0][k]
            nonzero_nonepi = nonepi[nonepi[k]!=0][k]

            
            prob = prob = stats.mannwhitneyu(nonzero_epi,nonzero_nonepi, alternative='greater')
            if prob[1] <= 0.05:
                enriched_epi.append(k)
            else:
                impt = False

            prob = prob = stats.mannwhitneyu(nonzero_nonepi,nonzero_epi, alternative='greater')
            if prob[1] <= 0.05:
                enriched_nonepi.append(k)
                impt = True

    elif test == "exact":

        enriched_epi = []
        enriched_nonepi = []

        for sf in sfs:

            row1 = []
            row2 = []

            row1.append((epi[sf] != 0).sum())
            row2.append((epi[sf] == 0).sum())

            row1.append((nonepi[sf] != 0).sum())
            row2.append((nonepi[sf] == 0).sum())


            table = [row1, row2]

            res = stats.barnard_exact(table, alternative='greater', pooled=False)
            if res.pvalue <= 0.05:
                enriched_epi.append(sf)
                

        for sf in sfs:
            row1 = []
            row2 = []

            row1.append((nonepi[sf] != 0).sum())
            row2.append((nonepi[sf] == 0).sum())

            row1.append((epi[sf] != 0).sum())
            row2.append((epi[sf] == 0).sum())


            table = [row1, row2]

            res = stats.barnard_exact(table, alternative='greater', pooled=False)
            if res.pvalue <= 0.05:
                enriched_nonepi.append(sf)

    elif test == "fisher":

        enriched_epi = []
        enriched_nonepi = []

        for sf in sfs:

            row1 = []
            row2 = []

            row1.append((epi[sf] != 0).sum())
            row2.append((epi[sf] == 0).sum())

            row1.append((nonepi[sf] != 0).sum())
            row2.append((nonepi[sf] == 0).sum())


            table = [row1, row2]

            res = stats.fisher_exact(table, alternative='greater')
            if res[1] <= 0.05:
                enriched_epi.append(sf)
                

        for sf in sfs:
            row1 = []
            row2 = []

            row1.append((nonepi[sf] != 0).sum())
            row2.append((nonepi[sf] == 0).sum())

            row1.append((epi[sf] != 0).sum())
            row2.append((epi[sf] == 0).sum())


            table = [row1, row2]

            res = stats.fisher_exact(table, alternative='greater')
            if res[1] <= 0.05:
                enriched_nonepi.append(sf)
               
    elif test == "gtest":
        enriched_epi = []
        enriched_nonepi = []

        for sf in sfs:

            row1 = []
            row2 = []

            row1.append((nonepi[sf] != 0).sum())
            row2.append((nonepi[sf] == 0).sum())

            row1.append((epi[sf] != 0).sum())
            row2.append((epi[sf] == 0).sum())


            table = [row1, row2]
            chi, pvalue, dof, ex = stats.chi2_contingency(table, lambda_="log-likelihood")
            
            if pvalue <= 0.05:
                enriched_nonepi.append(sf)


            row1 = []
            row2 = []

            row1.append((epi[sf] != 0).sum())
            row2.append((epi[sf] == 0).sum())

            row1.append((nonepi[sf] != 0).sum())
            row2.append((nonepi[sf] == 0).sum())


            table = [row1, row2]
            chi, pvalue, dof, ex = stats.chi2_contingency(table, lambda_="log-likelihood")
            
            if pvalue <= 0.05:
                enriched_epi.append(sf)


    elif test == "chi":

        enriched_epi = []
        enriched_nonepi = []

        for sf in sfs:

            row1 = []
            row2 = []

            row1.append((nonepi[sf] != 0).sum())
            row1.append((nonepi[sf] == 0).sum())

            row2.append((epi[sf] != 0).sum())
            row2.append((epi[sf] == 0).sum())


            table = [row1, row2]
            chi, pvalue, dof, ex = stats.chi2_contingency(table, correction=True)
            

            if pvalue <= 0.05:
                    enriched_nonepi.append(sf)
            

            row1 = []
            row2 = []

            row1.append((epi[sf] != 0).sum())
            row2.append((epi[sf] == 0).sum())

            row1.append((nonepi[sf] != 0).sum())
            row2.append((nonepi[sf] == 0).sum())


            table = [row1, row2]
            chi, pvalue, dof, ex = stats.chi2_contingency(table)
            
            if pvalue <= 0.05:
                enriched_epi.append(sf)

    print('RBPS enriched in epigene flanks: ' ,len(enriched_epi))
    print('RBPs enriched in non-epigene flanks ', len(enriched_nonepi))
    print('###############################')

    with open(f'{proc}_0_Files/enriched_epi.txt', 'w') as f:
        for rbp in enriched_epi:
            f.write(f"{rbp}\n")

    with open(f'{proc}_0_Files/enriched_nonepi.txt', 'w') as f:
        for rbp in enriched_nonepi:
            f.write(f"{rbp}\n")
            
    return enriched_epi, enriched_nonepi


if __name__ == "__main__":
    
    significane_test("welch")
