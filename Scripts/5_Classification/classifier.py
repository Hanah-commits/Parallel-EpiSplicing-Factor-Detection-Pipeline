import pandas as pd
import numpy as np
import sys
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve
from collections import defaultdict
from functools import reduce


def stratified_classifier(output_dir):

    features = pd.read_csv('0_Files/all_features.csv', delimiter='\t')


    sf = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'PTBP1', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']

    sf_data = features[sf]
    features['label'] = features['label'].map({'epigene': 1, 'non-epigene': 0}).astype(int)
    X, y = sf_data.values, features['label'].values

    clf = RandomForestClassifier(n_estimators=100, criterion='gini')
    kf = StratifiedKFold(n_splits=10)

    tprs = []
    aucs= []
    base_fpr = np.linspace(0, 1, 101)

    plt.figure(figsize=(5, 5))
    plt.axes().set_aspect('equal', 'datalim')
    gini_scores = []

    for i, (train, test) in enumerate(kf.split(X, y)):
        model = clf.fit(X[train], y[train])
        y_score = model.predict_proba(X[test])
        fpr, tpr, _ = roc_curve(y[test], y_score[:, 1])
        roc_auc = metrics.auc(fpr, tpr)
        aucs.append(roc_auc)

        plt.plot(fpr, tpr, 'b', alpha=0.15)
        tpr = np.interp(base_fpr, fpr, tpr)
        tpr[0] = 0.0
        tprs.append(tpr)
        gini_scores.append(dict(zip(sf,model.feature_importances_)))

    # obtain range of AUCs
    auc_range = np.percentile(aucs, (2.5, 97.5))
    ci = float("%.2f" % (auc_range[1] - auc_range[0]))/2
    mean_auc = "%.2f" % (auc_range[1] - ci)

    # mean gini_impurity across all folds
    def foo(r, d):
        for k in d:
            r[k].append(d[k])
    gini_scores = reduce(lambda r, d: foo(r, d) or r, gini_scores, defaultdict(list))
    gini_scores = pd.DataFrame(gini_scores)
    mean_gini = gini_scores.mean(axis=0).sort_values(ascending=False)
    print('Gini Scores: ')
    print(mean_gini)
    mean_gini.to_csv('0_Files/impt_features.csv', sep='\t')
    

    tprs = np.array(tprs)
    mean_tprs = tprs.mean(axis=0)
    std = tprs.std(axis=0)
    # aucs = np.array(aucs)
    # mean_auc = aucs.mean(axis=0)

    tprs_upper = np.minimum(mean_tprs + std, 1)
    tprs_lower = mean_tprs - std

    plt.plot(base_fpr, mean_tprs, 'b', label = 'Mean AUC = ' + str(mean_auc) + ' '+ r'$\pm$' + ' ' + str(ci))
    plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.3)
    plt.legend(loc = 'lower right')
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.title="Receiver Operating Characteristic"
    plt.savefig(output_dir +'ROC.png')


if __name__ == "__main__":
    stratified_classifier(output_dir=sys.argv[1])
