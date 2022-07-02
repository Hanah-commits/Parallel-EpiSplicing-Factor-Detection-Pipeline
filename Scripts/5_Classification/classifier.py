import pandas as pd
import numpy as np
import sys
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import resample
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve
from collections import defaultdict
from functools import reduce


def classifier(output_dir, method='LinearRegression', predictors='all', alpha=1.0, l1_ratio=0.5, folds=5):

    features = pd.read_csv('0_Files/all_features.csv', delimiter='\t')
    sf = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'PTBP1', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']

    sf_data = features[sf]

    # check if imbalanced data  TODO: check if under or oversamping is needed
    counts = features['label'].value_counts().sort_values(ascending=False)
    classes = counts.index.tolist()
    counts = counts.tolist()
    imbalance = False
    majority = ''
    minority = ''
    if counts[0] > counts[1]:
        imbalance = True
        majority = classes[0]
        minority = classes[1]

    data = features.values

    X, y = [[]], [[]]
    feature_names = []

    if predictors == 'all':
        X, y = data[:, 1:], data[:, :1]
        feature_names = features.columns.tolist()[1:]
    elif predictors == 'peaks':
        X, y = data[:, 1:5], data[:, :1]
        feature_names = features.columns.tolist()[1:5]
    elif predictors == 'sf':
        # X, y = data[:, 5:], data[:, :1]
        # feature_names = features.columns.tolist()[5:]
        X, y = sf_data.values, data[:, :1]
        feature_names = sf

    # Split dataset into training set and test set
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3) # 70% training and 30% test

    if imbalance:
        train_data = np.column_stack((y_train, X_train))

        train_df = pd.DataFrame(data=train_data)

        # Separate majority and minority classes
        train_majority = train_df[train_df[0] == majority]
        train_minority = train_df[train_df[0] == minority]

        # Downsample majority class
        train_majority_downsampled = resample(train_majority,
                                           replace=False,  # sample without replacement
                                           n_samples=len(train_minority),  # to match minority class
                                           random_state=123)  # reproducible results

        # Combine minority class with downsampled majority class
        train_df = pd.concat([train_majority_downsampled, train_minority])

        train_data = train_df.values
        X_train, y_train = train_data[:, 1:], train_data[:, :1]


    # Create a Gaussian Classifier
    clf = RandomForestClassifier(n_estimators=100, criterion='gini')

    # Train the model using the training sets y_pred=clf.predict(X_test)
    clf.fit(X_train, y_train)


    ax = plt.gca()
    display = RocCurveDisplay.from_estimator(clf, X_test, y_test, ax=ax)
    display.plot()
    plt.savefig(output_dir +'ROC.png')


    feature_imp = pd.Series(clf.feature_importances_, index=feature_names).sort_values(ascending=False)
    feature_imp.tail(10).plot(kind='barh').set(xlabel="Feature Importance")
    # plt.savefig('../Output Files/Impt_features.png')
    feature_imp.to_csv('0_Files/impt_features.csv', sep='\t')

    y_pred = clf.predict(X_test)

    # Model Accuracy, how often is the classifier correct?
    print("Accuracy:", metrics.accuracy_score(y_test, y_pred))


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
    aucs = np.array(aucs)
    mean_auc = aucs.mean(axis=0)

    tprs_upper = np.minimum(mean_tprs + std, 1)
    tprs_lower = mean_tprs - std

    plt.plot(base_fpr, mean_tprs, 'b', label = 'AUC = %0.2f' % mean_auc)
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
    # classifier(predictors='sf', output_dir=sys.argv[1])
    stratified_classifier(output_dir=sys.argv[1])
