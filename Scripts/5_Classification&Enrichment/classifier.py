import pandas as pd
import numpy as np
import sys
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import resample
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay


def classifier(output_dir, method='LinearRegression', predictors='all', alpha=1.0, l1_ratio=0.5, folds=5):

    features = pd.read_csv('0_Files/all_features.csv', delimiter='\t')
    sf = ['HNRNPL', 'HNRNPH1', 'HNRNPK', 'SFPQ', 'HNRNPA1', 'HNRNPA2B1', 'PTBP1', 'HNRNPF', 'HNRNPH2', 'HNRNPM', 'FUS',
          'YBX1', 'PCBP1', 'HNRNPC', 'HuR', 'TARDBP', 'HNRNPU', 'PCBP2', 'SRSF9', 'SRSF1', 'SRSF7', 'U2AF2', 'SRSF10',
          'SRSF2', 'RALY', 'MBNL1', 'RBM4', 'PABPN1', 'RBM3', 'TIA1', 'KHDRBS1', 'RBM28', 'PABPC1', 'RBM5', 'SART3',
          'SNRNP70', 'FXR2', 'FXR1', 'ESRP2', 'HNRNPA1L2', 'ZNF638', 'SNRPA', 'RBM8A', 'FMR1', 'DAZAP1', 'RBM42',
          'ZCRB1', 'KHDRBS3', 'QKI', 'KHDRBS2', 'ZC3H10', 'RBM24', 'RBFOX1', 'BRUNOL5', 'BRUNOL4', 'BRUNOL6']

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


if __name__ == "__main__":
    classifier(predictors='sf', output_dir=sys.argv[1])
