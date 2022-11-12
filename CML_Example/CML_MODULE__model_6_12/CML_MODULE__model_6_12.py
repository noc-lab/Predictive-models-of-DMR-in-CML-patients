def model_6_12():
    
    import numpy as np
    import pandas as pd
    import csv
    from sklearn import preprocessing
    from sklearn.metrics import roc_curve, auc, f1_score, precision_score, recall_score, precision_recall_curve, plot_precision_recall_curve 
    from sklearn.ensemble import RandomForestClassifier                                             
    from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score, KFold, StratifiedKFold 
    from sklearn.linear_model import LogisticRegression
    from sklearn.svm import LinearSVC
    from sklearn.neural_network import MLPClassifier 
    import lightgbm as lgb
    from sklearn.feature_selection import RFE,RFECV
    import scipy.io
    from scipy import stats
    import sklearn
    import math
    import pickle
    import warnings
    warnings.filterwarnings("ignore")
    sklearn.__version__
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)


    ## input
    Mol_3m_IM = None
    while Mol_3m_IM is None:
        try:
            Mol_3m_IM = float(input("Enter BCR-ABL1 ratio (percantage) at 3 months after start of Imatinib: "))
        except ValueError:
            print('The input value for BCR-ABL1 ratio (percantage) at 3 months is invalid.')
            Mol_3m_IM = float(input("Enter BCR-ABL1 ratio (percantage) at 3 months after start of Imatinib: "))

    Mol_6m_IM = None
    while Mol_6m_IM is None:
        try:
            Mol_6m_IM = float(input("Enter BCR-ABL1 ratio (percantage) at 6 months after start of Imatinib: "))
        except ValueError:
            print('The input value for BCR-ABL1 ratio (percantage) at 6 months is invalid.')
            Mol_6m_IM = float(input("Enter BCR-ABL1 ratio (percantage) at 6 months after start of Imatinib: "))

    WBC_count_at_diagnosis = None
    while WBC_count_at_diagnosis is None:
        try:
            WBC_count_at_diagnosis = float(input("Enter White Blood Cell count at diagnosis: "))
        except ValueError:
            print('The input value for White Blood Cell count at diagnosis is invalid.')
            WBC_count_at_diagnosis = float(input("Enter White Blood Cell count at diagnosis: "))    


    ## dataframe        
    variable_dict = {

        'Mol_3m_IM' : Mol_3m_IM,
        'Mol_6m_IM': Mol_6m_IM,
        'WHITE_BLOOD_CELL_COUNT_AT_DIAGNOSIS': (WBC_count_at_diagnosis - 125161.634615) / 102512.604218, #standardization
    }

    df = pd.DataFrame(columns=list(list(variable_dict.keys())))
    df = df.append(variable_dict, ignore_index=True)

    cols_Mol = ['Mol_3m_IM','Mol_6m_IM']
    df['Mol_median'] = df[cols_Mol].median(axis=1)
    df['Mol_median'] = (df['Mol_median'] - 10.273899) / 16.618399 #standardization


    ## loading models
    my_seeds_end = 2040
    prob = []
    pred = []
    my_seeds=range(2020, my_seeds_end) # the random_state that controls the shuffling applied to the data before applying the split
    for seed in my_seeds:
        filename = f'clf_{seed}'
        clf = f'clf_{seed}'
        clf = pickle.load(open(filename, 'rb'))
        pred.append(clf.predict(df[['Mol_median', 'WHITE_BLOOD_CELL_COUNT_AT_DIAGNOSIS']].values))
        prob.append(clf.predict_proba(df[['Mol_median', 'WHITE_BLOOD_CELL_COUNT_AT_DIAGNOSIS']].values)[:,1])

    mean = sum(prob) / len(prob)
    print('The probability of Molecular Response achievement at 12 months after start of Imatinib tretment, using patient information up to 6 months after start of Imatinib, is: ' + str(round(mean[0],2)))