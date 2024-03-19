def model_3_18():

    import numpy as np
    import pandas as pd
    import csv
    import scipy.io
    from scipy import stats
    import matplotlib.pyplot as plt
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
        'WHITE_BLOOD_CELL_COUNT_AT_DIAGNOSIS': (WBC_count_at_diagnosis - 117775.2) / 101378.98, #standardization
    }

    df = pd.DataFrame(columns=list(variable_dict.keys()))
    df = pd.concat([df, pd.DataFrame(variable_dict, index=[0])], ignore_index=True)

    cols_Mol = ['Mol_3m_IM']
    df['Mol_median'] = df[cols_Mol].median(axis=1)
    df['Mol_median'] = (df['Mol_median'] - 11.24) / 19.28 #standardization         


    ## loading models
    my_seeds_end = 2025
    prob = []
    pred = []
    my_seeds=range(2020, my_seeds_end) # the random_state that controls the shuffling applied to the data before applying the split
    for seed in my_seeds:
        filename = f'clf_{seed}'
        clf = f'clf_{seed}'
        clf = pickle.load(open(filename, 'rb'))
        pred.append(clf.predict(df[['WHITE_BLOOD_CELL_COUNT_AT_DIAGNOSIS', 'Mol_median']].values))
        prob.append(clf.predict_proba(df[['WHITE_BLOOD_CELL_COUNT_AT_DIAGNOSIS', 'Mol_median']].values)[:,1])

    mean = sum(prob) / len(prob)
    print('The probability of Molecular Response achievement at 18 months after start of Imatinib tretment, using patient information up to 3 months after start of Imatinib, is: ' + str(round(mean[0],2)))