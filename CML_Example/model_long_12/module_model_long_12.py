def model_long_12():
    
    import numpy as np
    import pandas as pd
    import csv
    import scipy.io
    from scipy import stats
    import matplotlib.pyplot as plt
    from datetime import datetime
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
    #while Mol_3m_IM is None:
    try:
        user_input = input("Enter BCR-ABL1 ratio (percantage) at 3 months after start of Imatinib (or 'None' to represent None): ")
        if user_input.lower() == 'none':
            Mol_3m_IM = None
        else:
            Mol_3m_IM = float(user_input)
    except ValueError:
        print('The input value for BCR-ABL1 ratio (percantage) at 3 months is invalid.')
        Mol_3m_IM = float(input("Enter BCR-ABL1 ratio (percantage) at 3 months after start of Imatinib: "))


    Mol_6m_IM = None
    # while Mol_6m_IM is None:
    try:
        user_input = input("Enter BCR-ABL1 ratio (percantage) at 6 months after start of Imatinib (or 'None' to represent None): ")
        if user_input.lower() == 'none':
            Mol_6m_IM = None
        else:
            Mol_6m_IM = float(user_input)
    except ValueError:
        print('The input value for BCR-ABL1 ratio (percantage) at 6 months is invalid.')
        Mol_6m_IM = float(input("Enter BCR-ABL1 ratio (percantage) at 6 months after start of Imatinib: "))

    Mol_9m_IM = None
    # while Mol_9m_IM is None:
    try:
        user_input = input("Enter BCR-ABL1 ratio (percantage) at 9 months after start of Imatinib (or 'None' to represent None): ")
        if user_input.lower() == 'none':
            Mol_9m_IM = None
        else:
            Mol_9m_IM = float(user_input)
    except ValueError:
        print('The input value for BCR-ABL1 ratio (percantage) at 9 months is invalid.')
        Mol_9m_IM = float(input("Enter BCR-ABL1 ratio (percantage) at 9 months after start of Imatinib: "))


    Mol_12m_IM = None
    # while Mol_12m_IM is None:
    try:
        user_input = input("Enter BCR-ABL1 ratio (percantage) at 12 months after start of Imatinib (or 'None' to represent None): ")
        if user_input.lower() == 'none':
            Mol_12m_IM = None
        else:
            Mol_12m_IM = float(user_input)
    except ValueError:
        print('The input value for BCR-ABL1 ratio (percantage) at 12 months is invalid.')
        Mol_12m_IM = float(input("Enter BCR-ABL1 ratio (percantage) at 12 months after start of Imatinib: "))


    COMPLETE_MOLECULAR_RESPONSE = None
    while COMPLETE_MOLECULAR_RESPONSE is None:
        try:
            COMPLETE_MOLECULAR_RESPONSE = int(input("Enter 1 if patient achieved Complete Molecular Response (CMR), otherwise enter 0: "))
        except ValueError:
            print('The input value for Complete Molecular Response (CMR) is invalid.')
            COMPLETE_MOLECULAR_RESPONSE = int(input("Enter 1 if patient achieved Complete Molecular Response (CMR), otherwise enter 0: "))



    date_of_CMR = None
    date_of_start_of_IM = None
    loss_of_CMR = None
    date_of_loss_of_CMR = None

    if COMPLETE_MOLECULAR_RESPONSE == 1:


        while date_of_CMR is None: #or isinstance(date_of_CMR, datetime) == False:
            try:
                datetime_str = input("Enter date of the patient's CMR achievement, in the format of MM/DD/YYYY : ")
                date_of_CMR = datetime.strptime(datetime_str, '%m/%d/%Y')
            except ValueError:
                try:
                    print("The input value for date of the patient's CMR achievement is invalid.")
                    datetime_str = input("Enter date of the patient's CMR achievement, in the format of MM/DD/YYYY : ")
                    date_of_CMR = datetime.strptime(datetime_str, '%m/%d/%Y')
                except ValueError:
                    print("The input value for date of the patient's CMR achievement is invalid.")
                    datetime_str = input("Enter date of the patient's CMR achievement, in the format of MM/DD/YYYY : ")
                    date_of_CMR = datetime.strptime(datetime_str, '%m/%d/%Y')


        while date_of_start_of_IM is None:
            try:
                datetime_str = input("Enter start date of Imatinib treatment, in the format of MM/DD/YYYY : ")
                date_of_start_of_IM = datetime.strptime(datetime_str, '%m/%d/%Y')
            except ValueError:
                try: 
                    print("The input value for start date of Imatinib treatment is invalid.")
                    datetime_str = input("Enter start date of Imatinib treatment, in the format of MM/DD/YYYY : ")
                    date_of_start_of_IM = datetime.strptime(datetime_str, '%m/%d/%Y')
                except ValueError:
                    print("The input value for start date of Imatinib treatment is invalid.")
                    datetime_str = input("Enter start date of Imatinib treatment, in the format of MM/DD/YYYY : ")
                    date_of_start_of_IM = datetime.strptime(datetime_str, '%m/%d/%Y')


        while loss_of_CMR is None:
            try:
                loss_of_CMR = int(input("Enter 1 if loss of CMR is recorded for the patient, otherwise enter 0: "))
            except ValueError:
                print("The input value for loss of CMR is invalid.")
                loss_of_CMR = int(input("Enter 1 if loss of CMR is recorded for the patient, otherwise enter 0: "))


        if loss_of_CMR == 1:

            while date_of_loss_of_CMR is None:
                try:
                    datetime_str = input("Enter date of loss of CMR in the format of MM/DD/YYYY : ")
                    date_of_loss_of_CMR = datetime.strptime(datetime_str, '%m/%d/%Y')
                except ValueError:
                    try:
                        print("The input value for date of loss of CMR is invalid.")
                        datetime_str = input("Enter date of loss of CMR in the format of MM/DD/YYYY : ")
                        date_of_loss_of_CMR = datetime.strptime(datetime_str, '%m/%d/%Y')
                    except ValueError:
                        print("The input value for date of loss of CMR is invalid.")
                        datetime_str = input("Enter date of loss of CMR in the format of MM/DD/YYYY : ")
                        date_of_loss_of_CMR = datetime.strptime(datetime_str, '%m/%d/%Y')


    ## dataframe        
    variable_dict = {
        'Mol_3m_IM' : Mol_3m_IM,
        'Mol_6m_IM': Mol_6m_IM,  
        'Mol_9m_IM': Mol_9m_IM,
        'Mol_12m_IM': Mol_12m_IM,
        'COMPLETE_MOLECULAR_RESPONSE' : COMPLETE_MOLECULAR_RESPONSE,
        'date_of_MR' : date_of_CMR,
        'start_of_IM' : date_of_start_of_IM,
        'loss_of_MR' : loss_of_CMR,
        'date_loss_MR' : date_of_loss_of_CMR,
    }

    df = pd.DataFrame(columns=list(variable_dict.keys()))
    df = pd.concat([df, pd.DataFrame(variable_dict, index=[0])], ignore_index=True)

    # Compute Mol_median
    cols_Mol = ['Mol_3m_IM','Mol_6m_IM','Mol_9m_IM','Mol_12m_IM']
    df['Mol_median'] = df[cols_Mol].median(axis=1)
    df['Mol_median'] = (df['Mol_median'] - 3.37) / 8.95 #standardization



    # if COMPLETE_MOLECULAR_RESPONSE = 1
    df['MR_0_3'] = 0
    if df['COMPLETE_MOLECULAR_RESPONSE'].values == 1:

        df['date_of_MR'] = pd.to_datetime(df['date_of_MR'])
        df['start_of_IM'] = pd.to_datetime(df['start_of_IM'])

        date_of_MR = df['date_of_MR']
        start_of_IM = df['start_of_IM']

        df['3_months'] = 3
        # 3_mo
        df['start_of_IM_3months'] = start_of_IM + df['3_months'].values.astype("timedelta64[M]")
        start_of_IM_3months = df['start_of_IM_3months']


        date_of_MR_before_start_of_IM = (date_of_MR<=start_of_IM)
        date_of_MR_start_of_IM_to_3_mo = (start_of_IM<=date_of_MR) & (date_of_MR<=start_of_IM_3months)


        condition_0 = (df['COMPLETE_MOLECULAR_RESPONSE']==1)
        # date_of_MR 
        # before_start_of_IM:
        condition = condition_0 & date_of_MR_before_start_of_IM
        df.loc[condition, 'MR_0_3'] = 1
        # start_of_IM_to_3_mo:
        condition = condition_0 & date_of_MR_start_of_IM_to_3_mo
        df.loc[condition, 'MR_0_3'] = 1

    #     print(df[['MR_0_3']])



    # if df['loss_of_MR'] == 1
    if df['loss_of_MR'].values == 1:

        loss_of_MR = df['loss_of_MR']

        df['date_loss_MR'] = pd.to_datetime(df['date_loss_MR'])
        date_loss_MR = df['date_loss_MR']

        date_loss_MR_before_start_of_IM = (date_loss_MR<=start_of_IM)
        date_loss_MR_start_of_IM_to_3_mo = (start_of_IM<=date_loss_MR) & (date_loss_MR<=start_of_IM_3months)


        condition_ = (df['COMPLETE_MOLECULAR_RESPONSE']==1) & (loss_of_MR==1)

        condition_0 = condition_ & (date_of_MR<=start_of_IM)
        # loss_of_MR 
        # before_start_of_IM:
        condition = condition_0 & date_loss_MR_before_start_of_IM
        df.loc[condition, 'MR_0_3'] = 0
        # start_of_IM_to_3_mo:
        condition = condition_0 & date_loss_MR_start_of_IM_to_3_mo
        df.loc[condition, 'MR_0_3'] = 0


        condition_0 = condition_ & (start_of_IM<=date_of_MR) & (date_of_MR<=start_of_IM_3months)
        # loss_of_MR 
        # start_of_IM_to_3_mo:
        condition = condition_0 & date_loss_MR_start_of_IM_to_3_mo
        df.loc[condition, 'MR_0_3'] = 0




    #  standardization
    df['MR_0_3'] = (df['MR_0_3'] - 0.06) / 0.24


    # Calculate the probability
    my_seeds_end = 2025
    prob = []
    pred = []
    my_seeds=range(2020, my_seeds_end) # the random_state that controls the shuffling applied to the data before applying the split
    for seed in my_seeds:
        filename = f'clf_{seed}'
        clf = f'clf_{seed}'
        clf = pickle.load(open(filename, 'rb'))
        pred.append(clf.predict(df[['MR_0_3', 'Mol_median']].values))
        prob.append(clf.predict_proba(df[['MR_0_3', 'Mol_median']].values)[:,1])

    mean = sum(prob) / len(prob)
    print('The probability of Deep Molecular Response achievement, using patient information up to 12 months after start of Imatinib, is: ' + str(round(mean[0],2)))