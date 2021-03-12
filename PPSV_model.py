# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 19:06:31 2021

@author: Satanat
"""

import os
os.chdir('/D002277') # set working directoey here
os.getcwd()

percent_split = 0.2 # Spliting data to test set
nfolds = 5 # For n-folds cross validations

# Python 3 or Anaconda 3 is needed.
import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve, confusion_matrix
import os

def rf_param_selection(X_fullData, X_train, y_train, X_test, y_test, nfolds):
    param_grid = {'n_estimators': [50, 100, 150, 250, 300]}
    grid_search = GridSearchCV(RandomForestClassifier(), param_grid, scoring='f1', verbose = 1, n_jobs = -2, cv=nfolds)
    grid_search.fit(X_train, y_train)
    # Generate model by training set
    y_pred = grid_search.predict(X_train)
    results_train = classification_report(y_train, y_pred)
    cf_train = confusion_matrix(y_train, y_pred)
    # Calculate prediction score by test set
    pred_prob = grid_search.predict_proba(X_test)[:, 1]
    y_pre = np.copy(pred_prob)
    y_pre[y_pre >= 0.5] = 1
    y_pre[y_pre < 0.5] = 0
    y_pre = y_pre.astype(int)
    results_test_prob = classification_report(y_test, y_pre)
    cf_test_prob = confusion_matrix(y_test, y_pre)
    # fullData 
    fullData_pred = grid_search.predict_proba(X_fullData)[:, 1]
    return fullData_pred, results_test_prob, cf_test_prob, pred_prob, grid_search.best_params_, results_train, cf_train 


def machines(percent_split, nfolds, L1, L0, X_fullData):
    # 5 iterations
    performance_prob = []
    predScore_disease = []
    for i in range(1,6): # Iterations
        #print('iterations = ', i)
        # test set
        pos_test = L1.sample(int(percent_split*len(L1)), random_state=1) # sample positive test
        neg_test = L0.sample(len(pos_test), random_state=1) # sample negative test set 
        #pos_test.shape 
        #neg_test.shape 
        test_set = pd.concat([pos_test, neg_test])
        #test_set.shape 
        #os.mkdir('Iter'+ str(i))
        #test_set.to_csv(path + '_Iter_' + str(i) + '_test_set.csv', header = True) # save test set
        X_test = test_set.drop(['Label'], axis=1)
        y_test = test_set.loc[:, 'Label']
        y_test.value_counts()
        
        # training set 
        pos_train = L1.loc[set(L1.index).symmetric_difference(set(pos_test.index))].dropna() 
        #len(pos_train) 
        
        # 5 machines
        mat_result_prob = []
        for m in range(1,6): # Machines
            #print('disease = ', loop_dis)
            #print('machine = ', m)
            neg_train = (L0.loc[set(L0.index).symmetric_difference(set(neg_test.index))].dropna()).sample(len(pos_train), random_state=1) 
            #len(neg_train) 
            training_set = (pd.concat([pos_train, neg_train])) # rbind positive and negative train
            #training_set.shape 
            #training_set.to_csv(path + '_Iter_' + str(i) + '_training_set_m_'+ str(m) + '.csv', header = True) # save training set
            X_train = training_set.drop(['Label'], axis=1)
            y_train = training_set.loc[:, 'Label']
            #y_train.value_counts()
            
            # ML
            fullData_pred, results_test_prob, cf_test_prob, pred_prob, best_params_, results_train, cf_train = rf_param_selection(X_fullData, X_train, y_train, X_test, y_test, nfolds)
            
            # Append prediction score for each machine
            mat_result_prob.append(pred_prob)
            
            # Calculate prediction score for each drug-disease pair
            predScore_disease.append(fullData_pred)
            
            
        ## Average prediction score for each iteration
        result_prob = np.mean(mat_result_prob, axis=0)
        
        # Evaluation performance for each iteration
        AUPR = average_precision_score(y_test, result_prob)
        AUC = roc_auc_score(y_test, result_prob)
        Prob_prec, Prob_rec, Prob_thresholds_pr = precision_recall_curve(y_test, result_prob)
        Prob_fpr, Prob_tpr, Prob_thresholds_roc = roc_curve(y_test, result_prob)
        Prob_f1 = 2 * (Prob_prec * Prob_rec) / (Prob_prec + Prob_rec)
        Prob_max_index = np.argwhere(Prob_f1 == max(Prob_f1))[0]
        Prob_threshold = Prob_thresholds_pr[Prob_max_index]
        Prob_y_pre = np.copy(result_prob)
        Prob_y_pre[Prob_y_pre >= Prob_threshold] = 1
        Prob_y_pre[Prob_y_pre <Prob_threshold] = 0
        Prob_y_pre = Prob_y_pre.astype(int)
        PRE = precision_score(y_test, Prob_y_pre)
        REC = recall_score(y_test, Prob_y_pre)
        F1 = f1_score(y_test, Prob_y_pre)
        ACC = accuracy_score(y_test, Prob_y_pre)

        Prob_avgSchema_metric = {'Name':['AUPR', 'AUC', 'PRE', 'REC', 'ACC', 'F1'],
            'Value':[AUPR, AUC, PRE, REC, ACC, F1]}
        performance_prob.append(Prob_avgSchema_metric['Value'])
        #Prob_df = pd.DataFrame(Prob_avgSchema_metric)
        #Prob_cf_all = confusion_matrix(y_test, Prob_y_pre)
        
    # Average performance score for all iteractions
    df_per_prob = {'Name':['AUPR', 'AUC', 'PRE', 'REC', 'ACC', 'F1'],
            'Value':list(np.mean(performance_prob, axis=0))}
    df_per_prob = pd.DataFrame(df_per_prob)
    
    return str(df_per_prob), predScore_disease
    
    


if __name__ == "__main__":
    f = pd.read_table('./Drug-disease matrix.txt', sep = ' ') # Read drug-disease matrix for each disease
    #f.shape 
    X_fullData = f.drop(['Label'], axis=1)
    fullData_drugs = f.index.values
    L1 = f[f.Label == 1] # positive label
    L0 = f[f.Label == 0] # negative label
    #L1.shape 
    #L0.shape 
    Perf_score, Prediction_score = machines(percent_split, nfolds, L1, L0, X_fullData)
    # Average prediction score of each drug-disease pair
    avgPrediction_score = np.mean(Prediction_score, axis=0)   
    fullData_results = np.array((fullData_drugs, avgPrediction_score)) # Save first column is drug's name.
    fullData_results = np.transpose(fullData_results) # transpose 
    fullData_results = fullData_results[fullData_results[:, 1].argsort()[::-1]] # Sort of prediction score 
    
    print("**********************************************************************************************")
    print("The performance scores: ")
    print(str(Perf_score))
    pd.DataFrame(fullData_results) # Prediction score of each drug for a disease
    #pd.DataFrame(fullData_results, columns=['Drugs', 'Prediction_score']).to_csv(.'/Prediction_score.csv', index = False)







    
