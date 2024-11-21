#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday August  29 23:33:04 2023

@author: Pratik
"""

import pandas as pd
import keras
from keras.models import Sequential,Model
from keras.layers import Dense, Dropout, BatchNormalization,Input
from keras.wrappers.scikit_learn import KerasRegressor
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.utils import np_utils
import keras.backend as Kr
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.utils import class_weight
import numpy as np
from numpy import exp
# Library for Gaussian process
# import GPy
##Library for visualization
import matplotlib.pyplot as plt
# %matplotlib inline
# %config InlineBackend.figure_format = 'svg'
import matplotlib;matplotlib.rcParams['figure.figsize'] = (10,7)
import pylab 
import time
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from tqdm.keras import TqdmCallback
# import GPy

# Load nonGaussian datasets and do classification on them 
num_sim = 50

def model_function(df_train, phi, dummy_y, num_class, sim_iteration):
    print("##### Warning messages ######")
    class_weights = class_weight.compute_class_weight('balanced',np.unique(df_train["class"]),
                                             df_train["class"])
    class_weight_dict = dict(enumerate(class_weights))
    # DeepKriging model for continuous data
    model = Sequential()
    # model.add(Dense(100, input_dim = 2,  kernel_initializer='he_uniform', activation='relu'))
    model.add(Dense(100, input_dim = phi.shape[1],  
            kernel_initializer='he_uniform', activation='relu'))
    # model.add(Dropout(rate=0.5))
    # model.add(BatchNormalization())
    model.add(Dense(50, activation='relu'))
    model.add(Dense(50, activation='relu'))
    #     model.add(Dense(100, activation='relu'))
    # model.add(Dense(100, activation='relu'))
    model.add(Dense(50, activation='relu'))
    model.add(Dense(10, activation='relu'))
    model.add(Dense(10, activation='relu'))
    #     model.add(Dense(50, activation='relu'))
    #     model.add(Dense(50, activation='relu'))
    #     model.add(Dense(10, activation='relu'))
    #model.add(Dense(50, activation='relu'))
    #model.add(Dropout(rate=0.5))
    #     model.add(Dense(10, activation='relu'))
    #model.add(BatchNormalization())
    model.add(Dense(10, activation='relu'))
    model.add(Dense(num_class, activation='softmax'))
    NB_START_EPOCHS = 50 
    # NB_START_EPOCHS = 200  # Number of epochs we usually start to train with
#     BATCH_SIZE = 64  
#     fold_no = 1
    optimizer = keras.optimizers.Adam(lr=0.001)
    model.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])

    
    callbacks = [EarlyStopping(monitor='val_accuracy', patience=200),
                 ModelCheckpoint(filepath='indicator_kriging.h5', 
                                 monitor='val_accuracy', save_best_only=True),
                                 TqdmCallback(verbose=1)]
    print("##### End of warning messages ######")
    print('<<<<<<<<<<<<<<<< Fitting DNN-model for %4d-th simulation >>>>>>>>>>>>>>>>>'%(sim_iteration + 1))
    result = model.fit(phi, dummy_y, callbacks=callbacks, class_weight = class_weight_dict,
               validation_split = 0.1, epochs = 500, batch_size = 128, verbose = 0)

    model = keras.models.load_model('indicator_kriging.h5')
    return model



"""main"""
def main():
    for sim in range(num_sim):

        df_loc = pd.read_csv("synthetic_data_simulations/2D_nonGaussian_1200_projection_"+str(sim+1)+".csv", sep = ",")
        df_train,df_test = train_test_split(df_loc, test_size = 0.1, random_state=123)
        df_train.reset_index(drop=True, inplace=True)
        df_test.reset_index(drop=True, inplace=True)
        # Saving the training and testing datasets 

        df_train.to_csv("synthetic_data_simulations/training_data/2D_nonGaussian_1200_projection_"+str(sim+1)+"train.csv",
                        index = False)
        df_test.to_csv("synthetic_data_simulations/testing_data/2D_nonGaussian_1200_projection_"+str(sim+1)+"test.csv",
                        index = False)

        df_train1 = df_train.copy()
        # print(df_train1.head(1))
        df_train1["class"] = df_train1["class"] - 1
        dummy_y = np_utils.to_categorical(df_train1["class"])
        n = dummy_y.shape[1]
        print('Total number of classes %4d' %(n))
        N = len(df_train1)
        print('Training data size %4d' %(N))
        s = np.vstack((df_train1["x"],df_train1["y"])).T

        num_basis = [5**2,7**2,11**2]
        knots_1d = [np.linspace(0,1,int(np.sqrt(i))) for i in num_basis]
        ##Wendland kernel
        K = 0
        phi = np.zeros((N, sum(num_basis)))

        for res in range(len(num_basis)):
            theta = 1/np.sqrt(num_basis[res])*2.5
            knots_s1, knots_s2 = np.meshgrid(knots_1d[res],knots_1d[res])
            knots = np.column_stack((knots_s1.flatten(),knots_s2.flatten()))
            for i in range(num_basis[res]):
                d = np.linalg.norm(s-knots[i,:],axis=1)/theta
                for j in range(len(d)):
                    if d[j] >= 0 and d[j] <= 1:
                        phi[j,i + K] = (1-d[j])**6 * (35 * d[j]**2 + 18 * d[j] + 3)/3
                    else:
                        phi[j,i + K] = 0
            K = K + num_basis[res]



        # Training the model 

        model = model_function(df_train,phi,dummy_y,n,sim)

        # Basis functions for test set 

        N = len(df_test)
        s = np.vstack((df_test["x"],df_test["y"])).T

        knots_1d = [np.linspace(0,1,int(np.sqrt(i))) for i in num_basis]
        ##Wendland kernel
        K = 0
        phi_test = np.zeros((N, sum(num_basis)))

        for res in range(len(num_basis)):
            theta = 1/np.sqrt(num_basis[res])*2.5
            knots_s1, knots_s2 = np.meshgrid(knots_1d[res],knots_1d[res])
            knots = np.column_stack((knots_s1.flatten(),knots_s2.flatten()))
            for i in range(num_basis[res]):
                d = np.linalg.norm(s-knots[i,:],axis=1)/theta
                for j in range(len(d)):
                    if d[j] >= 0 and d[j] <= 1:
                        phi_test[j,i + K] = (1-d[j])**6 * (35 * d[j]**2 + 18 * d[j] + 3)/3
                    else:
                        phi_test[j,i + K] = 0
            K = K + num_basis[res]


        pred = model.predict(phi_test)
        pred_df = pd.DataFrame(pred)
        df_test_preds = pd.concat([df_test,pred_df], axis = 1)

        # Prediction probabilities saved in file 

        df_test_preds.to_csv("Results_DNN/2D_nonGaussian_1200_predictions_"+str(sim+1)+".csv",
                             index = False)

        
if __name__ == '__main__':
    main()    













