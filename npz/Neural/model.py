#@file model.py
# @brief 定义model
# @version 0.3
# @date 2021-09-30
#
# @copyright Copyright (c) 2021
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
import joblib
# import keras
from tensorflow.python.keras.backend import dropout
from tensorflow.python.keras.models import Sequential
import os

class Model:
    __scaler:StandardScaler
    __model:Sequential
    def __init__(self) -> None:
        self.__scaler=StandardScaler()
        self.__model  =tf.keras.models.Sequential([
                        tf.keras.layers.InputLayer(input_shape=(3,)),
                        #tf.keras.layers.Dense(8, activation='relu')               
                        tf.keras.layers.Dense(8, activation='relu'),
                        tf.keras.layers.Dense(16, activation='relu'),
                        #tf.keras.layers.Dense(128, activation='relu'),
                        #tf.keras.layers.Dense(256, activation='relu'),
                        #tf.keras.layers.Dense(16, activation='relu'),
                        #tf.keras.layers.Dense(16, activation='relu'),
                        tf.keras.layers.LayerNormalization(),
                        #tf.keras.layers.Dense(128, activation='relu'),
                        tf.keras.layers.Dense(8, activation='relu'),
                        tf.keras.layers.Dense(3)
                     ])
        loss_fn =tf.keras.losses.MeanAbsoluteError()
        opt=tf.keras.optimizers.Adam(learning_rate=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-07)
        self.__model.compile(optimizer=opt
                    ,loss=loss_fn
                    # ,metrics=['accuracy']
                    )
    def scaler_fit(self,data)-> None:
        self.__scaler.fit(data)
    
    def scaler_transform(self,data):
        return self.__scaler.transform(data)
        # return data
    
    def load(self, name:str)->None:
        self.__scaler=joblib.load('models/%s/%s.bin' % (name,name))
        self.__model=tf.keras.models.load_model('models/%s/%s.h5' % (name,name))

    def save(self,name:str)->None:
        try:
            os.mkdir("models")
        except:
            pass
        try:
            os.mkdir('models/%s' % (name,))
        except:
            pass
        joblib.dump(self.__scaler,'models/%s/%s.bin' % (name,name), compress=True)        
        self.__model.save('models/%s/%s.h5' % (name,name))
    
    def predict(self,data):
        data=self.scaler_transform(data)
        return self.__model.predict(data)
    
    def evaluate(self,x,y):
        x=self.scaler_transform(x)
        return self.__model.evaluate(x,y)
    
    def fit(self,x,y,verbose:int=0,initial_epoch:int=0,epochs:int=1,batch_size:int=64)->None:
        try:
            x=self.scaler_transform(x)
        except:
            self.scaler_fit(x)
            x=self.scaler_transform(x)
        self.__model.fit(x, y,initial_epoch=initial_epoch, epochs=epochs,verbose=verbose,batch_size=batch_size)
    
    def a(self,a):
        return self.__model.get_layer(index=a)

