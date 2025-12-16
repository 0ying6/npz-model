#%%
from time import time
import scipy,io as sio
import csv
import os
import os.path
import traceback
import time
import numpy as np
#%%直接读取mat文件转为csv  多了一个指标行，如何去掉？
#import pandas as pd
#import scipy
#from scipy import io
#features_struct = scipy.io.loadmat('/Users/cindy/Desktop/jy/datas/ODE45.mat')
#features = features_struct['P']
#dfdata = pd.DataFrame(features)
#datapath1 = '/Users/cindy/Desktop/jy/datas/ODE45.csv'
#dfdata.to_csv(datapath1, index=True)

#%% 去重复
filename=("/Users/xiaobao/Desktop/san/LimitCycle/shiyan9")
data=np.genfromtxt('%s.csv' % filename,delimiter=',')
m=dict()
l=[]
for i in range(data.shape[1]):
    t=tuple(np.transpose(data[2:,i]))
    if t in m:
        m[t].append(i)
    else:
        print(t,i)
        l.append(i)
        m[t]=[t,i]
    # print(i)
data=data[:,l]
np.savetxt('%s_free.csv' % filename,data,delimiter=',')
#%% 随机化
filename="/Users/xiaobao/Desktop/san/LimitCycle/shiyan9"
data=np.genfromtxt('%s.csv' % filename,delimiter=',')
arr=np.arange(data.shape[1])
np.random.shuffle(arr)
data=data[:,arr]
np.savetxt('%s_random.csv' % filename,data,delimiter=',')

# %%
