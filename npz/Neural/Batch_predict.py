#import os
#os.environ['TF_CPP_MIN_LOG_LEVEL']='2'
#%%
from model import Model
import numpy as np

import time

#%%读取 模型
model=Model()
name="best_model2025032415"
model.load(name)
#%% 读输入
data=np.genfromtxt('/Users/xiaobao/Desktop/san/LimitCycle/LimitCycle_1.csv',delimiter=',')
x=np.transpose(data)
print(x.shape)
#%% 计算输出并保存
ans=np.transpose(model.predict(np.array(x)))
np.savetxt('/Users/xiaobao/Desktop/san/Pathway/LimitCycle_1_ans.csv',ans,delimiter=',')


# %%
