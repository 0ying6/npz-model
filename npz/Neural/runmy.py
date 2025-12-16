# %%
from model import Model
import numpy as np
import time

# %% 初始化 模型
model = Model()
# %% 读取 模型
# name="best_model2022070321"
# model.load(name)
# %% 读取数据
data = np.genfromtxt('/Users/xiaobao/Desktop/san/LimitCycle/shiyan9_random.csv', delimiter=',')
# data=np.genfromtxt('datas/Effecive_Data_new.csv',delimiter=',')

x=np.transpose(data[:3])
y=np.transpose(data[3:])
y=np.expand_dims(y,axis=1)
n=int(len(x)*0.8)
x_train=x[0:n]
x_test=x[n:]
y_train=y[0:n]
y_test=y[n:]
# %%
# model.scaler_fit(x_train)

# %% 训练
# x_train=x[0:100]
# y_train=y[0:100]

es = 0.01
i = 0
maxi = i
best_less = (None, None)
# best_model=Model()
while True:  # maxi<=10000:
    maxi += 100
    epochs = 200
    # loss=model.evaluate(x_train, y_train)

    while True:
        i += 1
        print("#%d..." % i, end="")
        # arr=np.arange(x_train.shape[0])
        # np.random.shuffle(arr)
        # arr=arr[0:100]
        t = time.time()
        model.fit(x_train, y_train,initial_epoch=(i - 1) * epochs, epochs=i * epochs, batch_size=32, verbose=0)
        print("done. 用时%.2fs." % (time.time() - t))
        loss = model.evaluate(x_train, y_train)
        loss_te = model.evaluate(x_test, y_test)
        if (best_less[0] == None or (loss < best_less[0] and loss_te < best_less[1])):
            best_less = (loss, loss_te)
            # best_model.clone(model)
            print("新最佳模型", end="")
            model.save(time.strftime("best_model%Y%m%d%H"))
        if (loss < es or i >= maxi):
            break
    if (loss < es):
        print("成功")
        break
    # print("now model:")
    # model.evaluate(x_train, y_train)
    # model.evaluate(x_test, y_test)
    model.save(time.strftime("model%Y%m%d%H%M%S"))
# %%  测试

model.evaluate(x_train, y_train)
model.evaluate(x_test, y_test)
# %% 实验
#input_3d_data=np.random.randint(2,10,5)
#model.predict(input_3d_data)  # -17.383 ,  -8.9591

model.predict(np.array([[6.21,9.46,10.89]]))
#%% 保存模型
model.save(time.strftime("model%Y%m%d%H%M%S"))


