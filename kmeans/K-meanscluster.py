from PIL import Image
import matplotlib.pyplot as plt

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import os
import torch
import numpy as np
from sklearn.cluster import KMeans    #MiniBatchKMeans



xy = np.loadtxt('C:/Users/wan/Desktop/dataB/2221.csv', delimiter=',', dtype=np.float32)        #delimiter=','分隔的符号必须改成‘，’，（如‘\t’）不然会出错。
xylen = xy.shape[0]
xydata = torch.from_numpy(xy[:, :])
xydata = xydata.data.numpy()

#XX, YY = np.meshgrid(X, Y)
#k = 36  # 设置颜色聚类的类别个数（我们分别设置8，16，32，64，128进行对比）
clusterresult = KMeans(n_clusters=36, n_init=10, max_iter=2000, tol=0.0001)  # 构造聚类器
C = clusterresult.fit(xydata)
for i in C.labels_:
    with open('C:/Users/wan/Desktop/dataB/KMout123456.csv', 'a+') as f:
        f.write(str(i) + '\n')

print(C.labels_)
'''
xyre = np.loadtxt('C:/Users/wan/Desktop/dataB/KMout1.csv', delimiter='\n', dtype=np.float32)
xyre = torch.from_numpy(xyre[1:, :])
xyre = xyre.data.numpy()

X_tsne = TSNE(n_components=2, random_state=33).fit_transform(xydata)
X_pca = PCA(n_components=2).fit_transform(xydata)

#ckpt_dir = "images"
#if not os.path.exists(ckpt_dir):
    #os.makedirs(ckpt_dir)

plt.figure(figsize=(10, 5))
plt.subplot(121)
plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=xyre, label="t-SNE")
plt.legend()
plt.subplot(122)
plt.scatter(X_pca[:, 0], X_pca[:, 1], c=xyre, label="PCA")
plt.legend()
#plt.savefig('images/digits_tsne-pca.png', dpi=120)
plt.show()

'''




#plt.figure()


'''
ax.set_title('KMeans')
ax.set_xticks(())
ax.set_yticks(())
'''
#plt.text(-3.5, 1.8, 'train time: %.2fs\ninertia: %f' % (t_batch, k_means.inertia_))
#plt.scatter(C, marker='o',  c=C)
#plt.scatter(xydata[:, 0], xydata[:, 1], c=Kmeans.labels_)
#plt.show()


'''
if __name__ == "__main__":
    data = pd.read_csv(r"C:/Users/wan/Desktop/dataB/2221.csv")
    #data = StandardScaler().fit_transform(data)   #标准化
    n = 546
    ''''''
    x, y = data['Time (sec)'], data['Height (m HAE)']
    n = len(x)
    x = np.array(x)
    x = x.reshape(n, 1)  # reshape 为一列
    y = np.array(y)
    y = y.reshape(n, 1)  # reshape 为一列
    data = np.hstack((x, y))  # 水平合并为两列
    
    k = 36  # 设置颜色聚类的类别个数（我们分别设置8，16，32，64，128进行对比）
    cluster = KMeans(n_clusters=k)  # 构造聚类器
    C = cluster.fit_predict(data)
    # C_Image = cluster.fit_predict(data)

    plt.figure()
    plt.scatter(data[:, 0], data[:, 1], marker='o', s=2, c=C)
    plt.show()
'''

