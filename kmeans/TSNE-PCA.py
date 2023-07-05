import matplotlib.pyplot as plt
import numpy as np
import torch
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import pandas as pd
from sklearn.preprocessing import StandardScaler
import os



xy = np.loadtxt('C:/Users/wan/Desktop/dataB/2221.csv', delimiter=',', dtype=np.float32)        #delimiter=','分隔的符号必须改成‘，’，（如‘\t’）不然会出错。
xylen = xy.shape[0]
xydata = torch.from_numpy(xy[:, :])
xydata = xydata.data.numpy()

xyre = np.loadtxt('C:/Users/wan/Desktop/dataB/KMout123456.csv', delimiter=' ', dtype=np.float32)
#xyre = torch.from_numpy(xyre[:, 0])
#xyre = xyre.data.numpy()
xyre = xyre.reshape(546, 1)

X_tsne = TSNE(n_components=2).fit_transform(xydata)  #random_state=33
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
