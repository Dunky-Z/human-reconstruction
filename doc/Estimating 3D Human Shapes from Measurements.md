# Estimating 3D Human Shapes from Measurements阅读笔记


## 训练模型  
在PCA空间中，每个三维网格$X_i$都由一个向量$W_i$表示，通过特征分析可以得到一个线性映射关系,给定一个新的测量值就可以得到对应的网格在PCA中的权重：

$$W_{new} = BP_{new} \tag{1}$$  

对于每个训练集里的网格$X_i$进行特征分析，可以得到一个平均模型 $\mu$ 和矩阵$A$。有了这两个参数就可以通过一个新的模型的权重$W_{new}$求出对应的模型$X_{new}$：

$$X_{new} = AW_{new}+\mu \tag{2}$$  

由公式（1）（2）给定测量值$P_{new}$就可以得到对应的人体模型：

$$X_{new} = ABP_{new}+\mu \tag{3}$$

## 模型精调
测量的数据分为三类，欧式距离值，测地距离值和围长值，对于每个待求的网格$X_i$，在网格上对应的求出的数据尽量和真实值保持一致，这就是一个最小化能量函数的问题。 
  
$$
E_{e}=\sum_{d \in \mathcal{D}}\left(\left(\vec{p}_{i}-\vec p_{j}\right)^{2}-\left(l_{t}(d)\right)^{2}\right)^{2}
$$

$$
E_{g}=\sum_{e \in \mathcal{P}}\left((\vec {p}_{k}-\vec {p}_{l})^{2}-\left(l_{t}(e)\right)^{2}\right)^{2}
$$

$$
E_{c}=\sum_{e \in \mathcal{C}}\left(\left(\vec {q}_{i}-\vec {q}_{j}\right)^{2}-\left(l_{t}(e)\right)^{2}\right)^{2}
$$

其中$\vec p$是顶点坐标向量，$\vec q$是切平面与网格上三角面片的边的交点坐标向量。$l(d)$是实际测量线段的长度，$l(e)$是实际测量的围长。

在模型预测部分我们知道，只要给定一个全新的降维后的主成分$W_{new}$，就可以带入线性回归模型$X_{new} = AW_{new}+\mu$，得到一个全新的模型。那么接下来目标就是，如何获得一个更准确的$W_{new}$。我们知道$W_{new} = A^{+}(X^{init}_{new} - \mu)$,所以可以从如何获取一个比较好的$X^{init}_{new}$着手。
- Minimization with respect to $W_{new}$  
$X^{init}_{new}$可以通过学习的方法算出一个初始的模型，通过对三类尺寸的优化获得更精确的模型。优化方法采用拟牛顿法，需要对方程求导：$\nabla_{W_{new}}E = A^+\nabla_{p_i}E$

$$\nabla_{\mathrm{pi}} E_{e}=\sum_{d \in D\left(p_{i}\right)} 4\left(\left(\mathrm{p}_{\mathrm{i}}-\mathrm{p}_{\mathrm{j}}\right)^{2}-\left(l_{t}(d)\right)^{2}\right)\left(\mathrm{p}_{\mathrm{i}}-\mathrm{p}_{\mathrm{j}}\right)$$

$$\nabla_{\mathrm{pi}} E_{g}=\sum_{e \in P\left(p_{i}\right)} 4\left(\left(\mathrm{p}_{\mathrm{i}}-\mathrm{p}_{\mathrm{j}}\right)^{2}-\left(l_{t}(e)\right)^{2}\right)\left(\mathrm{p}_{\mathrm{i}}-\mathrm{p}_{\mathrm{j}}\right)$$

$$\nabla_{\mathrm{pi}} E_{c}=\sum_{e \in C\left(p_{i}\right)} 4\left(\left(\mathrm{p}_{\mathrm{i}}-\mathrm{p}_{\mathrm{j}}\right)^{2}-\left(l_{t}(e)\right)^{2}\right)\left(\mathrm{p}_{\mathrm{i}}-\mathrm{p}_{\mathrm{j}}\right)$$
计算出新的$W_{new}$后就可以通过线性回归模型得到新的模型$X_{new}^{pca}=AW_{new}+\mu$。但是这个模型仍然是数据集空间中的模型。
- Minimization with respect to $p_i$  
接下来就用网格优化的方法，得到一个数据集构建的空间无法描述的全新的结果。  
如果直接对$E_m = E_e+E_g+E_c$最小化能量处理，可能会导致网格不光滑。所以加入一个平滑项来保证得到模型在人体空间内。该平滑项的意义是让相邻的顶点都有相似的变形，$\Delta{\mathbf{p}_i}$表示变形前后的平移向量。

$$
E_s = \sum_{p_{i}\in X_{new}}\sum_{p_j\in {N(p_i)}}(\Delta{\mathbf{p}_i}-\Delta{\mathbf{p}_j})^2
$$

$$\nabla_{\mathbf{p}_{\mathbf{i}}} E_{s}=\sum_{p_{j} \in N\left(p_{i}\right)} 2\left(\Delta \mathbf{p}_{\mathbf{i}}-\Delta \mathbf{p}_{\mathbf{j}}\right)$$
这一步的整体能量函数可以表示为$E=(1-\lambda)E_m+\lambda E_s$。并且将顶点初始化为上一步中求得的$X_{new}^{pca}$。
## 实现细节

