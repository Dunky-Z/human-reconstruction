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


- Minimization with respect to $W_{new}$  
现在的目标就是求$E_m = E_e+E_g+E_c$最小化。首先用学习的方法得到一个结果，但是这个结果处于数据集所构建的空间，有先验条件的约束，不能表达数据集构建的空间以外的模型。

$$\nabla_{\mathrm{pi}} E_{e}=\sum_{d \in D\left(p_{i}\right)} 4\left(\left(\mathrm{p}_{\mathrm{i}}-\mathrm{p}_{\mathrm{j}}\right)^{2}-\left(l_{t}(d)\right)^{2}\right)\left(\mathrm{p}_{\mathrm{i}}-\mathrm{p}_{\mathrm{j}}\right)$$

- Minimization with respect to $p_i$  
接下来就用网格优化的方法，得到一个数据集构建的空间无法描述的全新的结果。  
如果直接对$E_m = E_e+E_g+E_c$最小化能量处理，可能会导致网格不光滑。所以加入一个平滑项来保证得到模型在人体空间内。

$$
E_s = \sum_{p_{i}\in X_{new}}\sum_{p_j\in {N(p_i)}}(\Delta{\vec{p}_i}-\Delta{\vec{p}_j})^2
$$

$$\nabla_{\mathbf{p}_{\mathbf{i}}} E_{s}=\sum_{p_{j} \in N\left(p_{i}\right)} 2\left(\Delta \mathbf{p}_{\mathbf{i}}-\Delta \mathbf{p}_{\mathbf{j}}\right)$$


## 实现细节
其实我还没看懂作者的两步优化能量函数什么意思。不知道为什么是对$W_{new}$最小化能量。

我现在的想法是，我已经能够用学习的方法生成一个初始模型了，我要做的就是如何进行进一步优化。让每个尺寸对应的顶点构成的边长之和逼近目标尺寸。如果能够求出每个尺寸对应的顶点，就可以用泊松变形的方法将整个模型进行一次变形，保证光滑。

所以我最近在做的是如何让三角网格边长逼近已知的长度。但是没有算出结果，问了其他三个伙伴也没能解决。


为了承接现有工作，暂时没有采用将尺寸信息分为三大类，仍然使用控制点之间的欧式距离作为尺寸信息。$\mathcal{E}$为所有控制点构成的边的集合。

为了避免向量与标量的混合运算，将$l(e)$替换成一个二范数为$l(e)$的辅助向量$\mathbf {d}$:

$$\mathbf {d_e} = \frac{\mathbf{q_i}-\mathbf{q_j}}{\|\mathbf{q_i}-\mathbf{q_j}\|}l_t(e)$$

其中$l_t(e)$暂取目标总长度的平均$l_t(e) = l_t/num(e)$。

能量函数可以改写为：

$$
E_{c}=\sum_{e \in \mathcal{E}}\left(\left(\vec {p}_{i}-\vec {p}_{j}\right)-\mathbf{d_e}\right)^{2}
$$

整理成矩阵形式：

$$A\mathbf {x} = \mathbf {d}$$

其中矩阵$A$大小为$3E\times 3V$,$\mathbf{d}$大小为$3E\times 1$，$E$为所有边的个数，$V$为模型顶点个数。求一把发现方程欠定，求不出来。接下来想把$V$改为尺寸相关的顶点，而不是所有顶点求解。因为我现在不需要知道整个模型的顶点位置，我只需要知道19个尺寸相关顶点的位置，然后用泊松变形进行进一步求解剩余顶点位置。

不知道这个方法可不可行。