### 生成数据集
#### 对主成分采样
使用三角形变形Deformation作为特征，所有特征构成矩阵$D_{9m \times n}$，其中$m$为三角面片数量，$n$为数据集模型数量。如果我们通过SVD找到了矩阵 $DD^T$ 最大的 $k$ 个特征向量张成的 $9m\times k$ 维矩阵 $U$ 。则我们如果进行如下处理：

$$D^{'}_{k\times n}= U^{T}_{k\times 9m}D_{9m \times n}$$

对降维后的结果$D'$每一维进行高斯采样，生成新主成分$D'$,由

$$D_{9m \times n} = U_{9m \times k}D^{'}_{k \times n}$$

即可获取随机的Deformation信息。

设原主成分$U$每一维的标准差为$Std.Dev$，为了避免潜在的畸形体型，对于新采样的$U'$每一维的标准差控制在$\pm 3\times Std.Dev$。
#### 高斯采样 Box-Muller

若随机变量$z$服从标准正态分布$N(0,1)$，令$x = \sigma \cdot x + \mu$，则$x$服从均值为$\mu$,方差为$\sigma^2$的高斯分布$N(\mu, \sigma^2)$。因此，任意高斯分布其实都可以通过标准正态分布通过拉伸和平移得到，所以可以通过标准正态分布采用为例来研究这个问题。

##### 逆变换法  

采样步骤：  
1. 产生[0,1]均匀分布的随机数$\mu$
2. $z = \sqrt{2}erf^{-1}(2\mu - 1)$,$z$服从标准正态分布。其中$erf()$是高斯误差函数，他是标准正态分布累积分布函数经过简单的平移和拉伸变换后的形式，定义如下：

$$\operatorname{erf}(x)=\frac{2}{\sqrt{\pi}} \int_{0}^{x} e^{-t^{2}} d t$$

##### Box-Muller算法
上述逆变换算法需要求解$erf()$的逆函数，因为其不是初等函数，所以没有显示解，计算复杂。为了避免这种非初等函数的求逆操作，Box-Muller算法提出了如下解决方案：既然单个高斯分布的累积分布不好求逆，那么两个独立的高斯分布的联合分布呢？  
假设$x,y$是服从标准正态分布的独立随机变量，它的联合概率密度为：  

$$p(x, y)=\frac{1}{2 \pi} \mathrm{e}^{-\frac{x^{2}+y^{2}}{2}}$$

考虑$(x,y)$在圆盘$(x,y)|x^2+y^2\leq R^2|$上的概率：

$$F(R)=\int_{x^{2}+y^{2} \leq R^{2}} \frac{1}{2 \pi} \mathrm{e}^{-\frac{x^{2}+y^{2}}{2}} \mathrm{d} x \mathrm{d} y$$

通过极坐标变换，求得二重积分：

$$F(R)=1-\mathrm{e}^{-\frac{R^{2}}{2}}, R \ge 0$$


这里$F(R)$可以看成是极坐标中 $r$ 的累积分布函数。由于$F(R)$ 的计算公式比较简单， 逆函数也很容易求得，所以可以利用逆变换法来对r进行采样；对于$\theta$，在$[0,2\pi]$ 上进行均匀采样即可。这样就得到了 $(r,\theta)$ ，经过坐标变换即可得到符合标准正态分布的$(x,y)$。具体采样过程如下：
1. 产生$[0,1]$上的两个独立均匀分布随机数$\mu_1, \mu_2$
2. 令

$$
x=\sqrt{-2 \ln \left(u_{1}\right)} \cos 2 \pi u_{2} \\
y=\sqrt{-2 \ln \left(u_{1}\right)} \sin 2 \pi u_{2}
$$
则$x,y0$服从正态分布，并且相互独立。

想要得到服从$z \sim N(\mu,\sigma^2)$的高斯分布，则只需要对$x \sim N(0,1）$做出如下变换：

$$y = \sigma x + \mu$$


#### 数据集结果
##### Male
![](https://gitee.com/dominic_z/markdown_picbed/raw/master/img/dertcfvyghjbksafdaf.jpg)

##### Female
![](https://gitee.com/dominic_z/markdown_picbed/raw/master/img/女性数据集.jpg)

#### 检查身高是否符合正态分布

利用python计算出所有模型（20000个）的身高，然后绘制出直方图，基本符合正态分布。
##### 男性身高分布  
因为数据集采集的样本是欧美国家的成年人，身高偏高，所有数据集因为没有什么问题。

<head>  

![](https://gitee.com/dominic_z/markdown_picbed/raw/master/img/身高正态分布直方图.png)

#### 女性身高分布
![](https://gitee.com/dominic_z/markdown_picbed/raw/master/img/女性身高分布图.png)

### 测量与预测尺寸
#### 测量尺寸
采用计算控制点之间的欧式距离之和。

#### 预测尺寸

三种已有的矩阵插补方法的比较。（MICE，KNN，SoftImpute）

##### 20000个模型

![](https://gitee.com//dominic_z/markdown_picbed/raw/master/img/predict.png)

##### Song等人的4000个模型

![](https://gitee.com//dominic_z/markdown_picbed/raw/master/img/20200527164745.png)

总体来看误差还要大一些。
我求平均误差的方式是，给定一个模型的原始尺寸。随机丢掉1-18个尺寸，然后用剩余尺寸预测。重复进行100次，然后取平均误差。Song原文里也没说是如何去平均误差的。如果只人为固定丢掉的尺寸，比如身高或者腰围这些尺寸较大的数据，误差就会降低很多。

#### 算法概述
测量的数据分为三类，欧式距离值，测地距离值和围长值，对于每个待求的网格$X_i$，在网格上对应的求出的数据尽量和真实值保持一致，这就是一个最小化能量函数的问题 。对于欧式距离，给定目标长度$l_t(d)$，表示线段$d$两个端点$v_i$,$v_j$的欧式距离。能量函数可以定义为：

$$
E_{\mathcal{D}} = \sum_{d \in \mathcal{D}}\left(\|v_i - v_j \|^{2}-l_{t}(d)^{2}\right)^{2}
$$

其中$\mathcal{D}$就是所有欧式距离尺寸的集合，如身高。

对于测地距离，给定目标尺寸$l_t(P)$表示顶点$v_k,v_l$之间测地路径$P$的测地距离。我们用$l_g(P)$来表示$P$的真实测地距离。我们假定在变形前后，路径$P$上的每个边长$e$相对长度保持不变，有了这个假设，我们就可以得到每个边的目标边长，就可以将问题化成一个优化每条边长的问题。能量函数可以表示为：

$$
E_{\mathcal{P}} = \sum_{e \in \mathcal{P}}\left(\|v_k - v_l \|^{2}-l_{t}(e)^{2}\right)^{2}
$$

加入拉普拉斯能量全局保形，

$$
E =  \sum\left\|L\mathbf{V}' -  L\mathbf{V} \right\|^{2}+
\sum_{d \in \mathcal{D}}\left(\|v_i - v_j \|^{2}-l_{t}(d)^{2}\right)^{2}
+
 \sum_{e \in \mathcal{P}}\left(\|v_k - v_l \|^{2}-l_{t}(e)^{2}\right)^{2}
$$

其中$\mathcal{D},\mathcal{P}$分别为欧式距离边集合，测地距离边集合。 $V$为原模型顶点坐标，$V'$为目标模型顶点坐标。假定变形前后长度比例不变，可以通过下式算出逼近的长度：

$$l_t(v) = \frac{l_t(P)}{l_g(P)}l_g(e)$$

将能量方程降为2次，设$\|\mathbf{d}\|=l_t(d)$，$\mathbf{v}_{ij}:=v_i - v_j$，重写能量函数$E_{\mathcal{D}}$:

$$
E_{\mathcal{D}} = \sum_{d \in \mathcal{D}}\left(\|\mathbf{v}_{ij}\|-\|\mathbf{d}\|\right)^{2}
$$

由反三角不等式可得：

$$
E_{\mathcal{D}} = \sum_{d \in \mathcal{D}}\left(\|\mathbf{v}_{ij}\|-\|\mathbf{d}\|\right)^{2} \leq \sum_{d \in \mathcal{D}}\|\mathbf{v}_{ij}-\mathbf{d}\|^{2}
$$

将$\mathbf{d}=l_t{d}(\mathbf{v}_{ij}/\|\mathbf{v}_{ij}\|)$替换不等式右边的$\mathbf{d}$:

$$
\sum_{d \in \mathcal{D}}\|\mathbf{v}_{ij}-l_t{d}\frac{\mathbf{v}_{ij}}{\|\mathbf{v}_{ij}\|}\|^{2}=\sum_{d \in \mathcal{D}}\|  \frac{\mathbf{v}_{ij}}{\|\mathbf{v}_{ij}\|}(\|\mathbf{v}_{ij}\| - l_t{d})  \|^{2} = \sum_{d \in \mathcal{D}}(\|\mathbf{v}_{ij}\| - l_t{d})^{2}
$$

因此当$\mathbf{d}=l_t{d}(\mathbf{v}_{ij}/\|\mathbf{v}_{ij}\|)$时，方程左右都有相同的极小值。因此重写能量函数为：

$$
E_{\mathcal{D}} = \sum_{d \in \mathcal{D}}\|  (v_{i} - v_{j}) - \mathbf{d}\|
$$

同理可以改写测地距离的能量函数，最终能量函数函数为：

$$
E =  \sum\left\|L\mathbf{V}' -  L\mathbf{V} \right\|^{2}+
\sum_{d \in \mathcal{D}}\|  (v_{i} - v_{j}) - \mathbf{d}_d\|
+
\sum_{e \in \mathcal{P}}\|  (v_{k} - v_{l}) - \mathbf{d}_e\|
$$

改写为矩阵形式为：

$$
\left[\begin{array}{c}
L\\
C_1 \\
C_2
\end{array}\right] V
=\left[\begin{array}{c}
\Delta\\
\mathbf{d}_d \\
\mathbf{d}_e
\end{array}\right]
$$

$C$代表标记边相邻顶点的系数矩阵。
设

$$
A = \left[\begin{array}{c}
\mathcal{L}\\
C_1 \\
C_2
\end{array}\right],
b=\left[\begin{array}{c}
\Delta\\
\mathbf{d}_d \\
\mathbf{d}_e
\end{array}\right]
$$

A大小为：$(3|V| + 3|M|) \times (3|V|)$，V大小为$3|V|\times 1$，b大小为：$(3|V| + 3|M|) \times 1$。$|M|$为尺寸相关的所有边的个数。
在分解矩阵$A^TA$时，内存爆了，占用高达20G，无法分解计算。
