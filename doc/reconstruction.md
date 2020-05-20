### 生成数据集
#### 对主成分采样
使用三角形变形Deformation作为特征，所有特征构成矩阵$D_{9m \times n}$，其中$m$为三角面片数量，$n$为数据集模型数量。如果我们通过SVD找到了矩阵$DD^T$最大的$k$个特征向量张成的$9m\times k$维矩阵$U$。则我们如果进行如下处理：$D'_{k\times n}= U^T_{k\times 9m}D_{9m \times n}$。可以得到一个$k\times n$的矩阵，就完成了降维操作。  
只要能获得不同的$D$即可生成不同身材的人模型，于是只要对$U_{9m \times k}$进行高斯采样，得到新采样的$U'$，
由

$$D_{9m\times n} = U'_{9m\times k}D'_{k\times n}$$

升维到Deformation信息。

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
![](https://gitee.com/dominic_z/markdown_picbed/raw/master/img/01.jpg)
![](https://gitee.com/dominic_z/markdown_picbed/raw/master/img/02.jpg)

##### Female
![](https://gitee.com/dominic_z/markdown_picbed/raw/master/img/001.jpg)
![](https://gitee.com/dominic_z/markdown_picbed/raw/master/img/002.jpg)