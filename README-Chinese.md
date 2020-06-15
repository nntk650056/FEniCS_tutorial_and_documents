# FEniCS_tutorial_and_documents

（FEniCS 最初源自芝加哥大学与查尔姆斯理工大学的一个合作项目，是一个用于求解偏微分方程的开源计算平台，其目标是使用户能够快速地将科学模型转换为高效的有限元代码）简单教程可以参考官网教程主页（https://fenicsproject.org/pub/tutorial/sphinx1/）

本项目中包含的文件包括：

- [FEniCS_documents]：本人硕士论文、FEniCS简单教程（[Solving PDEs in Python——The FEniCS Tutorial I.pdf]、[fenics-tutorial-chinese.pdf]）、FEniCS指导书[The FEniCS Book—Automated Solution of Differential Equations by the_Finite_Element_Method.pdf]
- [Generate_mesh_XML]： 生成.xml格式网格文件的python程序
- [python_programs]：官方案例程序[Example_from_FEniCS_website]、超弹性体拉伸python程序[hyperelastic]、各种水凝胶溶胀python程序[standard_gel_swelling_program]

## 参考

以下安装教程参考FEniCS官网安装教程

（https://fenicsproject.org/download/）

（https://fenics.readthedocs.io/en/latest/installation.html#debian-ubuntu-packages）

以及CSDN博客：

（https://blog.csdn.net/qq_40492373/article/details/105207726）（作者：仰望天空的蚂蚁）

# 安装过程（此教程适用于Linux 以及 MacOS）

FEniCS可以安装在普通账户下,也可以安装在root账户下。加 * 项为必选操作，其他为可选操作。

## 1.创建root账户

Ubuntu首次为root账户设置密码

```
sudo passwd root
```

跟随指令输入当前账户的密码，并设置root账户的密码，接着输入su，进入root账户下。

## 2.更换源为国内镜像站的源

这里推荐清华源、中科大源以及阿里源。以Ubuntu为例首先备份/etc/apt/sources.list文件

```
mv /etc/apt/sources.list /etc/apt/sourses.list.backup
```

通过上述指令使得 sources.list文件变为空白。

```
#修改sources.list文件*
sudo gedit  /etc/apt/sources.list
```

通过上述指令打开 sources.list文件 ，将以下内容粘贴入sources.list文件中，添加源。

(1)Ubuntu 添加源

```
#清华源
# 默认注释了源码镜像以提高 apt update 速度，如有需要可自行取消注释
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ xenial main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ xenial main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ xenial-updates main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ xenial-updates main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ xenial-backports main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ xenial-backports main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ xenial-security main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ xenial-security main restricted universe multiverse
# 预发布软件源，不建议启用
# deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ xenial-proposed main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ xenial-proposed main restricted universe multiverse

#中科大源
deb https://mirrors.ustc.edu.cn/ubuntu/ bionic main restricted universe multiverse
deb https://mirrors.ustc.edu.cn/ubuntu/ bionic-updates main restricted universe multiverse
deb https://mirrors.ustc.edu.cn/ubuntu/ bionic-backports main restricted universe multiverse
deb https://mirrors.ustc.edu.cn/ubuntu/ bionic-security main restricted universe multiverse
deb https://mirrors.ustc.edu.cn/ubuntu/ bionic-proposed main restricted universe multiverse
deb-src https://mirrors.ustc.edu.cn/ubuntu/ bionic main restricted universe multiverse
deb-src https://mirrors.ustc.edu.cn/ubuntu/ bionic-updates main restricted universe multiverse
deb-src https://mirrors.ustc.edu.cn/ubuntu/ bionic-backports main restricted universe multiverse
deb-src https://mirrors.ustc.edu.cn/ubuntu/ bionic-security main restricted universe multiverse
deb-src https://mirrors.ustc.edu.cn/ubuntu/ bionic-proposed main restricted universe multiverse

#阿里云源
deb http://mirrors.aliyun.com/ubuntu/ bionic main restricted universe multiverse
deb http://mirrors.aliyun.com/ubuntu/ bionic-security main restricted universe multiverse
deb http://mirrors.aliyun.com/ubuntu/ bionic-updates main restricted universe multiverse
deb http://mirrors.aliyun.com/ubuntu/ bionic-proposed main restricted universe multiverse
deb http://mirrors.aliyun.com/ubuntu/ bionic-backports main restricted universe multiverse
deb-src http://mirrors.aliyun.com/ubuntu/ bionic main restricted universe multiverse
deb-src http://mirrors.aliyun.com/ubuntu/ bionic-security main restricted universe multiverse
deb-src http://mirrors.aliyun.com/ubuntu/ bionic-updates main restricted universe multiverse
deb-src http://mirrors.aliyun.com/ubuntu/ bionic-proposed main restricted universe multiverse
deb-src http://mirrors.aliyun.com/ubuntu/ bionic-backports main restricted universe multiverse
```

如果遇到：无法锁定管理目录(/var/lib/dpkg/)，是否有其他进程正占用它？

尝试如下命令：

```linux
sudo rm /var/cache/apt/archives/lock
sudo rm /var/lib/dpkg/lock
```

（2）MacOS添加源（参考博客https://blog.csdn.net/xs18952904/article/details/87261603）

```
#清华源
# 替换brew.git:
git -C "$(brew --repo)" remote set-url origin https://mirrors.tuna.tsinghua.edu.cn/git/homebrew/brew.git
# 替换homebrew-core.git:
git -C "$(brew --repo homebrew/core)" remote set-url origin https://mirrors.tuna.tsinghua.edu.cn/git/homebrew/homebrew-core.git
# 应用生效
brew update
# 替换homebrew-bottles:
echo 'export HOMEBREW_BOTTLE_DOMAIN=https://mirrors.tuna.tsinghua.edu.cn/homebrew-bottles' >> ~/.bash_profile
source ~/.bash_profile

#阿里源
# 替换brew.git:
cd "$(brew --repo)"
git remote set-url origin https://mirrors.aliyun.com/homebrew/brew.git
# 替换homebrew-core.git:
cd "$(brew --repo)/Library/Taps/homebrew/homebrew-core"
git remote set-url origin https://mirrors.aliyun.com/homebrew/homebrew-core.git
# 应用生效
brew update
# 替换homebrew-bottles:
echo 'export HOMEBREW_BOTTLE_DOMAIN=https://mirrors.aliyun.com/homebrew/homebrew-bottles' >> ~/.bash_profile
source ~/.bash_profile

#中科大源
# 替换brew.git:
cd "$(brew --repo)"
git remote set-url origin https://mirrors.ustc.edu.cn/brew.git
# 替换homebrew-core.git:
cd "$(brew --repo)/Library/Taps/homebrew/homebrew-core"
git remote set-url origin https://mirrors.ustc.edu.cn/homebrew-core.git
# 应用生效
brew update
# 替换homebrew-bottles:
echo 'export HOMEBREW_BOTTLE_DOMAIN=https://mirrors.ustc.edu.cn/homebrew/homebrew-bottles' >> ~/.bash_profile
source ~/.bash_profile
```

## 3.安装Anaconda *

前往anaconda官网download页面：https://www.anaconda.com/distribution/

下载最新版的anaconda安装文件，Linux对应的文件后缀名是.sh，复制进入自己现在账户的根目录（~/home/username/），在终端中输入命令：（然后一步步选择 yes 或 no 即可）

```
bash Anaconda3-2020.02-Linux-x86_64.sh
```

安装完成后在终端输入：

```
sudo gedit ~/.bashrc
```


在弹出的文件最后添加：

export PATH="/home/xupp/anaconda3/bin:$PATH"
保存文件，关闭文件后在终端输入

```
source ~/.bashrc
```


此时可以在终端中输入python检查是否安装成功，如果出现如下所示的anaconda字样说明安装成功

```
xxxxx@xxxxxxxx:~$ python
Python 3.7.6 (default, Jan  8 2020, 19:59:22) 
[GCC 7.3.0] :: Anaconda, Inc. on linux
Type "help", "copyright", "credits" or "license" for more information.
```

Anaconda 的安装CSDN以及简书上有许多教程，可以自行搜索参考。

重启终端后，用户名前会出现（base）字样，这是因为anaconda的base环境自启动了，可以通过如下命令关闭base环境自启动

```
## 永远关闭base环境自启动
conda config --set auto_activate_base false
 
## 只关闭本次自启动的base环境（事实上可以关闭任何已经进入的环境）
conda deactivate
```

此时，anaconda安装完成。

## 4.安装FEniCS有限元框架主体*

FEniCS有限元框架主体可以选择安装在 base 环境下，也可以选择安装在新环境下，鉴于国内网络问题，强烈建议按照笔者一下的安装方式进行，如果出现错误，在自行参考博客：https://blog.csdn.net/qq_40492373/article/details/105207726 或者官网上的安装教程（https://fenicsproject.org/download/）（https://fenics.readthedocs.io/en/latest/installation.html#debian-ubuntu-packages）。

执行如下命令：

```
conda create -n fenicsproject -c conda-forge python=3.6 jupyter fenics
#或者选择以下代码
conda create -n fenicsproject -c conda-forge python=3.7 jupyter fenics
```

接着就可以安心等待环境安装完毕。

等待安装完成后，可以尝试启动该环境：

```
conda activate fenicsproject
#或者
source activate fenicsproject
```

每一次当你需要调用fenics时，都需要执行conda activate fenicsproject命令。

## 5.安装fenics所需要的库函数*

```
## 通过conda install 或者 pip install 安装python的各种库函数 numpy matplotlib pandas等
conda install numpy
conda install matplotlib
conda install pandas
```

通过 conda list 命令可以查看环境中已经安装的库函数

## 6.安装paraview后处理工具*

ParaView是一款开源、跨平台数据分析和可视化软件（https://zhuanlan.zhihu.com/p/37503100）

```
sudo apt install paraview
```

安装完成后在命令行输入 paraview 即可启动

打开软件后，可以通过点击“Getting Started Guide”打开官方说明手册（都是英文说明）

## 7.安装mshr包

注意此时应该在自己创建的fenicsproject环境中或者是在base环境中

mshr 就是 FEniCS 的网格生成器，当然FEniCS不借助mshr也能生成比较简单的网格。

```
## 进入anaconda3目录
cd anaconda3
 
## 从GitHub上克隆mshr包的源代码
git clone https://bitbucket.org/fenics-project/mshr
 
## 检查cmake版本
cmake --version
 
## 编译安装mshr包
mkdir mshr/build   && cd mshr/build   && cmake .. && make install && cd ../..
## 安装mshr库
cd mshr/python   && pip3 install . && cd ../..
```

其中cmake ..这一步容易出现问题，假如出现projram killed error，不要灰心，重新执行cmake ..命令即可。

如果不是在root用户下安装的anaconda，那么make install这一步一定会报错，报错大致如下

```
Install the project...
-- Install configuration: ""
-- Installing: /usr/local/lib/libmshr.so.2019.2.0.dev0
CMake Error at cmake_install.cmake:60 (file):
  file INSTALL cannot copy file
  "/home/username/anaconda3/mshr/build/libmshr.so.2019.2.0.dev0" to
  "/usr/local/lib/libmshr.so.2019.2.0.dev0".
```

这是因为make安装时需要管理员权限，可以采用在make install前面加sudo解决。

但是这也同时产生了问题，那就是sudo创建的文件，在本用户下是无法访问的，这也就导致build目录下的install_manifest.txt无法被访问，所以当调用的时候会报找不到模组的错误。解决方案如下：

打开build目录下的install_manifest.txt文件，将里面的内容复制出来到一个新的txt文件中，然后删除原本带锁的文件，将你新创建的文件重命名为install_manifest.txt即可，此时文件已经不带锁了。



————————————————
上述方法的贡献者为 CSDN博主「仰望天空的蚂蚁」，非常感谢！

在命令行中输入如下命令：

```
python

>>from mshr import *
```

如果没有报错，说明mshr包已经安装完成。

## 8.安装Python IDE：使用 conda install 或者 pip install 安装 Spyder、jupyter notebook 或者 jupyterLab *

```
## 激活fenics环境
conda activate fenicsproject
 
## 安装spyder
conda install spyder

## 安装jupyterLab
conda install jupyterLab

## 安装jupyter notebook
conda install jupyter notebook
 
## 如果使用conda安装报错的话，尝试如下
pip install spyder jupyterLab jupyter notebook
 
## 安装完成后，尝试运行
spyder
jupyter lab
jupyter notebook

```

## 9.FEniCS网格划分

UnitSquare(nx, ny): 用于定义大小为[0,1]×[0,1]的矩形的网格，nx、ny 分别 为 x、y 方向上划分的网格数量，网格类型为三角形网格。同类型的函数对象还 包括 UnitInterval, UnitCube, UnitCircle, UnitSphere, Interval, Rectangle, Box 等。使用方法如下：

```python
# 1D domains 
mesh = UnitInterval(20) # 20 cells, 21 vertices 
mesh = Interval(20, -1, 1) # domain [-1,1] 
# 2D domains (6x10 divisions, 120 cells, 77 vertices) 
mesh = UnitSquare(6, 10) # "right" diagonal is default 
# The diagonals can be right, left or crossed 
mesh = UnitSquare(6, 10, "left") 
mesh = UnitSquare(6, 10, "crossed") 
# Domain [0,3]x[0,2] with 6x10 divisions and left diagonals 
mesh = Rectangle(0, 0, 3, 2, 6, 10, "left") 
# 6x10x5 boxes in the unit cube, each box gets 6 tetrahedra: 
mesh = UnitCube(6, 10, 5) 
# Domain [-1,1]x[-1,0]x[-1,2] with 6x10x5 divisions 
mesh = Box(-1, -1, -1, 1, 0, 2, 6, 10, 5) 
# 10 divisions in radial directions 
mesh = UnitCircle(10)
mesh = UnitSphere(10)
```

FEniCS 支持外部导入网格文件，文件是后缀为.xml 或者.xml.gz 格式的文本文件，文件的格式为如下所示：

```xml
<?xml version="1.0" encoding="UTF-8"?>

<dolfin xmlns:dolfin="http://fenicsproject.org"> 
    <mesh celltype="tetrahedron" dim="3"> 
        <vertices size="624">
            <vertex index="0" x="0.0" y="0.0" z="0.0" /> 
            <vertex index="1" x="0.0" y="0.0" z="4.0" /> 
        </vertices> 
        <cells size="2217">
            <tetrahedron index="0" v0="432" v1="464" v2="465" v3="433" /> 
            <tetrahedron index="1" v0="465" v1="466" v2="467" v3="464" />
            ... 
        </cells>
    </mesh> 
</dolfin>
```

在上述.xml 文件中给定几何的维度，单元类型，节点信息以及单元信息即可。 可以通过在 ANSYS、Abaqus 中快速定义好需要使用的网格，从对应的工程文件 中读取相关的节点以及单元信息，通过简单的编程即可得到相应的.xml 文件，笔者生成.xml网格文件的 Python 函数(程序在上述[Generate_mesh_XML]文件夹中)：

最后祝大家安装顺利！

