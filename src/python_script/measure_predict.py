import os
import time
import random
import numpy as np

from fancyimpute import MICE
from fancyimpute import KNN
from fancyimpute import SoftImpute
from sklearn.feature_selection import RFE
from sklearn.linear_model import LinearRegression
from multiprocessing import Pool
import scipy
import scipy.sparse
import scipy.sparse.linalg

import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

V_NUM = 12500
F_NUM = 25000
M_NUM = 19
BASIS_NUM = 10

# DATA_DIR = "D://ITabc//ITabc//human-reconstruction//build//data//train_tmp/"
DATA_DIR = "D://ITabc//ITabc//human-reconstruction//build//data//train//male/"
#DATA_DIR = "D://ITabc//ITabc//BodyReconstruction//BodyReconstruction//build//data//train-1000"

rows = 3
cols = 25000
# facets = np.fromfile(DATA_DIR + "/facets", dtype=int)[4:]
mean_measure = np.fromfile(DATA_DIR + "/mean_measure", dtype=float)[2:]
std_measure = np.fromfile(DATA_DIR + "/std_measure", dtype=float)[2:]
# measure_list = np.fromfile(DATA_DIR + "/measure_list_py", dtype=float)[2:]

measure_list = np.fromfile(DATA_DIR + "/measure_list", dtype=float)[2:]
# facets = np.array(facets).reshape(cols, rows)
mean_measure = np.array(mean_measure).reshape(19, 1)
std_measure = np.array(std_measure).reshape(19, 1)
num_model = int(measure_list.shape[0] / 19)
measure_list = np.array(measure_list).reshape(num_model, M_NUM)
measure_list = np.transpose(measure_list)


# using imputation for missing data
def get_predict(flag, in_data, t_measure, method="MICE"):
    '''
    调用不同求解器，预测缺失信息
    Args:
        flag:标注哪些尺寸为缺失
        in_data:原始数据
        t_measure:尺寸数据集（经过中心化， 标准化）
        method:采用的方法。默认为MICE（MICE，KNN，SoftImpute）
    Returns:
        预测的完整尺寸信息

    '''
    output = in_data.copy()
    output.shape = (M_NUM, 1)
    output[~flag] = np.nan
    solver = MICE()
    if method == "SoftImpute":
        solver = SoftImpute()
    elif method == "KNN":
        print("KNN")
        solver = KNN()
    tmp = t_measure.copy()
    tmp = np.column_stack((tmp, output)).transpose()
    tmp = solver.complete(tmp)
    output = np.array(tmp[-1, :]).reshape(M_NUM, 1)
    return output


def predict(data, method):
    '''
    输入原始尺寸，和预测方法名称， 获得预测尺寸
    Args:
        data:原始尺寸信息
        method:预测方法名（MICE，KNN，SoftImpute）

    Returns:
        预测的完整尺寸信息

    '''
    data = np.array(data).reshape(M_NUM, 1)
    mask = np.zeros((M_NUM, 1), dtype=bool)
    for i in range(0, M_NUM):
        if data[i, 0] != 0:
            data[i, 0] -= mean_measure[i, 0]
            data[i, 0] /= std_measure[i, 0]
            mask[i, 0] = 1
    input_data = get_predict(mask, data, measure_list, method)
    # output = input
    return input_data


def predict_with_mask(data, mask_list, method):
    '''
    给定一组mask信息，通过mask信息，完成一系列预测
    Args:
        data:原始尺寸信息
        mask_list:mask列表，包含18个子列表，分别为缺失1-18个数据的mask
        method:预测方法名（MICE，KNN，SoftImpute）

    Returns:
        18组完整的预测尺寸

    '''
    data = np.array(data).reshape(M_NUM, 1)
    input_data_list = []
    """ 生成18组尺寸信息，分别是缺失1-18个尺寸的预测结果"""
    for j in range(0, M_NUM):
        data[j, 0] -= mean_measure[j, 0]
        data[j, 0] /= std_measure[j, 0]
    for i in range(18):
        input_data = get_predict(mask_list[i], data, measure_list, method)
        input_data_list.append(input_data)
    # output = input
    return input_data_list


def generate_random_mask_list():
    '''
    随机生成缺失1-18个数据的mask列表
    Returns:
        mask_list:分别为缺失1个数据，2个数据，。。。，18个数据的mask
        random_index_list:记录随机生成的缺失数据下标，方便后面原始信息与预测信息做差值时查找。

    '''
    mask_list = []
    random_index_list = []

    for i in range(1, M_NUM):
        mask = np.ones((M_NUM, 1), dtype=bool)
        list_t = random.sample(range(0, 18), i)
        random_index_list.append(list_t)
        for j in range(i):
            mask[list_t[j], 0] = 0
        mask_list.append(mask)
    return mask_list, random_index_list


def get_mean_error_list(ori_data, input_data_list, random_index_list):
    '''
    计算预测一次的平均误差
    Args:
        ori_data:原始尺寸信息
        input_data_list:预测尺寸信息
        random_index_list:随机生成的下标信息

    Returns:
        预测一次的平均误差

    '''
    error_mean = []
    for i in range(len(input_data_list)):
        for j in range(0, M_NUM):
            input_data_list[i][j] *= std_measure[j, 0]
            input_data_list[i][j] += mean_measure[j, 0]
    for i in range(len(input_data_list)):
        error = 0
        for j in range(len(random_index_list[i])):
            xx = float(input_data_list[i][random_index_list[i][j]])
            error = error + abs(ori_data[random_index_list[i][j]] - xx)
        error_mean.append(error / len(random_index_list[i]))
    return error_mean


def draw_error_diagram(mean_error):
    '''
    画折线图
    Args:
        mean_error:平均误差列表，包含三种方法分别预测N次后的平均误差

    Returns:None

    '''
    x_label_mice = range(1, 19)
    x_label_knn = range(1, 19)
    x_label_soft = range(1, 19)
    y_label_mice = mean_error[0]
    y_label_knn = mean_error[1]
    y_label_soft = mean_error[2]
    # 把x轴的刻度间隔设置为1，并存在变量里
    x_major_locator = MultipleLocator(1)
    # ax为两条坐标轴的实例
    ax = plt.gca()
    # 把x轴的主刻度设置为1的倍数
    ax.xaxis.set_major_locator(x_major_locator)
    # 把x轴的刻度范围设置为-0.5到11，因为0.5不满一个刻度间隔，所以数字不会显示出来，但是能看到一点空白
    plt.xlim(-0.5, 19)
    # 设置图表标题，并给坐标轴加上标签
    plt.xlabel("The Number of Missing Measurements", fontsize=14)
    plt.ylabel("Mean Value", fontsize=14)
    plt.plot(x_label_mice, y_label_mice, label="MICE", marker='o', color='b')
    plt.plot(x_label_knn, y_label_knn, label="KNN", marker='s', color='g')
    plt.plot(x_label_soft, y_label_soft, label="SoftImpute", marker='^', color='r')
    # 左上方增加图例
    plt.legend(loc='upper left')
    # 绘制图表中虚线网格
    plt.grid(x_label_mice)
    return


def get_mean_error(measure, times):
    '''
    输入一个原始尺寸数据，调用预测函数，计算出预测times次后数据的平均误差
    Args:
        measure:原始尺寸数据
        times:预测次数

    Returns:
        返回一个shape为[3,18]的列表，分别为3种方法缺失1-18个数据的预测times次平均误差
    '''
    mean_error_list_mice = []
    mean_error_list_knn = []
    mean_error_list_soft = []
    for i in range(times):
        mask_list, random_index_list = generate_random_mask_list()
        input_data_list_mice = predict_with_mask(measure, mask_list, "MICE")
        input_data_list_knn = predict_with_mask(measure, mask_list, "KNN")
        input_data_list_soft = predict_with_mask(measure, mask_list, "SoftImpute")
        mean_error_list_mice.append(get_mean_error_list(measure, input_data_list_mice, random_index_list))
        mean_error_list_knn.append(get_mean_error_list(measure, input_data_list_knn, random_index_list))
        mean_error_list_soft.append(get_mean_error_list(measure, input_data_list_soft, random_index_list))
    mean_error_mice = [float(sum(col)) / len(col) for col in zip(*mean_error_list_mice)]
    mean_error_knn = [float(sum(col)) / len(col) for col in zip(*mean_error_list_knn)]
    mean_error_soft = [float(sum(col)) / len(col) for col in zip(*mean_error_list_soft)]
    return [mean_error_mice, mean_error_knn, mean_error_soft]


def add():
    '''
    无用的测试函数
    Returns:

    '''
    a = 1
    b = 2
    print(a + b)


if __name__ == "__main__":
    # 4698.85,
    measure = [4698.85, 1795.61, 460.47, 1212.81, 1098.78, 1134.35, 890.41, 823.41, 419.05, 824.58, 1126.354, 1199.55,
              1336.46, 649.92, 423.889, 204.25, 1313.27, 442.89, 726.47]

    mean_error = get_mean_error(measure, 50)
    draw_error_diagram(mean_error)
    add()
