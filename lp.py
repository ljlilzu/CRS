# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 17:00:14 2016

for  Link Prediction

@author: CNNVD
"""

import igraph as ig
import numpy as np
#import math
import random
import datetime
import sim2
'''
Link Prediction
@param graph_file   网络文件
@param out_file     输出文件，文件格式，已打开
@param sim_method   计算相似性的方法
@param t            独立运行的实验的次数
@param p            测试集的比例
'''
def LP(graph_file, out_file, sim_method, t, p):  
    G = ig.Graph.Read_Pajek(graph_file)#网络为.net文件
#    G = ig.Graph.Read_GML(graph_file)#网络为.gml文件 

    node_num = ig.Graph.vcount(G)#节点数
    edge_num = ig.Graph.ecount(G)#边数


    # 列出所有不存在的链接，存放到non_edge_list中
    non_edge_list=[]
    gvs=G.vs()
    for u in range(len(gvs)):
        for v in range(u+1,len(gvs)):
            if u not in G.neighbors(v):                
                u, v = sim2.pair(u, v)
                non_edge_list.append((u, v))
    non_edge_num = len(non_edge_list)

    # for debug
    print("V: %d\tE: %d\tNon: %d" % (node_num, edge_num, non_edge_num))

    # 执行t次独立的实验，每次从G中选择p*100%的链接作为测试集，剩余的链接作为训练集
    test_num = int(edge_num * p)
    pre_num = 0

    for l in range(10, 101, 10):
        if l < test_num:
            pre_num += 1
        else:
            break
        # end if
    # end for
    pre_num += 1

    # for debug
    print('test_edge_num: %d' % test_num)

    # 定义数组存放性能值
    auc_list = []
    rs_list = []
    time_list = []
#    pre_matrix = [[0 for it in range(t)] for num in range(pre_num)]


    # 迭代t次进行测试
    for it in range(t):
        if it % 10 == 0:
            print('turn: %d' % it)
        # end if

        # 首先产生一批随机数
        random.seed(a=None)
        rand_set = set(random.sample(range(edge_num), test_num))

        # 遍历G中链接，根据rand_set中的值分成训练集和测试集
        training_graph = ig.Graph()
        training_graph.add_vertices(range(node_num))
        test_edge_list = []

        r = 0
        for u, v in G.get_edgelist():
            u, v = sim2.pair(u, v)
            if r in rand_set:  # 测试链接
                test_edge_list.append((u, v))
                
            else:
                training_graph.add_edge(u, v)  # 训练网络
            # end if
            r += 1
        # end for         
        training_graph.to_undirected()

        # 计算相似度        
        start = datetime.datetime.now()
        sim_dict = sim2.similarities(training_graph, sim_method,test_edge_list)
        end = datetime.datetime.now()
        
        # 0. 计算时间
        time_list.append((end - start).microseconds)

        # 1. 计算AUC
        auc_value = AUC(sim_dict, test_edge_list, non_edge_list)
        auc_list.append(auc_value)
        # for debug

        # 创建一个数组，存放顶点对的相似度
        sim_list = [((u, v), s) for (u, v), s in sim_dict.items()]

        # sim_dict不在需要
        sim_dict.clear()

        # 对sim_list按照相似度降序排列
        sim_list.sort(key=lambda x: (x[1], x[0]), reverse=True)

        # 2. 计算Ranking Score
        rank_score = Ranking_score(sim_list, test_edge_list, non_edge_num)
        rs_list.append(rank_score)
        # for debug

        # 3. 计算精度列表
#        pre_list = Precision(sim_list, test_edge_list, test_num)
#
#        for num in range(pre_num):
#            pre_matrix[num][it] = pre_list[num]
        # end for
    # end for

    # 计算平均值和方差，并将结果输出到文件
    auc_avg, auc_std = stats(auc_list)
    
    print('AUC: %.4f(%.4f)' % (auc_avg, auc_std))
    out_file.write('%.4f(%.4f)\t' % (auc_avg, auc_std))

    rs_avg, rs_std = stats(rs_list)

    print('Ranking_Score: %.4f(%.4f)' % (rs_avg, rs_std))
    out_file.write('%.4f(%.4f)\t' % (rs_avg, rs_std))

    time_avg, time_std = stats(time_list)

    print('Time: %.4f(%.4f)' % (time_avg, time_std))
    out_file.write('%.4f(%.4f)\t' % (time_avg, time_std))

#    pre_avg_list = []
#    pre_std_list = []
#    for num in range(pre_num):
#        pre_avg, pre_std = stats(pre_matrix[num])
#        pre_avg_list.append(pre_avg)
#        pre_std_list.append(pre_std)
    # end for

#    print('Precision: ')
#    for num in range(pre_num):
#        print('%.4f(%.4f)\t' % (pre_avg_list[num], pre_std_list[num]))
#        out_file.write('%.4f(%.4f)\t' % (pre_avg_list[num], pre_std_list[num]))
    # end for    
    out_file.write('%d\n' % test_num)   
# end def

# 输入列表，计算平均值和方差
def stats(value_list):
    value_array = np.array(value_list)
    avg = np.mean(value_array)
    std = np.std(value_array)

    return avg, std
# end def

###############################################################################
"""
精度计算的函数
"""
# @param sim_dict 存放顶点对相似度的字典
# @param node_num 顶点个数
# @param missing_edge_list 测试集，丢失的链接  $E^p$
# @param non_edge_list 不存在的链接  $U - E$


def AUC(sim_dict, missing_edge_list, non_edge_list):
    if len(missing_edge_list) * len(non_edge_list) <= 10000:
        return auc1(sim_dict, missing_edge_list, non_edge_list)
    else:
        return auc2(sim_dict, missing_edge_list, non_edge_list)
    # end if
# end AUC


###############################################################################
# 计算ACU值，该方法中将测试集中的边与不存在的边进行两两比较
# @param sim_dict 存放顶点对的相似度，字典
# @param missing_edge_list 测试集，丢失的链接
# @param non_edge_list 不存在的链接
# @return auc值
def auc1(sim_dict, missing_edge_list, non_edge_list):

    n1 = 0
    n2 = 0
    for (u, v) in missing_edge_list:
        try:
            m_s = int(sim_dict[(u, v)] * 1000000)
        except KeyError:
            m_s = 0
        # end try
        for (x, y) in non_edge_list:
            try:
                n_s = int(sim_dict[(x, y)] * 1000000)
            except KeyError:
                n_s = 0
            # end try

            if m_s > n_s:
                n1 += 1
            elif m_s == n_s:
                n2 += 1
    
            # end if
        # end for
    # end for

    n = len(missing_edge_list) * len(non_edge_list)
    return (n1 + 0.5 * n2) / n        
# end def

# 计算ACU值，该方法进行10000次比较
def auc2(sim_dict, missing_edge_list, non_edge_list):

    n = 10000
    n1 = 0
    n2 = 0
    
    m_num = len(missing_edge_list)
    n_num = len(non_edge_list)
    
    for i in range(n):
        r1 = random.randint(0, m_num - 1)
        r2 = random.randint(0, n_num - 1)
        
        (u, v) = missing_edge_list[r1]
        (x, y) = non_edge_list[r2]

        try:
            m_s = int(sim_dict[(u, v)] * 1000000)
        except KeyError:
            m_s = 0
        # end try

        try:
            n_s = int(sim_dict[(x, y)] * 1000000)
        except KeyError:
            n_s = 0
        # end try

        if m_s > n_s:
            n1 += 1
        elif m_s == n_s:
            n2 += 1
        # end if
    # end for
    return (n1 + 0.5 * n2) / n        
# end def

###############################################################################
#计算Precision
def Precision(sim_list, missing_edge_list, missing_edge_num):

    # 将missing_edge_list转换成set
    missing_edge_set = set(missing_edge_list)

    pre_list = []

    count = 0
    ll = len(sim_list)
    for l in range(missing_edge_num):
        if l < ll:
            (u, v) = sim_list[l][0]
            if (u, v) in missing_edge_set:
                # (u, v) 是一条丢失的边
                count += 1
            # end if
        # end if

        if (l + 1) % 10 == 0 and l < 100:   # 输出top-(l+1)
            pre_list.append(count / (l + 1))
        # end if
    # end for
    pre_list.append(count / missing_edge_num)

    return pre_list
# end def
    
###############################################################################

def Ranking_score(sim_list, missing_edge_list, non_edge_num):

    missing_edge_num = len(missing_edge_list)
    
    H = missing_edge_num + non_edge_num
     
    # 定义rank_dict，存放预测的每条链接的rank
    rank_dict = {}
    
    # 变量sim_list,得到rank值
    for r in range(len(sim_list)):
        (u, v) = sim_list[r][0]
        rank_dict[(u, v)] = r + 1
    # end for
    
    rr = H - 1
    
    sum_rank = 0
    for (u, v) in missing_edge_list:
        try:
            rank = rank_dict[(u, v)]
        except KeyError:
            rank = rr   # 没有相似度的边
        # end try
        sum_rank += rank
    # end for
            
    return sum_rank / (missing_edge_num * H)
# end ranking_score
