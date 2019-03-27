# -*- coding: utf-8 -*-
"""
Personalized Link Prediction

Created on Fri Oct 21 09:33:16 2016

@author: Longjie Li
"""
import lp

'''----------------------------------------------------------------------------
参数设置
'''
t = 20		# 独立实验的次数
p = 0.1		# 测试数据的比例，测试集的大小 

suf = str(p * 10)
path = r'./Networks/'  # 网络数据的根目录

# 网络文件
networks = [
    'CE/celegansneural.net', 		# 0
    'karate/karate.gml',  		# 1
	'jazz/jazz.net',				# 2
    'netscience/netscience.net',  # 3
    'USAir/USAir97.net',  		# 4
    'yeast/yeast.net',  			# 5 
    'PB/polblogs.net',  			# 6
    'email/email.net',  			# 7
    'word/adjnoun.gml',           # 8
    'football/football.gml',       #9
    'polbooks/polbooks.gml',     #10
    'Facebook/FaceBook.net',     # 11
]

#结果文件
results = [
    r'./results/CE-' + suf,
    r'./results/Karate-' + suf,
    r'./results/Jazz-' + suf,
    r'./results/NS-' + suf,
    r'./results/USAir-' + suf,
    r'./results/Yeast-' + suf,
    r'./results/PB-' + suf,
    r'./results/Email-' + suf,
    r'./results/Word-' + suf,
    r'./results/football-' + suf,
    r'./results/polbooks-' + suf,
    r'./results/Facebook-' + suf,    
]

# 实验中可能只使用部分网络， 下面数组中指定相应网络的id
net_ids = [0,2,3,4,5,6,7,11]    #G = ig.Graph.Read_Pajek(graph_file) 网络为.net文件
#net_ids = [1,8,9,10]    #G = ig.Graph.Read_GML(graph_file) 网络为.gml文件

graph_file_list = []  # 网络文件列表
result_file_list = []  # 结果文件列表

for i in net_ids:
    graph_file_list.append(path + networks[i])
    result_file_list.append(results[i]+ 'maxE'+ suf)
# end for

# 相似性方法
sim_methods = [               
    'CN',  		# 0
    'AA',  		# 1
    'RA',  		# 2
    
    'CAR',		# 3
    'CAA',      # 4
    'CRA',      # 5
    
    'FASTQ_WIC',#6
    'LOUVAIN_WIC',#7

    'FASTQ_BRIDEG_CN', #8
    'FASTQ_BRIDEG_AA',#9
    'FASTQ_BRIDEG_RA',#10
    'LOUVAIN_BRIDEG_CN',#11
    'LOUVAIN_BRIDEG_AA',#12
    'LOUVAIN_BRIDEG_RA',#13
    
    'FASTQ_YANCN',#14
    'FASTQ_YANAA',#15
    'FASTQ_YANRA',#16
    'LOUVAIN_YANCN',#17
    'LOUVAIN_YANAA',#18
    'LOUVAIN_YANRA',#19

    'FASTQ_CRCN',#20
    'FASTQ_CRAA',#21
    'FASTQ_CRRA',#22  
    'LOUVAIN_CRCN',#23
    'LOUVAIN_CRAA',#24
    'LOUVAIN_CRRA',#25
    
    'FASTQ_CRSCN',#26
    'FASTQ_CRSAA',#27
    'FASTQ_CRSRA',#28    
    'LOUVAIN_CRSCN',#29
    'LOUVAIN_CRSAA',#30
    'LOUVAIN_CRSRA',#31
        
]

# 实验中使用的方法的id
method_ids = [0] #,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
sim_method_list = [sim_methods[i] for i in method_ids]

# 按照数据集，分别计算
for i in range(len(graph_file_list)):
    graph_file = graph_file_list[i]
    result_file = result_file_list[i]
    out_file = open(result_file, 'w')  # 打开结果文件
    print(graph_file)
    # 输出标题
#    out_file.write('Method\tAUC\tRanking_Score\ttime (ms)\tPrecision (10)\n')
    out_file.write('Method\tAUC\tRanking_Score\ttime (ms)\n')
    
    # 按照不同的相似度方法分别计算
    for method in sim_method_list:
        print(method)
        out_file.write(method + '\t')
        lp.LP(graph_file, out_file, method, t, p)
        out_file.flush()
    # end for
    out_file.close()
# end for

###############################################################################
