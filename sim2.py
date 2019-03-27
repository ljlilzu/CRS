# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 16:57:47 2016

@author: Longjie Li
"""
import igraph as ig
import math

#计算顶点间的相似度
def similarities(graph, method,test_edge_list, flag=False):

    #CN,AA,RA
    if method == 'CN':
        return common_neighbors_index(graph, flag)
    if method == 'AA':
        return adamic_adar_index(graph, flag)        
    if method == 'RA':
        return resource_allocation_index(graph, flag)

    #CAR,CAA,CRA
    if method == 'CAR':
        return CAR(graph)
    if method == 'CAA':
        return CAA(graph)
    if method == 'CRA':
        return CRA(graph)

    #WIC        
    if method == 'FASTQ_WIC':
        return FASTQ_WIC(graph) 
    if method == 'LOUVAIN_WIC':
        return LOUVAIN_WIC(graph) 

    #YAN
    if method == 'FASTQ_YANCN':
        return FASTQ_YANCN(graph)         
    if method == 'FASTQ_YANAA':
        return FASTQ_YANAA(graph)       
    if method == 'FASTQ_YANRA':
        return FASTQ_YANRA(graph)         
    if method == 'LOUVAIN_YANCN':
        return LOUVAIN_YANCN(graph)          
    if method == 'LOUVAIN_YANAA':
        return LOUVAIN_YANAA(graph)       
    if method == 'LOUVAIN_YANRA':
        return LOUVAIN_YANRA(graph)                 

    #Bridge                  
    if method == 'FASTQ_BRIDEG_CN':
        return FASTQ_BRIDEG_CN(graph)                  
    if method == 'FASTQ_BRIDEG_AA':
        return FASTQ_BRIDEG_AA(graph)          
    if method == 'FASTQ_BRIDEG_RA':
        return FASTQ_BRIDEG_RA(graph)        
    if method == 'LOUVAIN_BRIDEG_CN':
        return LOUVAIN_BRIDEG_CN(graph)                   
    if method == 'LOUVAIN_BRIDEG_AA':
        return LOUVAIN_BRIDEG_AA(graph)          
    if method == 'LOUVAIN_BRIDEG_RA':
        return LOUVAIN_BRIDEG_RA(graph)

    #series of CRCN
    if method == 'FASTQ_CRCN':
        return FASTQ_CRCN(graph)                 
    if method == 'FASTQ_CRAA':
        return FASTQ_CRAA(graph)          
    if method == 'FASTQ_CRRA':
        return FASTQ_CRRA(graph)
    if method == 'LOUVAIN_CRCN':
        return LOUVAIN_CRCN(graph)                  
    if method == 'LOUVAIN_CRAA':
        return LOUVAIN_CRAA(graph)          
    if method == 'LOUVAIN_CRRA':
        return LOUVAIN_CRRA(graph) 

    #our methods
    if method == 'FASTQ_CRSCN':
        return FASTQ_CRSCN(graph)         
    if method == 'FASTQ_CRSAA':
        return FASTQ_CRSAA(graph)       
    if method == 'FASTQ_CRSRA':
        return FASTQ_CRSRA(graph)        
    if method == 'LOUVAIN_CRSCN':
        return LOUVAIN_CRSCN(graph)           
    if method == 'LOUVAIN_CRSAA':
        return LOUVAIN_CRSAA(graph)       
    if method == 'LOUVAIN_CRSRA':
        return LOUVAIN_CRSRA(graph)
        
    else:
        raise Exception('方法错误', method)

###############################################################################

#节点对
def pair(x, y):
    if (x < y):
        return (x, y)
    else:
        return (y, x)


###############################################################################


"""
Link prediction algorithms.
"""
# CN
def common_neighbors_index(G):
    neighbor={}#邻居
    gvs=G.vs()#节点id
    for i in range(len(gvs)):
        neighbor[i] = G.neighbors(i)          
    sim_dict={}#相似性 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):#i不为j的邻居
                s=len(set(neighbor[i]) & set(neighbor[j]))#共同邻居
                if s>0:
                    sim_dict[(i,j)]=s
    return sim_dict
# end def

# AA
def adamic_adar_index(G):
    neighbor={}#邻居
    gvs=G.vs()#节点id
    for i in range(len(gvs)):
        neighbor[i] = G.neighbors(i) 
    sim_dict={}#相似性   
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):#i,j不相连
                s =0
                for z in (set(neighbor[i]) & set(neighbor[j])):
                    s+=1/math.log2(G.degree(z))
                if s>0:
                    sim_dict[(i,j)]=s
    return sim_dict
# end def

# RA
def resource_allocation_index(G):
    neighbor={}#邻居
    gvs=G.vs()#节点id
    for i in range(len(gvs)):
        neighbor[i] = G.neighbors(i) 
    sim_dict={} #相似性  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):#i,j不相连
                s =0
                for z in (set(neighbor[i]) & set(neighbor[j])):
                    s+=1/(G.degree(z))
                if s>0:
                    sim_dict[(i,j)]=s
    return sim_dict
# end def

# CAR
def CAR(G):
    neighbor={}#邻居
    gvs=G.vs()#节点id
    for i in range(len(gvs)):
        neighbor[i] = G.neighbors(i) 
    sim_dict = {}#相似性 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):#i,j不相连
                s = 0
                t = 0
                for w in (set(neighbor[i]) & set(neighbor[j])):
                    s += 1
                    t += len((set(neighbor[i]) & set(neighbor[w])) & (set(neighbor[w]) & set(neighbor[j])))
                if s > 0:
                    t = t / 2.0
                    sim_dict[(i, j)] = s * t
    return sim_dict
# end def

# CAA
def CAA(G):
    neighbor={}#邻居
    gvs=G.vs()#节点id
    for i in range(len(gvs)):
        neighbor[i] = G.neighbors(i) 
    sim_dict={}#相似性 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):#i,j不相连
                s = 0
                for w in (set(neighbor[i]) & set(neighbor[j])):
                    t = len((set(neighbor[i]) & set(neighbor[w])) & (set(neighbor[w]) & set(neighbor[j])))
                    s += t / math.log2(G.degree(w))
                if s > 0:
                    sim_dict[(i,j)] = s        

    return sim_dict
# end def

# CRA
def CRA(G):
    neighbor={}#邻居
    gvs=G.vs()#节点id
    for i in range(len(gvs)):
        neighbor[i] = G.neighbors(i) 
    #相似性
    sim_dict={}
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):#i,j不相连
                s = 0
                for w in (set(neighbor[i]) & set(neighbor[j])):
                    t = len((set(neighbor[i]) & set(neighbor[w])) & (set(neighbor[w]) & set(neighbor[j])))
                    s += t / (G.degree(w))
                if s > 0:
                    sim_dict[(i,j)] = s        
    return sim_dict
# end def

#my method
def FASTQ_CRSCN(G):   
    dendr=G.community_fastgreedy()#fastq
    cs=dendr.as_clustering()#社团   

    #节点：所属社团编号，节点邻居
    comnum_dict = {}
    neighbor={}#邻居
    gvs=G.vs()#节点id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]#节点：所属社团编号
        neighbor[i] = G.neighbors(i)

    #社团邻居
    com_nei={}
    neii={}#社团外邻居
    C_set={}#社团节点
    for i in range(len(cs)):
        neii[i]=set()
        A=set()
        for j in cs[i]: 
            for z in neighbor[j]:
                A.add(z)
                if comnum_dict[j] == comnum_dict[z]:
                    A.remove(z)
                neii[i]=neii[i] | set(A)
                C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i]     
           
    #求CRS
    community_CRS={}
    com_length_min={}#两社团节点和邻居最小值
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            com_length_min[(i,k)]=min(len(j),len(v))   
            community_CRS[(i,k)]=len(com_nei[i] & com_nei[k])/((com_length_min[(i,k)])+0.0)     
            
    #CN
    CN={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                CN[(i,j)]=len(set(neighbor[i]) & set(neighbor[j]))

    #求相似性            
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=(community_CRS[(comnum_dict[i],comnum_dict[j])])*(CN[(i,j)])
                if s>0:
                    sim_dict[(i,j)]=s                                      
    return sim_dict
#end def    

def FASTQ_CRSAA(G):
    dendr=G.community_fastgreedy()#fastq
    cs=dendr.as_clustering()#社团
    #节点：所属社团编号，节点邻居
    comnum_dict = {} #节点：所属社团编号
    neighbor={}#节点邻居
    gvs=G.vs()#节点id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #社团邻居
    com_nei={}
    neii={}#社团外邻居
    C_set={}#社团节点
    for i in range(len(cs)):
        neii[i]=set()
        A=set()
        for j in cs[i]: 
            for z in neighbor[j]:
                A.add(z)
                if comnum_dict[j] == comnum_dict[z]:
                    A.remove(z)
                neii[i]=neii[i] | set(A)
                C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i]     
           
    #求CRS
    community_CRS={}
    com_length_min={}
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            com_length_min[(i,k)]=min(len(j),len(v))   
            community_CRS[(i,k)]=len(com_nei[i] & com_nei[k])/((com_length_min[(i,k)])+0.0)          
                                     
    #AA
    AA={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                AA[(i,j)]=0
                for z in (set(neighbor[i]) & set(neighbor[j])):
                    AA[(i,j)]+=1/math.log2(G.degree(z))

    #相似性               
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=(community_CRS[(comnum_dict[i],comnum_dict[j])])*(AA[(i,j)])
                if s>0:
                    sim_dict[(i,j)]=s                                         
    return sim_dict
#end def 
    
def FASTQ_CRSRA(G):
    dendr=G.community_fastgreedy()#fastq
    cs=dendr.as_clustering()#社团划分
    
    #节点：所属社团编号，节点邻居
    comnum_dict = {}#节点：所属社团编号
    neighbor={}#节点邻居
    gvs=G.vs()#节点id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #社团邻居
    com_nei={}
    neii={}#社团外邻居
    C_set={}#社团节点
    for i in range(len(cs)):
        neii[i]=set()
        A=set()
        for j in cs[i]: 
            for z in neighbor[j]:
                A.add(z)
                if comnum_dict[j] == comnum_dict[z]:
                    A.remove(z)
                neii[i]=neii[i] | set(A)
                C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i]     
           
    #求CRS
    community_CRS={}
    com_length_min={}
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            com_length_min[(i,k)]=min(len(j),len(v))   
            community_CRS[(i,k)]=len(com_nei[i] & com_nei[k])/((com_length_min[(i,k)])+0.0)          

    #RA                                     
    RA={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                RA[(i,j)]=0
                for z in (set(neighbor[i]) & set(neighbor[j])):
                    RA[(i,j)]+=1/(G.degree(z))
                    
    #相似性
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=(community_CRS[(comnum_dict[i],comnum_dict[j])])*(RA[(i,j)])
                if s>0:
                    sim_dict[(i,j)]=s                                
    return sim_dict
#end def
    
def LOUVAIN_CRSCN(G):
    cs=G.community_multilevel()#Louvain，社团划分
    #节点：所属社团编号
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #社团邻居
    com_nei={}
    neii={}#社团外邻居
    C_set={}#社团节点
    for i in range(len(cs)):
        neii[i]=set()
        A=set()
        for j in cs[i]: 
            for z in neighbor[j]:
                A.add(z)
                if comnum_dict[j] == comnum_dict[z]:
                    A.remove(z)
                neii[i]=neii[i] | set(A)
                C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i]     
           
    #求CRS
    community_CRS={}
    com_length_min={}
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            com_length_min[(i,k)]=min(len(j),len(v))   
            community_CRS[(i,k)]=len(com_nei[i] & com_nei[k])/((com_length_min[(i,k)])+0.0)          

    #CN                                     
    CN={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                CN[(i,j)]=len(set(neighbor[i]) & set(neighbor[j]))

    #相似性            
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=(community_CRS[(comnum_dict[i],comnum_dict[j])])*(CN[(i,j)])
                if s>0:
                    sim_dict[(i,j)]=s 
    return sim_dict
#end def    
 
def LOUVAIN_CRSAA(G):
    cs=G.community_multilevel()#Louvain，社团划分
    #节点：所属社团编号
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #社团邻居
    com_nei={}
    neii={}#社团外邻居
    C_set={}#社团节点
    for i in range(len(cs)):
        neii[i]=set()
        A=set()
        for j in cs[i]: 
            for z in neighbor[j]:
                A.add(z)
                if comnum_dict[j] == comnum_dict[z]:
                    A.remove(z)
                neii[i]=neii[i] | set(A)
                C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i]     
           
    #求CRS
    community_CRS={}
    com_length_min={}
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            com_length_min[(i,k)]=min(len(j),len(v))   
            community_CRS[(i,k)]=len(com_nei[i] & com_nei[k])/((com_length_min[(i,k)])+0.0)          

    #AA                                     
    AA={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                AA[(i,j)]=0
                for z in (set(neighbor[i]) & set(neighbor[j])):
                    AA[(i,j)]+=1/math.log2(G.degree(z))

    #相似性
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=(community_CRS[(comnum_dict[i],comnum_dict[j])])*(AA[(i,j)])
                if s>0:
                    sim_dict[(i,j)]=s 
    return sim_dict
#end def
    
def LOUVAIN_CRSRA(G):
    cs=G.community_multilevel()#Louvain，社团划分
    #节点：所属社团编号
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #社团邻居
    com_nei={}
    neii={}#社团外邻居
    C_set={}#社团节点
    for i in range(len(cs)):
        neii[i]=set()
        A=set()
        for j in cs[i]: 
            for z in neighbor[j]:
                A.add(z)
                if comnum_dict[j] == comnum_dict[z]:
                    A.remove(z)
                neii[i]=neii[i] | set(A)
                C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i]     
           
    #求CRS
    community_CRS={}
    com_length_min={}
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            com_length_min[(i,k)]=min(len(j),len(v))   
            community_CRS[(i,k)]=len(com_nei[i] & com_nei[k])/((com_length_min[(i,k)])+0.0)          

    #RA                                     
    RA={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                RA[(i,j)]=0
                for z in (set(neighbor[i]) & set(neighbor[j])):
                    RA[(i,j)]+=1/(G.degree(z))

    #相似性                
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=(community_CRS[(comnum_dict[i],comnum_dict[j])])*(RA[(i,j)])
                if s>0:
                    sim_dict[(i,j)]=s 
    return sim_dict
#end def
    
#WIC    
def FASTQ_WIC(G):
    dendr=G.community_fastgreedy()#fastq
    cs=dendr.as_clustering()#社团划分

    #节点：所属社团编号
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)

    #根据社团求邻居                           
    CN={}
    CN_W={}#节点i,j,z在同一社团
    CN_IC={}#其他情况
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                if i!=j:
                    CN[(i,j)]=len(set(neighbor[i]) & set(neighbor[j]))
                    neighbor_ij=[]
                    CN_W[(i,j)]=0
                    CN_IC[(i,j)]=0
                    for z in set(neighbor[i]) & set(neighbor[j]):
                        if comnum_dict[i]==comnum_dict[j]==comnum_dict[z]:
                            neighbor_ij.append(z)
                            CN_W[(i,j)]=len(set(neighbor_ij))
                        CN_IC[(i,j)]=CN[(i,j)]-CN_W[(i,j)]                               
    #相似性
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=CN_W[(i,j)]/(CN_IC[(i,j)]+0.001)
                if s>0:
                    sim_dict[(i,j)]=s  
    return sim_dict
#end def
    
def LOUVAIN_WIC(G):
    cs=G.community_multilevel()#LOUVAIN，社团划分
    #节点：所属社团编号
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #根据社团求邻居                       
    CN={}
    CN_W={}#节点i,j,z在同一社团
    CN_IC={}#其他情况
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                if i!=j:
                    CN[(i,j)]=len(set(neighbor[i]) & set(neighbor[j]))
                    neighbor_ij=[]
                    CN_W[(i,j)]=0
                    CN_IC[(i,j)]=0
                    for z in set(neighbor[i]) & set(neighbor[j]):
                        if comnum_dict[i]==comnum_dict[j]==comnum_dict[z]:
                            neighbor_ij.append(z)
                            CN_W[(i,j)]=len(set(neighbor_ij))
                        CN_IC[(i,j)]=CN[(i,j)]-CN_W[(i,j)]                               
    #相似性
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=CN_W[(i,j)]/(CN_IC[(i,j)]+0.001)
                if s>0:
                    sim_dict[(i,j)]=s  
    return sim_dict
#end def    
    
#YAN  
def FASTQ_YANCN(G):
    dendr=G.community_fastgreedy()#fastq
    cs=dendr.as_clustering()#社团划分
    #节点：所属社团编号    
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #根据社团求邻居 
    CN_W={}#i,j在同一社团
    CN_O={}#i,j不在同一社团
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j): 
                if comnum_dict[i]==comnum_dict[j]:
                    CN_W[(i,j)]=len(set(neighbor[i]) & set(neighbor[j]))
                else:                    
                    CN_O[(i,j)]=len(set(neighbor[i]) & set(neighbor[j]))
    CNmax_node=max(CN_O.items(), key=lambda x: x[1])[1] # i,j不在同一社团的最大值
    #相似性
    sim_dict={}
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0          
                if len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i]==comnum_dict[j]:
                        s=CN_W[(i,j)]+CNmax_node
                    else:                    
                        s=CN_O[(i,j)]
                if s>0:
                    sim_dict[(i,j)]=s             
    return sim_dict
#end def
    
def FASTQ_YANAA(G):
    dendr=G.community_fastgreedy()#fastq
    cs=dendr.as_clustering()#社团划分
    #节点：所属社团编号 
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #根据社团求AA
    AA_W={}#i,j在同一社团
    AA_O={}#i,j不在同一社团
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                AA_O[(i,j)]=0 
                AA_W[(i,j)]=0 
                for z in (set(neighbor[i]) & set(neighbor[j])):                   
                    if comnum_dict[i]==comnum_dict[j]:
                        AA_W[(i,j)]+=1/math.log2(G.degree(z))
                    else:                    
                        AA_O[(i,j)]+=1/math.log2(G.degree(z))
    AAmax_node=max(AA_O.items(), key=lambda x: x[1])[1]# i,j不在同一社团的最大值
    #相似性
    sim_dict={}
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0          
                if len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i]==comnum_dict[j]:
                        s=AA_W[(i,j)]+AAmax_node
                    else:                    
                        s=AA_O[(i,j)]
                if s>0:
                    sim_dict[(i,j)]=s                  
    return sim_dict
#end def
    
def FASTQ_YANRA(G):
    dendr=G.community_fastgreedy()#fastq
    cs=dendr.as_clustering()#社团划分
    #节点：所属社团编号 
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #根据社团求RA
    RA_W={}#i,j在同一社团
    RA_O={}#i,j不在同一社团
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                RA_O[(i,j)]=0 
                RA_W[(i,j)]=0 
                for z in (set(neighbor[i]) & set(neighbor[j])):                   
                    if comnum_dict[i]==comnum_dict[j]:
                        RA_W[(i,j)]+=1/(G.degree(z)+0.0)
                    else:                    
                        RA_O[(i,j)]+=1/(G.degree(z)+0.0)
    RAmax_node=max(RA_O.items(), key=lambda x: x[1])[1]# i,j不在同一社团的最大值
    #相似性           
    sim_dict={}
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0          
                if len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i]==comnum_dict[j]:
                        s=RA_W[(i,j)]+RAmax_node
                    else:                    
                        s=RA_O[(i,j)]
                if s>0:
                    sim_dict[(i,j)]=s           
    return sim_dict
    
def LOUVAIN_YANCN(G):
    cs=G.community_multilevel()#LOUVAIN，划分社团
    #节点：所属社团编号 
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #根据社团求邻居
    CN_W={}#节点i,j在同一社团
    CN_O={}#节点i,j不在同一社团
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                if comnum_dict[i]==comnum_dict[j]:
                    CN_W[(i,j)]=len(set(neighbor[i]) & set(neighbor[j]))
                else:                    
                    CN_O[(i,j)]=len(set(neighbor[i]) & set(neighbor[j]))
    CNmax_node=max(CN_O.items(), key=lambda x: x[1])[1]# i,j不在同一社团的最大值
    #相似性             
    sim_dict={}
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0          
                if len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i]==comnum_dict[j]:
                        s=CN_W[(i,j)]+CNmax_node
                    else:                    
                        s=CN_O[(i,j)]
                if s>0:
                    sim_dict[(i,j)]=s    
    return sim_dict  
#end def
    
def LOUVAIN_YANAA(G):
    cs=G.community_multilevel()#LOUVAIN，划分社团
    #节点：所属社团编号 
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #根据社团求AA
    AA_W={}#节点i,j在同一社团
    AA_O={}#节点i,j不在同一社团
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                AA_O[(i,j)]=0 
                AA_W[(i,j)]=0 
                for z in (set(neighbor[i]) & set(neighbor[j])):                   
                    if comnum_dict[i]==comnum_dict[j]:
                        AA_W[(i,j)]+=1/math.log2(G.degree(z))
                    else:                    
                        AA_O[(i,j)]+=1/math.log2(G.degree(z))
    AAmax_node=max(AA_O.items(), key=lambda x: x[1])[1]#节点i,j不在同一社团最大值
    #相似性
    sim_dict={}
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0          
                if len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i]==comnum_dict[j]:
                        s=AA_W[(i,j)]+AAmax_node
                    else:                    
                        s=AA_O[(i,j)]
                if s>0:
                    sim_dict[(i,j)]=s                  
    return sim_dict
#end def
   
def LOUVAIN_YANRA(G):
    cs=G.community_multilevel()#LOUVAIN，划分社团
    #节点：所属社团编号 
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #根据社团求RA
    RA_W={}#节点i,j在同一社团
    RA_O={}#节点i,j不在同一社团
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                RA_O[(i,j)]=0 
                RA_W[(i,j)]=0 
                for z in (set(neighbor[i]) & set(neighbor[j])):                   
                    if comnum_dict[i]==comnum_dict[j]:
                        RA_W[(i,j)]+=1/(G.degree(z)+0.0)
                    else:                    
                        RA_O[(i,j)]+=1/(G.degree(z)+0.0)
    RAmax_node=max(RA_O.items(), key=lambda x: x[1])[1]#节点i,j不在同一社团最大值
    #相似性          
    sim_dict={}
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0          
                if len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i]==comnum_dict[j]:
                        s=RA_W[(i,j)]+RAmax_node
                    else:                    
                        s=RA_O[(i,j)]
                if s>0:
                    sim_dict[(i,j)]=s                 
    return sim_dict  
#end def
               
#Bridge   
def FASTQ_BRIDEG_CN(G):
    dendr=G.community_fastgreedy()#fastq
    P=dendr.as_clustering()#社团划分
    #节点：所属社团编号
    comnum_dict = {}
    neighbor={}#邻居
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=P.membership[i]
        neighbor[i] = G.neighbors(i)
    #找桥节点    
    betweenEdges=[]#跨社团的边
    for i,j in G.get_edgelist():
        if comnum_dict[i]!=comnum_dict[j]:
            betweenEdges.append((i,j))    
    bridgeNode = list()#初始桥节点
    for u,v in betweenEdges:
        if u not in bridgeNode:
            bridgeNode.append(u)
        if v not in bridgeNode:
            bridgeNode.append(v)

    #设置度阈值
    N = ig.Graph.vcount(G)#顶点数
    M = ig.Graph.ecount(G)#边数  
    averagedegree=(2*M)/N#平均度   
    for u in bridgeNode:
        degreeu=G.degree(u)
        if degreeu < averagedegree:
            bridgeNode.remove(u)#移除小于平均度的点
           
    community_length = len(P)#社团数
    LinkNum = {}#节点在某社团连接的边数
    for i in bridgeNode:
        for j in range(community_length):
            LinkNum[(i,j)] = 0
        for w in G.neighbors(i):
            for j in range(community_length):
                if comnum_dict[w] == j:
                    LinkNum[(i,j)] = LinkNum[(i,j)] + 1              
    
    MCDR={}#社团主导率
    max_comnode={}#与一社团连接的最大边数
    for (u,v) in LinkNum:
        m_list=[]#节点与不同社团连接的边
        for (x,y) in LinkNum:
            if u==x :
                m_list.append(LinkNum[(u,v)])
                m_list.append(LinkNum[(x,y)])
        max_comnode[u]=max(m_list)
        MCDR[u]=max_comnode[u]/(G.degree(u)+0.0) 
    #筛选桥节点        
    leavebridgeV=[]
    for u in bridgeNode:
        if MCDR[u]<0.7:
            leavebridgeV.append(u)    
    #相似性
    sim_dict={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=len(set(neighbor[i]) & set(neighbor[j])) 
                if i in leavebridgeV or j in leavebridgeV:
                    s=2*s       
                if s>0:
                    sim_dict[(i,j)]=s 
    return sim_dict    
        
def FASTQ_BRIDEG_AA(G):
    dendr=G.community_fastgreedy()#fastq
    P=dendr.as_clustering()#社团划分 
    #节点：所属社团编号   
    comnum_dict = {}
    neighbor={}#邻居
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=P.membership[i]
        neighbor[i] = G.neighbors(i)
    #找桥节点    
    betweenEdges=[]#跨社团的边
    for i,j in G.get_edgelist():
        if comnum_dict[i]!=comnum_dict[j]:
            betweenEdges.append((i,j))    
    bridgeNode = list()#初始桥节点
    for u,v in betweenEdges:
        if u not in bridgeNode:
            bridgeNode.append(u)
        if v not in bridgeNode:
            bridgeNode.append(v)
    #设置度阈值
    N = ig.Graph.vcount(G)
    M = ig.Graph.ecount(G) 
    averagedegree=(2*M)/N#平均度   
    for u in bridgeNode:
        degreeu=G.degree(u)
        if degreeu < averagedegree:
            bridgeNode.remove(u) #移除小于平均度的点      

    community_length = len(P)#社团数
    LinkNum = {}#节点与某社团连接的边数
    for i in bridgeNode:
        for j in range(community_length):
            LinkNum[(i,j)] = 0
        for w in G.neighbors(i):
            for j in range(community_length):
                if comnum_dict[w] == j:
                    LinkNum[(i,j)] = LinkNum[(i,j)] + 1              
    
    MCDR={}#社团主导率
    max_comnode={}#与一社团连接的最大边数
    for (u,v) in LinkNum:
        m_list=[]#节点与不同社团连接的边
        for (x,y) in LinkNum:
            if u==x :
                m_list.append(LinkNum[(u,v)])
                m_list.append(LinkNum[(x,y)])
        max_comnode[u]=max(m_list)
        MCDR[u]=max_comnode[u]/(G.degree(u)+0.0) 
    #筛选桥节点    
    leavebridgeV=[]
    for u in bridgeNode:
        if MCDR[u]<0.7:
            leavebridgeV.append(u)
    #相似性        
    sim_dict={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0
                for z in set(neighbor[i]) & set(neighbor[j]):
                    s+=1/(math.log2(G.degree(z))+0.0) 
                if i in leavebridgeV or j in leavebridgeV:
                    s=2*s       
                if s>0:
                    sim_dict[(i,j)]=s         
    return sim_dict     
    
def FASTQ_BRIDEG_RA(G):
    dendr=G.community_fastgreedy()#fastq
    P=dendr.as_clustering()#社团划分   
    #节点：所属社团编号
    comnum_dict = {}
    neighbor={}#邻居
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=P.membership[i]
        neighbor[i] = G.neighbors(i)
    #找桥节点    
    betweenEdges=[]
    for i,j in G.get_edgelist():
        if comnum_dict[i]!=comnum_dict[j]:
            betweenEdges.append((i,j))    
    bridgeNode = list()#初始桥节点
    for u,v in betweenEdges:
        if u not in bridgeNode:
            bridgeNode.append(u)
        if v not in bridgeNode:
            bridgeNode.append(v)
                          
    #设置度阈值
    N = ig.Graph.vcount(G)
    M = ig.Graph.ecount(G)  
    averagedegree=(2*M)/N#平均度   
    for u in bridgeNode:
        degreeu=G.degree(u)
        if degreeu < averagedegree:
            bridgeNode.remove(u)#移除小于平均度的点 
     
    community_length = len(P)#社团数
    LinkNum = {}#节点与某社团连接的边数
    for i in bridgeNode:
        for j in range(community_length):
            LinkNum[(i,j)] = 0
        for w in G.neighbors(i):
            for j in range(community_length):
                if comnum_dict[w] == j:
                    LinkNum[(i,j)] = LinkNum[(i,j)] + 1              
    MCDR={}#社团主导率
    max_comnode={}#与一社团连接的最大边数
    for (u,v) in LinkNum:
        m_list=[]#节点与不同社团连接的边
        for (x,y) in LinkNum:
            if u==x :
                m_list.append(LinkNum[(u,v)])
                m_list.append(LinkNum[(x,y)])
        max_comnode[u]=max(m_list)
        MCDR[u]=max_comnode[u]/(G.degree(u)+0.0) 
    #筛选桥节点    
    leavebridgeV=[]
    for u in bridgeNode:
        if MCDR[u]<0.7:
            leavebridgeV.append(u)
    #相似性        
    sim_dict={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0
                for z in set(neighbor[i]) & set(neighbor[j]):
                    s+=1/(G.degree(z)+0.0) 
                if i in leavebridgeV or j in leavebridgeV:
                    s=2*s       
                if s>0:
                    sim_dict[(i,j)]=s  
    return sim_dict    


def LOUVAIN_BRIDEG_CN(G):
    P=G.community_multilevel()#LOUVAIN，社团划分
    #节点：所属社团编号
    comnum_dict = {}
    neighbor={}#邻居
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=P.membership[i]
        neighbor[i] = G.neighbors(i)
    #找桥节点    
    betweenEdges=[]
    for i,j in G.get_edgelist():
        if comnum_dict[i]!=comnum_dict[j]:
            betweenEdges.append((i,j))    
    bridgeNode = list()
    for u,v in betweenEdges:
        if u not in bridgeNode:
            bridgeNode.append(u)
        if v not in bridgeNode:
            bridgeNode.append(v)
                    
    #设置度阈值
    N = ig.Graph.vcount(G)
    M = ig.Graph.ecount(G)  
    averagedegree=(2*M)/N   
    for u in bridgeNode:
        degreeu=G.degree(u)
        if degreeu < averagedegree:
            bridgeNode.remove(u)
        
    community_length = len(P)
    LinkNum = {}#节点与某社团连接的边数
    for i in bridgeNode:
        for j in range(community_length):
            LinkNum[(i,j)] = 0
        for w in G.neighbors(i):
            for j in range(community_length):
                if comnum_dict[w] == j:
                    LinkNum[(i,j)] = LinkNum[(i,j)] + 1              
    #社团主导率
    MCDR={}
    max_comnode={}
    for (u,v) in LinkNum:
        m_list=[]#节点与不同社团连接的边
        for (x,y) in LinkNum:
            if u==x :
                m_list.append(LinkNum[(u,v)])
                m_list.append(LinkNum[(x,y)])
        max_comnode[u]=max(m_list)
        MCDR[u]=max_comnode[u]/(G.degree(u)+0.0) 
    #筛选桥节点    
    leavebridgeV=[]
    for u in bridgeNode:
        if MCDR[u]<0.7:
            leavebridgeV.append(u)
    #相似性    
    sim_dict={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=len(set(neighbor[i]) & set(neighbor[j])) 
                if i in leavebridgeV or j in leavebridgeV:
                    s=2*s       
                if s>0:
                    sim_dict[(i,j)]=s
    return sim_dict    
       
def LOUVAIN_BRIDEG_AA(G):
    P=G.community_multilevel()#LOUVAIN，社团划分
    #节点：所属社团编号
    comnum_dict = {}
    neighbor={}#邻居
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=P.membership[i]
        neighbor[i] = G.neighbors(i)
    #找桥节点    
    betweenEdges=[]
    for i,j in G.get_edgelist():
        if comnum_dict[i]!=comnum_dict[j]:
            betweenEdges.append((i,j))    
    bridgeNode = list()
    for u,v in betweenEdges:
        if u not in bridgeNode:
            bridgeNode.append(u)
        if v not in bridgeNode:
            bridgeNode.append(v)
                     
    #设置度阈值
    N = ig.Graph.vcount(G)
    M = ig.Graph.ecount(G) 
    averagedegree=(2*M)/N   
    for u in bridgeNode:
        degreeu=G.degree(u)
        if degreeu < averagedegree:
            bridgeNode.remove(u)
       
    community_length = len(P)
    LinkNum = {}#节点与某社团连接的边数
    for i in bridgeNode:
        for j in range(community_length):
            LinkNum[(i,j)] = 0
        for w in G.neighbors(i):
            for j in range(community_length):
                if comnum_dict[w] == j:
                    LinkNum[(i,j)] = LinkNum[(i,j)] + 1              
    #社团主导率
    MCDR={}
    max_comnode={}
    for (u,v) in LinkNum:
        m_list=[]
        for (x,y) in LinkNum:
            if u==x :
                m_list.append(LinkNum[(u,v)])
                m_list.append(LinkNum[(x,y)])
        max_comnode[u]=max(m_list)
        MCDR[u]=max_comnode[u]/(G.degree(u)+0.0) 
    #筛选桥节点    
    leavebridgeV=[]
    for u in bridgeNode:
        if MCDR[u]<0.7:
            leavebridgeV.append(u)
    #相似性        
    sim_dict={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0
                for z in set(neighbor[i]) & set(neighbor[j]):
                    s+=1/(math.log2(G.degree(z))+0.0) 
                if i in leavebridgeV or j in leavebridgeV:
                    s=2*s       
                if s>0:
                    sim_dict[(i,j)]=s
    return sim_dict     
    
def LOUVAIN_BRIDEG_RA(G):
    P=G.community_multilevel()#LOUVAIN，划分社团
    #节点：所属社团编号
    comnum_dict = {}
    neighbor={}#邻居
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=P.membership[i]
        neighbor[i] = G.neighbors(i)
    #找桥节点    
    betweenEdges=[]
    for i,j in G.get_edgelist():
        if comnum_dict[i]!=comnum_dict[j]:
            betweenEdges.append((i,j))    
    bridgeNode = list()
    for u,v in betweenEdges:
        if u not in bridgeNode:
            bridgeNode.append(u)
        if v not in bridgeNode:
            bridgeNode.append(v)
                       
    #设置度阈值
    N = ig.Graph.vcount(G)
    M = ig.Graph.ecount(G)  
    averagedegree=(2*M)/N   
    for u in bridgeNode:
        degreeu=G.degree(u)
        if degreeu < averagedegree:
            bridgeNode.remove(u)
       
    community_length = len(P)
    LinkNum = {}#节点与某社团连接的边数
    for i in bridgeNode:
        for j in range(community_length):
            LinkNum[(i,j)] = 0
        for w in G.neighbors(i):
            for j in range(community_length):
                if comnum_dict[w] == j:
                    LinkNum[(i,j)] = LinkNum[(i,j)] + 1              
    #社团主导率
    MCDR={}
    max_comnode={}
    for (u,v) in LinkNum:
        m_list=[]
        for (x,y) in LinkNum:
            if u==x :
                m_list.append(LinkNum[(u,v)])
                m_list.append(LinkNum[(x,y)])
        max_comnode[u]=max(m_list)
        MCDR[u]=max_comnode[u]/(G.degree(u)+0.0) 
    #筛选桥节点    
    leavebridgeV=[]
    for u in bridgeNode:
        if MCDR[u]<0.7:#0.8,0.9
            leavebridgeV.append(u)
    #相似性        
    sim_dict={}  
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0
                for z in set(neighbor[i]) & set(neighbor[j]):
                    s+=1/(G.degree(z)+0.0) 
                if i in leavebridgeV or j in leavebridgeV:
                    s=2*s       
                if s>0:
                    sim_dict[(i,j)]=s  
    return sim_dict    
    
#CRCN,CRAA,CRRA    
def FASTQ_CRCN(G):
    dendr=G.community_fastgreedy()#fastq
    cs=dendr.as_clustering()#社团

    #节点：所属社团编号
    comnum_dict = {}
    neighbor={}
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #社团和社团邻居
    com_nei={}
    neii={}
    C_set={}
    for i in range(len(cs)):
        neii[i]=set()
        for j in cs[i]:
            neii[i]=neii[i] | set(G.neighbors(j))
            C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i] 
    
    #求CRCN
    community_CRCN={}
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            if i!=k:
                community_CRCN[(i,k)]=len(com_nei[i] & com_nei[k])       
    #相似性                          
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0
                if i in comnum_dict and j in comnum_dict and len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i] == comnum_dict[j]:
                        s=1                       
                    else:
                        s=(community_CRCN[(comnum_dict[i],comnum_dict[j])])
                elif i in comnum_dict and j not in comnum_dict:
                    s = ig.Graph.degree(G,i)/ig.Graph.ecount(G)
                else:
                    s=0
                if s>0:
                    sim_dict[(i,j)]=s 
    return sim_dict
#end def

def FASTQ_CRAA(G):
    dendr=G.community_fastgreedy()#fastq
    cs=dendr.as_clustering()#社团划分
    
    comnum_dict = {}#节点：所属社团编号
    neighbor={}#邻居
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #社团和社团邻居
    com_nei={}
    neii={}
    C_set={}
    for i in range(len(cs)):
        neii[i]=set()
        for j in cs[i]:
            neii[i]=neii[i] | set(G.neighbors(j))
            C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i] 
    
    #求CRAA   
    comni_dict = {}#与节点有连接的社团数量
    for i in range(len(gvs)):
        list_1=[]
        for w in G.neighbors(i):
            for x in range(len(cs)):
                if w in cs[x]:
                    list_1.append(x)
        comni_dict[i] = len(list(set(list_1)))
    community_CRAA={}
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            if i!=k:
                community_CRAA[(i,k)]=0 
                for z in (com_nei[i] & com_nei[k]):                
                    community_CRAA[(i,k)]+=1/math.log2(comni_dict[z])         
    #相似性                                  
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0
                if i in comnum_dict and j in comnum_dict and len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i] == comnum_dict[j]:
                        s=1                       
                    else:
                        s=(community_CRAA[(comnum_dict[i],comnum_dict[j])])
                elif i in comnum_dict and j not in comnum_dict:
                    s = ig.Graph.degree(G,i)/ig.Graph.ecount(G)
                else:
                    s=0
                if s>0:
                    sim_dict[(i,j)]=s 
    return sim_dict
#end def

def FASTQ_CRRA(G):
    dendr=G.community_fastgreedy()#fastq
    cs=dendr.as_clustering()#社团划分
    
    comnum_dict = {}#节点：所属社团编号
    neighbor={}#邻居
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #社团和社团邻居
    com_nei={}
    neii={}
    C_set={}
    for i in range(len(cs)):
        neii[i]=set()
        for j in cs[i]:
            neii[i]=neii[i] | set(G.neighbors(j))
            C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i] 

    #求CRRA    
    comni_dict = {}#与节点有连接的社团数量
    for i in range(len(gvs)):
        list_1=[]
        for w in G.neighbors(i):
            for x in range(len(cs)):
                if w in cs[x]:
                    list_1.append(x)
        comni_dict[i] = len(list(set(list_1)))
    community_CRRA={}
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            if i!=k:
                community_CRRA[(i,k)]=0 
                for z in (com_nei[i] & com_nei[k]): 
                    community_CRRA[(i,k)]+=1/(comni_dict[z]) 
    #相似性
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0
                if i in comnum_dict and j in comnum_dict and len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i] == comnum_dict[j]:
                        s=1                       
                    else:
                        s=(community_CRRA[(comnum_dict[i],comnum_dict[j])])
                elif i in comnum_dict and j not in comnum_dict:
                    s = ig.Graph.degree(G,i)/ig.Graph.ecount(G)
                else:
                    s=0
                if s>0:
                    sim_dict[(i,j)]=s 
    return sim_dict    
#end def
    
def LOUVAIN_CRCN(G):
    cs=G.community_multilevel()#LOUVAIN，社团划分
    
    comnum_dict = {}#节点：所属社团编号
    neighbor={}#邻居
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #社团和社团邻居
    com_nei={}
    neii={}
    C_set={}
    for i in range(len(cs)):
        neii[i]=set()
        for j in cs[i]:
            neii[i]=neii[i] | set(G.neighbors(j))
            C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i] 
    
    #求CRCN
    community_CRCN={}
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            if i!=k:
                community_CRCN[(i,k)]=len(com_nei[i] & com_nei[k])    
    #相似性                                  
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0
                if i in comnum_dict and j in comnum_dict and len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i] == comnum_dict[j]:
                        s=1                       
                    else:
                        s=(community_CRCN[(comnum_dict[i],comnum_dict[j])])
                elif i in comnum_dict and j not in comnum_dict:
                    s = ig.Graph.degree(G,i)/ig.Graph.ecount(G)
                else:
                    s=0
                if s>0:
                    sim_dict[(i,j)]=s 
    return sim_dict
#end def

def LOUVAIN_CRAA(G):
    cs=G.community_multilevel()#LOUVAIN社团划分
    
    comnum_dict = {}#节点：所属社团编号
    neighbor={}#邻居
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #社团和社团邻居
    com_nei={}
    neii={}
    C_set={}
    for i in range(len(cs)):
        neii[i]=set()
        for j in cs[i]:
            neii[i]=neii[i] | set(G.neighbors(j))
            C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i] 
    #求CRAA    
    comni_dict = {}#与节点有连接的社团数量
    for i in range(len(gvs)):
        list_1=[]
        for w in G.neighbors(i):
            for x in range(len(cs)):
                if w in cs[x]:
                    list_1.append(x)
        comni_dict[i] = len(list(set(list_1)))
    community_CRAA={}
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            if i!=k:
                community_CRAA[(i,k)]=0 
                for z in (com_nei[i] & com_nei[k]):
                    community_CRAA[(i,k)]+=1/math.log2(comni_dict[z])         
    #相似性                                 
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0
                if i in comnum_dict and j in comnum_dict and len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i] == comnum_dict[j]:
                        s=1                       
                    else:
                        s=(community_CRAA[(comnum_dict[i],comnum_dict[j])])
                elif i in comnum_dict and j not in comnum_dict:
                    s = ig.Graph.degree(G,i)/ig.Graph.ecount(G)
                else:
                    s=0
                if s>0:
                    sim_dict[(i,j)]=s 
    return sim_dict
#end def

def LOUVAIN_CRRA(G):
    cs=G.community_multilevel()#LOUVAIN，划分社团
    
    comnum_dict = {}#节点：所属社团编号
    neighbor={}#邻居
    gvs=G.vs()#id
    for i in range(len(gvs)):
        comnum_dict[i]=cs.membership[i]
        neighbor[i] = G.neighbors(i)
    #社团和社团邻居
    com_nei={}
    neii={}
    C_set={}
    for i in range(len(cs)):
        neii[i]=set()
        for j in cs[i]:
            neii[i]=neii[i] | set(G.neighbors(j))
            C_set[i]=set(cs[i])
        com_nei[i]=set(cs[i]) | neii[i] 

    #求CRRA    
    comni_dict = {}#与节点有连接的社团数量
    for i in range(len(gvs)):
        list_1=[]
        for w in G.neighbors(i):
            for x in range(len(cs)):
                if w in cs[x]:
                    list_1.append(x)
        comni_dict[i] = len(list(set(list_1)))
    community_CRRA={}
    for i,j in com_nei.items():
        for k,v in com_nei.items():           
            if i!=k:
                community_CRRA[(i,k)]=0 
                for z in (com_nei[i] & com_nei[k]): 
                    community_CRRA[(i,k)]+=1/(comni_dict[z]) 
    #相似性
    sim_dict={} 
    for i in range(len(gvs)):
        for j in range(i+1,len(gvs)):
            if i not in G.neighbors(j):
                s=0
                if i in comnum_dict and j in comnum_dict and len(set(neighbor[i]) & set(neighbor[j]))>0:
                    if comnum_dict[i] == comnum_dict[j]:
                        s=1                       
                    else:
                        s=(community_CRRA[(comnum_dict[i],comnum_dict[j])])
                elif i in comnum_dict and j not in comnum_dict:
                    s = ig.Graph.degree(G,i)/ig.Graph.ecount(G)
                else:
                    s=0
                if s>0:
                    sim_dict[(i,j)]=s 
    return sim_dict    
#end def    
