"""
@Author:weidong wu
@Author_Email:weidongwu404@gmail.com
@Time:2022/10/12 下午2:10
"""

"""
muation rate about time:
22 aa, length S:
S,N,M,E 
1274,420,223,76

:param  
    input : (22* N) array 
    time  : day/h/week un 

"""
import numpy as np
import torch
from  tqdm import tqdm
from script.db import query_pango_data
import pandas as pd
import pangolin
AMINO_ACID =['M', 'F', 'V', 'L', 'P', 'S', 'Q', 'C', 'N', 'T', 'R', 'A', 'Y', 'G', 'D', 'K', 'H', 'W', 'I', 'E', 'X','-','*']
AminoAcid = 23

def poosition_entory(p_array):
    #error = e-7
    #p = torch.tensor(p_array)
    indx = torch.nonzero(p_array>0,as_tuple=False)
    p = torch.index_select(p_array,index=indx.squeeze(),dim=0)
    #22
    H_sum = -(p*torch.log2(p)).sum()
    return H_sum

def computer_rate(H:torch.Tensor,T):
    #H.sum()
    Kaa = torch.mean(H)/(2*abs(T))

    return Kaa
def rate(a:torch.Tensor,T:int):
    '''

    :param a: torch ,size(22,n)
    :param T:
    :return:
    '''
    h_list = []
    for n in range(a.size(1)):
        # print(a[:,n])
        # h = torch.cat(h,)
        h_list.append(poosition_entory(a[:, n]))
    h = torch.tensor(h_list)

    rate_entroy = computer_rate(h, T)

    return rate_entroy

# computer done write csv

#def counts_aa()

def computer_all_rate(pango,N_parent,M_parent,S_parent,E_parent,parent_total,T,area:list):
    #a = torch.full([22,76],1/22)
    
    data = query_pango_data(pango)


    #N, M, S, E
    #gap suan , X bussuan

    AminoAcid = 23 # 20 + X + gap
    N = np.zeros((AminoAcid,420))
    M = np.zeros((AminoAcid,223))
    S = np.zeros((AminoAcid,1274))
    E = np.zeros((AminoAcid,76))

    total = len(data)
    for sigal_data in data:
        falg = False

        for a in area:
            if str(a) in str(sigal_data[4]):
                falg=True
        if not falg:
            continue

        i = 0
        for nn in sigal_data[0]:
            index_aa = AMINO_ACID.index(nn)
            N[index_aa,i] = N[index_aa,i] + 1
            i = i + 1
        i = 0
        for nn in sigal_data[1]:
            index_aa = AMINO_ACID.index(nn)
            M[index_aa,i] = M[index_aa,i] + 1
            i = i + 1

        i = 0
        for nn in sigal_data[2]:
            index_aa = AMINO_ACID.index(nn)
            S[index_aa,i] = S[index_aa,i] + 1
            i = i + 1
        i = 0
        for nn in sigal_data[3]:
            index_aa = AMINO_ACID.index(nn)
            E[index_aa,i] = E[index_aa,i] + 1
            i = i + 1
    #update_alias = pangolin.decompress(Lineage)
    # parent_data = query_pango_data(pango_parent)
    # N_parent = np.zeros((AminoAcid, 420))
    # M_parent = np.zeros((AminoAcid,223))
    # S_parent = np.zeros((AminoAcid,1274))
    # E_parent = np.zeros((AminoAcid,76))
    # parent_total = len(parent_data)
    # for sigal_data in parent_data:
    #     i = 0
    #     for nn in sigal_data[0]:
    #         index_aa = AMINO_ACID.index(nn)
    #         N_parent[index_aa,i] = N_parent[index_aa,i] + 1
    #         i = i + 1
    #     i = 0
    #     for nn in sigal_data[1]:
    #         index_aa = AMINO_ACID.index(nn)
    #         M_parent[index_aa,i] = M_parent[index_aa,i] + 1
    #         i = i + 1
    #     i = 0
    #     for nn in sigal_data[2]:
    #         index_aa = AMINO_ACID.index(nn)
    #         S_parent[index_aa,i] = S_parent[index_aa,i] + 1
    #         i = i + 1
    #     i = 0
    #     for nn in sigal_data[3]:
    #         index_aa = AMINO_ACID.index(nn)
    #         E_parent[index_aa,i] = E_parent[index_aa,i] + 1
    #         i = i + 1
    N_abs_input = abs((N/total)-(N_parent/parent_total))

    M_abs_input = abs((M / total) - (M_parent / parent_total))
    S_abs_input = abs((S / total) - (S_parent / parent_total))
    E_abs_input = abs((E / total) - (E_parent / parent_total))
    #print(abs_input)
    N_rate = rate(torch.tensor(N_abs_input),T)
    M_rate = rate(torch.tensor(M_abs_input), T)
    S_rate = rate(torch.tensor(S_abs_input), T)
    E_rate = rate(torch.tensor(E_abs_input), T)
    #rate_entroy = rate(a,20)
    del data
    return N_rate,M_rate,S_rate,E_rate
    #print(N_rate,M_rate,S_rate,E_rate)

# if __name__ == '__main__':
#
#
#     computer()
    # data = pd.read_csv('../data/findl_gap_lineage.csv')
    #
    # val_list = []
    # for i in tqdm(range(len(data))):
    #     T = data.iloc[i]['gap']
    #     if T ==0:
    #         continue
    #     N_rate, M_rate, S_rate, E_rate = computer_all_rate(data.iloc[i]['Lineage'],data.iloc[i]['parent'],T)
    #     #tensor.item
    #     val_list.append([data.iloc[i]['Lineage'],pangolin.compress(data.iloc[i]['parent']),data.iloc[i]['MIN(collection_date)'],data.iloc[i]['parent_date'],T,'{:.10f}'.format(N_rate), '{:.10f}'.format(M_rate), '{:.10f}'.format(S_rate),'{:.10f}'.format( E_rate)])
    #
    # val = np.array(val_list)
    # new_data = pd.DataFrame(val,columns=['Lineage','parent','MIN(collection_date)','parent_date','gap','N','M','S','E'])
    # new_data.to_csv('../data/final_rate.csv')
    # print('ok')
def get_parent(pango_parent,area):
    parent_data = query_pango_data(pango_parent)
    N_parent = np.zeros((AminoAcid, 420))
    M_parent = np.zeros((AminoAcid,223))
    S_parent = np.zeros((AminoAcid,1274))
    E_parent = np.zeros((AminoAcid,76))
    parent_total = len(parent_data)
    for sigal_data in parent_data:
        falg = False
        for a in area:
            if str(a) in str(sigal_data[4]):
                falg=True
        if not falg:

            continue


        i = 0
        for nn in sigal_data[0]:
            index_aa = AMINO_ACID.index(nn)
            N_parent[index_aa,i] = N_parent[index_aa,i] + 1
            i = i + 1
        i = 0
        for nn in sigal_data[1]:
            index_aa = AMINO_ACID.index(nn)
            M_parent[index_aa,i] = M_parent[index_aa,i] + 1
            i = i + 1
        i = 0
        for nn in sigal_data[2]:
            index_aa = AMINO_ACID.index(nn)
            S_parent[index_aa,i] = S_parent[index_aa,i] + 1
            i = i + 1
        i = 0
        for nn in sigal_data[3]:
            index_aa = AMINO_ACID.index(nn)
            E_parent[index_aa,i] = E_parent[index_aa,i] + 1
            i = i + 1
    del parent_data
    return N_parent,M_parent,S_parent,E_parent,parent_total
def computer_location():
    Area = ['Europe','Africa','North America','Oceania','Asia','South America','China','United Kingdom','India','Japan','Australia']

    data = pd.read_csv('../data/findl_gap_lineage.csv')
    for area in Area:

        val_list = []
        parent_list = set(data['parent'].to_list())
        f = 0
        for parent in tqdm(parent_list):
            # f = f + 1
            # if f<110:
            #     continue

            N_parent,M_parent,S_parent,E_parent,parent_total = get_parent(parent,area)
            current_data = data[data['parent']==parent]
            for i in range(len(current_data)):
                T = current_data.iloc[i]['gap']
                if T ==0:
                    continue

                N_rate, M_rate, S_rate, E_rate = computer_all_rate(current_data.iloc[i]['Lineage'],N_parent,M_parent,S_parent,E_parent,parent_total,T,[area])
                #tensor.item
                # except:
                #     print(current_data.iloc[i]['Lineage'],parent)
                val_list.append([current_data.iloc[i]['Lineage'],pangolin.compress(current_data.iloc[i]['parent']),current_data.iloc[i]['MIN(collection_date)'],current_data.iloc[i]['parent_date'],T,'{:.10f}'.format(N_rate), '{:.10f}'.format(M_rate), '{:.10f}'.format(S_rate),'{:.10f}'.format( E_rate)])

        val = np.array(val_list)
        new_data = pd.DataFrame(val,columns=['Lineage','parent','MIN(collection_date)','parent_date','gap','N','M','S','E'])
        new_data.to_csv(f'../data/final_rate_{area}.csv')
        print(f'-----{area}')
    #East_Asia = ['China、'Japan'、'South Korea'、朝鲜,'Mongolia]
    #North_Asia = [不丹、尼泊尔、印度、巴基斯坦、孟加拉国、斯里兰卡，马尔代]
    #Country = ['China','United Kingdom','India','Japan','Australia']
    pass



if __name__ == '__main__':
        computer_location()