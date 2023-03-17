"""
@Author:weidong wu
@Author_Email:weidongwu404@gmail.com
@Time:2022/10/12 下午1:51
"""

import pymysql
import numpy as np
conn = pymysql.connect(
    host='127.0.0.1',
    user='sars',
    password='sarspassword',
    database='sars',
    charset='utf8',
    # autocommit=True,    # 如果插入数据，， 是否自动提交? 和conn.commit()功能一致。
)

def query_pango_data(pango:str):

    # cursor = conn.cursor()
    # sql = 'select N,M,S,E,location from NSME_final where pango_lineage="{}";'.format(pango)
    # #print(sql)
    # cursor.execute(sql)
    # rest = cursor.fetchall()
    try:

        rest = np.load(f'/media/wvdon/sdata/covid_evolution/npy/{pango}.npy')
    except:
        print('errorp,',pango)
    return rest

def test():

    data = query_pango_data('B.1.142')

    print(len(data))