import math

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy as sp
from scipy.stats import ranksums



def grade_calc(ser):
    grade=ser.sum()/12
    return grade


def virus_prot_grade_calc(base_file,trf,translator):
    rows=base_file.shape[0]
    dict={'PA':[0,0,0],'M1':[0,0,0],'M2':[0,0,0],'HA':[0,0,0],'NA':[0,0,0],'NP':[0,0,0],'NS1':[0,0,0],'NS2':[0,0,0]
        ,'PB1':[0,0,0],'PB2':[0,0,0]}
    problem=[]
    flag=0
    for i in range(rows):
        flag=0
        row=base_file.iloc[i]
        name=row.loc['Host symbol']
        if name in translator.keys():
            homologues=translator[name]
            for protein in homologues:
                if protein in trf.index:
                    flag=1
                    connections = row.loc['PR8 viral protein']
                    neighbors=connections.split(",")
                    for gene in neighbors:
                        r_grade = trf.loc[protein,'average_r_grade']
                        t_grade = trf.loc[protein, 'average_t_grade']
                        dict[gene][0] = dict[gene][0] + r_grade
                        dict[gene][1] = dict[gene][1] + t_grade
                        dict[gene][2] = dict[gene][2]+1
            if flag==0:
                problem.append(name)
        else:
            problem.append(name)
    for key in dict.keys():
        if dict[key][2]!=0:
            dict[key][0] = dict[key][0] / dict[key][2]
            dict[key][1] = dict[key][1] / dict[key][2]
    return dict


def create_dict(file):
    n=file.shape[0]
    human_to_mouse={}
    mouse_to_human={}
    for i in range(n):
        human=file.iloc[i,0]
        mouse=file.iloc[i,2]
        if human in human_to_mouse.keys():
            human_to_mouse[human].append(mouse)
        else:
            human_to_mouse[human]=[mouse]
        if mouse in mouse_to_human.keys():
            mouse_to_human[mouse].append(human)
        else:
            mouse_to_human[mouse]=[human]
    return (human_to_mouse, mouse_to_human)

def translate(gene):
    return translate_full(gene,full_corr,trans_dict)

def translate_full(gene,file,dictt):
    name=gene
    if gene in dictt.keys():
        for option in dictt[gene]:
            if option in file.index:
                name=option
                break
    elif type(name)==type("str"):
        name=name.lower()
        name=name.capitalize()
    if gene=="SKP1":
        name="Skp1a"
    return name

def vector_creator(db,prog):
    base=db.join(full_positions)
    base=base.dropna()
    if prog=="t":
        vector1 = base['T_position']
        vector1 = vector1.to_numpy()
        vector2 = background2['T_position']
        vector2 = vector2.to_numpy()
    if prog=="r":
        vector1 = base['R_position']
        vector1 = vector1.to_numpy()
        vector2 = background2['R_position']
        vector2 = vector2.to_numpy()
    return (vector1,vector2,base)

def NaN_fix(prot):
    if pd.isna(prot):
        return "Na"
    else:
        return prot

def keep_first(db):
    rep = db.index.duplicated(keep='first')
    not_rep = ~rep
    final = db[not_rep]
    return final

def pearson_calc(db):
    db1=no_infection=db[db['IAV infection']=="with"]
    db2=no_infection=db[db['IAV infection']=="without"]
    db1 = keep_first(db1)
    db2 = keep_first(db2)
    pearson_wt = sp.stats.pearsonr(db1['SAINT score'], db1['T_correlation_CC'])
    pearson_wr = sp.stats.pearsonr(db1['SAINT score'], db1['R_correlation_CC'])
    pearson_nt = sp.stats.pearsonr(db2['SAINT score'], db2['T_correlation_CC'])
    pearson_nr = sp.stats.pearsonr(db2['SAINT score'], db2['R_correlation_CC'])
    return ([pearson_wt,pearson_wr,pearson_nt,pearson_nr])



if __name__ == '__main__':
    #total_file=pd.read_excel(r'C:\Users\omerr\PycharmProjects\lab\venv\R_T_gene_correlation.xlsx')
    #total_file.rename(columns={'Unnamed: 0': 'Gene','Unnamed: 1': 't_pos','Unnamed: 2': 'r_pos'}, inplace=True)
    xls=pd.ExcelFile(r'C:\Users\omerr\PycharmProjects\lab\venv\Cell2009_Shapira_et al_data.xls')
    xls2=pd.ExcelFile(r'C:\Users\omerr\PycharmProjects\lab\venv\unified_gene_states_corr_with_healthy_sep.igv8.xlsx')
    full_corr=pd.read_excel(xls2,'Sheet1')
    full_corr.rename(columns={'Unnamed: 0': 'Gene'},inplace=True)
    full_corr=full_corr.set_index('Gene')
    full_positions=full_corr.iloc[:,:2]
    full_positions=full_positions.dropna()
    main_curr=full_corr.iloc[:,2:4]
    main_curr=main_curr.dropna()
    interactions=pd.read_excel(xls,'4.viral-host PR8')
    bank= pd.read_excel(xls,'1.Gene expression')
    #corr=total_file.dropna()
    #corr = corr.set_index('Gene')
    #R_file=corr.iloc[:,2:14]
    #T_file=corr.iloc[:,14:26]
    #positions=corr.iloc[:,:2]
    #num=R_file.shape[0]
    #t_grades=[grade_calc(T_file.iloc[i]) for i in range(num)]
    #r_grades=[grade_calc(R_file.iloc[i]) for i in range(num)]
    #corr['average_t_grade']=t_grades
    #corr['average_r_grade']=r_grades
    translation=pd.read_csv(r'HMD_HumanPhenotype.rpt.txt',delimiter="\t")
    translation.columns=['1','2','3','4','5','6']
    translation=translation.dropna(axis=1)
    dicts=create_dict(translation)
    trans_dict=dicts[0]
    #result=virus_prot_grade_calc(interactions,corr,trans_dict)
    bank['symbol']=bank['symbol'].apply(translate)
    bank=bank.set_index('symbol')
    background=bank.join(full_positions)
    background=background.dropna()
    #ax1=background.plot.scatter(x='r_pos',y='t_pos',c='gray')
    '''interactions['Host symbol'] = interactions['Host symbol'].apply(translate)
    interactions=interactions.set_index('Host symbol')
    active=interactions.join(corr)
    active = active.dropna()
    #active.plot.scatter(x='r_pos',y='t_pos',c='red',ax=ax1)
    sns.set_theme()
    graph1=sns.relplot(data=background,x='r_pos',y='t_pos',color="grey")
    sns.scatterplot(data=active,x='r_pos',y='t_pos',color="red")'''


    #second_db
    all_pr8 = pd.read_excel(r"C:\Users\omerr\PycharmProjects\lab\venv\PR8 complete.xlsx")
    all_pr8 = all_pr8.rename(columns=all_pr8.iloc[1])
    all_pr8 = all_pr8.drop([0, 1, 2])
    all_pr8 = all_pr8.rename({"Spectral count": "sc rep1", np.nan: "sc rep2"}, axis=1)
    all_pr8['Gene'] = all_pr8['Gene'].apply(translate)
    all_pr8 = all_pr8.set_index('Gene')
    is_dup = all_pr8.index.duplicated(keep='first')
    not_dup = ~is_dup
    single_pr8 = all_pr8[not_dup]
    background2 = single_pr8.join(full_positions)
    background2 = background2.dropna()
    #saint correlations
    corr_back = all_pr8.join(main_curr)
    corr_back['Vrial gene ']=corr_back['Vrial gene '].apply(NaN_fix)
    corr_back = corr_back.dropna()
    ha_prot_corr = corr_back[corr_back['Vrial gene '] == "HA"]
    m1_prot_corr = corr_back[corr_back['Vrial gene '] == "M1"]
    m2_prot_corr = corr_back[corr_back['Vrial gene '] == "M2"]
    pb1_prot_corr = corr_back[corr_back['Vrial gene '] == "PB1"]
    pb1f2_prot_corr = corr_back[corr_back['Vrial gene '] == "PB1F2"]
    pa_prot_corr = corr_back[corr_back['Vrial gene '] == "PA"]
    pb2_prot_corr = corr_back[corr_back['Vrial gene '] == "PB2"]
    na_prot_corr = corr_back[corr_back['Vrial gene '] == "Na"]
    np_prot_corr = corr_back[corr_back['Vrial gene '] == "NP"]
    ns1_prot_corr = corr_back[corr_back['Vrial gene '] == "NS1"]
    ns2_prot_corr = corr_back[corr_back['Vrial gene '] == "NS2"]
    pearson_calc(ns2_prot_corr)
    #enrichment analysis - pr8 db
    inf_inter = pd.read_excel(r'C:\Users\omerr\PycharmProjects\lab\venv\PR8 SAINT TOP.xlsx')
    inf_inter = inf_inter.rename(columns=inf_inter.iloc[1])
    inf_inter = inf_inter.drop([0, 1])
    inf_inter['Gene ID']=inf_inter['Gene ID'].apply(translate)
    inf_inter = inf_inter.set_index('Gene ID')
    inf_inter['Vrial gene ']=inf_inter['Vrial gene '].apply(NaN_fix)
    ha_prot=inf_inter[inf_inter['Vrial gene ']=="HA"]
    m1_prot = inf_inter[inf_inter['Vrial gene '] == "M1"]
    m2_prot = inf_inter[inf_inter['Vrial gene '] == "M2"]
    pb1_prot = inf_inter[inf_inter['Vrial gene '] == "PB1"]
    pb1f2_prot = inf_inter[inf_inter['Vrial gene '] == "PB1F2"]
    pa_prot = inf_inter[inf_inter['Vrial gene '] == "PA"]
    pb2_prot = inf_inter[inf_inter['Vrial gene '] == "PB2"]
    na_prot = inf_inter[inf_inter['Vrial gene '] == "Na"]
    np_prot = inf_inter[inf_inter['Vrial gene '] == "NP"]
    ns1_prot = inf_inter[inf_inter['Vrial gene '] == "NS1"]
    ns2_prot = inf_inter[inf_inter['Vrial gene '] == "NS2"]

    #enrichment analysis - other strains db
    other_inter= pd.read_excel(r'C:\Users\omerr\PycharmProjects\lab\venv\PR8 SAINT TOP.xlsx')





    vectors=vector_creator(ns2_prot,"r")
    results=ranksums(vectors[0],vectors[1],alternative="less")
    pval=results[1]
    statval=results[0]
    results = ranksums(vectors[0], vectors[1])
    pval = results[1]
    statval = results[0]
    results = ranksums(vectors[0], vectors[1], alternative="greater")
    pval = results[1]
    statval = results[0]
    log_pval=math.log(pval)*-1
    act=vectors[2]
    rep = act.index.duplicated(keep='first')
    not_rep = ~rep
    act = act[not_rep]
    ax2=sns.relplot(data=background2, x='R_position', y='T_position', color='blue')
    sns.scatterplot(data=act, x='R_position', y='T_position', color="red")
    plt.scatter(data=background2, x='R_position', y='T_position',c='SAINT score',cmap="hot")
    plt.show()












