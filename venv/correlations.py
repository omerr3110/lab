#%%
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy as sp
from scipy.stats import ranksums
#%%
#1

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

def vector_creator(db,prog,background):
    base=db.join(full_positions)
    base=base.dropna()
    base=keep_first(base)
    if prog=="t":
        vector1 = base['T_position']
        vector1 = vector1.to_numpy()
        vector2 = background['T_position']
        vector2 = vector2.to_numpy()
    if prog=="r":
        vector1 = base['R_position']
        vector1 = vector1.to_numpy()
        vector2 = background['R_position']
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

def df_creator(df,prot):
    res=df[df['Vrial gene ']== prot]
    return res


def rep_rank(lists):
    dict={}
    for list in lists:
        for prot in list:
            if prot in dict:
                dict[prot]=dict[prot]+1
            else:
                dict[prot]=1
    return dict

def strain_repair(name):
    if pd.isna(name):
        return "PR8"
    return name

    #%%
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
    translation=pd.read_csv(r'C:\Users\omerr\PycharmProjects\lab\venv\HMD_HumanPhenotype.rpt.txt',delimiter="\t")
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
#3
    # %%
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
    '''corr_back = all_pr8.join(main_curr)
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
    pearson_calc(ns2_prot_corr)'''
#%%
    #enrichment analysis - pr8 db
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
    background2=keep_first(background2)
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
    #%%
    #enrichment analysis - other strains db
    all_other = pd.read_excel(r"C:\Users\omerr\PycharmProjects\lab\venv\other strains complete.xlsx")
    all_other = all_other.rename(columns=all_other.iloc[1])
    all_other = all_other.drop([0, 1, 2])
    all_other = all_other.rename({"Spectral count": "sc rep1", np.nan: "sc rep2"}, axis=1)
    all_other['Gene'] = all_other['Gene'].apply(translate)
    all_other = all_other.set_index('Gene')
    background3 = all_other.join(full_positions)
    background3 = background3.dropna()
    aichi_all = background3[background3['Strain'] == "Aichi"]
    single_aichi = keep_first(aichi_all)
    NY2009_all = background3[background3['Strain'] == "NY/2009"]
    single_NY2009 = keep_first(NY2009_all)
    WSN33_all = background3[background3['Strain'] == "WSN/33"]
    single_WSN33 = keep_first(WSN33_all)
    H5N1_all = background3[background3['Strain'] == "H5N1"]
    single_H5N1 = keep_first(H5N1_all)
    other_inter= pd.read_excel(r'C:\Users\omerr\PycharmProjects\lab\venv\other strains.xlsx')
    other_inter=other_inter.dropna()
    header=other_inter.iloc[0]
    other_inter.columns=header
    other_inter=other_inter.drop(other_inter.index[:1])
    other_inter['Gene ID'] = other_inter['Gene ID'].apply(translate)
    other_inter = other_inter.set_index('Gene ID')
    other_inter['Vrial gene '] = other_inter['Vrial gene '].apply(NaN_fix)
    aichi_inter=other_inter[other_inter['Strain'] == "Aichi"]
    NY2009_inter = other_inter[other_inter['Strain'] == "NY/2009"]
    WSN33_inter = other_inter[other_inter['Strain'] == "WSN/33"]
    H5N1_inter = other_inter[other_inter['Strain'] == "H5N1"]
    #%%
    #5
    prot_check=df_creator(aichi_inter,"NS1")
    vectors=vector_creator(prot_check,"r",single_aichi)
    results=ranksums(vectors[0],vectors[1],alternative="less")
    pval=results[1]
    statval=results[0]
    print(pval,statval)
    results = ranksums(vectors[0], vectors[1])
    pval = results[1]
    statval = results[0]
    print(pval, statval)
    results = ranksums(vectors[0], vectors[1], alternative="greater")
    pval = results[1]
    statval = results[0]
    print(pval, statval)
    log_pval=math.log(pval)*-1
    act=vectors[2]
    rep = act.index.duplicated(keep='first')
    not_rep = ~rep
    act = act[not_rep]
    ax2=sns.scatterplot(data=single_aichi, x='T_position', y='R_position', color='gray')
    sns.scatterplot(data=act, x='T_position', y='R_position', color="red")
    ax2.set(title='Aichi - NS1')
    plt.show()
    #plt.scatter(data=background2, x='R_position', y='T_position',c='SAINT score',cmap="hot")
    # %%
    #6
    np_inter = other_inter[other_inter['Vrial gene '] == "NP"]
    np_total = np_inter.append(np_prot)
    np_single = keep_first(np_total)
    np_single.columns.name = None
    count = np_total.index.value_counts(sort=False)
    np_single = np_single.assign(count=count)
    np_single['Strain'] = np_single['Strain'].apply(strain_repair)

    m2_inter = other_inter[other_inter['Vrial gene '] == "M2"]
    m2_total = m2_inter.append(m2_prot)
    m2_single = keep_first(m2_total)
    m2_single.columns.name = None
    count = m2_total.index.value_counts(sort=False)
    m2_single = m2_single.assign(count=count)
    m2_single['Strain'] = m2_single['Strain'].apply(strain_repair)
    m2_single.loc['Gm11127','count']=4

    ns1_inter = other_inter[other_inter['Vrial gene '] == "NS1"]
    ns1_total = ns1_inter.append(ns1_prot)
    ns1_single = keep_first(ns1_total)
    ns1_single.columns.name = None
    count = ns1_total.index.value_counts(sort=False)
    ns1_single = ns1_single.assign(count=count)
    ns1_single['Strain'] = ns1_single['Strain'].apply(strain_repair)

    pb1_inter = other_inter[other_inter['Vrial gene '] == "PB1"]
    pb1_total = pb1_inter.append(pb1_prot)
    pb1_single = keep_first(pb1_total)
    pb1_single.columns.name = None
    count = pb1_total.index.value_counts(sort=False)
    pb1_single = pb1_single.assign(count=count)
    pb1_single['Strain'] = pb1_single['Strain'].apply(strain_repair)

    pb2_inter = other_inter[other_inter['Vrial gene '] == "PB2"]
    pb2_total = pb2_inter.append(pb2_prot)
    pb2_single = keep_first(pb2_total)
    pb2_single.columns.name = None
    count = pb2_total.index.value_counts(sort=False)
    pb2_single = pb2_single.assign(count=count)
    pb2_single['Strain'] = pb2_single['Strain'].apply(strain_repair)

    full_background=background3.append(background2)
    full_background['Strain']=full_background['Strain'].apply(strain_repair)
    full_single_background = keep_first(full_background)
    vectors=vector_creator(np_single,'t',full_single_background)
    col_dict={0:"mistyrose",1:"rosybrown",2:"salmon",3:"indianred",4:"firebrick",5:"darkred"}
    a=sns.scatterplot(data=vectors[2], x='T_position', y='R_position',hue='count',palette=col_dict,size='count')
    a.axis('equal')
    a.set(title='np interactions count')

    plt.show()

    #combine p values






















