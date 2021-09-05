import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


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
    return translate_full(gene,file1,trans_dict)

def translate_full(gene,file,dictt):
    name=gene
    if gene in dictt.keys():
        for option in dictt[gene]:
            if option in file.index:
                name=option
    return name

if __name__ == '__main__':
    total_file=pd.read_excel(r'C:\Users\omerr\PycharmProjects\lab\venv\R_T_gene_correlation.xlsx')
    total_file.rename(columns={'Unnamed: 0': 'Gene','Unnamed: 1': 't_pos','Unnamed: 2': 'r_pos'}, inplace=True)
    xls=pd.ExcelFile(r'C:\Users\omerr\PycharmProjects\lab\venv\Cell2009_Shapira_et al_data.xls')
    interactions=pd.read_excel(xls,'4.viral-host PR8')
    bank= pd.read_excel(xls,'1.Gene expression')
    file1=total_file.dropna()
    file1 = file1.set_index('Gene')
    R_file=file1.iloc[:,2:14]
    T_file=file1.iloc[:,14:26]
    positions=file1.iloc[:,:2]
    num=R_file.shape[0]
    t_grades=[grade_calc(T_file.iloc[i]) for i in range(num)]
    r_grades=[grade_calc(R_file.iloc[i]) for i in range(num)]
    file1['average_t_grade']=t_grades
    file1['average_r_grade']=r_grades
    translation=pd.read_csv(r'HMD_HumanPhenotype.rpt.txt',delimiter="\t")
    translation.columns=['1','2','3','4','5','6']
    translation=translation.dropna(axis=1)
    dicts=create_dict(translation)
    trans_dict=dicts[0]
    #result=virus_prot_grade_calc(interactions,file1,trans_dict)
    bank['symbol']=bank['symbol'].apply(translate)
    bank=bank.set_index('symbol')
    background=bank.join(file1)
    background=background.dropna()
    ax1=background.plot.scatter(x='r_pos',y='t_pos',c='gray')
    interactions['Host symbol'] = interactions['Host symbol'].apply(translate)
    interactions=interactions.set_index('Host symbol')
    active=interactions.join(file1)
    active = active.dropna()
    active.plot.scatter(x='r_pos',y='t_pos',c='red',ax=ax1)
    all_pr8 = pd.read_excel(r"C:\Users\omerr\PycharmProjects\lab\venv\PR8 complete.xlsx")
    all_pr8 = all_pr8.rename(columns=all_pr8.iloc[1])
    all_pr8 = all_pr8.drop([0, 1, 2])
    all_pr8 = all_pr8.rename({"Spectral count": "sc rep1", np.nan: "sc rep2"}, axis=1)
    all_pr8['Gene'] = all_pr8['Gene'].apply(translate)
    all_pr8 = all_pr8.set_index('Gene')
    is_dup = all_pr8.index.duplicated(keep='first')
    not_dup = ~is_dup
    single_pr8 = all_pr8[not_dup]
    background2 = single_pr8.join(file1)
    background2 = background2.dropna()
    ax2 = background2.plot.scatter(x='r_pos', y='t_pos', c='gray')
    inf_inter = pd.read_excel(r'C:\Users\omerr\PycharmProjects\lab\venv\PR8 SAINT TOP.xlsx')
    inf_inter = inf_inter.rename(columns=inf_inter.iloc[1])
    inf_inter = inf_inter.drop([0, 1, 2])
    inf_inter['Gene ID']=inf_inter['Gene ID'].apply(translate)
    inf_inter = inf_inter.set_index('Gene ID')
    is_dup2 = inf_inter.index.duplicated(keep='first')
    not_dup2 = ~is_dup2
    single_inf_inter = inf_inter[not_dup2]
    act=inf_inter.join(file1)
    act=act.dropna()
    act.plot.scatter(x='r_pos',y='t_pos',c='red',ax=ax2)



    plt.show()













