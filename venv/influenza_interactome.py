import numpy as np
import pandas as pd
import matplotlib.pyplot as plt




if __name__ == '__main__':
    all_pr8 = pd.read_excel(r"C:\Users\omerr\PycharmProjects\lab\venv\PR8 complete.xlsx")
    all_pr8 = all_pr8.rename(columns=all_pr8.iloc[1])
    all_pr8 = all_pr8.drop([0, 1, 2])
    all_pr8 = all_pr8.rename({"Spectral count": "sc rep1", np.nan: "sc rep2"}, axis=1)
    all_pr8=all_pr8.set_index('Gene')
    is_dup=all_pr8.index.duplicated(keep='first')
    not_dup=~is_dup
    single_pr8=all_pr8[not_dup]
    inf_inter=pd.read_excel(r'C:\Users\omerr\PycharmProjects\lab\venv\PR8 SAINT TOP.xlsx')
    inf_inter = inf_inter.rename(columns=inf_inter.iloc[1])
    inf_inter = inf_inter.drop([0, 1, 2])
    inf_inter = inf_inter.set_index(['Gene ID'])
    is_dup2 = inf_inter.index.duplicated(keep='first')
    not_dup2 = ~is_dup2
    single_inf_inter = inf_inter[not_dup2]
    print(single_inf_inter)





