import numpy as np
import pandas as pd
import matplotlib.pyplot as plt




if __name__ == '__main__':
    all_pr8 = pd.read_excel(r"C:\Users\omerr\PycharmProjects\lab\venv\PR8 complete.xlsx")
    all_pr8=all_pr8.rename(columns=all_pr8.iloc[1])
    all_pr8=all_pr8.drop([0,1,2])
    all_pr8=all_pr8.rename({"Spectral count":"sc rep1",np.nan:"sc rep2"},axis=1)
    print(all_pr8)


