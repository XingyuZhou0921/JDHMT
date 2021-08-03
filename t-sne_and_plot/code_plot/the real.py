import matplotlib.pyplot as plt  
import seaborn as sns
import pandas as pd
import numpy as np
tips = pd.read_csv('out_DG_out1.txt.csv')
area = np.pi * 4
fig = sns.scatterplot( x="c1", y="c2",data=tips[:86],alpha = 0.1,s = area,color = 'yellowgreen')
scatter_fig = fig.get_figure()
scatter_fig.savefig('test.png', dpi = 400)
plt.show()



