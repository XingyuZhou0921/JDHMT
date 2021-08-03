#1.#C9DAA7 and #F9BE8E 2.#BED296 3.#AEC87B
#调色后:1.#AEC87B & #F9B984  2.#AEC87B & #f97306 3.#ff796c
#4.#F89D54
import pandas as pd
import matplotlib.pyplot as plt
data = pd.read_csv('CE.csv')
gene_data=data[data['label']==1]# 绿色
disease_data=data[data['label']==0]# 橙色

f,ax=plt.subplots()
ax.scatter(gene_data['c1'],gene_data['c2'],s=15,marker='.',facecolors='none',edgecolors='#AEC87B',linewidths=0.5)
#ax.scatter(disease_data['c1'],disease_data['c2'],s=15,marker='.',facecolors='none',edgecolors='#F89D54',linewidths=0.5)
#plt.show()
scatter_fig = ax.get_figure()
scatter_fig.savefig('CE_gene.png', dpi = 400)
#plt.show()
