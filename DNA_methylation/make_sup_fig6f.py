import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl

df = pd.read_table("~/haplotype_methylation_matrix.tsv", index_col=0)
#mark 2 to "None"
df = df.replace("None", 2).astype(float)

#devide data frame between HP1 and HP2
even = []
odd = []
for i in list(range(len(df.columns))):
    if i % 2 == 0:
        even.append(i)
    else:
        odd.append(i)

df_HP1 = df.iloc[:,even]
for i in list(df_HP1.columns):
    i_new = i.split("_")[0]
    df_HP1 = df_HP1.rename(columns={i : i_new})

df_HP2 = df.iloc[:,odd]
for i in list(df_HP2.columns):
    i_new = i.split("_")[0]
    df_HP2 = df_HP2.rename(columns={i : i_new})

#get sample names for x axis label and gene names for y axis label
x_sample_list = list(df_HP1.columns)
y_gene_list = list(df_HP1.index)
y_gene_list_rev = y_gene_list[::-1]
#make list of label location 
x_label_loc = []
a = 0.5
for i in range(60):
    x_label_loc.append(a)
    a += 1

y_label_loc = []
b = 0.5
for i in range(29):
    y_label_loc.append(b)
    b += 1


#change pandas object to numpy array, then prepare mask array
mat_HP1 = df_HP1.to_numpy()
mat_HP1_flipud = np.flipud(mat_HP1)
mat_HP1_val = mat_HP1_flipud.ravel()
Mask_HP1 = np.zeros(np.shape(mat_HP1_flipud))
for i in range(29):
    for j in range(60):
        if(mat_HP1_flipud[i][j] == 2.0):
            Mask_HP1[i][j] = True
        else:
            Mask_HP1[i][j] = False
Mask_HP1 = Mask_HP1.ravel()
            
mat_HP2 = df_HP2.to_numpy()
mat_HP2_flipud = np.flipud(mat_HP2)
mat_HP2_val = mat_HP2_flipud.ravel()
Mask_HP2 = np.zeros(np.shape(mat_HP2_flipud))
for i in range(29):
    for j in range(60):
        if(mat_HP2_flipud[i][j] == 2.0):
            Mask_HP2[i][j] = True
        else:
            Mask_HP2[i][j] = False
Mask_HP2 = Mask_HP2.ravel()


list_NA_HP1 = []
for i in range(len(Mask_HP1)):
    if Mask_HP1[i] == 1:
        list_NA_HP1.append(i)
list_NA_HP2 = []
for i in range(len(Mask_HP2)):
    if Mask_HP2[i] == 1:
        list_NA_HP2.append(i)
        

#make figure
M = 60
N = 29
x = np.arange(M + 1)
y = np.arange(N + 1)
xs, ys = np.meshgrid(x, y)
triangles1 = [(i + j*(M+1), i + (j+1)*(M+1), i+1 + (j+1)*(M+1)) for j in range(N) for i in range(M)]
triangles2 = [(i + j*(M+1), i+1 + j*(M+1), i+1 + (j+1)*(M+1)) for j in range(N) for i in range(M)]

#get location of "None" in HP1
triangles3 = []
for i in list_NA_HP1:
    triangles3.append(triangles1[i])

mat_HP1_grey = np.empty(len(list_NA_HP1))
for i in range(len(list_NA_HP1)):
    mat_HP1_grey[i] = 0.8

#get location of "None" in HP2
triangles4 = []
for i in list_NA_HP2:
    triangles4.append(triangles2[i])

mat_HP2_grey = np.empty(len(list_NA_HP2))
for i in range(len(list_NA_HP2)):
    mat_HP2_grey[i] = 0.8


#make color map
cdict = {'red':   [(0.0,  0.0, 0.0),
                   (0.5,  1.0, 1.0),
                   (1.0,  1.0, 1.0)],

         'green': [(0.0,  0.4375, 0.4375),
                   (0.5,  1.0, 1.0),
                   (1.0,  0.0, 0.0)],

         'blue':  [(0.0,  0.75, 0.75),
                   (0.5,  1.0, 1.0),
                   (1.0,  0.0, 0.0)]}
newcmp = LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)



triang1 = Triangulation(xs.ravel(), ys.ravel(), triangles1, mask=Mask_HP1)
triang2 = Triangulation(xs.ravel(), ys.ravel(), triangles2, mask=Mask_HP2)
triang3 = Triangulation(xs.ravel(), ys.ravel(), triangles3)
triang4 = Triangulation(xs.ravel(), ys.ravel(), triangles4)
img1 = plt.tripcolor(triang1, mat_HP1_val, cmap=newcmp, vmax=1, vmin=0, edgecolors='0.3', linewidth=0.7)
img2 = plt.tripcolor(triang2, mat_HP2_val, cmap=newcmp, vmax=1, vmin=0, edgecolors='0.3', linewidth=0.7)
img3 = plt.tripcolor(triang3, mat_HP1_grey, cmap=plt.get_cmap('gray'), vmax=1, vmin=0, edgecolors='0.3', linewidth=0.7)
img4 = plt.tripcolor(triang4, mat_HP2_grey, cmap=plt.get_cmap('gray'), vmax=1, vmin=0, edgecolors='0.3', linewidth=0.7)
for i in range(29):
    plt.axhline(y=i, color="white", linewidth=1.2)
for i in range(60):
    plt.axvline(x=i, color="white", linewidth=1.2)
plt.colorbar(img1, ticks=None)
plt.xlim(x[0], x[-1])
plt.ylim(y[0], y[-1])
plt.xticks(x_label_loc, x_sample_list, rotation=90)
plt.yticks(y_label_loc, y_gene_list_rev)

plt.show()