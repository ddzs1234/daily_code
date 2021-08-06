#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 11:24:01 2017

@author: ashley
"""

"""
nameg3=[]
numg3=[]
redg3=[]
namel3=[]
numl3=[]
redl3=[]

for i in range(0,len(name3),1):
    if num_3[i]>3:
        nameg3.append(name3[i])
        numg3.append(num_3[i])
        redg3.append(red_3[i])
    if num_3[i]>=2 and num_3[i]<=3:
        namel3.append(name3[i])
        numl3.append(num_3[i])
        redl3.append(red_3[i])

print('len_nameg3',len(nameg3))
print('len_namel3',len(namel3))
"""
'''
nameg3.extend(namel3)
#print('type',type(allname))
a=list(set(name_2).difference(set(nameg3)))
print(a)
print('test:',namel3[0])
b=name2.index('CXO J121749.7+301151')
print('index_b',b)
'''
['CXO J121749.7+301151',
 'CXO J121345.5+364412',
 'CXO J141756.6+254326',
 'CXO J225715.9+074251', 
 'CXO J093037.5+495025', 
 'CXO J124203.5+063829',
 'CXO J105329.4+572104',
 'CXO J103223.5+532319',
 'CXO J121049.6+392821', 
 'CXO J122911.6+020313',
 'CXO J141953.7+541920', 
 'CXO J165258.8+022403',
 'CXO J103140.2+505436', 
 'CXO J141756.8+254430', 
 'CXO J122915.4+020529',
 'CXO J124155.0+063327',
 'CXO J121032.5+392420']
#lost data
"""
#========num<2,3>
indexl3=[]
varl3=[]
varl3_real=[]
#print('single:',name3.index(namel3[0]))
for i in range(0,len(namel3),1):
    fluxtmp=[]
    index_l3=name3.index(namel3[i])
    #index_tmp=[n for n,v in enumerate(name3) if v==namel3[i]]
    #index_l3=np.min(index_tmp)
    #print('index_l3:',index_l3)
    indexl3.append(index_l3)
    a_num=numl3[i]
    for j in range(index_l3,a_num+index_l3,1):
        if flux3[j]==0.0:
            flux_tmp=flux3_hilim[j]
        else:
            flux_tmp=flux3[j]
        fluxtmp.append(flux_tmp)
    #print('fluxtmp:',fluxtmp)
    flux_tmpmin=np.min(fluxtmp)
    flux_tmpmax=np.max(fluxtmp)
    #print
    vary=(flux_tmpmax)/(flux_tmpmin)
    vary_real=np.var(fluxtmp)/(np.mean(flux_tmp)**2)
    #print('vary_real',vary_real)

    #print('vary',vary)
    varl3.append(vary)
    varl3_real.append(vary_real)
#index=[i for i,v in enumerate(varl3) if v=='nan']
#print('index_nan:',index) 
#print('len_varl3',len(varl3))
#print('var:',var)
ymax1=np.max(varl3_real)*1.1
Dll3=[i*3E8*308568E16/67.73999999999998 for i in redl3]
Dl3=[i*3.24077929e-25 for i in Dll3]
plt.figure()
#plt.ylim(0,ymax1)
grey=plt.scatter(Dl3,varl3_real,s=4,facecolors='none',edgecolor='c',label='3>=acis_num>=2')
plt.legend(loc=0)          
plt.xscale('log')
plt.xlabel('Distance'+'  '+'Mpc')
plt.ylabel('var_real')

plt.title('var_real23')
plt.savefig('/home/ashley/Link_xray/result/var_real23',format='png',dpi=200)


#===========num>3:
indexg3=[]
varg3=[]
varg3_real=[]
#index=name3.index(nameg3[828])
#print(index)
for i in range(0,len(nameg3),1):
    fluxtmp=[]
    index_g3=name3.index(nameg3[i])
    #index_tmp=[n for n,v in enumerate(name3) if v==namel3[i]]
    #index_l3=np.min(index_tmp)
    #print('index_l3:',index_l3)
    indexg3.append(index_g3)
    #print('i',i)
    a_num=numg3[i]
    #print(index_g3)
    #print(a_num)
    for j in range(index_g3,a_num+index_g3,1):
        
        if flux3[j]==0.0:
            flux_tmp=flux3_hilim[j]
        else:
            flux_tmp=flux3[j]
        fluxtmp.append(flux_tmp)
    #print('fluxtmp:',fluxtmp)
    flux_tmpmin=np.min(fluxtmp)
    flux_tmpmax=np.max(fluxtmp)
    #print
    vary=(flux_tmpmax)/(flux_tmpmin)
    vary_real=np.var(fluxtmp)/(np.mean(flux_tmp)**2)
    #print('vary',vary)
    varg3.append(vary)
    varg3_real.append(vary_real)
#print('len_varg3',len(varg3))
#print('var:',var)
#index1=[i for i,v in enumerate(varg3) if v=='nan']
#print('index_nan:',index1) 
ymax=np.max(varg3_real)*1.1
Dgg3=[i*3E8*308568E16/67.73999999999998 for i in redg3]
Dg3=[i*3.24077929e-25 for i in Dgg3]
plt.figure()
#plt.ylim(0,ymax)
glow=plt.scatter(Dg3,varg3_real,s=4,facecolors='none',edgecolor='m',label='acis_num>3')
plt.legend(loc=0)          
plt.xscale('log')
plt.xlabel('Distance'+'  '+'Mpc')
plt.ylabel('var_real')

plt.title('var3_real')
plt.savefig('/home/ashley/Link_xray/result/var_real',format='png',dpi=200)

file3.close()

 
"""





#variance.py



"""
mean_sigma7=np.sum(sigma7)/len(sigma7)
for n in range(0,len(sigma7),1):
    tmp_err7+=(sigma7[n]-mean_sigma7)/len(sigma7)*(len(sigma7)-1)
sqr_err7=[np.sqrt(n) for n in tmp_err7]
mean_ml7=np.mean(meanl7)

mean_sigma8=np.sum(sigma8)/len(sigma8)
for n in range(0,len(sigma8),1):
    tmp_err8+=(sigma8[n]-mean_sigma8)/len(sigma8)*(len(sigma8)-1)
sqr_err8=[np.sqrt(n) for n in tmp_err8]
mean_ml8=np.mean(meanl8)

mean_sigma9=np.sum(sigma9)/len(sigma9)
for n in range(0,len(sigma9),1):
    tmp_err9+=(sigma9[n]-mean_sigma9)/len(sigma9)*(len(sigma9)-1)
sqr_err9=[np.sqrt(n) for n in tmp_err9]
mean_ml9=np.mean(meanl9)

mean_sigma10=np.sum(sigma10)/len(sigma10)
for n in range(0,len(sigma10),1):
    tmp_err10+=(sigma10[n]-mean_sigma10)/len(sigma10)*(len(sigma10)-1)
sqr_err10=[np.sqrt(n) for n in tmp_err10]
mean_ml10=np.mean(meanl10)

mean_sigma11=np.sum(sigma11)/len(sigma11)
for n in range(0,len(sigma11),1):
    tmp_err11+=(sigma11[n]-mean_sigma11)/len(sigma11)*(len(sigma11)-1)
sqr_err11=[np.sqrt(n) for n in tmp_err11]
mean_ml11=np.mean(meanl11)
"""


#collect


"""
collect info 
output=sys.stdout
f=open('/home/ashley/Link_xray/result/info_variance_acisnum>2','a+')
sys.stdout=f
"""