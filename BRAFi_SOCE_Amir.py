


import pylab
import os
import boolean2
from boolean2 import Model
import numpy as np
from boolean2 import util


#main_folder= "C:\users\amir\desktop\ Bi_weekley_melanoma_meeting\boolean_modeling"

n_delays=50

Initial_conditions='''
BRAFi=False
BRAF= True
MEK= False
ERK=False
Gene_exp=False
Ca_ext=True
Ca_channel= True
Ca_cyt_1=False
pumpi=False
Ca_pump_ER=True
Ca_ER=False
Ca_cyt_2=False
Ca_cyt_3=False
'''
for i in range(n_delays):
    Initial_conditions+='Delay%d=False\n'%i

rules='''
BRAF*= not BRAFi
MEK*= BRAF or(Ca_cyt_1 and Ca_cyt_2 and Ca_cyt_3)
ERK*=MEK
Ca_channel*=ERK or Gene_exp
Ca_cyt_1*=Ca_cyt_1 or (Ca_ext and Ca_channel)or Ca_ER
Ca_pump_ER*= not pumpi
Ca_ER*=Ca_cyt_1 and Ca_pump_ER
Ca_cyt_2*= Ca_cyt_1 and not Ca_pump_ER
Ca_cyt_3*= Ca_cyt_2 and Ca_ext and Ca_channel
'''
#Delay1*= not ERK
#Delay2*= Delay1 and not ERK
#Delay3*= Delay2 and not ERK
#Gene_exp*= Delay3 and not ERK

for i in range(n_delays):
    if i == 0:
        rules+='Delay0*= not ERK\n'
    else:
        rules+= 'Delay%d*= Delay%d and not ERK\n'%(i,i-1)
    if i== n_delays-1:
        rules+= 'Gene_exp*= Delay%d and not ERK\n'%(n_delays-1)



text=Initial_conditions+rules



species_to_plot= ["BRAF","MEK","ERK","Gene_exp","Ca_channel","Ca_pump_ER","Ca_ER"]
marker=["o", "^","*","|","+","d","H"]
colors=["green","grey","black","yellow","red","purple","brown"]





"""
model = Model(text=text, mode='sync')
model.initialize()
model.iterates(steps=10)

for node in model.data:
    print (node,model.data[node])

p1= pylab.plot(model.data["Ca_cyt_1"],'sb-',label= 'Ca_cyt_1')
p2= pylab.plot(model.data["ERK"], 'go--', label= 'ERK')
p3= pylab.plot(model.data["BRAFi"], 'og-', label='BRAFi')

pylab.legend('best')
pylab.ylim((-0.1,1.1))
pylab.show()
"""


coll= util.Collector()

for i in range(1000):
    print(i)

    # step 1 , equilibration
    model = boolean2.Model(text, mode='async')
    model.initialize()
    model.iterate(steps=10)
    coll.collect(states=model.states, nodes=model.nodes)
    #print (coll.store)


    # step 2 , add BRAF inhibitor
    #turnon= []
    #turnoff=[]
    Initial_conditions=""
    for species in model.states[-1].keys():
        if species == "BRAFi" or model.states[-1][species] is True:
            #turnon.append(species)
            Initial_conditions+='%s=True\n'%species
        else:
            #turnoff.append(species)
            Initial_conditions += '%s=False\n' % species
    #new_txt = boolean2.modify_states(text, turnon=turnon, turnoff=turnoff)
    new_txt= Initial_conditions+rules
    #print(new_txt)
    #quit()

    model = boolean2.Model(new_txt, mode='async')
    model.initialize()
    model.iterate(steps=10)
    coll_tmp = util.Collector()
    coll_tmp.collect(states=model.states, nodes=model.nodes)
    for species in coll.store.keys():
        coll.store[species][i] += coll_tmp.store[species][0][1:]

    #print (coll.store)
    #quit()

    #step3= remove external calcium AND add pump inhibitor

    Initial_conditions = ""
    for species in model.states[-1].keys():
        if (species != "Ca_ext" and model.states[-1][species] is True) \
                or species=="pumpi":
            Initial_conditions+= '%s=True\n' % species
        else:
            Initial_conditions += '%s=False\n' % species

    new_txt = Initial_conditions + rules
    #print(new_txt)
    #quit()

    model = boolean2.Model(new_txt, mode='async')
    model.initialize()
    model.iterate(steps=10)
    coll_tmp = util.Collector()
    coll_tmp.collect(states=model.states, nodes=model.nodes)
    for species in coll.store.keys():
        coll.store[species][i] += coll_tmp.store[species][0][1:]

    #step4  add external Calcium

    Initial_conditions = ""
    for species in model.states[-1].keys():
        if species == "Ca_ext" or model.states[-1][species] is True:
            Initial_conditions += '%s=True\n' % species
        else:
            Initial_conditions += '%s=False\n' % species

    new_txt = Initial_conditions + rules
    #print(new_txt)
    # quit()

    model = boolean2.Model(new_txt, mode='async')
    model.initialize()
    model.iterate(steps=10)
    coll_tmp = util.Collector()
    coll_tmp.collect(states=model.states, nodes=model.nodes)
    for species in coll.store.keys():
        coll.store[species][i] += coll_tmp.store[species][0][1:]

avgs = coll.get_averages(normalize=True)

Ca_cyt_1=np.array(avgs['Ca_cyt_1'])
Ca_cyt_2=np.array(avgs['Ca_cyt_2'])
Ca_cyt_3=np.array(avgs['Ca_cyt_3'])

Ca_cyt_avg= (Ca_cyt_1+Ca_cyt_2+Ca_cyt_3)/3.0


#pylab.plot(avgs['Ca_cyt_1'], 'sb-', label='Ca_cyt_1')
pylab.plot(Ca_cyt_avg, 'sb-', label='Ca_cyt')
for specie,m,c in zip(species_to_plot,marker,colors):
    pylab.plot(avgs[specie],marker= m , color=c, label=specie)

#pylab.plot(avgs['BRAF'], 'vr-', label='BRAF')
#pylab.plot(avgs['ERK'], 'og-', label='ERK')
#pylab.plot(avgs['Ca_cyt_2'], '*y-', label='Ca_cyt_2')


pylab.legend(loc="best", ncol=2)
pylab.ylim((-0.1, 1.1))

pylab.show()

