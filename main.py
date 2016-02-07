import matplotlib.pyplot as pl
import matplotlib.animation as anim
import numpy as np

import objects as obj

pl.close('all')

#Traçé :

fig = pl.figure()

ax = fig.add_subplot(111)
#energy = fig.add_subplot(212, xlim=(-1,15), ylim=(-1,5))

parts = [
    obj.Part(function=lambda t:10, param='POLAR', xlim=(-np.pi/2,np.pi/2), origin=(10,-20.1)),
    obj.Part(function=lambda x:0, xlim=(-30,30)),
    obj.Part(function=lambda x:-30, xlim=(-30,10)),
    obj.Part(function=lambda x:-x-40, xlim=(-60,-20)),
    obj.Part(function=lambda x:1/2*x-50, xlim=(-60,50)),
    obj.Part(function=lambda x:-1/2*x-50, xlim=(-60,50)),

]
"""
parts = [
    obj.Part(function=lambda x:1/100*x**2, xlim=(-60,50)),
    obj.Part(function=lambda x:-30, xlim=(-60,50)),
]
"""
M = None#obj.Point((2,0.15),10, (1,0))

#Ep, = energy.plot([],[],label='$E_p$')
#Ec, = pl.plot([],[],label='$E_c$')

Ep_val = []

def init():
    global M
    M = obj.Point((2,8),10, (-20,1))
    for o in parts:
        X,Y = o.get_points(600)
        ax.plot(X,Y, c='blue')

def update(i,clear=True):
    #print('\n')
    if clear or i==0:
        ax.clear()
        for o in parts:
            X,Y = o.get_points(600)
            ax.plot(X,Y, c='blue')
    M.draw(ax)
    M.move(parts, i/500,draw_debug=ax)
    #Ep_val.append(M.mass*9.81*M.pos[1])
    #   Ep.set_data(np.linspace(0,i/30, len(Ep_val)), Ep_val)
    #return Ep

pl.axis('equal')

init()
from time import time
t0 = time()
#update(1)
t1 = time()
interval = 1000 * 1/30 - (t1 - t0)
ani = anim.FuncAnimation(fig, update, 190,interval=5)

#ani.save('test.mp4')
#for i in range(1,40):
#    update(i,False)
#init()
pl.show()