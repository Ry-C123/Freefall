import matplotlib.pyplot as plt
from CONFIG import *



for i in range(OUTPUT_int, n_steps, OUTPUT_int*5):
    x = []
    y = []
    for lin in open(runname+'_'+str(i)+'.simo').readlines():
        if '#' not in lin:
            tmp = lin.split(',')
            x.append(float(tmp[1]))
            y.append(float(tmp[2]))

    plt.xlim(-1*AU, 2*AU)
    plt.ylim(-1*AU, 1*AU)
    rot_s = dt*i*omega*57.2958 + 270
    plt.plot(0, 0, marker=(3, 0, rot_s), markersize=8, linestyle='None', c='k')
    TIME_YEARS = round((i*dt*3.17098e-8)/0.0836, 4)
    plt.text(2.0*AU/7, 0.9*AU/7,str(TIME_YEARS)+' orbits')
    plt.scatter(x,y, s=4, color='r')
    plt.savefig(runname+'_'+str(i)+'_plot.png')
    plt.clf()
