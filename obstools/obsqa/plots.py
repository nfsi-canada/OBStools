from ObsQA.imports import *

def overplot_taup(arrivals):
        ii=4
        ang = 30
        for a in arrivals:
                ii+=1
                y = plt.gca().get_ylim()
                xl = plt.gca().get_xlim()
                x = [a[1],a[1]]
                alg = ['enter','left','right']
                e = -1
                if a[0][-1].upper()=='S':
                        c = 'b'
                else:
                        c = 'r'
                if (len([m.start() for m in re.finditer('S',a[0].upper())])<=1) or (len([m.start() for m in re.finditer('I',a[0].upper())])>=1) or (len([m.start() for m in re.finditer('C',a[0].upper())])<1):
                        if x[0]<np.max(xl):
                                if len(a[0])<2:
                                        plt.plot(x,y,c,linewidth=1,linestyle='-.')
                                        if a[0][-1].upper()=='S':
                                                plt.text(x[0],(y[1]*(0.1)),a[0],rotation=0,color='b',fontsize=20,horizontalalignment='right',fontweight='bold')
                                        else:
                                                plt.text(x[0],(y[1]*(0.1)),a[0],rotation=0,color='r',fontsize=20,horizontalalignment='right',fontweight='bold')
                                elif len(a[0])<3:
                                        plt.plot(x,y,c,linewidth=0.3,alpha=0.2,linestyle='-.')
                                        if a[0][-1].upper()=='S':
                                                plt.text(x[0],(y[1]*(0.4)) * (-1),a[0],rotation=-0,color='b',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',alpha=0.5)
                                        else:
                                                plt.text(x[0],(y[1]*(0.4)) * (1),a[0],rotation=-0,color='r',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',alpha=0.5)
                                elif len(a[0])==3:
                                        plt.plot(x,y,c,linewidth=0.3,alpha=0.2,linestyle='-.')
                                        if a[0][-1].upper()=='S':
                                                plt.text(x[0],(y[1]*(0.5)) * (-1),a[0],rotation=-0,color='b',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',alpha=0.5,verticalalignment='bottom')
                                        else:
                                                plt.text(x[0],(y[1]*(0.5)) * (1),a[0],rotation=-0,color='r',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',alpha=0.5,verticalalignment='top')
                                elif (len([m.start() for m in re.finditer('I',a[0].upper())])>=1):
                                        plt.plot(x,y,'gray',linewidth=0.3,alpha=0.9,linestyle='-.')
                                        if a[0][-1].upper()=='S':
                                                plt.text(x[0],(y[1]*(1)) * (-1),a[0],rotation=90,color='gray',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',verticalalignment='bottom',alpha=0.5)
                                        else:
                                                plt.text(x[0],(y[1]*(1)) * (1),a[0],rotation=90,color='gray',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',verticalalignment='top',alpha=0.5)
                                elif (len([m.start() for m in re.finditer('K',a[0].upper())])>=1):
                                        plt.plot(x,y,'gray',linewidth=0.3,alpha=0.9,linestyle='-.')
                                        if a[0][-1].upper()=='S':
                                                plt.text(x[0],(y[1]*(0.7)) * (-1),a[0],rotation=90,color='gray',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',verticalalignment='bottom')
                                        else:
                                                plt.text(x[0],(y[1]*(0.6)) * (1),a[0],rotation=90,color='gray',horizontalalignment=alg[e**np.abs(ii)],fontweight='bold',verticalalignment='top')