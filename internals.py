import numpy
import random
import time
import sys
from matplotlib import pyplot as plt
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
rcParams['legend.loc'] = 'best'
rcParams['legend.fontsize'] = 8.5
rcParams['legend.isaxes'] = False
rcParams['figure.facecolor'] = 'white'
rcParams['figure.edgecolor'] = 'white'

NUM_ITERATIONS = 1000000
ANGLE_RESOLUTION = 7.5

def val(x):
    return round(x/ANGLE_RESOLUTION)*ANGLE_RESOLUTION

def analyse():
    """Perform some more analysis on saved output files after the simulation is done"""
    alice = numpy.load('Alice.npy') # angle, outcome 
    bob = numpy.load('Bob.npy')  # angle, outcome 
    alice_hv = numpy.load('Alice-hv.npy') # hidden variables
    bob_hv = numpy.load('Bob-hv.npy') # hidden variables

    ab = alice[:,0] - bob[:,0]        
    abdeg = numpy.round(numpy.degrees(ab)/ANGLE_RESOLUTION)*ANGLE_RESOLUTION
    adeg  = numpy.round(numpy.degrees(alice[:,0])/ANGLE_RESOLUTION)*ANGLE_RESOLUTION
    bdeg  = numpy.round(numpy.degrees(bob[:,0])/ANGLE_RESOLUTION)*ANGLE_RESOLUTION
    
    sl_pp = ((alice[:,-1] == +1.0) & (bob[:,-1] == +1.0))
    sl_mm = ((alice[:,-1] == -1.0) & (bob[:,-1] == -1.0))
    sl_pm = ((alice[:,-1] == +1.0) & (bob[:,-1] == -1.0))
    sl_mp = ((alice[:,-1] == -1.0) & (bob[:,-1] == +1.0))
    sl_ap = (alice[:,-1] == +1.0)
    sl_bp = (bob[:,-1] == +1.0)
    sl_am = (alice[:,-1] == -1.0)
    sl_bm = (bob[:,-1] == -1.0)
    sl_sd = (numpy.abs(alice[:,-1])+numpy.abs(bob[:,-1]) > 0.0) # single detections
    sl_dd = (alice[:,-1]*bob[:,-1] != 0.0)                      # double detections (coincidences)
    
    pp = abdeg[sl_pp]
    mm = abdeg[sl_mm]
    pm = abdeg[sl_pm]
    mp = abdeg[sl_mp]
    bp = abdeg[sl_bp]
    ap = abdeg[sl_ap]
    am = abdeg[sl_am]
    bm = abdeg[sl_bm]
    
    x = numpy.arange(0.0, 360.0, ANGLE_RESOLUTION)
    ypp = numpy.zeros_like(x)
    ymm = numpy.zeros_like(x)
    ypm = numpy.zeros_like(x)
    ymp = numpy.zeros_like(x)
    yap = numpy.zeros_like(x)
    ybp = numpy.zeros_like(x)
    Eab = numpy.zeros_like(x)
    coinc = numpy.zeros_like(x)
    geff = numpy.zeros_like(x)
    
    for i, a in enumerate(x):
        lpp = (pp==a).sum() # ++
        lmm = (mm==a).sum() # --
        lpm = (pm==a).sum() # -+
        lmp = (mp==a).sum() # +-
        lap = (ap==a).sum() # A+
        lbp = (bp==a).sum() # B+
        Na = (ap==a).sum() + (am==a).sum()
        Nb = (bp==a).sum() + (bm==a).sum()
        tot = float(lpp + lmm + lpm + lmp)
        sel = (abdeg == a)
        coinc[i] = float((sel & sl_dd).sum())/(sel & sl_sd).sum()
        geff[i] = float((sel & sl_sd).sum())/(sel).sum()
        if tot != 0.0:
            ypp[i] = lpp/tot
            ymm[i] = lmm/tot
            ypm[i] = lpm/tot
            ymp[i] = lmp/tot
        
        Eab[i] = ypp[i] + ymm[i] - ypm[i] - ymp[i]
        # single sided +/- outcome probabilities
        if Na != 0.0:
            yap[i] = float(lap)/Na
        if Nb != 0.0:
            ybp[i] = float(lbp)/Nb
    
    # Display results   
    a = val(0.0)
    ap = val(45.0)
    b = val(22.5)
    bp = val(67.5)
    
    sl_same = (adeg == bdeg) & sl_dd
    sl_opp = (numpy.abs(adeg - bdeg) == val(180.0)) & sl_dd
    ysame = (alice[sl_same, 1] * bob[sl_same,1]).mean()
    yopp = (alice[sl_opp, 1] *  bob[sl_opp, 1]).mean()
    
    DESIG = {0 : '<a1b1>', 1: '<a2d2>', 2: '<c3b3>', 3: '<c4d4>'}    
    CHSH = []
    QM = []
    for k,(i,j) in enumerate([(a,b),(a,bp), (ap,b), (ap, bp)]):
        sel = (adeg==i) & (bdeg==j) & (alice[:,1] != 0.0) & (bob[:,1] != 0.0)
        Ai = alice[sel, 1] 
        Bj = bob[sel, 1]
        print "%s: E(%5.1f,%5.1f), AB=%+0.2f, QM=%+0.2f, N=%d" % (DESIG[k],i, j, (Ai*Bj).mean(),
             -numpy.cos(numpy.radians(j-i)), len(Ai))
        CHSH.append( (Ai*Bj).mean())
        QM.append( -numpy.cos(numpy.radians(j-i)) )
   
    print
    print "Same Angle <AB> = %+0.2f, QM = -1.00" % (ysame)
    print "Oppo Angle <AB> = %+0.2f, QM = +1.00" % (yopp)
    print "CHSH: < 2.0, Sim: %0.3f, QM: %0.3f" % (abs(CHSH[0]-CHSH[1]+CHSH[2]+CHSH[3]), abs(QM[0]-QM[1]+QM[2]+QM[3]))
    print "Detector Efficiency:  %0.1f %%" % (100.0*sl_sd.mean())  
    print "%% Coincidences:  %0.1f %%" % (100.0*sl_dd.sum()/sl_sd.sum())  
            
    gs = gridspec.GridSpec(3,2)
    ax1 = plt.subplot(gs[0,:])    
    ax1.plot(x, Eab, 'm-', label='Model: E(a,b)')
    ax1.plot(x, -numpy.cos(numpy.radians(x)), 'b.:', label='QM: -cos(ab)')
    ax1.plot([0.0, 180.0, 360.0], [-1.0, 1.0, -1.0], 'r--')
    ax1.legend()
    
    ax2 = plt.subplot(gs[1,0])       
    ax2.plot(x, ypp, label='P++')
    ax2.plot(x, ymm, label='P--')
    ax2.plot(x, ypm, label='P+-')
    ax2.plot(x, ymp, label='P-+')
    ax2.plot(x, yap, label='P+(A)')
    ax2.plot(x, ybp, label='P+(B)')
    ax2.legend()
    
    ax3 =  plt.subplot(gs[1,1])       
    ax3.plot(x, 100*coinc, 'b--', label='% Coincidence')
    ax3.plot(x, 100*geff, 'r--', label='% Detector Efficiency')
    ax3.set_ylim(0,105)
    ax3.legend()
    
    ax4 = plt.subplot(gs[2,:])
    ax4.hist(alice_hv[:,1], 50, normed=1, histtype='step', label='emitted', color='black')    
    ax4.hist(alice_hv[sl_sd,1], 50, normed=1, histtype='step', label='detected one side', color='green')    
    ax4.hist(alice_hv[sl_dd,1], 50, normed=1, histtype='step', label='detected both sides', color='red')
    ax4.set_xlabel('histogram of p values')
    ax4.legend()

    for ax in [ax1, ax2, ax3]:  ax.set_xlim(0, 360)

    plt.savefig('internals.png', dpi=72)
    plt.show()
           
if __name__ == '__main__':
    analyse()        
