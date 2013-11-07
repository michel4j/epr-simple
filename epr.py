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

NUM_ITERATIONS = 50000000
MAX_ANGLE = 360.0
ANGLE_RESOLUTION = 7.5

class Source(object):
    """Generate and emit two particles with hidden variables"""
    def __init__(self, spin=1):
        self.spin = spin
        self.phase = self.spin*2*numpy.pi
        self.n = 2*self.spin
        self.angles = numpy.radians(numpy.arange(0, 360.0, ANGLE_RESOLUTION))  # pre-calculate angles to choose from 
        #self.ps = numpy.linspace(0, numpy.pi/4, 1000)**(numpy.pi/2) # pre-calculate p-values to choose from
        self.ps = 0.5*numpy.sin(numpy.linspace(0, numpy.pi/2, 1000))**2 # alternate p-distribution
        
    def emit(self):
        e = numpy.random.choice(self.angles) # pick one of the angles randomly 
        p = numpy.random.choice(self.ps)  # pick one of the p values randomly
        lp = (e, p, self.n)  # Left particle        
        rp = (e+self.phase, p, self.n)  # Right particle
        return lp, rp

class Station(object):
    """Detect a particle with a given/random setting"""
    def __init__(self, name='Alice'):
        self.name = name
        # create arrays to store results for this station
        self.results = numpy.empty((NUM_ITERATIONS, 2))  # columns: angle, +1/-1/0 outcome        
        self.properties = numpy.empty((NUM_ITERATIONS, 2)) # columns: e, p
        self.count = 0
        self.angles = numpy.radians(numpy.arange(0, MAX_ANGLE, ANGLE_RESOLUTION)) # pre-calculate p-values to choose from
        
    def get_setting(self):
        """Return the current detector setting"""
        return numpy.random.choice(self.angles) # pick one angle randomly from pre-calculated set
        
    def detect(self, particle):
        """Calculate the station outcome for the given `particle`"""
        a = self.get_setting()
        e, p, n = particle        
        C = ((-1)**n)*numpy.cos(n*(a-e))
        out =  p < abs(C) and numpy.sign(C) or 0.0  # if |C| > p then particle is detected by sign(C) channel
        self.results[self.count] = numpy.array([a, out]) # save angle, outcome
        self.properties[self.count] = numpy.array([e, p]) # save e, p
        self.count += 1 

    def save(self):
        """Save the results"""
        numpy.save("%s" % self.name, self.results)
        numpy.save("%s-hv" % self.name, self.properties)
        

class Simulation(object):
    def __init__(self):
        self.source = Source(spin=0.5)
        self.alice = Station('Alice')
        self.bob = Station('Bob')
        
    def run(self):
        """generate and detect particles"""
        progress = ProgressMeter(total=NUM_ITERATIONS)
        for i in range(NUM_ITERATIONS):
            left_particle, right_particle = self.source.emit() # generate two particles from source
            self.alice.detect(left_particle)  # send left particle to Alice
            self.bob.detect(right_particle)  # send right particle to Bob
            progress.update(1)
            
        # save results in separate files each is a list of angles and +/- outcomes
        self.alice.save()
        self.bob.save()

    def analyse(self, load=False):
        """select coincidences and calculate probabilities for each angle diff"""
        if not load:
            alice = self.alice.results
            bob = self.bob.results
            alice_hv = self.alice.properties
            bob_hv = self.bob.properties
        else:
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
        
        x = numpy.arange(0.0, MAX_ANGLE, ANGLE_RESOLUTION)
        ypp = numpy.zeros_like(x)
        ymm = numpy.zeros_like(x)
        ypm = numpy.zeros_like(x)
        ymp = numpy.zeros_like(x)
        yap = numpy.zeros_like(x)
        ybp = numpy.zeros_like(x)
        Eab = numpy.zeros_like(x)
        coinc = numpy.zeros_like(x)
        
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
            print "%s: E(%5.1f,%5.1f), AB=%+0.2f, QM=%+0.2f" % (DESIG[k],i, j, (Ai*Bj).mean(),
                 -numpy.cos(numpy.radians(j-i)))
            CHSH.append( (Ai*Bj).mean())
            QM.append( -numpy.cos(numpy.radians(j-i)) )
       
        print
        print "Same Angle <AB> = %+0.2f, QM = -1.00" % (ysame)
        print "Oppo Angle <AB> = %+0.2f, QM = +1.00" % (yopp)
        print "CHSH: < 2.0, Sim: %0.3f, QM: %0.3f" % (abs(CHSH[0]-CHSH[1]+CHSH[2]+CHSH[3]), abs(QM[0]-QM[1]+QM[2]+QM[3]))
        print "Detector Efficiency:  %0.1f %%" % (100.0*sl_sd.mean())  
        print "%% Coincidences:  %0.1f %%" % (100.0*sl_dd.sum()/sl_sd.sum())  
          
        gs = gridspec.GridSpec(2,1)
        ax1 = plt.subplot(gs[0])    
        ax1.plot(x, Eab, 'm-', label='Sim: E(a,b)')
        ax1.plot(x, -numpy.cos(numpy.radians(x)), 'b.', label='QM: -cos(ab)')
        ax1.plot([0.0, 180.0, 360.0], [-1.0, 1.0, -1.0], 'r--')
        ax1.legend()
        
        ax2 = plt.subplot(gs[1])       
        ax2.plot(x, ypp, label='P++')
        ax2.plot(x, ymm, label='P--')
        ax2.plot(x, ypm, label='P+-')
        ax2.plot(x, ymp, label='P-+')
        ax2.plot(x, yap, label='P+(A)')
        ax2.plot(x, ybp, label='P+(B)')
        ax2.legend()
        
        for ax in [ax1, ax2]:  ax.set_xlim(0, MAX_ANGLE)
        ax1.set_title('%d emitted particle pairs' % NUM_ITERATIONS)
        plt.savefig('epr.png', dpi=100)
        plt.show()

def val(x):
    """Round value to the nearest ANGLE_RESOLUTION"""
    return round(x/ANGLE_RESOLUTION)*ANGLE_RESOLUTION

class ProgressMeter(object):
    """Displays a progress bar so we know how long it will take"""
    def __init__(self, **kw):
        self.total = int(kw.get('total', 100)) # Number of units to process
        self.count = int(kw.get('count', 0))   # Number of units already processed
        self.refresh_rate = float(kw.get('rate_refresh', .5)) # Refresh rate in seconds
        self.meter_ticks = int(kw.get('ticks', 50)) # Number of ticks in meter
        self.meter_division = float(self.total) / self.meter_ticks
        self.meter_value = int(self.count / self.meter_division)
        self.current_rate = 0.0
        self.last_refresh = 0
        self.start_time = time.time()
        sys.stdout.write(chr(27) + '[s')

    def update(self, count, **kw):
        """Increment progres count by `count` amount"""
        now = time.time()
        self.count += count 
        self.current_rate = self.count/(now - self.start_time)
        self.meter_value = max(int(self.count / self.meter_division), self.meter_value)
        if (now - self.last_refresh) > self.refresh_rate or (self.count >= self.total):
            sys.stdout.write(chr(27) + '[2K') # clear the line and reset cursor
            sys.stdout.write(chr(27) + '[u')
            sys.stdout.write(chr(27) + '[s')
            meter = '[%s>%s] %d%%  %10.2e/sec' % (
                    '-' * self.meter_value, ' ' * (self.meter_ticks - self.meter_value), 
                    (float(self.count) / self.total) * 100, self.current_rate)
            sys.stdout.write(meter)
            if self.count >= self.total:
                sys.stdout.write('\n')
            sys.stdout.flush()
            self.last_refresh = time.time()
           
if __name__ == '__main__':
    sim = Simulation()
    if len(sys.argv) == 1:
        sim.run()
        load = False
    else:
        load = True
    sim.analyse(load)        
