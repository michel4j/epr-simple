import numpy
import random
import time
import sys

NUM_ITERATIONS = 10000000

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

class Source(object):
    """Generate and emit two particles with hidden variables"""
    def __init__(self, spin=1, phase=numpy.pi):
        self.spin = spin
        self.phase = phase
        
    def emit(self):
        e = numpy.random.uniform(0.0,2*numpy.pi)
        #p1 = numpy.random.normal(loc=0.5, scale=0.25)
        p1 = numpy.random.uniform(0.0, 1.0) 
        p2 = 1 - p1                 # complementary entangled particles
        lp = (e, p1, self.spin)  # Left particle        
        rp = (e+self.phase, p2, -self.spin)  # Right particle
        return lp, rp

class Station(object):
    """Detect a particle with a given/random setting"""
    def __init__(self):
        """Initialize with fixed_angle, otherwise random fast switching will be used"""
        self.results = numpy.empty((NUM_ITERATIONS, 2))  # columns: angle, +1/-1/0 outcome
        self.results.fill(numpy.nan)
        self.count = 0
        
    def get_setting(self):
        """Return the current detector setting"""
        return numpy.random.uniform(0.0,numpy.pi*2)
        
    def detect(self, particle):
        """Calculate the station outcome for the given `particle`"""
        a = self.get_setting()
        e, p, s = particle
        g = numpy.sign(s)
        n = 2*abs(s)
        C = ((-1)**n)*numpy.cos(n*(a-e))/2.0
        Cd = (abs(C) + 0.5)**2
        if g < 0:  Cd = 1 - Cd # Flip for the other side
         
        if (p > Cd and g > 0) or (p < Cd and g < 0):
            out = 0 # not detected
        else:
            out = numpy.sign(C)
                        
        self.results[self.count] = numpy.array([a, out]) # save angle, outcome pair
        self.count += 1


class Simulation(object):
    def __init__(self):
        self.source = Source(spin=0.5)
        self.alice = Station()
        self.bob = Station()
        
    def run(self):
        """generate and detect particles"""
        progress = ProgressMeter(total=NUM_ITERATIONS)
        for i in range(NUM_ITERATIONS):
            left_particle, right_particle = self.source.emit()
            self.alice.detect(left_particle)
            self.bob.detect(right_particle)
            progress.update(1)
            
        # save results in separate files each is a list of angles and +/- outcomes
        numpy.save("Alice", self.alice.results)
        numpy.save("Bob", self.bob.results)

    def analyse(self, load=False):
        """select coincidences and calculate probabilities for each angle diff"""
        if not load:
            alice = self.alice.results
            bob = self.bob.results
        else:
            alice = numpy.load('Alice.npy') # angle, outcome 
            bob = numpy.load('Bob.npy')  # angle, outcome 
        ab = alice[:,0] - bob[:,0]
        abdeg = numpy.round(numpy.degrees(ab))
        sl_pp = ((alice[:,-1] == +1.0) & (bob[:,-1] == +1.0))
        sl_mm = ((alice[:,-1] == -1.0) & (bob[:,-1] == -1.0))
        sl_pm = ((alice[:,-1] == +1.0) & (bob[:,-1] == -1.0))
        sl_mp = ((alice[:,-1] == -1.0) & (bob[:,-1] == +1.0))
        sl_ap = (alice[:,-1] == +1.0)
        sl_bp = (bob[:,-1] == +1.0)
        sl_am = (alice[:,-1] == -1.0)
        sl_bm = (bob[:,-1] == -1.0)
        sl_nd = ((alice[:,-1] != 0.0) | (bob[:,-1] != 0.0))
        sl_nc = ((alice[:,-1] != 0.0) & (bob[:,-1] != 0.0))
        
        pp = abdeg[sl_pp]
        mm = abdeg[sl_mm]
        pm = abdeg[sl_pm]
        mp = abdeg[sl_mp]
        bp = abdeg[sl_bp]
        ap = abdeg[sl_ap]
        am = abdeg[sl_am]
        bm = abdeg[sl_bm]
        
        x = numpy.arange(360.)
        ypp = numpy.zeros_like(x)
        ymm = numpy.zeros_like(x)
        ypm = numpy.zeros_like(x)
        ymp = numpy.zeros_like(x)
        yap = numpy.zeros_like(x)
        ybp = numpy.zeros_like(x)
        Eab = numpy.zeros_like(x)

        sl_same = (abdeg == 0.0) & ((alice[:,-1] != 0.0) & (bob[:,-1] != 0.0))
        sl_opp = (abdeg == 180.0) & ((alice[:,-1] != 0.0) & (bob[:,-1] != 0.0))
        ysame = (alice[sl_same,1] *  bob[sl_same,1]).mean()
        yopp = (alice[sl_opp,1] *  bob[sl_opp,1]).mean()
        ceff = 100.0*(1 - len(abdeg[sl_nc])/float(len(abdeg[sl_nd])))
        aeff = 100.0*(1 - len(abdeg[sl_nd])/float(len(abdeg)))
        print "Absolute Efficiency %0.4f %%" % (aeff)
        print "Conditional efficiency %0.4f %%" % (ceff)
        print "Same Angle <AB> = %0.4f, QM = -1.0" % (ysame)
        print "Opposite Angle <AB> = %0.4f, QM = 1.0" % (yopp)
        
        for a in x:
            i = int(a)
            lpp = len(pp[(pp==a)]) # ++
            lmm = len(mm[(mm==a)]) # --
            lpm = len(pm[(pm==a)]) # -+
            lmp = len(mp[(mp==a)]) # +-
            lap = len(ap[ap==a])
            lbp = len(bp[bp==a])
            Na = len(ap[ap==a]) +len(am[am==a])
            Nb = len(bp[bp==a]) + len(bm[bm==a])
            tot = lpp + lmm + lpm + lmp
            
            if tot != 0.0:
                ypp[i] = float(lpp)/tot
                ymm[i] = float(lmm)/tot
                ypm[i] = float(lpm)/tot
                ymp[i] = float(lmp)/tot
            
            Eab[i] = ypp[i] + ymm[i] - ypm[i] - ymp[i]
            # single sided +/- outcome probabilities
            if Na != 0.0:
                yap[i] = float(lap)/Na
            if Nb != 0.0:
                ybp[i] = float(lbp)/Nb
        
        # Plot results   
        from matplotlib import pyplot as plt
        from matplotlib import rcParams
        
        CHSH = Eab[23] - Eab[68] + Eab[360+23-45] + Eab[68-45]
        QM = (-numpy.cos(numpy.radians(23)) +
              numpy.cos(numpy.radians(68)) -
              numpy.cos(numpy.radians(23-45)) -
              numpy.cos(numpy.radians(68-45))
             )
        rcParams['legend.loc'] = 'best'
        rcParams['legend.fontsize'] = 8.5
        rcParams['legend.isaxes'] = False
        rcParams['figure.facecolor'] = 'white'
        rcParams['figure.edgecolor'] = 'white'

        plt.plot(x, ypp, label='++')
        plt.plot(x, ymm, label='--')
        plt.plot(x, ypm, label='+-')
        plt.plot(x, ymp, label='-+')
        plt.plot(x, yap, label='A+')
        plt.plot(x, ybp, label='B+')
        plt.plot(x, Eab, label='E(a,b)')
        plt.plot(x, -numpy.cos(numpy.radians(x)), 'r:', label='-cos(ab)')
        plt.plot([0.0, 180.0, 360.0], [-1.0, 1.0, -1.0], 'r--')
        plt.legend()
        plt.savefig('epr.png')
        print "CHSH = %0.4f, Classical <= 2, QM=%0.4f" % (abs(CHSH), abs(QM))
        plt.show()
           
if __name__ == '__main__':
    sim = Simulation()
    if len(sys.argv) == 1:
        sim.run()
        load = False
    else:
        load = True
    sim.analyse(load)        
