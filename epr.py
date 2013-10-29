import numpy
import random
import time
import sys

NUM_ITERATIONS = 10000000

class ProgressMeter(object):
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
    def __init__(self, variant=1.0):
        self.variant= variant
        
    def emit(self):
        e = numpy.random.uniform(0.0,numpy.pi*2)
        p1 = numpy.random.uniform(0.0, 1.0)
        p2 = numpy.random.uniform(0.0, 1.0)
        lp = (e, p1, self.variant)  # Left particle        
        rp = (e, p2, self.variant)  # Right particle
        return lp, rp

class Station(object):
    def __init__(self):
        self.results = numpy.empty((NUM_ITERATIONS, 2))  # angle, channel
        self.results.fill(numpy.nan)
        self.count = 0
        
    def get_setting(self):
        return numpy.random.uniform(0.0,numpy.pi*2)
        
    def detect(self, particle):
        a = self.get_setting()
        e, p, s = particle
        threshold = numpy.random.uniform(0, 1)
        C = numpy.cos(s*(e - a))**2 - p
        if abs(C) > threshold:
            out = numpy.sign(C)
        else:
            out = 0.0             
        self.results[self.count] = numpy.array([a, out])
        self.count += 1


class Simulation(object):
    def __init__(self):
        self.source = Source(0.5)
        self.alice = Station()
        self.bob = Station()
        
    def run(self):
        # generate and detect particles
        progress = ProgressMeter(total=NUM_ITERATIONS)
        for i in range(NUM_ITERATIONS):
            left_particle, right_particle = self.source.emit()
            self.alice.detect(left_particle)
            self.bob.detect(right_particle)
            progress.update(1)
            
        # save results in separate files each is a list of angles and +/- outcomes
        numpy.save("Alice", self.alice.results)
        numpy.save("Bob", self.bob.results)

    def analyse(self):
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

        pp = abdeg[sl_pp]
        mm = abdeg[sl_mm]
        pm = abdeg[sl_pm]
        mp = abdeg[sl_mp]
        bp = abdeg[sl_bp]
        ap = abdeg[sl_ap]
        am = abdeg[sl_am]
        bm = abdeg[sl_bm]
        
        x = numpy.arange(180.)
        ypp = numpy.zeros_like(x)
        ymm = numpy.zeros_like(x)
        ypm = numpy.zeros_like(x)
        ymp = numpy.zeros_like(x)
        yap = numpy.zeros_like(x)
        ybp = numpy.zeros_like(x)
        Eab = numpy.zeros_like(x)
        
        # select coincidences and calculate probabilities for each angle diff
        for a in x:
            i = int(a)
            lpp = len(pp[(pp==a)]) # ++
            lmm = len(mm[(mm==a)]) # --
            lpm = len(pm[(pm==a)]) # -+
            lmp = len(mp[(mp==a)]) # +-
            tot = lpp + lmm + lpm + lmp
            
            ypp[i] = float(lpp)/tot
            ymm[i] = float(lmm)/tot
            ypm[i] = float(lpm)/tot
            ymp[i] = float(lmp)/tot
            
            Eab[i] = ypp[i] + ymm[i] - ypm[i] - ymp[i]
            
            # single sided +/- outcome probabilities
            yap[i] = float(len(ap[ap==a]))/(len(ap[ap==a]) + len(am[am==a]))
            ybp[i] = float(len(bp[bp==a]))/(len(bp[bp==a]) + len(bm[bm==a]))
        
        # Plot results   
        from matplotlib import pyplot as plt
        from matplotlib import rcParams
        
        rcParams['legend.loc'] = 'best'
        rcParams['legend.fontsize'] = 8.5
        rcParams['legend.isaxes'] = False
        rcParams['figure.facecolor'] = 'white'
        rcParams['figure.edgecolor'] = 'white'
        rcParams['font.sans-serif'] = 'Cantarell'

        plt.plot(x, ypp, label='++')
        plt.plot(x, ymm, label='--')
        plt.plot(x, ypm, label='+-')
        plt.plot(x, ymp, label='-+')
        plt.plot(x, yap, label='A+')
        plt.plot(x, ybp, label='B+')
        plt.plot(x, Eab, label='E(a,b)')
        plt.legend()
        plt.savefig('epr.png')
        plt.show()
           
if __name__ == '__main__':
    sim = Simulation()
    sim.run()
    sim.analyse()        
