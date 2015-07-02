#!/usr/bin/env python
#    Simple EPR simulation
#    Copyright (C) 2015  Michel Fodje
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

from __future__ import division
import numpy
import sys
import time
import gzip
import os

class Source(object):
    """Generate and emit two particles with hidden variables"""
    def __init__(self, spin=1.0):
        self.left = []
        self.right = []
        self.n = 2*spin
        self.phase = self.n*numpy.pi
        self.angles = numpy.linspace(0, 2*numpy.pi, 33)
        self.ps = 0.5*numpy.sin(numpy.linspace(0, numpy.pi/2, 1000))**2
        
    def emit(self):
        e = numpy.random.choice(self.angles)
        p = numpy.random.choice(self.ps)
        self.left.append(numpy.array([e, p, self.n]))  
        self.right.append(numpy.array([e+self.phase, p, self.n]))  

    def save(self, fname, a):
        f = gzip.open(fname, 'wb')
        numpy.save(f, a)
        f.close()

    def run(self, duration=60.0):
        start_t = time.time()
        print "Generating particle spin-{0} particle pairs".format(self.n*0.5)
        count = 0
        while time.time() - start_t <= duration:
            self.emit()
            count += 1
            sys.stdout.write("\rETA: %4ds [%8d pairs generated]" % (duration - time.time() + start_t, count))
            sys.stdout.flush()
        self.save('SrcLeft.npy.gz', numpy.array(self.left))
        self.save('SrcRight.npy.gz', numpy.array(self.right))
        print
        print "%d particles in 'SrcLeft.npy.gz'" % (len(self.left))
        print "%d particles in 'SrcRight.npy.gz'" % (len(self.right))

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: \n\t source.py <duration in seconds> <spin>\n"
    else:
        if len(sys.argv) == 3:
            spin = float(sys.argv[2])
        else:
            spin = 1.0
        duration = float(sys.argv[1])
        source = Source(spin=spin)
        source.run(duration=duration)
