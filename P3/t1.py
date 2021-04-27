#!/usr/bin/env python
#
#  Sample demonstration of particle filters usage for localization purposes.
#  Copyright 2014 Stanislav Vechet.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#

from pylab import *
from random import *
from math import *

__samples = 30
__steps = 9
#
# believe
__sigma_f = 0.1
__sigma_d = 0.02
#
# movement
__alpha = 0
__dist = 0.5
#
# box
__x1, __y1 = 0.0, 0.1
__x2, __y2 = 1.5, -0.5
#
# functions

def motion_mdl(state):
    x,y,f = state
    f = gauss(f+__alpha,__sigma_f)
    dist = gauss(__dist, __sigma_d)
    x += dist * cos(f)
    y += dist * sin(f)
    return x, y, f

def particles(p):
    return [motion_mdl(state) for i in range(__samples) for state in p]

def weight(m,state):
    x,y,f = state
    if x < __x1 or x > __x2 or y > __y1 or y < __y2:
        return 0.0
    elif m == True:
        if abs(x-__x1) < 0.02 or \
           abs(x-__x2) < 0.02 or \
           abs(y-__y1) < 0.02 or \
           abs(y-__y2) < 0.02:
               return 1.0
        else:
            return 0.0
    return 1.0

def resampling(state_space, state_weights):
    size, count = 0, 0;
    result_space, result_weights = [],[]
    while (size < __samples):
        count = randrange(len(state_weights))
        if (random() < state_weights[count]):
            result_space.append(state_space[count])
            result_weights.append(1. / __samples)
            del state_space[count]
            del state_weights[count]
            size+=1
    return result_space, result_weights

def estimate(states):
    ex, ey, ef = 0,0,0
    for x,y,f in states:
        ex += x
        ey += y
        ef += ef
    return ex/len(states),ey/len(states),ef/len(states)
#
# initial settings
x,y,f = 0,0,0
states, allstates = [],[]
pos = [(0,0),(0.5,0),(1.0,0),(1.5,0),(1.5,-0.5),(1.0,-0.5),(0.5,-0.5),(0.0,-0.5),(0,0)]
pos = pos[0:__steps]
estim_pos = []
#
# know position
p = [(x,y,f)]
states.append((x,y,f))
allstates.append((x,y,f))
#
# run the simulation for given number of time steps
print 'Simulation started. Sometimes happends that it newer ends, \
because its simple and particular RANDOM:-) So if it doesnt finish \
until you read this, try it again ...',
for t in range(1,__steps):
    #
    # turning in specific moments
    if t == 4:
        __alpha = -3.14/2
    if t == 6:
        __alpha = 0
    if t == 8:
        __alpha = -3.14/2
    #
    # generate set of new particles
    p = particles(p)
    #
    # simulate a measurement in specific steps
    m = False
    if t == 3 and t == 4 and t == 7 and t == 8:
        m = True
    #
    # incorporate measurement into weigths
    w = [weight(m,state) for state in p]
    #
    # normalize
    w = [we/sum(w) for we in w]
    allstates += p # add before resamplig to store the source distribution
    #
    # resampling except in the first step,
    # because in the first step the number of new particles is the same as resampled
    # due to the fact that the first set is sampled from initial single position
    if t != 1:
        p, w = resampling(p,w)
    states += p
    estim_pos.append(estimate(p))
#=============
#
print 'OK, this time it worked!'
print 'Plotting results. Wait for awhile ...',
#
# show results
#
# plot box
plot([__x1,__x2],[__y1,__y1],'k-')
plot([__x2,__x2],[__y1,__y2],'k-')
plot([__x1,__x2],[__y2,__y2],'k-')
plot([__x1,__x1],[__y1,__y2],'k-')
#
# plot all particles created during whole simulation
for x,y,f in allstates:
    plot([x,x+0.02*cos(f)],[y,y+0.02*sin(f)],'r-')
#
# plot resampled particles
for x,y,f in states:
    plot([x,x+0.02*cos(f)],[y,y+0.02*sin(f)],'k-')
#
# estimated position
estim_pos = zip(*estim_pos)
x, y = estim_pos[0], estim_pos[1]
plot(x,y,'o-')
#
# true position
pos = zip(*pos)
x, y = pos[0], pos[1]
plot(x,y,'o-')
#
# show the miracle
print 'Done!'
show()