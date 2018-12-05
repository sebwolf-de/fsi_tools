#! /usr/bin/env python

import math

n = 64  
for k in range(0,n):
    print('Point('+str(k+1)+') = { '+str(math.cos(2*math.pi*k/n))+', '+str(math.sin(2*math.pi*k/n))+', 0, 1.0};')

for k in range(0,n):
    print('Line('+str(k+1)+') = {'+str(k+1)+', '+str((k+1)%n+1)+'};')

s = 'Line Loop(1) = {'
for k in range(0,n-1):
    s += str(k+1)+', '
s += str(n) + '};'
print(s)

print('Plane Surface(1) = {1};')
