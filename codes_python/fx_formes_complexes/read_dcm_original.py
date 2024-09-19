#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 20:33:44 2021

@author: sarvi
"""
import numpy as np
from scipy.spatial import ConvexHull, convex_hull_plot_2d

# rng = np.random.default_rng()
# points = rng.random((30, 2))   # 30 random points in 2-D
# hull = ConvexHull(points)

import matplotlib.pyplot as plt
# plt.plot(points[:,0], points[:,1], 'o')
# for simplex in hull.simplices:
#     plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

import os
import pydicom
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import math
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from collections import Counter
from shapely.geometry import point, Polygon

planPath = 'Tic/'
planName = 'RP.1.3.6.1.4.1.33868.20210210162030.730598'
plan = pydicom.read_file(planPath+planName, force=True)

dummyX, dummyY, centroid = [], [], []
polypoint = []
beamSequence = plan[(0x300a, 0x00b0)].value
beamContours = []
for i in range(len(beamSequence)):
    blockSequence = beamSequence[i][(0x300a, 0x00f4)].value
    blockData = blockSequence[0][(0x300a, 0x0106)].value
    dummy = []
    for entry in blockData:
        dummy.append(entry)
    dummyArray = np.array(dummy).astype(float)
    dummyX.append(dummyArray[::2])
    dummyX[i] = dummyX[i].reshape(dummyX[i].shape[0], 1)
    dummyY.append(dummyArray[1::2])
    dummyY[i] = dummyY[i].reshape(dummyY[i].shape[0], 1)
    dummyX = dummyX
    dummyY = dummyY
    centroid.append((sum(dummyX[i][:, 0])/len(dummyX[i][:, 0]),
                    sum(dummyY[i][:, 0])/len(dummyY[i][:, 0])))
    beamContours.append(np.concatenate((dummyX[i], dummyY[i]), axis=1))
    polypoint.append(list(zip(list(dummyX[i][:, 0]), list(dummyY[i][:, 0]))))
    poly = Polygon(polypoint[i])
    
#@@@@@@@@@@@

counter = 6
planes_list=[]
beamNumber,isoCenter=[],[]
for beamdirection in plan[(0x300a,0x00b0)].value:
    beamNumber.append(beamdirection[(0x300a,0x00c0)].value)
    ctrlpoint=beamdirection[(0X300A,0X0111)].value[0]
    ctrlpoint[(0x300a,0x012c)].value[0]=round(ctrlpoint[(0x300a,0x012c)].value[0],3)
    ctrlpoint[(0x300a,0x012c)].value[1]=round(ctrlpoint[(0x300a,0x012c)].value[1],3)
    ctrlpoint[(0x300a,0x012c)].value[2]=round(ctrlpoint[(0x300a,0x012c)].value[2],3)
    isoCenter.append(ctrlpoint[0x300a,0x012c].value)
    
isoCenterPenMRT=[float(isoCenter[0][0]),float(isoCenter[0][2])*-1,float(isoCenter[0][1])*-1]
print('isoCenter PenMRT=', isoCenterPenMRT)
print('isoCenter IsoGray=', isoCenter[0])

points = beamContours[1][:-1]
def is_col(pt1, pt2, pt3):
    [a,b], [m,n], [x,y] = pt1, pt2, pt3
    return (n-b)*(x-m)==(y-n)*(m-a)

def is_col_planes(plane1, plane2):
    return is_col(plane1[0], plane1[1], plane2[0]) and is_col(plane1[0], plane1[1], plane2[1])

hull = ConvexHull(points, qhull_options='Qc')
hull.simplices = np.array([np.sort(simplex)
                           if hull.vertices.max() not in simplex
                           or hull.vertices.min() not in simplex else
                           np.sort(simplex)[::-1]
                           for simplex in hull.simplices])
hull.simplices = hull.simplices[np.argsort(hull.simplices[:, 0])]
plt.plot(points[:, 0], points[:, 1], 'o')
for simplex in hull.simplices:
    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

class Concave():
    def __init__(self, hull):
        self.hull = hull
        self.children = []
        self.plan_indices = []


def build_concave(hull: ConvexHull, parent_plane = None):
    global counter
    global planes_list
    concave = Concave(hull)

    for simplex in hull.simplices[:-1]:
        is_coplanar = False
        simplex_plane = hull.points[simplex]
        for i, plane in enumerate(planes_list):
            if is_col_planes(simplex_plane, plane):
                is_coplanar = True
                concave.plan_indices.append(i)
                print('co_p detected')
                continue
                
        if not is_coplanar:
            planes_list.append(simplex_plane)
            counter += 1
            concave.plan_indices.append(counter)

    if not parent_plane:
        counter += 1
        concave.plan_indices.append(counter)
    else:
        concave.plan_indices.append(parent_plane)

    if hull.coplanar.size != 0:
        outer_points = np.concatenate((hull.vertices, hull.coplanar[:,0]))
        print('outers:',outer_points)
    else:
        outer_points = hull.vertices

    if len(outer_points) == len(hull.points):
        return concave

    first_hull_point = min(outer_points)
    inner_pts = [i for i in range(hull.npoints) if i not in outer_points]
    print('inners:',inner_pts)
    while inner_pts:

        first_inner_pt = min(inner_pts)
        print("first inner pt:", first_inner_pt)
        last_inner_pt = min([i for i in outer_points if i >
                            first_inner_pt], default=first_hull_point)
        print("last inner pt:", last_inner_pt)

        if last_inner_pt >= first_inner_pt:
            sub_points = hull.points[first_inner_pt - 1:last_inner_pt+1]
            inner_pts = [i for i in range(
                hull.npoints) if i not in outer_points and i > last_inner_pt]
        else:
            print("1st hull pt:", first_hull_point)
            sub_points = np.concatenate(
                (hull.points[first_inner_pt-1:], hull.points[:first_hull_point+1]))
            inner_pts = []
        print("sub points =", sub_points)
        sub_hull = ConvexHull(sub_points, qhull_options='Qc')
        sub_hull.simplices = np.array([np.sort(simplex)
                                       if sub_hull.vertices.max() not in simplex
                                       or sub_hull.vertices.min() not in simplex else
                                       np.sort(simplex)[::-1]
                                       for simplex in sub_hull.simplices])
        sub_hull.simplices = sub_hull.simplices[np.argsort(
            sub_hull.simplices[:, 0])]

        concave.children.append(build_concave(sub_hull, concave.plan_indices[list(hull.simplices[:,0]).index(first_inner_pt - 1)]))

    return concave


concave = build_concave(hull)


for k in range(len(concave.children)):
    for j in range(concave.children[k].hull.points.shape[0]):
        if j != concave.children[k].hull.points.shape[0]-1:
            if concave.children[k].hull.points[j, 1] == concave.children[k].hull.points[j+1, 1]:
                np.delete(concave.children[k].hull.points, j, 0)


def plot_concave(concave):
    points = concave.hull.points
    plt.plot(points[:, 0], points[:, 1], 'x')
    for simplex in concave.hull.simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
    for child in concave.children:
        points = child.hull.points
        for simplex in child.hull.simplices:
            plt.plot(points[simplex, 0], points[simplex, 1], 'r-')
        for grandchild in child.children:
            points = grandchild.hull.points
            for simplex in grandchild.hull.simplices:
                plt.plot(points[simplex, 0], points[simplex, 1], 'b--')


plot_concave(concave)


def angle_of_line(x1, y1, x2, y2):
    return math.degrees(math.atan2(y2-y1, x2-x1))


angle2 = [angle_of_line(simplex[0, 0], simplex[0, 1], simplex[1, 0], simplex[1, 1])
          for simplex in hull.points[hull.simplices]]


geofile = open("mask.geo", 'w', encoding='utf-8')
geofile.write(str('0000000000000000000000000000000000000000000000000000000000000000')
              + '\n' + str('SURFACE (  13)    Plane Y=  -0.75')
              + '\n'+str('INDICES=( 0, 0, 0, 0, 0)')+'\n' +
              str('     AY=( 1.000000000000000E+00,   0)')
              + '\n'+str('     A0=( 7.500000000000000E-01,   0)')+'\n'
              + str('0000000000000000000000000000000000000000000000000000000000000000')
              + '\n' + str('SURFACE (  14)    Plane Y=  0.75')
              + '\n'+str('INDICES=( 0, 0, 0, 0, 0)')+'\n' +
              str('     AY=( 1.000000000000000E+00,   0)')
              + '\n'+str('     A0=(-7.500000000000000E-01,   0)')+'\n'
              + str('0000000000000000000000000000000000000000000000000000000000000000')
              + '\n' + str('SURFACE (  15)    Plane X= -6')
              + '\n'+str('INDICES=( 0, 0, 0, 0, 0)')+'\n' +
              str('     AX=( 1.000000000000000E+00,   0)')
              + '\n'+str('     A0=( 6.000000000000000E+00,   0)')+'\n'
              + str('0000000000000000000000000000000000000000000000000000000000000000')
              + '\n' + str('SURFACE (  16)    Plane X=  6')
              + '\n'+str('INDICES=( 0, 0, 0, 0, 0)')+'\n' +
              str('     AX=( 1.000000000000000E+00,   0)')
              + '\n'+str('     A0=(-6.000000000000000E+00,   0)')+'\n'
              + str('0000000000000000000000000000000000000000000000000000000000000000')
              + '\n' + str('SURFACE (  17)    Plane Z= -6')
              + '\n'+str('INDICES=( 0, 0, 0, 0, 0)')+'\n' +
              str('     AZ=( 1.000000000000000E+00,   0)')
              + '\n'+str('     A0=( 6.000000000000000E+00,   0)')+'\n'
              + str('0000000000000000000000000000000000000000000000000000000000000000')
              + '\n' + str('SURFACE (  18)    Plane Z=  6')
              + '\n'+str('INDICES=( 0, 0, 0, 0, 0)')+'\n' +
              str('     AZ=( 1.000000000000000E+00,   0)')
              + '\n'+str('     A0=(-6.000000000000000E+00,   0)')+'\n')

def write_planes(concave, top = False):
    for j, simplex in enumerate(concave.hull.points[concave.hull.simplices]):
        if top or j < len(concave.hull.simplices) - 1:
            angle = angle_of_line(simplex[0, 0], simplex[0, 1], simplex[1, 0], simplex[1, 1])
            geofile.write(str('0000000000000000000000000000000000000000000000000000000000000000')
                          + '\n' +'SURFACE ({:4d})'.format(concave.plan_indices[j]+12)+
                          str('  Plane X,Z=   ') + str(
                              simplex[0, 0])+str(',')+str(simplex[0, 1])
                          + '\n'+str('INDICES=( 0, 0, 0, 1, 0)')+'\n'+str(
                              '  THETA=(')+str('{0:22.15E}'.format(-1*angle))
                          + str(',   0)')
                          + '\n'+str('X-SHIFT=(') + str('{0:22.15E}'.format(
                              round((simplex[0, 0]/10), 15)))
                          + str(',   0)')
                          + '\n'+str('Z-SHIFT=(') + str('{0:22.15E}'.format(
                              round((simplex[0, 1]/10), 15)))
                          + str(',   0)')+'\n')
    for child in concave.children:
        write_planes(child)
        
write_planes(concave, True)


print("Producing the modules")
module_counter = 1
def write_modules(concave: Concave, material, top = False):
    global module_counter
    children_modules = []
    for child in concave.children:
        children_modules.append(write_modules(child, 3 - material))
    geofile.write('0000000000000000000000000000000000000000000000000000000000000000\n'
                  + 'MODULE  ({:4d})\n'.format(module_counter+1)
                  + 'MATERIAL({:4d})\n'.format(material))
    for i, surface_index in enumerate(concave.plan_indices):
        if material == 2:
            side_pointer =-1
        else:
            side_pointer = 1
        if i == len(concave.plan_indices) - 1 and not top:
            side_pointer *= -1
        geofile.write('SURFACE ({:4d}), SIDE POINTER=({:+d})\n'.format(surface_index+12, side_pointer))
    geofile.write('SURFACE (  13), SIDE POINTER=(+1)\n'
                  +'SURFACE (  14), SIDE POINTER=(-1)\n')
    for child_module in children_modules:
        geofile.write('MODULE  ({:4d})\n'.format(child_module+1))
    module_counter+=1
    return module_counter - 1

write_modules(concave, 1, True)


geofile.write('0000000000000000000000000000000000000000000000000000000000000000\n'
              + 'MODULE  ({:4d})\n'.format(module_counter+1)
              + 'MATERIAL(   2)\n'
              + 'SURFACE (  13), SIDE POINTER=(+1)\n'
              + 'SURFACE (  14), SIDE POINTER=(-1)\n'
              + 'SURFACE (  15), SIDE POINTER=(+1)\n'
              + 'SURFACE (  16), SIDE POINTER=(-1)\n'
              + 'SURFACE (  17), SIDE POINTER=(+1)\n'
              + 'SURFACE (  18), SIDE POINTER=(-1)\n'
              + 'MODULE  ({:4d})\n'.format(module_counter - 0)
              +'0000000000000000000000000000000000000000000000000000000000000000\n'
              +'END     00000000000000000000000000000000000000000000000000000000\n')

geofile.close()

import re

# Open the file in read mode
with open('mask.geo', 'r') as file:
    # Read the content
    content = file.read()

# Define a function to replace MATERIAL numbers based on conditions
def replace_material(match):
    material_number = int(match.group(1))
    if material_number == 1:
        return 'MATERIAL(   2)'
    elif material_number == 2:
        return 'MATERIAL(   4)'

# Use regular expression to find and replace MATERIAL numbers
content = re.sub(r'MATERIAL\((   \d+)\)', replace_material, content)

# Open the file in write mode and overwrite the content
with open('mask.geo', 'w') as file:
    file.write(content)
