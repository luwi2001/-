# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 21:30:16 2023

@author: MICROSOFT
"""

import math

class Point3D:
    def __init__(self, x=0, y=0, z=0, ix=0, iy=0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.ix = int(ix)
        self.iy = int(iy)

def dist(p1, p2):
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    return math.sqrt(dx * dx + dy * dy)

def midpoint(p1, p2):
    return Point3D(0.5*p1.x + 0.5*p2.x, 0.5*p1.y + 0.5*p2.y, 0.5*p1.z + 0.5*p2.z, 0, 0)

def _centroid(p1, p2, p3):
    return Point3D((p1.x + p2.x + p3.x) / 3.0, (p1.y + p2.y + p3.y) / 3.0, (p1.z + p2.z + p3.z) / 3.0, 0, 0)

class Triangle:
    def __init__(self, pt1=0, pt2=0, pt3=0):
        self.pt1 = pt1
        self.pt2 = pt2
        self.pt3 = pt3

class Mesh:
    def __init__(self, nodes, nodeArray, triangles, width, height):
        self.nodes = nodes
        self.nodeArray = nodeArray
        self.triangles = triangles
        self.width = width
        self.height = height

def triangle_area(mesh, index):
    triangle = mesh.triangles[index]
    pt1 = mesh.nodes[triangle.pt1]
    pt2 = mesh.nodes[triangle.pt2]
    pt3 = mesh.nodes[triangle.pt3]
    return _triangle_area(pt1, pt2, pt3)

def _triangle_area(p1, p2, p3):
    a = dist(p1, p2)
    b = dist(p2, p3)
    c = dist(p3, p1)
    s = (a + b + c) / 2.0
    return math.sqrt(s * (s - a) * (s - b) * (s - c))



