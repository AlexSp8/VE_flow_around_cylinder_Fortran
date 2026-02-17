#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###


#sys module: for system operations
import sys
#salome module: control salome
import salome

#initiate salome
salome.salome_init()
#import notebook
import salome_notebook
#add study to notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/aspyridakis/Projects/Stability/2D_Cylinder/2D_VE_Cylinder_T_full')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS

#geomBuilder: class which will use the study to create geompy object
geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

# geompy.addToStudy(O, 'O')
# geompy.addToStudy(OX, 'OX')
# geompy.addToStudy(OY, 'OY')
# geompy.addToStudy(OZ, 'OZ')

M_ref = 0.125*2#*8*3   #Mesh density
Tetra = False
Quadratic = False
Domain_Mult = 1

print('M = ', M_ref)

BR = 0.5
L  = 2*(1.0/BR)  #x
H  = 2*(1.0/BR)  #y
W = 0.5*Domain_Mult/2   #z

#outer box
L1 = 12.5*(H/2.0)
L2 = 12.5*(H/2.0)
H2 = H
# W_outer = 100.0

print('BR = ', BR)
print('H = ', H/2.0, '/ L = ', L1, L2, '/ W = ', W)

e  = 0.0001

Lines  = []
Groups = []
IDs    = []

r = [1.0]

#make Face_1 (box): first x, then y
Face_1 = geompy.MakeFaceHW(L, H, 1)
# geompy.addToStudy(Face_1, 'Face_1')

#make Face_2 (channel): first x, then y
Face_2 = geompy.MakeFaceHW(L1+L2, H, 1)
#For unequal Lengths
if (abs(L1-L2) > e):
  Face_2 = geompy.MakeTranslation(Face_2, L1, 0, 0)
# geompy.addToStudy(Face_2, 'Face_2')

# #make Face_3 (channel): first x, then y
# Face_3 = geompy.MakeFaceHW(L1+L2, H2, 1)
# geompy.addToStudy(Face_3, 'Face_3')

#make disks
ndisks = 1           #number of disks
Disk = geompy.MakeDiskR(1.0, 1)
# geompy.addToStudy(Disk, 'Disk')

angle = math.atan(H/L) #in rad
geompy.Rotate(Disk, OZ, angle)
# geompy.addToStudy(Disk, 'Disk')

#################################################################################
#make diagonal lines
#Line 1
vertex_1 = geompy.MakeVertex(L/2, H/2,  0.0)
vertex_2 = geompy.MakeVertex(-L/2, -H/2, 0.0)
Line = geompy.MakeLineTwoPnt(vertex_1, vertex_2)
Lines.append(Line)

#Line 2
vertex_1 = geompy.MakeVertex(L/2, -H/2, 0.0)
vertex_2 = geompy.MakeVertex(-L/2, H/2, 0.0)
Line = geompy.MakeLineTwoPnt(vertex_1, vertex_2)
Lines.append(Line)

#Line 3
vertex_1 = geompy.MakeVertex(-L1, 0, 0.0)
vertex_2 = geompy.MakeVertex(+L2, 0, 0.0)
Line = geompy.MakeLineTwoPnt(vertex_1, vertex_2)
Lines.append(Line)

#Line 4
vertex_1 = geompy.MakeVertex(0, H/2, 0.0)
vertex_2 = geompy.MakeVertex(0, -H/2, 0.0)
Line = geompy.MakeLineTwoPnt(vertex_1, vertex_2)
Lines.append(Line)

#################################################################################
#partition Face, disks, connecting lines, vertices
Partition_1 = geompy.MakePartition([Face_1, Face_2, Disk], [], [], [], 
  geompy.ShapeType["FACE"], 0, [], 0)

for i in range(0,len(Lines)):
  Partition_1 = geompy.MakePartition([Partition_1, Lines[i]], [], [], [], 
    geompy.ShapeType["FACE"], 0, [], 0)

# geompy.addToStudy(Partition_1, 'Partition_1')
# Face_4 = geompy.MakeTranslation(Face_2, 0, -H/2, z)
# geompy.addToStudy(Face_4, 'Face_4')
# Cut_1 = geompy.MakeCutList(Partition_1, [Face_4, Disk], True)
Cut_1 = geompy.MakeCutList(Partition_1, [Disk], True)
# geompy.addToStudy(Cut_1, 'Cut_1')

Partition_1 = geompy.MakeTranslation(Cut_1, 0, 0, +W)
# Extrusion_1 = geompy.MakePrismVecH(Partition_1, OZ, -2*W)
# Partition_1 = Extrusion_1
geompy.addToStudy(Partition_1, 'Partition_1')
########################################################################
# make groups in Partition_1
j = 1  #group counter

r2  = [1.0, L/2]
#order: right, middle, left
#1
Group = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
for i in range(0,len(r2)):
  #up
  #right
  vertex = geompy.MakeVertex(+e, r2[i], W)
  #get sub-shape (edge)
  edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
  #edge ID
  ID = geompy.GetSubShapeID(Partition_1, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])

  #middle
  vertex = geompy.MakeVertex(-e, r2[i], W)
  #get sub-shape (edge)
  edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
  #edge ID
  ID = geompy.GetSubShapeID(Partition_1, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])

  #left
  vertex = geompy.MakeVertex(-r2[i], +e, W)
  #get sub-shape (edge)
  edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
  #edge ID
  ID = geompy.GetSubShapeID(Partition_1, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])

  #down
  #right
  vertex = geompy.MakeVertex(+e, -r2[i], W)
  #get sub-shape (edge)
  edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
  #edge ID
  ID = geompy.GetSubShapeID(Partition_1, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])

  #middle
  vertex = geompy.MakeVertex(-e, -r2[i], W)
  #get sub-shape (edge)
  edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
  #edge ID
  ID = geompy.GetSubShapeID(Partition_1, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])

  #left
  vertex = geompy.MakeVertex(-r2[i], -e, W)
  #get sub-shape (edge)
  edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
  #edge ID
  ID = geompy.GetSubShapeID(Partition_1, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])

#up left
vertex = geompy.MakeVertex(-L1, +e, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
#print(ID)
geompy.UnionIDs(Group, [ID])

#down left
vertex = geompy.MakeVertex(-L1, -e, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
#print(ID)
geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(Partition_1, Group, 'Group_'+str(j))
Groups.append(Group)
j = j+1
########################################################################
#2
Group = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
for i in range(0,len(r2)):
  #up right
  vertex = geompy.MakeVertex(r2[i], +e, W)
  #get sub-shape (edge)
  edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
  #edge ID
  ID = geompy.GetSubShapeID(Partition_1, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])

  #down right
  vertex = geompy.MakeVertex(r2[i], -e, W)
  #get sub-shape (edge)
  edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
  #edge ID
  ID = geompy.GetSubShapeID(Partition_1, edge)
  IDs.append(ID)
  geompy.UnionIDs(Group, [ID])

#up right out
vertex = geompy.MakeVertex(L2, +e, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
#print(ID)
geompy.UnionIDs(Group, [ID])

#down right out
vertex = geompy.MakeVertex(L2, -e, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
#print(ID)
geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(Partition_1, Group, 'Group_'+str(j))
Groups.append(Group)
j = j+1

########################################################################
#3: channel right
Group = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
#up
vertex = geompy.MakeVertex(+L/2+e, +H/2, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
geompy.UnionIDs(Group, [ID])

#middle
vertex = geompy.MakeVertex(+L/2+e, 0, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
geompy.UnionIDs(Group, [ID])

#down
vertex = geompy.MakeVertex(+L/2+e, -H/2, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(Partition_1, Group, 'Group_'+str(j))
Groups.append(Group)
j = j+1

#4: channel left
Group = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
#up
vertex = geompy.MakeVertex(-L/2-e, +H/2, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
geompy.UnionIDs(Group, [ID])

#middle
vertex = geompy.MakeVertex(-L/2-e, 0, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
geompy.UnionIDs(Group, [ID])

#down
vertex = geompy.MakeVertex(-L/2-e, -H/2, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(Partition_1, Group, 'Group_'+str(j))
Groups.append(Group)
j = j+1
  
########################################################################
#5: diagonal line
angle = math.atan(H/L) #in rad

Group = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
x = r2[0]*math.cos(angle)
y = r2[0]*math.sin(angle)
vertex = geompy.MakeVertex(x+e, y+e, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
#print(ID)
geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(Partition_1, Group, 'Group_'+str(j))
Groups.append(Group)
j = j+1

########################################################################
#Inflow 1D
Group = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])

vertex = geompy.MakeVertex(-L1, e, +W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
geompy.UnionIDs(Group, [ID])

vertex = geompy.MakeVertex(-L1, -e, +W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(Partition_1, Group, 'Inflow_1D')
Groups.append(Group)
j = j+1

########################################################################
#Outflow 1D
Group = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])

vertex = geompy.MakeVertex(+L1, e, +W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
geompy.UnionIDs(Group, [ID])

vertex = geompy.MakeVertex(+L1, -e, +W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID = geompy.GetSubShapeID(Partition_1, edge)
IDs.append(ID)
geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(Partition_1, Group, 'Outflow_1D')
Groups.append(Group)
j = j+1

########################################################################
#Cylinder
Group = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
x = [1.0,    +e,    -e, -1.0, -1.0,     -e,     +e, 1.0]
y = [ +e, 1.0-e, 1.0-e,   +e,   -e, -1.0+e, -1.0+e,  -e]

for i in range(0,len(x)):
  #up
  vertex = geompy.MakeVertex(x[i], y[i], +W)
  face = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
  # geompy.addToStudyInFather(Partition_1, face, 'face')
  ID = geompy.GetSubShapeID(Partition_1, face)
  geompy.UnionIDs(Group, [ID])

  #down
  vertex = geompy.MakeVertex(x[i], -y[i], +W)
  face = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
  # geompy.addToStudyInFather(Partition_1, face, 'face')
  ID = geompy.GetSubShapeID(Partition_1, face)
  geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(Partition_1, Group, 'Cylinder')
Groups.append(Group)
j = j+1

########################################################################
#Wall faces
x = [-L1/2, -e, +e, +L1/2]
y = [+H/2, -H/2]
name = ['Top_wall', 'Bottom_wall']
for k in range(len(y)):
  Group = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
  for i in range(len(x)):
    vertex = geompy.MakeVertex(x[i], y[k], +W)
    #get sub-shape (edge)
    edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
    #edge ID
    ID = geompy.GetSubShapeID(Partition_1, edge)
    IDs.append(ID)
    geompy.UnionIDs(Group, [ID])

  geompy.addToStudyInFather(Partition_1, Group, name[k])
  Groups.append(Group)
  j = j+1

########################################################################
#Cylinder_2D: front
x = [-L1/2, -L/4-e,    -e,    +e, +L/4+e, +L1/2]
y = [   +e,     +e, 2.0-e, 2.0-e,     +e,    +e]
Group = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
for k in range(1,3):
  for i in range(len(x)):
    yi = y[i]*((-1)**k)
    vertex = geompy.MakeVertex(x[i], yi, +W)
    # geompy.addToStudyInFather(Partition_1, vertex, 'vertex'+str(i))
    #get sub-shape (edge)
    edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["FACE"])
    #edge ID
    ID = geompy.GetSubShapeID(Partition_1, edge)
    IDs.append(ID)
    geompy.UnionIDs(Group, [ID])

geompy.addToStudyInFather(Partition_1, Group, 'Front')
Groups.append(Group)
j = j+1

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish(False) 
# Set to False to avoid publish in study if not needed or in some particular situations:
# multiples meshes built in parallel, complex and numerous mesh edition (performance)

Sub_meshes = []

Mesh_1 = smesh.Mesh(Partition_1,'Mesh_1')
# Mesh_1 = smesh.Mesh(Partition_1)
# smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')

Regular_1D = Mesh_1.Segment()
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')

#Default number of segments if it is left unspecified somewhere
Number_of_Segments_1 = Regular_1D.NumberOfSegments(5)
smesh.SetName(Number_of_Segments_1, 'Number_of_Segments_1')

Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
# Quadrangle_2D = Mesh_1.Quadrangle()
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
Quadrangle_Parameters_1 = Quadrangle_2D.QuadrangleParameters(smeshBuilder.QUAD_TRIANGLE_PREF,-1,[],[])
smesh.SetName(Quadrangle_Parameters_1, 'Quadrangle Parameters_1')

# for i in range(0,len(Groups)):
#   Group_ = Mesh_1.GroupOnGeom(Groups[i],'Group_'+str(i+1))
#   #smesh.SetName(Group_, 'Group_'+str(i+1))
#   # Group_ = Mesh_1.GetGroups()

size = len(Groups)
Inflow_1D = Mesh_1.GroupOnGeom(Groups[size-6],'Inflow_1D')
Outflow_1D = Mesh_1.GroupOnGeom(Groups[size-5],'Outflow_1D')
Cylinder_1D = Mesh_1.GroupOnGeom(Groups[size-4],'Cylinder')
Top_wall_1D = Mesh_1.GroupOnGeom(Groups[size-3],'Top_wall')
Bottom_wall_1D = Mesh_1.GroupOnGeom(Groups[size-2],'Bottom_wall')
Front = Mesh_1.GroupOnGeom(Groups[size-1],'Front')

########################################################################
if BR == 0.9:
  Factors = [1.0, 0.2, 800.0, 1.0/800.0, 1.0/1.0, 1.0]
elif BR == 0.8:
  Factors = [1.0, 0.2, 250.0, 1.0/250.0, 1.0/2.4, 1.0]
elif BR == 0.7:
  Factors = [1.0, 0.2, 120.0, 1.0/120.0, 1.0/5.0, 1.0]
elif BR == 0.6:
  Factors = [1.0, 0.2, 70.0, 1.0/70.0, 1.0/10.0, 1.0]
elif BR == 0.5:
  Factors = [1.0, 1.0/10.0, 50.0, 1.0/50.0, 1.0/18.0, 1.0]
elif BR == 0.4:
  Factors = [1.0, 0.2, 34.0, 1.0/34.0, 1.0/30.0, 1.0]
elif BR == 0.3:
  Factors = [1.0, 0.2, 24.0, 1.0/24.0, 1.0/52.0, 1.0]
elif BR == 0.25:
  Factors = [1.0, 0.2, 22.0, 1.0/22.0, 1.0/70.0, 1.0]
elif BR == 0.2:
  Factors = [1.0, 0.2, 18.0, 1.0/18.0, 1.0/100.0, 1.0]
elif BR == 0.1:
  Factors = [1.0, 0.2, 14.0, 1.0/14.0, 1.0/250.0, 1.0]
elif BR == 0.05:
  Factors = [1.0, 0.2, 12.0, 1.0/12.0, 1.0/580.0, 1.0]
else:
  print("BR not defined")

j = 0       #Sub-mesh counter
k = 0       #scale factor

#1: upstream
Nelem = int(M_ref*(32))
print('Upstream: ', Nelem, Factors[j])
Regular_1D = Mesh_1.Segment(geom=Groups[j])
Number_of_Segments = Regular_1D.NumberOfSegments(Nelem)
smesh.SetName(Number_of_Segments, 'Number of Segments_'+str(j))
j = j+1
Sub_mesh = Regular_1D.GetSubMesh()
Sub_meshes.append(Sub_mesh)
smesh.SetName(Sub_mesh, 'Sub-mesh_'+str(j))

#2: downstream up
Nelem = int(M_ref*(80))
print('Downstream: ', Nelem, Factors[j])
Regular_1D = Mesh_1.Segment(geom=Groups[j])
Number_of_Segments = Regular_1D.NumberOfSegments(Nelem)
smesh.SetName(Number_of_Segments, 'Number of Segments_'+str(j))
Number_of_Segments.SetScaleFactor(Factors[j])
smesh.SetName(Number_of_Segments, 'Number of Scaled Segments_'+str(k))
Propagation_of_Node = Regular_1D.PropagationOfDistribution()
smesh.SetName(Propagation_of_Node, 'Propagation of Node Distribution on Opposite Edges')

#reversed edges
vertex = geompy.MakeVertex(L/2, +e, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID1 = geompy.GetSubShapeID(Partition_1, edge)

# #down cylinder
# vertex = geompy.MakeVertex(1.0, -e, W)
# #get sub-shape (edge)
# edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
# #edge ID
# ID2 = geompy.GetSubShapeID(Partition_1, edge)

#out
vertex = geompy.MakeVertex(L2, +e, W)
#get sub-shape (edge)
edge = geompy.GetShapesNearPoint(Partition_1, vertex, geompy.ShapeType["EDGE"])
#edge ID
ID3 = geompy.GetSubShapeID(Partition_1, edge)
# Number_of_Segments.SetReversedEdges([ID1, ID2, ID3])
Number_of_Segments.SetReversedEdges([ID1, ID3])

j = j+1
Sub_mesh = Regular_1D.GetSubMesh()
Sub_meshes.append(Sub_mesh)
smesh.SetName(Sub_mesh, 'Sub-mesh_'+str(j))

########################################################################
#3, 4: outer box
Nelem = [32, 32] #3D

for i in range(0,len(Nelem)):

  Nelem[i] = int(M_ref*Nelem[i])

  print('Out: ', Nelem[i], Factors[j])

  Regular_1D = Mesh_1.Segment(geom=Groups[j])
  Number_of_Segments = Regular_1D.NumberOfSegments(Nelem[i])
  Number_of_Segments.SetScaleFactor(Factors[j])
  smesh.SetName(Number_of_Segments, 'Number of Scaled Segments_'+str(k))
  j = j+1
  Propagation_of_Node = Regular_1D.PropagationOfDistribution()
  smesh.SetName(Propagation_of_Node, 'Propagation of Node Distribution on Opposite Edges')
  smesh.SetName(Number_of_Segments, 'Number of Segments_'+str(j))
  Sub_mesh = Regular_1D.GetSubMesh()
  Sub_meshes.append(Sub_mesh)
  smesh.SetName(Sub_mesh, 'Sub-mesh_'+str(j))

########################################################################
#Scale Distribution Sub-meshes
#5: diagonal
Nelem = int(M_ref*64)
print('Diag: ', Nelem, Factors[j])

# if BR == 0.9:
#   Nelem = int(0.375*Nelem)
# elif BR == 0.8:
#   Nelem = int(0.625*Nelem)
# elif BR == 0.7:
#   Nelem = int(0.75*Nelem)
# elif BR == 0.7:
#   Nelem = int(0.875*Nelem)

Regular_1D = Mesh_1.Segment(geom=Groups[j])
Number_of_Segments = Regular_1D.NumberOfSegments(Nelem)
Number_of_Segments.SetScaleFactor(Factors[j])
smesh.SetName(Number_of_Segments, 'Number of Scaled Segments_'+str(k))
Propagation_of_Node = Regular_1D.PropagationOfDistribution()
smesh.SetName(Propagation_of_Node, 'Propagation of Node Distribution on Opposite Edges')
# Number_of_Segments.SetObjectEntry('Partition_1')
Sub_mesh = Regular_1D.GetSubMesh()
Sub_meshes.append(Sub_mesh)
j = j+1
k = k+1
smesh.SetName(Sub_mesh, 'Sub-mesh_'+str(j))

#compute the mesh
isDone = Mesh_1.Compute()

# if (Quadratic):
#   # Mesh_1.ConvertToQuadratic(0) #medium node on the edge
#   Mesh_1.ConvertToQuadratic(0, Mesh_1, True) #bi-quadratic
#   # Mesh_1.ConvertToQuadratic(1) #medium node in the middle
#   # Mesh_1.ConvertToQuadratic(1, Mesh_1, True) #bi-quadratic

########################################################################
#6: z
Nelem = int(M_ref*4*8*Domain_Mult) #3D
Nelem = Nelem*2
# Nelem = 1 #2D

dw = 2*W/Nelem

print('z: ', Nelem)

Groups_extruded = Mesh_1.ExtrusionSweepObjects( [Mesh_1], [Mesh_1], [Mesh_1], 
  [0, 0, -dw], Nelem, 1, [ ], 0, [ ], [ ], 0 )

if (Tetra):
  # Mesh_1.QuadTo4Tri(Mesh_1)
  Mesh_1.SplitVolumesIntoTetra(Mesh_1, 1) #5 tetra
  # Mesh_1.SplitVolumesIntoTetra(Mesh_1, 1, 1) #avoid over-constrained volumes
  # Mesh_1.SplitVolumesIntoTetra(Mesh_1, 2) #6 tetra

Ngroups = len(Groups_extruded)

Inflow_2D = Groups_extruded[0]
Outflow_2D = Groups_extruded[1]
Cylinder_2D = Groups_extruded[2]
Top_wall_2D = Groups_extruded[3]
Bottom_wall_2D = Groups_extruded[4]
# Front = Groups_extruded[5]
Mesh_3D = Groups_extruded[5]
Back = Groups_extruded[Ngroups-1]

try:
  # Mesh_1.ExportUNV(r'/home/aspyridakis/Downloads/Cylinder_3D_MESH.unv')

  Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Cylinder_3D_MESH.DTA')
  # Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Cylinder_3D_MESH.DTA', Mesh_3D)
  Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Inflow_1D_MESH.DTA', Inflow_1D)
  # Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Outflow_1D.DTA', Outflow_1D)
  # Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Cylinder_1D.DTA', Cylinder_1D)
  # Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Top_wall_1D.DTA', Top_wall_1D)
  # Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Bottom_wall_1D.DTA', Bottom_wall_1D)
  Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Inflow_2D_MESH.DTA', Inflow_2D)
  # Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Outflow_2D.DTA', Outflow_2D)
  # Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Cylinder_2D.DTA', Cylinder_2D)
  # Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Top_wall_2D.DTA', Top_wall_2D)
  # Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Bottom_wall_2D.DTA', Bottom_wall_2D)
  Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Cylinder_2D_MESH.DTA', Front)
  # Mesh_1.ExportDAT(r'/home/aspyridakis/Downloads/Back.DTA', Back)

  # #Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Cylinder_3D_MESH.DTA')
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Cylinder_3D_MESH.DTA', Mesh_3D)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Inflow_1D_MESH.DTA', Inflow_1D)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Outflow_1D.DTA', Outflow_1D)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Cylinder_1D.DTA', Cylinder_1D)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Top_wall_1D.DTA', Top_wall_1D)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Bottom_wall_1D.DTA', Bottom_wall_1D)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Inflow_2D_MESH.DTA', Inflow_2D)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Outflow_2D.DTA', Outflow_2D)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Cylinder.DTA', Cylinder)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Top_wall.DTA', Top_wall)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Bottom_wall.DTA', Bottom_wall)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Cylinder_2D_MESH.DTA', Front)
  # Mesh_1.ExportDAT(r'C:/Users/alex_/Downloads/Back.DTA', Back)

  pass
except:
  print('ExportPartToDAT() failed. Invalid file name?')

########################################################################
if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
