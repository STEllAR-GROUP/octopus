# Copyright (c) 2012 Zach Byerly 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

from os.path import join

OpenDatabase(join(DATA_DIRECTORY, "3d_torus.silo")) #opens data file

HideToolbars()
InvertBackgroundColor() #sets background to black
Source(join(VISUALIZATION_DIRECTORY, "clear_annotation.py")) #clears annotations

AddPlot("Contour", "rho_level_0", 1, 1)
AddPlot("Mesh", "mesh_level_0", 1, 1)
SetActivePlots(1) #making the mesh the active plot
AddOperator("Slice", 0) #slicing the mesh
SliceAtts = SliceAttributes()
SliceAtts.normal = (0,0,1) #slice along the z-plane
SliceAtts.axisType = SliceAtts.ZAxis 
SliceAtts.project2d = 0 #don't project to 2D
SliceAtts.meshName = "mesh_level_0"
SetOperatorOptions(SliceAtts,0)

SetActivePlots(0) #making the density the active plot
AddOperator("Transform", 0) #raising the density above the mesh 
TransformAtts = TransformAttributes()
TransformAtts.doTranslate = 1
TransformAtts.translateZ = 5e-05 #raising it by this much
SetOperatorOptions(TransformAtts, 0)

#setting the contours
ContourAtts = ContourAttributes()
ContourAtts.contourNLevels = 3
SetPlotOptions(ContourAtts)

AddOperator("Reflect", 0) #reflecting because we use z-reflect in the code
ReflectAtts = ReflectAttributes() 
ReflectAtts.reflections = (1, 0, 0, 0, 1, 0, 0, 0)
SetOperatorOptions(ReflectAtts, 0)

AddOperator("Clip", 0)#clipping planes
ClipAtts = ClipAttributes() 
ClipAtts.plane1Status = 1
ClipAtts.plane2Status = 1
ClipAtts.plane1Normal = (-1, 0, 0)
ClipAtts.plane2Normal = (0, -1, 0)
SetOperatorOptions(ClipAtts, 0)

DrawPlots() 
Source(join(VISUALIZATION_DIRECTORY, "set_view.py")) #sets view according to file

MeshAtts = MeshAttributes()
MeshAtts.meshColor = (102, 102, 153, 255) #makes mesh purple
MeshAtts.meshColorSource = MeshAtts.MeshCustom
SetPlotOptions(MeshAtts)


