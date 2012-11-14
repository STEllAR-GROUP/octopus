from os.path import join

OpenDatabase(join(DATA_DIRECTORY, "3d_torus.silo")) #opens data file

HideToolbars()
#InvertBackgroundColor() #sets background to black
Source(join(VISUALIZATION_DIRECTORY, "clear_annotation.py")) #clears annotations


AddPlot("Pseudocolor", "rho_level_0", 1, 1)
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.minFlag = 1
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.min = 0
PseudocolorAtts.max = 14
PseudocolorAtts.colorTableName = "orangehot"
SetPlotOptions(PseudocolorAtts)

AddOperator("Isovolume", 1)
IsovolumeAtts = IsovolumeAttributes()
IsovolumeAtts.lbound = 1
IsovolumeAtts.ubound = 1e+37
IsovolumeAtts.variable = "default"
SetOperatorOptions(IsovolumeAtts, 1)

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
#Source("/home/zbyerly/SC12/scripts/set_view.py") #sets view according to file
Source(join(VISUALIZATION_DIRECTORY, "set_view.py")) #sets view according to file
Source(join(VISUALIZATION_DIRECTORY, "set_view_angle.py")) #sets view according to file
Source(join(VISUALIZATION_DIRECTORY, "set_zoom.py")) #sets view according to file
set_view_angle(70,135)
set_zoom(2.5)

MeshAtts = MeshAttributes()
MeshAtts.meshColor = (102, 102, 153, 255) #makes mesh purple
MeshAtts.meshColorSource = MeshAtts.MeshCustom
SetPlotOptions(MeshAtts)
