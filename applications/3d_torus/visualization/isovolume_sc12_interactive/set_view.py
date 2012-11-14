# I want to turn this into a function that takes an angle from 0 to 90 degrees
# and translates that into the appropriate view coordinates.

def set_view(view):
    
   #create some control points
    c0 = View3DAttributes()
    c0.viewNormal = (0, 0, 1)
    c0.focus = (0, 0, 0)
    c0.viewUp = (0.5, 0.86, 0)
    c0.viewAngle = 30
    c0.parallelScale = 0.0004
    c0.nearPlane = -0.00075
    c0.farPlane = 0.00075
    c0.perspective = 1
 
    c1 = View3DAttributes()
    c1.viewNormal = (-0.25, -0.4, 0.9)
    c1.focus = (0, 0, 0)
    c1.viewUp = (0.5, 0.86, 0.45)
    c1.viewAngle = 30
    c1.parallelScale = 0.0004
    c1.nearPlane = -0.00075
    c1.farPlane = 0.00075
    c1.perspective = 1
 
    c2 = View3DAttributes()
    c2.viewNormal = (-0.4, -0.75, 0.53)
    c2.focus = (0, 0, 0)
    c2.viewUp = (0.26, 0.45, 0.85)
    c2.viewAngle = 30
    c2.parallelScale = 0.0004
    c2.nearPlane = -0.00075
    c2.farPlane = 0.00075
    c2.perspective = 1
 
    c3 = View3DAttributes()
    c3.viewNormal = (-0.438771, 0.523661, -0.730246)
    c3.focus = (0, 0, 0)
    c3.viewUp = (-0.0199911, 0.80676, 0.590541)
    c3.viewAngle = 30
    c3.parallelScale = 8.28257
    c3.nearPlane = 3.5905
    c3.farPlane = 48.2315
    c3.perspective = 1
 
    c4 = View3DAttributes()
    c4.viewNormal = (0.286142, -0.342802, -0.894768)
    c4.focus = (0, 0, 0)
    c4.viewUp = (-0.0382056, 0.928989, -0.36813)
    c4.viewAngle = 30
    c4.parallelScale = 10.4152
    c4.nearPlane = 1.5495
    c4.farPlane = 56.1905
    c4.perspective = 1
 
    c5 = View3DAttributes()
    c5.viewNormal = (0.974296, -0.223599, -0.0274086)
    c5.focus = (0, 0, 0)
    c5.viewUp = (0.222245, 0.97394, -0.0452541)
    c5.viewAngle = 30
    c5.parallelScale = 1.1052
    c5.nearPlane = 24.1248
    c5.farPlane = 58.7658
    c5.perspective = 1
 
    c6 = c0

    if view == 0:
        c=c0
    elif view == 1:
        c=c1
    elif view == 2:
        c=c2
    elif view == 3:
        c=c3
    elif view == 4:
        c=c4
    elif view == 5:
        c=c5
    else:
        print "ERROR views go from 0 to 5"

    SetView3D(c)

    
        
        



#View3DAtts = View3DAttributes()
#View3DAtts.viewNormal = (-0.454822, -0.852912, 0.256277)
#View3DAtts.focus = (0, 0, 0.00015)
#View3DAtts.viewUp = (0.110004, 0.231757, 0.966534)
#View3DAtts.viewAngle = 30
#View3DAtts.parallelScale = 0.000259808
#View3DAtts.nearPlane = -0.000519615
#View3DAtts.farPlane = 0.000519615
#View3DAtts.imagePan = (-0.0158611, 0.296352)
#View3DAtts.imageZoom = 1.21
#View3DAtts.perspective = 1
#View3DAtts.eyeAngle = 2
#View3DAtts.centerOfRotationSet = 0
#View3DAtts.centerOfRotation = (0, 0, 0.00015)
#View3DAtts.axis3DScaleFlag = 0
#View3DAtts.axis3DScales = (1, 1, 1)
#View3DAtts.shear = (0, 0, 1)
#SetView3D(View3DAtts)
