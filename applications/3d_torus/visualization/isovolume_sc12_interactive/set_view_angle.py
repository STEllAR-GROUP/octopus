# I want to turn this into a function that takes two angles from 0 to 90 degrees
# and translates that into the appropriate view coordinates.

# Phi is the azimuthal angle; theta is the polar angle.
# 0 <= phi <= 360     0 <= theta <= 180

import math

def set_view_angle(phi,theta):    
    
    while phi <= 0:
        phi = phi + 360
    while phi > 360:
        phi = phi -360
    
    while theta <= 0:
        theta = theta + 180
    while theta > 180:
        theta = theta - 180

    up_theta = theta + 90
#    while up_theta > 180:
#        up_theta = up_theta - 180    

    print "phi =", phi
    print "theta =", theta
    print "up_theta =", up_theta

    #convert deg to radians
    phi = phi*math.pi/180
    theta = theta*math.pi/180
    up_theta = up_theta*math.pi/180
    

    normal_x = -math.sin(theta)*math.cos(phi)
    normal_y = -math.sin(theta)*math.sin(phi)
    normal_z = -math.cos(theta)

    up_x = -math.sin(up_theta)*math.cos(phi) 
    up_y = -math.sin(up_theta)*math.sin(phi)
    up_z = -math.cos(up_theta) 


    c0 = GetView3D()
    c0.viewNormal = (normal_x,normal_y,normal_z)
    c0.focus = (0, 0, 0)    
    c0.viewUp = (up_x,up_y,up_z)    
    c0.viewAngle = 30
    c0.parallelScale = 0.0004
    c0.nearPlane = -0.00075
    c0.farPlane = 0.00075
    c0.perspective = 1

    SetView3D(c0)
    print "done setting view"
