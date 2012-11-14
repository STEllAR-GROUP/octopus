#function to set zoom level.
#SC12 recommended default zoom = 2.5, min=1 max=5

def set_zoom(zoom):
    c = GetView3D()
    c.imageZoom = zoom
    SetView3D(c)
