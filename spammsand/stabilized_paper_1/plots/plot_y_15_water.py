# Recorded script from Mayavi2
from numpy import array
try:
    engine = mayavi.engine
except NameError:
    from mayavi.api import Engine
    engine = Engine()
    engine.start()
if len(engine.scenes) == 0:
    engine.new_scene()
# ------------------------------------------- 

from mayavi.modules.iso_surface import IsoSurface


vtk_file_reader2 = engine.open(u'/home/m/spammsand/spammsand/water_6311gss_duals/y_15.vtk')
##vtk_file_reader2 = engine.open(u'/home/matcha/Desktop/RESEARCH/spammsand_may_10_2015/spammsand/350_6311gss/y_12.vtk')
iso_surface2 = IsoSurface()
engine.add_module(iso_surface2, obj=None)
iso_surface2.actor.mapper.scalar_mode = 'use_field_data'

#iso_surface2.actor.property.specular_color = (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
#iso_surface2.actor.property.diffuse_color = (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
#iso_surface2.actor.property.ambient_color = (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
#iso_surface2.actor.property.color = (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)

iso_surface2.actor.property.specular_color = (1.0, 1.0, 1.0)
iso_surface2.actor.property.diffuse_color  = (1.0, 1.0, 1.0)
iso_surface2.actor.property.ambient_color  = (1.0, 1.0, 1.0)
iso_surface2.actor.property.color          = (1.0, 1.0, 1.0)

iso_surface2.actor.property.opacity = 0.3
iso_surface2.contour.contours[0:1] = [0.0001]

scene = engine.scenes[0]

#from mayavi.modules.axes import Axes
#axes = Axes()
#engine.add_module(axes, obj=None)

from mayavi.modules.outline import Outline
outline1 = Outline()
engine.add_module(outline1, obj=None)

outline1.actor.mapper.scalar_range = array([ 0.,  1.])
outline1.outline_mode = 'full'
outline1.actor.property.specular_color = (0.0, 0.0, 0.0)
outline1.actor.property.diffuse_color = (0.0, 0.0, 0.0)
outline1.actor.property.ambient_color = (0.0, 0.0, 0.0)
outline1.actor.property.color = (0.0, 0.0, 0.0)
outline1.actor.property.line_width = 4.
outline1.actor.property.line_width = 4.

#scene.scene.background = (0.7529411764705882, 0.7529411764705882, 0.7529411764705882)
scene.scene.background = (1.0, 1.0, 1.0)
scene.scene.jpeg_quality = 100

from mayavi.modules.axes import Axes
axes = Axes()
engine.add_module(axes, obj=None)

axes.axes.x_label = 'i'
axes.axes.y_label = 'j'
axes.axes.y_label = 'k'
axes.axes.label_format = ''
axes.property.display_location = 'background'


scene.scene.isometric_view()

camera_light = engine.scenes[0].scene.light_manager.lights[0]
camera_light.activate  = True
camera_light.azimuth   = -20.0
camera_light.elevation = 35.0 
camera_light.color = (1.0, 0.0, 1.0)

camera_light1 = engine.scenes[0].scene.light_manager.lights[1]
camera_light1.activate = False
camera_light2= engine.scenes[0].scene.light_manager.lights[2]
camera_light2.activate = False

camera_light3 = engine.scenes[0].scene.light_manager.lights[3]
camera_light3.activate = True
camera_light3.elevation = -2.0
camera_light3.azimuth = -10.0
camera_light3.color = (0.0, 1.0, 0.0)


scene.scene.camera.position = [5760.0510962263688, 8264.1602192001847, 8166.9237003172129]
scene.scene.camera.focal_point = [1545.0000000000143, 1544.9999999999961, 1544.9999999999879]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.31733327151558843, -0.55718686597616895, 0.7673606656409151]
scene.scene.camera.clipping_range = [5008.7568762710653, 17059.25614172122]

scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
#scene.scene.save(u'/home/m/spammsand/spammsand/y_15_tube_5_36.png',size=(512,512))
scene.scene.save(u'/home/m/spammsand/spammsand/y_15_water.png',size=(1024,1024))



