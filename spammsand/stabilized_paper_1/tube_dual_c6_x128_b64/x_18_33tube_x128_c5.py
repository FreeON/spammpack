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


#vtk_file_reader2 = engine.open(u'/home/m/spammsand/spammsand/water_to_duals/z_15.vtk')

vtk_file_reader2 = engine.open(u'x_19.vtk')

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

iso_surface2.actor.property.opacity = 0.4
iso_surface2.contour.contours[0:1] = [0.03]

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
axes.axes.z_label = 'k'
#axes.axes.label_format = ''
axes.property.display_location = 'background'

scene.scene.isometric_view()

camera_light = engine.scenes[0].scene.light_manager.lights[0]
camera_light.activate  = True
camera_light.azimuth   = -20.0
camera_light.elevation = 35.0 
#camera_light.color = (1.0, 0.0, 1.0)
camera_light.color = (1.0, 1.0, 0.0)

camera_light1 = engine.scenes[0].scene.light_manager.lights[1]
camera_light1.activate = False
camera_light2= engine.scenes[0].scene.light_manager.lights[2]
camera_light2.activate = False

camera_light3 = engine.scenes[0].scene.light_manager.lights[3]
camera_light3.activate = True
camera_light3.elevation = -2.0
camera_light3.azimuth = -10.0
camera_light3.color = (0.0, 1.0, 0.0)

scene.scene.camera.position = [19371.507490174925, 19371.507490174925, 19371.507490174925]
scene.scene.camera.focal_point = [6913.0, 6913.0, 6913.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0, 0.0, 1.0]
scene.scene.camera.clipping_range = [51.892147109552056, 51892.147109552054]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()

scene.scene.save(u'x_19_scene1.png')
#,size=(1024,1024))

scene.scene.camera.position = [9690.9622653846054, 21883.988892956138, 20095.618502138386]
scene.scene.camera.focal_point = [6730.4512952534214, 4402.5403496259396, 5728.5184666368496]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.033795933371882338, -0.63114808074832296, 0.77492576099600119]
scene.scene.camera.clipping_range = [46.840525535698553, 46840.525535698551]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()

scene.scene.save(u'x_19_scene2.png')
