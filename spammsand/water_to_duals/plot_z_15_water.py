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
vtk_file_reader2 = engine.open(u'/home/matcha/Desktop/RESEARCH/spammsand_may_21_2015/spammsand/water_to_duals/z_14.vtk')

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
iso_surface2.contour.contours[0:1] = [0.01]

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


scene.scene.camera.position = [16403.433503518248, 16131.551555667136, 24041.649958155223]
scene.scene.camera.focal_point = [4609.0, 4609.0, 4609.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.53662464383091613, -0.54176353604187599, 0.64693605761987416]
scene.scene.camera.clipping_range = [9848.773230415045, 45230.169453467308]


#scene.scene.camera.position = [5760.0510962263688, 8264.1602192001847, 8166.9237003172129]
#scene.scene.camera.focal_point = [1545.0000000000143, 1544.9999999999961, 1544.9999999999879]
#scene.scene.camera.view_angle = 30.0
#scene.scene.camera.view_up = [-0.31733327151558843, -0.55718686597616895, 0.7673606656409151]
#scene.scene.camera.clipping_range = [5008.7568762710653, 17059.25614172122]

scene.scene.camera.compute_view_plane_normal()
scene.scene.render()

scene.scene.save(u'/home/matcha/Desktop/RESEARCH/spammsand_may_21_2015/spammsand/stabilized_paper_1/z_water_to_duals_scn1.png')
#!,size=(1024,1024))


scene.scene.camera.position = [5138.678723236244, 5354.3588072568082, 30077.993327509143]
scene.scene.camera.focal_point = [4609.0, 4609.0, 4609.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.72145710637153859, -0.69156169586610861, 0.035242935133142715]
scene.scene.camera.clipping_range = [15607.760545286075, 37980.790949069298]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()

scene.scene.save(u'/home/matcha/Desktop/RESEARCH/spammsand_may_21_2015/spammsand/stabilized_paper_1/z_water_to_duals_scn2.png')

