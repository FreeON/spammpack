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

#vtk_file_reader  = engine.open(u'/home/m/spammsand/spammsand/64/z_0.vtk')
#iso_surface = IsoSurface()
#engine.add_module(iso_surface, obj=None)
#iso_surface.actor.mapper.scalar_mode = 'use_field_data'
#iso_surface.actor.property.specular_color = (1.0, 0.0, 0.0)
#iso_surface.actor.property.diffuse_color  = (1.0, 0.0, 0.0)
#iso_surface.actor.property.ambient_color  = (1.0, 0.0, 0.0)
#iso_surface.actor.property.color          = (1.0, 0.0, 0.0)
#iso_surface.contour.contours[0:1]         = [0.01]



#vtk_file_reader1 = engine.open(u'/home/m/spammsand/spammsand/64/z_6.vtk')
#vtk_file_reader1 = engine.open(u'/home/matcha/Desktop/RESEARCH/spammsand_may_10_2015/spammsand/350_6311gss/z_6.vtk')
#iso_surface1 = IsoSurface()
#engine.add_module(iso_surface1, obj=None)

#iso_surface1.actor.mapper.scalar_mode = 'use_field_data'
#iso_surface1.actor.property.specular_color = (1.0, 0.0, 0.0)
#iso_surface1.actor.property.diffuse_color  = (1.0, 0.0, 0.0)
#iso_surface1.actor.property.ambient_color  = (1.0, 0.0, 0.0)
#iso_surface1.actor.property.color          = (1.0, 0.0, 0.0)
#iso_surface1.contour.contours[0:1]         = [0.01]
#iso_surface1.actor.property.opacity = 1.0



#vtk_file_reader2 = engine.open(u'/home/m/spammsand/spammsand/64/z_11.vtk')
vtk_file_reader2 = engine.open(u'/home/matcha/Desktop/RESEARCH/spammsand_may_10_2015/spammsand/350_6311gss/z_12.vtk')
iso_surface2 = IsoSurface()
engine.add_module(iso_surface2, obj=None)
iso_surface2.actor.mapper.scalar_mode = 'use_field_data'
iso_surface2.actor.property.specular_color = (0.0, 1.0, 0.0)
iso_surface2.actor.property.diffuse_color  = (0.0, 1.0, 0.0)
iso_surface2.actor.property.ambient_color  = (0.0, 1.0, 0.0)
iso_surface2.actor.property.color          = (0.0, 1.0, 0.0)
iso_surface2.actor.property.opacity = .20
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
outline1.actor.property.line_width = 2.
outline1.actor.property.line_width = 2.

scene.scene.camera.position = [11196.156762432653, 7613.1210378714022, 12209.119939168879]
scene.scene.camera.focal_point = [3447.3178710937505, 3449.5173339843691, 3447.4838867187505]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [-0.64427143769366901, -0.29008441174688149, 0.707647757456772]
scene.scene.camera.clipping_range = [824.16725951679109, 27002.18051922306]
scene.scene.camera.compute_view_plane_normal()

scene.scene.render()
scene.scene.save(u'/home/m/spammsand/spammsand/snapshot.png')

scene.scene.disable_render = True

