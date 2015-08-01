#    $ ipython --gui=wx
#    In [1]: %run visualization/plot-2.py
#    In [2]: plot("z_dual_1.norms")


from mayavi import mlab
import numpy
import re

def mycolor(x):
    """Returns a color vector (a triple of floats) based on x, where x is
    in the range [0, 1].

    """

    lut = [
        [  0,   0,   0, 255],
        [  1,   0,   0, 255],
        [  2,   0,   0, 255],
        [  4,   0,   0, 255],
        [  5,   0,   0, 255],
        [  6,   0,   0, 255],
        [  8,   0,   0, 255],
        [  9,   0,   0, 255],
        [ 10,   0,   0, 255],
        [ 12,   0,   0, 255],
        [ 14,   0,   0, 255],
        [ 16,   0,   0, 255],
        [ 17,   0,   0, 255],
        [ 18,   0,   0, 255],
        [ 20,   0,   0, 255],
        [ 21,   0,   0, 255],
        [ 23,   0,   0, 255],
        [ 24,   0,   0, 255],
        [ 26,   0,   0, 255],
        [ 27,   0,   0, 255],
        [ 28,   0,   0, 255],
        [ 29,   0,   0, 255],
        [ 31,   0,   0, 255],
        [ 32,   0,   0, 255],
        [ 33,   0,   0, 255],
        [ 35,   0,   0, 255],
        [ 36,   0,   0, 255],
        [ 37,   0,   0, 255],
        [ 39,   0,   0, 255],
        [ 40,   0,   0, 255],
        [ 42,   0,   0, 255],
        [ 43,   0,   0, 255],
        [ 46,   0,   0, 255],
        [ 47,   0,   0, 255],
        [ 48,   0,   0, 255],
        [ 50,   0,   0, 255],
        [ 51,   0,   0, 255],
        [ 53,   0,   0, 255],
        [ 54,   0,   0, 255],
        [ 55,   0,   0, 255],
        [ 56,   0,   0, 255],
        [ 58,   0,   0, 255],
        [ 59,   0,   0, 255],
        [ 60,   0,   0, 255],
        [ 62,   0,   0, 255],
        [ 63,   0,   0, 255],
        [ 65,   0,   0, 255],
        [ 66,   0,   0, 255],
        [ 68,   0,   0, 255],
        [ 69,   0,   0, 255],
        [ 70,   0,   0, 255],
        [ 71,   0,   0, 255],
        [ 73,   0,   0, 255],
        [ 74,   0,   0, 255],
        [ 77,   0,   0, 255],
        [ 78,   0,   0, 255],
        [ 80,   0,   0, 255],
        [ 81,   0,   0, 255],
        [ 82,   0,   0, 255],
        [ 84,   0,   0, 255],
        [ 85,   0,   0, 255],
        [ 86,   0,   0, 255],
        [ 88,   0,   0, 255],
        [ 89,   0,   0, 255],
        [ 91,   0,   0, 255],
        [ 93,   0,   0, 255],
        [ 95,   0,   0, 255],
        [ 96,   0,   0, 255],
        [ 97,   0,   0, 255],
        [ 98,   0,   0, 255],
        [100,   0,   0, 255],
        [101,   0,   0, 255],
        [102,   0,   0, 255],
        [104,   0,   0, 255],
        [105,   0,   0, 255],
        [108,   0,   0, 255],
        [110,   0,   0, 255],
        [111,   0,   0, 255],
        [113,   0,   0, 255],
        [114,   0,   0, 255],
        [115,   0,   0, 255],
        [116,   0,   0, 255],
        [118,   0,   0, 255],
        [119,   0,   0, 255],
        [120,   0,   0, 255],
        [122,   0,   0, 255],
        [123,   0,   0, 255],
        [124,   0,   0, 255],
        [126,   0,   0, 255],
        [127,   0,   0, 255],
        [128,   0,   0, 255],
        [130,   0,   0, 255],
        [131,   0,   0, 255],
        [133,   0,   0, 255],
        [134,   0,   0, 255],
        [135,   0,   0, 255],
        [138,   0,   0, 255],
        [140,   0,   0, 255],
        [140,   0,   0, 255],
        [142,   0,   0, 255],
        [143,   0,   0, 255],
        [145,   0,   0, 255],
        [146,   0,   0, 255],
        [147,   0,   0, 255],
        [149,   0,   0, 255],
        [150,   0,   0, 255],
        [152,   0,   0, 255],
        [153,   0,   0, 255],
        [155,   0,   0, 255],
        [156,   0,   0, 255],
        [157,   0,   0, 255],
        [158,   0,   0, 255],
        [160,   0,   0, 255],
        [161,   0,   0, 255],
        [162,   0,   0, 255],
        [164,   0,   0, 255],
        [165,   0,   0, 255],
        [167,   0,   0, 255],
        [169,   0,   0, 255],
        [170,   0,   0, 255],
        [172,   0,   0, 255],
        [173,   0,   0, 255],
        [175,   1,   0, 255],
        [176,   3,   0, 255],
        [177,   4,   0, 255],
        [179,   6,   0, 255],
        [180,   8,   0, 255],
        [182,  10,   0, 255],
        [183,  13,   0, 255],
        [185,  16,   0, 255],
        [187,  17,   0, 255],
        [188,  19,   0, 255],
        [189,  20,   0, 255],
        [191,  22,   0, 255],
        [192,  24,   0, 255],
        [194,  26,   0, 255],
        [195,  28,   0, 255],
        [197,  30,   0, 255],
        [198,  32,   0, 255],
        [200,  34,   0, 255],
        [202,  36,   0, 255],
        [203,  38,   0, 255],
        [205,  40,   0, 255],
        [206,  42,   0, 255],
        [207,  44,   0, 255],
        [209,  46,   0, 255],
        [210,  48,   0, 255],
        [211,  49,   0, 255],
        [212,  51,   0, 255],
        [214,  52,   0, 255],
        [215,  54,   0, 255],
        [217,  56,   0, 255],
        [218,  58,   0, 255],
        [220,  60,   0, 255],
        [221,  61,   0, 255],
        [222,  63,   0, 255],
        [224,  65,   0, 255],
        [225,  67,   0, 255],
        [226,  68,   0, 255],
        [227,  70,   0, 255],
        [229,  72,   0, 255],
        [232,  76,   0, 255],
        [233,  77,   0, 255],
        [234,  79,   0, 255],
        [236,  81,   0, 255],
        [237,  83,   0, 255],
        [239,  85,   0, 255],
        [240,  86,   0, 255],
        [241,  88,   0, 255],
        [242,  89,   0, 255],
        [244,  91,   0, 255],
        [245,  93,   0, 255],
        [247,  95,   0, 255],
        [248,  97,   0, 255],
        [249,  99,   0, 255],
        [251, 101,   0, 255],
        [252, 102,   0, 255],
        [253, 103,   0, 255],
        [255, 105,   0, 255],
        [255, 107,   0, 255],
        [255, 109,   0, 255],
        [255, 111,   0, 255],
        [255, 114,   0, 255],
        [255, 117,   0, 255],
        [255, 118,   0, 255],
        [255, 120,   0, 255],
        [255, 121,   0, 255],
        [255, 123,   0, 255],
        [255, 125,   0, 255],
        [255, 127,   0, 255],
        [255, 129,   0, 255],
        [255, 131,   0, 255],
        [255, 133,   1, 255],
        [255, 136,   8, 255],
        [255, 137,  11, 255],
        [255, 139,  15, 255],
        [255, 141,  19, 255],
        [255, 143,  22, 255],
        [255, 145,  26, 255],
        [255, 146,  30, 255],
        [255, 148,  34, 255],
        [255, 150,  37, 255],
        [255, 152,  41, 255],
        [255, 154,  47, 255],
        [255, 157,  52, 255],
        [255, 159,  55, 255],
        [255, 161,  59, 255],
        [255, 162,  63, 255],
        [255, 164,  67, 255],
        [255, 166,  70, 255],
        [255, 168,  74, 255],
        [255, 170,  78, 255],
        [255, 171,  81, 255],
        [255, 173,  85, 255],
        [255, 174,  89, 255],
        [255, 176,  93, 255],
        [255, 178,  96, 255],
        [255, 180, 100, 255],
        [255, 182, 103, 255],
        [255, 184, 107, 255],
        [255, 186, 110, 255],
        [255, 187, 114, 255],
        [255, 188, 118, 255],
        [255, 190, 122, 255],
        [255, 192, 126, 255],
        [255, 196, 133, 255],
        [255, 198, 137, 255],
        [255, 200, 140, 255],
        [255, 202, 144, 255],
        [255, 203, 148, 255],
        [255, 205, 152, 255],
        [255, 206, 155, 255],
        [255, 208, 158, 255],
        [255, 210, 162, 255],
        [255, 212, 166, 255],
        [255, 214, 169, 255],
        [255, 216, 173, 255],
        [255, 217, 177, 255],
        [255, 219, 181, 255],
        [255, 221, 184, 255],
        [255, 222, 188, 255],
        [255, 224, 192, 255],
        [255, 226, 195, 255],
        [255, 228, 199, 255],
        [255, 229, 203, 255],
        [255, 231, 206, 255],
        [255, 234, 212, 255],
        [255, 237, 217, 255],
        [255, 238, 221, 255],
        [255, 240, 225, 255],
        [255, 242, 228, 255],
        [255, 244, 232, 255],
        [255, 245, 236, 255],
        [255, 247, 240, 255],
        [255, 249, 243, 255],
        [255, 251, 247, 255]
    ]

    if(x < 0 or x > 1):
        raise Exception("illegal scale")

    i = int(x*255)
    return (lut[i][0]/255., lut[i][1]/255., lut[i][2]/255.)

def stratify(number_bins, norms, *args):
    """Stratify the lists by norms into number_bins strata.

    The args contain the centers and the widths, i.e. args <- center_i,
    center_j, width_i, width_j.  The function returns a tuple (norms,
    centers) of the stratified result.
    """

    length = len(norms)
    for i in range(len(args)):
        if len(args[i]) != length:
            print(norms)
            print(i)
            print(args[i])
            raise Exception("All lengths have to match")

    min_norm = numpy.amin(norms)
    max_norm = numpy.amax(norms)

    def bound(i):
        return min_norm+i*(max_norm-min_norm)/float(number_bins)

    args_stratified = [ [ [] for j in range(number_bins) ] for i in range(len(args)) ]
    norms_stratified = [ [] for i in range(number_bins) ]

    print("stratifying into {:d} bins".format(number_bins))
    for i in range(len(norms)):
        found_bin = False
        for j in range(number_bins):
            if norms[i] >= bound(j) and norms[i] < bound(j+1):
                for k in range(len(args)):
                    args_stratified[k][j].append(args[k][i])
                norms_stratified[j].append(norms[i])
                found_bin = True
                break
        if not found_bin:
            for k in range(len(args)):
                args_stratified[k][number_bins-1].append(args[k][i])
            norms_stratified[j].append(norms[i])

    # for i in range(number_bins):
    #     print("{:d} norm [{:1.2f},{:1.2f})".format(
    #         len(args_stratified[0][i]), bound(i), bound(i+1)))

    result = [ norms_stratified ]
    for arg in args_stratified:
        result.append(arg)

    return result

def read_squares(fd, start, end=None):
    """Reads a norms file from a call to spamm_tree_print_leaves_2d_symm().
    """

    # Use readline() with a length argument so we can tell whether the
    # file as reached EOF.
    LINE_LENGTH = 1000

    i = []
    j = []
    width_i = []
    width_j = []
    norm = []

    re_matrix_square = re.compile("^\s*([0-9.eEdD+-]+)"
                                  + "\s+([0-9.eEdD+-]+)"
                                  + "\s+([0-9]+)"
                                  + "\s+([0-9]+)"
                                  + "\s+([0-9.eEdD+-]+)$")

    while True:
        line = fd.readline(LINE_LENGTH)
        if len(line) == 0:
            return None
        if start.search(line):
            matrix_name = line.rstrip()
            break

    line = fd.readline()
    block_size = int(line)

    while True:
        old_position = fd.tell()
        line = fd.readline(LINE_LENGTH)
        if len(line) == 0:
            break
        if end != None:
            if end.search(line):
                fd.seek(old_position)
                break
        result = re_matrix_square.search(line)
        i.append(float(result.group(1)))
        j.append(float(result.group(2)))
        width_i.append(int(result.group(3)))
        width_j.append(int(result.group(4)))
        norm.append(float(result.group(5)))

    print("loaded {:d} matrix squares from {:s}".format(len(i), matrix_name))
    result = (block_size, i, j, width_i, width_j, norm)
    #print(result)
    return result

def read_cubes(fd, start, end=None):
    """Reads a norms file from a call to spamm_tree_print_leaves_2d_symm().
    """

    # Use readline() with a length argument so we can tell whether the
    # file as reached EOF.
    LINE_LENGTH = 1000

    i = []
    j = []
    k = []
    width_i = []
    width_j = []
    width_k = []
    norm = []

    re_product_cube = re.compile("^\s*([0-9.eEdD+-]+)"
                                  + "\s+([0-9.eEdD+-]+)"
                                  + "\s+([0-9.eEdD+-]+)"
                                  + "\s+([0-9]+)"
                                  + "\s+([0-9]+)"
                                  + "\s+([0-9]+)"
                                  + "\s+([0-9.eEdD+-]+)$")

    while True:
        line = fd.readline(LINE_LENGTH)
        if len(line) == 0:
            return None
        if start.search(line):
            matrix_name = line.rstrip()
            break

    line = fd.readline()
    block_size = int(line)

    while True:
        old_position = fd.tell()
        line = fd.readline(LINE_LENGTH)
        if len(line) == 0:
            break
        if end != None:
            if end.search(line):
                fd.seek(old_position)
                break
        result = re_product_cube.search(line)
        i.append(float(result.group(1)))
        j.append(float(result.group(2)))
        k.append(float(result.group(3)))
        width_i.append(int(result.group(4)))
        width_j.append(int(result.group(5)))
        width_k.append(int(result.group(6)))
        norm.append(float(result.group(7)))

    print("loaded {:d} product cubes from {:s}".format(len(i), matrix_name))
    return (i, j, k, width_i, width_j, width_k, norm)

@mlab.show
def plot(filename, number_bins=2):
    """Plot the cubes from a file.

    The cubes are stratified into number_bins norm bins. The
    transparency of the cubes is set depending on which norm bin the
    cube is in.
    """

    re_matrix_A = re.compile("^\s*Matrix A$")
    re_matrix_B = re.compile("^\s*Matrix B$")
    re_matrix_C = re.compile("^\s*Matrix C$")
    re_product_space = re.compile("^\s*Product Space$")

    fd = open(filename)
    (block_size, A_i, A_j,
     A_width_i, A_width_j, A_norm) = read_squares(fd, re_matrix_A, end=re_matrix_B)
    (block_size, B_i, B_j,
     B_width_i, B_width_j, B_norm) = read_squares(fd, re_matrix_B, end=re_matrix_C)
    (block_size, C_i, C_j,
     C_width_i, C_width_j, C_norm) = read_squares(fd, re_matrix_C, end=re_product_space)

    (prod_i, prod_j, prod_k,
     prod_width_i, prod_width_j, prod_width_k, prod_norm) = read_cubes(fd, re_product_space)

    # Get the current figure.
    figure = mlab.gcf()

    # Get the engine.
    engine = mlab.get_engine()

    # Clean the figure.
    mlab.clf()

    # Turn off rendering (for performance).
    figure.scene.disable_render = True

    # Tune background color.
    figure.scene.background = (1, 1, 1)

    # Stratify matrix squares.
    (norms_stratified,
     A_i_stratified, A_j_stratified,
     A_width_i_stratified, A_width_j_stratified) = stratify(number_bins, A_norm,
                                                            A_i, A_j, A_width_i, A_width_j)

    # Add matrices.
    print("Plotting matrix A")
    for i in range(number_bins):
        if len(A_i_stratified[i]) > 0:
            points = mlab.points3d(A_i_stratified[i],
                                   [1 for j in range(len(A_i_stratified[i]))],
                                   A_j_stratified[i],
                                  mode='cube',
                                 color=(0.0, 0.5019607843137255, 0.5019607843137255),
                                scale_factor=1,
                                   opacity=0.5*(i+1)/float(number_bins))

            points.glyph.glyph_source.glyph_source.x_length = block_size
            points.glyph.glyph_source.glyph_source.y_length = 0
            points.glyph.glyph_source.glyph_source.z_length = block_size

    (norms_stratified,
     B_i_stratified, B_j_stratified,
     B_width_i_stratified, B_width_j_stratified) = stratify(number_bins, B_norm,
                                                            B_i, B_j, B_width_i, B_width_j)

    # Add matrices.
    print("Plotting matrix B")
    for i in range(number_bins):
        if len(B_i_stratified[i]) > 0:
            points = mlab.points3d([1 for j in range(len(B_i_stratified[i]))],
                                   B_j_stratified[i],
                                   B_i_stratified[i],
                                   mode='cube',
                                   color=(0.5019607843137255, 0.0, 0.0),
                                   scale_factor=1,
                                   opacity=0.5*(i+1)/float(number_bins))
            points.glyph.glyph_source.glyph_source.x_length = 0
            points.glyph.glyph_source.glyph_source.y_length = block_size
            points.glyph.glyph_source.glyph_source.z_length = block_size

    (norms_stratified,
     C_i_stratified, C_j_stratified,
     C_width_i_stratified, C_width_j_stratified) = stratify(number_bins, C_norm,
                                                            C_i, C_j, C_width_i, C_width_j)

    # Add matrices.
    print("Plotting matrix C")
    for i in range(number_bins):
        if len(C_i_stratified[i]) > 0:
            points = mlab.points3d(C_i_stratified[i],
                                   C_j_stratified[i],
                                   [1 for j in range(len(C_i_stratified[i]))],
                                   mode='cube',
                                   color=(0.5019607843137255, 0.0, 0.5019607843137255),
                                   scale_factor=1,
                                   opacity=0.5*(i+1)/float(number_bins))
            points.glyph.glyph_source.glyph_source.x_length = block_size
            points.glyph.glyph_source.glyph_source.y_length = block_size
            points.glyph.glyph_source.glyph_source.z_length = 0

    # Stratify cubes by norm.
    (norms_stratified, prod_i_stratified, prod_j_stratified, prod_k_stratified) = stratify(
        number_bins, prod_norm, prod_i, prod_j, prod_k)

    # Add cubes.
    print("Plotting product cubes")
    for i in range(number_bins):
        if len(prod_i_stratified[i]) > 0:
            points = mlab.points3d(prod_i_stratified[i],
                                   prod_j_stratified[i],
                                   prod_k_stratified[i],
                                   mode='cube',
                                   color=(0.3,.3,.3),
                                   scale_factor=1,
                                   opacity=0.5*(i+1)/float(number_bins))
            points.glyph.glyph_source.glyph_source.x_length = block_size
            points.glyph.glyph_source.glyph_source.y_length = block_size
            points.glyph.glyph_source.glyph_source.z_length = block_size

#                                   color=(0.0,0.0,0.0), 

#                                   colormap='gist_heat',
#                                   color=mycolor((i+1)/float(number_bins)),

            # This is how to the colormap values:
            #
            # lut = points.module_manager.scalar_lut_manager.lut.table.to_array()
            # for j in range(lut.shape[0]):
            #     print(lut[j,:])

    i_max = max(numpy.amax(prod_i), numpy.amax(prod_j), numpy.amax(prod_k))+block_size/2
    print("i_max = {:e}".format(i_max))

    # Insert fake invisible data-set for axes.
    mlab.points3d([1, i_max], [1, i_max], [1, i_max], mode='cube', scale_factor=0)

    #mlab.axes(xlabel="i", ylabel="j", zlabel="k", extent=[1, xmax, 1, xmax, 1, xmax])

    # Box around the whole thing.
    mlab.outline(extent=[1, i_max, 1, i_max, 1, i_max])

    outline = engine.scenes[0].children[-1].children[0].children[1]
    outline.actor.property.color = (0, 0, 0)
    outline.actor.property.line_width = 2

    # Add axes.
    from mayavi.modules.axes import Axes
    axes = Axes()
    engine.add_module(axes, obj=None)

    axes.axes.label_format = '%-3.0f'
    axes.axes.width = 2
    axes.axes.x_label = 'i'
    axes.axes.y_label = 'j'
    axes.axes.z_label = 'k'
    axes.label_text_property.color = (0, 0, 0)
    axes.label_text_property.opacity = 0.0
    axes.label_text_property.shadow = True
    axes.label_text_property.shadow_offset = numpy.array([ 1, -1])
    axes.property.color = (0, 0, 0)
    axes.property.display_location = 'background'
    axes.title_text_property.color = (0, 0, 0)
    axes.title_text_property.shadow_offset = numpy.array([ 1, -1])

    # Fix camera position.
    #print(figure.scene.camera)
#    figure.scene.camera.position = [7000, 9000, 8500]
#    figure.scene.camera.focal_point = [1500, 1500, 1500]
#    figure.scene.camera.view_angle = 30.0
#    figure.scene.camera.view_up = [0, 0, 1]

    figure.scene.disable_render = False
    figure.scene.camera.compute_view_plane_normal()

    import os.path

    figure.scene.camera.position = [1975, 1975., 2400.]
    figure.scene.camera.focal_point = [440.0, 440.0, 440.0]
    figure.scene.camera.view_angle = 30.0
    figure.scene.camera.view_up = [-0.476, -0.476, 0.740]
    figure.scene.camera.clipping_range = [1406.7446099663789, 4872.2717825401005]

    png_filename = os.path.splitext(filename)[0] + "_cant1.png"
    print("Saving image to " + png_filename)
    figure.scene.save(png_filename)

    # For Y: CLOSE UP  ALONG I=K & CUBE DIAGONAL
    figure.scene.camera.position = [1878.3518173107655, 2210.072124516224, 1963.6777932457967]
    figure.scene.camera.focal_point = [904.5, 904.5, 904.5]
    figure.scene.camera.view_angle = 30.0
    figure.scene.camera.view_up = [-0.40500836669875068, -0.37455739560584284, 0.83407132806551876]
    figure.scene.camera.clipping_range = [5.8611765479361146, 5861.1765479361147]
    figure.scene.camera.compute_view_plane_normal()

    png_filename = os.path.splitext(filename)[0] + "_zoom_cube_diag_ylense.png"
    print("Saving image to " + png_filename)
    figure.scene.save(png_filename)

    # For X: 
    figure.scene.camera.position = [4725.8633422234443, 1998.1802839889251, 4736.2307754209651]
    figure.scene.camera.focal_point = [904.5, 904.5, 904.5]
    figure.scene.camera.view_angle = 30.0
    figure.scene.camera.view_up = [-0.13532282092230383, 0.980162372840087, -0.14480834577509852]
    figure.scene.camera.clipping_range = [2617.2771154002421, 9189.439101845117]
    figure.scene.camera.compute_view_plane_normal()
    figure.scene.render()

    png_filename = os.path.splitext(filename)[0] + "_xresolve.png"
    print("Saving image to " + png_filename)
    figure.scene.save(png_filename)




    figure.scene.isometric_view()

    png_filename = os.path.splitext(filename)[0] + "_isov.png"
    print("Saving image to " + png_filename)
    figure.scene.save(png_filename)


    # Turn rendering back on.
    #    figure.scene.disable_render = False

    # Save the figure to file.
    #    import os.path
    #    png_filename = os.path.splitext(filename)[0] + ".png"
    #    print("Saving image to " + png_filename)
    #    figure.scene.save(png_filename)
