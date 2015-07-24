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

@mlab.show
def plot(filename, number_bins=100):
    """Plot the cubes from a file.

    The cubes are stratify them into number_bins norm bins. The transparency
    of the cubes is set depending on which norm bin the cube is in.
    """

    parser = re.compile("^\s*([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9.eEdD+-]+)$")

    x = []
    y = []
    z = []
    width = []
    norm = []

    fd = open(filename)
    for line in fd:
        result = parser.search(line)
        x.append(int(result.group(1)))
        y.append(int(result.group(2)))
        z.append(int(result.group(3)))
        width.append(int(result.group(4)))
        norm.append(float(result.group(5)))
    fd.close()

    print("loaded {:d} cubes".format(len(x)))

    # Stratify cubes by norm.
    min_norm = numpy.amin(norm)
    max_norm = numpy.amax(norm)

    def bound(i):
        return min_norm+i*(max_norm-min_norm)/float(number_bins)

    print("norms in interval [{:1.2f}, {:1.2f}]".format(min_norm, max_norm))

    for i in range(number_bins):
        print("bin {:d} in [{:1.2f}, {:1.2f})".format(
            i, bound(i), bound(i+1)))

    x_stratified = [ [] for i in range(number_bins)]
    y_stratified = [ [] for i in range(number_bins)]
    z_stratified = [ [] for i in range(number_bins)]
    norm_stratified = [ [] for i in range(number_bins)]

    print("stratifying into {:d} bins".format(number_bins))
    for i in range(len(x)):
        found_bin = False
        for j in range(number_bins):
            if norm[i] >= bound(j) and norm[i] < bound(j+1):
                x_stratified[j].append(x[i])
                y_stratified[j].append(y[i])
                z_stratified[j].append(z[i])
                norm_stratified[j].append(norm[i])
                found_bin = True
                break
        if not found_bin:
            x_stratified[number_bins-1].append(x[i])
            y_stratified[number_bins-1].append(y[i])
            z_stratified[number_bins-1].append(z[i])
            norm_stratified[j].append(norm[i])
            print("norm {:1.2f} outside".format(norm[i]))

    # Get the current figure.
    figure = mlab.gcf()

    # Get the engine.
    engine = mlab.get_engine()

    # Clean the figure.
    mlab.clf()

    # Turn off rendering (for performance).
    figure.scene.disable_render = True

    # Tune background color.
    #mlab.figure(bgcolor=(1, 1, 1))

    # Add cubes.
    for i in range(number_bins):
        print("{:d} cubes with norm [{:e},{:e})".format(
            len(x_stratified[i]), bound(i), bound(i+1)))
        if len(x_stratified[i]) > 0:
            points = mlab.points3d(x_stratified[i],
                                   y_stratified[i],
                                   z_stratified[i],
                                   mode='cube',
                                   colormap='gist_heat',
                                   color=mycolor((i+1)/float(number_bins)),
                                   scale_factor=width[0],
                                   opacity=(i+1)/float(number_bins))
            # lut = points.module_manager.scalar_lut_manager.lut.table.to_array()
            # for j in range(lut.shape[0]):
            #     print(lut[j,:])

    # Add axes.
    from mayavi.modules.axes import Axes
    axes = Axes()
    engine.add_module(axes, obj=None)

    xmax = max(numpy.amax(x), numpy.amax(y), numpy.amax(z))+width[0]/2-2
    print("xmax = {:d}".format(xmax))

    axes.axes.x_label = 'i'
    axes.axes.y_label = 'j'
    axes.axes.z_label = 'k'
    axes.axes.label_format = '%-3.0f'
    axes.property.display_location = 'background'
    axes.axes.ranges = [1, xmax,
                        1, xmax,
                        1, xmax]

    # Box around the whole thing.
    mlab.outline(extent=[1, xmax,
                         1, xmax,
                         1, xmax])

    # Turn rendering back on.
    figure.scene.disable_render = False
