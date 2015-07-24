from mayavi import mlab
import numpy
import re

@mlab.show
def plot(filename, number_bins=6):
    """Plot the cubes from a file.

    The cubes are stratify them into number_of_bins norm bins. The
    transparency of the cubes is set depending on which norm bin the
    cube is in.
    """

    parser = re.compile("^\s*([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9.eEdD+-]+)$")

    x = []
    y = []
    z = []
    norm = []

    fd = open(filename)
    for line in fd:
        result = parser.search(line)
        x.append(int(result.group(1)))
        y.append(int(result.group(2)))
        z.append(int(result.group(3)))
        norm.append(float(result.group(4)))
    fd.close()

    print("loaded {:d} cubes".format(len(x)))

    # Stratify cubes by norm.
    min_norm = numpy.amin(norm)
    max_norm = numpy.amax(norm)

    x_stratified = [ [] for i in range(number_bins)]
    y_stratified = [ [] for i in range(number_bins)]
    z_stratified = [ [] for i in range(number_bins)]

    def bound(i):
        return min_norm+i*(max_norm-min_norm)

    for i in range(len(x)):
        for j in range(number_bins):
            if norm[i] >= bound(j) and norm[i] < bound(j+1):
                x_stratified[j].append(x[i])
                y_stratified[j].append(y[i])
                z_stratified[j].append(z[i])
                break

    # Get the current figure.
    f = mlab.gcf()

    # Clean the figure.
    mlab.clf()

    # Turn off rendering (for performance).
    f.scene.disable_render = True

    # Add cubes.
    for i in range(number_bins):
        if len(x_stratified[i]) > 0:
            mlab.points3d(x_stratified[i],
                          y_stratified[i],
                          z_stratified[i],
                          mode='cube',
                          opacity=(i+1)/float(number_bins))

    # Turn rendering back on.
    f.scene.disable_render = False
