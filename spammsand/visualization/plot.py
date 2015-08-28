from mayavi import mlab
from mayavi import tools

@mlab.show
def main(vtk_file):
    plot = tools.pipeline.open(vtk_file)
    tools.pipeline.glyph(plot, scale_factor=0.95, mode='cube')
    mlab.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "VTK",
        help="The VTK data file")
    options = parser.parse_args()
    main(options.VTK)
