#!/usr/bin/env python

def main ():
    """
    The main function.
    """

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
            "SOURCE",
            help = "The source file for the unit test",
            nargs = "+"
            )

    parser.add_argument(
            "--output",
            help = "The output file")

    parser.add_argument(
            "--spammpack-lib",
            help = "The spammpack library to link against",
            choices = [
                "spammpack_serial_shared",
                "spammpack_serial_static"
                ],
            default = "spammpack_serial_shared"
            )

    options = parser.parse_args()

    if options.output:
        fd = open(options.output, "w")
    else:
        import sys
        fd = sys.stdout
        print("writing to %s" % (options.output))

    for source in options.SOURCE:
        import os.path
        testbasename = os.path.splitext(os.path.basename(source))[0]
        testexename = "unit-test-" + testbasename
        fd.write("# Unit test from %s\n" % (source))
        fd.write("add_executable( %s %s )\n" % (testexename, source))
        fd.write("target_link_libraries( %s %s spammtests)\n" % (testexename, options.spammpack_lib))
        fd.write("target_include_directories( %s PRIVATE ${CMAKE_BINARY_DIR}/src )\n" % (testexename))
        fd.write("add_test( %s %s )\n" % (testbasename, testexename))
        fd.write("\n")
        fd.write("# Unit test using valgrind from %s\n" % (source))
        fd.write("if( VALGRIND )\n")
        fd.write("  add_test( valgrind-%s ${VALGRIND} --error-exitcode=1 ${CMAKE_CURRENT_BINARY_DIR}/%s )\n" % (testbasename, testexename))
        fd.write("endif()\n")

    fd.close()

if __name__ == "__main__":
    main()
