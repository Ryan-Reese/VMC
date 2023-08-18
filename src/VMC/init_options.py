"""Module containing the initialisation options for PyMol.

The options are chosen to assist various functions in VMC.
All can be changed to suit your needs, but the final `undo_disable` option is highly, highly recommended.

"""

from pymol import cmd

def init_options() -> None:

    # sets background colour to grey
    # assists in visualising the community colours
    cmd.bg_color("grey")  # Can be changed

    # gives ray-traced images an opaque (default: grey) background
    cmd.set("ray_opaque_background", "on")

    # puts PyMol in maximum-quality mode
    cmd.util.performance(0)

    # Disabling undo greatly reduces memory cost.
    # PyMol `cmd.set` operations apparently have a very high memory overhead.
    # From rough tests, this option will increase the rate of atom/bond colouring by ~100x.
    cmd.undo_disable()

