# Using nMOLDYN in the Python 3 era

nMOLDYN has been developed over many years. It was originally written
in Python version 2.2, based on libraries that went back to the early
days of the scientific Python ecosystem in the 1990s. Porting these
libraries to Python 3 would have represented a major effort for which
we could not gather the required resources. After [Python 2 was sunset
in 2020](https://www.python.org/doc/sunset-python-2/), it became
increasingly difficult to install and use on modern operating systems.
While in 2022 it is still possible to follow the original nMOLDYN
installation instructions (see README), starting from source code, it
is probable that this will become more and more cumbersome.

We therefore recommended installing nMOLDYN via the [Guix package
manager](https://guix.gnu.org/), which allows using conflicting
versions of libraries in parallel and has excellent support for
reproducible computations. [A special Guix
channel](https://gitlab.inria.fr/guix-hpc/guix-past) provides much of
the Python 2 ecosystem and makes it easy to install on any modern
Linux system. We have added nMOLDYN to this collection.

To install nMOLDYN via Guix, you must be running the GNU/Linux
operating system. The steps are:

  1. Install Guix as described in [its manual](https://guix.gnu.org/en/manual/devel/en/guix.html#Binary-Installation)
  2. Create the file `~/.config/guix/channels.scm` containing the following lines:
```
(use-modules (guix ci))

(list (channel
        (name 'nmoldyn)
        (url "https://github.com/khinsen/guix-past")
        (branch "nmoldyn"))
      (channel-with-substitutes-available
          %default-guix-channel
          "https://ci.guix.gnu.org"))
```
  2. Run `guix pull` to update Guix' package definitions.
  3. Run `guix install nmoldyn`.

You can then run `nMOLDYNStart.py`.
