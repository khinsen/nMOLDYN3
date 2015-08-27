# -*- coding: iso-8859-1 -*-

aboutAuthors = """
nMoldyn program v2.1 written by:

Tomasz Róg (email: tomekr@mol.uj.edu.pl)
Krzysztof Murzyn (email: murzyn@mol.uj.edu.pl)
Konrad Hinsen (email: hinsen@cnrs-orleans.fr)

nMoldyn project supervisor:
Gerald Kneller (email: kneller@cnrs-orleans.fr)
"""

aboutHistory = """
 HISTORY

 30.06.1999 nMoldyn v0.5
            Tomasz Róg releases the first version of
            the program. nMoldyn has fully functional
            * scattering functions * and routines for
            calculation (auto)correlation functions.
            It also have a simple GUI. * Dynamics *
            menu partially works (there are
            problems with RCF and AVACF, MSD and DOS.
            Routines in * Results * menu work and is
            almost complete. Atom selection is very
            primitive yet, group selection does not exist.


 01.12.2000 nMoldyn v1.2
            Krzysztof Murzyn works on the existing code,
            brushing it up and patching in a few places.
            The code is rewritten to form a Python site
            package and it is organised in a bunch of
            modules. The functions of crucial importance
            are grouped in * nMoldyn.core * module,
            auxilliary functions and classes in
            * nMoldyn.calc * and so on ... (for a full
            description see the program documentation).
            In * Dynamics * menu, AVACF and RCF are now
            fully functional. Konrad Hinsen fixes a bug
            in MSD. Basic changes and new features in
            this version:
            (1) rewritten GUI (xMoldyn program),
            (2) a new way to run calculation handled by
                pMoldyn program (command line switches..),
                input files are now Python programs
                executed in the nMoldyn environment,
            (3) flexible and powerful group and atom
                selection (including PDB files processing
                (reference structure, selection), graphical
                interface, selection via external user
                defined Python program),
            (4) parallelism of a few functions
                (thanks to Konrad)
            (5) new core functions: SAVACF, MPCF, DF, AT,
                RBT, FT (for abreviations see program
                documentation) and lots of a new code in
                nMoldyn.calc and nMoldyn.misc modules
 12.08.2001 nMoldyn v1.3
            A few buglets were fixed. A GUI was slightly 
            reorgenized. Many minor improvements...
"""               


logo = """\
R0lGODlheABTAMYAABoabI6OjFJWlL7K1KaytEJCRG52fLK+zHqCrCYqfEpKTJKevGZqbHZ6hKq2
vDY6hL6+vJqmvB4idC4yfD5ChF5eXEZKjF5mnP7+/LrGzIKKrHJ6pE5SjGZunI6WtMbS3KKuxKam
pFZelMbO3GpypLrC1LK6zCYmdDY2hEZGjHp+rKqyxF5inBoedMLO1LbCzIKGrC4ufJqivK66zD4+
hJ6qxDI2fEJGjE5OjGZqnHZ+rFpinHJ2pIaOtH6GrCoufJaivKq2zDo+hL7G1J6mvCImdEJChEpO
jGJqnIaKtHZ6pFJSlGpupJKavKayxFpelG52pBoadFZalMLK1La+zH6CrCoqfJaevH5+fDo6hCIi
dDIyfEpKjGJmnLrG1IKKtE5SlGZupI6WvKauxKKqxG5ypCYmfKqyzB4edMLO3LbC1DI2hIqStK62
zMbGxJ6mxEJCjHZ6rFpenP///////////////////////////////////////////////////yH+
Dk1hZGUgd2l0aCBHSU1QACwAAAAAeABTAAAH/oAYgoOEhYaHiImKi4yNjo+QkZKTlJWWl5iZmo1u
IZ6foKGio6SlpqeoqaqibocKBg2xsrO0tba3uLm6u7y9swwVri8fxMXGx8jJysvMzc7P0MYOwYYK
w8ku2drb3N3e397R4uPM064Z2x/a5OzILurw7/Lx9PP29fj3+sQjBNSFCjIUezcQXDqDCLvlW7iv
IcOH2YiZq4aO2DqCzjC2MxZxo8dnK/4RCghvoLqEFhOqNPmx2YABGRz00EBCh0ZlEwEK3JfxZruI
DoNCVKdhh5AJLQAoVTrsXsoP/oQRnHrSW4erV5lgxepEpYsIV0lgJSFWrA+haOEN2LI0ylKm/vU6
ShQ5iCS9ku7eun2rFAk0AXwD22hZEAGSH3uj7BW47F3Ic1SBdkuTJozRpIozu43hJeEBK0r3vrXA
AgFhkyOaKA4Nt1nOkTvFjRgARvTbFguc6QCwmnfosybT8jwGArNbt7E5zntMUaM9yXjZ+D6+VMDT
ePJutPXNm4zFoUH1jdHC2m3TdMde12VMFW9kYl8UY+bt9seMqfPUOZEQpQV1t1qAgNdHNxW32l7X
uAMPczph5Bx+T0kXxQNmBGZaSRohwBsLornVQg2nBUfMeNOZ54x6gpDkIIbYDdSDWzfsoJlSNwxw
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"""






