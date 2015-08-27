Additional atomic parameters for nMOLDYN (by Pellegrini Eric -- pellegrini@ill.fr).

You will find in Database\Atoms some extra atomic parameters (i.e. coherent and incoherent 
scattering lengths) that are not provided by the standard version of MMTK.
Almost all the periodical elements classification has been registered. When available the 
datas for isotopes are also included.


The source for the coherent and incoherent scattering lengths is:

	Neutron Data Booklet
	Edited by Albert-José Dianoux and Gerry Lander
	ILL


If during some file conversions or trajectory readings this kind or error message occurs:


	raise IOError("Database entry %s/%s not found" % (directory, filename))
	IOError: Database entry Atoms/os not found


this may indicate that some atoms of your system are not recognized by MMTK (Osmium atom, os, in 
the above case).

To solve this kind of trouble, just copy the file corresponding to the atom symbol that is missing
(i.e. for osmium Database\Atoms\os) and paste it in:

	$YOUR_PYTHON_PATH\Lib\site-packages\MMTK\Database\Atoms\

Where $YOUR_PYTHON_PATH is the path where your version of python has been installed.
Everything should work find afterwards.

Have fun