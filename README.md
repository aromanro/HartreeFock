# HartreeFock
A program implementing the Hartree–Fock/self-consistent field method (and more, see the README at the end)

[![CodeFactor](https://www.codefactor.io/repository/github/aromanro/hartreefock/badge)](https://www.codefactor.io/repository/github/aromanro/hartreefock)

Description is available here: http://compphys.go.ro/the-hartree-fock-program/

Some Hartree-Fock theory here: http://compphys.go.ro/the-hartree-fock-method/

Some things more general (Schrödinger equation, Born-Oppenheimer approximation, variational principle), here: http://compphys.go.ro/how-to-solve-a-quantum-many-body-problem/

### PROGRAM IN ACTION

[![Program video](https://img.youtube.com/vi/tmCRJxOAIH8/0.jpg)](https://youtu.be/tmCRJxOAIH8)

### HOW TO COMPUTE

Here are some things about usage:

Using the classes should be easy. Here is how to grab some atoms from the 'basis':

```c++
	Systems::AtomWithShells H1, H2, O, N, C, He, Li, Ne, Ar;

	for (auto &atom : basis.atoms)
	{
		if (atom.Z == 1) H1 = H2 = atom;
		else if (atom.Z == 2) He = atom;
		else if (atom.Z == 3) Li = atom;
		else if (atom.Z == 8) O = atom;
		else if (atom.Z == 6) C = atom;
		else if (atom.Z == 7) N = atom;
		else if (atom.Z == 10) Ne = atom;
		else if (atom.Z == 18) Ar = atom;
	}
```

Here is how to set the H2O molecule with the coordinates from the 'Mathematica Journal' (referenced in the code):

```c++
	H1.position.X = H2.position.X = O.position.X = 0;

	H1.position.Y = 1.43233673;
	H1.position.Z = -0.96104039;
	H2.position.Y = -1.43233673;
	H2.position.Z = -0.96104039;

	O.position.Y = 0;
	O.position.Z = 0.24026010;

	Systems::Molecule H2O;

	H2O.atoms.push_back(H1);
	H2O.atoms.push_back(H2);
	H2O.atoms.push_back(O);

	H2O.Init();
  ```
  
  And here is how you calculate:
  
  ```c++
  
  HartreeFock::RestrictedHartreeFock HartreeFockAlgorithm;
  HartreeFockAlgorithm.alpha = 0.5;
  HartreeFockAlgorithm.initGuess = 0;
  
  HartreeFockAlgorithm.Init(&H2O);
  double result = HartreeFockAlgorithm.Calculate();
  ```
  
You can do computation for a single atom, too, for now by putting it into a dummy molecule with a single atom in it. For example for He:

```c++
  Systems::Molecule Heatom;
  Heatom.atoms.push_back(He);
  Heatom.Init();
```

### LATEST ADDITIONS

I added more basis sets, it's not limited to STO-nG. 
While doing that, I found and fixed a bug in the integrals repository that manifested for orbitals having L > 1.
Now it seems to work fine.

Basis sets added:
- Split valence orbitals: 3-21G, 6-21G, 6-31G
- Besides 'split valence', polarization on heavy atoms: 6-31G*
- 'Split valence', polarization on heavy atoms and hydrogen and diffusion functions on heavy atoms: 6-31+G**

You may add more of them, the parsed format is nwchem.

### Tutorials from Crawford Group

Most of those are already implemented:

https://github.com/CrawfordGroup/ProgrammingProjects

Projects for MP2, DIIS, CCSD(T), CIS and TDHF/RPA are implemented.
DIIS and MP2 are implemented for both the restricted and unrestricted methods.

Projects 1 and 2 are not implemented here, maybe I'll add a python workbook for those in the python repository, 3 was already implemented obviously before finding this.
Project 7 is also not implemented, also 9 and 13.

