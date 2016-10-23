# HartreeFock
A program implementing the Hartree–Fock/self-consistent field method

Description comming soon on http://compphys.go.ro

I'll add here some things, too.

For now I suspect an issue with unrestricted Hartree-Fock, but I committed it anyway, I'll fix it if I'll find the cause.
The restricted method seems to work fine.

Using the classes should be easy. Here is how to grab some atoms from the 'basis':

```c++
  Systems::AtomWithShells H1, H2, O, N, C, He, Li, Ne, Ar;
  // Titanium has Z = 22, could be used to test up to d orbitals

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

	H1.SetCenterForShells();
	H2.SetCenterForShells();
	O.SetCenterForShells();

	Systems::Molecule H2O;

	H2O.atoms.push_back(H1);
	H2O.atoms.push_back(H2);
	H2O.atoms.push_back(O);

	H2O.SetIDs();
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
  Heatom.SetIDs();
```
