
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
