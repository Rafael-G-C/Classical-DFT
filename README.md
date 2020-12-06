## A python code based on the article by Jeanmairet et al.(https://doi.org/10.1021/ed500049m)

on lines

```python
distances = np.linspace(1,10,100)# angstroms, starting distance, ending distance, points between the start and the end
densities = np.linspace(0.01,0.1,100)# molecules/angstrom^3 starting density, ending density, points between the start and the end
```
In _distances_ we can change the maximum distance and how many points we want to calculate while in _densities_ we can change how many densities should the numerical aproximation try for the minimization

to keep things as simple as possible the code outputs the data

```
distances,numerical_density,analytical density
1.0,0.01,0.0
1.0909090909090908,0.01,0.0
1.1818181818181819,0.01,0.0
1.2727272727272727,0.01,0.0
1.3636363636363638,0.01,0.0
1.4545454545454546,0.01,0.0
1.5454545454545454,0.01,0.0
```
