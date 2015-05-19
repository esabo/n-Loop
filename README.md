# n-Loop
Computation of the Dimofte-Garoufalidis n-loop invariant.

For a more complete description of this project, please visit http://people.math.gatech.edu/~esabo3, or more specifically, http://people.math.gatech.edu/~esabo3/knot-data.html.

GitHub does not support folders with a large number of files, so the data, pre-computed Neumann-Zagier datum, and SageMath object file results need to be downloaded separately under the "release" header. Human readable results may be found http://people.math.gatech.edu/~esabo3/knot-data.html, while text data file results may be found http://people.math.gatech.edu/~stavros/publications/nloop.data/index.html.

A brief usage guide is given here, while a more detailed explanation may be found in http://arxiv.org/abs/1503.02554. Note that the following code assumes SnapPy is installed inside of SageMath. The Manifold class used in the code refers to SnapPy and not SageManifolds.

- COMPUTATION OF EXACT N-LOOP INVARIANTS FROM MANIFOLD
```python
attach('nloop_exact.py')
all_diagrams = load('6diagrams.sobj')
M = Manifold('6_2')
nloop_from_manifold(M, 2, all_diagrams, engine="retrieve")
```

- COMPUTATION OF EXACT N-LOOP INVARIANTS FROM PRECOMPUTED NZ-DATA
```python
attach('nloop_exact.py')
all_diagrams = load('6diagrams.sobj')
nz62 = load('nzdata/nz_exact_K5_19.sobj')
nloop_from_nzdatum(nz62, 2, all_diagrams)
```

- MANUAL COMPUTATION OF EXACT N-LOOP INVARIANTS FROM MANIFOLD
```python
attach('nloop_exact.py')
all_diagrams = load('6diagrams.sobj')
M = Manifold('6_2')
D = NeumannZagierDatum(M, engine="retrieve")
D.generate_nz_data()
E = nloop(D.nz, 2, all_diagrams)
E.nloop_invariant()
E.one_loop()
```

- GENERATE EXACT NZ-DATA
```python
attach('nloop_exact.py')
M = Manifold('6_2')
D = NeumannZagierDatum(M, engine="retrieve")
D.generate_nz_data()
D.nz
```
