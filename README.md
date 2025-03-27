# bemol

Blade Element Momentum (BEM) academic library.

![bemol logo](images/logo-bemol.svg)

**bemol**: (B)lade (E)lement (M)omentum (O)pen-source (L)ibrary.

From musical notation, bémol means a note a semitone lower. This package is
a semitone lower to commercial or fully packed BEM frameworks, a playground
for students and academics.

## Scope

The aim is to have a very modular basic BEM solver for testing different
models and different aspects of the method. There is no intention to provide
a complete or fast implementation.

### Limitations

- Model is always rigid.
- Blade and rotor forces and moments (thus power) are not calculated, only
  sectional forces and induction factors are calculated.
- Complex geometries (pre-bend or pre-swept) not considered.
- Non-uniform inflow not taken into account.
- Not exhaustively validated.

## Implemented models

### Ning

Implementation of the coupled and uncoupled BEM models of Ning et al.:

Andrew Ning, Gregory Hayman, Rick Damiani and Jason M. Jonkman.
"Development and Validation of a New Blade Element Momentum Skewed-Wake Model
within AeroDyn," AIAA 2015-0215. 33rd Wind Energy Symposium. January 2015. doi:
[doi.org/10.2514/6.2015-0215](10.2514/6.2015-0215), available at
[https://www.nrel.gov/docs/fy15osti/63217.pdf](https://www.nrel.gov/docs/fy15osti/63217.pdf).

## Dependencies

The library is based on [Numpy](https://numpy.org/), [Scipy](https://scipy.org/),
[pandas](https://pandas.pydata.org/) and [PyYAML](https://pyyaml.org/).
[matplotlib](https://matplotlib.org/) is only used in examples, but highly recommended.
[pytest](https://pytest.org/) is as well optional, used for the evaluation
of the tests inside the `tests` sub-folder.

The list of considered/tested versions are available in environment.yaml.
For use with conda:

```
conda env create -f environment.yaml
```