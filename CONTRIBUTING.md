# Contributing to atommovr

Interested in building out *atommovr* or adding a custom feature? You've come to the right place. 

*atommovr* was intentionally designed to be customizable to every user's needs. Instead of being sucked into the black hole of trying to code up *every imaginable feature and use case* (and inevitably failing), we tried to keep *atommovr* very modular so that users could easily add and/or amend features. Because after all, who knows your use case better than *you*?

In particular, the `Algorithm` and `ErrorModel` classes are very general, and their source files describe the requirements for implementing new algorithms or error models.

## Specific opportunities for contribution

*Code maintenance and performance*
- Adding more unit tests to `atommover.tests/`

*New features*
- Extending framework to support general array/lattice shapes
  - i.e. by building out the `ArrayGeometry` class in `atommover.utils.core.py` and building support for animations in `atommover.utils.animation.py`
- Adding automatic plotting to the benchmarking module
  - i.e. by building out the `BenchmarkingFigure` class in `atommover.utils.benchmarking.py`

*Growing the library*
- Adding more algorithms (see `atommover.algorithms.Algorithm.py` for a template)
- Adding more error models (see `atommover.utils.ErrorModel.py` for a template)

Thinking of something that's not on this list? Feel free to contact [Nikhil](mailto:nikhil.harle@colorado.edu).