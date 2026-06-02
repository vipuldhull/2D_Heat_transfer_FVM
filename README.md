# 2D Finite Volume Heat Conduction Solver

A Python implementation of a 2D steady-state heat conduction solver using the Finite Volume Method (FVM).

Developed as part of the Computational Thermofluid Dynamics course at the Technical University of Munich.

## Features

- Structured body-fitted mesh generation
- Finite Volume discretization
- Dirichlet boundary conditions
- Neumann boundary conditions
- Robin boundary conditions
- Arbitrary domain geometries
- 2D contour visualization
- 3D temperature visualization

## Governing Equation

The solver computes the steady-state heat conduction equation:

∇·(k∇T)=0

using a cell-centered finite volume formulation.

## Domain Shapes

Implemented geometries:

- Rectangular
- Linear
- Quadratic
- Crazy geometry

## Example Results

Insert screenshots here.

## How to Run

```bash
pip install -r requirements.txt
jupyter notebook
