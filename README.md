\# Composite Material Properties Calculator



This project was developed in July 2021 as part of an engineering project at École Polytechnique de Thiès.



It is a Python tool for calculating mechanical properties, stresses, and strength criteria of composite laminates using Classical Lamination Theory.



\## Features



\- Calculation of ply mechanical properties

\- Laminate stiffness matrices A, B, D and ABD

\- Mid-plane strains and curvatures

\- Stress calculation in each ply

\- Stress transformation between global and principal material axes

\- Tsai-Hill failure criterion

\- Input data read from an Excel configuration file



\## Project structure



```text

main.py              Main execution file

Configuration.py     Reads laminate configuration from Excel

Ply.py               Defines ply mechanical properties

Laminate.py          Performs laminate calculations

PlyStrength.py       Computes ply strength values

materials.py         Defines fiber and matrix materials

settings.xlsx        Input data file

docs/                Project report

## Requirements

This project requires Python 3 and the following Python libraries:

- numpy
- pandas
- sympy
- openpyxl

You can install them with:

```bash
pip install -r requirements.txt

## How to run

Make sure that `settings.xlsx` is in the same folder as the Python files.

Then run:

```bash
python main.py

## Author

Moussa Kantibo DIEDHIOU
