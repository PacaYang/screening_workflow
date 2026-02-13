The tree should look like this
```
External FSx
/shared/
├───── /task_1
├───── /task_2
└───── /task_3
       ├── README
       ├── /Input
       │    ├── compounds_smiles.csv
       │    ├── sequences.csv
       │    └── /protein_file
       │         ├── /protein1
       │         └── /protein2
       │              ├── pdb file
       │              ├── gro file
       │              ├── itp file
       │              └── system_EM.top
       ├── /protein1
       ├── …
       └── /proteinN
            ├── summary.csv
            ├── /inital_screen
            │    ├── input.csv
            │    ├── /GraphDTA
            │    ├── /HMSA
            │    ├── /ColdDTA
            │    ├── summary.csv
            │    └── selected.csv
            └── /fine_screen
                 ├── /PBSA
                 │    ├── /DiffDock
                 │    │    ├── /output 1
                 │    │    └── /output n
                 │    │
                 │    ├── /PBSA
                 │    │    ├── /MD
                 │    │    └── /PBSA
                 │    └── summary.csv
                 ├── /Vina
                 │    ├── /output
                 │    └── summary.csv
                 ├── /AF3
                 │    ├── /output
                 │    └── summary.csv
                 └── /Boltz2
                      ├── /output
                      └── summary.csv
```
