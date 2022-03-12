# Sequenzalignierungstool
Autoren: 
Al-Juboori, Faisal;
Kandhari, Neha;
Li, Huije

- [Kurzbeschreibung](#kurzbeschreibung)
- [Verwendung](#verwendung)

## Kurzbeschreibung
Sequenzalignierungstool (für RNA-Seq-Data) mit Referenzgenom inklusive Read-Aufbereitung (trimmen, Adaptoren / Chimären entfernen, Fehlerkorrektur) und anschließend Ausgabe der Scaffolds.

## Benötigte Software:

1. Julia 
2. Cutadapt
3. Fastqc
4. Vsearch
5. Hisat2
6. samtools

# Verwendung
1. Open Julia
2. Change to working dir. ,Where your .jl File is using ``` cd("Path to .jl File") ```
3. Now Excute the .jl file using ``` include(yourFile.jl) ```
4. follow the instructions (Input is your Fastq File)
5. finished 
