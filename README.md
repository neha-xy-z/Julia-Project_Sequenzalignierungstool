# Sequenzalignierungstool
Für das B.Sc. Modul BI2024 Einführung in die Programmiersprache Julia WiSe21/22 <br>

- [Kurzbeschreibung](#kurzbeschreibung)
- [Benötigte Software](#benötigte-software)
- [Verwendung](#verwendung)

## Kurzbeschreibung
Sequenzalignierungstool (für RNA-Seq-Data) mit Referenzgenom inklusive Read-Aufbereitung (trimmen, Adaptoren / Chimären entfernen, Fehlerkorrektur) und anschließend Ausgabe der Scaffolds.

#### Screenshot
<img src="Screenshot.png" width="700" height="600">

## Benötigte Software

1. Julia 
2. Cutadapt
3. FastQC
4. VSearch
5. Hisat2
6. Samtools

Für Linux müssen diese über `sudo apt-get <Software>` vorinstalliert werden.<br>
Für Windows muss dieser Ordner heruntergeladen werden. <br>
**Hinweis**: Da Windows kein Hisat2 und Samtools verwenden kann, wird auch keine Alignment stattfinden. Die Ausführung mit Windows endet mit der Readaufbereitung.

_P.S.: Wir haben keinen STAR Aligner verwendet, weil diese Software **30 GB RAM** braucht!_

## Verwendung
1. Öffne Julia.
2. Wechsele ins Working Directory, wo sich auch die ProjektJulia.jl befindet. <br>
` cd("Path to .jl File") `
3. Führe nun die .jl-Datei mit aus `include("ProjektJulia.jl")`
4. Befolge nun die Anweisungen.

Autoren: Neha Kandhari, Faisal Al-Juboori, Huijie Li
