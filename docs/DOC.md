# Sequenzalignierungstool: Documentation
Autoren: 
Al-Juboori, Faisal;
Kandhari, Neha;
Li, Huije

- [File Input](#file-input)
- [FastQC](#fastqc)
- [Cutadapt](#cutadapt)
- [VSearch](#vsearch)
- [Trimming](#trimming)
- [Error Correction](#error-correction)

## File Input
Pfad zur FASTQ Datei wird mit der Funktion `askForFileAndDoesFileExist()` verlangt und es wird auf Exitenz überprüft.

## FastQC
Verwendung von FastQC um mit _Overrepresented Sequences_ eine mögliche Source einer Adaptersequenz einlesen zu können.
`fastqc(input_path::String, output_path::String)`
Nach dieser Funktion wird ein Zip-Archive ausgegeben.
Mit `function read_adapter(fastqc_result_path::String)` wird die Adaptersequenz gelesen.

## Cutadapt
Anschließend wird durch die Ausgabe aus der vorherigen Funktion in die `function cutadapt(Fastq_1, adapt_Seq_Fwd, output_1, Fastq_2=nothing, dapt_Seq_Rev=nothing, output_2=nothing)`als Parameter, zusammen mit dem File-Input und einer selbstbestimmten Output-Datei, übergeben.
Mit dieser Funktion werden mittels Cutadapt die Adaptersequenz aus dem File-Input entfernt.

## VSearch
Die Funktion `function vsearch(input_path::String, output_path::String)` lässt Chimäre entfernen und gibt eine Datei aus, in dem die Chimären gesammelt worden sind.

## Trimming

Mit `function trimming(org_records::Vector{Any}, threshold::Int, min_length::Int)` werden die Reads getrimmt, d.h. es werden schlechte Reads, kurze Reads und Reads mit 'N' entfernt.

## Error Correction

