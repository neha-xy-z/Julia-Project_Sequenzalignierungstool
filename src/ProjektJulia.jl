# Needed Modules
using FASTX
using ZipFile

# Set Path For Linux
function get_path_variable(file_path::String)

    if Sys.islinux()
        file_fastx = string(split(file_path, "/")[end])     # xxx.fastq
        file_name = string(split(file_fastx, ".")[1])       # xxx
        work_path = string(replace(file_path, "/"*file_fastx => ""))    # xx/xx/xx
    else
        file_fastx = string(split(file_path, "\\")[end])
        file_name = string(split(file_fastx, ".")[1])
        work_path = string(replace(file_path, "\\"*file_fastx => ""))
    end
    return file_fastx, file_name, work_path
end

# Set Path to used Softwares for Windows   
function setPathToSoftwareForWindows()
    curr_dir = pwd()
    curr_dir = replace(curr_dir, "\\src" => "")

    path_to_cutadapt = curr_dir*"\\cutadapt"
    path_to_fastqc = curr_dir*"\\fastqc\\FastQC"
    path_to_vsearch = curr_dir*"\\vsearch\\vsearch\\bin"
   
    return path_to_cutadapt, path_to_fastqc, path_to_vsearch
end

path_to_cutadapt, path_to_fastqc, path_to_vsearch = setPathToSoftwareForWindows()

# FastQC Unix Function
function fastqcUnix(input_file::String)
    run(`fastqc $input_file`)
end

# FastQC Windows Function
function fastqcWindows(input_fastqc_path::String, input_path::String, output_path::String)
    cd(input_fastqc_path)
    run(`java -Xmx250m -classpath .";"./sam-1.103.jar";"./jbzip2-0.9.jar uk.ac.babraham.FastQC.FastQCApplication $input_path`)
    curr_dir = replace(input_fastqc_path, "\\fastqc\\FastQC" => "")
    cd(curr_dir*"\\src")
end

# Cutadapt Run for Linux
function cutadaptUnix(adapter_seqs_str::String, input_file::String)
    run(`cutadapt -j 0 $adapter_seqs_str -o without_adapters.fastq $input_file`)
end

# Cutadapt Run for Windows
function cutadaptWindows(input_cutadapt_path::String, adapter_seqs_str::String, file::String)
    cd(input_cutadapt_path)
    run(`cutadapt-3.4.exe $adapter_seqs_str -o without_adapters.fastq $file`) # Windows
    curr_dir = replace(input_cutadapt_path, "\\cutadapt" => "")
    cd(curr_dir*"\\src")
end 

# Function for reading adapter from FastQC File
function read_adapter(fastqc_result_path::String)
            
    global fastqc_res = ZipFile.Reader(fastqc_result_path)

    flag = false
    seqs = String[]
    counts = String[]
    percentages = String[]
    sources = String[]
    adapter_seq_str = ""

    for file in fastqc_res.files
        if contains(file.name, "fastqc_data.txt")
            for line in eachline(file)
                if contains(line, ">>Overrepresented")
                    flag = true
                    continue
                end

                if flag == true && contains(line, ">>END_MODULE")
                    flag = false
                    break
                end

                if flag && !contains(line, "#")
                    seq, count, percentage, source = split(line, "\t")
                    if string(source) != "No Hit"
                        push!(seqs, seq)
                        push!(counts, count)
                        push!(percentages, percentage)
                        push!(sources, source)
                    end
                end
            end
        end
    end 

    if length(seqs) == 0
        #TODO
        adapter_seq_str = "no adapter seq"
    else

        for adapter_seq in seqs
            adapter_seq_str *= "-a $adapter_seq "
        end
    end

    return seqs, sources, string(strip(adapter_seq_str))
end

# Read adapterfree fastq File and store into records
function readFastq(input_file::String)
    # println("Reading file without adapters...")
    reader = FASTQ.Reader(open(input_file, "r"))
    records = []
    for record in reader
        push!(records, record)
    end
    close(reader)
    return records
end

# Read adapterfree fastq File and store into records
function readFasta(input_file::String)
    reader = FASTA.Reader(open(input_file, "r"))
    records = []
    for record in reader
        push!(records, record)
    end
    close(reader)
    return records
end

# Collect garbage and free memory 
function free_memory(variable::Any)
    empty!(variable)
    GC.gc()
end

# Get threshold parameter for the trimming function
function getCheckParameter(parameter::String, min::Int, max::Int)

    while true
        print("\n Enter Parameter for: ", parameter, "\nMinimum: ", min, ", Maximum: ", max, "\n>>> ")

        input_num = parse(Int, string(strip(readline())))

        if input_num >= min && input_num <= max
            return input_num
        else
            print("You can only enter a threshold size of 3, 4 or 5. Please enter again.")
        end
    end
end

# Get threshold parameter for the trimming function (Overwrite)
function getCheckParameter(parameter::String, min::Int)

    while true
        print("\n Enter Parameter for: ", parameter, "\nMinimum: ", min, "\n>>> ")

        input_num = parse(Int, string(strip(readline())))

        if input_num >= min
            return input_num
        else
            print("You can only enter a threshold size of 3, 4 or 5. Please enter again.")
        end
    end
end

# Trimming
function trimming(input_records::Vector{Any}, threshold::Int=3, min_read_length::Int=30)

    file = FASTQ.Writer(open("trimmed.fastq", "w"), true)
    
    for record in input_records
        quality_list = quality(record)
        seq_list = sequence(String, record)
        count = 0
        position = 0

        # 'N' und short-reads skip
        if length(seq_list) < min_read_length || contains(seq_list, 'N')
            continue
        end
        
        for base_index = 1 : length(quality_list)    
            if Int64(quality_list[base_index]) <= 20
                count += 1
                if count == threshold
                    position = base_index
                    break
                end
            end
        end

        # low quality reads 
        if count == threshold 

            # only work with long-reads 
            if (position-threshold) >= min_read_length

                # New Sequence with new Quality
                new_seq = sequence(String, record, (1:position-threshold))
                new_quality = quality(record, 33, (1:position-threshold))
                # Get the same identifiers and description from the old Record and give it to the new Record
                ident = identifier(record)
                desc = description(record)
                new_record = FASTQ.Record(ident, desc, new_seq, new_quality)
                write(file, new_record)
            end
        else
            write(file, record) 
        end
    end
    free_memory(input_records)
    close(file)
end

# VSearch Run for different OS Linux, Windows
function vsearchrun(input_vsearch_path::String, input_path::String)
    if Sys.islinux()
        run(`vsearch --uchime3_denovo trimmed.fastq --nonchimeras nonchimeras.fasta`)
    else
        # println(pwd())
        cd(input_vsearch_path)
        run(`vsearch --uchime3_denovo $input_path --nonchimeras nonchimeras.fasta`)
        curr_dir = replace(input_vsearch_path, "\\vsearch\\vsearch\\bin" => "")
        cd(curr_dir*"\\src")
    end
end

# Ask for kmer size
function kmer_input()
    
    println("\nEnter a kmer-size:\n")
    n = parse(UInt8, readline())
    if n != UInt8
        println("\nWrong input.\n")
        kmer_input()
    else
        return n
    end
end

# Count frequent kmers
function kmers_Histo(k, records)
    kmer_dict = Dict{String,Integer}()

    for record in records
        rec_seq = string(FASTA.sequence(record))

        for i = 1:((length(rec_seq) - k) + 1)
            sequence_k = rec_seq[i:((i+k)-1)]
            if !haskey(kmer_dict, sequence_k)
                kmer_dict[sequence_k] = 1
            else
                kmer_dict[sequence_k] += 1
            end
        end
    end
    return kmer_dict
end

# Find neighbours from a kmer
function hammingDistance(kmer)

    alphabet = ['A','C','G','T']
    neighbors = []

    for all_pos = reverse(1:length(kmer))

        for character in alphabet

            new_kmer = collect(kmer)

            if character != new_kmer[all_pos]
                new_kmer[all_pos] = character
                push!(neighbors, join(new_kmer))
            end
        end

    end

    return neighbors
end

# Correct Reads into free-error Reads
function newReads(kmersHisto, recs, t, k)

    file = FASTA.Writer(open("new_errorfree_reads.fasta", "w"), 1000)
    
    for rec in recs

        seq = string(FASTA.sequence(rec))
        n_Kmers = (length(seq) - k) + 1

        for i = 1:n_Kmers

            kmer = seq[i:(i+k)-1]

            if haskey(kmersHisto, kmer)

                if kmersHisto[kmer] < t

                    all_neighbors = hammingDistance(kmer)

                    for nk in all_neighbors
                        if haskey(kmersHisto, nk)

                            if kmersHisto[nk] >= t

                                first_part = seq[1:i-1]
                                kmer_part = nk
                                second_part = seq[k+i:end]
    
                                new_seq = string(first_part, kmer_part, second_part)
                                id = identifier(rec)
    
                                new_record = FASTA.Record(id,new_seq)
                                write(file, new_record)

                                break
                            end
                        end
                    end
                    free_memory(all_neighbors)
                end
            end
        end
    end
    close(file)
end 

# Hisat2 build index
function hisat2_build_Unix(input_file::String)
    run(`hisat2-build -p 4 $input_file index`)
end

# Hisat2 Alignment
function hisat2Unix()
    run(`hisat2 -p 4 --dta -x index -f nonchimeras.fasta -S alignment`)
end

# Hisat2 Alignment (Overwrite) 
function hisat2Unix(input_index::String)
    run(`hisat2 -p 4 --dta -x $input_index -f nonchimeras.fasta -S alignment.sam`)
end

# Samtool Output to Scaffolds
function samtoolsUnix()
    run(`samtools sort alignment.sam -o sorted_alignment.bam`)
    run(`samtools index sorted_alignment.bam`)
    write("scaffolds.txt", read(`samtools idxstats sorted_alignment.bam`))
end

# FastQC and Cutadapt Run for different OS Linux, Windows
function communication()
    org_fastq_file_path = ""
    work_path = ""
    file_name = ""
    with_adapter = true

    # input file 
    while true
        print("\nEnter absolute path of your FASTQ file. Or type 'q' to quit.\n>>> ")

        input = string(strip(readline()))
        if input == "q"
            println("Bye.")
            exit()
        elseif isfile(input)
            if endswith(input, ".fastq")
                org_fastq_file_path = input
                break
            else
                println("ERROR:\nFile input is not a FASTQ File. Check your path again.")
            end
        else
            println("ERROR:\nFile not found. Check your path again and enter on request.")
        end
    end

    # Get Path
    file_name, work_path = get_path_variable(org_fastq_file_path)

    # Open Path Directory
    cd(work_path)

    # fastqc and cutadapt
    while true 
        print("\nDo you want to run FastQC to get a possible Source for an Adapter Sequence, which could be found in your input-file?
                 \nIf yes, type 'y'.
                 \nIf no, type 'n'. Then you can enter an individual Adapter sequence afterwards on request.
                 \nOr type 'q' to quit.
                 \n>>> ")

        input = readline()

        if input == "q"
            println("Bye.")
            exit()

        elseif input == "y"
        
            if Sys.islinux()
                # FastQC Run for OS Linux
                fastqcUnix(file_name*".fastq")
            else 
                # FastQC Run for Windows
                #fastqcWindows(path_to_fastqc, org_fastq_file_path, path_to_fastqc)
            end

            println("\nFastQC successful.\n")

            println("\nReading adapter sequences...\n")

            if Sys.islinux()
                adapter_seqs, adapter_sources, adapter_seq_str = read_adapter(file_name*"_fastqc.zip")
            else
                adapter_seqs, adapter_sources, adapter_seq_str = read_adapter(work_path*"\\"*file_name*"_fastqc.zip")
            end
            
            println("\nAdapters found:\n")

            if startswith(adapter_seq_str, "no")
                with_adapter = false
                println("\t"*adapter_seq_str)
                break
            else
                for i = 1 : length(adapter_seqs)
                    println(adapter_seqs[i], adapter_sources[i])
                end

                # Run Cutadapt
                if Sys.islinux()
                    cutadaptUnix(adapter_seq_str, file_name*".fastq")
                    break
                else
                    #  Windows
                    #cutadaptWindows(path_to_cutadapt, adapter_seq_str, work_path*"\\"*file_name*".fastq")
                    break
                end
            end

        elseif input == "n"
            print("\nPlease type in your adapter sequence:\n>>> ")

            adapter_seq = string(strip(readline()))

            if adapter_seq == ""
                with_adapter = false
                #TODO
                println("no adapter seq")
                break
            elseif contains(adapter_seq, r"[^atcgATCG]")
                #TODO
                println("\n warning check again")
            else
                # Run Cutadapt
                if Sys.islinux()
                    cutadaptUnix("-a "*adapter_seq, file_name*".fastq")
                    break
                else
                    cutadaptWindows(path_to_cutadapt, adapter_seq_str, work_path*"\\"*file_name*".fastq")
                    break
                end
            end

        else
            #TODO
            println("ERROR:\nWrong input.")
        end
    end

    if Sys.islinux()
        # Read fastq after cutadapt
        if with_adapter
            records = readFastq("without_adapters.fastq")
        else
            records = readFastq(file_name*".fastq")
        end
    else
        if with_adapter
            # find file without adapters 
            file_without_adapters_file = path_to_cutadapt*"\\without_adapters.fastq"
    
            # read fastq into records
            records = readFastq(file_without_adapters_file)
        else
            records = readFastq(work_path*"\\"*file_name*".fastq")
        end
    end

    # input: threshould and min_length
    trimm_threshold = getCheckParameter("threshold", 3, 5)
    trimm_min_length = getCheckParameter("Minimum length", 30)

    # Trimming
    trimming(records, trimm_threshold, trimm_min_length)

    # vsearch and trimming
    if Sys.islinux()
        vsearchUnix()
    else
        curr_dir = pwd()
        fileForVSearch = curr_dir*"\\trimmed.fastq"
        vsearchrun(path_to_vsearch, fileForVSearch)
    end

    # read fasta
    if Sys.islinux()
        records = readFasta("nonchimeras.fasta")
    else
        # find vsearch file and read into records
        vsearch_file = path_to_vsearch*"\\nonchimeras.fasta"
        records = readFasta(vsearch_file)
    end

    # input: kmer
    kmer = getCheckParameter("Kmer", 8, 55)

    kmers_Histo_Dict = kmers_Histo(kmer, records)

    # input: newRead threshold
    newReads_threshold = getCheckParameter("newRead threshold", 6)

    newReads(kmers_Histo_Dict, records, newReads_threshold, kmer)

    #### Only for Linux Operating System

    if Sys.islinux()
        # Hisat2 1. Step: build index,  2. Step Alignment
        while true
            println("\nDo you have indexed files for the required genome? Type 'y', if yes. Else 'n' for no.\n Or type 'q' to quit.")
            input = string(strip(readline()))
            if input == "q"
                exit()
            elseif input == "n"  # ohne index

                println("\nType in the absoulte path of your genome FASTA-file.")

                input_file = string(strip(readline())) #hg38.fa

                hisat2_build_Unix(input_file)
                hisat2Unix()
                break

            elseif input == "y"
                println("\nType in the absolute path for the index-file, which should only contain the prefix without extension.\n (Example: 'your/path/to/indexFiles/genome')")

                input_file = string(strip(readline())) # index

                hisat2Unix(input_file)
                break
            else
                println("Wrong input.")
            end
        end
        # samtools
        samtoolsUnix()
    else
        println("FASTA File has been succesfully created and saved.")
    end

end

communication()