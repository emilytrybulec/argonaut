process OUTPUT_COMBINE {
    label 'process_low'

    input:
    path(input_files)

    output:
    path("all_assemblyStats.tsv")       , emit: assemblyStats
   
    script: 
    def prefix
    """
    file_array_1=(\$(ls ${input_files}))

    ls ${input_files} > file_names.txt

    # List the files and their line lengths

    awk '{ print length(), \$0 | "sort -n" }' file_names.txt > sorted_files_int.txt
    awk '{print \$2}' sorted_files_int.txt  > sorted_files.txt

    # Extract only the filenames from the sorted file
    file_array_2=(\$(cat sorted_files.txt))

    output_file="\${file_array_2[0]}"

    # Loop through the rest of the files and join them one by one
    for file in "\${file_array_2[@]:1}"; do
        tmp_file=\$(mktemp)
        join -t \$'\\t' "\$output_file" "\$file" > "\$tmp_file"
        mv "\$tmp_file" "\$output_file"
    done

    # Optionally save the final output to a file
    mv "\$output_file" all_assemblyStats.tsv

    """
}
