#!/bin/bash

# Enable debugging output by the shell
set -e -x -o pipefail
gatk_meg=$(head -n1 /proc/meminfo | awk '{print int(0.9*($2/1024/1024))}')

# Download inputs
dx-download-all-inputs --parallel

tar zxf  R_qualimap.tar.gz
cp ./R/bin/* /usr/bin/

# Prepare output folder; link all inputs there
mkdir -p out/out/
cd out/out/
ln -sf ~/in/*/* .

if [ -z "${output_prefix}" ]; then
    output_base=`dx describe "${in}" --name`
    output_base="${output_base%.bam}"
else
    output_base="${output_prefix}"
    name=`dx describe "${in}" --name`
    if [ "$name" != "${output_base}.bam" ]; then
        mv $name ${output_base}.bam
    fi
fi

#Command to run
cmd="qualimap $option -bam $output_base.bam --java-mem-size=${gatk_meg}G $extra_cmd"

#GTF file needed for comp-counts and rnaseq options
if [ "$ingtf" != "" ]; then
    if [ "$option" == "comp-counts" ] || [ "$option" == "rnaseq" ]; then
        dx download "$ingtf" -o - | zcat >  ${ingtf_name%.gz}
        cmd="$cmd -gtf ${ingtf_name%.gz}"
    fi
fi

#Specification of the outputs files generated
if [ "$option" == "comp-counts" ]; then
   cmd="$cmd -out ${output_base}.comp-counts.txt"
else
   cmd="$cmd -outfile ${output_base}.pdf -outdir ${output_base}" 
fi 


# Run command
#echo "$cmd -bam $output_base.bam -outfile ${output_base}.pdf -outdir ${output_base} --java-mem-size=${gatk_meg}G"
#eval "$cmd -bam $output_base.bam -outfile ${output_base}.pdf -outdir ${output_base} --java-mem-size=${gatk_meg}G"
echo $cmd
eval "$cmd"
if [ "$option" != "comp-counts" ]; then
   cd ${output_base}
   for file in $(ls *txt); do mv $file ${output_base}_$file; done
   cd ..
fi
if [ "$option" == "comp-counts" ] || [ "$option" == "rnaseq" ]; then
    rm *gtf
fi
# Remove input links from output folder
find . -type l -delete

# Upload outputs
dx-upload-all-outputs --parallel
