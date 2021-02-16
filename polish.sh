#!/bin/bash
#echo "$@"
map_q=0
base_q=13
threads=16
gpu_id=-1
coverage=40
config_file="./brnn/config/wtdbg2+racon.3x3.yaml"
model_file="./brnn/ctc.ecoli_yeast_chr21.wtdbg2+racon.layer3x3.model/model.chkpt"

usage(){
    echo "usage: ./polish.sh -i <draft_assembly> -r <basecalled_reads> -o <output_dir> [options ...]"
    echo "options:"
    echo "  -i STR   input file in FASTA format, containing draft assembly which will be polished"
    echo "  -r STR   input file in FASTA/FASTQ format, containing basecalled raw reads used for polishing"
    echo "  -o STR   output directory of files"
    echo "  -q INT   skip alignments with mapQ smaller than INT, default: 0"
    echo "  -Q INT   skip bases with baseQ/BAQ smaller than INT, default: 13"
    echo "  -t INT   number of threads, default: 16"
    echo "  -g INT   index of gpu device, default: None"
    echo "  -c INT   coverage of reads used for polishing, default: 40"
    echo "  -f STR   neural network config file, default: ./brnn/config/wtdbg2+racon.3x3.yaml"
    echo "  -m STR   neural network pre-trained file, default: ./brnn/ctc.ecoli_yeast_chr21.wtdbg2+racon.layer3x3.model/model.chkpt"
}

opt_cnt=0
while getopts ":i:r:o:q:Q:t:g:c:f:m:h" opt; do
    case ${opt} in
        i) draft_assembly="${OPTARG}" && ((opt_cnt++)) 
        ;;
        r) basecalled_reads="${OPTARG}" && ((opt_cnt++)) 
        ;;
        o) output_dir="${OPTARG}" && ((opt_cnt++))
        ;;
        q) map_q="${OPTARG}" && ((opt_cnt++)) 
        ;;
        Q) base_q="${OPTARG}" && ((opt_cnt++)) 
        ;;
        t) threads="${OPTARG}" && ((opt_cnt++)) 
        ;;
        g) gpu_id="${OPTARG}" && ((opt_cnt++)) 
        ;;
        c) coverage="${OPTARG}" && ((opt_cnt++))
        ;;
        f) config_file="${OPTARG}" && ((opt_cnt++))
        ;;
        m) model_file="${OPTARG}" && ((opt_cnt++))
        ;;
        h) usage && ((opt_cnt++)) && exit 1
    esac
done


if [ ${opt_cnt} -lt 3 ]; then
    usage && exit 1
fi

mkdir ${output_dir}

output_dir=`readlink -f ${output_dir}`
config_file=`readlink -f ${config_file}`
model_file=`readlink -f ${model_file}`

echo $config_file
echo $model_file


echo '---------- Racon prepolishing ------------'
minimap2 -x map-ont ${draft_assembly} ${basecalled_reads} -t ${threads} > ${output_dir}/aln.paf
racon ${basecalled_reads} ${output_dir}/aln.paf ${draft_assembly} -t ${threads} > ${output_dir}/assembly_racon.fasta

echo ''
echo '---------- Aligning ------------'
echo ''
minimap2 -ax map-ont ${output_dir}/assembly_racon.fasta ${basecalled_reads} -t ${threads} > ${output_dir}/aln.sam

echo ''
echo '---------- Filtering alignments ------------'
echo ''
python brnn/scan_filter_alignment.py --input ${output_dir}/aln.sam --output ${output_dir}/filter_aln.sam --max_indels 5000

echo ''
echo '---------- Samtools ------------'
echo ''
samtools view -bS -@ ${threads} ${output_dir}/filter_aln.sam -o ${output_dir}/filter_aln.bam
samtools sort -@ ${threads} ${output_dir}/filter_aln.bam -o ${output_dir}/filter_aln.sort.bam
samtools index -@ ${threads} ${output_dir}/filter_aln.sort.bam

echo ''
echo '---------- Generating features ------------'
echo ''
cd neuralpolishextract/
python run_samtools.py -r ${output_dir}/assembly_racon.fasta -b ${output_dir}/filter_aln.sort.bam -o ${output_dir}/pileups -t ${threads} -bq ${base_q} -mq ${map_q}
./extract -p ${output_dir}/pileups -w 64 -t ${threads} -o ${output_dir}/features -c ${coverage}

echo ''
echo '---------- Network polishing ------------'
echo ''
cd ../brnn
#CUDA_VISIBLE_DEVICES=${gpu_id} python ctc_polish.py -config config/human_guppy.yaml -model_path human.guppy.wtdbg2+racon.layer3x3.model/model.chkpt -data_path ${output_dir} -prepolished_path ${output_dir}/assembly_racon.fasta -output ${output_dir}/polished.fasta
if [ ${gpu_id} -ne -1 ]; then
    CUDA_VISIBLE_DEVICES=${gpu_id} python ctc_polish.py -config ${config_file} -model_path ${model_file} -data_path ${output_dir} -prepolished_path ${output_dir}/assembly_racon.fasta -output ${output_dir}/polished.fasta
else
    python ctc_polish.py -config ${config_file} -model_path ${model_file} -data_path ${output_dir} -prepolished_path ${output_dir}/assembly_racon.fasta -output ${output_dir}/polished.fasta
fi
cd ..

echo ''
echo '---------- Deleting files ---------------'
echo ''
cd ${output_dir}
rm -rf features/ pileups/ aln.sam aln.paf assembly_racon.fasta assembly_racon.fasta.fai filter_aln.sam filter_aln.bam filter_aln.sort.bam filter_aln.sort.bam.bai

echo 'polishing done.'
 

