#!/bin/bash
###### Part 1 ######
#SBATCH --partition=gpu
#SBATCH --qos=normal
#SBATCH --qos=normal
#SBATCH --account=higgsgpu
#SBATCH --job-name=higgs125charge
#SBATCH --ntasks=16
#SBATCH --output=logs/higgs125charge.log
#SBATCH --mem-per-cpu=24576
#SBATCH --gres=gpu:v100:1

###### Part 2 ######

srun -l hostname

/usr/bin/nvidia-smi -L

echo "Allocate GPU cards : ${CUDA_VISIBLE_DEVICES}"


set -x

source /hpcfs/cepc/higgsgpu/zhuyf/particle_transformer/env.sh

conda activate weaver

echo "args: $@"

# set the dataset dir via `DATADIR_JetClass`
DATADIR=${DATADIR_JetClass}
[[ -z $DATADIR ]] && DATADIR='./datasets/JetClass'

# set a comment via `COMMENT`
suffix=${COMMENT}

# set the number of gpus for DDP training via `DDP_NGPUS`
NGPUS=${DDP_NGPUS}
NGPUS=1
[[ -z $NGPUS ]] && NGPUS=1
if ((NGPUS > 1)); then
    CMD="torchrun --standalone --nnodes=1 --nproc_per_node=$NGPUS $(which weaver) --backend nccl"
else
    CMD="weaver"
fi

epochs=30
samples_per_epoch=$((5*1000 * 1024 ))
samples_per_epoch_val=$((2*1000 * 1024))
dataopts="--num-workers 8 --fetch-step 0.01"

# PN, PFN, PCNN, ParT
model=$1
if [[ "$model" == "ParT" ]]; then
    modelopts="networks/example_ParticleTransformer.py --use-amp"
    batchopts="--batch-size 512 --start-lr 1e-3"
elif [[ "$model" == "PN" ]]; then
    modelopts="networks/example_ParticleNet.py"
    batchopts="--batch-size 512 --start-lr 1e-2"
elif [[ "$model" == "PFN" ]]; then
    modelopts="networks/example_PFN.py"
    batchopts="--batch-size 4096 --start-lr 2e-2"
elif [[ "$model" == "PCNN" ]]; then
    modelopts="networks/example_PCNN.py"
    batchopts="--batch-size 4096 --start-lr 2e-2"
else
    echo "Invalid model $model!"
    exit 1
fi

# "kin", "kinpid", "full"
FEATURE_TYPE=$2
[[ -z ${FEATURE_TYPE} ]] && FEATURE_TYPE="full"

if ! [[ "${FEATURE_TYPE}" =~ ^(full|kin|kinpid)$ ]]; then
    echo "Invalid feature type ${FEATURE_TYPE}!"
    exit 1
fi

# currently only Pythia
SAMPLE_TYPE=Pythia/higgs125Charge

$CMD \
    --data-train \
    "B:${DATADIR}/${SAMPLE_TYPE}/bbtrain/*.root" \
    "C:${DATADIR}/${SAMPLE_TYPE}/cctrain/*.root" \
    "G:${DATADIR}/${SAMPLE_TYPE}/ggtrain/*.root" \
    "Bbar:${DATADIR}/${SAMPLE_TYPE}/bbbartrain/*.root" \
    "Cbar:${DATADIR}/${SAMPLE_TYPE}/ccbartrain/*.root" \
    "D:${DATADIR}/${SAMPLE_TYPE}/ddtrain/*.root" \
    "Dbar:${DATADIR}/${SAMPLE_TYPE}/ddbartrain/*.root" \
    "U:${DATADIR}/${SAMPLE_TYPE}/uutrain/*.root" \
    "Ubar:${DATADIR}/${SAMPLE_TYPE}/uubartrain/*.root" \
    "S:${DATADIR}/${SAMPLE_TYPE}/sstrain/*.root" \
    "Sbar:${DATADIR}/${SAMPLE_TYPE}/ssbartrain/*.root" \
    --data-val "${DATADIR}/${SAMPLE_TYPE}/val/*.root" \
    --data-test \
    "B:${DATADIR}/${SAMPLE_TYPE}/bbtest/*.root" \
    "C:${DATADIR}/${SAMPLE_TYPE}/cctest/*.root" \
    "G:${DATADIR}/${SAMPLE_TYPE}/ggtest/*.root" \
    "Bbar:${DATADIR}/${SAMPLE_TYPE}/bbbartest/*.root" \
    "Cbar:${DATADIR}/${SAMPLE_TYPE}/ccbartest/*.root" \
    "D:${DATADIR}/${SAMPLE_TYPE}/ddtest/*.root" \
    "Dbar:${DATADIR}/${SAMPLE_TYPE}/ddbartest/*.root" \
    "U:${DATADIR}/${SAMPLE_TYPE}/uutest/*.root" \
    "Ubar:${DATADIR}/${SAMPLE_TYPE}/uubartest/*.root" \
    "S:${DATADIR}/${SAMPLE_TYPE}/sstest/*.root" \
    "Sbar:${DATADIR}/${SAMPLE_TYPE}/ssbartest/*.root" \
    --data-config  data/JetClass/JetClass_M11.yaml --network-config $modelopts \
    --model-prefix saved_model/JetClass/${SAMPLE_TYPE}/${FEATURE_TYPE}/${model}3/{auto}${suffix}/net \
    $dataopts $batchopts \
    --samples-per-epoch ${samples_per_epoch} --samples-per-epoch-val ${samples_per_epoch_val} --num-epochs $epochs --gpus '0' \
    --optimizer ranger --log logs/JetClass_${SAMPLE_TYPE}_${FEATURE_TYPE}_${model}3_{auto}${suffix}.log --predict-output higgs125/pred.root \
    --tensorboard JetClass_${SAMPLE_TYPE}_${FEATURE_TYPE}_${model}3${suffix} \
    "${@:3}"
