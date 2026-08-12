#!/bin/zsh

echo -e "\nAllocating node non-exclusively (oversubscribe) of partition ''hzg'' for 7 days using following command:"
echo -e "salloc --partition=hzg --time=7-00:00:00 --oversubscribe --mem=500G"
salloc --partition=hzg --time=7-00:00:00 --oversubscribe 

echo -e "\nLog in to allocated node using:\n ssh -Y USER@NODE"
#ssh -Y  $USER@$SLURM_NODELIST # not working
