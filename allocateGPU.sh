#!/bin/zsh


salloc --partition=hzg,allgpu --mem=500G -c 40 --oversubscribe --time=7-00:00:00
