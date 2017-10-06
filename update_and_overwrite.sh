#!/bin/sh

echo -e '\nUpdate repository:'

echo -e '\ngit fetch origin master'
git fetch origin master

echo -e '\ngit reset --hard origin/master'
git reset --hard origin/master'
