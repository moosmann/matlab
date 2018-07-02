#!/bin/sh

echo -e "\nUpdate repository and overwrite locally modified files (new unstaged files won't be deleted)[y/n]?"
read ans

if [ $ans = y -o $ans = Y -o $ans = yes -o $ans = Yes -o $ans = YES ]
then

echo -e '\ngit fetch origin master'
git fetch origin master

echo -e '\ngit reset --hard origin/master'
git reset --hard origin/master

fi

if [ $ans = n -o $ans = N -o $ans = no -o $ans = No -o $ans = NO ]
then
echo "Not updating."
fi