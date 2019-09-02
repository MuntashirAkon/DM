#!/bin/bash

path='../ADACT/exec/'
os=$(uname -s)
file='dm'

if ! [ -e "$path" ]; then
  >&2 echo "ADACT path doesn't exists!"
  exit 1
fi

if [ "$os" == "Darwin" ]; then
  file='dm_mac'
fi

if ! [ -e "./dm" ]; then
  >&2 echo "dm not found. Did you build it yet?"
  exit 1
fi

mv "./dm" "$path/$file"
if [ $? -eq 0 ]; then
  echo "Success!"
fi