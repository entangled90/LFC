#!/bin/bash

./heat | ffmpeg -f rawvideo -pix_fmt rgb24 -s 500x500 -r 60  -vcodec libx264 -vpre libx264 -i - -an -threads 0 -b:v 2100k -f mp4 


#