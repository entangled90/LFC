#!/bin/bash

./schrodinger3D | ffmpeg  -f rawvideo -pix_fmt rgb24 -s 640x480 -r 30 -i - -an -vcodec libx264 -preset libx264-max -threads 0 -b:v 2100k -f mp4 foo.mp4 


#
