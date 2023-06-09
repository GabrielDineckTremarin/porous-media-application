python3 video.py
ffmpeg -framerate 12 -start_number 101 -i transient/t%d.png output/animation.mp4
xdg-open output/animation.mp4
