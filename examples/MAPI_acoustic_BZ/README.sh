
# VMD render with Tacyon; you turned off occlusion as your laptop wouldn't take it :^)
ffmpeg -framerate 30 -i conf.%05d.ppm -s:v 1280x720 -c:v libx264 -profile:v high -crf 23 -pix_fmt yuv420p -r 30 toney.mp4

# A bit too small
ffmpeg -i toney.mp4 -vf scale=300:-2 small-toney.mp4

# A bit crummy, and 10x bigger than the MP4
ffmpeg  -i toney.mp4 -filter_complex "[0:v] split [a][b];[a] palettegen [p];[b][p] paletteuse" toney.gif

