# Spatially-structured models of phage infection

The scripts in this repository implement two different models of phage infecting bacteria on a plate. See `latex/notes.pdf` for a complete description of each model.

## Dependencies

The Python scripts were written for Python3. They both require NumPy, a module for efficient numerical operations. Assuming you have Python3 installed, you can install it with `pip3` using the following command:

```
pip3 install numpy
```

For making heatmap animations, matplotlib is required. It can also be installed with `pip3`:

```
pip3 install matplotlib
```

## Converting `matplotlib`-generated PNGs to an animated gif

If you have a directory of numbered PNGs (e.g. 1.png, 2.png, 3.png, etc.), navigate to inside of that directory, and use the following ImageMagick command to convert the PNGs to an animated GIF:

```
convert -delay 15 -loop 0 `ls | sort -n` ../test.gif
```
