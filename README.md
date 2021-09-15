`axi-art`
=========
Procgen line art to plot on an AxiDraw

`crosswords`
============

Every 7x7 crossword grid with the following constraints:

- "American" style with black squares (not bars)
- At least one black square
- All entries are at least 3 letters long
- 180 degree rotational symmetry

This piece was never fully plotted on the AxiDraw; it has *way too many* line segments.

`three_views`
=============

300 randomly sampled points visualized 3 ways. Points sampled with Poisson-disc sampling.

The original idea was that I could come up with new ways of visualizing these points over time, and cycle through which 3 panels I used.

My favorites are `#1, #2, #4`, displayed in that order (left-to-right).

### #1

Points connected with both their minimum-spanning tree as well as their Steiner tree.

### #2

Voronoi diagram of the source points.

### #4

The cells of the Voronoi diagram are 4-colored with the following four styles:

* Blank
* Polygon outline of Voronoi ridge verticies
* Star polygon ({n/2}) of Voronoi ridge verticies
* Line segments connecting each point to all of its Voronoi ridge veritices

Generated With
--------------

    $ python three_views.py -N 300 -p 1 -H 12 -W 8.5 -R 32 -S -r three_views/1.png -o three_views/1.axi
    $ python three_views.py -N 300 -p 2 -H 12 -W 8.5 -R 32 -S -r three_views/2.png -o three_views/2.axi
    $ python three_views.py -N 300 -p 4 -H 12 -W 8.5 -R 32 -S -r three_views/3.png -o three_views/3.axi
