# 3 views of a set of points
# Generate points with even distribution?
# Minimum spanning tree or Steiner tree?
# Voronoi boundaries or Denaulay triangularization
# Partition graph by 4-coloring; overlay 4 delaunay triangularizations

import axi
import axi.paths
import numpy
import random
import shapely.geometry
import subprocess
from scipy.spatial import Voronoi, Delaunay
import networkx
import claspy as cy

def distance((x1, y1), (x2, y2)):
    return ((x1 - x2) ** 2. + (y1 - y2) ** 2.) ** 0.5

def generate_points(n=256, min_distance=None, k=30):
    # Ref: Robert Bridson, "Fast Poisson Disk Sampling in Arbitrary Dimensions"
    # Bastardized O(N^2) version for simplicity, I only need ~300 points

    # This seems to approximately yield `n` points
    if min_distance is None:
        min_distance = (n ** -0.5) / 1.25

    xs = [(random.random(), random.random())]
    active_list = [xs[0]]

    while active_list:
        idx = random.randrange(len(active_list))
        x = active_list[idx]
        for i in range(k):
            while True:
                point = (random.random(), random.random())
                if min_distance < distance(point, x) <= 2.0 * min_distance:
                    break
            if all((distance(point, q) > min_distance for q in xs)):
                xs.append(point)
                active_list.append(point)
                break
        else:
            active_list.pop(idx)
    return xs

def steiner_points(points):
    # Get list of steiner points
    # The min spanning tree of (points + steiner points) is the steiner tree of (points)
    # XXX: I think this could theoretically deadlock, but shouldn't in practice
    efst = subprocess.Popen(["./geosteiner-3.1/efst"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    bb = subprocess.Popen(["./geosteiner-3.1/bb"], stdin=efst.stdout, stdout=subprocess.PIPE)
    efst.stdout.close()

    for x, y in points:
        efst.stdin.write("{} {}\n".format(x, y))
    efst.stdin.close()
    
    sps = []
    for line in bb.stdout:
        if line.startswith(" % @C"):
            row = line.split()
            sps.append((float(row[2]), float(row[3])))

    return sps

def mst(points):
    # min spanning tree
    # O(N^2) naive greedy
    paths = []
    mst_points = [points.pop()]
    while points:
        # Find min edge from node in points to note in mst_points
        min_distance = None
        min_edge = None
        for to_point in mst_points:
            for from_point in points:
                dist = distance(from_point, to_point)
                if min_distance is None or dist < min_distance:
                    min_distance = dist
                    min_edge = from_point, to_point
        from_point, to_point = min_edge
        paths.append([from_point, to_point])
        points.remove(from_point)
        mst_points.append(from_point)
    return paths

def voronoi_finite_polygons_2d(vor, radius=2.):
    """
    https://stackoverflow.com/questions/20515554/colorize-voronoi-diagram
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= -0.5 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= numpy.linalg.norm(t)
            n = numpy.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = numpy.sign(numpy.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = numpy.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = numpy.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = numpy.array(new_region)[numpy.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, numpy.asarray(new_vertices)

def voronoi(points, wtf=True):
    vor = Voronoi(points)
    regions, vertices = voronoi_finite_polygons_2d(vor)
    paths = []
    for region in regions:
        polygon = vertices[region].tolist()
        paths.append(polygon + [polygon[0]])
    return paths

def colored_voronoi(vor):
    edges = {}

    colors = [cy.MultiVar(*range(4)) for p in vor.points]

    # Special case filtering
    for c, pr in zip(colors, vor.point_region):
        region = vor.regions[pr]
        if any(r < -0.5 for r in region):
            cy.require(c != 0)
        if len(region) <= 4:
            cy.require(c != 3)

    for ridge_pts in vor.ridge_points:
        left, right = ridge_pts
        cy.require(colors[left] != colors[right])

    cy.solve()
    print colors
    d = {i: c.value() for i, c in enumerate(colors)}

    subgraphs = [[] for i in range(4)]
    for i, c in enumerate(colors):
        subgraphs[c.value()].append(vor.points[i].tolist())

    raw_colors = [c.value() for c in colors]

    return subgraphs, raw_colors

def delaunay(points):
    tri = Delaunay(points)
    paths = []
    for simp in tri.simplices:
        paths.append(tri.points[simp])
    return paths

# 

def panel_one(points):
    sps = steiner_points(points)
    steiner_tree = mst(sps + points)
    plain_tree = mst(points)
    d = axi.Drawing(steiner_tree + plain_tree)
    return d

def panel_two(points, ratio):
    paths = voronoi(points)
    bbox = shapely.geometry.box(0, 0, 1, ratio)
    spath = axi.paths.paths_to_shapely(paths)
    spath = spath.intersection(bbox)
    paths = axi.paths.shapely_to_paths(spath)
    d = axi.Drawing(paths)
    return d

def panel_three(points):
    vor = Voronoi(points)
    subgraphs, _ = colored_voronoi(vor)

    paths = []
    for sg in subgraphs:
        paths.extend(delaunay(sg))
        #paths.extend(mst(sg))
        #sps = steiner_points(sg)
        #paths.extend(mst(sps + sg))

    d = axi.Drawing(paths)
    return d


def panel_four(points):
    # 4-colored voronoi diagram, then draw each color with a different style
    paths = []
    def _a(center, vertices):
        pass
    def _b(center, vertices):
        paths.append(vertices + [vertices[0]])
    def _c(center, vertices):
        paths.extend([[center, v] for v in vertices])
    def _d(center, vertices):
        paths.extend([[a, b] for a, b in zip(vertices, vertices[2:] + vertices)])
    draw_fns = [_a, _b, _c, _d]

    vor = Voronoi(points)
    subgraphs, colors = colored_voronoi(vor)
    regions, vertices = vor.regions, vor.vertices

    for point, pr, c in zip(vor.points, vor.point_region, colors):
        region = regions[pr]
        while -1 in region:
            region.remove(-1)
        polygon = vertices[region].tolist()
        draw_fns[c](point, polygon)

    bbox = shapely.geometry.box(0, 0, 1, ratio)
    spath = axi.paths.paths_to_shapely(paths)
    spath = spath.intersection(bbox)
    paths = axi.paths.shapely_to_paths(spath)
    d = axi.Drawing(paths)
    return d

def panel_five(points):
    # MSTs of voronoi cell centers with the same number of edges
    vor = Voronoi(points)

    sides = {}
    for point, pr in zip(vor.points, vor.point_region):
        region = vor.regions[pr]
        l = len(region)
        sides.setdefault(l, []).append(point.tolist())

    paths = []
    for n, points in sides.items():
        print n, len(points)
        paths.extend(mst(points))
    
    d = axi.Drawing(paths)
    return d

def panel_six(points):
    # Flow, with points pushing flow away

    paths = []
    for x in numpy.linspace(0, 1, 101):
        path = []
        for y in numpy.linspace(0, 1, 2001):
            cx, cy = x, y
            for px, py in points:
                dx = px - x
                dy = py - y
                r = (dx ** 2 + dy ** 2) ** 0.5
                k = 1e-10
                scale = k / (r ** 5)
                scale = -min(scale, 0.01)
                cx += scale * dx / r
                cy += scale * dy / r
            path.append((cx, cy))
            #print (x, y), (cx, cy)
        paths.append(path)

    paths.append([
        (0, 0),
        (1, 0),
        (1, 1),
        (0, 1),
        (0, 0),
    ])

    d = axi.Drawing(paths)
    return d

#
def find_seed_with_n_points(n, ratio, target, start=0):
    for i in range(start, 1000):
        random.seed(i)
        points = generate_points(n=n/ratio)
        points = [p for p in points if p[1] < ratio]
        print i, len(points)
        if len(points) != target:
            continue
        print "done"
        return points

if __name__ == "__main__":
    cli = axi.cli()
    cli.add_argument("-N", "--points", type=int, help="Approximate number of points to generate", default=256.)
    cli.add_argument("-p", "--panel", help="panel to generate")
    cli.add_argument("-R", "--random", type=int, help="random seed", default=0)
    args = cli.parse_args()

    ratio = min(args.width / args.height, args.height / args.width)

    #random.seed(args.random)
    #points = generate_points(n=args.points / ratio)
    #points = [p for p in points if p[1] < ratio]

    points = find_seed_with_n_points(args.points, ratio, args.points, start=args.random)
    print "generated %d points" % len(points)

    if args.panel == "1":
        d = panel_one(points)
    elif args.panel == "2":
        d = panel_two(points, ratio=ratio)
    elif args.panel == "3":
        d = panel_three(points)
    elif args.panel == "4":
        d = panel_four(points)
    elif args.panel == "5":
        d = panel_five(points)
    elif args.panel == "6":
        d = panel_six(points)
    else:
        cli.error("invalid panel %s" % args.panel)

    cli.draw(d)
