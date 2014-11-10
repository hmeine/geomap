import vigra, numpy, cellimage

def map_from_boundaries(boundaries_filename, cornerType = cellimage.CellType.Line):
    edgeImage = vigra.readImage(boundaries_filename) > 0
    gm = cellimage.GeoMap(edgeImage.astype(numpy.float32), 1, cornerType)
    print boundaries_filename, len(gm.nodes), len(gm.edges), len(gm.faces)
    return gm

def test_simple():
    gm = map_from_boundaries('simple_open_edges.png')
    assert len(gm.edges) == 4
    assert len(gm.faces) == 2
    assert len(gm.nodes) == 9

def test_closed_contours():
    gm = map_from_boundaries('closed_contours.png')
    assert len(gm.edges) == 4
    assert len(gm.faces) == 4
    assert len(gm.nodes) == 3

def test_bridges():
    gm = map_from_boundaries('bridges.png')
    assert len(gm.edges) == 6
    assert len(gm.faces) == 4
    assert len(gm.nodes) == 5

def test_four_connected():
    gm = map_from_boundaries('four_connected.png')
    assert len(gm.edges) == 7
    assert len(gm.faces) == 5
    assert len(gm.nodes) == 5

def test_touching_boundary():
    gm = map_from_boundaries('touching_boundary.png')
    assert len(gm.faces) == 3
    # FIXME: at the moment, there will always be a node at (-1, -1):
    assert len(gm.edges) == 3+1
    assert len(gm.nodes) == 2+1

def test_self_loop():
    gm = map_from_boundaries('self_loop.png')
    assert len(gm.faces) == 3
    assert len(gm.edges) == 2
    assert len(gm.nodes) == 2

def test_regression():
    gm = map_from_boundaries('failing.png')
    assert len(gm.faces) == 3
    assert len(gm.edges) == 2
    assert len(gm.nodes) == 3

